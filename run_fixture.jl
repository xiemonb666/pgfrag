using Peridynamics
using LinearAlgebra
using TOML

function require(condition, message)
    condition || error(message)
end

function find_root!(parent::Vector{Int}, value::Int)
    root = value
    while parent[root] != root
        parent[root] = parent[parent[root]]
        root = parent[root]
    end
    while parent[value] != value
        next_value = parent[value]
        parent[value] = root
        value = next_value
    end
    return root
end

function join_roots!(parent::Vector{Int}, i::Int, j::Int)
    ri = find_root!(parent, i)
    rj = find_root!(parent, j)
    ri == rj || (parent[max(ri, rj)] = min(ri, rj))
    return nothing
end

function write_rows(path::String, header::String, rows)
    open(path, "w") do io
        println(io, header)
        for row in rows
            println(io, join(string.(row), ","))
        end
    end
end

function run_fixture(args)
    require(length(args) == 17, "expected 17 runner arguments")
    outdir = args[1]
    dims_m = parse.(Float64, args[2:4]) .* 1e-3
    dx_m = parse(Float64, args[5]) * 1e-3
    horizon_m = parse(Float64, args[6]) * 1e-3
    material_name = args[7]
    rho = parse(Float64, args[8])
    young = parse(Float64, args[9])
    nu = parse(Float64, args[10])
    gc = parse(Float64, args[11])
    load_kind = args[12]
    speed = parse(Float64, args[13])
    requested_steps = parse(Int, args[14])
    dt_arg = args[15]
    safety = parse(Float64, args[16])
    max_points = parse(Int, args[17])

    require(material_name in ("BB", "OSB"), "unsupported material model")
    require(load_kind in ("zero_load", "symmetric_axial_velocity"), "unsupported loading kind")
    require(all(>(0), dims_m) && dx_m > 0 && horizon_m > dx_m, "invalid geometry or horizon")
    require(rho > 0 && young > 0 && 0 < nu < 0.5 && gc > 0, "invalid material parameters")
    require(material_name != "BB" || isapprox(nu, 0.25; atol=1e-12), "BBMaterial in this 3D backend requires nu=0.25")
    require(requested_steps > 0 && 0 < safety < 1 && max_points > 0, "invalid solver/resource settings")

    mkpath(outdir)
    position, volume = uniform_box(dims_m[1], dims_m[2], dims_m[3], dx_m)
    n = size(position, 2)
    require(0 < n <= max_points, "generated point count exceeds max_points")
    require(all(isfinite, position) && all(isfinite, volume) && all(>(0), volume), "fixture point set invalid")
    material = material_name == "BB" ? BBMaterial() : OSBMaterial()
    body = Body(material, position, volume)
    zmin = minimum(@view position[3, :])
    zmax = maximum(@view position[3, :])
    tol = dx_m * 1e-8
    bottom_ids = findall(abs.(@view(position[3, :]) .- zmin) .<= tol)
    top_ids = findall(abs.(@view(position[3, :]) .- zmax) .<= tol)
    require(!isempty(bottom_ids) && !isempty(top_ids) && minimum(top_ids) > maximum(bottom_ids),
            "fixture must have distinct top and bottom loading layers")
    point_set!(body, :bottom, bottom_ids)
    point_set!(body, :top, top_ids)
    material!(body; horizon=horizon_m, rho=rho, E=young, nu=nu, Gc=gc)
    load_speed = load_kind == "zero_load" ? 0.0 : speed
    velocity_bc!(t -> -load_speed, body, :bottom, 3)
    velocity_bc!(t -> load_speed, body, :top, 3)

    solver = dt_arg == "AUTO" ?
             VelocityVerlet(steps=requested_steps, safety_factor=safety) :
             VelocityVerlet(steps=requested_steps, stepsize=parse(Float64, dt_arg), safety_factor=safety)
    # Preflight the exact native bound for this fixture before integration.
    # The one-chunk exporter does not claim this is a general OSB/contact guarantee.
    stability_dh = Peridynamics.threads_data_handler(body, solver, Threads.nthreads())
    native_stability_bound = Peridynamics.calc_stable_timestep(stability_dh, Float64(safety))
    require(isfinite(native_stability_bound) && native_stability_bound > 0,
            "STOP_UNVERIFIED_DT: upstream native stable timestep is not finite and positive")
    requested_dt = dt_arg == "AUTO" ? nothing : parse(Float64, dt_arg)
    dt_tolerance = 128 * eps(Float64) * max(native_stability_bound, requested_dt === nothing ? 0.0 : requested_dt)
    if requested_dt !== nothing
        require(isfinite(requested_dt) && requested_dt > 0, "STOP_UNVERIFIED_DT: explicit dt must be finite and positive")
        require(requested_dt <= native_stability_bound + dt_tolerance,
                "STOP_UNVERIFIED_DT: explicit dt $(requested_dt) exceeds upstream native bound $(native_stability_bound)")
    end
    job = Job(body, solver)
    wall_time_s = @elapsed dh = submit(job; quiet=true)
    require(length(dh.chunks) == 1, "local_smoke requires a single CPU solver chunk")
    if requested_dt === nothing
        require(isapprox(solver.Δt, native_stability_bound; rtol=1e-12, atol=dt_tolerance),
                "STOP_UNVERIFIED_DT: auto-accepted dt does not match native calc_stable_timestep")
    else
        require(isapprox(solver.Δt, requested_dt; rtol=1e-12, atol=dt_tolerance),
                "STOP_UNVERIFIED_DT: accepted dt does not match explicit request")
    end
    chunk = dh.chunks[1]
    storage = chunk.storage
    system = chunk.system
    point_ids = system.chunk_handler.point_ids
    require(length(point_ids) == n && length(unique(point_ids)) == n && Set(point_ids) == Set(1:n),
            "solver returned incomplete or duplicate global point IDs")
    now = storage.position
    velocity = storage.velocity
    damage = storage.damage
    active_count = storage.n_active_bonds
    sed = Peridynamics.export_field(Val(:strain_energy_density), chunk.mat, system, storage,
                                    chunk.paramsetup, solver.n_steps * solver.Δt)
    require(length(damage) == n && length(sed) == n, "unexpected state field length")
    require(all(isfinite, damage) && all(isfinite, sed) && all(>=(0), sed), "non-finite state or strain energy")

    directed = Dict{Tuple{Int,Int},Vector{Bool}}()
    for local_i in 1:system.chunk_handler.n_loc_points
        global_i = point_ids[local_i]
        for bond_id in Peridynamics.each_bond_idx(system, local_i)
            bond = system.bonds[bond_id]
            global_j = point_ids[bond.neighbor]
            require(global_i != global_j, "self-bond detected")
            edge = global_i < global_j ? (global_i, global_j) : (global_j, global_i)
            push!(get!(directed, edge, Bool[]), storage.bond_active[bond_id])
        end
    end
    for (edge, states) in directed
        require(length(states) == 2, "expected exactly two directed entries for bond $edge")
        require(states[1] == states[2], "reciprocal bond states disagree for $edge")
    end
    edges = sort!(collect(keys(directed)))
    broken = [edge for edge in edges if !directed[edge][1]]

    initial_parent = collect(1:n)
    for (i, j) in edges
        join_roots!(initial_parent, i, j)
    end
    initial_roots = unique(find_root!(initial_parent, i) for i in 1:n)
    intact_parent = collect(1:n)
    for (i, j) in edges
        directed[(i, j)][1] && join_roots!(intact_parent, i, j)
    end
    roots = [find_root!(intact_parent, global_id) for global_id in point_ids]
    ordered_roots = sort!(unique(roots))
    root_to_component = Dict(root => component for (component, root) in enumerate(ordered_roots))
    point_component = [root_to_component[root] for root in roots]
    component_count = length(ordered_roots)
    component_sizes = [count(==(component), point_component) for component in 1:component_count]
    component_masses = zeros(Float64, component_count)
    for local_i in 1:n
        global_i = point_ids[local_i]
        component_masses[point_component[local_i]] += rho * volume[global_i]
    end
    active_degree = zeros(Int, n)
    for (i, j) in edges
        if directed[(i, j)][1]
            active_degree[i] += 1
            active_degree[j] += 1
        end
    end
    isolated_ids = findall(==(0), active_degree)

    point_rows = Vector{Tuple}(undef, n)
    for local_i in 1:n
        global_i = point_ids[local_i]
        point_rows[local_i] = (
            global_i, position[1, global_i], position[2, global_i], position[3, global_i],
            now[1, local_i], now[2, local_i], now[3, local_i],
            velocity[1, local_i], velocity[2, local_i], velocity[3, local_i],
            damage[local_i], active_count[local_i], sed[local_i], volume[global_i], rho * volume[global_i],
            point_component[local_i],
        )
    end
    sort!(point_rows; by=row -> row[1])
    write_rows(joinpath(outdir, "points.csv"),
               "point_id,x_ref_m,y_ref_m,z_ref_m,x_current_m,y_current_m,z_current_m,vx_mps,vy_mps,vz_mps,damage,n_active_bonds,strain_energy_density_J_m3,volume_m3,mass_kg,intact_bond_component_id",
               point_rows)

    bond_rows = [(i, j, directed[(i, j)][1], norm(system.position[:, system.chunk_handler.localizer[i]] - system.position[:, system.chunk_handler.localizer[j]]),
                  norm(now[:, system.chunk_handler.localizer[i]] - now[:, system.chunk_handler.localizer[j]])) for (i, j) in edges]
    write_rows(joinpath(outdir, "bond_states.csv"), "point_i_global,point_j_global,active,reference_length_m,current_length_m", bond_rows)
    broken_rows = filter(row -> !row[3], bond_rows)
    write_rows(joinpath(outdir, "broken_bonds.csv"), "point_i_global,point_j_global,reference_length_m,current_length_m",
               [(row[1], row[2], row[4], row[5]) for row in broken_rows])
    component_rows = [(point_ids[i], point_component[i], rho * volume[point_ids[i]], damage[i]) for i in 1:n]
    write_rows(joinpath(outdir, "components.csv"), "point_id,component_id,point_mass_kg,damage", component_rows)
    isolated_rows = []
    for global_i in isolated_ids
        local_i = system.chunk_handler.localizer[global_i]
        push!(isolated_rows, (global_i, system.position[1, local_i], system.position[2, local_i],
                              system.position[3, local_i], rho * volume[global_i]))
    end
    write_rows(joinpath(outdir, "isolated_points.csv"), "point_id,x_ref_m,y_ref_m,z_ref_m,mass_kg", isolated_rows)

    point_center_span_mm = [1e3 * (maximum(@view(position[k, :])) - minimum(@view(position[k, :]))) for k in 1:3]
    full_support_mm = point_center_span_mm .+ (dx_m * 1e3)
    volume_sum_mm3 = sum(volume) * 1e9
    mass_sum_kg = rho * sum(volume)
    kinetic_energy_j = 0.0
    for local_i in 1:n
        global_i = point_ids[local_i]
        kinetic_energy_j += 0.5 * rho * volume[global_i] * sum(abs2, @view velocity[:, local_i])
    end
    strain_energy_j = sum(sed .* volume[point_ids])
    time_s = solver.n_steps * solver.Δt
    summary = Dict{String,Any}(
        "status" => "COMPLETED",
        "solver_model" => "Peridynamics.jl 0.5.4+pgfrag.1 $(material_name)Material / CriticalStretch",
        "backend_id" => "peridynamics_jl_headless_0_5_4",
        "backend_version" => string(Base.pkgversion(Peridynamics)),
        "units" => "SI internally; fixture request dimensions and spacing are mm",
        "calibration_status" => "UNVALIDATED_PROTOTYPE",
        "rho_kg_m3" => rho, "E_pa" => young, "nu" => nu, "Gc_j_m2" => gc,
        "dimensions_requested_mm" => parse.(Float64, args[2:4]),
        "dx_m" => dx_m, "horizon_m" => horizon_m,
        "n_points" => n, "point_center_span_mm" => point_center_span_mm,
        "implied_full_point_cell_support_mm" => full_support_mm,
        "point_volume_sum_mm3" => volume_sum_mm3,
        "point_mass_kg" => mass_sum_kg,
        "load_kind" => load_kind, "speed_m_s" => speed,
        "n_steps" => solver.n_steps, "dt_s" => solver.Δt, "time_s" => time_s,
        "requested_dt_mode" => requested_dt === nothing ? "AUTO" : "EXPLICIT",
        "native_stability_bound_dt_s" => native_stability_bound,
        "accepted_dt_s" => solver.Δt,
        "native_stability_check" => requested_dt === nothing ? "PASS_NATIVE_AUTO_DT" : "PASS_EXPLICIT_DT_WITHIN_NATIVE_BOUND",
        "native_stability_scope" => "single-body fixture only; not a general OSB/contact stability guarantee",
        "max_damage" => maximum(damage), "mean_damage" => sum(damage) / length(damage),
        "strain_energy_j" => strain_energy_j, "kinetic_energy_j" => kinetic_energy_j,
        "initial_graph_component_count" => length(initial_roots),
        "intact_bond_connected_component_count" => component_count,
        "component_sizes_points" => component_sizes, "component_masses_kg" => component_masses,
        "canonical_undirected_bond_count" => length(edges), "broken_undirected_bond_count" => length(broken),
        "isolated_point_ids" => isolated_ids,
        "reciprocal_bond_state_consistent" => true,
        "all_global_point_ids_complete_unique" => true,
        "total_energy_conservation" => "UNKNOWN: imposed boundary velocities perform external work; no complete ledger",
        "contact" => "UNKNOWN: contact not exercised",
        "surface_mapping" => "UNKNOWN: point fixture has no geometry surface reconstruction",
        "wall_time_s" => wall_time_s,
    )
    open(joinpath(outdir, "summary.toml"), "w") do io
        TOML.print(io, summary)
    end
    println("PGFRAG_SOLVER_COMPLETED points=$n bonds=$(length(edges)) broken=$(length(broken)) components=$component_count time_s=$time_s")
    return nothing
end

try
    run_fixture(ARGS)
catch err
    showerror(stderr, err, catch_backtrace())
    println(stderr)
    exit(1)
end
