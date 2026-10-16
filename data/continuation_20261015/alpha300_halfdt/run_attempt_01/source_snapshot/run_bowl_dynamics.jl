using Peridynamics
using LinearAlgebra
using TOML

function require(condition, message)
    condition || error(message)
end

"""Small type-specific adapter: preserve the native Body API but filter the analytic cavity."""
struct GeometryFilteredBody{M<:Peridynamics.AbstractMaterial} <: Peridynamics.AbstractBody{M}
    inner::Peridynamics.Body{M}
    center_m::NTuple{3,Float64}
    inner_radius_m::Float64
    outer_radius_m::Float64
    excluded::Dict{Tuple{Int,Int},Float64}
end

function GeometryFilteredBody(inner::Peridynamics.Body{M}, center_m, inner_radius_m, outer_radius_m) where {M}
    return GeometryFilteredBody{M}(inner, Tuple(Float64.(center_m)), Float64(inner_radius_m),
                                   Float64(outer_radius_m), Dict{Tuple{Int,Int},Float64}())
end

function Base.getproperty(body::GeometryFilteredBody, name::Symbol)
    name in (:inner, :center_m, :inner_radius_m, :outer_radius_m, :excluded) && return getfield(body, name)
    return getproperty(getfield(body, :inner), name)
end

function Base.propertynames(body::GeometryFilteredBody, private::Bool=false)
    native = propertynames(getfield(body, :inner), private)
    return Tuple(unique((native..., :inner, :center_m, :inner_radius_m, :outer_radius_m, :excluded)))
end

function segment_minimum_radius(a, b, center)
    dx = b .- a
    d2 = dot(dx, dx)
    t = d2 == 0.0 ? 0.0 : clamp(-dot(a .- center, dx) / d2, 0.0, 1.0)
    return norm(a .+ t .* dx .- center)
end

function Peridynamics.find_bonds(body::GeometryFilteredBody, loc_points::AbstractVector{Int})
    bonds, native_neighbors = Peridynamics.find_bonds(body.inner, loc_points)
    ranges = Peridynamics.find_bond_ids(native_neighbors)
    keep = trues(length(bonds))
    neighbors = copy(native_neighbors)
    center = collect(body.center_m)
    for (local_i, point_i) in enumerate(loc_points)
        for bond_id in ranges[local_i]
            bond = bonds[bond_id]
            radius = segment_minimum_radius(@view(body.position[:, point_i]),
                                            @view(body.position[:, bond.neighbor]), center)
            if radius < body.inner_radius_m
                keep[bond_id] = false
                neighbors[local_i] -= 1
                i, j = minmax(point_i, bond.neighbor)
                old = get(body.excluded, (i, j), Inf)
                body.excluded[(i, j)] = min(old, radius)
            end
        end
    end
    any(keep) || error("cavity filter removed every native bond")
    return bonds[keep], neighbors
end

function Peridynamics.pre_submission_check(body::GeometryFilteredBody,
                                            solver::Peridynamics.AbstractTimeSolver;
                                            body_in_multibody_setup::Bool=false)
    return Peridynamics.pre_submission_check(body.inner, solver;
                                              body_in_multibody_setup=body_in_multibody_setup)
end

Peridynamics.parameter_setup_type(body::GeometryFilteredBody, spec::Peridynamics.AbstractParamSpec) =
    Peridynamics.parameter_setup_type(body.inner, spec)

function write_rows(path, header, rows)
    open(path, "w") do io
        println(io, header)
        for row in rows
            println(io, join(string.(row), ","))
        end
    end
end

function json_escape(s::AbstractString)
    io = IOBuffer()
    for c in s
        if c == '"'
            print(io, "\\\"")
        elseif c == '\\'
            print(io, "\\\\")
        elseif c == '\n'
            print(io, "\\n")
        elseif c == '\r'
            print(io, "\\r")
        elseif c == '\t'
            print(io, "\\t")
        else
            print(io, c)
        end
    end
    return String(take!(io))
end

function json_value(x)
    if x isa AbstractString || x isa Symbol
        return "\"$(json_escape(string(x)))\""
    elseif x isa Bool
        return x ? "true" : "false"
    elseif x === nothing
        return "null"
    elseif x isa Real
        isfinite(x) || error("non-finite value cannot be written to JSON")
        return repr(x)
    elseif x isa AbstractDict
        return "{" * join((json_value(string(k)) * ":" * json_value(v) for (k, v) in x), ",") * "}"
    elseif x isa Tuple || x isa AbstractVector
        return "[" * join(json_value.(collect(x)), ",") * "]"
    end
    error("unsupported JSON value type $(typeof(x))")
end

function write_json(path, object)
    open(path, "w") do io
        println(io, json_value(object))
    end
end

function read_points(path)
    lines = readlines(path)
    length(lines) >= 2 || error("point CSV is empty")
    position = Matrix{Float64}(undef, 3, length(lines) - 1)
    for (i, line) in enumerate(lines[2:end])
        fields = split(strip(line), ',')
        length(fields) == 3 || error("point CSV row $i must have three coordinates")
        position[:, i] .= parse.(Float64, fields)
    end
    all(isfinite, position) || error("point CSV includes non-finite coordinates")
    return position
end

function write_little_endian(path, values::Array{T}) where {T<:Union{Float64,Float32,Int32,UInt8}}
    Base.ENDIAN_BOM == 0x04030201 || error("native result writer requires a little-endian host")
    open(path, "w") do io
        write(io, reinterpret(UInt8, vec(values)))
    end
end

function stable_root!(parent::Vector{Int}, i::Int)
    root = i
    while parent[root] != root
        parent[root] = parent[parent[root]]
        root = parent[root]
    end
    while parent[i] != i
        next = parent[i]
        parent[i] = root
        i = next
    end
    return root
end

function component_labels(n_points, edges, active)
    parent = collect(1:n_points)
    for bond_index in eachindex(edges)
        active[bond_index] || continue
        i, j = edges[bond_index]
        ri, rj = stable_root!(parent, i), stable_root!(parent, j)
        ri == rj || (parent[max(ri, rj)] = min(ri, rj))
    end
    return Int32[stable_root!(parent, i) - 1 for i in 1:n_points]
end

function pack_active(active::Vector{Bool})
    packed = zeros(UInt8, cld(length(active), 8))
    for k in eachindex(active)
        active[k] && (packed[(k - 1) ÷ 8 + 1] |= UInt8(1) << ((k - 1) % 8))
    end
    return packed
end

function state_active(system, storage, point_ids, edge_to_bond_ids, edges)
    active = Vector{Bool}(undef, length(edges))
    for k in eachindex(edges)
        ids = edge_to_bond_ids[edges[k]]
        length(ids) == 2 || error("canonical pair $(edges[k]) lacks two directed native bonds")
        a, b = storage.bond_active[ids[1]], storage.bond_active[ids[2]]
        a == b || error("reciprocal active flags disagree for pair $(edges[k])")
        active[k] = a
    end
    return active
end

function save_state!(state_index, position_history, velocity_history, damage_history,
                     component_history, active_history, sed_history, chunk, edges,
                     edge_to_bond_ids, dt)
    system, storage = chunk.system, chunk.storage
    point_ids = system.chunk_handler.point_ids
    n_points = length(point_ids)
    for local_i in eachindex(point_ids)
        global_i = point_ids[local_i]
        for d in 1:3
            position_history[state_index, global_i, d] = storage.position[d, local_i]
            velocity_history[state_index, global_i, d] = storage.velocity[d, local_i]
        end
        damage_history[state_index, global_i] = Float32(storage.damage[local_i])
    end
    active = state_active(system, storage, point_ids, edge_to_bond_ids, edges)
    labels = component_labels(n_points, edges, active)
    component_history[state_index, :] .= labels
    active_history[state_index, :] .= pack_active(active)
    sed = Peridynamics.export_field(Val(:strain_energy_density), chunk.mat, system,
                                    storage, chunk.paramsetup, (state_index - 1) * dt)
    for local_i in eachindex(point_ids)
        sed_history[state_index, point_ids[local_i]] = Float32(sed[local_i])
    end
    return active
end

function compute_state_ledger(chunk, edge_pairs, active, rest_length, epsilon_c,
                              point_mass, point_volume, bond_constant, time_s)
    system, storage = chunk.system, chunk.storage
    point_ids = system.chunk_handler.point_ids
    x = storage.position
    v = storage.velocity
    kinetic = 0.0
    momentum = zeros(3)
    com = zeros(3)
    for local_i in eachindex(point_ids)
        global_i = point_ids[local_i]
        m = point_mass[global_i]
        kinetic += 0.5 * m * sum(abs2, @view v[:, local_i])
        momentum .+= m .* @view(v[:, local_i])
        com .+= m .* @view(x[:, local_i])
    end
    com ./= sum(point_mass)
    angular = zeros(3)
    force = zeros(3)
    for local_i in eachindex(point_ids)
        global_i = point_ids[local_i]
        p_i = point_mass[global_i] .* @view(v[:, local_i])
        angular .+= cross(@view(x[:, local_i]) .- com, p_i)
        force .+= point_volume[global_i] .* @view(storage.b_int[:, local_i])
    end
    pair_energy = 0.0
    for k in eachindex(edge_pairs)
        active[k] || continue
        i, j = edge_pairs[k]
        li = system.chunk_handler.localizer[i]
        lj = system.chunk_handler.localizer[j]
        stretch = (norm(@view(x[:, lj]) .- @view(x[:, li])) - rest_length[k]) / rest_length[k]
        pair_energy += 0.5 * bond_constant * point_volume[i] * point_volume[j] * rest_length[k] * stretch^2
    end
    native_sed = Peridynamics.export_field(Val(:strain_energy_density), chunk.mat, system,
                                           storage, chunk.paramsetup, time_s)
    native_energy = sum(native_sed[local_i] * system.volume[local_i] for local_i in eachindex(point_ids))
    return (kinetic, pair_energy, native_energy, force, momentum, angular)
end

function setup_velocity!(body, point_matrix, volume, rho, alpha)
    masses = rho .* volume
    total_mass = sum(masses)
    com = zeros(3)
    for i in eachindex(masses)
        com .+= masses[i] .* @view(point_matrix[:, i])
    end
    com ./= total_mass
    A = Diagonal([alpha, -0.25 * alpha, 0.0])
    velocity = Matrix{Float64}(undef, 3, length(masses))
    for i in eachindex(masses)
        velocity[:, i] = A * (@view(point_matrix[:, i]) .- com)
    end
    mean_velocity = zeros(3)
    for i in eachindex(masses)
        mean_velocity .+= masses[i] .* @view(velocity[:, i])
    end
    mean_velocity ./= total_mass
    velocity .-= mean_velocity
    inertia = zeros(3, 3)
    angular = zeros(3)
    for i in eachindex(masses)
        r = @view(point_matrix[:, i]) .- com
        p = masses[i] .* @view(velocity[:, i])
        inertia .+= masses[i] .* (dot(r, r) .* Matrix{Float64}(I, 3, 3) .- r * transpose(r))
        angular .+= cross(r, p)
    end
    omega = inertia \ angular
    for i in eachindex(masses)
        r = @view(point_matrix[:, i]) .- com
        velocity[:, i] .-= cross(omega, r)
    end
    # Native ICs retain the exact pointwise field, with the point function argument named p.
    for d in 1:3
        f = p -> begin
            r = p .- com
            (A * r .- mean_velocity .- cross(omega, r))[d]
        end
        Peridynamics.velocity_ic!(f, body, :all_points, d)
    end
    return (com, mean_velocity, omega, velocity)
end

function run_bowl(case_path, outdir)
    Threads.nthreads() == 1 || error("native bowl run requires exactly one Julia thread")
    mkpath(outdir)
    case = TOML.parsefile(case_path)
    case["schema"] == "pgfrag.native_bowl_case.v1" || error("case schema mismatch")
    geom, mat, load, sol = case["geometry"], case["material"], case["loading"], case["solver"]
    n_steps, n_points = Int(sol["steps"]), Int(geom["n_points"])
    0 < n_steps <= Int(sol["max_steps"]) || error("step budget exceeded")
    0 < n_points <= Int(sol["max_points"]) || error("point budget exceeded")
    mat["rho_kg_m3"] == 2300.0 && mat["E_pa"] == 70.0e9 && mat["nu"] == 0.25 && mat["Gc_j_m2"] == 7.6 ||
        error("material values differ from the source-bound illustrative prototype case")
    mat["calibration_status"] == "UNVALIDATED_PROTOTYPE" || error("calibration status must remain unvalidated")
    delta, h = Float64(geom["delta_m"]), Float64(geom["spacing_m"])
    0 < delta && h > 0 || error("invalid point spacing or horizon")
    position = read_points(joinpath(dirname(case_path), case["points_csv"]))
    size(position, 2) == n_points || error("point count differs from frozen case")
    volume = fill(h^3, n_points)
    material = BBMaterial{NoCorrection}(dmgmodel=CriticalStretch())
    inner = Body(material, position, volume)
    material!(inner; horizon=delta, rho=Float64(mat["rho_kg_m3"]), E=Float64(mat["E_pa"]),
               nu=Float64(mat["nu"]), Gc=Float64(mat["Gc_j_m2"]))
    center = (Float64(geom["center_x_m"]), Float64(geom["center_y_m"]), Float64(geom["center_z_m"]))
    body = GeometryFilteredBody(inner, center, Float64(geom["inner_radius_m"]), Float64(geom["outer_radius_m"]))
    com, translation_velocity, omega, initial_velocity = setup_velocity!(body, position, volume,
        Float64(mat["rho_kg_m3"]), Float64(load["alpha_s_1"]))
    solver = VelocityVerlet(steps=n_steps, safety_factor=Float64(sol["safety_factor"]))
    job = Job(body, solver)
    dh = Peridynamics.threads_data_handler(body, solver, 1)
    length(dh.chunks) == 1 || error("native run must have one chunk")
    chunk = dh.chunks[1]
    system, storage = chunk.system, chunk.storage
    point_ids = system.chunk_handler.point_ids
    sort(point_ids) == collect(1:n_points) || error("native point IDs are not a complete one-chunk permutation")
    any(==(0), system.n_neighbors) && error("cavity filter left a point with zero initial native neighbors")
    length(system.bonds) <= 2 * Int(sol["max_pairs"]) || error("directed native bond budget exceeded")
    native_stable_dt = Peridynamics.calc_stable_timestep(dh, 1.0)
    isfinite(native_stable_dt) && native_stable_dt > 0 || error("native stable timestep is invalid")
    Peridynamics.init_time_solver!(solver, dh)
    Peridynamics.initialize!(dh, solver)
    dt = solver.Δt
    dt <= native_stable_dt * Float64(sol["safety_factor"]) * (1 + 1e-12) || error("accepted dt exceeds native bound")
    # Initialize t=0 acceleration from the native force density without changing the IC velocity.
    for local_i in eachindex(point_ids), d in 1:3
        storage.acceleration[d, local_i] = (storage.b_int[d, local_i] + storage.b_ext[d, local_i]) /
                                            Float64(mat["rho_kg_m3"])
    end
    epsilon_c = Float64(inner.point_params[1].εc)
    bond_constant = Float64(inner.point_params[1].bc)
    expected_epsilon = sqrt(5 * Float64(mat["Gc_j_m2"]) /
                            (9 * (Float64(mat["E_pa"]) / (3 * (1 - 2 * Float64(mat["nu"])))) * delta))
    isapprox(epsilon_c, expected_epsilon; rtol=1e-12, atol=0.0) || error("native critical stretch differs from source formula")
    expected_c = 18 * (Float64(mat["E_pa"]) / (3 * (1 - 2 * Float64(mat["nu"])))) / (π * delta^4)
    isapprox(bond_constant, expected_c; rtol=1e-12, atol=0.0) || error("native bond constant differs from BBMaterial formula")

    directed = Dict{Tuple{Int,Int},Vector{Int}}()
    for local_i in Peridynamics.each_point_idx(system)
        global_i = point_ids[local_i]
        for bond_id in Peridynamics.each_bond_idx(system, local_i)
            global_j = point_ids[system.bonds[bond_id].neighbor]
            global_i == global_j && error("native self-bond detected")
            edge = minmax(global_i, global_j)
            push!(get!(directed, edge, Int[]), bond_id)
        end
    end
    edges = sort!(collect(keys(directed)))
    length(edges) <= Int(sol["max_pairs"]) || error("canonical native pair budget exceeded")
    for edge in edges
        length(directed[edge]) == 2 || error("canonical pair $edge does not have two reciprocal native bonds")
    end
    edge_to_bond_ids = directed
    edge_pairs = [(edge[1], edge[2]) for edge in edges]
    rest_length = Vector{Float64}(undef, length(edges))
    bond_index_by_edge = Dict{Tuple{Int,Int},Int}()
    for (k, edge) in enumerate(edges)
        bond_index_by_edge[edge] = k
        bondid = directed[edge][1]
        local_i = system.chunk_handler.localizer[edge[1]]
        rest_length[k] = system.bonds[bondid].length
    end
    point_volume = fill(h^3, n_points)
    point_mass = Float64(mat["rho_kg_m3"]) .* point_volume
    epsilon_by_bond = fill(epsilon_c, length(edges))
    write_rows(joinpath(outdir, "bonds.csv"), "bond_i,bond_j,rest_length_m,epsilon_c",
        [(edge[1] - 1, edge[2] - 1, rest_length[k], epsilon_c) for (k, edge) in enumerate(edges)])
    excluded_edges = sort!(collect(body.excluded); by=first)
    write_rows(joinpath(outdir, "excluded_bonds.csv"), "bond_i,bond_j,closest_radius_m",
        [(edge[1] - 1, edge[2] - 1, radius) for (edge, radius) in excluded_edges])
    length(body.excluded) == length(excluded_edges) || error("duplicate geometric exclusion ledger entries")

    n_states = n_steps + 1
    packed_width = cld(length(edges), 8)
    # Raw histories and a worst-case equally sized NPZ must fit the allocation.
    # Check this before allocating or integrating, rather than after a long run.
    checkpoint_count_upper = 6
    raw_history_bytes = n_states * n_points * 60 + n_states * packed_width +
                        checkpoint_count_upper * n_points * 3 * 8 * 2
    2 * raw_history_bytes <= Int(sol["max_output_bytes"]) ||
        error("full trajectory plus packaging exceeds the frozen output budget before integration")
    position_history = Array{Float64}(undef, n_states, n_points, 3)
    velocity_history = Array{Float64}(undef, n_states, n_points, 3)
    damage_history = Array{Float32}(undef, n_states, n_points)
    component_history = Array{Int32}(undef, n_states, n_points)
    active_history = Array{UInt8}(undef, n_states, packed_width)
    sed_history = Array{Float32}(undef, n_states, n_points)
    checkpoint_steps = sort!(unique!(Int[0, 1, n_steps ÷ 4, n_steps ÷ 2, (3 * n_steps) ÷ 4, n_steps]))
    checkpoint_acceleration = Array{Float64}(undef, length(checkpoint_steps), n_points, 3)
    checkpoint_force = Array{Float64}(undef, length(checkpoint_steps), n_points, 3)
    checkpoint_idx = Dict(step => idx for (idx, step) in enumerate(checkpoint_steps))
    event_bond = Int32[]
    event_step = Int32[]
    event_time = Float64[]
    event_stretch = Float64[]
    event_energy = Float64[]
    ledger_rows = Vector{Tuple}()
    previous_active = save_state!(1, position_history, velocity_history, damage_history,
        component_history, active_history, sed_history, chunk, edges, edge_to_bond_ids, dt)
    function record_state!(step)
        state_index = step + 1
        active = state_active(system, storage, point_ids, edge_to_bond_ids, edges)
        if step > 0
            current = storage.position
            for k in eachindex(edges)
                if previous_active[k] && !active[k]
                    i, j = edge_pairs[k]
                    li = system.chunk_handler.localizer[i]
                    lj = system.chunk_handler.localizer[j]
                    stretch = (norm(@view(current[:, lj]) .- @view(current[:, li])) - rest_length[k]) / rest_length[k]
                    stretch > epsilon_c || error("native first-break does not satisfy strict stretch > epsilon_c")
                    U = 0.5 * bond_constant * point_volume[i] * point_volume[j] * rest_length[k] * stretch^2
                    push!(event_bond, Int32(k - 1)); push!(event_step, Int32(step))
                    push!(event_time, step * dt); push!(event_stretch, stretch); push!(event_energy, U)
                end
            end
        end
        previous_active .= active
        if state_index > 1
            save_state!(state_index, position_history, velocity_history, damage_history,
                component_history, active_history, sed_history, chunk, edges, edge_to_bond_ids, dt)
        end
        if haskey(checkpoint_idx, step)
            ci = checkpoint_idx[step]
            for local_i in eachindex(point_ids)
                global_i = point_ids[local_i]
                for d in 1:3
                    checkpoint_acceleration[ci, global_i, d] = storage.acceleration[d, local_i]
                    checkpoint_force[ci, global_i, d] = storage.b_int[d, local_i] * point_volume[global_i]
                end
            end
        end
        t = step * dt
        kinetic, pair_u, native_u, force, momentum, angular = compute_state_ledger(
            chunk, edge_pairs, active, rest_length, epsilon_c, point_mass, point_volume, bond_constant, t)
        broken_u = isempty(event_energy) ? 0.0 : sum(event_energy)
        push!(ledger_rows, (step, t, kinetic, pair_u, broken_u, native_u,
            force[1], force[2], force[3], momentum[1], momentum[2], momentum[3],
            angular[1], angular[2], angular[3]))
        return nothing
    end
    record_state!(0)
    Δt_half = 0.5 * dt
    for step in 1:n_steps
        Peridynamics.verlet_timestep!(dh, job.options, dt, Δt_half, step)
        record_state!(step)
        if time() - start_time > Float64(sol["max_runtime_s"])
            error("native solver runtime budget exceeded")
        end
    end
    write_rows(joinpath(outdir, "events.csv"), "bond_index,step,time_s,stretch,pair_energy_J",
        zip(event_bond, event_step, event_time, event_stretch, event_energy))
    write_rows(joinpath(outdir, "energy.csv"),
        "step,time_s,kinetic_energy_J,pair_potential_energy_J,break_dissipation_cumulative_J,native_strain_energy_J,total_force_x_N,total_force_y_N,total_force_z_N,momentum_x_kg_m_s,momentum_y_kg_m_s,momentum_z_kg_m_s,angular_momentum_x_kg_m2_s,angular_momentum_y_kg_m2_s,angular_momentum_z_kg_m2_s",
        ledger_rows)
    write_rows(joinpath(outdir, "checkpoints.csv"), "step,time_s",
        [(step, step * dt) for step in checkpoint_steps])
    write_little_endian(joinpath(outdir, "position.f64le"), position_history)
    write_little_endian(joinpath(outdir, "velocity.f64le"), velocity_history)
    write_little_endian(joinpath(outdir, "damage.f32le"), damage_history)
    write_little_endian(joinpath(outdir, "component.i32le"), component_history)
    write_little_endian(joinpath(outdir, "native_sed.f32le"), sed_history)
    write_little_endian(joinpath(outdir, "bond_active_packed.u8"), active_history)
    write_little_endian(joinpath(outdir, "checkpoint_acceleration.f64le"), checkpoint_acceleration)
    write_little_endian(joinpath(outdir, "checkpoint_force.f64le"), checkpoint_force)
    status = Dict{String,Any}(
        "status" => "NATIVE_COMPLETE",
        "solver_model" => "Peridynamics.jl 0.5.4 BBMaterial{NoCorrection} + CriticalStretch + VelocityVerlet",
        "material_model" => "BBMaterial{NoCorrection}", "damage_model" => "CriticalStretch",
        "surface_correction" => "NoCorrection", "time_integrator" => "VelocityVerlet",
        "n_points" => n_points, "canonical_pair_count" => length(edges),
        "initial_geometric_exclusion_count" => length(body.excluded),
        "stable_dt_s" => native_stable_dt, "dt_s" => dt, "time_s" => n_steps * dt,
        "steps" => n_steps, "threads" => Threads.nthreads(), "safety_factor" => Float64(sol["safety_factor"]),
        "calibration_status" => mat["calibration_status"],
        "rho_kg_m3" => Float64(mat["rho_kg_m3"]), "E_pa" => Float64(mat["E_pa"]),
        "nu" => Float64(mat["nu"]), "Gc_j_m2" => Float64(mat["Gc_j_m2"]),
        "point_mass_total_kg" => sum(point_mass), "mass_center_m" => collect(com),
        "loading_alpha_s_1" => Float64(load["alpha_s_1"]),
        "initial_translation_velocity_removed_mps" => collect(translation_velocity),
        "initial_angular_velocity_removed_rad_s" => collect(omega),
        "checkpoint_steps" => checkpoint_steps,
        "initial_velocity_max_mps" => maximum(norm(@view(initial_velocity[:, i])) for i in 1:n_points),
        "scope" => "free-release native bowl evolution; no impact, drop, contact, grip, or precrack",
    )
    write_json(joinpath(dirname(outdir), "native_status.json"), status)
    println("NATIVE_COMPLETE steps=$n_steps points=$n_points pairs=$(length(edges)) dt_s=$dt")
    return nothing
end

const start_time = time()
length(ARGS) == 2 || error("usage: run_bowl_dynamics.jl case.toml output-directory")
run_bowl(ARGS[1], joinpath(ARGS[2], "native_raw"))
