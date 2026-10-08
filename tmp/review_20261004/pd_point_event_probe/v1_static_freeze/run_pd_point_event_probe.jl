using Peridynamics
using LinearAlgebra
using TOML

function require(condition, message)
    condition || error(message)
end

mutable struct EventRecorder
    point_rows::Vector{Vector{Any}}
    directed_state_rows::Vector{Vector{Any}}
    canonical_state_rows::Vector{Vector{Any}}
    directed_transition_rows::Vector{Vector{Any}}
    canonical_transition_rows::Vector{Vector{Any}}
    previous_directed::Union{Nothing,Dict{Tuple{Int,Int},Bool}}
    previous_canonical::Union{Nothing,Dict{Tuple{Int,Int},Bool}}
    last_step::Int
    sample_count::Int
    initial_directed_count::Int
    initial_canonical_count::Int
end

EventRecorder() = EventRecorder(Vector{Vector{Any}}(), Vector{Vector{Any}}(),
                                Vector{Vector{Any}}(), Vector{Vector{Any}}(),
                                Vector{Vector{Any}}(), nothing, nothing, -1, 0, 0, 0)

mutable struct EventProbeOptions <: Peridynamics.AbstractJobOptions
    export_allowed::Bool
    recorder::EventRecorder
    solver::Peridynamics.VelocityVerlet
    point_ids_by_native_index::Vector{Int}
    rest_position::Matrix{Float64}
    volumes::Vector{Float64}
    rho_kg_m3::Float64
end

# The pinned headless package's logger only specializes JobOptions. This option subtype
# deliberately disables VTK/logfile export while observing the native result hooks below.
Peridynamics.log_job_options(::EventProbeOptions) = nothing

function canonical_pair(i::Int, j::Int)
    return i < j ? (i, j) : (j, i)
end

function map_ids_to_native_indices(point_ids::Vector{Int}, selected_ids, label::String)
    id_to_index = Dict(point_id => index for (index, point_id) in enumerate(point_ids))
    selected = Int.(selected_ids)
    require(!isempty(selected), "$label point-id set must be non-empty")
    require(length(unique(selected)) == length(selected), "$label point-id set has duplicates")
    require(all(haskey(id_to_index, point_id) for point_id in selected),
            "$label point-id set references an unknown point ID")
    return [id_to_index[point_id] for point_id in selected]
end

function load_point_input(points_path::String, metadata_path::String, max_points::Int)
    metadata = TOML.parsefile(metadata_path)
    require(get(metadata, "schema", "") == "pgfrag.pd.pointset.v1",
            "unsupported point-input schema")
    require(get(metadata, "coordinate_units", "") == "m",
            "point-input coordinates must use SI metres")
    require(get(metadata, "volume_units", "") == "m^3",
            "point-input volumes must use SI cubic metres")
    require(!isempty(strip(string(get(metadata, "coordinate_frame", "")))),
            "point-input coordinate frame must be declared")
    require(!isempty(strip(string(get(metadata, "point_id_semantics", "")))),
            "point-id semantics must be declared")

    source = get(metadata, "source", Dict{String,Any}())
    require(source isa AbstractDict, "point-input source provenance is required")
    for key in ("kind", "description", "reference", "source_sha256")
        require(!isempty(strip(string(get(source, key, "")))),
                "point-input source.$key must be declared")
    end
    source_hash = string(source["source_sha256"])
    require(length(source_hash) == 64 && all(isxdigit, source_hash),
            "point-input source.sha256 must be a 64-character hex digest")
    quadrature = get(metadata, "quadrature", Dict{String,Any}())
    require(quadrature isa AbstractDict, "point-input quadrature declaration is required")
    require(!isempty(strip(string(get(quadrature, "method", "")))),
            "quadrature method must be declared")
    require(!isempty(strip(string(get(quadrature, "volume_semantics", "")))),
            "quadrature volume semantics must be declared")
    require(get(quadrature, "voronoi_claim", true) == false &&
            get(quadrature, "real_shell_claim", true) == false,
            "input may not claim exact Voronoi volumes or a reconstructed real shell")

    expected_header = "point_id,x_m,y_m,z_m,volume_m3"
    lines = readlines(points_path)
    require(!isempty(lines) && strip(first(lines)) == expected_header,
            "point CSV header must be $expected_header")
    records = NamedTuple[]
    for (row_number, line) in enumerate(lines[2:end], 2)
        isempty(strip(line)) && continue
        fields = split(strip(line), ',')
        require(length(fields) == 5, "point CSV row $row_number must have 5 numeric fields")
        point_id = parse(Int, fields[1])
        xyzv = parse.(Float64, fields[2:5])
        require(point_id > 0, "point IDs must be positive integers")
        require(all(isfinite, xyzv), "point CSV row $row_number contains non-finite values")
        require(xyzv[4] > 0, "point CSV row $row_number has a non-positive volume")
        push!(records, (point_id=point_id, x_m=xyzv[1], y_m=xyzv[2], z_m=xyzv[3],
                        volume_m3=xyzv[4]))
    end

    n = length(records)
    require(0 < n <= max_points <= 512,
            "point count must be in 1:max_points and max_points may not exceed 512")
    require(n == Int(get(metadata, "point_count", -1)),
            "point count disagrees with point-input metadata")
    point_ids = [row.point_id for row in records]
    require(length(unique(point_ids)) == n, "point IDs must be unique")
    position = Matrix{Float64}(undef, 3, n)
    volume = Vector{Float64}(undef, n)
    for (index, row) in enumerate(records)
        position[:, index] .= (row.x_m, row.y_m, row.z_m)
        volume[index] = row.volume_m3
    end

    volume_sum = sum(volume)
    declared_volume_sum = Float64(get(metadata, "declared_volume_sum_m3", NaN))
    require(isfinite(declared_volume_sum) &&
            isapprox(volume_sum, declared_volume_sum; rtol=1e-11, atol=1e-18),
            "input volume sum does not match the declared SI volume sum")
    declared_com = Float64.(get(metadata, "declared_com_m", Float64[]))
    require(length(declared_com) == 3 && all(isfinite, declared_com),
            "declared_com_m must contain three finite SI coordinates")
    measured_com = vec(position * volume) ./ volume_sum
    require(all(isapprox.(measured_com, declared_com; rtol=1e-10, atol=1e-12)),
            "volume-weighted input COM does not match declared_com_m")

    point_data_hash = string(get(metadata, "point_data_sha256", ""))
    require(length(point_data_hash) == 64 && all(isxdigit, point_data_hash),
            "point_data_sha256 must be a 64-character hex digest; launcher pins its file bytes")
    return (; metadata, point_ids, position, volume, volume_sum, measured_com, point_data_hash)
end

function collect_native_state!(options::EventProbeOptions,
                               dh::Peridynamics.ThreadsBodyDataHandler,
                               step::Int, time_s::Float64)
    require(length(dh.chunks) == 1, "observer requires exactly one native solver chunk")
    require(step == options.recorder.last_step + 1,
            "native output hook steps are missing, duplicated, or out of order")
    dt_s = options.solver.Δt
    require(isfinite(dt_s) && dt_s > 0, "native solver did not provide a positive accepted dt")
    require(isapprox(time_s, step * dt_s; rtol=1e-12, atol=128 * eps(Float64) * dt_s),
            "native output-hook time does not match step × accepted dt")

    chunk = dh.chunks[1]
    system = chunk.system
    storage = chunk.storage
    native_point_ids = system.chunk_handler.point_ids
    n = length(options.point_ids_by_native_index)
    require(length(native_point_ids) == n &&
            length(unique(native_point_ids)) == n && Set(native_point_ids) == Set(1:n),
            "native chunk point IDs do not map one-to-one to frozen input row order")
    require(system.chunk_handler.n_loc_points == n,
            "single chunk unexpectedly contains halo or missing material points")

    for local_i in 1:n
        native_i = native_point_ids[local_i]
        source_i = options.point_ids_by_native_index[native_i]
        mass_kg = options.rho_kg_m3 * options.volumes[native_i]
        push!(options.recorder.point_rows, Any[
            step, time_s, dt_s, native_i, source_i,
            options.rest_position[1, native_i], options.rest_position[2, native_i],
            options.rest_position[3, native_i],
            storage.position[1, local_i], storage.position[2, local_i],
            storage.position[3, local_i],
            storage.displacement[1, local_i], storage.displacement[2, local_i],
            storage.displacement[3, local_i],
            storage.velocity[1, local_i], storage.velocity[2, local_i],
            storage.velocity[3, local_i],
            storage.acceleration[1, local_i], storage.acceleration[2, local_i],
            storage.acceleration[3, local_i],
            storage.b_int[1, local_i], storage.b_int[2, local_i], storage.b_int[3, local_i],
            storage.b_ext[1, local_i], storage.b_ext[2, local_i], storage.b_ext[3, local_i],
            storage.damage[local_i], storage.n_active_bonds[local_i],
            options.volumes[native_i], mass_kg,
        ])
    end

    directed = Dict{Tuple{Int,Int},Bool}()
    directions_by_pair = Dict{Tuple{Int,Int},Vector{Tuple{Int,Int}}}()
    for local_i in 1:system.chunk_handler.n_loc_points
        native_i = native_point_ids[local_i]
        source_i = options.point_ids_by_native_index[native_i]
        for bond_id in Peridynamics.each_bond_idx(system, local_i)
            bond = system.bonds[bond_id]
            require(1 <= bond.neighbor <= length(native_point_ids),
                    "native bond neighbor index is outside the one-chunk point map")
            native_j = native_point_ids[bond.neighbor]
            source_j = options.point_ids_by_native_index[native_j]
            require(source_i != source_j, "native solver produced a self-bond")
            direction = (source_i, source_j)
            require(!haskey(directed, direction), "duplicate directed native bond entry")
            directed[direction] = storage.bond_active[bond_id]
            pair = canonical_pair(source_i, source_j)
            push!(get!(directions_by_pair, pair, Tuple{Int,Int}[]), direction)
        end
    end

    canonical = Dict{Tuple{Int,Int},Bool}()
    for pair in sort!(collect(keys(directions_by_pair)))
        directions = sort!(directions_by_pair[pair])
        require(length(directions) == 2 && directions[1] == reverse(directions[2]),
                "canonical bond $pair does not have exactly two reciprocal directed entries")
        active_a = directed[directions[1]]
        active_b = directed[directions[2]]
        require(active_a == active_b, "reciprocal native bond states disagree for $pair")
        canonical[pair] = active_a
    end
    require(length(directed) == 2 * length(canonical),
            "directed/canonical bond counts are inconsistent")

    if step == 0
        require(all(values(directed)), "initial state unexpectedly contains an inactive bond")
        options.recorder.initial_directed_count = length(directed)
        options.recorder.initial_canonical_count = length(canonical)
    else
        require(options.recorder.previous_directed !== nothing &&
                options.recorder.previous_canonical !== nothing,
                "native reference output hook did not capture step 0 before integration")
        require(Set(keys(directed)) == Set(keys(options.recorder.previous_directed)),
                "directed native bond graph changed during integration")
        require(Set(keys(canonical)) == Set(keys(options.recorder.previous_canonical)),
                "canonical native bond graph changed during integration")
        for direction in sort!(collect(keys(directed)))
            old_active = options.recorder.previous_directed[direction]
            new_active = directed[direction]
            if old_active != new_active
                pair = canonical_pair(direction...)
                push!(options.recorder.directed_transition_rows,
                      Any[step, time_s, dt_s, direction[1], direction[2], pair[1], pair[2],
                          old_active, new_active])
            end
        end
        for pair in sort!(collect(keys(canonical)))
            old_active = options.recorder.previous_canonical[pair]
            new_active = canonical[pair]
            if old_active != new_active
                push!(options.recorder.canonical_transition_rows,
                      Any[step, time_s, dt_s, pair[1], pair[2], old_active, new_active])
            end
        end
    end

    for direction in sort!(collect(keys(directed)))
        pair = canonical_pair(direction...)
        push!(options.recorder.directed_state_rows,
              Any[step, time_s, dt_s, direction[1], direction[2], pair[1], pair[2],
                  directed[direction], canonical[pair]])
    end
    for pair in sort!(collect(keys(canonical)))
        push!(options.recorder.canonical_state_rows,
              Any[step, time_s, dt_s, pair[1], pair[2], canonical[pair]])
    end

    options.recorder.previous_directed = directed
    options.recorder.previous_canonical = canonical
    options.recorder.last_step = step
    options.recorder.sample_count += 1
    return nothing
end

# Native solve! dispatches this once after initialize!(dh, solver), before step 1.
function Peridynamics.export_reference_results(dh::Peridynamics.ThreadsBodyDataHandler,
                                                options::EventProbeOptions)
    collect_native_state!(options, dh, 0, 0.0)
    return nothing
end

# Native velocity_verlet.jl dispatches this after each completed verlet_timestep! update.
function Peridynamics.export_results(dh::Peridynamics.ThreadsBodyDataHandler,
                                     options::EventProbeOptions, chunk_id::Int,
                                     timestep::Int, time::Float64)
    require(chunk_id == 1, "single-chunk observer received an unexpected chunk id")
    collect_native_state!(options, dh, timestep, time)
    return nothing
end

function write_csv(path::String, header::String, rows::Vector{Vector{Any}})
    open(path, "w") do io
        println(io, header)
        for row in rows
            println(io, join(string.(row), ','))
        end
    end
    return nothing
end

function run_probe(args)
    require(length(args) == 5,
            "expected arguments: snapshot_dir point_csv point_metadata.toml common.toml case.toml")
    snapshot_dir, points_path, metadata_path, common_path, case_path = args
    require(!ispath(snapshot_dir), "snapshot output directory already exists; refusing overwrite")
    mkpath(snapshot_dir)

    common = TOML.parsefile(common_path)
    case = TOML.parsefile(case_path)
    limits = common["limits"]
    material_config = common["material"]
    geometry = common["geometry"]
    solver_config = common["solver"]
    max_points = Int(limits["max_points"])
    require(max_points <= 512, "configured max_points exceeds 512")
    require(Int(limits["max_chunks"]) == 1,
            "configured max_chunks must remain one for event observer completeness")
    require(string(solver_config["dt_mode"]) == "AUTO_NATIVE_STABLE_BOUND",
            "probe requires native automatic stable timestep selection")
    require(Threads.nthreads() == 1, "this observer requires JULIA_NUM_THREADS=1")

    input = load_point_input(points_path, metadata_path, max_points)
    require(get(material_config, "model", "") == "OSB",
            "this experimental observer is pinned to native OSBMaterial")
    rho = Float64(material_config["rho_kg_m3"])
    young = Float64(material_config["E_pa"])
    nu = Float64(material_config["nu"])
    gc = Float64(material_config["Gc_j_m2"])
    require(rho > 0 && young > 0 && 0 < nu < 0.5 && gc > 0,
            "invalid OSB material values")
    horizon_m = Float64(geometry["horizon_m"])
    require(horizon_m > 0, "horizon_m must be positive")
    steps = Int(case["solver"]["steps"])
    safety = Float64(solver_config["safety_factor"])
    require(0 < steps <= Int(limits["max_steps"]), "solver steps exceed configured limits")
    require(0 < safety < 1, "solver safety factor must be in (0,1)")

    body = Body(OSBMaterial(), input.position, input.volume)
    load = case["load"]
    load_kind = string(load["kind"])
    load_speed = Float64(get(load, "speed_m_s", 0.0))
    load_axis = Int(get(load, "axis", 3))
    require(load_kind in ("zero_load", "symmetric_axial_velocity"),
            "unsupported probe load kind")
    require(abs(load_speed) < Inf && load_speed >= 0, "load speed must be finite and non-negative")
    if load_kind == "zero_load"
        require(iszero(load_speed), "zero_load case must specify speed_m_s = 0")
    else
        require(1 <= load_axis <= 3, "load axis must be 1, 2, or 3")
        load_sets = input.metadata["load_sets"]
        require(Int(load_sets["axis"]) == load_axis,
                "case load axis disagrees with frozen input load-set axis")
        bottom_indices = map_ids_to_native_indices(input.point_ids,
                                                   load_sets["bottom_point_ids"], "bottom")
        top_indices = map_ids_to_native_indices(input.point_ids,
                                                load_sets["top_point_ids"], "top")
        require(isempty(intersect(bottom_indices, top_indices)),
                "top and bottom point sets must be disjoint")
        point_set!(body, :probe_bottom, bottom_indices)
        point_set!(body, :probe_top, top_indices)
        velocity_bc!(t -> -load_speed, body, :probe_bottom, load_axis)
        velocity_bc!(t -> load_speed, body, :probe_top, load_axis)
    end
    material!(body; horizon=horizon_m, rho=rho, E=young, nu=nu, Gc=gc)

    solver = VelocityVerlet(steps=steps, safety_factor=safety)
    stability_dh = Peridynamics.threads_data_handler(body, solver, Threads.nthreads())
    require(length(stability_dh.chunks) == 1,
            "native stability handler must contain exactly one solver chunk")
    native_stability_bound = Peridynamics.calc_stable_timestep(stability_dh, Float64(safety))
    require(isfinite(native_stability_bound) && native_stability_bound > 0,
            "STOP_UNVERIFIED_DT: native stable timestep is not finite and positive")
    options = EventProbeOptions(false, EventRecorder(), solver, input.point_ids,
                                input.position, input.volume, rho)
    job = Job(body, solver, options)
    dh = submit(job; quiet=true)
    require(length(dh.chunks) == 1, "native submit returned more than one solver chunk")
    dt_tolerance = 128 * eps(Float64) * native_stability_bound
    require(isapprox(solver.Δt, native_stability_bound; rtol=1e-12, atol=dt_tolerance),
            "STOP_UNVERIFIED_DT: accepted dt does not match native stable timestep")
    require(options.recorder.sample_count == solver.n_steps + 1,
            "native reference/step output hooks did not produce steps 0 through n_steps")
    require(options.recorder.last_step == solver.n_steps,
            "native event observer did not finish on the requested final step")

    if load_kind == "zero_load"
        require(isempty(options.recorder.directed_transition_rows) &&
                isempty(options.recorder.canonical_transition_rows),
                "zero-load control unexpectedly changed native bond activity")
    end

    write_csv(joinpath(snapshot_dir, "point_state_history.csv"),
              "step,time_s,dt_s,native_body_index,point_id,x_ref_m,y_ref_m,z_ref_m,x_current_m,y_current_m,z_current_m,ux_m,uy_m,uz_m,vx_mps,vy_mps,vz_mps,ax_mps2,ay_mps2,az_mps2,b_int_x_N_m3,b_int_y_N_m3,b_int_z_N_m3,b_ext_x_N_m3,b_ext_y_N_m3,b_ext_z_N_m3,damage,n_active_bonds,volume_m3,mass_kg",
              options.recorder.point_rows)
    write_csv(joinpath(snapshot_dir, "directed_bond_state_history.csv"),
              "step,time_s,dt_s,directed_source_point_id,directed_neighbor_point_id,canonical_point_i_id,canonical_point_j_id,directed_active,canonical_active",
              options.recorder.directed_state_rows)
    write_csv(joinpath(snapshot_dir, "canonical_bond_state_history.csv"),
              "step,time_s,dt_s,point_i_id,point_j_id,active",
              options.recorder.canonical_state_rows)
    write_csv(joinpath(snapshot_dir, "directed_bond_transitions.csv"),
              "step,time_s,dt_s,directed_source_point_id,directed_neighbor_point_id,canonical_point_i_id,canonical_point_j_id,previous_active,current_active",
              options.recorder.directed_transition_rows)
    write_csv(joinpath(snapshot_dir, "canonical_bond_transitions.csv"),
              "step,time_s,dt_s,point_i_id,point_j_id,previous_active,current_active",
              options.recorder.canonical_transition_rows)

    summary = Dict{String,Any}(
        "status" => "COMPLETED",
        "experiment_kind" => "isolated_native_pd_point_event_probe",
        "input_id" => string(input.metadata["input_id"]),
        "point_input_sha256" => input.point_data_hash,
        "source_point_count" => length(input.point_ids),
        "native_body_index_semantics" => "row order in the frozen point CSV; point_id is the preserved source identity",
        "point_volume_sum_m3" => input.volume_sum,
        "volume_weighted_com_m" => input.measured_com,
        "point_mass_sum_kg" => rho * input.volume_sum,
        "native_material" => "OSBMaterial / CriticalStretch",
        "rho_kg_m3" => rho,
        "E_pa" => young,
        "nu" => nu,
        "Gc_j_m2" => gc,
        "calibration_status" => "UNVALIDATED_PROTOTYPE",
        "horizon_m" => horizon_m,
        "load_kind" => load_kind,
        "load_speed_m_s" => load_speed,
        "load_axis" => load_axis,
        "n_steps" => solver.n_steps,
        "dt_s" => solver.Δt,
        "time_s" => solver.n_steps * solver.Δt,
        "native_stability_bound_dt_s" => native_stability_bound,
        "native_stability_check" => "PASS_NATIVE_SINGLE_CHUNK_AUTO_DT",
        "native_thread_count" => Threads.nthreads(),
        "native_chunk_count" => length(dh.chunks),
        "sample_count_including_step_zero" => options.recorder.sample_count,
        "directed_bond_state_rows" => length(options.recorder.directed_state_rows),
        "canonical_bond_state_rows" => length(options.recorder.canonical_state_rows),
        "directed_transition_count" => length(options.recorder.directed_transition_rows),
        "canonical_transition_count" => length(options.recorder.canonical_transition_rows),
        "initial_directed_bond_count" => options.recorder.initial_directed_count,
        "initial_canonical_undirected_bond_count" => options.recorder.initial_canonical_count,
        "reciprocal_state_consistency" => "PASS_AT_EVERY_CAPTURED_STEP",
        "input_volume_and_com_validation" => "PASS",
        "contact" => "UNKNOWN: contact is outside this single-body observer probe",
        "energy_ledger" => "UNKNOWN: no complete work/energy ledger is emitted; missing quantities are not replaced with zero",
        "physical_fragment_count" => "UNKNOWN: no physical fragment or surface reconstruction is claimed",
        "R2_or_training_eligibility" => "UNKNOWN: native graph events do not establish verified physical fragments or R2 lineage",
    )
    open(joinpath(snapshot_dir, "summary.toml"), "w") do io
        TOML.print(io, summary)
    end
    println("PD_POINT_EVENT_PROBE_COMPLETED points=$(length(input.point_ids)) steps=$(solver.n_steps) directed_transitions=$(length(options.recorder.directed_transition_rows)) canonical_transitions=$(length(options.recorder.canonical_transition_rows)) dt_s=$(solver.Δt)")
    return nothing
end

try
    run_probe(ARGS)
catch err
    showerror(stderr, err, catch_backtrace())
    println(stderr)
    exit(1)
end
