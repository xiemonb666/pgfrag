using Peridynamics
using LinearAlgebra
using SHA

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const DEFAULT_OUT = joinpath(ROOT, "relicweave-synlab", "tmp", "review_20261003",
                             "implementation", "pg", "contact_spike")
const DX = 2e-3
const DELTA = 3e-3
const RHO = 1000.0
const YOUNG = 1e6
const NU = 0.25
const GC = 1.0
const CONTACT_RADIUS = 2e-3
const PENALTY_FACTOR = 1e6
const SPEED = 0.2
const MAX_POINTS = 128
const WALL_LIMIT_S = 120.0

function require_ok(ok, message)
    ok || error(message)
    return nothing
end

function csv_cell(x)
    s = string(x)
    if occursin(',', s) || occursin('"', s) || occursin('\n', s) || occursin('\r', s)
        return "\"" * replace(s, "\"" => "\"\"") * "\""
    end
    return s
end

function write_csv(path, header, rows)
    open(path, "w") do io
        println(io, join(header, ','))
        for row in rows
            require_ok(length(row) == length(header),
                       "STOP_CSV_SCHEMA: row width $(length(row)) differs from header width $(length(header)) for $path")
            println(io, join(csv_cell.(row), ','))
        end
    end
    return nothing
end

function append_csv(path, rows)
    expected_width = length(split(first(readlines(path)), ','))
    open(path, "a") do io
        for row in rows
            require_ok(length(row) == expected_width,
                       "STOP_CSV_SCHEMA: row width $(length(row)) differs from existing header width $expected_width for $path")
            println(io, join(csv_cell.(row), ','))
        end
    end
    return nothing
end

function json_string(s::AbstractString)
    io = IOBuffer()
    print(io, '"')
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
        elseif Int(c) < 0x20
            print(io, "\\u", string(Int(c); base=16, pad=4))
        else
            print(io, c)
        end
    end
    print(io, '"')
    return String(take!(io))
end

function json_value(x)
    if x === nothing
        return "null"
    elseif x isa Bool
        return x ? "true" : "false"
    elseif x isa AbstractString || x isa Symbol
        return json_string(string(x))
    elseif x isa Integer
        return string(x)
    elseif x isa AbstractFloat
        return isfinite(x) ? repr(x) : "null"
    elseif x isa AbstractDict
        keys_sorted = sort!(collect(keys(x)); by=string)
        return "{" * join((json_string(string(k)) * ":" * json_value(x[k]) for k in keys_sorted), ",") * "}"
    elseif x isa Tuple || x isa AbstractVector
        return "[" * join(json_value.(collect(x)), ",") * "]"
    end
    return json_string(string(x))
end

function write_json(path, value)
    open(path, "w") do io
        println(io, json_value(value))
    end
    return nothing
end

function sha256_file(path)
    return bytes2hex(open(sha256, path))
end

function all_point_states(dh)
    rows = NamedTuple[]
    for body_id in Peridynamics.each_body_idx(dh)
        body_name = dh.body_names[body_id]
        body_dh = Peridynamics.get_body_dh(dh, body_id)
        for chunk in body_dh.chunks
            system, storage = chunk.system, chunk.storage
            body_local_ids = system.chunk_handler.point_ids
            n_local = system.chunk_handler.n_loc_points
            for local_i in 1:n_local
                body_local_id = body_local_ids[local_i]
                push!(rows, (body=string(body_name), point_id=body_local_id,
                    export_point_id=(body_name === :a ? body_local_id : body_local_id + 8),
                    position=copy(storage.position[:, local_i]),
                    velocity=copy(storage.velocity[:, local_i]),
                    total_b_int=copy(storage.b_int[:, local_i]),
                    volume=system.volume[local_i],
                    mass=RHO * system.volume[local_i]))
            end
        end
    end
    sort!(rows; by=x -> (x.body, x.point_id))
    return rows
end

function point_positions(dh, body_id::Symbol)
    body_dh = Peridynamics.get_body_dh(dh, body_id)
    positions = Vector{Vector{Float64}}()
    for chunk in body_dh.chunks
        n_local = chunk.system.chunk_handler.n_loc_points
        for i in 1:n_local
            push!(positions, copy(chunk.storage.position[:, i]))
        end
    end
    return positions
end

function verify_native_point_mapping(dh)
    require_ok(dh.body_names == [:a, :b],
               "native body order must match explicit names :a and :b")
    require_ok(length(dh.body_dhs) == 2 && all(length(body_dh.chunks) == 1 for body_dh in dh.body_dhs),
               "expected two bodies and one native chunk per body")
    native_local_ids = [[chunk.system.chunk_handler.point_ids[i]
                         for chunk in body_dh.chunks for i in 1:chunk.system.chunk_handler.n_loc_points]
                        for body_dh in dh.body_dhs]
    require_ok(native_local_ids == [collect(1:8), collect(1:8)],
               "native body-local point IDs must each be complete 1:8")
    export_ids = vcat(native_local_ids[1], 8 .+ native_local_ids[2])
    require_ok(length(export_ids) == 16 && length(unique(export_ids)) == 16 &&
               Set(export_ids) == Set(1:16),
               "combined export point IDs must be complete and unique for all 16 points")
    return export_ids
end

function capture_total_b_int(dh)
    total = Dict{Tuple{Symbol,Int},Vector{Float64}}()
    for body_id in (:a, :b)
        body_dh = Peridynamics.get_body_dh(dh, body_id)
        for chunk in body_dh.chunks
            for local_i in 1:chunk.system.chunk_handler.n_loc_points
                local_id = chunk.system.chunk_handler.point_ids[local_i]
                total[(body_id, local_id)] = copy(chunk.storage.b_int[:, local_i])
            end
        end
    end
    require_ok(length(total) == 16 && length(unique(keys(total))) == 16,
               "native total b_int samples must cover 16 unique (body, local_id) pairs")
    return total
end

function minimum_cross_distance(positions_a, positions_b)
    min_distance = Inf
    for p in positions_a, q in positions_b
        d = norm(p - q)
        require_ok(isfinite(d) && d > 0, "STOP_INVALID_CONTACT_DISTANCE: cross-body L must be finite and > 0")
        min_distance = min(min_distance, d)
    end
    return min_distance
end

function predicted_minimum_cross_distance(dh, dt)
    predicted = Dict{Symbol,Vector{Vector{Float64}}}()
    for body_id in (:a, :b)
        body_dh = Peridynamics.get_body_dh(dh, body_id)
        positions = Vector{Vector{Float64}}()
        for chunk in body_dh.chunks
            storage = chunk.storage
            n_local = chunk.system.chunk_handler.n_loc_points
            for i in 1:n_local
                p = storage.position[:, i] .+ (storage.velocity[:, i] .+ 0.5dt .* storage.acceleration[:, i]) .* dt
                require_ok(all(isfinite, p), "STOP_INVALID_CONTACT_DISTANCE: predicted position is not finite")
                push!(positions, p)
            end
        end
        predicted[body_id] = positions
    end
    return minimum_cross_distance(predicted[:a], predicted[:b])
end

function preflight_positive_bond_lengths(dh, dt)
    for body_dh in dh.body_dhs
        for chunk in body_dh.chunks
            storage, system = chunk.storage, chunk.system
            predicted = similar(storage.position)
            local_points = 1:system.chunk_handler.n_loc_points
            for i in local_points
                @views predicted[:, i] .= storage.position[:, i] .+
                    (storage.velocity[:, i] .+ 0.5 * dt .* storage.acceleration[:, i]) .* dt
            end
            for i in local_points
                for bond_id in Peridynamics.each_bond_idx(system, i)
                    bond = system.bonds[bond_id]
                    L = norm(@view(predicted[:, i]) .- @view(predicted[:, bond.neighbor]))
                    require_ok(isfinite(L) && L > 0,
                               "STOP_INVALID_BOND_DISTANCE: native bond L must be finite and > 0 before force evaluation")
                end
            end
        end
    end
    return nothing
end

function snapshot_storage_arrays(storage)
    snapshots = Dict{Symbol,Any}()
    for field in fieldnames(typeof(storage))
        value = getfield(storage, field)
        value isa AbstractArray && (snapshots[field] = copy(value))
    end
    return snapshots
end

function restore_storage!(storage, snapshots)
    for (field, saved) in snapshots
        copyto!(getfield(storage, field), saved)
    end
    return nothing
end

function unchanged_storage_fields(storage, snapshots; ignore=())
    changed = Symbol[]
    for (field, saved) in snapshots
        field in ignore && continue
        isequal(getfield(storage, field), saved) || push!(changed, field)
    end
    return changed
end

function split_native_bond_and_contact_density(dh, total_b_int, t, dt)
    chunks = [(chunk, snapshot_storage_arrays(chunk.storage))
              for body_dh in dh.body_dhs for chunk in body_dh.chunks]
    bond_b_int = Dict{Tuple{Symbol,Int},Vector{Float64}}()
    try
        for body_id in (:a, :b)
            body_dh = Peridynamics.get_body_dh(dh, body_id)
            Peridynamics.calc_force_density!(body_dh, t, dt)
            for chunk in body_dh.chunks
                for local_i in 1:chunk.system.chunk_handler.n_loc_points
                    local_id = chunk.system.chunk_handler.point_ids[local_i]
                    bond_b_int[(body_id, local_id)] = copy(chunk.storage.b_int[:, local_i])
                end
            end
        end
    finally
        changed_by_diagnostic = Symbol[]
        for (chunk, snapshot) in chunks
            append!(changed_by_diagnostic,
                    [Symbol(string(chunk.body_name, ".", field))
                     for field in unchanged_storage_fields(chunk.storage, snapshot; ignore=(:b_int,))])
        end
        for (chunk, snapshot) in chunks
            restore_storage!(chunk.storage, snapshot)
        end
        restored_mismatch = Symbol[]
        for (chunk, snapshot) in chunks
            append!(restored_mismatch,
                    [Symbol(string(chunk.body_name, ".", field))
                     for field in unchanged_storage_fields(chunk.storage, snapshot)])
        end
        require_ok(isempty(changed_by_diagnostic),
                   "STOP_DIAGNOSTIC_STATE_MUTATION: diagnostic bond pass changed native fields $(changed_by_diagnostic)")
        require_ok(isempty(restored_mismatch),
                   "STOP_DIAGNOSTIC_STATE_RESTORE: native arrays were not restored $(restored_mismatch)")
    end
    contact_b_int = Dict{Tuple{Symbol,Int},Vector{Float64}}()
    for (key, total) in total_b_int
        haskey(bond_b_int, key) || error("missing native BB bond density for point $key")
        contact_b_int[key] = total .- bond_b_int[key]
    end
    return bond_b_int, contact_b_int, true
end

function point_csv_header()
    return ["case_id", "step", "time_s", "export_point_id", "body", "body_local_point_id",
        "x_m", "y_m", "z_m", "vx_mps", "vy_mps", "vz_mps", "mass_kg",
        "native_b_int_total_x_N_m3", "native_b_int_total_y_N_m3", "native_b_int_total_z_N_m3",
        "native_bb_bond_b_int_x_N_m3", "native_bb_bond_b_int_y_N_m3", "native_bb_bond_b_int_z_N_m3",
        "native_contact_b_int_increment_x_N_m3", "native_contact_b_int_increment_y_N_m3", "native_contact_b_int_increment_z_N_m3",
        "native_contact_force_x_N", "native_contact_force_y_N", "native_contact_force_z_N",
        "native_contact_impulse_trapezoid_step_x_Ns", "native_contact_impulse_trapezoid_step_y_Ns", "native_contact_impulse_trapezoid_step_z_Ns"]
end

function body_momenta(rows)
    momenta = Dict(:a => zeros(3), :b => zeros(3))
    kinetic = 0.0
    total_momentum = zeros(3)
    for p in rows
        momenta[Symbol(p.body)] .+= p.mass .* p.velocity
        total_momentum .+= p.mass .* p.velocity
        kinetic += 0.5 * p.mass * sum(abs2, p.velocity)
    end
    return momenta, total_momentum, kinetic
end

function neighbor_metrics(dh)
    isempty(dh.srf_contacts) && return (count=0, min_contact_distance=Inf)
    contact = first(dh.srf_contacts)
    a_idx = dh.body_idxs[contact.body_id_a]
    b_idx = dh.body_idxs[contact.body_id_b]
    pos_a = dh.position_caches[a_idx]
    pos_b = dh.position_caches[b_idx]
    count = 0
    min_distance = Inf
    for i in axes(pos_a, 2)
        Peridynamics.foreach_neighbor(pos_a, pos_b, contact.nhs, i) do _, _, _, L
            require_ok(isfinite(L) && L > 0,
                       "STOP_INVALID_CONTACT_DISTANCE: native contact neighbor L must be finite and > 0")
            count += 1
            min_distance = min(min_distance, L)
        end
    end
    return (count=count, min_contact_distance=min_distance)
end

function check_free_motion_setup(body_a, body_b, initial_a, initial_b)
    require_ok(isempty(body_a.single_dim_bcs) && isempty(body_b.single_dim_bcs),
               "unexpected velocity or other single-dimension boundary condition")
    require_ok(isempty(body_a.posdep_single_dim_bcs) && isempty(body_b.posdep_single_dim_bcs),
               "unexpected position-dependent boundary condition")
    require_ok(isempty(body_a.pos_single_dim_bcs) && isempty(body_b.pos_single_dim_bcs),
               "unexpected displacement boundary condition")
    require_ok(isempty(body_a.data_bcs) && isempty(body_b.data_bcs),
               "unexpected data boundary condition")
    require_ok(size(initial_a, 2) == 8 && size(initial_b, 2) == 8,
               "fixed 4 mm / 2 mm geometry must produce exactly 8 points per body")
    export_ids = vcat(collect(1:size(initial_a, 2)),
                      size(initial_a, 2) .+ collect(1:size(initial_b, 2)))
    require_ok(length(export_ids) == 16 && length(unique(export_ids)) == 16 &&
               Set(export_ids) == Set(1:16),
               "combined export point IDs must be complete and unique for all 16 points")
    require_ok(isapprox(sum(initial_a[3, :]) / size(initial_a, 2), -2.25e-3; atol=1.0e-12, rtol=0),
               "body A initial center offset changed")
    require_ok(isapprox(sum(initial_b[3, :]) / size(initial_b, 2), 2.25e-3; atol=1.0e-12, rtol=0),
               "body B initial center offset changed")
    require_ok(all(isapprox.(body_a.volume, 8e-9; rtol=0, atol=1e-20)) &&
               all(isapprox.(body_b.volume, 8e-9; rtol=0, atol=1e-20)),
               "point volumes must all be dx^3 = 8e-9 m^3")
    return nothing
end

function build_setup(with_contact::Bool)
    pos_a, volume_a = uniform_box(4e-3, 4e-3, 4e-3, DX; center=(0.0, 0.0, -2.25e-3))
    pos_b, volume_b = uniform_box(4e-3, 4e-3, 4e-3, DX; center=(0.0, 0.0, 2.25e-3))
    body_a = Body(BBMaterial(), pos_a, volume_a)
    body_b = Body(BBMaterial(), pos_b, volume_b)
    for body in (body_a, body_b)
        material!(body; horizon=DELTA, rho=RHO, E=YOUNG, nu=NU, Gc=GC)
        no_failure!(body)
    end
    velocity_ic!(body_a, :all_points, :z, SPEED)
    velocity_ic!(body_b, :all_points, :z, -SPEED)
    check_free_motion_setup(body_a, body_b, pos_a, pos_b)
    ms = MultibodySetup(:a => body_a, :b => body_b)
    if with_contact
        contact!(ms, :a, :b; radius=CONTACT_RADIUS, penalty_factor=PENALTY_FACTOR)
    end
    return ms, body_a, body_b, pos_a, pos_b, volume_a, volume_b
end

function native_stability_bound(ms)
    probe = VelocityVerlet(steps=1)
    dh = Peridynamics.threads_data_handler(ms, probe, Threads.nthreads())
    bound = Peridynamics.calc_stable_timestep(dh, 1.0)
    require_ok(isfinite(bound) && bound > 0,
               "STOP_UNVERIFIED_DT: native bond-based stability bound is not positive finite")
    return bound
end

function make_step_rows(dh, case_id, step, time_s, contact_b_int, bond_b_int, total_b_int,
                        previous_contact_force, dt)
    rows = Vector{Tuple}(undef, 16)
    k = 0
    for body_id in (:a, :b)
        body_dh = Peridynamics.get_body_dh(dh, body_id)
        for chunk in body_dh.chunks
            system, storage = chunk.system, chunk.storage
            body_local_ids = system.chunk_handler.point_ids
            for local_i in 1:system.chunk_handler.n_loc_points
                k += 1
                local_id = body_local_ids[local_i]
                export_id = body_id === :a ? local_id : local_id + 8
                key = (body_id, local_id)
                b_total = total_b_int[key]
                b_bond = bond_b_int[key]
                b_contact = contact_b_int[key]
                volume = system.volume[local_i]
                force = b_contact .* volume
                impulse = 0.5 .* (previous_contact_force[key] .+ force) .* dt
                rows[k] = (case_id, step, time_s, export_id, string(body_id), local_id,
                    storage.position[1, local_i], storage.position[2, local_i], storage.position[3, local_i],
                    storage.velocity[1, local_i], storage.velocity[2, local_i], storage.velocity[3, local_i],
                    RHO * volume,
                    b_total[1], b_total[2], b_total[3],
                    b_bond[1], b_bond[2], b_bond[3],
                    b_contact[1], b_contact[2], b_contact[3],
                    force[1], force[2], force[3], impulse[1], impulse[2], impulse[3])
            end
        end
    end
    require_ok(k == 16 && sort([row[4] for row in rows]) == collect(1:16),
               "native point rows must uniquely export combined IDs 1:16")
    expected_local_ids = Set(vcat([(string(:a), i) for i in 1:8],
                                  [(string(:b), i) for i in 1:8]))
    require_ok(Set((row[5], row[6]) for row in rows) == expected_local_ids,
               "native point rows must preserve complete per-body local IDs 1:8")
    return rows
end

function point_contact_forces(dh, contact_by_body)
    result = Dict{Tuple{Symbol,Int},Vector{Float64}}()
    for body_id in (:a, :b)
        body_dh = Peridynamics.get_body_dh(dh, body_id)
        for chunk in body_dh.chunks
            for local_i in 1:chunk.system.chunk_handler.n_loc_points
                local_id = chunk.system.chunk_handler.point_ids[local_i]
                key = (body_id, local_id)
                result[key] = contact_by_body[key] .* chunk.system.volume[local_i]
            end
        end
    end
    return result
end

function sum_point_contact_forces(force_by_point)
    forces = Dict(:a => zeros(3), :b => zeros(3))
    for ((body_id, _), force) in force_by_point
        forces[body_id] .+= force
    end
    return forces
end

function run_native_preflight(outdir, native_bound)
    preflight_start = time()
    preflight_dir = joinpath(outdir, "native_interface_preflight")
    mkpath(preflight_dir)
    point_path = joinpath(preflight_dir, "one_step_point_states.csv")
    metric_path = joinpath(preflight_dir, "one_step_metrics.csv")
    write_csv(point_path, point_csv_header(), Tuple[])
    metric_header = ["step", "time_s", "min_cross_body_point_distance_m", "radius_penetration_m",
        "native_contact_neighbor_pair_count", "min_native_contact_neighbor_distance_m",
        "total_momentum_x_kg_mps", "total_momentum_y_kg_mps", "total_momentum_z_kg_mps",
        "body_a_contact_force_x_N", "body_a_contact_force_y_N", "body_a_contact_force_z_N",
        "body_b_contact_force_x_N", "body_b_contact_force_y_N", "body_b_contact_force_z_N",
        "action_reaction_force_error_norm_N", "action_reaction_relative_error",
        "body_a_delta_momentum_x_kg_mps", "body_a_delta_momentum_y_kg_mps", "body_a_delta_momentum_z_kg_mps",
        "body_b_delta_momentum_x_kg_mps", "body_b_delta_momentum_y_kg_mps", "body_b_delta_momentum_z_kg_mps",
        "body_a_contact_impulse_trapezoid_x_Ns", "body_a_contact_impulse_trapezoid_y_Ns", "body_a_contact_impulse_trapezoid_z_Ns",
        "body_b_contact_impulse_trapezoid_x_Ns", "body_b_contact_impulse_trapezoid_y_Ns", "body_b_contact_impulse_trapezoid_z_Ns",
        "body_a_momentum_impulse_balance_error_norm_kg_mps", "body_b_momentum_impulse_balance_error_norm_kg_mps",
        "diagnostic_state_unchanged"]
    write_csv(metric_path, metric_header, Tuple[])
    progress_path = joinpath(preflight_dir, "progress.json")
    write_json(progress_path, Dict("status" => "RUNNING", "native_steps_completed" => 0,
                                   "fixed_cases_started" => false))

    ms, body_a, body_b, _, _, _, _ = build_setup(true)
    solver = VelocityVerlet(steps=1, stepsize=0.5 * native_bound)
    job = Job(ms, solver)
    dh = Peridynamics.threads_data_handler(ms, solver, Threads.nthreads())
    export_ids = verify_native_point_mapping(dh)
    require_ok(length(export_ids) == 16 && n_points(ms) == 16 && n_points(ms) <= MAX_POINTS,
               "native preflight must contain exactly 16 mapped points within the point cap")
    checked_bound = Peridynamics.calc_stable_timestep(dh, 1.0)
    require_ok(isapprox(checked_bound, native_bound; rtol=1e-12, atol=0.0),
               "native bond bound changed between preflight setups")
    require_ok(solver.Δt == 0.5 * checked_bound,
               "native preflight did not use the fixed 0.5 raw-bound step")
    Peridynamics.init_time_solver!(solver, dh)
    Peridynamics.initialize!(dh, solver)
    initial_min_distance = minimum_cross_distance(point_positions(dh, :a), point_positions(dh, :b))
    require_ok(isapprox(initial_min_distance, 2.5e-3; atol=1e-12, rtol=0.0),
               "native preflight initial minimum point gap must be 2.5 mm")
    require_ok(isempty(body_a.single_dim_bcs) && isempty(body_b.single_dim_bcs) &&
               isempty(body_a.posdep_single_dim_bcs) && isempty(body_b.posdep_single_dim_bcs) &&
               isempty(body_a.pos_single_dim_bcs) && isempty(body_b.pos_single_dim_bcs) &&
               isempty(body_a.data_bcs) && isempty(body_b.data_bcs),
               "native preflight must have no velocity or other boundary conditions")
    require_ok(Threads.nthreads() == 1, "native preflight must use one Julia thread")
    require_ok(all(iszero, reduce(vcat, [chunk.storage.b_ext for body_dh in dh.body_dhs for chunk in body_dh.chunks])),
               "native preflight external force density must be zero")

    initial_state = all_point_states(dh)
    initial_body_momentum, initial_momentum, _ = body_momenta(initial_state)
    initial_total = Dict{Tuple{Symbol,Int},Vector{Float64}}()
    initial_bond = Dict{Tuple{Symbol,Int},Vector{Float64}}()
    zero_contact_density = Dict{Tuple{Symbol,Int},Vector{Float64}}()
    initial_contact_force = Dict{Tuple{Symbol,Int},Vector{Float64}}()
    for p in initial_state
        key = (Symbol(p.body), p.point_id)
        initial_total[key] = copy(p.total_b_int)
        initial_bond[key] = copy(p.total_b_int)
        zero_contact_density[key] = zeros(3)
        initial_contact_force[key] = zeros(3)
    end
    require_ok(length(initial_total) == 16 && length(unique(keys(initial_total))) == 16,
               "native preflight initial (body, local_id) map must cover 16 unique points")
    initial_rows = make_step_rows(dh, "native_interface_preflight", 0, 0.0,
        zero_contact_density, initial_bond, initial_total, initial_contact_force, solver.Δt)
    require_ok(length(initial_rows) == 16, "native preflight must export all 16 initial point states")
    append_csv(point_path, initial_rows)
    initial_neighbors = neighbor_metrics(dh)
    append_csv(metric_path, ((0, 0.0, initial_min_distance,
        max(0.0, CONTACT_RADIUS - initial_min_distance), initial_neighbors.count,
        initial_neighbors.min_contact_distance,
        initial_momentum[1], initial_momentum[2], initial_momentum[3],
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, true),))

    predicted_minimum_cross_distance(dh, solver.Δt)
    preflight_positive_bond_lengths(dh, solver.Δt)
    Peridynamics.verlet_timestep!(dh, job.options, solver.Δt, 0.5 * solver.Δt, 1)
    write_json(progress_path, Dict("status" => "RUNNING", "native_steps_completed" => 1,
                                   "fixed_cases_started" => false))
    total_b_int = capture_total_b_int(dh)
    bond_b_int, contact_b_int, diagnostic_state_unchanged =
        split_native_bond_and_contact_density(dh, total_b_int, solver.Δt, solver.Δt)
    current_state = all_point_states(dh)
    require_ok(all(isfinite(p) for row in current_state for p in row.position) &&
               all(isfinite(v) for row in current_state for v in row.velocity),
               "native preflight produced non-finite positions or velocities")
    point_forces = point_contact_forces(dh, contact_b_int)
    forces = sum_point_contact_forces(point_forces)
    force_error = norm(forces[:a] + forces[:b])
    force_scale = max(norm(forces[:a]), norm(forces[:b]), eps())
    current_body_momentum, current_momentum, _ = body_momenta(current_state)
    delta_momentum = Dict(body_id => current_body_momentum[body_id] .- initial_body_momentum[body_id]
                          for body_id in (:a, :b))
    trapezoid_impulse = Dict(body_id => 0.5 .* forces[body_id] .* solver.Δt for body_id in (:a, :b))
    balance_residual = Dict(body_id => delta_momentum[body_id] .- trapezoid_impulse[body_id]
                            for body_id in (:a, :b))
    final_min_distance = minimum_cross_distance(point_positions(dh, :a), point_positions(dh, :b))
    final_neighbors = neighbor_metrics(dh)
    final_rows = make_step_rows(dh, "native_interface_preflight", 1, solver.Δt,
        contact_b_int, bond_b_int, total_b_int, initial_contact_force, solver.Δt)
    require_ok(length(final_rows) == 16, "native preflight must export all 16 post-step point states")
    append_csv(point_path, final_rows)
    append_csv(metric_path, ((1, solver.Δt, final_min_distance,
        max(0.0, CONTACT_RADIUS - final_min_distance), final_neighbors.count,
        final_neighbors.min_contact_distance,
        current_momentum[1], current_momentum[2], current_momentum[3],
        forces[:a][1], forces[:a][2], forces[:a][3], forces[:b][1], forces[:b][2], forces[:b][3],
        force_error, force_error / force_scale,
        delta_momentum[:a][1], delta_momentum[:a][2], delta_momentum[:a][3],
        delta_momentum[:b][1], delta_momentum[:b][2], delta_momentum[:b][3],
        trapezoid_impulse[:a][1], trapezoid_impulse[:a][2], trapezoid_impulse[:a][3],
        trapezoid_impulse[:b][1], trapezoid_impulse[:b][2], trapezoid_impulse[:b][3],
        norm(balance_residual[:a]), norm(balance_residual[:b]), diagnostic_state_unchanged),))
    require_ok(time() - preflight_start <= WALL_LIMIT_S,
               "STOP_WALL_LIMIT: native preflight exceeded 120 seconds")
    # Diagnostic re-evaluation snapshots/restores the native total b_int. The saved
    # total and point-force maps therefore remain the actual post-contact Verlet values.
    preflight_summary = Dict{String,Any}(
        "status" => "PASS",
        "preflight" => "one native MultibodySetup/contact! VelocityVerlet timestep",
        "native_backend_version" => string(Base.pkgversion(Peridynamics)),
        "native_contact_api" => "Peridynamics.contact! and Peridynamics.verlet_timestep!",
        "points_total" => n_points(ms),
        "points_per_body" => 8,
        "body_local_ids" => Dict("a" => collect(1:8), "b" => collect(1:8)),
        "unique_export_point_ids" => collect(1:16),
        "threads" => Threads.nthreads(),
        "native_bond_stability_bound_dt_s" => checked_bound,
        "native_bond_bound_scope" => "BB bond-only; not a penalty-contact stability proof",
        "dt_s" => solver.Δt,
        "steps" => 1,
        "initial_min_point_distance_m" => initial_min_distance,
        "final_min_point_distance_m" => final_min_distance,
        "initial_native_contact_neighbor_pairs" => initial_neighbors.count,
        "final_native_contact_neighbor_pairs" => final_neighbors.count,
        "native_contact_force_on_a_N" => collect(forces[:a]),
        "native_contact_force_on_b_N" => collect(forces[:b]),
        "action_reaction_force_error_norm_N" => force_error,
        "body_a_delta_momentum_kg_mps" => collect(delta_momentum[:a]),
        "body_b_delta_momentum_kg_mps" => collect(delta_momentum[:b]),
        "body_a_trapezoid_contact_impulse_Ns" => collect(trapezoid_impulse[:a]),
        "body_b_trapezoid_contact_impulse_Ns" => collect(trapezoid_impulse[:b]),
        "body_a_delta_momentum_minus_impulse_kg_mps" => collect(balance_residual[:a]),
        "body_b_delta_momentum_minus_impulse_kg_mps" => collect(balance_residual[:b]),
        "native_total_b_int_captured_before_diagnostic" => true,
        "diagnostic_state_unchanged_and_restored" => diagnostic_state_unchanged,
        "point_csv_rows" => 32,
        "metric_csv_rows" => 2,
        "wall_time_s" => time() - preflight_start)
    write_json(joinpath(preflight_dir, "summary.json"), preflight_summary)
    write_json(progress_path, Dict("status" => "PASS", "native_steps_completed" => 1,
                                   "fixed_cases_started" => false))
    println("NATIVE_INTERFACE_PREFLIGHT_PASSED steps=1 dt_s=$(solver.Δt) points=16 wall_time_s=$(preflight_summary[\"wall_time_s\"])")
    return preflight_summary
end

function run_case(case_id, dt_factor, n_steps, with_contact, outdir, native_bound)
    case_start = time()
    case_dir = joinpath(outdir, case_id)
    mkpath(case_dir)
    point_path = joinpath(case_dir, "point_steps.csv")
    metric_path = joinpath(case_dir, "step_metrics.csv")
    point_header = point_csv_header()
    metric_header = ["case_id", "step", "time_s", "min_cross_body_point_distance_m",
        "radius_penetration_m", "native_contact_neighbor_pair_count", "min_native_contact_neighbor_distance_m",
        "total_momentum_x_kg_mps", "total_momentum_y_kg_mps", "total_momentum_z_kg_mps",
        "total_momentum_drift_norm_kg_mps", "native_contact_force_on_a_x_N", "native_contact_force_on_a_y_N", "native_contact_force_on_a_z_N",
        "native_contact_force_on_b_x_N", "native_contact_force_on_b_y_N", "native_contact_force_on_b_z_N",
        "action_reaction_force_error_norm_N", "action_reaction_relative_error",
        "cumulative_native_contact_impulse_on_a_x_Ns", "cumulative_native_contact_impulse_on_a_y_Ns", "cumulative_native_contact_impulse_on_a_z_Ns",
        "cumulative_native_contact_impulse_on_b_x_Ns", "cumulative_native_contact_impulse_on_b_y_Ns", "cumulative_native_contact_impulse_on_b_z_Ns",
        "kinetic_energy_J",
        "body_a_delta_momentum_x_kg_mps", "body_a_delta_momentum_y_kg_mps", "body_a_delta_momentum_z_kg_mps",
        "body_b_delta_momentum_x_kg_mps", "body_b_delta_momentum_y_kg_mps", "body_b_delta_momentum_z_kg_mps",
        "body_a_contact_impulse_trapezoid_step_x_Ns", "body_a_contact_impulse_trapezoid_step_y_Ns", "body_a_contact_impulse_trapezoid_step_z_Ns",
        "body_b_contact_impulse_trapezoid_step_x_Ns", "body_b_contact_impulse_trapezoid_step_y_Ns", "body_b_contact_impulse_trapezoid_step_z_Ns",
        "body_a_delta_momentum_minus_contact_impulse_x_kg_mps", "body_a_delta_momentum_minus_contact_impulse_y_kg_mps", "body_a_delta_momentum_minus_contact_impulse_z_kg_mps",
        "body_a_momentum_impulse_balance_error_norm_kg_mps",
        "body_b_delta_momentum_minus_contact_impulse_x_kg_mps", "body_b_delta_momentum_minus_contact_impulse_y_kg_mps", "body_b_delta_momentum_minus_contact_impulse_z_kg_mps",
        "body_b_momentum_impulse_balance_error_norm_kg_mps"]
    write_csv(point_path, point_header, Tuple[])
    write_csv(metric_path, metric_header, Tuple[])
    ms, body_a, body_b, _, _, _, _ = build_setup(with_contact)
    total_points = n_points(ms)
    require_ok(total_points <= MAX_POINTS, "point cap exceeded: $total_points > $MAX_POINTS")
    solver = VelocityVerlet(steps=n_steps, stepsize=dt_factor * native_bound)
    job = Job(ms, solver)
    dh = Peridynamics.threads_data_handler(ms, solver, Threads.nthreads())
    verify_native_point_mapping(dh)
    # This is the exact native BB bond bound: calc_stable_timestep multiplies the native
    # pointwise bound by its safety argument, here 1.0. It says nothing about contact.
    checked_bound = Peridynamics.calc_stable_timestep(dh, 1.0)
    require_ok(isapprox(checked_bound, native_bound; rtol=1e-12, atol=0.0),
               "native bond bound changed between case setup and run")
    require_ok(solver.Δt == dt_factor * checked_bound,
               "fixed dt factor was not applied exactly")
    Peridynamics.init_time_solver!(solver, dh)
    Peridynamics.initialize!(dh, solver)
    initial_min_distance = minimum_cross_distance(point_positions(dh, :a), point_positions(dh, :b))
    require_ok(isapprox(initial_min_distance, 2.5e-3; atol=1e-12, rtol=0.0),
               "initial minimum point gap must be 2.5 mm")
    require_ok(all(iszero, reduce(vcat, [chunk.storage.b_ext for body_dh in dh.body_dhs for chunk in body_dh.chunks])),
               "external force density must be zero")
    require_ok(Threads.nthreads() == 1, "native trial must use exactly one Julia thread")

    initial_state = all_point_states(dh)
    previous_body_momentum, initial_momentum, initial_kinetic = body_momenta(initial_state)
    cumulative_impulse = Dict(:a => zeros(3), :b => zeros(3))
    point_rows = Tuple[]
    metric_rows = Tuple[]
    no_contact = Dict{Tuple{Symbol,Int},Vector{Float64}}()
    initial_bond = Dict{Tuple{Symbol,Int},Vector{Float64}}()
    initial_total = Dict{Tuple{Symbol,Int},Vector{Float64}}()
    for p in initial_state
        body_id = Symbol(p.body)
        key = (body_id, p.point_id)
        initial_bond[key] = copy(p.total_b_int)
        initial_total[key] = copy(p.total_b_int)
    end
    for body_id in (:a, :b), i in 1:8
        zero = zeros(3)
        no_contact[(body_id, i)] = zero
    end
    previous_contact_force = deepcopy(no_contact)
    previous_contact_force_by_body = Dict(:a => zeros(3), :b => zeros(3))
    diagnostic_state_unchanged_all_steps = true
    initial_point_rows = make_step_rows(dh, case_id, 0, 0.0, no_contact, initial_bond, initial_total,
                                       previous_contact_force, solver.Δt)
    append!(point_rows, initial_point_rows)
    append_csv(point_path, initial_point_rows)

    initial_contact_metrics = neighbor_metrics(dh)
    initial_metric = (case_id, 0, 0.0, initial_min_distance,
        max(0.0, CONTACT_RADIUS - initial_min_distance), initial_contact_metrics.count,
        initial_contact_metrics.min_contact_distance,
        initial_momentum[1], initial_momentum[2], initial_momentum[3], 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, initial_kinetic,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    push!(metric_rows, initial_metric)
    append_csv(metric_path, (initial_metric,))
    for step in 1:n_steps
        time() - case_start <= WALL_LIMIT_S || error("STOP_WALL_LIMIT: case exceeded 120 seconds")
        # Predict the exact native Verlet position update to reject L<=0 before native
        # contact! evaluates its (r-L)/L term. No coordinates or forces are altered.
        predicted_minimum_cross_distance(dh, solver.Δt)
        preflight_positive_bond_lengths(dh, solver.Δt)
        Peridynamics.verlet_timestep!(dh, job.options, solver.Δt, 0.5 * solver.Δt, step)

        total_b_int = capture_total_b_int(dh)
        # Recover native BB-only b_int diagnostically, then restore and verify every
        # native storage array before any integration or reporting proceeds.
        bond_b_int, contact_b_int, diagnostic_state_unchanged =
            split_native_bond_and_contact_density(dh, total_b_int, step * solver.Δt, solver.Δt)
        diagnostic_state_unchanged_all_steps &= diagnostic_state_unchanged
        for body_id in (:a, :b)
            body_dh = Peridynamics.get_body_dh(dh, body_id)
            for chunk in body_dh.chunks
                require_ok(all(isfinite, chunk.storage.position) && all(isfinite, chunk.storage.velocity),
                           "STOP_NONFINITE_STATE: native solver produced a nonfinite position or velocity")
                require_ok(all(isfinite, chunk.storage.b_ext) && all(iszero, chunk.storage.b_ext),
                           "STOP_NONZERO_EXTERNAL_FORCE: external force density changed during the free-motion trial")
            end
        end

        point_forces = point_contact_forces(dh, contact_b_int)
        forces = sum_point_contact_forces(point_forces)
        force_error = norm(forces[:a] + forces[:b])
        force_scale = max(norm(forces[:a]), norm(forces[:b]), eps())
        current_rows = all_point_states(dh)
        current_body_momentum, momentum, kinetic = body_momenta(current_rows)
        step_contact_impulse = Dict{Symbol,Vector{Float64}}()
        for body_id in (:a, :b)
            step_contact_impulse[body_id] = 0.5 .* (previous_contact_force_by_body[body_id] .+ forces[body_id]) .* solver.Δt
            cumulative_impulse[body_id] .+= step_contact_impulse[body_id]
        end
        delta_momentum = Dict(body_id => current_body_momentum[body_id] .- previous_body_momentum[body_id]
                              for body_id in (:a, :b))
        balance_residual = Dict(body_id => delta_momentum[body_id] .- step_contact_impulse[body_id]
                                for body_id in (:a, :b))
        min_distance = minimum_cross_distance(point_positions(dh, :a), point_positions(dh, :b))
        contact_neighbors = neighbor_metrics(dh)
        radius_penetration = max(0.0, CONTACT_RADIUS - min_distance)
        step_point_rows = make_step_rows(dh, case_id, step, step * solver.Δt, contact_b_int,
            bond_b_int, total_b_int, previous_contact_force, solver.Δt)
        append!(point_rows, step_point_rows)
        append_csv(point_path, step_point_rows)
        metric_row = (case_id, step, step * solver.Δt, min_distance, radius_penetration,
            contact_neighbors.count, contact_neighbors.min_contact_distance,
            momentum[1], momentum[2], momentum[3], norm(momentum - initial_momentum),
            forces[:a][1], forces[:a][2], forces[:a][3], forces[:b][1], forces[:b][2], forces[:b][3],
            force_error, force_error / force_scale,
            cumulative_impulse[:a][1], cumulative_impulse[:a][2], cumulative_impulse[:a][3],
            cumulative_impulse[:b][1], cumulative_impulse[:b][2], cumulative_impulse[:b][3], kinetic,
            delta_momentum[:a][1], delta_momentum[:a][2], delta_momentum[:a][3],
            delta_momentum[:b][1], delta_momentum[:b][2], delta_momentum[:b][3],
            step_contact_impulse[:a][1], step_contact_impulse[:a][2], step_contact_impulse[:a][3],
            step_contact_impulse[:b][1], step_contact_impulse[:b][2], step_contact_impulse[:b][3],
            balance_residual[:a][1], balance_residual[:a][2], balance_residual[:a][3], norm(balance_residual[:a]),
            balance_residual[:b][1], balance_residual[:b][2], balance_residual[:b][3], norm(balance_residual[:b]))
        push!(metric_rows, metric_row)
        append_csv(metric_path, (metric_row,))
        if step % 10 == 0
            open(point_path, "a") do io; flush(io); end
            open(metric_path, "a") do io; flush(io); end
        end
        previous_body_momentum = current_body_momentum
        previous_contact_force_by_body = forces
        previous_contact_force = point_forces
    end

    last_metric = isempty(metric_rows) ? nothing : last(metric_rows)
    final_state = all_point_states(dh)
    final_momentum = isempty(metric_rows) ? initial_momentum : collect(last_metric[8:10])
    summary = Dict{String,Any}(
        "status" => "COMPLETED",
        "case_id" => case_id,
        "native_backend_version" => string(Base.pkgversion(Peridynamics)),
        "native_solver" => "Peridynamics.VelocityVerlet / Peridynamics.verlet_timestep!",
        "native_multibody_setup" => "Peridynamics.MultibodySetup",
        "native_contact_api" => (with_contact ? "Peridynamics.contact!" : "not called; no-contact control"),
        "point_count_total" => total_points,
        "point_count_per_body" => 8,
        "threads" => Threads.nthreads(),
        "dimensions_m" => [0.004, 0.004, 0.004],
        "dx_m" => DX,
        "horizon_m" => DELTA,
        "point_volume_m3_each" => 8e-9,
        "rho_kg_m3" => RHO,
        "youngs_modulus_Pa" => YOUNG,
        "poisson_ratio" => NU,
        "fracture_energy_J_m2" => GC,
        "failure_mode" => "no_failure! on both bodies",
        "initial_center_a_z_m" => -2.25e-3,
        "initial_center_b_z_m" => 2.25e-3,
        "initial_velocity_a_z_mps" => SPEED,
        "initial_velocity_b_z_mps" => -SPEED,
        "velocity_boundary_conditions" => "none",
        "external_force_density" => "zero; checked at every step",
        "initial_min_point_distance_m" => initial_min_distance,
        "initial_radius_penetration_m" => max(0.0, CONTACT_RADIUS - initial_min_distance),
        "contact_radius_m" => CONTACT_RADIUS,
        "contact_penalty_factor" => PENALTY_FACTOR,
        "native_contact_force_density_formula" => "C=9*penalty_factor/(pi*delta_a^4); native b_int is force density; force computed as native contact b_int increment times own point volume",
        "native_contact_source_sha256" => sha256_file(joinpath(ROOT, "pgfrag", "vendor", "Peridynamics.jl", "src", "physics", "short_range_force_contact.jl")),
        "native_bond_stability_bound_dt_s" => native_bound,
        "native_bond_stability_bound_safety_multiplier" => 1.0,
        "native_bond_bound_scope" => "BB bond-only bound; does not prove penalty-contact stability",
        "dt_factor_of_native_bond_bound" => dt_factor,
        "dt_s" => solver.Δt,
        "steps" => n_steps,
        "end_time_s" => solver.end_time,
        "point_steps_rows" => length(point_rows),
        "metrics_rows" => length(metric_rows),
        "native_contact_neighbor_pairs_at_initialization" => initial_contact_metrics.count,
        "initial_native_contact_min_neighbor_distance_m" => initial_contact_metrics.min_contact_distance,
        "final_min_cross_body_point_distance_m" => isempty(metric_rows) ? initial_min_distance : last_metric[4],
        "final_radius_penetration_m" => isempty(metric_rows) ? 0.0 : last_metric[5],
        "final_total_momentum_kg_mps" => final_momentum,
        "final_total_momentum_drift_norm_kg_mps" => isempty(metric_rows) ? 0.0 : last_metric[11],
        "max_action_reaction_force_error_norm_N" => isempty(metric_rows) ? 0.0 : maximum(r[18] for r in metric_rows),
        "max_action_reaction_relative_error" => isempty(metric_rows) ? 0.0 : maximum(r[19] for r in metric_rows),
        "cumulative_contact_impulse_on_a_Ns" => collect(cumulative_impulse[:a]),
        "cumulative_contact_impulse_on_b_Ns" => collect(cumulative_impulse[:b]),
        "energy" => "UNKNOWN: the native short-range contact API exposes no contact potential energy; no total-energy conservation claim",
        "point_semantics" => "discrete point-to-point penalty contact; not CAD surface penetration and not fragment self-contact",
        "point_id_semantics" => "body_local_point_id is the native per-body point ID; export_point_id is a stable unique 1:16 row key assigned A=1:8 and B=9:16",
        "native_total_b_int_captured_before_bond_diagnostic" => true,
        "diagnostic_state_unchanged_and_restored_every_step" => diagnostic_state_unchanged_all_steps,
        "wall_time_s" => time() - case_start,
    )
    write_json(joinpath(case_dir, "summary.json"), summary)
    println("CASE_COMPLETED case=$case_id points=$total_points steps=$n_steps dt_s=$(solver.Δt) end_time_s=$(solver.end_time) min_distance_m=$(summary["final_min_cross_body_point_distance_m"]) max_AR_error_N=$(summary["max_action_reaction_force_error_norm_N"]) wall_time_s=$(summary["wall_time_s"]) ")
    return (summary=summary, point_rows=point_rows, metric_rows=metric_rows,
            initial_state=initial_state, final_state=final_state)
end

function parse_outdir(args)
    length(args) == 2 && args[1] == "--outdir" || error("usage: julia run_contact_spike.jl --outdir PATH")
    return abspath(args[2])
end

function main(args)
    outdir = parse_outdir(args)
    mkpath(outdir)
    require_ok(Threads.nthreads() == 1, "STOP_RESOURCE_CAP: expected one Julia thread")
    require_ok(Base.pkgversion(Peridynamics) == v"0.5.4+pgfrag.1",
               "STOP_BACKEND_VERSION: expected local headless fork 0.5.4+pgfrag.1")
    println("CONTACT_SPIKE_START julia=$(VERSION) backend=$(Base.pkgversion(Peridynamics)) threads=$(Threads.nthreads()) max_points=$MAX_POINTS max_steps=500 wall_limit_s=$WALL_LIMIT_S")
    println("OUTDIR $outdir")

    baseline_setup, = build_setup(true)
    raw_bound = native_stability_bound(baseline_setup)
    println("NATIVE_BB_BOND_BOUND dt_s=$raw_bound source=Peridynamics.calc_stable_timestep(dh, 1.0) contact_bound=false")
    cases = [
        ("no_contact_control", 0.5, 240, false),
        ("contact_baseline", 0.5, 240, true),
        ("contact_half_dt", 0.25, 480, true),
    ]
    results = Dict{String,Any}[]

    source_contact = joinpath(ROOT, "pgfrag", "vendor", "Peridynamics.jl", "src", "physics", "short_range_force_contact.jl")
    source_verlet = joinpath(ROOT, "pgfrag", "vendor", "Peridynamics.jl", "src", "time_solvers", "velocity_verlet.jl")
    source_multibody = joinpath(ROOT, "pgfrag", "vendor", "Peridynamics.jl", "src", "discretization", "multibody_setup.jl")
    preflight = try
        run_native_preflight(outdir, raw_bound)
    catch err
        preflight_dir = joinpath(outdir, "native_interface_preflight")
        mkpath(preflight_dir)
        message = sprint(showerror, err, catch_backtrace())
        open(joinpath(preflight_dir, "failure.log"), "w") do io
            println(io, message)
        end
        progress_path = joinpath(preflight_dir, "progress.json")
        progress_text = isfile(progress_path) ? read(progress_path, String) : ""
        progress_match = match(r"\"native_steps_completed\"\s*:\s*(\d+)", progress_text)
        completed_native_steps = progress_match === nothing ? 0 : parse(Int, progress_match.captures[1])
        failed_preflight = Dict{String,Any}("status" => "FAILED",
            "steps_completed" => completed_native_steps,
            "native_bond_stability_bound_dt_s" => raw_bound,
            "error" => sprint(showerror, err), "stacktrace" => message,
            "fixed_cases_started" => false)
        write_json(joinpath(preflight_dir, "summary.json"), failed_preflight)
        aggregate = Dict{String,Any}(
            "experiment" => "RelicWeave PG-Frag bounded native point-contact spike",
            "status" => "PREFLIGHT_FAILED",
            "native_backend_version" => string(Base.pkgversion(Peridynamics)),
            "julia_version" => string(VERSION), "threads" => Threads.nthreads(),
            "resource_caps" => Dict("max_points" => MAX_POINTS, "max_steps_per_case" => 500,
                                     "max_wall_seconds_per_case" => WALL_LIMIT_S),
            "native_bond_stability_bound_dt_s" => raw_bound,
            "preflight" => failed_preflight, "cases" => Any[],
            "source_hashes_sha256" => Dict(
                "native_short_range_force_contact_jl" => sha256_file(source_contact),
                "native_velocity_verlet_jl" => sha256_file(source_verlet),
                "native_multibody_setup_jl" => sha256_file(source_multibody),
                "experiment_runner_jl" => sha256_file(@__FILE__)))
        write_json(joinpath(outdir, "report.json"), aggregate)
        println("NATIVE_INTERFACE_PREFLIGHT_FAILED fixed_cases_started=false error=$(sprint(showerror, err))")
        return 2
    end

    consecutive_failures = 0
    for (case_id, factor, steps, with_contact) in cases
        try
            result = run_case(case_id, factor, steps, with_contact, outdir, raw_bound)
            push!(results, result.summary)
            consecutive_failures = 0
        catch err
            consecutive_failures += 1
            backtrace = catch_backtrace()
            case_dir = joinpath(outdir, case_id)
            mkpath(case_dir)
            message = sprint(showerror, err, backtrace)
            open(joinpath(case_dir, "failure.log"), "w") do io
                println(io, message)
            end
            fail_summary = Dict{String,Any}("case_id" => case_id, "status" => "FAILED",
                "error" => sprint(showerror, err), "failure_count_same_root_consecutive" => consecutive_failures,
                "native_bond_stability_bound_dt_s" => raw_bound,
                "dt_factor_of_native_bond_bound" => factor, "steps_requested" => steps,
                "contact" => with_contact, "wall_limit_s" => WALL_LIMIT_S)
            write_json(joinpath(case_dir, "summary.json"), fail_summary)
            push!(results, fail_summary)
            println("CASE_FAILED case=$case_id error=$(sprint(showerror, err)) failures=$consecutive_failures")
            consecutive_failures >= 3 && break
        end
    end

    aggregate = Dict{String,Any}(
        "experiment" => "RelicWeave PG-Frag bounded native point-contact spike",
        "status" => all(x["status"] == "COMPLETED" for x in results) ? "COMPLETED" : "PARTIAL_CASE_FAILURES",
        "native_backend_version" => string(Base.pkgversion(Peridynamics)),
        "julia_version" => string(VERSION),
        "threads" => Threads.nthreads(),
        "resource_caps" => Dict("max_points" => MAX_POINTS, "max_steps_per_case" => 500, "max_wall_seconds_per_case" => WALL_LIMIT_S),
        "native_bond_stability_bound_dt_s" => raw_bound,
        "case_end_time_s" => 120.0 * raw_bound,
        "fixed_contact_parameters" => Dict("radius_m" => CONTACT_RADIUS, "penalty_factor" => PENALTY_FACTOR),
        "native_interface_preflight" => preflight,
        "contact_model_scope" => "native discrete cross-body point-to-point penalty contact; no CAD surfaces and no self-contact",
        "energy_status" => "UNKNOWN: contact potential energy is unavailable; no total-energy conservation claim",
        "source_hashes_sha256" => Dict(
            "native_short_range_force_contact_jl" => sha256_file(source_contact),
            "native_velocity_verlet_jl" => sha256_file(source_verlet),
            "native_multibody_setup_jl" => sha256_file(source_multibody),
            "experiment_runner_jl" => sha256_file(@__FILE__)),
        "cases" => results,
    )
    write_json(joinpath(outdir, "report.json"), aggregate)
    println("CONTACT_SPIKE_DONE status=$(aggregate["status"]) report=$(joinpath(outdir, "report.json"))")
    return all(x["status"] == "COMPLETED" for x in results) ? 0 : 2
end

exit_code = try
    main(ARGS)
catch err
    message = sprint(showerror, err, catch_backtrace())
    println(stderr, message)
    1
end
exit(exit_code)
