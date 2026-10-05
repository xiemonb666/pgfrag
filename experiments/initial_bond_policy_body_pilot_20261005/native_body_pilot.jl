"""Source-only driver for the fixed GEOMETRY_FILTER_FIXTURE_ONLY native graph pilot.

Expected pins are explicit CLI inputs from the just-produced fixture-only transport receipt;
they are not production authorization. The script constructs Body/BondSystem only, never a
Job, BodyChunk, ThreadsBodyDataHandler, time integrator, force evaluation, or CAD workflow.
"""

using Peridynamics

const PILOT_ROOT = @__DIR__
const POINTS_CSV = joinpath(PILOT_ROOT, "fixedcase_points.csv")
const COMPONENTS_CSV = joinpath(PILOT_ROOT, "fixedcase_components.csv")
const HORIZON_M = 0.0018
const VOLUME_M3 = 1.0e-9
const MAX_POINTS = 125
const MAX_DIRECTED_PAIRS = 20_000
const MAX_WALL_SECONDS = 60.0
const MAX_RESULT_BYTES = 20 * 1024 * 1024
const EXTERNAL_ISOLATED_ID = 103

function require(condition::Bool, message::AbstractString)
    condition || error("PILOT_ASSERTION_FAILED: " * String(message))
    return nothing
end

function parse_fixedcase()
    lines = split(chomp(String(read(POINTS_CSV))), '\n')
    require(first(lines) == "point_id,x_mm,y_mm,z_mm", "ordered point CSV header changed")
    rows = Tuple{Int,NTuple{3,Float64}}[]
    for line in lines[2:end]
        fields = split(strip(line), ',')
        require(length(fields) == 4, "fixed point row must have four columns")
        push!(rows, (parse(Int, fields[1]),
                     (parse(Float64, fields[2]), parse(Float64, fields[3]),
                      parse(Float64, fields[4]))))
    end
    component_lines = split(chomp(String(read(COMPONENTS_CSV))), '\n')
    require(first(component_lines) == "point_id,component", "component fixture header changed")
    components = Dict{Int,String}()
    for line in component_lines[2:end]
        fields = split(strip(line), ',')
        require(length(fields) == 2, "component row must have two columns")
        components[parse(Int, fields[1])] = fields[2]
    end
    require(1 <= length(rows) <= MAX_POINTS, "fixture point resource cap failed")
    require(length(components) == length(rows), "component mapping does not cover all points")
    ids = first.(rows)
    require(length(unique(ids)) == length(ids), "point IDs are not unique")
    return rows, components
end

function make_body(rows, policy)
    n = length(rows)
    position = Matrix{Float64}(undef, 3, n)
    for (native_index, (_, xyz_mm)) in enumerate(rows), axis in 1:3
        position[axis, native_index] = xyz_mm[axis] * 0.001
    end
    body = Peridynamics.Body(Peridynamics.OSBMaterial(), position,
                             fill(VOLUME_M3, n); initial_bond_policy=policy)
    Peridynamics.material!(body; horizon=HORIZON_M, rho=1000.0, E=1.0e6, nu=0.25)
    require(body.mat isa Peridynamics.OSBMaterial, "the fixture must use native OSBMaterial")
    require(body.mat.kernel === Peridynamics.linear_kernel,
            "weighted-volume fixture must use the pinned default linear kernel")
    require(isempty(body.point_sets_precracks), "the fixed fixture must have no precracks")
    require(body.n_points == n && size(body.position) == (3, n), "native Body shape mismatch")
    return body
end

function expected_decisions(rows, components; isolate_id::Union{Nothing,Int}=nothing)
    ids = first.(rows)
    xyz = [collect(row[2]) .* 0.001 for row in rows]
    decisions = Dict{Tuple{Int,Int},String}()
    directed = Set{Tuple{Int,Int}}()
    source_degree = zeros(Int, length(rows))
    for i in eachindex(rows), j in eachindex(rows)
        i == j && continue
        distance2 = sum(abs2, xyz[i] .- xyz[j])
        # PointNeighbors GridNeighborhoodSearch v0.6.7 uses distance2 <= radius^2.
        # Fixture coordinates deliberately place no pair exactly at the horizon.
        if distance2 <= HORIZON_M^2
            push!(directed, (ids[i], ids[j]))
            source_degree[i] += 1
            if i < j
                pair = ids[i] < ids[j] ? (ids[i], ids[j]) : (ids[j], ids[i])
                state = components[pair[1]] == components[pair[2]] ? "IN_MATERIAL" : "OUT_OR_VOID"
                if isolate_id !== nothing && isolate_id in pair && state == "IN_MATERIAL"
                    state = "OUT_OR_VOID"
                end
                decisions[pair] = state
            end
        end
    end
    require(length(directed) <= MAX_DIRECTED_PAIRS, "directed candidate resource cap failed")
    return decisions, directed, source_degree
end

function directed_raw_native(bonds, degrees, native_to_external)
    length(degrees) == length(native_to_external) || error("raw degree vector has wrong length")
    observed = Set{Tuple{Int,Int}}()
    offset = 0
    for source_native in eachindex(degrees)
        count = degrees[source_native]
        for bond_index in (offset + 1):(offset + count)
            target_native = bonds[bond_index].neighbor
            require(1 <= target_native <= length(native_to_external), "raw target index is invalid")
            require(source_native != target_native, "raw native roster contains a self-bond")
            direction = (native_to_external[source_native], native_to_external[target_native])
            require(!(direction in observed), "raw native roster has a duplicate direction")
            push!(observed, direction)
        end
        offset += count
    end
    require(offset == length(bonds), "raw neighbor counts do not cover the native ledger")
    return observed
end

function read_native_graph(system, native_to_external)
    n = length(native_to_external)
    require(system.chunk_handler.n_loc_points == n, "native graph is not a complete one-chunk body")
    require(system.chunk_handler.point_ids == collect(1:n), "one chunk point order changed")
    for native_id in 1:n
        require(system.chunk_handler.localizer[native_id] == native_id,
                "one-chunk native localization is not identity")
        require(length(system.bond_ids[native_id]) == system.n_neighbors[native_id],
                "bond ID range length disagrees with native degree")
    end
    directed = Set{Tuple{Int,Int}}()
    pair_ids = Dict{Tuple{Int,Int},Int}()
    for source_native in 1:n
        for bond_id in system.bond_ids[source_native]
            bond = system.bonds[bond_id]
            target_native = bond.neighbor
            source_external = native_to_external[source_native]
            target_external = native_to_external[target_native]
            direction = (source_external, target_external)
            require(!(direction in directed), "post-mask native graph repeats a directed edge")
            push!(directed, direction)
            pair_ids[direction] = bond_id
        end
    end
    require(length(pair_ids) == length(system.bonds), "native pair IDs are not a complete ledger")
    require(sort(collect(values(pair_ids))) == collect(1:length(system.bonds)),
            "native pair IDs are not a unique contiguous ledger")
    return directed, pair_ids
end

function expected_undirected(states, wanted_state)
    return Set(pair for (pair, state) in states if state == wanted_state)
end

function make_readonly_osb_state(system)
    # The native OSB weighted-volume method reads only `bond_active` from its storage argument.
    # This fresh all-active scratch value invokes that method without constructing a solver,
    # BodyChunk, force field, or physical event state.
    n = size(system.position, 2)
    nloc = system.chunk_handler.n_loc_points
    nbonds = length(system.bonds)
    return Peridynamics.OSBStorage(
        copy(system.position),
        zeros(3, n), zeros(3, n), zeros(3, n), zeros(3, n), zeros(3, n),
        zeros(3, n), zeros(3, n), zeros(3, n), zeros(3, n),
        zeros(nloc), zeros(Int, nloc), zeros(nloc),
        zeros(nbonds), trues(nbonds), zeros(nloc),
        zeros(3, nloc), zeros(3, nloc), zeros(3), zeros(3), zeros(3), zeros(3),
    )
end

function native_weighted_volume(system, body, storage)
    return [Peridynamics.calc_weighted_volume(
                storage, system, body.mat, Peridynamics.get_point_param(body, native_i), native_i)
            for native_i in 1:body.n_points]
end

function independent_weighted_volume(rows, included_pairs, body)
    expected_values = zeros(Float64, length(rows))
    included = Set(included_pairs)
    for (source_native, source_row) in enumerate(rows)
        source_id = source_row[1]
        delta = Peridynamics.get_point_param(body, source_native).δ
        for (target_native, target_row) in enumerate(rows)
            source_native == target_native && continue
            target_id = target_row[1]
            pair = source_id < target_id ? (source_id, target_id) : (target_id, source_id)
            pair in included || continue
            dx = body.position[:, source_native] .- body.position[:, target_native]
            distance2 = sum(abs2, dx)
            distance = sqrt(distance2)
            expected_values[source_native] += body.mat.kernel(delta, distance) * distance2 *
                                             body.volume[target_native]
        end
    end
    return expected_values
end

function validate_good_policy(json_path, toml_path, json_sha, toml_sha, rows, components)
    policy = Peridynamics.InitialBondPolicy(json_path, toml_path;
        expected_json_sha256=json_sha, expected_toml_sha256=toml_sha)
    model = Peridynamics._initial_policy_current_model(policy)
    require(model.eligibility == "CANDIDATE_MASK_POSITIVE_DEGREE_NATIVE_UNVERIFIED",
            "positive fixture policy has a blocked or verified status")
    require(isempty(model.zero_degree_ids), "positive fixture policy has a zero-degree row")
    expected_states, expected_raw, expected_raw_degrees = expected_decisions(rows, components)
    require(model.decision_by_pair == expected_states,
            "transport ledger differs from the independently reconstructed fixture rule")
    require(length(expected_raw) == 120 && length(expected_states) == 60,
            "fixedcase expected raw roster changed from 120 directed / 60 canonical pairs")
    body = make_body(rows, policy)
    for i in 1:body.n_points
        δ = Peridynamics.get_point_param(body, :δ, i)
        require(isfinite(δ) && δ == model.horizon_m && δ == HORIZON_M,
                "native material horizons are not exactly uniform")
        for axis in 1:3
            expected_m = rows[i][2][axis] * 0.001
            require(isapprox(body.position[axis, i], expected_m; atol=1.0e-12, rtol=0.0),
                    "native Body position does not match ordered CSV mm-to-m mapping")
        end
    end

    raw_bonds, raw_degrees = Peridynamics.find_bonds(body, collect(1:body.n_points))
    raw_native_to_external = first.(rows)
    raw_observed = directed_raw_native(raw_bonds, raw_degrees, raw_native_to_external)
    require(raw_observed == expected_raw,
            "actual native unmasked candidate roster differs from independent Euclidean rule")
    require(raw_degrees == expected_raw_degrees,
            "actual native raw candidate degrees differ from independent Euclidean degrees")

    point_decomposition = Peridynamics.PointDecomposition(body, 1)
    require(point_decomposition.n_chunks == 1 && point_decomposition.decomp == [1:body.n_points],
            "native point decomposition is not one ordered chunk")
    system = Peridynamics.BondSystem(body, point_decomposition, 1)
    require(system.correction isa Peridynamics.NoCorrection,
            "weighted-volume fixture must use native NoCorrection")
    receipt = system.initial_bond_admission
    require(receipt !== nothing, "native BondSystem omitted the admission receipt")
    require(receipt.raw_candidate_set_reproduced && receipt.filtered_bond_list_reproduced,
            "native seam did not reproduce raw candidates and filtered roster")
    require(receipt.geometry_exclusions_are_physics_events == false,
            "geometry exclusions were mislabeled as physics events")
    require(receipt.raw_directed_bond_count == length(raw_bonds),
            "admission raw count differs from actual unmasked native roster")
    require(receipt.included_directed_bond_count == 112 && length(system.bonds) == 112,
            "native post-mask roster is not 112 directed edges")
    require(receipt.raw_candidate_pair_count == 60 && receipt.included_pair_count == 56 &&
            receipt.geometry_excluded_pair_count == 4 && receipt.unknown_pair_count == 0,
            "native admission receipt pair counts changed")

    graph_directed, pair_ids = read_native_graph(system, raw_native_to_external)
    included_pairs = expected_undirected(expected_states, "IN_MATERIAL")
    excluded_pairs = expected_undirected(expected_states, "OUT_OR_VOID")
    expected_graph = Set{Tuple{Int,Int}}()
    for (left, right) in included_pairs
        push!(expected_graph, (left, right))
        push!(expected_graph, (right, left))
    end
    require(graph_directed == expected_graph,
            "native graph readback differs from complete bidirectional included mask")
    for (left, right) in excluded_pairs
        require(!((left, right) in graph_directed) && !((right, left) in graph_directed),
                "a geometry-excluded pair remains in one or both native directions")
    end
    degrees_from_graph = zeros(Int, body.n_points)
    for (source_external, _) in graph_directed
        degrees_from_graph[model.external_to_native[source_external]] += 1
    end
    require(degrees_from_graph == system.n_neighbors && all(==(7), system.n_neighbors),
            "native post-mask degree readback differs from K8+K8 fixture graph")

    readonly_storage = make_readonly_osb_state(system)
    native_wvol = native_weighted_volume(system, body, readonly_storage)
    expected_wvol = independent_weighted_volume(rows, collect(included_pairs), body)
    require(all(value -> isfinite(value) && value > 0.0, native_wvol),
            "native read-only OSB weighted volumes must be finite and positive")
    for i in eachindex(native_wvol)
        require(isapprox(native_wvol[i], expected_wvol[i]; atol=0.0, rtol=1.0e-12),
                "native OSB weighted-volume method differs from independent post-mask formula")
    end
    require(receipt.native_graph_readback == "NOT_RUN",
            "static policy transport was mutated to claim readback evidence")
    return (policy, system, receipt, raw_bonds, graph_directed, pair_ids,
            native_wvol, expected_wvol, model, included_pairs, excluded_pairs)
end

function validate_zero_degree_block(json_path, toml_path, json_sha, toml_sha, rows, components)
    policy = Peridynamics.InitialBondPolicy(json_path, toml_path;
        expected_json_sha256=json_sha, expected_toml_sha256=toml_sha)
    model = Peridynamics._initial_policy_current_model(policy)
    require(model.eligibility == "BLOCKED_ZERO_DEGREE" &&
            model.zero_degree_ids == [EXTERNAL_ISOLATED_ID],
            "negative fixture policy did not declare exactly one zero-degree row")
    expected_states, _, _ = expected_decisions(rows, components; isolate_id=EXTERNAL_ISOLATED_ID)
    require(model.decision_by_pair == expected_states,
            "zero-degree transport state ledger differs from the probe definition")
    body = make_body(rows, policy)
    point_decomposition = Peridynamics.PointDecomposition(body, 1)
    rejected = false
    try
        Peridynamics.BondSystem(body, point_decomposition, 1)
    catch err
        message = lowercase(sprint(showerror, err))
        rejected = err isa ArgumentError && occursin("blocked", message)
    end
    require(rejected, "blocked zero-degree policy was not rejected before native BondSystem admission")
    return model
end

function main()
    started = time()
    length(ARGS) == 9 || error(
        "usage: native_body_pilot.jl good.json good.toml good-json-sha good-toml-sha " *
        "blocked.json blocked.toml blocked-json-sha blocked-toml-sha result.txt")
    good_json, good_toml, good_json_sha, good_toml_sha = ARGS[1:4]
    blocked_json, blocked_toml, blocked_json_sha, blocked_toml_sha = ARGS[5:8]
    result_path = abspath(ARGS[9])
    result_dir = dirname(abspath(good_json))
    require(dirname(result_path) == result_dir,
            "result receipt must remain beside the disposable positive fixture bundle")
    require(!ispath(result_path), "refusing to overwrite a prior native pilot result")

    rows, components = parse_fixedcase()
    require(rows[1][1] == 100 && rows[8][1] == EXTERNAL_ISOLATED_ID,
            "ordered external/native ID mapping changed")
    good = validate_good_policy(good_json, good_toml, good_json_sha, good_toml_sha,
                                rows, components)
    blocked_model = validate_zero_degree_block(blocked_json, blocked_toml,
                                               blocked_json_sha, blocked_toml_sha,
                                               rows, components)
    elapsed = time() - started
    require(elapsed <= MAX_WALL_SECONDS, "native graph fixture exceeded the wall-time budget")

    _, system, admission, raw_bonds, graph_directed, pair_ids, native_wvol,
    expected_wvol, model, included_pairs, excluded_pairs = good
    result = Dict{String,Any}(
        "label" => "GEOMETRY_FILTER_FIXTURE_ONLY",
        "source_artifact" => "fixed synthetic point cloud; no CAD input or accuracy claim",
        "native_body" => "OSBMaterial single Body constructed",
        "native_graph" => "BondSystem read back; one PointDecomposition chunk; ThreadsBodyDataHandler not constructed",
        "native_candidates" => "raw native roster exactly matched independent Euclidean <= rule",
        "raw_canonical_pair_count" => admission.raw_candidate_pair_count,
        "raw_directed_bond_count" => length(raw_bonds),
        "ordered_point_payload_sha256" => model.points["payload_sha256"],
        "audit_contract_json_sha256" => admission.audit_contract_json_sha256,
        "native_transport_toml_sha256" => admission.native_transport_toml_sha256,
        "source_marker_sha256" => admission.source_cad_sha256,
        "included_canonical_pair_count" => length(included_pairs),
        "geometry_excluded_canonical_pair_count" => length(excluded_pairs),
        "included_directed_bond_count" => length(graph_directed),
        "zero_degree_probe" => "rejected before BondSystem admission",
        "pair_id_ledger_count" => length(pair_ids),
        "weighted_volume" => "native calc_weighted_volume direct read-only formula matched independent post-mask recomputation",
        "native_weighted_volume_readback" => native_wvol,
        "independent_weighted_volume_recalculation" => expected_wvol,
        "policy_transport_native_graph_readback" => admission.native_graph_readback,
        "geometry_exclusions_are_physics_events" => admission.geometry_exclusions_are_physics_events,
        "physical_failure_or_fracture" => "NOT_RUN",
        "solver_force_or_time_step" => "NOT_CREATED_OR_RUN",
        "body_chunk_or_threads_handler" => "NOT_CREATED",
        "CAD_verification" => "NOT_RUN",
        "blocked_probe_eligibility" => blocked_model.eligibility,
        "elapsed_seconds" => elapsed,
        "limits" => "<=125 points; <=20000 directed raw edges; <=60 s; output <=20 MiB",
    )
    mkpath(result_dir)
    open(result_path, "w") do stream
        for key in sort(collect(keys(result)))
            value = result[key]
            rendered = value isa AbstractVector ? join(string.(value), ",") : string(value)
            println(stream, key, "=", rendered)
        end
    end
    require(filesize(result_path) <= MAX_RESULT_BYTES, "result receipt exceeds the output byte cap")
    println("PILOT_RESULT_TEXT=" * result_path)
    println("PILOT_LABEL=GEOMETRY_FILTER_FIXTURE_ONLY")
    println("PILOT_RAW_DIRECTED=" * string(length(raw_bonds)))
    println("PILOT_INCLUDED_DIRECTED=" * string(length(graph_directed)))
    println("PILOT_GRAPH_READBACK=PASS")
    println("PILOT_PHYSICAL_FAILURE=NOT_RUN")
    println("PILOT_CAD_VERIFICATION=NOT_RUN")
    return nothing
end

main()
