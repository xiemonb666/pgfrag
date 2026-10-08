const INITIAL_POLICY_AUDIT_MAX_BYTES = 64 * 1024 * 1024
const INITIAL_POLICY_TOML_MAX_BYTES = 128 * 1024 * 1024
const INITIAL_POLICY_TOTAL_INPUT_MAX_BYTES = 256 * 1024 * 1024
const INITIAL_POLICY_ARTIFACT_TOTAL_MAX_BYTES = 64 * 1024 * 1024
const INITIAL_POLICY_CAD_MAX_BYTES = 16 * 1024 * 1024
const INITIAL_POLICY_RESULT_MAX_BYTES = 8 * 1024 * 1024
const INITIAL_POLICY_CSV_MAX_BYTES = 2 * 1024 * 1024
const INITIAL_POLICY_SOURCE_MAX_BYTES = 2 * 1024 * 1024
const INITIAL_POLICY_SIDECAR_MAX_BYTES = 32 * 1024 * 1024
const INITIAL_POLICY_MAX_POINTS = 5000
const INITIAL_POLICY_MAX_PAIRS = 250_000
const INITIAL_POLICY_MAX_EXTERNAL_ID = typemax(Int64)
const INITIAL_POLICY_MAX_REASON_BYTES = 4096
const INITIAL_POLICY_MAX_PATH_BYTES = 16_384
const INITIAL_POLICY_POSITION_ATOL_M = 1.0e-12
const INITIAL_POLICY_STATES = ("IN_MATERIAL", "OUT_OR_VOID", "UNKNOWN")
const INITIAL_POLICY_ELIGIBILITIES = (
    "BLOCKED_UNKNOWN",
    "BLOCKED_EMPTY_CANDIDATE_SET",
    "BLOCKED_EMPTY_MASK",
    "BLOCKED_ZERO_DEGREE",
    "CANDIDATE_MASK_POSITIVE_DEGREE_NATIVE_UNVERIFIED",
)

struct _InitialBondPolicyConstructionToken end
const _INITIAL_BOND_POLICY_TOKEN = _InitialBondPolicyConstructionToken()

"""
An immutable, externally pinned TOML projection of one strict JSON mask contract.

The expected digest keywords must come from the trusted job/compiler receipt. A digest
embedded only in the TOML is not sufficient to establish the identity of that TOML.
"""
struct InitialBondPolicy
    audit_json_path::String
    native_toml_path::String
    native_toml_bytes::String
    audit_json_sha256::String
    native_toml_sha256::String

    function InitialBondPolicy(audit_json_path::String, native_toml_path::String,
                               native_toml_bytes::String, audit_json_sha256::String,
                               native_toml_sha256::String,
                               ::_InitialBondPolicyConstructionToken)
        return new(audit_json_path, native_toml_path, native_toml_bytes,
                   audit_json_sha256, native_toml_sha256)
    end
end

"""
Admission evidence for the geometry filter at the raw-bond construction seam.

This receipt does not claim physical fracture events, correction/OSB initialization,
full BodyChunk readback, or authorization to start a simulation.
"""
struct InitialBondAdmissionReceipt
    schema::String
    status::String
    audit_contract_json_sha256::String
    native_transport_toml_sha256::String
    source_cad_sha256::String
    point_payload_sha256::String
    classifier_result_sha256::String
    classifier_source_sha256::String
    raw_hit_interval_sidecar_sha256::String
    pair_map_sha256::String
    included_mask_sha256::String
    raw_candidate_pair_count::Int
    raw_directed_bond_count::Int
    included_pair_count::Int
    included_directed_bond_count::Int
    geometry_excluded_pair_count::Int
    unknown_pair_count::Int
    raw_candidate_set_reproduced::Bool
    filtered_bond_list_reproduced::Bool
    geometry_exclusions_are_physics_events::Bool
    native_graph_readback::String
    initialization_authorized::Bool
end

initial_bond_policy(::AbstractBody) = nothing

function _initial_policy_require(condition::Bool, message::AbstractString)
    condition || throw(ArgumentError(String(message)))
    return nothing
end

function _initial_policy_file_metadata(path::String)
    return (filesize(path), mtime(path))
end

function _read_initial_policy_bundle(
    specs::Vector{Tuple{String,String,Int}};
    max_total_bytes::Int,
)
    isempty(specs) && throw(ArgumentError("initial bond policy input bundle is empty"))
    prepared = Dict{String,Tuple{String,Int,Tuple{Int64,Float64}}}()
    resolved_paths = Dict{String,String}()
    seen_paths = Set{String}()
    total_bytes = Int64(0)

    for (label, raw_path, limit) in specs
        resolved = try
            realpath(raw_path)
        catch
            throw(ArgumentError("cannot resolve initial bond policy input: $(label)"))
        end
        isfile(resolved) || throw(ArgumentError("$(label) must be a regular file"))
        metadata = try
            _initial_policy_file_metadata(resolved)
        catch
            throw(ArgumentError("cannot stat initial bond policy input: $(label)"))
        end
        size = metadata[1]
        (size >= 0 && size <= limit) ||
            throw(ArgumentError("$(label) exceeds its byte limit"))
        resolved in seen_paths &&
            throw(ArgumentError("initial bond policy roles must name distinct resolved paths"))
        push!(seen_paths, resolved)
        total_bytes += size
        total_bytes <= max_total_bytes ||
            throw(ArgumentError("initial bond policy inputs exceed the aggregate byte limit"))
        prepared[label] = (resolved, Int(size), metadata)
        resolved_paths[label] = resolved
    end

    function read_once(path::String, size::Int, expected, label::String)
        try
            isfile(path) || throw(ArgumentError("$(label) ceased to be a regular file"))
            _initial_policy_file_metadata(path) == expected ||
                throw(ArgumentError("$(label) changed after size preflight"))
            bytes = open(path, "r") do stream
                read(stream, size + 1)
            end
            length(bytes) == size ||
                throw(ArgumentError("$(label) changed size during bounded read"))
            _initial_policy_file_metadata(path) == expected ||
                throw(ArgumentError("$(label) changed during bounded read"))
            return bytes
        catch err
            err isa ArgumentError && rethrow()
            throw(ArgumentError("cannot bounded-read initial bond policy input: $(label)"))
        end
    end

    captured = Dict{String,Vector{UInt8}}()
    for (label, (path, size, metadata)) in prepared
        first = read_once(path, size, metadata, label)
        second = read_once(path, size, metadata, label)
        first == second || throw(ArgumentError("$(label) changed between bounded reads"))
        captured[label] = first
    end
    return captured, resolved_paths
end

function _initial_policy_sha256(bytes::AbstractVector{UInt8})
    return bytes2hex(SHA.sha256(bytes))
end

function _initial_policy_sha256_ok(value)
    return value isa String && ncodeunits(value) == 64 &&
           all(c -> ('0' <= c <= '9') || ('a' <= c <= 'f'), value)
end

function _initial_policy_sha(value, label::AbstractString)
    _initial_policy_require(_initial_policy_sha256_ok(value),
                            "$(label) must be a lowercase SHA-256 digest")
    return value::String
end

function _initial_policy_string(value, label::AbstractString; max_bytes::Int=INITIAL_POLICY_MAX_PATH_BYTES)
    _initial_policy_require(value isa String && !isempty(value), "$(label) must be a non-empty string")
    _initial_policy_require(ncodeunits(value) <= max_bytes, "$(label) exceeds its string byte limit")
    return value::String
end

function _initial_policy_int(value, label::AbstractString;
                             minimum::Int64=0, maximum::Int64=INITIAL_POLICY_MAX_EXTERNAL_ID)
    _initial_policy_require(value isa Integer && !(value isa Bool),
                            "$(label) must be an integer")
    _initial_policy_require(minimum <= value <= maximum,
                            "$(label) is outside its allowed integer range")
    return Int(value)
end

function _initial_policy_float(value, label::AbstractString;
                               minimum::Float64=-Inf, strict_minimum::Bool=false)
    _initial_policy_require(value isa Float64 && isfinite(value),
                            "$(label) must be a finite TOML float")
    _initial_policy_require(strict_minimum ? value > minimum : value >= minimum,
                            "$(label) is outside its allowed numeric range")
    return value::Float64
end

function _initial_policy_bool(value, label::AbstractString)
    _initial_policy_require(value isa Bool, "$(label) must be a TOML boolean")
    return value::Bool
end

function _initial_policy_array(value, label::AbstractString; max_items::Int=INITIAL_POLICY_MAX_PAIRS)
    _initial_policy_require(value isa AbstractVector, "$(label) must be a TOML array")
    _initial_policy_require(length(value) <= max_items, "$(label) exceeds its item cap")
    return value
end

function _initial_policy_table(value, required::Tuple, label::AbstractString)
    _initial_policy_require(value isa AbstractDict, "$(label) must be a TOML table")
    table = value::AbstractDict
    all(key -> key isa String, keys(table)) ||
        throw(ArgumentError("$(label) contains a non-string TOML key"))
    Set(keys(table)) == Set(required) ||
        throw(ArgumentError("$(label) has missing or unsupported keys"))
    return table
end

function _initial_policy_value_tree(value, label::AbstractString)
    if value isa AbstractDict
        for (key, child) in value
            key isa String || throw(ArgumentError("$(label) has a non-string key"))
            _initial_policy_value_tree(child, "$(label).$(key)")
        end
    elseif value isa AbstractVector
        for (index, child) in enumerate(value)
            _initial_policy_value_tree(child, "$(label)[$(index)]")
        end
    elseif value isa AbstractFloat
        isfinite(value) || throw(ArgumentError("$(label) contains NaN or infinity"))
    elseif value isa Integer || value isa Bool || value isa String
        return nothing
    else
        throw(ArgumentError("$(label) contains an unsupported TOML value"))
    end
    return nothing
end

function _initial_policy_transform(value)
    rows = _initial_policy_array(value, "coordinates.native_frame_from_reference_frame_m"; max_items=4)
    length(rows) == 4 || throw(ArgumentError("native frame transform must be 4x4"))
    matrix = Matrix{Float64}(undef, 4, 4)
    for row_index in 1:4
        row = _initial_policy_array(rows[row_index], "native frame transform row"; max_items=4)
        length(row) == 4 || throw(ArgumentError("native frame transform must be 4x4"))
        for column_index in 1:4
            matrix[row_index, column_index] =
                _initial_policy_float(row[column_index], "native frame transform element")
        end
    end
    matrix[4, :] == [0.0, 0.0, 0.0, 1.0] ||
        throw(ArgumentError("native frame transform must be homogeneous affine form"))
    for i in 1:3, j in 1:3
        dot = sum(matrix[k, i] * matrix[k, j] for k in 1:3)
        expected = i == j ? 1.0 : 0.0
        isapprox(dot, expected; atol=1.0e-12, rtol=0.0) ||
            throw(ArgumentError("native frame transform rotation must be orthonormal"))
    end
    rotation = matrix[1:3, 1:3]
    determinant = (
        rotation[1,1] * (rotation[2,2] * rotation[3,3] - rotation[2,3] * rotation[3,2]) -
        rotation[1,2] * (rotation[2,1] * rotation[3,3] - rotation[2,3] * rotation[3,1]) +
        rotation[1,3] * (rotation[2,1] * rotation[3,2] - rotation[2,2] * rotation[3,1])
    )
    isapprox(determinant, 1.0; atol=1.0e-12, rtol=0.0) ||
        throw(ArgumentError("native frame transform must preserve orientation"))
    return matrix
end

function _initial_policy_point_order_sha(point_ids::Vector{Int})
    io = IOBuffer()
    write(io, "pgfrag.point-ids.v1\n")
    for point_id in point_ids
        print(io, point_id, '\n')
    end
    return _initial_policy_sha256(take!(io))
end

function _initial_policy_pair_map_sha(decisions::Vector{Tuple{Int,Int,String}})
    io = IOBuffer()
    write(io, "pgfrag.native-pair-map.v1\n")
    for (left, right, state) in decisions
        print(io, left, ',', right, ',', state, '\n')
    end
    return _initial_policy_sha256(take!(io))
end

function _initial_policy_mask_sha(pairs::Vector{Tuple{Int,Int}})
    io = IOBuffer()
    write(io, "pgfrag.native-included-mask.v1\n")
    for (left, right) in pairs
        print(io, left, ',', right, '\n')
    end
    return _initial_policy_sha256(take!(io))
end

function _initial_policy_pair(row, label::AbstractString, point_id_set::Set{Int})
    left = _initial_policy_int(row["point_i"], "$(label).point_i")
    right = _initial_policy_int(row["point_j"], "$(label).point_j")
    left < right && left in point_id_set && right in point_id_set ||
        throw(ArgumentError("$(label) must use canonical external IDs i<j present in the point rows"))
    return (left, right)
end

function _initial_policy_validate_transport(model, audit_json_sha256::String,
                                            audit_json_path::String)
    _initial_policy_value_tree(model, "transport")
    root = _initial_policy_table(model,
        ("transport", "source", "coordinates", "points", "candidate_generation",
         "raw_euclidean_candidate_pairs", "counts", "hashes", "mask", "evidence",
         "native_gate", "initial_physics_events"), "root")

    transport = _initial_policy_table(root["transport"],
        ("schema", "audit_contract_json_path", "audit_contract_json_sha256",
         "static_contract_status", "native_mask_eligibility", "qualification_status",
         "initialization_authorized"), "transport")
    transport["schema"] == "pgfrag.native-material-visibility-toml.v1" ||
        throw(ArgumentError("unsupported native TOML transport schema"))
    transport["static_contract_status"] == "VALID" ||
        throw(ArgumentError("static audit contract is not VALID"))
    transport["qualification_status"] == "TRANSPORT_ONLY_NATIVE_UNVERIFIED" ||
        throw(ArgumentError("native transport cannot claim native verification"))
    _initial_policy_bool(transport["initialization_authorized"], "transport.initialization_authorized") == false ||
        throw(ArgumentError("native transport cannot authorize initialization"))
    eligibility = _initial_policy_string(transport["native_mask_eligibility"], "transport.native_mask_eligibility"; max_bytes=128)
    eligibility in INITIAL_POLICY_ELIGIBILITIES ||
        throw(ArgumentError("unsupported native mask eligibility state"))
    embedded_json_path = _initial_policy_string(transport["audit_contract_json_path"], "transport.audit_contract_json_path")
    realpath(embedded_json_path) == audit_json_path ||
        throw(ArgumentError("audit JSON path does not match the path embedded in TOML"))
    embedded_json_sha = _initial_policy_sha(transport["audit_contract_json_sha256"], "transport.audit_contract_json_sha256")
    embedded_json_sha == audit_json_sha256 ||
        throw(ArgumentError("TOML audit JSON SHA does not match the actual pinned JSON bytes"))

    source = _initial_policy_table(root["source"],
        ("source_cad_path", "source_cad_sha256", "classifier_source_path",
         "classifier_source_sha256", "source_hash_status", "source_status",
         "source_accuracy_mm", "source_accuracy_basis"), "source")
    _initial_policy_string(source["source_cad_path"], "source.source_cad_path")
    _initial_policy_sha(source["source_cad_sha256"], "source.source_cad_sha256")
    _initial_policy_string(source["classifier_source_path"], "source.classifier_source_path")
    _initial_policy_sha(source["classifier_source_sha256"], "source.classifier_source_sha256")
    source["source_hash_status"] == "MATCHED_ACTUAL_FILE_BYTES" ||
        throw(ArgumentError("source CAD is not byte-SHA-bound"))
    source["source_status"] == "READY" ||
        throw(ArgumentError("source CAD is not READY"))
    _initial_policy_float(source["source_accuracy_mm"], "source.source_accuracy_mm"; minimum=0.0)
    _initial_policy_string(source["source_accuracy_basis"], "source.source_accuracy_basis"; max_bytes=4096)

    coordinates = _initial_policy_table(root["coordinates"],
        ("input_unit", "native_unit", "m_per_input_unit", "reference_frame",
         "native_frame", "native_frame_from_reference_frame_m"), "coordinates")
    coordinates["input_unit"] == "mm" && coordinates["native_unit"] == "m" &&
        coordinates["m_per_input_unit"] == 0.001 ||
        throw(ArgumentError("coordinates must specify exactly one mm-to-m conversion"))
    _initial_policy_string(coordinates["reference_frame"], "coordinates.reference_frame"; max_bytes=512)
    _initial_policy_string(coordinates["native_frame"], "coordinates.native_frame"; max_bytes=512)
    transform = _initial_policy_transform(coordinates["native_frame_from_reference_frame_m"])

    points = _initial_policy_table(root["points"],
        ("payload_path", "payload_sha256", "order", "point_ids_order_sha256", "count", "rows"),
        "points")
    point_path = _initial_policy_string(points["payload_path"], "points.payload_path")
    point_sha = _initial_policy_sha(points["payload_sha256"], "points.payload_sha256")
    points["order"] == "exact_payload_row_order" ||
        throw(ArgumentError("point rows are not declared in payload order"))
    point_order_sha = _initial_policy_sha(points["point_ids_order_sha256"], "points.point_ids_order_sha256")
    point_count = _initial_policy_int(points["count"], "points.count";
                                      minimum=1, maximum=INITIAL_POLICY_MAX_POINTS)
    point_rows = _initial_policy_array(points["rows"], "points.rows"; max_items=INITIAL_POLICY_MAX_POINTS)
    length(point_rows) == point_count ||
        throw(ArgumentError("points.count disagrees with points.rows"))
    point_ids = Int[]
    external_to_native = Dict{Int,Int}()
    native_to_external = Vector{Int}(undef, point_count)
    positions_mm = Vector{NTuple{3,Float64}}(undef, point_count)
    for (row_index, raw_row) in enumerate(point_rows)
        row = _initial_policy_table(raw_row, ("point_id", "native_index", "position_mm"),
                                    "points.rows[$(row_index)]")
        point_id = _initial_policy_int(row["point_id"], "points.rows[$(row_index)].point_id")
        native_index = _initial_policy_int(row["native_index"], "points.rows[$(row_index)].native_index";
                                           minimum=1, maximum=INITIAL_POLICY_MAX_POINTS)
        native_index == row_index ||
            throw(ArgumentError("native_index must exactly equal the one-based payload row"))
        haskey(external_to_native, point_id) &&
            throw(ArgumentError("point IDs must be unique"))
        position = _initial_policy_array(row["position_mm"], "points.rows.position_mm"; max_items=3)
        length(position) == 3 || throw(ArgumentError("each point position must have three values"))
        xyz = ntuple(i -> _initial_policy_float(position[i], "point coordinate"), 3)
        external_to_native[point_id] = native_index
        native_to_external[native_index] = point_id
        push!(point_ids, point_id)
        positions_mm[native_index] = xyz
    end
    _initial_policy_point_order_sha(point_ids) == point_order_sha ||
        throw(ArgumentError("point_ids_order_sha256 does not match points.rows order"))
    point_id_set = Set(point_ids)

    generation = _initial_policy_table(root["candidate_generation"],
        ("metric", "horizon_m", "candidate_rule", "scope"), "candidate_generation")
    generation["metric"] == "euclidean_reference_distance" &&
        generation["candidate_rule"] == "exact_pinned_Peridynamics_find_bonds_convention" &&
        generation["scope"] == "one_body_one_uniform_horizon" ||
        throw(ArgumentError("unsupported candidate-generation convention"))
    horizon_m = _initial_policy_float(generation["horizon_m"], "candidate_generation.horizon_m";
                                      minimum=0.0, strict_minimum=true)

    raw_candidates = _initial_policy_array(root["raw_euclidean_candidate_pairs"],
        "raw_euclidean_candidate_pairs"; max_items=INITIAL_POLICY_MAX_PAIRS)
    decisions = Tuple{Int,Int,String}[]
    decision_by_pair = Dict{Tuple{Int,Int},String}()
    previous_pair = nothing
    for (index, raw_row) in enumerate(raw_candidates)
        row = _initial_policy_table(raw_row,
            ("point_i", "point_j", "native_i", "native_j", "classification", "reason"),
            "raw_euclidean_candidate_pairs[$(index)]")
        pair = _initial_policy_pair(row, "raw candidate", point_id_set)
        (previous_pair === nothing || previous_pair < pair) ||
            throw(ArgumentError("raw external candidate pairs must be unique and sorted"))
        previous_pair = pair
        native_i = _initial_policy_int(row["native_i"], "raw candidate.native_i";
                                       minimum=1, maximum=INITIAL_POLICY_MAX_POINTS)
        native_j = _initial_policy_int(row["native_j"], "raw candidate.native_j";
                                       minimum=1, maximum=INITIAL_POLICY_MAX_POINTS)
        native_i == external_to_native[pair[1]] && native_j == external_to_native[pair[2]] ||
            throw(ArgumentError("raw candidate native rows do not map to its external IDs"))
        state = _initial_policy_string(row["classification"], "raw candidate.classification"; max_bytes=32)
        state in INITIAL_POLICY_STATES || throw(ArgumentError("invalid raw candidate state"))
        _initial_policy_string(row["reason"], "raw candidate.reason"; max_bytes=INITIAL_POLICY_MAX_REASON_BYTES)
        push!(decisions, (pair[1], pair[2], state))
        decision_by_pair[pair] = state
    end

    counts = _initial_policy_table(root["counts"],
        ("raw_euclidean_candidate_pairs", "initial_material_bond_pairs",
         "geometry_excluded_pairs", "unknown_pairs"), "counts")
    state_counts = Dict(state => count(row -> row[3] == state, decisions) for state in INITIAL_POLICY_STATES)
    expected_counts = (
        raw_euclidean_candidate_pairs=length(decisions),
        initial_material_bond_pairs=state_counts["IN_MATERIAL"],
        geometry_excluded_pairs=state_counts["OUT_OR_VOID"],
        unknown_pairs=state_counts["UNKNOWN"],
    )
    for (name, expected) in pairs(expected_counts)
        actual = _initial_policy_int(counts[string(name)], "counts.$(name)";
                                     maximum=INITIAL_POLICY_MAX_PAIRS)
        actual == expected || throw(ArgumentError("counts.$(name) disagrees with raw candidate rows"))
    end

    hashes = _initial_policy_table(root["hashes"], ("pair_map_sha256", "mask_sha256"), "hashes")
    pair_map_sha = _initial_policy_sha(hashes["pair_map_sha256"], "hashes.pair_map_sha256")
    mask_sha = _initial_policy_sha(hashes["mask_sha256"], "hashes.mask_sha256")
    _initial_policy_pair_map_sha(decisions) == pair_map_sha ||
        throw(ArgumentError("pair_map_sha256 does not match the complete state ledger"))

    mask = _initial_policy_table(root["mask"],
        ("included_pair_count", "included_pairs", "geometry_excluded_pair_count",
         "geometry_excluded_pairs", "unknown_pair_count", "unknown_pairs"), "mask")
    groups = (
        ("included_pairs", "included_pair_count", "IN_MATERIAL"),
        ("geometry_excluded_pairs", "geometry_excluded_pair_count", "OUT_OR_VOID"),
        ("unknown_pairs", "unknown_pair_count", "UNKNOWN"),
    )
    expected_mask_pairs = Dict{String,Vector{Tuple{Int,Int}}}()
    for (field, count_field, state) in groups
        expected = [(row[1], row[2]) for row in decisions if row[3] == state]
        raw_rows = _initial_policy_array(mask[field], "mask.$(field)"; max_items=INITIAL_POLICY_MAX_PAIRS)
        actual = Tuple{Int,Int}[]
        for (index, raw_row) in enumerate(raw_rows)
            row = _initial_policy_table(raw_row, ("point_i", "point_j", "native_i", "native_j"),
                                        "mask.$(field)[$(index)]")
            pair = _initial_policy_pair(row, "mask.$(field) pair", point_id_set)
            native_i = _initial_policy_int(row["native_i"], "mask.$(field).native_i";
                                           minimum=1, maximum=INITIAL_POLICY_MAX_POINTS)
            native_j = _initial_policy_int(row["native_j"], "mask.$(field).native_j";
                                           minimum=1, maximum=INITIAL_POLICY_MAX_POINTS)
            native_i == external_to_native[pair[1]] && native_j == external_to_native[pair[2]] ||
                throw(ArgumentError("mask native rows do not map to canonical external IDs"))
            push!(actual, pair)
        end
        actual == expected ||
            throw(ArgumentError("mask.$(field) is not the exact state partition of the full ledger"))
        _initial_policy_int(mask[count_field], "mask.$(count_field)";
                            maximum=INITIAL_POLICY_MAX_PAIRS) == length(expected) ||
            throw(ArgumentError("mask.$(count_field) disagrees with pair rows"))
        expected_mask_pairs[field] = expected
    end
    _initial_policy_mask_sha(expected_mask_pairs["included_pairs"]) == mask_sha ||
        throw(ArgumentError("mask_sha256 does not match the included external-ID pairs"))

    degrees = Dict(point_id => 0 for point_id in point_ids)
    for (left, right) in expected_mask_pairs["included_pairs"]
        degrees[left] += 1
        degrees[right] += 1
    end
    zero_degree_ids = [point_id for point_id in point_ids if degrees[point_id] == 0]
    gate = _initial_policy_table(root["native_gate"],
        ("status", "raw_candidate_set_reproduced", "native_graph_readback",
         "initialization_authorized", "zero_degree_external_point_ids"), "native_gate")
    gate["status"] == "UNVERIFIED" ||
        throw(ArgumentError("audit contract cannot claim native verification"))
    _initial_policy_bool(gate["raw_candidate_set_reproduced"], "native_gate.raw_candidate_set_reproduced") == false ||
        throw(ArgumentError("static transport cannot claim raw native candidates were reproduced"))
    gate["native_graph_readback"] == "NOT_RUN" ||
        throw(ArgumentError("static transport cannot claim native graph readback"))
    _initial_policy_bool(gate["initialization_authorized"], "native_gate.initialization_authorized") == false ||
        throw(ArgumentError("static transport cannot authorize initialization"))
    reported_zero = _initial_policy_array(gate["zero_degree_external_point_ids"],
        "native_gate.zero_degree_external_point_ids"; max_items=INITIAL_POLICY_MAX_POINTS)
    parsed_zero = [_initial_policy_int(id, "native gate zero-degree ID") for id in reported_zero]
    length(parsed_zero) == length(Set(parsed_zero)) && parsed_zero == zero_degree_ids ||
        throw(ArgumentError("native gate zero-degree IDs disagree with included mask degrees"))

    expected_eligibility = if state_counts["UNKNOWN"] > 0
        "BLOCKED_UNKNOWN"
    elseif isempty(decisions)
        "BLOCKED_EMPTY_CANDIDATE_SET"
    elseif isempty(expected_mask_pairs["included_pairs"])
        "BLOCKED_EMPTY_MASK"
    elseif !isempty(zero_degree_ids)
        "BLOCKED_ZERO_DEGREE"
    else
        "CANDIDATE_MASK_POSITIVE_DEGREE_NATIVE_UNVERIFIED"
    end
    eligibility == expected_eligibility ||
        throw(ArgumentError("native_mask_eligibility disagrees with candidate, unknown, and degree states"))

    evidence = _initial_policy_table(root["evidence"],
        ("classifier_result_path", "classifier_result_sha256", "raw_hit_interval_sidecar_path",
         "raw_hit_interval_sidecar_sha256", "sidecar_test_nodeid", "sidecar_call_index",
         "visibility_algorithm"), "evidence")
    _initial_policy_string(evidence["classifier_result_path"], "evidence.classifier_result_path")
    _initial_policy_sha(evidence["classifier_result_sha256"], "evidence.classifier_result_sha256")
    _initial_policy_string(evidence["raw_hit_interval_sidecar_path"], "evidence.raw_hit_interval_sidecar_path")
    _initial_policy_sha(evidence["raw_hit_interval_sidecar_sha256"], "evidence.raw_hit_interval_sidecar_sha256")
    _initial_policy_string(evidence["sidecar_test_nodeid"], "evidence.sidecar_test_nodeid"; max_bytes=4096)
    _initial_policy_int(evidence["sidecar_call_index"], "evidence.sidecar_call_index";
                        maximum=Int64(typemax(Int32)))
    _initial_policy_string(evidence["visibility_algorithm"], "evidence.visibility_algorithm"; max_bytes=512)

    physics = _initial_policy_table(root["initial_physics_events"],
        ("evidence_status", "planned_physical_broken_bond_pairs", "geometry_exclusions_are_physics_events"),
        "initial_physics_events")
    physics["evidence_status"] == "PLANNED_NOT_OBSERVED_NATIVE_NOT_RUN" ||
        throw(ArgumentError("initial physics event evidence must remain planned"))
    _initial_policy_int(physics["planned_physical_broken_bond_pairs"], "planned physical broken pairs") == 0 ||
        throw(ArgumentError("static transport cannot report observed physical broken bonds"))
    _initial_policy_bool(physics["geometry_exclusions_are_physics_events"],
                         "geometry_exclusions_are_physics_events") == false ||
        throw(ArgumentError("geometry exclusions cannot be physical events"))

    return (
        root=root,
        transport=transport,
        source=source,
        coordinates=coordinates,
        transform=transform,
        points=points,
        point_rows=point_rows,
        point_ids=point_ids,
        external_to_native=external_to_native,
        native_to_external=native_to_external,
        positions_mm=positions_mm,
        horizon_m=horizon_m,
        decisions=decisions,
        decision_by_pair=decision_by_pair,
        expected_counts=expected_counts,
        expected_mask_pairs=expected_mask_pairs,
        pair_map_sha256=pair_map_sha,
        mask_sha256=mask_sha,
        evidence=evidence,
        zero_degree_ids=zero_degree_ids,
        eligibility=eligibility,
    )
end

function _initial_policy_artifact_specs(model)
    source = model.source
    points = model.points
    evidence = model.evidence
    return [
        ("source CAD", _initial_policy_string(source["source_cad_path"], "source CAD path"), INITIAL_POLICY_CAD_MAX_BYTES),
        ("ordered point CSV", _initial_policy_string(points["payload_path"], "point CSV path"), INITIAL_POLICY_CSV_MAX_BYTES),
        ("classifier result JSON", _initial_policy_string(evidence["classifier_result_path"], "classifier result path"), INITIAL_POLICY_RESULT_MAX_BYTES),
        ("classifier source", _initial_policy_string(source["classifier_source_path"], "classifier source path"), INITIAL_POLICY_SOURCE_MAX_BYTES),
        ("raw-hit interval sidecar", _initial_policy_string(evidence["raw_hit_interval_sidecar_path"], "sidecar path"), INITIAL_POLICY_SIDECAR_MAX_BYTES),
    ]
end

function _initial_policy_csv_row(bytes::Vector{UInt8}, index::Int)
    fields = String[]
    n = length(bytes)
    while true
        field = IOBuffer()
        if index <= n && bytes[index] == UInt8('"')
            index += 1
            closed = false
            while index <= n
                byte = bytes[index]
                if byte == UInt8('"')
                    if index < n && bytes[index + 1] == UInt8('"')
                        write(field, UInt8('"'))
                        index += 2
                    else
                        index += 1
                        closed = true
                        break
                    end
                else
                    write(field, byte)
                    index += 1
                end
            end
            closed || throw(ArgumentError("ordered point CSV has an unterminated quoted field"))
            if index <= n && bytes[index] != UInt8(',') &&
               bytes[index] != UInt8('\n') && bytes[index] != UInt8('\r')
                throw(ArgumentError("ordered point CSV has characters after a quoted field"))
            end
        else
            while index <= n && bytes[index] != UInt8(',') &&
                  bytes[index] != UInt8('\n') && bytes[index] != UInt8('\r')
                bytes[index] == UInt8('"') &&
                    throw(ArgumentError("ordered point CSV has a quote inside an unquoted field"))
                write(field, bytes[index])
                index += 1
            end
        end
        push!(fields, String(take!(field)))
        if index > n
            return fields, index, true
        elseif bytes[index] == UInt8(',')
            index += 1
            continue
        else
            if bytes[index] == UInt8('\r') && index < n && bytes[index + 1] == UInt8('\n')
                index += 2
            else
                index += 1
            end
            return fields, index, false
        end
    end
end

function _initial_policy_parse_points_csv(bytes::Vector{UInt8})
    csv_text = String(copy(bytes))
    isvalid(csv_text) || throw(ArgumentError("ordered point CSV is not valid UTF-8"))
    isempty(bytes) && throw(ArgumentError("ordered point CSV is empty"))
    rows = Vector{Tuple{Int,NTuple{3,Float64}}}()
    index = 1
    header, index, at_eof = _initial_policy_csv_row(bytes, index)
    header == ["point_id", "x_mm", "y_mm", "z_mm"] ||
        throw(ArgumentError("ordered point CSV header must be point_id,x_mm,y_mm,z_mm"))
    at_eof && throw(ArgumentError("ordered point CSV must contain at least one point"))
    seen = Set{Int}()
    while index <= length(bytes)
        fields, next_index, ended_at_eof = _initial_policy_csv_row(bytes, index)
        length(fields) == 4 || throw(ArgumentError("ordered point CSV row must contain four columns"))
        id_text = fields[1]
        canonical_id = !isempty(id_text) &&
            all(char -> '0' <= char <= '9', id_text) &&
            (id_text == "0" || first(id_text) != '0')
        canonical_id || throw(ArgumentError("ordered point CSV ID is not canonical non-negative decimal"))
        point_id = tryparse(Int, id_text)
        point_id === nothing && throw(ArgumentError("ordered point CSV ID is outside native integer range"))
        point_id in seen && throw(ArgumentError("ordered point CSV contains a duplicate ID"))
        push!(seen, point_id)
        coordinates = ntuple(axis -> begin
            value = tryparse(Float64, strip(fields[axis + 1]))
            value !== nothing && isfinite(value) ||
                throw(ArgumentError("ordered point CSV contains an invalid or non-finite coordinate"))
            value
        end, 3)
        push!(rows, (point_id, coordinates))
        length(rows) <= INITIAL_POLICY_MAX_POINTS ||
            throw(ArgumentError("ordered point CSV exceeds the point cap"))
        index = next_index
        ended_at_eof && break
    end
    isempty(rows) && throw(ArgumentError("ordered point CSV must contain at least one point"))
    return rows
end

function _initial_policy_match_points_csv(bytes::Vector{UInt8}, point_ids::Vector{Int},
                                          positions_mm::Vector{NTuple{3,Float64}})
    csv_rows = _initial_policy_parse_points_csv(bytes)
    length(csv_rows) == length(point_ids) == length(positions_mm) ||
        throw(ArgumentError("ordered point CSV row count differs from TOML point rows"))
    for (native_index, (point_id, position_mm)) in enumerate(csv_rows)
        point_id == point_ids[native_index] &&
            all(isequal(position_mm[axis], positions_mm[native_index][axis]) for axis in 1:3) ||
            throw(ArgumentError("ordered point CSV IDs/coordinates differ from TOML point rows"))
    end
    return nothing
end

function _verify_initial_policy_artifact_pins(model; reserved_paths::Vector{String}=String[])
    specs = _initial_policy_artifact_specs(model)
    captured, resolved = _read_initial_policy_bundle(specs;
        max_total_bytes=INITIAL_POLICY_ARTIFACT_TOTAL_MAX_BYTES)
    any(path -> path in reserved_paths, values(resolved)) &&
        throw(ArgumentError("source evidence paths must be distinct from policy JSON/TOML paths"))
    pin_by_label = Dict(
        "source CAD" => model.source["source_cad_sha256"],
        "ordered point CSV" => model.points["payload_sha256"],
        "classifier result JSON" => model.evidence["classifier_result_sha256"],
        "classifier source" => model.source["classifier_source_sha256"],
        "raw-hit interval sidecar" => model.evidence["raw_hit_interval_sidecar_sha256"],
    )
    for (label, _, _) in specs
        _initial_policy_sha256(captured[label]) == pin_by_label[label] ||
            throw(ArgumentError("actual $(label) bytes do not match their TOML SHA pin"))
    end
    _initial_policy_match_points_csv(captured["ordered point CSV"], model.point_ids,
                                    model.positions_mm)
    return resolved
end

function InitialBondPolicy(audit_json_path::AbstractString,
                           native_toml_path::AbstractString;
                           expected_json_sha256::AbstractString,
                           expected_toml_sha256::AbstractString)
    expected_json = _initial_policy_sha(String(expected_json_sha256), "expected_json_sha256")
    expected_toml = _initial_policy_sha(String(expected_toml_sha256), "expected_toml_sha256")
    audit_path = String(audit_json_path)
    toml_path = String(native_toml_path)
    captured, resolved = _read_initial_policy_bundle([
        ("audit JSON", audit_path, INITIAL_POLICY_AUDIT_MAX_BYTES),
        ("native TOML", toml_path, INITIAL_POLICY_TOML_MAX_BYTES),
    ]; max_total_bytes=INITIAL_POLICY_TOTAL_INPUT_MAX_BYTES)
    actual_json = _initial_policy_sha256(captured["audit JSON"])
    actual_toml = _initial_policy_sha256(captured["native TOML"])
    actual_json == expected_json ||
        throw(ArgumentError("audit JSON bytes do not match the externally pinned receipt SHA"))
    actual_toml == expected_toml ||
        throw(ArgumentError("native TOML bytes do not match the externally pinned receipt SHA"))
    resolved_audit = resolved["audit JSON"]
    resolved_toml = resolved["native TOML"]
    toml_text = String(copy(captured["native TOML"]))
    model = try
        TOML.parse(toml_text)
    catch
        throw(ArgumentError("native mask transport is not valid TOML"))
    end
    parsed = _initial_policy_validate_transport(model, actual_json, resolved_audit)
    _verify_initial_policy_artifact_pins(parsed;
        reserved_paths=[resolved_audit, resolved_toml])
    return InitialBondPolicy(resolved_audit, resolved_toml, toml_text, actual_json, actual_toml,
                             _INITIAL_BOND_POLICY_TOKEN)
end

function _initial_policy_current_model(policy::InitialBondPolicy)
    captured, resolved = _read_initial_policy_bundle([
        ("audit JSON", policy.audit_json_path, INITIAL_POLICY_AUDIT_MAX_BYTES),
        ("native TOML", policy.native_toml_path, INITIAL_POLICY_TOML_MAX_BYTES),
    ]; max_total_bytes=INITIAL_POLICY_TOTAL_INPUT_MAX_BYTES)
    actual_json = _initial_policy_sha256(captured["audit JSON"])
    actual_toml = _initial_policy_sha256(captured["native TOML"])
    actual_json == policy.audit_json_sha256 ||
        throw(ArgumentError("audit JSON changed after InitialBondPolicy construction"))
    actual_toml == policy.native_toml_sha256 ||
        throw(ArgumentError("native TOML changed after InitialBondPolicy construction"))
    resolved["audit JSON"] == policy.audit_json_path &&
        resolved["native TOML"] == policy.native_toml_path ||
        throw(ArgumentError("policy artifact paths changed after construction"))
    captured["native TOML"] == collect(codeunits(policy.native_toml_bytes)) ||
        throw(ArgumentError("native TOML bytes differ from immutable InitialBondPolicy bytes"))
    model = try
        TOML.parse(policy.native_toml_bytes)
    catch
        throw(ArgumentError("immutable native mask transport is not valid TOML"))
    end
    parsed = _initial_policy_validate_transport(model, actual_json, policy.audit_json_path)
    _verify_initial_policy_artifact_pins(parsed;
        reserved_paths=[policy.audit_json_path, policy.native_toml_path])
    return parsed
end

function _validate_initial_policy_scope(body::AbstractBody, pd,
                                        chunk_id::Int, policy::InitialBondPolicy)
    model = _initial_policy_current_model(policy)
    model.eligibility == "CANDIDATE_MASK_POSITIVE_DEGREE_NATIVE_UNVERIFIED" ||
        throw(ArgumentError("blocked or unverified-empty mask cannot construct a native body"))
    model.expected_counts.unknown_pairs == 0 ||
        throw(ArgumentError("UNKNOWN candidate pairs block the complete native asset"))
    isempty(model.decisions) && throw(ArgumentError("empty Euclidean candidate set is blocked"))
    isempty(model.expected_mask_pairs["included_pairs"]) &&
        throw(ArgumentError("empty included mask is blocked"))
    isempty(model.zero_degree_ids) || throw(ArgumentError("zero-degree native points are blocked"))
    body.mat isa OSBMaterial ||
        throw(ArgumentError("initial bond policy v1 supports only OSBMaterial"))
    pd.n_chunks == 1 && length(pd.decomp) == 1 && chunk_id == 1 ||
        throw(ArgumentError("initial bond policy v1 requires exactly one chunk"))
    pd.decomp[1] == (1:body.n_points) ||
        throw(ArgumentError("one native chunk must contain every body row in exact order"))

    point_count = model.points["count"]
    body.n_points == point_count && size(body.position) == (3, point_count) ||
        throw(ArgumentError("native Body point count/shape differs from the pinned point rows"))
    scale = model.coordinates["m_per_input_unit"]
    matrix = model.transform
    for row in model.point_rows
        native_index = _initial_policy_int(row["native_index"], "point native_index";
                                           minimum=1, maximum=INITIAL_POLICY_MAX_POINTS)
        xyz_mm = row["position_mm"]
        for native_axis in 1:3
            expected_m = matrix[native_axis, 4] +
                         sum(matrix[native_axis, reference_axis] *
                             (Float64(xyz_mm[reference_axis]) * scale)
                             for reference_axis in 1:3)
            actual_m = body.position[native_axis, native_index]
            isfinite(actual_m) &&
                isapprox(actual_m, expected_m; atol=INITIAL_POLICY_POSITION_ATOL_M, rtol=0.0) ||
                throw(ArgumentError("Body.position does not match the pinned point row/frame transform"))
        end
    end

    horizon = model.horizon_m
    for native_index in 1:body.n_points
        delta = get_point_param(body, :δ, native_index)
        isfinite(delta) && delta == horizon ||
            throw(ArgumentError("initial bond policy v1 requires the exact pinned uniform horizon"))
    end
    return model
end

function _reject_initial_policy_in_multibody(ms::AbstractMultibodySetup,
                                             handler_name::AbstractString)
    for body in each_body(ms)
        if initial_bond_policy(body) !== nothing
            throw(ArgumentError("initial bond policy v1 requires one body; $(handler_name) multibody setup is rejected"))
        end
    end
    return nothing
end

function _reject_initial_policy_in_mpi_body(body::AbstractBody)
    initial_bond_policy(body) === nothing ||
        throw(ArgumentError("initial bond policy v1 is supported only by the Threads CPU handler"))
    return nothing
end

function _initial_policy_admission_receipt(policy::InitialBondPolicy, model,
                                           raw_directed::Int, included_directed::Int)
    return InitialBondAdmissionReceipt(
        "pgfrag.native-bond-admission.v1",
        "RAW_CANDIDATES_REPRODUCED_MASK_APPLIED_BEFORE_PRECRACK",
        policy.audit_json_sha256,
        policy.native_toml_sha256,
        model.source["source_cad_sha256"],
        model.points["payload_sha256"],
        model.evidence["classifier_result_sha256"],
        model.source["classifier_source_sha256"],
        model.evidence["raw_hit_interval_sidecar_sha256"],
        model.pair_map_sha256,
        model.mask_sha256,
        length(model.decisions),
        raw_directed,
        length(model.expected_mask_pairs["included_pairs"]),
        included_directed,
        length(model.expected_mask_pairs["geometry_excluded_pairs"]),
        0,
        true,
        true,
        false,
        "NOT_RUN",
        false,
    )
end

function _filter_initial_policy_candidates!(bonds::AbstractVector,
                                            n_neighbors::Vector{Int},
                                            loc_points::AbstractVector{Int},
                                            body::AbstractBody,
                                            model,
                                            policy::InitialBondPolicy)
    length(loc_points) == body.n_points && length(n_neighbors) == body.n_points ||
        throw(ArgumentError("native raw search is not a complete one-chunk point list"))
    native_to_external = model.native_to_external
    decision_by_pair = model.decision_by_pair
    observed_directions = Dict{Tuple{Int,Int},Int}()
    offset = 0
    raw_directed_count = 0
    for (local_index, source_native) in enumerate(loc_points)
        1 <= source_native <= body.n_points ||
            throw(ArgumentError("raw native candidate source is outside the Body rows"))
        count = n_neighbors[local_index]
        count >= 0 || throw(ArgumentError("raw native neighbor count is negative"))
        last_index = offset + count
        last_index <= length(bonds) ||
            throw(ArgumentError("raw native neighbor counts exceed the bond ledger"))
        source_external = native_to_external[source_native]
        for bond_index in (offset + 1):last_index
            target_native = bonds[bond_index].neighbor
            1 <= target_native <= body.n_points && target_native != source_native ||
                throw(ArgumentError("raw native candidate has an invalid endpoint"))
            target_external = native_to_external[target_native]
            left, right = source_external < target_external ?
                          (source_external, target_external) : (target_external, source_external)
            haskey(decision_by_pair, (left, right)) ||
                throw(ArgumentError("native Euclidean search produced an extra external candidate"))
            direction = (source_external, target_external)
            observed_directions[direction] = get(observed_directions, direction, 0) + 1
            raw_directed_count += 1
        end
        offset = last_index
    end
    offset == length(bonds) ||
        throw(ArgumentError("raw native bond vector contains uncounted directed entries"))

    expected_directions = Set{Tuple{Int,Int}}()
    for (left, right, _) in model.decisions
        push!(expected_directions, (left, right))
        push!(expected_directions, (right, left))
    end
    Set(keys(observed_directions)) == expected_directions ||
        throw(ArgumentError("native Euclidean candidate set is incomplete or lacks reciprocal directions"))
    all(count == 1 for count in values(observed_directions)) ||
        throw(ArgumentError("native Euclidean search contains duplicate directed candidates"))
    raw_directed_count == 2 * length(model.decisions) ||
        throw(ArgumentError("native directed candidate count is not exactly two per canonical pair"))

    kept_bonds = Bond[]
    sizehint!(kept_bonds, raw_directed_count)
    kept_neighbors = zeros(Int, length(n_neighbors))
    included_directions = Set{Tuple{Int,Int}}()
    offset = 0
    included_directed_count = 0
    for (local_index, source_native) in enumerate(loc_points)
        count = n_neighbors[local_index]
        last_index = offset + count
        source_external = native_to_external[source_native]
        for bond_index in (offset + 1):last_index
            bond = bonds[bond_index]
            target_external = native_to_external[bond.neighbor]
            left, right = source_external < target_external ?
                          (source_external, target_external) : (target_external, source_external)
            if decision_by_pair[(left, right)] == "IN_MATERIAL"
                direction = (source_external, target_external)
                direction in included_directions &&
                    throw(ArgumentError("included native graph contains a duplicate directed candidate"))
                push!(included_directions, direction)
                push!(kept_bonds, bond)
                kept_neighbors[local_index] += 1
                included_directed_count += 1
            end
        end
        offset = last_index
    end
    all(count > 0 for count in kept_neighbors) ||
        throw(ArgumentError("geometry mask would leave a native point with zero included degree"))
    expected_included_directions = Set{Tuple{Int,Int}}()
    for (left, right) in model.expected_mask_pairs["included_pairs"]
        push!(expected_included_directions, (left, right))
        push!(expected_included_directions, (right, left))
    end
    included_directions == expected_included_directions ||
        throw(ArgumentError("filtered directed native graph does not equal the complete reciprocal included mask"))
    included_directed_count == 2 * length(model.expected_mask_pairs["included_pairs"]) ||
        throw(ArgumentError("filtered directed graph does not match included canonical mask"))

    empty!(bonds)
    append!(bonds, kept_bonds)
    n_neighbors .= kept_neighbors
    receipt = _initial_policy_admission_receipt(policy, model, raw_directed_count,
                                                included_directed_count)
    return receipt
end
