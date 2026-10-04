using Peridynamics
using LinearAlgebra
using SHA
using TOML

function require(condition, message)
    condition || error(message)
end

function read_fixture(points_path::String, metadata_path::String)
    metadata = TOML.parsefile(metadata_path)
    actual_hash = bytes2hex(SHA.sha256(read(points_path)))
    require(actual_hash == lowercase(String(metadata["point_data_sha256"])),
            "point CSV hash does not match frozen metadata")
    lines = readlines(points_path)
    require(!isempty(lines) && strip(first(lines)) == "point_id,x_m,y_m,z_m,volume_m3",
            "unexpected point CSV header")
    n = length(lines) - 1
    require(n == Int(metadata["point_count"]), "point count disagrees with metadata")
    point_ids = Vector{Int}(undef, n)
    position = Matrix{Float64}(undef, 3, n)
    volume = Vector{Float64}(undef, n)
    seen = Set{Int}()
    for (i, line) in enumerate(lines[2:end])
        fields = split(strip(line), ',')
        require(length(fields) == 5, "invalid point CSV row $i")
        point_ids[i] = parse(Int, fields[1])
        require(point_ids[i] > 0 && !(point_ids[i] in seen), "invalid or duplicate point id")
        push!(seen, point_ids[i])
        position[:, i] .= parse.(Float64, fields[2:4])
        volume[i] = parse(Float64, fields[5])
        require(isfinite(volume[i]) && volume[i] > 0 && all(isfinite, position[:, i]),
                "invalid point geometry")
    end
    return metadata, point_ids, position, volume
end

function write_bond_csv(path::String, rows)
    open(path, "w") do io
        println(io, "state_id,eta_over_native_epsilon_c,eta,native_epsilon_c,source_point_id,neighbor_point_id,canonical_i_id,canonical_j_id,reference_length_m,current_length_m,actual_stretch,fail_permit,active_before,active_after,transition")
        for row in rows
            println(io, join(string.(row), ','))
        end
    end
end

function write_pair_csv(path::String, rows)
    open(path, "w") do io
        println(io, "state_id,eta_over_native_epsilon_c,eta,native_epsilon_c,canonical_i_id,canonical_j_id,reference_length_m,current_length_m,actual_stretch,fail_permit,active_before,active_after,transition")
        for row in rows
            println(io, join(string.(row), ','))
        end
    end
end

function build_native_handler(reference_position, volume, material_cfg,
                              geometry_cfg, safety_factor::Float64)
    body = Body(OSBMaterial(), reference_position, volume)
    material!(body; horizon=Float64(geometry_cfg["horizon_m"]),
              rho=Float64(material_cfg["rho_kg_m3"]), E=Float64(material_cfg["E_pa"]),
              nu=Float64(material_cfg["nu"]), Gc=Float64(material_cfg["Gc_j_m2"]))
    solver = VelocityVerlet(steps=1, safety_factor=safety_factor)
    dh = Peridynamics.threads_data_handler(body, solver, Threads.nthreads())
    require(length(dh.chunks) == 1, "native diagnostic requires one chunk")
    native_dt = Peridynamics.calc_stable_timestep(dh, safety_factor)
    require(isfinite(native_dt) && native_dt > 0, "native stability reference is invalid")
    solver.Δt = native_dt
    return solver, dh
end

function main(args)
    require(length(args) == 5,
            "expected arguments: points.csv metadata.toml common.toml case.toml output_dir")
    points_path, metadata_path, common_path, case_path, output_dir = args
    require(!ispath(output_dir), "output directory already exists; refusing overwrite")
    case = TOML.parsefile(case_path)
    common = TOML.parsefile(common_path)
    fixture, point_ids, reference_position, volume = read_fixture(points_path, metadata_path)
    require(Threads.nthreads() == 1, "diagnostic requires exactly one Julia thread")
    require(length(point_ids) <= Int(common["limits"]["max_points"]), "point limit exceeded")
    require(length(point_ids) == 125, "unexpected diagnostic point count")

    material = common["material"]
    geometry = common["geometry"]
    solver_cfg = common["solver"]
    expected = case["expected_graph"]
    expected_canonical = Int(expected["canonical_bonds"])
    state_specs = case["states"]
    require(length(state_specs) == 3, "exactly three affine states are frozen")
    require([Float64(s["epsilon_c_multiplier"]) for s in state_specs] == [0.0, 0.5, 1.5],
            "affine threshold multipliers changed")

    rows = Vector{Vector{Any}}()
    canonical_rows = Vector{Vector{Any}}()
    state_summaries = Vector{NamedTuple}()
    epsilon_c = NaN
    native_dt = NaN
    for (state_index, spec) in enumerate(state_specs)
        solver, dh = build_native_handler(reference_position, volume, material, geometry,
                                          Float64(solver_cfg["safety_factor"]))
        chunk = dh.chunks[1]
        system = chunk.system
        storage = chunk.storage
        require(size(system.position) == size(storage.position), "native position shape mismatch")
        require(size(storage.position, 2) == length(point_ids),
                "one-chunk native point order does not match frozen CSV row order")
        require(length(system.bonds) == Int(expected["directed_bonds"]),
                "native directed graph differs from frozen run01 reference")
        native_epsilon_by_point = [Peridynamics.get_params(chunk, i).εc
                                   for i in eachindex(point_ids)]
        state_epsilon_c = first(native_epsilon_by_point)
        require(all(isapprox(value, state_epsilon_c; rtol=1e-12, atol=0.0)
                    for value in native_epsilon_by_point),
                "uniform-material fixture has inconsistent native epsilon_c")
        if state_index == 1
            epsilon_c = state_epsilon_c
            native_dt = solver.Δt
        else
            require(isapprox(state_epsilon_c, epsilon_c; rtol=1e-12, atol=0.0),
                    "native epsilon_c changed between pristine state branches")
            require(isapprox(solver.Δt, native_dt; rtol=1e-12, atol=0.0),
                    "native reference dt changed between identical pristine state branches")
        end
        multiplier = Float64(spec["epsilon_c_multiplier"])
        eta = multiplier * epsilon_c
        pristine = copy(system.position)
        storage.position .= pristine
        @views storage.position[3, :] .= (1 + eta) .* pristine[3, :]
        active_before = copy(storage.bond_active)
        Peridynamics.calc_force_density!(dh, 0.0, solver.Δt)
        active_after = copy(storage.bond_active)
        require(length(active_before) == length(active_after) == length(system.bonds),
                "native bond state shape changed during force evaluation")
        pair_rows = Dict{Tuple{Int,Int},Vector{NamedTuple}}()
        pair_changed = Set{Tuple{Int,Int}}()
        directed_changes = 0
        for i in Peridynamics.each_point_idx(system)
            eps_i = Peridynamics.get_params(chunk, i).εc
            for bond_id in Peridynamics.each_bond_idx(system, i)
                bond = system.bonds[bond_id]
                j = bond.neighbor
                source_id = point_ids[i]
                neighbor_id = point_ids[j]
                pair = minmax(source_id, neighbor_id)
                before = active_before[bond_id]
                after = active_after[bond_id]
                transition = before != after
                directed_changes += transition
                transition && push!(pair_changed, pair)
                reference_length = bond.length
                current_length = storage.bond_length[bond_id]
                stretch = (current_length - reference_length) / reference_length
                require(isfinite(stretch) && isfinite(eps_i), "non-finite native bond measure")
                measure = (before=before, after=after,
                           reference_length=reference_length,
                           current_length=current_length, stretch=stretch,
                           epsilon_c=eps_i, fail_permit=bond.fail_permit)
                push!(get!(pair_rows, pair, NamedTuple[]), measure)
                push!(rows, Any[String(spec["state_id"]), multiplier, eta, eps_i,
                                 source_id, neighbor_id, pair[1], pair[2], reference_length,
                                 current_length, stretch, bond.fail_permit, before, after, transition])
            end
        end
        require(length(pair_rows) == expected_canonical, "canonical graph size changed")
        reciprocal_ok = true
        for (pair, directions) in pair_rows
            if length(directions) != 2
                reciprocal_ok = false
                continue
            end
            forward, reverse = directions
            pair_consistent = forward.before == reverse.before &&
                              forward.after == reverse.after &&
                              forward.fail_permit == reverse.fail_permit &&
                              isapprox(forward.reference_length, reverse.reference_length; rtol=1e-12, atol=0.0) &&
                              isapprox(forward.current_length, reverse.current_length; rtol=1e-12, atol=0.0) &&
                              isapprox(forward.stretch, reverse.stretch; rtol=1e-12, atol=1e-15) &&
                              isapprox(forward.epsilon_c, reverse.epsilon_c; rtol=1e-12, atol=0.0)
            reciprocal_ok &= pair_consistent
            transition = forward.before != forward.after
            push!(canonical_rows, Any[String(spec["state_id"]), multiplier, eta,
                forward.epsilon_c, pair[1], pair[2], forward.reference_length,
                forward.current_length, forward.stretch, forward.fail_permit,
                forward.before, forward.after, transition])
        end
        push!(state_summaries, (state_id=String(spec["state_id"]),
                                eta_over_epsilon_c=multiplier,
                                eta=eta,
                                directed_state_rows=length(system.bonds),
                                canonical_pairs=length(pair_rows),
                                directed_transition_count=directed_changes,
                                canonical_transition_count=length(pair_changed),
                                reciprocal_status=reciprocal_ok ? "PASS" : "FAIL"))
    end

    low_ok = state_summaries[1].canonical_transition_count == 0 &&
             state_summaries[2].canonical_transition_count == 0
    high_ok = state_summaries[3].canonical_transition_count > 0
    directed_bond_counts = [state.directed_state_rows for state in state_summaries]
    all_graph_counts_match = length(directed_bond_counts) == 3 &&
        all(count -> count == Int(expected["directed_bonds"]), directed_bond_counts)
    require(all_graph_counts_match, "three pristine state branches do not share the expected native graph")
    reciprocity_ok = all(s.reciprocal_status == "PASS" for s in state_summaries)
    gate_status = low_ok && high_ok && reciprocity_ok ? "PASS_NATIVE_BRANCH_DIAGNOSTIC_ONLY" : "FAIL"
    mkpath(output_dir)
    write_bond_csv(joinpath(output_dir, "affine_native_bond_states.csv"), rows)
    write_pair_csv(joinpath(output_dir, "affine_canonical_pair_states.csv"), canonical_rows)
    open(joinpath(output_dir, "summary.toml"), "w") do io
        TOML.print(io, Dict{String,Any}(
            "schema" => "pgfrag.pd.affine_critical_stretch_result.v1",
            "status" => gate_status,
            "experiment_kind" => "native constitutive branch diagnostic; no time integration",
            "physical_qualification" => "NONE; UNVALIDATED_PROTOTYPE parameters",
            "input_id" => String(fixture["input_id"]),
            "point_csv_sha256" => bytes2hex(SHA.sha256(read(points_path))),
            "point_count" => length(point_ids),
            "native_directed_bonds" => state_summaries[1].directed_state_rows,
            "native_canonical_bonds" => expected_canonical,
            "native_epsilon_c" => epsilon_c,
            "native_stable_dt_reference_s_not_integrated" => native_dt,
            "native_time_integration" => false,
            "manual_bond_activity_edits" => false,
            "below_threshold_controls_clear" => low_ok,
            "above_threshold_native_event_observed" => high_ok,
            "reciprocal_pair_states_match" => reciprocity_ok,
            "energy_ledger" => "UNKNOWN_NOT_A_GATE",
            "contact" => "NOT_RUN",
            "free_fracture_or_fragment_claim" => "NONE"))
        for (i, state) in enumerate(state_summaries)
            println(io)
            println(io, "[state_$(i)]")
            for (key, value) in pairs(state)
                println(io, key, " = ", value isa String ? repr(value) : value)
            end
        end
    end
    println("AFFINE_NATIVE_DIAGNOSTIC_STATUS=", gate_status)
    gate_status == "PASS_NATIVE_BRANCH_DIAGNOSTIC_ONLY" || exit(2)
end

main(ARGS)
