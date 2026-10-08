using Test
using LinearAlgebra
using Peridynamics

const PB = Peridynamics

struct NativePolicyProbeBody <: PB.AbstractBody{PB.AbstractMaterial}
    n_points::Int
    policy::Union{Nothing,PB.InitialBondPolicy}
end

NativePolicyProbeBody(n_points::Int) = NativePolicyProbeBody(n_points, nothing)
PB.initial_bond_policy(body::NativePolicyProbeBody) = body.policy

function probe_policy()
    digest = repeat("a", 64)
    return PB.InitialBondPolicy("unused.json", "unused.toml", "", digest, digest,
                                PB._INITIAL_BOND_POLICY_TOKEN)
end

@testset "MPI body handler excludes policy-bearing bodies" begin
    @test PB._reject_initial_policy_in_mpi_body(NativePolicyProbeBody(1)) === nothing
    policy_body = NativePolicyProbeBody(1, probe_policy())
    @test_throws ArgumentError PB._reject_initial_policy_in_mpi_body(policy_body)
end

function probe_model(native_to_external, decisions, included, excluded)
    by_pair = Dict((left, right) => state for (left, right, state) in decisions)
    return (
        native_to_external=native_to_external,
        decision_by_pair=by_pair,
        decisions=decisions,
        expected_mask_pairs=Dict(
            "included_pairs" => included,
            "geometry_excluded_pairs" => excluded,
            "unknown_pairs" => Tuple{Int,Int}[],
        ),
        source=Dict("source_cad_sha256" => repeat("b", 64),
                    "classifier_source_sha256" => repeat("c", 64)),
        points=Dict("payload_sha256" => repeat("d", 64)),
        evidence=Dict("classifier_result_sha256" => repeat("e", 64),
                      "raw_hit_interval_sidecar_sha256" => repeat("f", 64)),
        pair_map_sha256=repeat("1", 64),
        mask_sha256=repeat("2", 64),
    )
end

@testset "initial bond policy transport helpers" begin
    csv = collect(codeunits("point_id,x_mm,y_mm,z_mm\r\n50,1.0,-0.0,3e1\r\n0,4,5,6"))
    parsed = PB._initial_policy_parse_points_csv(csv)
    @test parsed == [(50, (1.0, -0.0, 30.0)), (0, (4.0, 5.0, 6.0))]
    @test PB._initial_policy_match_points_csv(
        csv, [50, 0], [(1.0, -0.0, 30.0), (4.0, 5.0, 6.0)]) === nothing
    @test_throws ArgumentError PB._initial_policy_match_points_csv(
        csv, [0, 50], [(4.0, 5.0, 6.0), (1.0, -0.0, 30.0)])
    @test_throws ArgumentError PB._initial_policy_match_points_csv(
        csv, [50, 0], [(1.0, -0.0, 31.0), (4.0, 5.0, 6.0)])
    @test_throws ArgumentError PB._initial_policy_parse_points_csv(
        collect(codeunits("point_id,x_mm,y_mm,z_mm\n01,1,2,3\n")))
    @test_throws ArgumentError PB._initial_policy_parse_points_csv(
        collect(codeunits("point_id,x_mm,y_mm,z_mm\n0,1,2,3\n0,4,5,6\n")))
    @test_throws ArgumentError PB._initial_policy_parse_points_csv(
        collect(codeunits("point_id,x_mm,y_mm,z_mm\n0,nan,2,3\n")))

    identity = [[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]
    @test PB._initial_policy_transform(identity) == Matrix{Float64}(I, 4, 4)
    reflection = [[-1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0],
                  [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]
    @test_throws ArgumentError PB._initial_policy_transform(reflection)
end

@testset "full directed candidate comparison and canonical external IDs" begin
    body = NativePolicyProbeBody(2)
    model = probe_model([100, 7], [(7, 100, "IN_MATERIAL")], [(7, 100)],
                        Tuple{Int,Int}[])
    bonds = [PB.Bond(2, 1.0, true), PB.Bond(1, 1.0, true)]
    neighbors = [1, 1]
    receipt = PB._filter_initial_policy_candidates!(bonds, neighbors, [1, 2], body,
                                                    model, probe_policy())
    @test [bond.neighbor for bond in bonds] == [2, 1]
    @test neighbors == [1, 1]
    @test receipt.raw_candidate_pair_count == 1
    @test receipt.raw_directed_bond_count == 2
    @test receipt.included_pair_count == 1
    @test receipt.included_directed_bond_count == 2
    @test receipt.raw_candidate_set_reproduced
    @test receipt.filtered_bond_list_reproduced
    @test !receipt.geometry_exclusions_are_physics_events
    @test receipt.native_graph_readback == "NOT_RUN"
    @test !receipt.initialization_authorized

    missing_reciprocal = [PB.Bond(2, 1.0, true)]
    @test_throws ArgumentError PB._filter_initial_policy_candidates!(
        missing_reciprocal, [1, 0], [1, 2], body, model, probe_policy())
end

@testset "geometry exclusions are removed before later graph setup" begin
    body = NativePolicyProbeBody(3)
    decisions = [(0, 20, "IN_MATERIAL"), (0, 50, "IN_MATERIAL"),
                 (20, 50, "OUT_OR_VOID")]
    model = probe_model([50, 0, 20], decisions, [(0, 20), (0, 50)], [(20, 50)])
    bonds = [PB.Bond(2, 1.0, true), PB.Bond(3, 1.0, true),
             PB.Bond(1, 1.0, true), PB.Bond(3, 1.0, true),
             PB.Bond(1, 1.0, true), PB.Bond(2, 1.0, true)]
    neighbors = [2, 2, 2]
    receipt = PB._filter_initial_policy_candidates!(bonds, neighbors, [1, 2, 3], body,
                                                    model, probe_policy())
    @test [bond.neighbor for bond in bonds] == [2, 1, 3, 2]
    @test neighbors == [1, 2, 1]
    @test receipt.geometry_excluded_pair_count == 1
    @test receipt.included_directed_bond_count == 4
end
