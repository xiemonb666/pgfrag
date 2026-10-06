"""SOURCE_READY-only static checker; never imports or executes numeric modules."""
from __future__ import annotations
import argparse
import ast
import hashlib
import itertools
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
PGFRAG = HERE.parents[1]
WORKSPACE = PGFRAG.parent
REVIEW = PGFRAG / "tmp" / "review_20261006" / "mixed93_actual_cohesive_bridge_20261006"
CASE_PATH = HERE / "case.json"
PINS_PATH = REVIEW / "SOURCE_PINS.json"
MAX_DRIVER_LINES = 420
ALLOWED_IMPORT_ROOTS = {"__future__", "argparse", "ast", "hashlib", "itertools", "json", "pathlib", "sys"}
EXPECTED_CACHE_NPZ_SHA256 = "3e1e6fd17751e4e06113010f20ea33c7b69914ca90431e26bbf76124951b5d41"
EXPECTED_AUDIT_SHA256 = "b622b771f092768fc1dd73ef5fc4a7360b0f17329ff592639db48e5600ecd495"
EXPECTED_AUDIT_RECEIPT_SHA256 = "ea61fbb508a2f5ef052cc31e996f575b9bf5897ec6ca459bbb29ea2656ede703"
EXPECTED_UPSTREAM_RESULT_SHA256 = "2510eabd3e3406f243d9fe900ac5af3cfdf17071cdbe50383093de0a13adff4b"
EXPECTED_FIXED_CASE_SHA256 = "475e63bb059a8464de20f8f530ecf90d3fa3a94163f04342b8fa75b6ab230ed6"


def require(ok: bool, message: str) -> None:
    if not ok:
        raise SystemExit(f"BLOCKED: {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def imported_roots(tree: ast.AST) -> set[str]:
    roots: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            roots.update(item.name.split(".", 1)[0] for item in node.names)
        elif isinstance(node, ast.ImportFrom):
            roots.add((node.module or "").split(".", 1)[0])
    return roots


def semantic_id(scope: str, owner: str, component: str, index: int) -> str:
    return json.dumps([scope, owner, component, index], separators=(",", ":"), ensure_ascii=True)


def declared_names(tree: ast.AST) -> set[str]:
    return {node.name for node in ast.walk(tree)
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))}


def validate_case(case: dict[str, object]) -> None:
    require(case.get("schema") == "mixed93_actual_cohesive_bridge.v1", "unexpected case schema")
    require(case.get("status") == "SOURCE_READY_EXECUTION_DISABLED", "case is not SOURCE_READY")
    auth = case.get("execution_authorization")
    require(isinstance(auth, dict) and auth.get("authorized") is False,
            "numerical execution must remain disabled")
    binding, registry = case.get("source_binding"), case.get("registry")
    require(isinstance(binding, dict) and isinstance(registry, dict), "source or registry contract missing")
    axis = binding.get("axis_coordinates_m")
    expected_ids = [f"node-{i:03d}" for i in range(27)]
    require(axis == [-0.005, 0.0, 0.005], "source axis differs from frozen 3x3x3 grid")
    expected_positions = [list(point) for point in itertools.product(axis, axis, axis)]
    require(binding.get("source_ids") == expected_ids, "27 source IDs are not in frozen canonical order")
    require(binding.get("source_positions_xyz_m") == expected_positions,
            "explicit source coordinates differ from frozen product-grid order")
    require("compare the cached source_ids and source_positions_xyz_m arrays exactly"
            in binding.get("accepted_cache_identity_gate", ""),
            "future numerical path must bind the actual cached source cloud before MLS")
    components = ["ux", "uy", "uz"]
    expected_semantics = [
        semantic_id("SOURCE_BASE", source_id, component, 0)
        for source_id in expected_ids for component in components
    ] + [
        semantic_id("GLOBAL_FRONT", "global", component, branch)
        for branch in range(1, 5) for component in components
    ]
    require(registry.get("semantic_coefficient_ids_in_order") == expected_semantics,
            "93 SOURCE_BASE/GLOBAL_FRONT tokens or order diverge from frozen registry")
    require(registry.get("operator_layout") == "MIXED93_BASE81_GLOBAL12",
            "explicit mixed operator-layout label missing")
    require(registry.get("trace_family") == "GLOBALFRONT",
            "trace tag must describe enrichment family, not the whole DOF layout")
    require(registry.get("dof_count") == 93 and len(expected_semantics) == 93,
            "93-coefficient contract mismatch")
    mls, qrule, limits = case.get("mls"), case.get("interface_quadrature"), case.get("execution_budget")
    require(isinstance(mls, dict) and isinstance(qrule, dict) and isinstance(limits, dict),
            "MLS, quadrature or execution budget missing")
    require(mls.get("base_row_only") is True and mls.get("visibility_mask") is None,
            "direct continuous BASE row contract is required")
    require(mls.get("new_unique_queries") == 16 and mls.get("expected_support_count_each_query") == 27,
            "continuous MLS query/support budget changed")
    require(qrule.get("point_count") == 16 and qrule.get("area_counted_once") is True,
            "finite 4x4 one-side interface quadrature contract changed")
    require(qrule.get("one_side_total_area_m2") == 1.0e-4, "mapped face area must equal R*L")
    require(limits.get("attempts") == 1 and limits.get("coupled_kkt_solves") == 1
            and limits.get("primitive_assembly_calls_max") == 3,
            "single-run or primitive-trial budget changed")
    require(limits.get("execution_authorization") is False, "execution budget must remain disabled")
    coupled = case.get("coupled_problem")
    require(isinstance(coupled, dict)
            and coupled.get("mechanical_change_roundoff_floors") == {
                "q_change_m": "128*float64_eps*norm2(q0_m)",
                "cohesive_internal_force_n": "128*float64_eps*norm2(f_ext_scaled_N)"
            }, "mechanical-change roundoff floors must be predeclared")
    family = case.get("cohesive_law")
    require(isinstance(family, dict) and family.get("family_tag") == "GLOBALFRONT"
            and family.get("operator_layout") == "MIXED93_BASE81_GLOBAL12",
            "primitive enrichment-family/layout metadata is inconsistent")
    calls = family.get("primitive_call_budget")
    require(isinstance(calls, dict) and len(calls.get("calls", [])) == 3
            and [call.get("name") for call in calls["calls"]]
            == ["baseline", "zero_tangent_control", "coupled_final"]
            and calls.get("qps_per_call") == 16,
            "baseline, zero-tangent control and coupled-final calls must be explicit")
    zero_control = calls["calls"][1]
    require("fully-damaged" in zero_control.get("history", "")
            and "discarded" in zero_control.get("purpose", ""),
            "zero-tangent control must be synthetic and excluded from physical updates")
    cache = case.get("load_and_cache_binding")
    require(isinstance(cache, dict)
            and cache.get("accepted_arrays_sha256") == EXPECTED_CACHE_NPZ_SHA256
            and cache.get("independent_audit_sha256") == EXPECTED_AUDIT_SHA256
            and cache.get("recovery_receipt_sha256") == EXPECTED_AUDIT_RECEIPT_SHA256,
            "recovered-array, audit or receipt identity differs from root-accepted pins")
    keys = cache.get("actual_array_keys", {})
    require(keys.get("q_solution") == "pure_I_q_solution" and keys.get("q_exact") == "pure_I_q_exact",
            "actual solution and exact manufactured coefficients must remain separate")
    material = family.get("material")
    require(isinstance(material, dict) and material.get("k0_pa_per_m") == 1.0e8
            and material.get("tpeak_pa") == 1000.0 and material.get("gc_j_per_m2") == 0.01
            and material.get("delta0_m") == 1.0e-5 and material.get("deltaf_m") == 2.0e-5,
            "root-approved toy elastic-branch parameters changed")
    interface = case.get("interface_mechanics")
    coupled = case.get("coupled_problem")
    require(isinstance(interface, dict)
            and interface.get("expected_KG_f1_uy_n_per_m") == 20000.0
            and interface.get("only_nonzero_entry") == "zero-based [82,82]"
            and interface.get("analytic_matrix_relative_error_max") == 1.0e-12,
            "analytic elastic interface tangent contract changed")
    require(isinstance(coupled, dict) and coupled.get("kkt_order") == 99
            and coupled.get("solve_count") == 1,
            "exactly one 99-order coupled KKT solve is permitted")
    require(any("no real material" in claim for claim in case.get("nonclaims", [])),
            "material nonclaim missing")


def validate_recovery_metadata(case: dict[str, object], pins: dict[str, object],
                                paths: dict[str, Path]) -> None:
    cache = case["load_and_cache_binding"]
    refs = pins.get("references")
    require(isinstance(refs, dict), "source pin references are missing")
    audit = json.loads(paths["recovery_audit"].read_text(encoding="utf-8"))
    receipt = json.loads(paths["recovery_run_receipt"].read_text(encoding="utf-8"))
    require(audit.get("status") == "RECOVERED_CACHE_INDEPENDENTLY_CONSISTENT_PHYSICS_UNQUALIFIED",
            "accepted recovery audit status changed")
    require(audit.get("pins_post", {}).get("failed_arrays") == cache["accepted_arrays_sha256"],
            "recovery audit no longer binds accepted array digest")
    post_pins = audit.get("pins_post", {})
    require(post_pins.get("frozen_entry") == cache["producer_frozen_driver_sha256"]
            and post_pins.get("case") == cache["cached_source_case_sha256"],
            "recovered arrays no longer bind the reviewed producer source/case")
    spectrum = audit.get("matrix_layout_spectrum_rigid_audit", {}).get("spectrum_counts", {})
    require(spectrum.get("zero") == 6 and spectrum.get("positive") == 87
            and spectrum.get("negative") == 0, "recovery audit spectrum gate changed")
    gates = audit.get("consistency_gates", {})
    require(all(gates.get(key) is True for key in
                ("matrix_and_layout", "spectrum_and_rigid_modes", "four_case_algebra_and_fields")),
            "independent recovery gates did not pass")
    pure_i = audit.get("case_audits", {}).get("pure_I", {})
    require(pure_i.get("KKT_saved_algebra_gate") is True
            and pure_i.get("saved_solution_field_reconstruction_gate") is True,
            "actual pure-I q_solution is not independently algebra-checked")
    require(receipt.get("process", {}).get("outer_return_code") == 0,
            "independent recovery process did not complete")
    outputs = receipt.get("generated_outputs", [])
    require(any(item.get("sha256") == cache["independent_audit_sha256"]
                and str(item.get("path", "")).replace("\\", "/").endswith("/audit.json")
                for item in outputs), "recovery receipt does not bind the audit JSON")


def validate_upstream_metadata(case: dict[str, object], paths: dict[str, Path]) -> None:
    binding = case["load_and_cache_binding"]
    result = json.loads(paths["upstream_result"].read_text(encoding="utf-8"))
    fixed = json.loads(paths["upstream_fixed_case"].read_text(encoding="utf-8"))
    source = case["source_binding"]
    require(result.get("registry", {}).get("source_ids") == source["source_ids"],
            "upstream result source IDs differ from the explicit new case source order")
    require(result.get("study_case_sha256") == "9d1b2153a1a3db36387d838f584714dae8595c722f542155e83bfe1f037e1c5d",
            "upstream result does not match the pinned tip-aligned case")
    geometry = fixed.get("geometry", {})
    sources = fixed.get("sources", {})
    mls = fixed.get("mls", {})
    require(geometry.get("radius_m") == 0.01 and geometry.get("length_m") == 0.01
            and geometry.get("z_bounds_m") == [-0.005, 0.005],
            "upstream coordinate frame differs from the frozen centered cylinder")
    require(sources.get("axis_coordinates_m") == source["axis_coordinates_m"]
            and sources.get("source_count") == 27
            and sources.get("source_id_format") == "node-{row:03d}",
            "upstream 27-source geometry differs from the explicit new case")
    expected_mls = case["mls"]
    require(mls.get("reference_origin_m") == expected_mls["reference_origin_m"]
            and mls.get("reference_scale_m") == expected_mls["reference_scale_m"]
            and mls.get("support_radius_m") == expected_mls["support_radius_m"]
            and mls.get("max_condition_number") == expected_mls["max_condition_number"]
            and mls.get("visibility_mask") is None,
            "continuous MLS configuration differs from its frozen provider contract")
    require(binding.get("upstream_result_sha256") == EXPECTED_UPSTREAM_RESULT_SHA256
            and binding.get("upstream_fixed_case_sha256") == EXPECTED_FIXED_CASE_SHA256,
            "upstream source/result identity differs from accepted metadata pins")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="AST and metadata-only SOURCE_READY check.")
    parser.add_argument("--execute", action="store_true", help="always refused by this source-only preparation")
    args = parser.parse_args(argv)
    if args.execute:
        raise SystemExit("execution disabled: root review must authorize a separate numerical implementation")
    driver_path = Path(__file__).resolve()
    source_text = driver_path.read_text(encoding="utf-8")
    require(len(source_text.splitlines()) <= MAX_DRIVER_LINES, f"driver exceeds {MAX_DRIVER_LINES} lines")
    own_tree = ast.parse(source_text, filename=str(driver_path))
    require(imported_roots(own_tree) <= ALLOWED_IMPORT_ROOTS,
            "SOURCE_READY driver imports a numerical or non-stdlib module")
    case = json.loads(CASE_PATH.read_text(encoding="utf-8"))
    pins = json.loads(PINS_PATH.read_text(encoding="utf-8"))
    require(pins.get("schema") == "mixed93_actual_cohesive_bridge_source_pins.v1",
            "unexpected SOURCE_PINS schema")
    require(pins.get("execution_authorized") is False, "pins must keep execution disabled")
    validate_case(case)
    refs = pins.get("references")
    require(isinstance(refs, dict) and refs, "reference pin map is empty")
    accepted = pins.get("accepted_cache_binding")
    binding = case["load_and_cache_binding"]
    require(isinstance(accepted, dict)
            and accepted.get("arrays_path") == binding["failed_serialization_arrays_path"]
            and accepted.get("arrays_sha256") == binding["accepted_arrays_sha256"]
            and accepted.get("audit_sha256") == binding["independent_audit_sha256"]
            and accepted.get("audit_receipt_sha256") == binding["recovery_receipt_sha256"]
            and pins.get("producer_frozen_driver_sha256") == binding["producer_frozen_driver_sha256"]
            and accepted.get("npz_opened_by_static_check") is False,
            "SOURCE_PINS cache binding differs from case metadata")
    require(refs.get("recovery_audit", {}).get("path") == binding["independent_audit_path"]
            and refs.get("recovery_run_receipt", {}).get("path") == binding["recovery_receipt_path"]
            and refs.get("upstream_result", {}).get("path") == binding["upstream_result_path"]
            and refs.get("upstream_fixed_case", {}).get("path") == binding["upstream_fixed_case_path"],
            "case paths and source-pin paths diverge")
    resolved: dict[str, Path] = {}
    trees: dict[str, ast.AST] = {}
    for name, record in refs.items():
        require(isinstance(record, dict) and isinstance(record.get("path"), str),
                f"malformed source pin record: {name}")
        path = WORKSPACE / record["path"]
        require(path.is_file(), f"pinned source or metadata is absent: {record['path']}")
        require(path.suffix.lower() != ".npz", "binary cached arrays must not be opened in SOURCE_READY")
        require(sha256(path) == record.get("sha256"), f"SHA-256 mismatch for {name}")
        resolved[name] = path
        if path.suffix.lower() == ".py":
            trees[name] = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    required_symbols = {
        "mls3d_mls": {"MLS3DConfig", "MLSStatus"},
        "mls3d_shape": {"ShapeFunctions3DResult", "shape_functions_3d"},
        "mixed_lefm_field": {"branch_functions_local"},
        "mixed_enriched_operator": {"semantic_coefficient_ids", "rigid_body_mode_coefficients"},
        "trace_primitive": {"TraceFamily", "EnrichedFacePoint", "RigidTraceCandidates", "assemble_enriched_trial"},
        "cohesive_scalar_law": {"BilinearCohesiveMaterial", "initial_state", "update_state"},
    }
    for source_name, names in required_symbols.items():
        require(source_name in trees and names <= declared_names(trees[source_name]),
                f"required source API symbols missing from {source_name}")
    require(resolved.get("entry_driver") == driver_path, "entry-driver pin targets a different file")
    require(refs["case"]["sha256"] == sha256(CASE_PATH), "case pin differs from SOURCE_PINS")
    require(refs["entry_driver"]["sha256"] == sha256(driver_path), "driver pin differs from SOURCE_PINS")
    require(refs["recovery_audit"]["sha256"] == EXPECTED_AUDIT_SHA256
            and refs["recovery_run_receipt"]["sha256"] == EXPECTED_AUDIT_RECEIPT_SHA256,
            "accepted audit/receipt pins differ from root-reviewed values")
    validate_recovery_metadata(case, pins, resolved)
    validate_upstream_metadata(case, resolved)
    require(case["load_and_cache_binding"]["accepted_arrays_sha256"] == EXPECTED_CACHE_NPZ_SHA256,
            "accepted array digest changed; array contents are intentionally not read")
    print(json.dumps({
        "status": "SOURCE_READY",
        "scope": "AST parsing and metadata SHA-256 only",
        "numeric_modules_imported": False,
        "npz_opened_or_hashed": False,
        "provider_called": False,
        "solver_called": False,
        "source_count": len(case["source_binding"]["source_ids"]),
        "dof_count": case["registry"]["dof_count"],
        "new_mls_queries": case["mls"]["new_unique_queries"],
        "interface_qps": case["interface_quadrature"]["point_count"],
        "primitive_trial_assembly_max": case["cohesive_law"]["primitive_call_budget"]["assemble_enriched_trial_max"],
        "coupled_kkt_solves": case["coupled_problem"]["solve_count"],
        "execution_authorized": False
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


