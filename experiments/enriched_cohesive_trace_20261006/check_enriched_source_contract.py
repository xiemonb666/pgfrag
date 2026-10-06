"""Static AST/JSON contract check only; never imports or runs any cohesive module."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
PGFRAG = HERE.parents[1]
SOURCE = HERE / "enriched_cohesive_trace.py"
MICRO = HERE / "check_enriched_numerics_future.py"
SCHEMA = HERE / "enriched_trace_contract.schema.json"
CASE = HERE / "case_globalfront_toy.json"
SCALAR = HERE.parent / "cohesive_scalar_law_20261006" / "cohesive_scalar_law.py"
DERIVATIVE = HERE.parent / "cohesive_interface_assembly_20261006" / "cohesive_interface_assembly.py"
RECEIPT = PGFRAG / "tmp" / "review_20261006" / "enriched_cohesive_trace_20261006" / "source_check_receipt.json"
SCALAR_PIN = "78433638608ad4216979607412fe0d4966f1423706c47f734f50dd2a8121fa17"
DERIVATIVE_PIN = "d2bf0dd79507e5ab26efe69b5489c3284b7c9118656ab49fffe942303ba11854"
MAX_SOURCE_LINES, MAX_MICRO_LINES, MAX_CHECKER_LINES = 330, 170, 150
STDLIB = {"__future__", "ast", "collections", "dataclasses", "enum", "hashlib", "importlib", "json", "math", "numbers", "pathlib", "sys", "time", "unittest"}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"FAIL: {message}")


def roots(tree: ast.AST) -> set[str]:
    out = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import): out.update(x.name.split(".")[0] for x in node.names)
        elif isinstance(node, ast.ImportFrom): out.add((node.module or "").split(".")[0])
    return out


def declared(tree: ast.AST) -> set[str]:
    return {x.name for x in ast.walk(tree) if isinstance(x, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))}


def main() -> None:
    source_text = SOURCE.read_text(encoding="utf-8")
    micro_text = MICRO.read_text(encoding="utf-8")
    checker_text = Path(__file__).read_text(encoding="utf-8")
    source_tree = ast.parse(source_text, filename=str(SOURCE))
    micro_tree = ast.parse(micro_text, filename=str(MICRO))
    schema = json.loads(SCHEMA.read_text(encoding="utf-8"))
    case = json.loads(CASE.read_text(encoding="utf-8"))
    source_lines, micro_lines, checker_lines = len(source_text.splitlines()), len(micro_text.splitlines()), len(checker_text.splitlines())
    require(source_lines <= MAX_SOURCE_LINES, f"implementation {source_lines} lines exceeds {MAX_SOURCE_LINES}")
    require(micro_lines <= MAX_MICRO_LINES, f"future micro-checker {micro_lines} lines exceeds {MAX_MICRO_LINES}")
    require(checker_lines <= MAX_CHECKER_LINES, f"static checker {checker_lines} lines exceeds {MAX_CHECKER_LINES}")
    require("MAX_COEFFICIENTS, MAX_QUADRATURE_POINTS = 93, 16" in source_text, "enriched limits must remain 93 coefficients / 16 QPs")
    require(roots(source_tree) <= STDLIB and roots(micro_tree) <= STDLIB, "implementation/future checker must use only standard library")
    required = {"RigidTraceCandidates", "EnrichedFacePoint", "PointTrialResult", "EnrichedTrialResult", "assemble_enriched_trial", "_validate_rigid_trace", "_allocate_tangent"}
    require(required <= declared(source_tree), "enriched trace API or rigid evidence validator is missing")
    api = next(x for x in source_tree.body if isinstance(x, ast.FunctionDef) and x.name == "assemble_enriched_trial")
    require(api.args.args[-1].arg == "rigid_candidates" and not api.args.defaults, "rigid evidence must be mandatory at the assembly API")
    calls = [x for x in ast.walk(api) if isinstance(x, ast.Call) and isinstance(x.func, ast.Name) and x.func.id == "_allocate_tangent"]
    qp_cap = [x for x in ast.walk(api) if isinstance(x, ast.If) and "MAX_QUADRATURE_POINTS" in ast.unparse(x.test)]
    require(len(calls) == 1 and qp_cap and calls[0].lineno > qp_cap[0].lineno, "16-QP cap must precede dense tangent allocation")
    for pin in (SCALAR_PIN, DERIVATIVE_PIN): require(pin in source_text, "runtime dependency SHA pins must be embedded")
    require("_algorithmic_tangent" in source_text and "_scalar_update" in source_text and "_scalar_traction" in source_text, "shared scalar law/derivative must replace any duplicate law")
    required_sections = {"status", "units", "limits", "point_fields", "rigid_evidence", "history", "result_fields", "dependencies", "qualification"}
    require(schema.get("$schema") == "https://json-schema.org/draft/2020-12/schema", "JSON schema dialect missing")
    require(set(schema.get("required", ())) == required_sections, "schema must describe all nine interface contract sections")
    require(case.get("status") == "SOURCE_ONLY_NOT_COUPLED" and case.get("basis_family") == "GLOBALFRONT", "case manifest must state scope")
    require(len(case.get("coefficient_ids_in_order", ())) == 24 and case["point"]["one_side_area_weight_m2"] == 2.0, "future toy case must be 24 coefficients and one-side area")
    require(len(case["point"]["rigid_candidate_sparse_entries"]) == 6, "toy case must declare six rigid modes")
    require(case["point"]["positive_reference_point_m"] == [0.2, -0.1, 0.3] and case["point"]["opening_gradient_a"] == [0.5, -0.25, 1.0], "fixture must exercise non-origin rigid traces and a general opening gradient")
    require(all(token in micro_text for token in ("directional_fd", "range(6,9)", "CALLS == 21", "MAX_SECONDS, CALLS = 15.0, 0", "expected_area_t_a_n", "expected_area_k_a_outer_a_n_per_m", "action_reaction_max_residual_n", "future_run01/numeric_result.json")), "future micro-checker tangent/force/JSON budget plan incomplete")
    scalar_sha = hashlib.sha256(SCALAR.read_bytes()).hexdigest()
    derivative_sha = hashlib.sha256(DERIVATIVE.read_bytes()).hexdigest()
    require(scalar_sha == SCALAR_PIN and derivative_sha == DERIVATIVE_PIN, "frozen cohesive dependency source changed")
    source_sha = hashlib.sha256(source_text.encode()).hexdigest()
    receipt = {"status": "PASS", "check_kind": "AST parse and JSON/schema structure only", "cohesive_source_executed": False, "scalar_law_executed": False, "nodal_derivative_executed": False, "future_micro_checker_executed": False, "source_lines": source_lines, "source_line_limit": MAX_SOURCE_LINES, "future_micro_checker_lines": micro_lines, "future_micro_checker_line_limit": MAX_MICRO_LINES, "static_checker_lines": checker_lines, "static_checker_line_limit": MAX_CHECKER_LINES, "source_sha256": source_sha, "future_micro_checker_sha256": hashlib.sha256(micro_text.encode()).hexdigest(), "scalar_source_sha256": scalar_sha, "scalar_source_pin": SCALAR_PIN, "nodal_derivative_sha256": derivative_sha, "nodal_derivative_pin": DERIVATIVE_PIN, "case_coefficients": len(case["coefficient_ids_in_order"]), "case_qps": 1, "schema_sections": sorted(required_sections)}
    RECEIPT.parent.mkdir(parents=True, exist_ok=True)
    RECEIPT.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__": main()
