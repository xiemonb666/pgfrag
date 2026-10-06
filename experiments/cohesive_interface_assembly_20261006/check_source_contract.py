"""stdlib-only AST/schema check; it never imports or runs the assembler."""

from __future__ import annotations

import ast
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
PGFRAG = HERE.parents[1]
SOURCE = HERE / "cohesive_interface_assembly.py"
SCHEMA = HERE / "interface_contract.schema.json"
NUMERIC = HERE / "check_interface_numerics.py"
SCALAR = HERE.parent / "cohesive_scalar_law_20261006" / "cohesive_scalar_law.py"
SCALAR_PIN = "78433638608AD4216979607412FE0D4966F1423706C47F734F50DD2A8121FA17"
RECEIPT = PGFRAG / "tmp" / "review_20261006" / "cohesive_interface_assembly_20261006" / "source_check_receipt.json"
MAX_SOURCE_LINES, MAX_CHECKER_LINES, MAX_NUMERIC_LINES = 300, 150, 185


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(f"FAIL: {message}")


def names(tree: ast.AST) -> set[str]:
    return {node.name for node in ast.walk(tree) if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))}


def main() -> None:
    source_text = SOURCE.read_text(encoding="utf-8")
    source_tree = ast.parse(source_text, filename=str(SOURCE))
    numeric_text = NUMERIC.read_text(encoding="utf-8")
    numeric_tree = ast.parse(numeric_text, filename=str(NUMERIC))
    schema = json.loads(SCHEMA.read_text(encoding="utf-8"))
    source_lines, checker_lines = len(source_text.splitlines()), len(Path(__file__).read_text(encoding="utf-8").splitlines())
    numeric_lines = len(numeric_text.splitlines())
    source_sha = hashlib.sha256(source_text.encode()).hexdigest()
    require(source_lines <= MAX_SOURCE_LINES, f"assembler is {source_lines} lines; limit is {MAX_SOURCE_LINES}")
    require(checker_lines <= MAX_CHECKER_LINES, f"checker is {checker_lines} lines; limit is {MAX_CHECKER_LINES}")
    require(numeric_lines <= MAX_NUMERIC_LINES, f"numeric checker is {numeric_lines} lines; limit is {MAX_NUMERIC_LINES}")
    required_names = {"CohesiveQuadraturePoint", "PointTrialResult", "InterfaceTrialResult", "assemble_trial", "_validate_geometry", "_algorithmic_tangent"}
    require(required_names <= names(source_tree), "required interface API/validation declarations are absent")
    constants = {}
    for node in ast.walk(source_tree):
        if isinstance(node, ast.Assign) and isinstance(node.value, (ast.Tuple, ast.List)) and all(isinstance(x, ast.Constant) for x in node.value.elts):
            for target in node.targets:
                if isinstance(target, ast.Tuple):
                    constants.update((name.id, value.value) for name, value in zip(target.elts, node.value.elts) if isinstance(name, ast.Name))
        elif isinstance(node, ast.Assign) and isinstance(node.value, ast.Constant):
            for target in node.targets:
                if isinstance(target, ast.Name): constants[target.id] = node.value.value
    require((constants.get("MAX_NODES"), constants.get("MAX_DOF"), constants.get("MAX_QUADRATURE_POINTS")) == (12, 36, 16), "source caps must be 12 nodes / 36 DOFs / 16 QPs")
    require((constants.get("_FORCE_ABS_TOL_N"), constants.get("_FORCE_REL_TOL")) == (1.0e-12, 1.0e-10), "source force-balance tolerances must retain N and dimensionless units")
    imported_roots = set()
    for node in ast.walk(source_tree):
        if isinstance(node, ast.Import):
            imported_roots.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            imported_roots.add((node.module or "").split(".")[0])
    require(imported_roots <= {"__future__", "collections", "dataclasses", "enum", "importlib", "math", "numbers", "pathlib", "sys"}, "assembler imports a non-stdlib dependency")
    numeric_imports = set()
    for node in ast.walk(numeric_tree):
        if isinstance(node, ast.Import): numeric_imports.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom): numeric_imports.add((node.module or "").split(".")[0])
    require(numeric_imports <= {"__future__", "hashlib", "importlib", "json", "math", "pathlib", "sys", "time", "unittest"}, "numeric checker imports a non-stdlib dependency")
    require({"directional_fd", "resultant_and_torque", "main"} <= names(numeric_tree), "focused numerical checker declarations are incomplete")
    require("CALLS == 17" in numeric_text and "MAX_SECONDS = 10.0" in numeric_text, "checker call count or time budget is not pinned")
    require(f'ASSEMBLY_PIN = "{source_sha}"' in numeric_text, "numeric checker must pin the reviewed interface source before import")
    classes = [node for node in source_tree.body if isinstance(node, ast.ClassDef)]
    frozen = {node.name for node in classes if any(isinstance(d, ast.Call) and getattr(d.func, "id", None) == "dataclass" and any(k.arg == "frozen" and isinstance(k.value, ast.Constant) and k.value.value is True for k in d.keywords) for d in node.decorator_list)}
    require({"CohesiveQuadraturePoint", "PointTrialResult", "InterfaceTrialResult"} <= frozen, "interface input/results must be frozen dataclasses")
    api = next(node for node in source_tree.body if isinstance(node, ast.FunctionDef) and node.name == "assemble_trial")
    api_calls = {node.func.id for node in ast.walk(api) if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)}
    require({"_scalar_update", "_scalar_traction", "_validate_geometry", "_force_close_n", "_allocate_tangent"} <= api_calls, "trial assembly must use scalar history, force balance, geometry, and bounded tangent allocation")
    require("commit" not in names(source_tree), "trial module must not define a history commit operation")
    require(schema.get("$schema") == "https://json-schema.org/draft/2020-12/schema" and schema.get("type") == "object", "contract schema header is invalid")
    required_schema = {"units", "point_fields", "history", "result_fields", "caps", "limits"}
    require(set(schema.get("required", ())) == required_schema, "schema must describe all six contract sections")
    properties = schema.get("properties", {})
    for section in required_schema:
        require(section in properties and properties[section].get("type") == "object", f"schema section {section} is missing")
        require(bool(properties[section].get("required")), f"schema section {section} has no required fields")
    field_properties = properties["point_fields"].get("properties", {})
    require("J = H_positive - H_negative" in field_properties["jump_operator"].get("const", ""), "schema must pin the two-face jump operator")
    require("one physical face side once" in properties["units"]["properties"]["area_weight"].get("const", ""), "schema must pin one-side area convention")
    caps = properties["caps"]["properties"]
    require((caps["max_nodes"].get("const"), caps["max_dof"].get("const"), caps["max_quadrature_points"].get("const")) == (12, 36, 16), "schema caps must match the nodal primitive")
    require("162-enriched" in caps["coupling_scope"].get("const", ""), "schema must exclude 162-enriched coupling")
    require(properties["units"]["properties"]["force_balance_abs_tolerance"].get("const") == "1e-12 N", "force tolerance must carry Newton units")
    require(properties["units"]["properties"]["force_balance_rel_tolerance"].get("const") == "1e-10 dimensionless", "relative force tolerance must be dimensionless")
    require(SCALAR.is_file(), "existing scalar law source is not available for the declared import")
    scalar_sha = hashlib.sha256(SCALAR.read_bytes()).hexdigest().upper()
    require(scalar_sha == SCALAR_PIN, "frozen scalar-law source hash changed")
    require(SCALAR_PIN in numeric_text, "numeric checker must pin the frozen scalar-law source before import")
    receipt = {"status": "PASS", "check_kind": "stdlib AST and JSON-schema structure only", "assembler_executed": False, "numeric_checker_executed": False, "source_lines": source_lines, "source_line_limit": MAX_SOURCE_LINES, "checker_lines": checker_lines, "checker_line_limit": MAX_CHECKER_LINES, "numeric_checker_lines": numeric_lines, "numeric_checker_line_limit": MAX_NUMERIC_LINES, "source_sha256": source_sha, "numeric_source_sha256": hashlib.sha256(numeric_text.encode()).hexdigest(), "scalar_source_sha256": scalar_sha, "scalar_source_pin": SCALAR_PIN, "numeric_assembly_source_pin": source_sha, "schema_sections": sorted(required_schema), "checks": ["AST parse only", "stdlib-only imports", "frozen interface dataclasses", "trial uses existing scalar law", "no commit API", "unit/area/J schema fields", "caps and Newton-force tolerance", "reviewed assembly and scalar source pins", "numerical checker source present but not run"]}
    RECEIPT.parent.mkdir(parents=True, exist_ok=True)
    RECEIPT.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
