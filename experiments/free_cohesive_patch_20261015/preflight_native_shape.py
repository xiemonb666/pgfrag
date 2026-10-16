from __future__ import annotations

import ast
import hashlib
import json
import os
import shutil
import sys
import time
import traceback
from pathlib import Path
from typing import Any

import numpy as np


EXP = Path(__file__).resolve().parent
CASE_PATH = EXP / "case.json"
MESH_PATH = EXP / "block_2x2x4mm.msh"
MATERIAL_PATH = EXP / "material.dat"
OUT = Path("/home/xiemo/free_cohesive_patch_20261015/preflight_native_shape_20261015_02")
WALL_CAP_S = 60.0
OUTPUT_CAP_BYTES = 32 * 1024 * 1024
SOURCE_NAMES = (
    "case.json",
    "block_2x2x4mm.msh",
    "material.dat",
    "preflight_native_shape.py",
    "run_native_free_patch.py",
    "analyze_native_free_patch.py",
    "free_patch_contract.py",
    "test_source_contract.py",
    "launch_free_patch.sh",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def output_size_bytes(root: Path) -> int:
    return sum(path.stat().st_size for path in root.rglob("*") if path.is_file())


def called_solve_step(source: str) -> bool:
    tree = ast.parse(source)
    return any(
        isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr == "solveStep"
        for node in ast.walk(tree)
    )


def run() -> dict[str, Any]:
    started = time.monotonic()
    OUT.mkdir(parents=True, exist_ok=False)
    snapshot_dir = OUT / "source_bytes"
    snapshot_dir.mkdir()
    snapshot_rows = []
    for name in SOURCE_NAMES:
        source_path = EXP / name
        snapshot_path = snapshot_dir / name
        shutil.copyfile(source_path, snapshot_path)
        snapshot_rows.append({
            "path": name,
            "bytes": snapshot_path.stat().st_size,
            "sha256": sha256(snapshot_path),
        })

    driver_text = (EXP / "preflight_native_shape.py").read_text(encoding="utf-8")
    if called_solve_step(driver_text):
        raise RuntimeError("preflight source contains a solveStep() call")
    case = json.loads(CASE_PATH.read_text(encoding="utf-8"))
    if sha256(MESH_PATH) != case["geometry"]["mesh_sha256"]:
        raise RuntimeError("frozen mesh SHA-256 mismatch")
    if sha256(MATERIAL_PATH) != case["material"]["material_sha256"]:
        raise RuntimeError("frozen material SHA-256 mismatch")

    import akantu as aka

    aka.parseInput(str(MATERIAL_PATH))
    mesh = aka.Mesh(3)
    mesh.read(str(MESH_PATH))
    input_nodes = np.asarray(mesh.getNodes()).copy()
    input_tets = np.asarray(mesh.getConnectivity(aka._tetrahedron_4)).copy()
    model = aka.SolidMechanicsModelCohesive(mesh)
    model.initFull(_analysis_method=aka._explicit_lumped_mass, _is_extrinsic=False)
    model.assembleMassLumped()

    native_nodes = np.asarray(mesh.getNodes()).copy()
    native_tets = np.asarray(mesh.getConnectivity(aka._tetrahedron_4)).copy()
    native_cohesive = np.asarray(mesh.getConnectivity(aka._cohesive_3d_6)).copy()
    mass = np.asarray(model.getMass()).copy()
    bulk_stable_step = float(model.getStableTimeStep())
    fields = {
        name: np.asarray(model.flattenInternalReal(name, aka._ek_cohesive, aka._not_ghost)(aka._cohesive_3d_6, aka._not_ghost)).copy()
        for name in ("opening", "tractions", "damage", "delta max")
    }
    if native_cohesive.ndim != 2 or native_cohesive.shape[1] != 6:
        raise RuntimeError(f"native cohesive connectivity shape unexpected: {native_cohesive.shape}")
    if not fields["opening"].ndim == 2 or fields["opening"].shape[1] != 3:
        raise RuntimeError(f"native opening shape unexpected: {fields['opening'].shape}")
    if fields["tractions"].shape != fields["opening"].shape:
        raise RuntimeError(f"native traction/opening shape mismatch: {fields['tractions'].shape}, {fields['opening'].shape}")
    if fields["damage"].shape != (len(fields["opening"]), 1) or fields["delta max"].shape != (len(fields["opening"]), 1):
        raise RuntimeError("native scalar cohesive field shapes do not match the opening point count")
    if len(native_cohesive) == 0 or len(fields["opening"]) % len(native_cohesive) != 0:
        raise RuntimeError("native cohesive point rows do not divide evenly by element count")
    points_per_element = len(fields["opening"]) // len(native_cohesive)
    if mass.shape != (len(native_nodes), 3):
        raise RuntimeError(f"native lumped mass shape unexpected: {mass.shape}")
    if not all(np.isfinite(value).all() for value in (*fields.values(), mass)):
        raise RuntimeError("non-finite native initial cohesive field or mass")

    npz_path = OUT / "native_zero_step_state.npz"
    np.savez_compressed(
        npz_path,
        input_nodes=input_nodes,
        input_bulk_tets=input_tets,
        native_nodes=native_nodes,
        native_bulk_tets=native_tets,
        native_cohesive_connectivity=native_cohesive,
        mass_kg=mass,
        opening=fields["opening"],
        tractions=fields["tractions"],
        damage=fields["damage"],
        delta_max=fields["delta max"],
    )
    if time.monotonic() - started > WALL_CAP_S:
        raise TimeoutError(f"zero-step preflight exceeded {WALL_CAP_S} seconds")

    elapsed = time.monotonic() - started
    result = {
        "status": "PREFLIGHT_PASS_ZERO_SOLVER_STEPS",
        "actual_argv": sys.argv,
        "solver_steps": 0,
        "solve_step_calls_in_preflight_source": False,
        "elapsed_wall_seconds": elapsed,
        "wall_cap_seconds": WALL_CAP_S,
        "output_cap_bytes": OUTPUT_CAP_BYTES,
        "case_id": case["case_id"],
        "counts": {
            "input_nodes": int(len(input_nodes)),
            "input_bulk_tetrahedra": int(len(input_tets)),
            "native_nodes": int(len(native_nodes)),
            "native_bulk_tetrahedra": int(len(native_tets)),
            "native_cohesive_3d6": int(len(native_cohesive)),
            "native_quadrature_points_per_cohesive_element": int(points_per_element),
            "native_total_cohesive_quadrature_points": int(len(fields["opening"])),
        },
        "native_field_shapes": {name: list(value.shape) for name, value in fields.items()},
        "native_initial_max_abs": {
            "opening_m": float(np.max(np.abs(fields["opening"]), initial=0.0)),
            "traction_Pa": float(np.max(np.abs(fields["tractions"]), initial=0.0)),
            "damage": float(np.max(np.abs(fields["damage"]), initial=0.0)),
            "delta_max_m": float(np.max(np.abs(fields["delta max"]), initial=0.0)),
        },
        "mass": {
            "shape": list(mass.shape),
            "component_sums_kg": [float(mass[:, component].sum()) for component in range(3)],
            "total_component_sum_kg": float(mass[:, 0].sum()),
        },
        "bulk_only_get_stable_time_step_s": bulk_stable_step,
        "total_cohesive_plus_bulk_stability": "NOT_ESTABLISHED; native getStableTimeStep excludes cohesive elements",
        "source_snapshot": snapshot_rows,
        "native_state_archive": {
            "path": npz_path.name,
            "bytes": npz_path.stat().st_size,
            "sha256": sha256(npz_path),
        },
        "runtime": {
            "python": sys.version.split()[0],
            "akantu_version": str(getattr(aka, "__version__", "unreported")),
            "backend": "isolated workspace WSL runtime",
        },
    }
    if output_size_bytes(OUT) > OUTPUT_CAP_BYTES:
        raise RuntimeError(f"preflight outputs exceed {OUTPUT_CAP_BYTES} bytes")
    write_json(OUT / "preflight_result.json", result)
    receipt = {
        "receipt_status": result["status"],
        "actual_argv": sys.argv,
        "started_unix_time_s": time.time() - elapsed,
        "elapsed_wall_seconds": elapsed,
        "solver_steps": 0,
        "source_snapshot_sha256s": {row["path"]: row["sha256"] for row in snapshot_rows},
        "artifact_sha256s": {
            "native_zero_step_state.npz": sha256(npz_path),
            "preflight_result.json": sha256(OUT / "preflight_result.json"),
        },
        "output_size_bytes_including_receipt": output_size_bytes(OUT),
        "output_cap_bytes": OUTPUT_CAP_BYTES,
        "wall_cap_seconds": WALL_CAP_S,
    }
    write_json(OUT / "receipt.json", receipt)
    if output_size_bytes(OUT) > OUTPUT_CAP_BYTES:
        raise RuntimeError(f"preflight outputs including receipt exceed {OUTPUT_CAP_BYTES} bytes")
    return {"status": result["status"], "output_directory": str(OUT), "receipt": receipt}


def main() -> int:
    try:
        summary = run()
    except Exception as error:
        failure = {
            "status": "FAIL_NO_SOLVER_STEPS",
            "actual_argv": sys.argv,
            "error": repr(error),
            "traceback": traceback.format_exc(),
        }
        if OUT.is_dir():
            write_json(OUT / "preflight_failure.json", failure)
        print(json.dumps(failure, sort_keys=True))
        return 2
    print(json.dumps(summary, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
