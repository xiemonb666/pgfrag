"""Repeat the declared tet-stratified point-IN check on a saved ROI MSH without remeshing."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

import cadquery as cq
import numpy as np
from OCP.BRepClass3d import BRepClass3d_SolidClassifier

import mesh_audit as ma
import run_roi_trial as trial


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    mesh_path = args.mesh.resolve(strict=True)
    output_path = args.output.resolve()
    mesh = ma.mesh_from_msh(mesh_path)
    count = min(trial.INTERIOR_SAMPLE_COUNT, int(mesh["tet_xyz"].shape[0]))
    selected = np.linspace(0, mesh["tet_xyz"].shape[0] - 1, count, dtype=np.int64)
    rng = np.random.default_rng(trial.INTERIOR_SAMPLE_SEED)
    barycentric_weights = rng.dirichlet(np.ones(4, dtype=np.float64), size=count)
    samples = np.einsum("ij,ijk->ik", barycentric_weights, mesh["tet_xyz"][selected])

    source_hash = sha256(trial.STEP)
    if source_hash.lower() != trial.EXPECTED_SHA.lower():
        raise ValueError(f"canonical source STEP SHA mismatch: {source_hash}")
    solid = cq.importers.importStep(str(trial.STEP)).val().Solids()[0]
    checker = BRepClass3d_SolidClassifier()
    checker.Load(solid.wrapped)
    labels = [trial.classify_one(checker, point) for point in samples]
    non_in_samples: list[dict[str, Any]] = []
    for sample_index, label in enumerate(labels):
        if label == "IN":
            continue
        cell_index = int(selected[sample_index])
        point = samples[sample_index]
        vertex = cq.Vertex.makeVertex(*map(float, point))
        non_in_samples.append({
            "sample_ordinal": int(sample_index),
            "tetra_cell_id": int(mesh["cell_ids"][cell_index]),
            "tetra_centroid_xyz_mm": mesh["centroids"][cell_index].tolist(),
            "barycentric_weights": barycentric_weights[sample_index].tolist(),
            "sample_xyz_mm": point.tolist(),
            "source_state": label,
            "distance_to_native_solid_mm": float(vertex.distance(solid)),
        })
    counts = {state: int(labels.count(state)) for state in ("IN", "OUT", "UNKNOWN")}
    result = {
        "experiment_id": "pg.wall-resolution-local-native-cad-pointin-audit/0.1.0",
        "python_executable": sys.executable,
        "argv": sys.argv,
        "method": "read-only repeat of the exact deterministic tet-stratified Dirichlet samples from the stored MSH; no meshing",
        "source_step": str(trial.STEP),
        "source_step_sha256": source_hash,
        "mesh_msh": str(mesh_path),
        "mesh_msh_sha256": sha256(mesh_path),
        "sample_seed": trial.INTERIOR_SAMPLE_SEED,
        "sample_count": int(count),
        "state_counts": counts,
        "non_IN_samples": non_in_samples,
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({
        "output": str(output_path),
        "output_sha256": sha256(output_path),
        "mesh_msh_sha256": result["mesh_msh_sha256"],
        "state_counts": counts,
        "non_IN_samples": non_in_samples,
    }, indent=2))
    return 0 if counts["IN"] == count else 1


if __name__ == "__main__":
    raise SystemExit(main())
