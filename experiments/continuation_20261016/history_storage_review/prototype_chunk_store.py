"""Tiny lossless chunk/manifest prototype; this is not PG solver code.

Run with the workspace NumPy environment. All generated files stay below this
experiment directory. The test verifies atomic chunk publication, manifest
coverage/shape checks, and that an incomplete temporary chunk is ignored.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
OUT = HERE / "prototype_output" / "attempt_03"
CHUNK_FRAMES = 5
CORE_FIELDS = (
    "opening_m", "traction_Pa", "damage", "delta_max_m",
    "contact_traction_Pa", "contact_opening_m", "contact_energy_check",
    "interface_work_check", "energy_J", "ledger", "patch_displacement_x_m",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, value: dict) -> None:
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="\n") as stream:
        json.dump(value, stream, indent=2, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def field_schema(core: dict[str, np.ndarray], nodes: np.ndarray) -> dict:
    result = {
        name: {"dtype": array.dtype.str, "tail_shape": list(array.shape[1:])}
        for name, array in core.items()
    }
    result["node_state"] = {"dtype": nodes.dtype.str, "tail_shape": list(nodes.shape[1:])}
    return result


def payload_for_chunk(
    core: dict[str, np.ndarray], node_steps: np.ndarray, node_state: np.ndarray,
    start: int, stop: int,
) -> dict[str, np.ndarray]:
    local = (node_steps >= start) & (node_steps < stop)
    payload = {name: np.asarray(array[start:stop]) for name, array in core.items()}
    payload["node_saved_steps"] = np.asarray(node_steps[local], dtype=np.int64)
    payload["node_state"] = np.asarray(node_state[local])
    return payload


def write_chunk(path: Path, payload: dict[str, np.ndarray]) -> None:
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("wb") as stream:
        np.savez_compressed(stream, **payload)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def read_committed(root: Path, *, allow_incomplete: bool) -> tuple[dict, dict[str, np.ndarray], np.ndarray, np.ndarray]:
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    if manifest.get("schema") != "pgfrag.lossless_history_chunks.prototype.v1":
        raise ValueError("unexpected manifest schema")
    if manifest.get("state") != "COMPLETE" and not allow_incomplete:
        raise ValueError("history chunks are not complete")
    cursor = 0
    cores: dict[str, list[np.ndarray]] = {name: [] for name in CORE_FIELDS}
    node_parts: list[np.ndarray] = []
    node_step_parts: list[np.ndarray] = []
    expected_fields = set(CORE_FIELDS) | {"node_state", "node_saved_steps"}
    for entry in manifest["chunks"]:
        start, stop = int(entry["frame_start"]), int(entry["frame_stop"])
        if start != cursor or stop <= start or stop > manifest["total_frames"]:
            raise ValueError("committed chunk frame ranges are not contiguous")
        chunk = root / entry["file"]
        if chunk.stat().st_size != int(entry["bytes"]) or sha256(chunk) != entry["sha256"]:
            raise ValueError("committed chunk hash/size mismatch")
        with np.load(chunk, allow_pickle=False) as archive:
            if set(archive.files) != expected_fields:
                raise ValueError("chunk member set mismatch")
            local_steps = np.asarray(archive["node_saved_steps"])
            local_nodes = np.asarray(archive["node_state"])
            expected_local_steps = np.asarray([
                step for step in manifest["node_saved_steps"] if start <= step < stop
            ], dtype=np.int64)
            if local_steps.dtype != np.dtype("int64") or not np.array_equal(local_steps, expected_local_steps):
                raise ValueError("chunk node-step index mismatch")
            node_spec = manifest["fields"]["node_state"]
            if local_nodes.dtype.str != node_spec["dtype"] or local_nodes.shape != (len(local_steps), *node_spec["tail_shape"]):
                raise ValueError("chunk node-state dtype/shape mismatch")
            for name in CORE_FIELDS:
                values = np.asarray(archive[name])
                spec = manifest["fields"][name]
                if values.dtype.str != spec["dtype"] or values.shape != (stop - start, *spec["tail_shape"]):
                    raise ValueError(f"chunk field dtype/shape mismatch: {name}")
                cores[name].append(values.copy())
            node_step_parts.append(local_steps.copy())
            node_parts.append(local_nodes.copy())
        cursor = stop
    if cursor != int(manifest["committed_frames"]):
        raise ValueError("manifest committed-frame count disagrees with chunk ranges")
    if manifest["state"] == "COMPLETE" and cursor != int(manifest["total_frames"]):
        raise ValueError("complete manifest lacks all expected frames")
    joined = {name: np.concatenate(parts, axis=0) if parts else np.empty((0,)) for name, parts in cores.items()}
    joined_steps = np.concatenate(node_step_parts) if node_step_parts else np.empty((0,), dtype=np.int64)
    joined_nodes = np.concatenate(node_parts, axis=0) if node_parts else np.empty((0,))
    return manifest, joined, joined_steps, joined_nodes


def main() -> None:
    if OUT.exists():
        raise FileExistsError(f"prototype output already exists: {OUT}")
    OUT.mkdir(parents=True)
    frames, qps, patch_nodes, nodes = 13, 5, 9, 7
    core_shapes = {
        "opening_m": (frames, qps, 3), "traction_Pa": (frames, qps, 3),
        "damage": (frames, qps), "delta_max_m": (frames, qps),
        "contact_traction_Pa": (frames, qps, 3), "contact_opening_m": (frames, qps, 3),
        "contact_energy_check": (frames, 3), "interface_work_check": (frames, 4),
        "energy_J": (frames, 5), "ledger": (frames, 28),
        "patch_displacement_x_m": (frames, patch_nodes),
    }
    core = {
        name: (np.arange(int(np.prod(shape)), dtype=np.float64).reshape(shape) + i / 7.0) / 31.0
        for i, (name, shape) in enumerate(core_shapes.items())
    }
    saved_steps = np.asarray([0, 3, 6, 9, 12], dtype=np.int64)
    node_state = (np.arange(len(saved_steps) * nodes * 4 * 3, dtype=np.float64)
                  .reshape(len(saved_steps), nodes, 4, 3) / 17.0)
    fields = field_schema(core, node_state)
    manifest_path = OUT / "manifest.json"
    manifest = {
        "schema": "pgfrag.lossless_history_chunks.prototype.v1",
        "state": "IN_PROGRESS", "total_frames": frames, "committed_frames": 0,
        "chunk_frames": CHUNK_FRAMES, "node_saved_steps": saved_steps.tolist(),
        "fields": fields, "chunks": [],
    }
    atomic_json(manifest_path, manifest)

    starts = list(range(0, frames, CHUNK_FRAMES))
    for chunk_index, start in enumerate(starts[:2]):
        stop = min(start + CHUNK_FRAMES, frames)
        filename = f"chunk_{chunk_index:06d}.npz"
        path = OUT / filename
        write_chunk(path, payload_for_chunk(core, saved_steps, node_state, start, stop))
        manifest["chunks"].append({
            "file": filename, "frame_start": start, "frame_stop": stop,
            "bytes": path.stat().st_size, "sha256": sha256(path),
        })
        manifest["committed_frames"] = stop
        atomic_json(manifest_path, manifest)

    # Simulate a process crash while writing the next chunk. The manifest still
    # commits only the first ten frames, so this partial ZIP must be ignored.
    incomplete = OUT / "chunk_000002.npz.tmp"
    incomplete.write_bytes(b"PK\x03\x04truncated-incomplete-prototype")
    partial_manifest, partial_core, partial_steps, partial_nodes = read_committed(OUT, allow_incomplete=True)
    assert partial_manifest["committed_frames"] == 10
    assert len(partial_core["opening_m"]) == 10 and len(partial_steps) == 4 and len(partial_nodes) == 4
    assert incomplete.exists() and "chunk_000002.npz" not in {x["file"] for x in partial_manifest["chunks"]}
    incomplete.unlink()

    for chunk_index, start in enumerate(starts[2:], start=2):
        stop = min(start + CHUNK_FRAMES, frames)
        filename = f"chunk_{chunk_index:06d}.npz"
        path = OUT / filename
        write_chunk(path, payload_for_chunk(core, saved_steps, node_state, start, stop))
        manifest["chunks"].append({
            "file": filename, "frame_start": start, "frame_stop": stop,
            "bytes": path.stat().st_size, "sha256": sha256(path),
        })
        manifest["committed_frames"] = stop
        manifest["state"] = "COMPLETE" if stop == frames else "IN_PROGRESS"
        atomic_json(manifest_path, manifest)

    complete_manifest, recovered_core, recovered_steps, recovered_nodes = read_committed(OUT, allow_incomplete=False)
    for name in CORE_FIELDS:
        assert np.array_equal(recovered_core[name], core[name]), f"not bitwise-equal: {name}"
    assert np.array_equal(recovered_steps, saved_steps)
    assert np.array_equal(recovered_nodes, node_state)
    assert complete_manifest["state"] == "COMPLETE"
    result = {
        "status": "TINY_LOSSLESS_CHUNK_STORE_PASS",
        "schema": complete_manifest["schema"],
        "numpy_version": np.__version__,
        "synthetic_frames": frames,
        "chunk_frame_ranges": [[x["frame_start"], x["frame_stop"]] for x in complete_manifest["chunks"]],
        "incomplete_temp_ignored_at_committed_prefix": True,
        "committed_prefix_frames_during_crash_simulation": 10,
        "all_core_fields_bitwise_equal_float64": True,
        "node_state_and_step_index_bitwise_equal": True,
        "manifest_shape_dtype_coverage_hash_checks": "PASS",
        "prototype_output_bytes_including_manifest_and_validation": 0,
        "output_cap_bytes": 64 * 1024 * 1024,
        "not_solver_or_physics_validation": True,
    }
    validation_path = OUT / "prototype_validation.json"
    total_bytes = 0
    for _ in range(5):
        result["prototype_output_bytes_including_manifest_and_validation"] = total_bytes
        validation_path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        final_bytes = sum(p.stat().st_size for p in OUT.iterdir() if p.is_file())
        if final_bytes == total_bytes:
            break
        total_bytes = final_bytes
    else:
        raise AssertionError("could not stabilize prototype output-byte count")
    assert final_bytes < 64 * 1024 * 1024
    print(json.dumps(result, separators=(",", ":")))


if __name__ == "__main__":
    main()
