"""Recompute exact logical history costs and a bounded chunk-size design.

Read-only inputs: the frozen zero-step cost forecast, root independent
zero-step readback, and the current PG source/audit/exporter implementations.
The generated report is written only to the matching history_storage_review
data directory. No solver or physical trajectory is executed.
"""
from __future__ import annotations

import hashlib
import io
import json
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[4]
PG = ROOT / "pgfrag"
FORECAST_PATH = PG / "data/continuation_20261016/reference_vessel_preflight/run_attempt_02/history_cost_forecast.json"
ROOT_MATH_PATH = ROOT / "experiments/continuation_20261016/root_preflight_math.json"
RUN_DIR = PG / "experiments/continuation_20261016/reference_vessel_preflight/run_attempt_02/native_preflight/run01"
FROZEN_INPUTS = PG / "experiments/continuation_20261016/reference_vessel_preflight/run_attempt_02/inputs"
SOURCE_PATH = PG / "src/pgfrag/cohesive_vessel.py"
AUDIT_PATH = PG / "experiments/cohesive_vessel_20261015/audit_native_vessel.py"
EXPORT_V5 = PG / "data/continuation_20261015/native_cohesive_fragments_v5/export_native_fragments.py"
EXPORT_V6 = PG / "data/continuation_20261015/native_cohesive_fragments_v6/export_native_fragments.py"
OUT_DIR = PG / "data/continuation_20261016/history_storage_review"
OUT_JSON = OUT_DIR / "history_storage_review.json"
OUT_MD = OUT_DIR / "README.md"


CORE_SHAPES = {
    "opening_m": lambda f, q, p, n: (f, q, 3),
    "traction_Pa": lambda f, q, p, n: (f, q, 3),
    "damage": lambda f, q, p, n: (f, q),
    "delta_max_m": lambda f, q, p, n: (f, q),
    "contact_traction_Pa": lambda f, q, p, n: (f, q, 3),
    "contact_opening_m": lambda f, q, p, n: (f, q, 3),
    "contact_energy_check": lambda f, q, p, n: (f, 3),
    "interface_work_check": lambda f, q, p, n: (f, 4),
    "energy_J": lambda f, q, p, n: (f, 5),
    "ledger": lambda f, q, p, n: (f, 28),
    "patch_displacement_x_m": lambda f, q, p, n: (f, p),
}
FIELD_USAGE = {
    "opening_m": "Every-step QP opening: cohesive/interface-work integration, contact/full-gap reconstruction, final fracture-state extraction.",
    "traction_Pa": "Every-step QP cohesive traction: work integration, full-damage/zero-traction failure gate, post-failure checks, final extraction.",
    "damage": "Every-step QP damage: bounds/irreversibility, first-full-damage timing, damaged-area/fracture-energy history, post-failure gate, final extraction.",
    "delta_max_m": "Every-step QP maximum opening: irreversible history-variable check.",
    "contact_traction_Pa": "Every-step QP contact traction: zero-contact check, failure/post-failure gate, final extraction.",
    "contact_opening_m": "Every-step QP contact opening: independent quadrature contact energy/work and full-gap reconstruction.",
    "contact_energy_check": "Every-step wrapper-generated three-column check against independently rebuilt contact energy; redundant if the audit explicitly recomputes it.",
    "interface_work_check": "Every-step wrapper-generated work decomposition/cross-work check; redundant if recomputed from raw QP opening/traction.",
    "energy_J": "Every-step native kinetic/potential/dissipated/reversible/contact energy: independent energy and fracture balance.",
    "ledger": "Every-step time/load/work/momentum/torque/force/recurrence/energy residual diagnostics (28 float64 columns).",
    "patch_displacement_x_m": "Every-step displacement on 3708 loaded nodes: independent reconstruction of external work.",
    "node_state": "Saved-step nodal displacement, velocity, acceleration, internal force: independent force/momentum/energy/deformation reconstruction; stride 1 also enables strict every-step recurrences and energy balance. Extraction checks saved states after joint failure and the final state.",
    "node_saved_steps": "Exact int64 step map for node_state; must retain every declared stride frame and final frame.",
    "qp_first_full_damage_step": "Small per-QP int64 summary written by native loop and compared with timing independently recomputed from the full damage history.",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def npy_header_bytes(shape: tuple[int, ...], dtype: np.dtype) -> int:
    stream = io.BytesIO()
    np.lib.format.write_array_header_1_0(stream, {
        "descr": np.dtype(dtype).str,
        "fortran_order": False,
        "shape": shape,
    })
    return stream.tell()


def npy_array_bytes(shape: tuple[int, ...], dtype: np.dtype) -> int:
    return int(np.prod(shape, dtype=np.int64)) * np.dtype(dtype).itemsize + npy_header_bytes(shape, dtype)


def npy_saved_steps_bytes(saved_steps: np.ndarray) -> int:
    stream = io.BytesIO()
    np.save(stream, saved_steps, allow_pickle=False)
    return len(stream.getvalue())


def main() -> None:
    forecast = json.loads(FORECAST_PATH.read_text(encoding="utf-8"))
    math_report = json.loads(ROOT_MATH_PATH.read_text(encoding="utf-8"))
    assert forecast["status"] == "COST_FORECAST_FROM_ZERO_STEP_NATIVE_PREFLIGHT"
    assert math_report["status"] == "INDEPENDENT_ZERO_STEP_READBACK_PASS"
    nodes, tets, facets = map(int, forecast["native_counts_nodes_tets_cohesive_facets"])
    qps = int(forecast["native_quadrature_points"])
    steps = int(forecast["stability_and_time"]["selected_steps_for_50us"])
    frames = steps + 1
    patch_nodes = int(forecast["balanced_patch"]["patch_node_count_total"])
    assert [nodes, tets, facets, qps, steps] == [18634, 56139, 541, 1623, 22504]
    assert math_report["solver_steps"] == 0 and math_report["forecast_steps"] == steps
    dtype = np.dtype(np.float64)
    scalar_dtype = np.dtype(np.int64)
    specs = {name: shape_fn(frames, qps, patch_nodes, nodes) for name, shape_fn in CORE_SHAPES.items()}
    specs["node_state"] = (0, nodes, 4, 3)  # completed below for each selected node stride
    per_frame = {
        "qp_float64_doubles": qps * 14,
        "qp_float64_bytes": qps * 14 * dtype.itemsize,
        "energy_and_check_ledger_float64_doubles": 3 + 4 + 5 + 28,
        "energy_and_check_ledger_bytes": (3 + 4 + 5 + 28) * dtype.itemsize,
        "patch_displacement_bytes": patch_nodes * dtype.itemsize,
        "all_non_node_bytes": qps * 14 * dtype.itemsize + (3 + 4 + 5 + 28) * dtype.itemsize + patch_nodes * dtype.itemsize,
        "node_snapshot_bytes": nodes * 4 * 3 * dtype.itemsize,
    }
    assert per_frame["all_non_node_bytes"] == 211760
    static_run_bytes = int(forecast["node_stride_forecasts"]["50"]["actual_zero_step_preflight_output_bytes_excluding_receipt"])
    source_paths = {
        "history_cost_forecast": FORECAST_PATH,
        "root_independent_math": ROOT_MATH_PATH,
        "native_wrapper": SOURCE_PATH,
        "cohesive_mesh_module": PG / "src/pgfrag/cohesive_mesh.py",
        "cohesive_stability_module": PG / "src/pgfrag/cohesive_stability.py",
        "independent_native_audit": AUDIT_PATH,
        "fragment_export_v5": EXPORT_V5,
        "fragment_export_v6": EXPORT_V6,
        "preflight_receipt": RUN_DIR / "receipt.json",
        "executed_frozen_wrapper": FROZEN_INPUTS / "cohesive_vessel.py",
        "executed_frozen_cohesive_mesh": FROZEN_INPUTS / "cohesive_mesh.py",
        "executed_frozen_cohesive_stability": FROZEN_INPUTS / "cohesive_stability.py",
    }
    source_hashes = {name: {"path": path.relative_to(ROOT).as_posix(), "sha256": sha256(path), "bytes": path.stat().st_size}
                     for name, path in source_paths.items()}

    stride_reports = {}
    for stride, raw_current in forecast["node_stride_forecasts"].items():
        stride = int(stride)
        saved = np.unique(np.r_[np.arange(0, frames, stride, dtype=np.int64), steps]).astype(np.int64)
        node_shape = (len(saved), nodes, 4, 3)
        arrays = {**specs, "node_state": node_shape}
        raw_numeric = sum(int(np.prod(shape, dtype=np.int64)) * dtype.itemsize for shape in arrays.values())
        assert raw_numeric == int(raw_current["raw_numeric_array_bytes"])
        full_headers = sum(npy_header_bytes(shape, dtype) for shape in arrays.values())
        assert raw_current["wrapper_raw_history_bytes_forecast"] == raw_numeric + 12 * 256
        saved_steps_file = npy_saved_steps_bytes(saved)

        # Choose each time chunk so its raw QP/diagnostic + nodal snapshot
        # payload stays below the existing 64 MiB output block budget.
        chunk_frames = 250 if stride == 50 else 25
        chunks: list[dict] = []
        chunk_npy_member_bytes = 0
        chunk_header_bytes = 0
        for chunk_index, start in enumerate(range(0, frames, chunk_frames)):
            stop = min(start + chunk_frames, frames)
            local_node_steps = saved[(saved >= start) & (saved < stop)]
            local_node_frames = len(local_node_steps)
            member_rows: dict[str, dict] = {}
            chunk_raw_bytes = 0
            for name, shape_fn in CORE_SHAPES.items():
                shape = shape_fn(stop - start, qps, patch_nodes, nodes)
                size = int(np.prod(shape, dtype=np.int64)) * dtype.itemsize
                header = npy_header_bytes(shape, dtype)
                chunk_raw_bytes += size
                chunk_header_bytes += header
                chunk_npy_member_bytes += size + header
                member_rows[name] = {
                    "path": f"{name}/chunk_{chunk_index:06d}.npy",
                    "shape": list(shape), "dtype": dtype.str, "bytes": size + header,
                    "sha256": "0" * 64,
                }
            node_local_shape = (local_node_frames, nodes, 4, 3)
            node_local_bytes = int(np.prod(node_local_shape, dtype=np.int64)) * dtype.itemsize
            node_header = npy_header_bytes(node_local_shape, dtype)
            chunk_raw_bytes += node_local_bytes
            chunk_header_bytes += node_header
            chunk_npy_member_bytes += node_local_bytes + node_header
            member_rows["node_state"] = {
                "path": f"node_state/chunk_{chunk_index:06d}.npy",
                "shape": list(node_local_shape), "dtype": dtype.str,
                "bytes": node_local_bytes + node_header, "sha256": "0" * 64,
            }
            chunks.append({
                "frame_start": start, "frame_stop": stop,
                "node_saved_steps": local_node_steps.tolist(),
                "numeric_payload_bytes": chunk_raw_bytes,
                "uncompressed_member_bytes": chunk_raw_bytes + sum(
                    row["bytes"] - int(np.prod(row["shape"], dtype=np.int64)) * dtype.itemsize
                    for row in member_rows.values()
                ),
                "arrays": member_rows,
            })

        layout = {
            "schema": "pgfrag.native_cohesive_vessel_history_chunks.v1",
            "state": "COMPLETE", "dtype_contract": "float64 histories; int64 saved step indices",
            "time_frames": frames, "solver_steps": steps, "chunk_frames": chunk_frames,
            "node_stride": stride, "qp_order": "native element-major, three quadrature points per cohesive_3d_6; no averaging",
            "fields": {name: {"dtype": dtype.str, "tail_shape": list(shape[1:])} for name, shape in arrays.items()},
            "node_saved_steps": {"path": "node_saved_steps.npy", "shape": [len(saved)],
                                 "dtype": scalar_dtype.str, "bytes": saved_steps_file, "sha256": "0" * 64},
            "chunks": chunks,
        }
        manifest_bytes = len((json.dumps(layout, indent=2, allow_nan=False) + "\n").encode("utf-8"))
        node_layout_bytes = int(raw_current["layout_json_bytes"])
        # Chunked NPY files keep the same numeric payload and add exact NPY
        # headers per chunk, one saved-step index, layout, and a verifiable
        # manifest. Directory allocation units/metadata are filesystem-specific.
        chunked_logical = chunk_npy_member_bytes + saved_steps_file + node_layout_bytes + manifest_bytes
        chunk_count = len(chunks)
        file_count = chunk_count * 12 + 3
        node_frame_distribution: dict[str, int] = {}
        for row in chunks:
            count = str(len(row["node_saved_steps"]))
            node_frame_distribution[count] = node_frame_distribution.get(count, 0) + 1
        max_chunk_raw = max(row["numeric_payload_bytes"] for row in chunks)
        max_chunk_members = max(row["uncompressed_member_bytes"] for row in chunks)
        assert max_chunk_members < 64 * 1024 * 1024
        current_folder = int(raw_current["full_history_folder_estimate_bytes"])
        current_run = int(raw_current["estimated_run_directory_bytes_before_safety_margin"])
        margin = int(raw_current["wrapper_2mb_safety_margin_bytes"])
        cap = int(raw_current["configured_output_cap_bytes"])
        stride_reports[str(stride)] = {
            "node_stride": stride, "steps": steps, "time_frames": frames,
            "node_saved_frames": len(saved), "node_saved_steps_dtype": scalar_dtype.str,
            "all_history_array_dtype": dtype.str,
            "array_shapes": {name: list(shape) for name, shape in arrays.items()},
            "numeric_payload_bytes_exact": raw_numeric,
            "single_npy_array_header_bytes_exact_numpy_2_5_3": full_headers,
            "current_wrapper_forecast_bytes_256_byte_header_allowance_each": int(raw_current["wrapper_raw_history_bytes_forecast"]),
            "single_array_npy_files_logical_bytes_if_exact_header_used": raw_numeric + full_headers,
            "node_saved_steps_npy_bytes_exact": saved_steps_file,
            "layout_json_bytes_current_contract": node_layout_bytes,
            "current_history_folder_estimate_bytes": current_folder,
            "current_full_run_estimate_before_margin_bytes": current_run,
            "current_run_allocation_with_2mb_plus_20pct_gib_rounded_up": int(raw_current["recommended_output_allocation_gib_rounded_up"]),
            "chunk_design": {
                "format_forecasted": "uncompressed per-field .npy chunk files with atomic manifest commit; exact raw-storage alternative",
                "chunk_frames": chunk_frames,
                "chunk_count": chunk_count,
                "node_frames_per_chunk_count": node_frame_distribution,
                "max_raw_combined_chunk_bytes": max_chunk_raw,
                "max_raw_combined_chunk_mib": max_chunk_raw / (1024 ** 2),
                "max_uncompressed_chunk_members_including_npy_headers_bytes": max_chunk_members,
                "chunk_npy_numeric_payload_plus_exact_headers_bytes": chunk_npy_member_bytes,
                "chunk_npy_header_overhead_bytes": chunk_header_bytes,
                "manifest_json_bytes_for_declared_schema_and_paths": manifest_bytes,
                "node_saved_steps_npy_bytes": saved_steps_file,
                "layout_json_bytes": node_layout_bytes,
                "logical_history_bytes_exact_excluding_filesystem_allocation_units": chunked_logical,
                "total_history_files_including_manifest_and_layout": file_count,
                "history_size_difference_from_current_folder_estimate_bytes": chunked_logical - current_folder,
                "uncompressed_npz_bundle_member_bytes_before_zip_container": chunk_npy_member_bytes,
                "compressed_npz_final_bytes": None,
                "compressed_size_reason": "No full trajectory exists; the zero-step preflight has no evolved history. Lossless deflate output cannot be forecast from mesh counts or synthetic data.",
                "existing_64mib_case_cap_satisfied_by_this_uncompressed_alternative": False,
            },
            "full_run_logical_bytes_with_chunked_npy_before_margin": static_run_bytes + chunked_logical,
            "filesystem_allocation_overhead_not_included": True,
            "existing_output_cap_bytes": cap,
            "current_history_forecast_bytes_over_case_cap": int(raw_current["bytes_over_configured_history_cap"]),
            "still_requires_storage_allocation_until_lossless_compression_is_measured": True,
        }

    report = {
        "schema": "pgfrag.native_cohesive_vessel_history_storage_review.v1",
        "status": "READ_ONLY_STORAGE_AUDIT_COMPLETE",
        "classification": "COST_ESTIMATE_UNCALIBRATED; NO SPECIMEN PHYSICAL DATA; ZERO SOLVER EVOLUTION",
        "solver_evolution_performed": False,
        "runtime_forecast": "NOT_ESTIMATED; no trajectory was run, and lossless compression adds data-dependent CPU/I/O cost.",
        "native_preflight_status": forecast["native_preflight_status"],
        "root_independent_math_status": math_report["status"],
        "root_independent_receipt_sha256": math_report["source_native_receipt_sha256"],
        "frozen_preflight_modules_match_current_sources": {
            "cohesive_vessel": sha256(SOURCE_PATH) == sha256(FROZEN_INPUTS / "cohesive_vessel.py"),
            "cohesive_mesh": sha256(PG / "src/pgfrag/cohesive_mesh.py") == sha256(FROZEN_INPUTS / "cohesive_mesh.py"),
            "cohesive_stability": sha256(PG / "src/pgfrag/cohesive_stability.py") == sha256(FROZEN_INPUTS / "cohesive_stability.py"),
        },
        "mesh_counts_nodes_tets_cohesive_facets": [nodes, tets, facets],
        "quadrature_points": qps, "balanced_patch_nodes": patch_nodes,
        "selected_dt_s": forecast["stability_and_time"]["selected_dt_s"],
        "steps_for_50us": steps, "frames_including_initial": frames,
        "per_frame_bytes": per_frame,
        "history_fields": [
            {"name": name, "dtype": dtype.str,
             "shape_stride_50": stride_reports["50"]["array_shapes"][name],
             "shape_stride_1": stride_reports["1"]["array_shapes"][name],
             "lifetime": "allocated/pre-sized before evolution; written at every solver frame" if name != "node_state" else "pre-sized for declared node_saved_steps; written only at saved stride frames, always including final frame",
             "consumer": FIELD_USAGE[name]}
            for name in (*CORE_SHAPES.keys(), "node_state")
        ] + [
            {"name": "node_saved_steps", "dtype": scalar_dtype.str,
             "shape_stride_50": [stride_reports["50"]["node_saved_frames"]],
             "shape_stride_1": [stride_reports["1"]["node_saved_frames"]],
             "lifetime": "written once before solver evolution; index map for node_state",
             "consumer": FIELD_USAGE["node_saved_steps"]},
            {"name": "qp_first_full_damage_step", "dtype": scalar_dtype.str,
             "shape_stride_50": [qps], "shape_stride_1": [qps],
             "lifetime": "computed in RAM across the trajectory; written as one array in finally",
             "consumer": FIELD_USAGE["qp_first_full_damage_step"]},
        ],
        "consumers_and_compatibility": {
            "independent_audit": {
                "path": AUDIT_PATH.relative_to(ROOT).as_posix(),
                "format_contract": "direct np.load(path, mmap_mode='r', allow_pickle=False) of standard NPY files; exact fixed shapes, finite values, complete completion.json required",
                "all_step_requirements": ["QP damage and delta_max for monotonicity/first failure", "QP cohesive/contact opening and traction for work, contact energy, failure and post-failure checks", "energy_J and 28-column ledger for balances/residuals", "patch displacement every step for independently rebuilt external work"],
                "node_stride_effect": "stride 1 enables strict full-step central-difference and internal/external work identities; stride 50 only permits sampled node-state reconstruction and is explicitly reported as such",
                "streaming_risk": "audit currently uses full-history NumPy expressions including np.diff, trapezoid integration, contact norms and comparisons; an NPZ reader alone would decompress/copy chunks and old audit expressions can still materialize large arrays. A chunk-aware audit iterator/reduction is required for bounded peak RAM.",
            },
            "fragment_extraction": {
                "paths": [EXPORT_V5.relative_to(ROOT).as_posix(), EXPORT_V6.relative_to(ROOT).as_posix()],
                "format_contract": "direct NPY memmaps for node_state, damage, traction, opening, contact_traction, and direct node_saved_steps NPY",
                "must_preserve": ["all-step damage and traction to establish first joint full-damage/zero-traction frame and continued failure gate", "all saved node snapshots after joint failure for gap/separation witness", "final node displacement/velocity and final QP damage/traction/opening/contact state for output extraction"],
                "compatibility_risk": "both exporters require existing individual NPY paths and receipt coverage; compressed chunks need a reader/exporter adaptation or a verified lazy compatibility layer. Reconstructing monolithic arrays would recreate peak-disk pressure.",
            },
            "receipt_and_checkpoint": {
                "current": "receipt.json inventories and hashes every file; audit rejects extra/missing files. Full-size NPY memmaps are pre-sized before step 0 and flushed every flush_stride frames (case stride 250), then checkpoint.json is written; finally flush occurs even on failure.",
                "risk": "a pre-sized NPY tail is zero-filled/unwritten after a crash; checkpoint/result/completion must be used to reject incomplete histories. write_json writes directly, so current checkpoint JSON is not itself an atomic commit record.",
                "chunk_design": "publish each chunk via temp file, flush+fsync, atomic rename, then atomically replace a manifest that records contiguous committed frame range, node steps, exact shapes/dtypes, byte count and SHA-256. Only a COMPLETE manifest with exact expected coverage can pass an independent audit; partial prefixes are diagnostic only and do not imply solver restartability.",
            },
        },
        "storage_forecasts": stride_reports,
        "forecast_scope": {
            "current_history_specs": "11 float64 arrays sampled every solver step plus node_state float64 at the declared stride; excludes the small QP first-full-damage summary and runtime JSON/control files.",
            "excluded_qp_failure_summary": {
                "path": "qp_first_full_damage_step.npy", "shape": [qps], "dtype": scalar_dtype.str,
                "logical_bytes_numpy_2_5_3": npy_array_bytes((qps,), scalar_dtype),
                "consumer": "independent auditor recomputes and compares the per-QP full-damage step",
            },
            "run_budget_note": "The prior 2,000,000-byte safety allowance covers this small summary and JSON/control outputs, but its forecast is not a complete per-file inventory of a completed evolution run. Filesystem allocation units are also excluded.",
        },
        "design_recommendation": {
            "bounded_next_step": "Keep float64 and every QP value at every time step. Prototype chunked NumPy storage behind a versioned manifest; first preserve the current NPY layout as a legacy read contract. The tiny NPZ prototype verifies lossless chunks/commit markers but does not estimate vessel-data compression.",
            "raw_npy_chunk_storage_effect": "recovery/commit granularity and bounded read blocks improve; logical numeric storage remains essentially the same as current full NPY arrays. It is not a fix for the 64 MiB cap.",
            "lossless_compressed_npz_effect": "may reduce bytes, but exact real size and runtime are unknowable before a representative trajectory. Retain raw-size allocation as the safe forecast until the first actual chunk is measured and the entire run can be budgeted.",
            "suggested_initial_chunk_sizes": {"node_stride_50_time_frames": 250, "node_stride_1_time_frames": 25,
                "reason": "keeps combined raw QP/diagnostic and node-state payload below 64 MiB per chunk (about 59.0 MiB and 47.7 MiB respectively), before small headers"},
            "prohibited_savings_for_this_review": ["diagnostic downsampling", "float precision reduction", "new physics parameters", "treating compression ratios from synthetic arrays as real-history estimates"],
            "no_solver_or_physical_qualification": True,
        },
        "tiny_prototype": {
            "script": "pgfrag/experiments/continuation_20261016/history_storage_review/prototype_chunk_store.py",
            "result": "pgfrag/experiments/continuation_20261016/history_storage_review/prototype_output/attempt_03/prototype_validation.json",
            "scope": "synthetic float64 arrays only; no native model or solver evolution",
        },
        "source_hashes": source_hashes,
    }
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    lines = [
        "# Native vessel history storage review",
        "",
        "This read-only audit covers the current history writer, independent audit, fragment exporters, and the actual zero-step full-size preflight. No solver evolution or storage change to PG source was made. Material values remain uncalibrated placeholders; this is a cost/layout review, not specimen qualification.",
        "",
        "The 50 μs forecast contains 22,505 frames, 1,623 interface quadrature points, 3,708 patch nodes, and 18,634 native nodes after seam insertion. All primary histories are float64. Each frame writes 211,760 bytes of QP/work/ledger/patch values; each saved nodal frame adds 1,788,864 bytes.",
        "",
        "| Node stride | Saved node frames | Raw numeric histories | Current wrapper history estimate | Current full-run allocation with margin |",
        "|---:|---:|---:|---:|---:|",
    ]
    for key in ("50", "1"):
        item = stride_reports[key]
        lines.append(f"| {key} | {item['node_saved_frames']:,} | {item['numeric_payload_bytes_exact']:,} B | {item['current_history_folder_estimate_bytes']:,} B | {item['current_run_allocation_with_2mb_plus_20pct_gib_rounded_up']} GiB |")
    lines += [
        "",
        "The raw chunked-NPY alternative preserves every numeric value and adds exact per-chunk NPY headers plus the declared manifest and saved-step index. It improves fail-closed prefix handling and chunked reads, but does not materially shrink the forecast. Compressed NPZ chunks are lossless and the synthetic test passed bitwise equality, shape/dtype/hash checks, and incomplete-temp exclusion; the actual vessel compression ratio is unknown because there are no evolved histories.",
        "",
        "The independent auditor expects individual NPY files and performs several full-history NumPy operations. Both v5/v6 fragment exporters directly memory-map named NPY files. A chunk format therefore needs a versioned audit/export reader and streaming reductions before it can replace that contract. All QP frames remain necessary for failure/work checks; node state remains exactly at the declared stride (stride 1 for strict full-step nodal identities, stride 50 for sampled reconstruction), including the final frame.",
        "",
        "See `history_storage_review.json` for field-by-field shapes, dtype/lifetime/consumer, exact chunk/header forecasts, compatibility constraints, source hashes, and migration gates. The prototype is under `pgfrag/experiments/continuation_20261016/history_storage_review/`.",
        "",
    ]
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps({"status": report["status"], "report": str(OUT_JSON),
                      "history_sizes": {k: stride_reports[k]["current_history_folder_estimate_bytes"] for k in ("50", "1")},
                      "chunk_logical_sizes": {k: stride_reports[k]["chunk_design"]["logical_history_bytes_exact_excluding_filesystem_allocation_units"] for k in ("50", "1")},
                      "chunk_counts": {k: stride_reports[k]["chunk_design"]["chunk_count"] for k in ("50", "1")}}, separators=(",", ":")))


if __name__ == "__main__":
    main()
