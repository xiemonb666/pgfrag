#!/usr/bin/env python3
"""Frozen, root-gated explicit refinement capture for the 96-T4 coupon.

`--prepare` validates hashes, mesh-derived interface areas, budgets, and writes
the protocol. It does not import Akantu or construct a model. Native execution
requires `--launch` plus the exact SHA-256 of the reviewed protocol file.
"""
from __future__ import annotations

import argparse
from collections import deque
import csv
from datetime import datetime, timezone
from decimal import Decimal
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import sys
import time
import traceback

import numpy as np


ROOT = Path(__file__).resolve().parents[4]
EXPERIMENT = Path(__file__).resolve().parent
DATA = ROOT / "pgfrag/data/continuation_20261016/extrinsic_onset_refinement"
CASE_ROOT = ROOT / "relicweave-synlab/tmp/review_20261004/implementation/cohesive_event_probe_20261004/coarse654_20261004_01"
INPUT = CASE_ROOT / "snapshot/input"
MANIFEST_PATH = CASE_ROOT / "snapshot/manifest.json"
OLD_DRIVER = ROOT / "pgfrag/experiments/cohesive_event_probe_20261004/probe_driver.py"
CAPTURE_PREFLIGHT = ROOT / "pgfrag/experiments/continuation_20261016/extrinsic_capture_preflight/capture_preflight.py"
CAPTURE_EVIDENCE = ROOT / "pgfrag/data/continuation_20261016/extrinsic_capture_preflight/extrinsic_capture.json"
INTRINSIC_EVIDENCE = ROOT / "pgfrag/data/continuation_20261016/extrinsic_capture_preflight/intrinsic_capability_capture.json"
SHAPE_PREFLIGHT_SCRIPT = EXPERIMENT / "shape_only_preflight.py"
SHAPE_PREFLIGHT_EVIDENCE = DATA / "shape_only_preflight/shape_only_preflight.json"
SHAPE_PREFLIGHT_ARRAYS = DATA / "shape_only_preflight/zero_step_native_arrays.npz"
SHAPE_PREFLIGHT_STATUS = DATA / "shape_only_preflight/preflight_status_final.json"
RUNTIME_WHEEL = ROOT / "third_party/akantu-runtime/downloads/akantu-5.0.7.post1-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl"
AKANTU_SOURCE = ROOT / "third_party/akantu-runtime/source-v5.0.7"
ROOT_MATH = ROOT / "docs/continuation/20261016/EXTRINSIC_EVENT_MATH.md"
PROTOCOL_PATH = DATA / "frozen_protocol.json"
API_EVIDENCE_PATH = DATA / "wheel_feengine_public_api.json"

LOCKED_HASHES = {
    "source_msh": "b98402c5a7b7af36b8e5e8a2a90c7db689cb1af86a089c8e1b9a0dbf07ba2638",
    "material_dat": "c300c3091ac323eba756f2d7bfbf55e888c68d6b10f5fcb9763c0e6b6194305f",
    "generated_geometry": "eec01312e181581c1d433b2619b6fd8a4a42677769c4dbc25053d0a2c355188b",
    "manifest": "21b563069550e2091474fe94877bd3ceaae094ac4c92aed3900ecdd94702dd68",
    "runtime_wheel": "f618c31bb4cdf23d86e74a10f7e94dac3dae7b0ae2cfd057d352c1b156c0654b",
}
SOURCE_HASHES = {
    "solid_mechanics_model.cc": ("src/model/solid_mechanics/solid_mechanics_model.cc", "fc7873c73f830381fd73b24bf80ad5fb07ca02258ba24d30241632835fcf0465"),
    "solid_mechanics_model_cohesive.cc": ("src/model/solid_mechanics/solid_mechanics_model_cohesive/solid_mechanics_model_cohesive.cc", "b9a3b89ef9b1f17f82fb3d80c99666d292cd5e8615d56d4fa9ced648030f3d64"),
    "material_cohesive.cc": ("src/model/solid_mechanics/solid_mechanics_model_cohesive/materials/material_cohesive.cc", "46bbb9639a975c04e1ca136bb2ecd26efccb80b21294ea98ade9aae9ac849b13"),
    "material_cohesive_linear.cc": ("src/model/solid_mechanics/solid_mechanics_model_cohesive/materials/constitutive_laws/material_cohesive_linear.cc", "90dcef30fc1f60ad668f33a3756875d0bdeb87f4a354d157674051eb8494fda2"),
    "material_cohesive_linear.hh": ("src/model/solid_mechanics/solid_mechanics_model_cohesive/materials/constitutive_laws/material_cohesive_linear.hh", "597e59fd02ca1853f1c861785684e69b9aa071dd380ab50c0a4c9d917e353bbc"),
}

DT_H = Decimal("5.8582308142940834e-8")
BASE_STEPS = 666
OUTPUT_CAP_BYTES = 256 * 1024 * 1024
WALL_CAP_PER_RUN_S = 300.0
WALL_CAP_CAMPAIGN_S = 20 * 60.0
MAX_NODES = 45 + 8 * 3  # 45 source nodes plus the conservative 3 new nodes per eligible Tri3.
MAX_COHESIVE_ELEMENTS = 8
QPS_PER_COHESIVE_ELEMENT = 3
MAX_COHESIVE_QPS = MAX_COHESIVE_ELEMENTS * QPS_PER_COHESIVE_ELEMENT
D_CHUNK_STEPS = 128
ROOT_REVIEW_REQUIRED = True

ENERGY_KEYS = {
    "kinetic_J": "kinetic",
    "bulk_potential_J": "potential",
    "cohesive_dissipated_J": "dissipated",
    "cohesive_reversible_J": "reversible",
    "cohesive_contact_native_key_J": "cohesive contact",
    "contact_alias_diagnostic_J": "contact",
}

COHESIVE_FIELDS = (
    ("sigma_c_eff", 1),
    ("delta_c_eff", 1),
    ("insertion_stress", 3),
    ("normal", 3),
    ("delta max", 1),
    ("damage", 1),
    ("opening", 3),
    ("tractions", 3),
    ("contact_tractions", 3),
    ("contact_opening", 3),
    ("reversible_energy", 1),
    ("total_energy", 1),
)

NODAL_ARRAY_NAMES = (
    "nodal_displacement",
    "nodal_velocity",
    "nodal_acceleration",
    "nodal_mass",
    "nodal_internal_force",
    "nodal_external_force",
    "nodal_blocked_dofs",
)
BARRIER_ARRAY_NAMES = (
    "nodes_rest",
    "bulk_tetra_connectivity",
    "cohesive_connectivity",
) + NODAL_ARRAY_NAMES
CHUNK_ARRAY_NAMES = NODAL_ARRAY_NAMES + ("bulk_stress_native",) + tuple(
    "cohesive_field__" + name.replace(" ", "_") for name, _ in COHESIVE_FIELDS
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relative(path: Path) -> str:
    return path.resolve().relative_to(ROOT.resolve()).as_posix()


def file_record(path: Path, expected: str | None = None) -> dict[str, object]:
    record: dict[str, object] = {"path": relative(path), "exists": path.is_file()}
    if path.is_file():
        record.update(byte_count=path.stat().st_size, sha256=sha256_file(path))
    else:
        record.update(byte_count=None, sha256=None)
    if expected is not None:
        record["expected_sha256"] = expected
        record["hash_matches_expected"] = record["sha256"] == expected
    return record


def jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        value = float(value)
    if isinstance(value, float) and not math.isfinite(value):
        return repr(value)
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    if isinstance(value, dict):
        return {str(key): jsonable(item) for key, item in value.items()}
    return repr(value)


def write_json(path: Path, value, *, exclusive: bool = False) -> None:
    payload = json.dumps(jsonable(value), indent=2, sort_keys=True) + "\n"
    if exclusive:
        try:
            with path.open("x", encoding="utf-8", newline="\n") as stream:
                stream.write(payload)
        except FileExistsError:
            if path.read_text(encoding="utf-8") != payload:
                raise RuntimeError(f"refusing to overwrite differing frozen file: {relative(path)}")
        return
    path.write_text(payload, encoding="utf-8", newline="\n")


def existing_bytes(root: Path) -> int:
    if not root.exists():
        return 0
    return sum(path.stat().st_size for path in root.rglob("*") if path.is_file())


def validate_hashes() -> dict[str, object]:
    frozen_paths = {
        "source_msh": INPUT / "block_2x2x4mm.msh",
        "material_dat": INPUT / "material.dat",
        "generated_geometry": INPUT / "generated_geometry.json",
        "manifest": MANIFEST_PATH,
        "runtime_wheel": RUNTIME_WHEEL,
    }
    records = {name: file_record(path, LOCKED_HASHES[name]) for name, path in frozen_paths.items()}
    mismatches = [name for name, item in records.items() if not item["exists"] or not item["hash_matches_expected"]]
    if mismatches:
        raise RuntimeError(f"frozen input/runtime hash mismatch: {mismatches}")
    for name, (relpath, expected) in SOURCE_HASHES.items():
        item = file_record(AKANTU_SOURCE / relpath, expected)
        records["cached_source_" + name] = item
        if not item["exists"] or not item["hash_matches_expected"]:
            raise RuntimeError(f"cached native source hash mismatch: {name}")

    ancillary = {
        "old_probe_driver": OLD_DRIVER,
        "root_math_derivation": ROOT_MATH,
        "capture_preflight_driver": CAPTURE_PREFLIGHT,
        "extrinsic_capture_preflight_result": CAPTURE_EVIDENCE,
        "intrinsic_capability_preflight_result": INTRINSIC_EVIDENCE,
        "shape_only_preflight_driver": SHAPE_PREFLIGHT_SCRIPT,
        "shape_only_preflight_result": SHAPE_PREFLIGHT_EVIDENCE,
        "shape_only_preflight_raw_arrays": SHAPE_PREFLIGHT_ARRAYS,
        "shape_only_preflight_status": SHAPE_PREFLIGHT_STATUS,
        "facet_ancestry_module": ROOT / "pgfrag/src/pgfrag/facet_ancestry.py",
        "cohesive_mesh_reader": ROOT / "pgfrag/src/pgfrag/cohesive_mesh.py",
        "wheel_feengine_api_evidence": API_EVIDENCE_PATH,
    }
    for name, path in ancillary.items():
        if not path.is_file():
            raise FileNotFoundError(f"required provenance file is absent: {relative(path)}")
        records[name] = file_record(path)
    if not API_EVIDENCE_PATH.is_file():
        raise RuntimeError("wheel FEEngine API inventory is missing; run --inspect-wheel-api before --prepare")
    api = json.loads(API_EVIDENCE_PATH.read_text(encoding="utf-8"))
    if api.get("wheel_sha256") != LOCKED_HASHES["runtime_wheel"]:
        raise RuntimeError("FEEngine API evidence belongs to a different wheel")
    if api.get("model_initialized") is not False or api.get("solver_steps") != 0:
        raise RuntimeError("FEEngine API inventory must be module-level only, with no model or solver steps")
    shape_probe = json.loads(SHAPE_PREFLIGHT_EVIDENCE.read_text(encoding="utf-8"))
    shape_status = json.loads(SHAPE_PREFLIGHT_STATUS.read_text(encoding="utf-8"))
    if shape_probe.get("status") != "VERIFIED_COUPON_AFFINE_TRI3":
        raise RuntimeError("zero-step shape preflight did not verify the coupon-specific Tri3 measure")
    if shape_status.get("status") != "COMPLETED_ZERO_STEP_SHAPE_ONLY" or any(
        shape_status.get(name) != 0 for name in ("solver_steps", "cohesive_checks", "observation_assemblies")
    ):
        raise RuntimeError("zero-step shape preflight status/counters are not clean")
    if shape_probe.get("frozen_input_hashes", {}).get("wheel_sha256") != LOCKED_HASHES["runtime_wheel"]:
        raise RuntimeError("zero-step shape preflight belongs to a different runtime wheel")
    if shape_probe.get("public_api_layout", {}).get("FEEngine.getShapes", {}).get("layout", {}).get("layout") != "FLAT_ELEMENT_MAJOR_REPEATED_QP_NODE":
        raise RuntimeError("native Tri3 shape layout differs from the frozen verified preflight layout")
    return records


def load_mesh_and_area_report() -> dict[str, object]:
    sys.path.insert(0, str(ROOT / "pgfrag/src"))
    from pgfrag.cohesive_mesh import read_gmsh22

    mesh = read_gmsh22(INPUT / "block_2x2x4mm.msh")
    interface = [triangle for triangle in mesh["triangles"] if triangle["physical_name"] == "interface"]
    areas: list[float] = []
    for triangle in interface:
        ids = [int(node_id) for node_id in triangle["node_tags"]]
        id_to_row = {int(node_id): row for row, node_id in enumerate(mesh["node_tags"])}
        points = np.asarray([mesh["nodes"][id_to_row[node_id]] for node_id in ids], dtype=np.float64)
        areas.append(float(np.linalg.norm(np.cross(points[1] - points[0], points[2] - points[0])) * 0.5))
    if len(interface) != 8 or len(areas) != 8:
        raise RuntimeError(f"expected 8 saved physical interface Tri3 elements; found {len(interface)}")
    if any(not math.isclose(area, 0.5e-6, rel_tol=0.0, abs_tol=1e-18) for area in areas):
        raise RuntimeError(f"actual saved interface Tri3 areas differ from frozen geometry: {areas}")
    return {
        "source": "recomputed from frozen Gmsh 2.2 node coordinates and physical interface triangles",
        "triangle_count": len(areas),
        "triangle_areas_m2": areas,
        "total_area_m2": math.fsum(areas),
        "method": "0.5 * norm(cross(x1-x0, x2-x0)); no nominal rectangle area substituted",
    }


def run_matrix() -> list[dict[str, object]]:
    result = []
    for divisor, name in ((1, "h"), (2, "h_over_2"), (4, "h_over_4"), (8, "h_over_8_reserved")):
        dt = DT_H / Decimal(divisor)
        steps = BASE_STEPS * divisor
        result.append({
            "resolution": name,
            "divisor": divisor,
            "dt_decimal_s": str(dt),
            "dt_binary64_s": float(dt),
            "steps_to_common_endpoint": steps,
            "common_endpoint_decimal_s": str(DT_H * BASE_STEPS),
            "status": "PLANNED_ROOT_REVIEW" if divisor <= 4 else "RESERVED_REQUIRES_SEPARATE_ROOT_CONFIRMATION",
        })
    return result


def output_forecast() -> dict[str, object]:
    # D chunk payload: six nodal float arrays, blocked mask, bulk stress, and
    # all twelve cohesive internal fields. Topology/rest coordinates are
    # captured and ancestry-checked every D in memory, but persisted at the
    # initial state and event phases only; per-D connectivity hashes are kept.
    nodal_float_arrays = 6 * MAX_NODES * 3 * 8
    blocked = MAX_NODES * 3
    bulk_stress = 96 * 9 * 8
    cohesive_vectors = 6 * MAX_COHESIVE_QPS * 3 * 8
    cohesive_scalars = 6 * MAX_COHESIVE_QPS * 8
    d_row_bytes = nodal_float_arrays + blocked + bulk_stress + cohesive_vectors + cohesive_scalars
    def estimate(divisors: tuple[int, ...], scope: str) -> dict[str, object]:
        d_steps = sum(BASE_STEPS * divisor for divisor in divisors)
        d_chunk_count = sum(math.ceil((BASE_STEPS * divisor) / D_CHUNK_STEPS) for divisor in divisors)
        d_payload = d_row_bytes * d_steps
        d_container_allowance = d_chunk_count * 64 * 1024
        scalar_trace_allowance = d_steps * 1024
        # A/B/C pre-energy snapshots, one C post-getter snapshot, one D
        # post-getter snapshot, up to four preceding B frames, and four
        # terminal B frames. D pre-energy arrays are already in dense chunks.
        event_archive_count = (8 * 9 + 4) * len(divisors)
        dense_event_allowance = event_archive_count * 64 * 1024
        lineage_and_phase_allowance = d_steps * 512
        logs_status_and_protocol = len(divisors) * 2 * 1024 * 1024
        total = d_payload + d_container_allowance + scalar_trace_allowance + dense_event_allowance + lineage_and_phase_allowance + logs_status_and_protocol
        return {
            "scope": scope,
            "resolution_divisors": list(divisors),
            "d_snapshots": d_steps,
            "d_chunk_count_at_128_steps": d_chunk_count,
            "d_array_payload_bytes": d_payload,
            "d_chunk_container_allowance_bytes": d_container_allowance,
            "scalar_trace_allowance_bytes": scalar_trace_allowance,
            "event_archive_allowance_bytes": dense_event_allowance,
            "ancestry_phase_metadata_allowance_bytes": lineage_and_phase_allowance,
            "logs_status_and_protocol_allowance_bytes": logs_status_and_protocol,
            "total_forecast_bytes": total,
            "total_forecast_mib": total / (1024 * 1024),
            "under_cap": total <= OUTPUT_CAP_BYTES,
        }

    h_forecast = estimate((1,), "single baseline h run upper bound")
    requested_forecast = estimate((1, 2, 4), "currently requested h/h2/h4 envelope; h/8 remains separately root-gated")
    four_run_forecast = estimate((1, 2, 4, 8), "maximum four-resolution envelope including separately gated h/8")
    per_resolution = {
        name: estimate((divisor,), f"single {name} run upper bound")
        for divisor, name in ((1, "h"), (2, "h_over_2"), (4, "h_over_4"), (8, "h_over_8_reserved"))
    }
    return {
        "output_cap_bytes": OUTPUT_CAP_BYTES,
        "output_cap_mib": OUTPUT_CAP_BYTES / (1024 * 1024),
        "worst_case_node_rows": MAX_NODES,
        "worst_case_cohesive_elements": MAX_COHESIVE_ELEMENTS,
        "worst_case_cohesive_qps": MAX_COHESIVE_QPS,
        "d_snapshot_uncompressed_array_bytes_per_step": d_row_bytes,
        "forecast_method": "uncompressed D array payload plus explicit conservative fixed allowances; no compression ratio assumed",
        "baseline_h_run": h_forecast,
        "requested_h_h2_h4": requested_forecast,
        "all_four_including_reserved_h8": four_run_forecast,
        "per_resolution": per_resolution,
        # Compatibility summary: top-level totals continue to denote the
        # maximum envelope, while each currently planned scope is explicit.
        "total_forecast_bytes": four_run_forecast["total_forecast_bytes"],
        "total_forecast_mib": four_run_forecast["total_forecast_mib"],
        "under_cap": four_run_forecast["under_cap"],
    }


def build_protocol() -> dict[str, object]:
    hashes = validate_hashes()
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    areas = load_mesh_and_area_report()
    forecast = output_forecast()
    if not forecast["under_cap"]:
        raise RuntimeError("four-resolution output forecast exceeds 256 MiB; campaign is not launchable")
    if manifest["geometry"]["bulk_tetrahedra"] != 96 or manifest["expected_insertion_facets"] != 8:
        raise RuntimeError("frozen manifest is not the accepted 96-T4/eight-interface case")
    return {
        "schema": "pgfrag.extrinsic_onset_refinement.protocol.v1",
        "status": "FROZEN_FOR_ROOT_REVIEW_NO_CAMPAIGN_LAUNCHED",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "driver": file_record(Path(__file__)),
        "hash_inventory": hashes,
        "runtime": {
            "distribution": "akantu 5.0.7.post1",
            "wheel_sha256": LOCKED_HASHES["runtime_wheel"],
            "cached_native_source_tag": "v5.0.7",
            "cached_source_is_exact_wheel_source": False,
            "python_executable_expected": "/mnt/c/Users/xiemo/Desktop/ShardForge/third_party/akantu-runtime/venv/bin/python",
            "thread_limits": {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"},
            "native_model_initialized_during_prepare": False,
            "solver_steps_during_prepare": 0,
        },
        "frozen_model": {
            "mesh_elements": 96,
            "preinsert_nodes": 45,
            "eligible_physical_surface": "interface",
            "eligible_interface_facets": 8,
            "initial_cohesive_elements": 0,
            "analysis_method": "_explicit_lumped_mass",
            "insertion_mode": "extrinsic",
            "material_input": "frozen material.dat; no edits",
            "material_values_from_frozen_manifest": manifest["material"],
            "actual_interface_geometry": areas,
            "root_derivation": file_record(ROOT_MATH),
        },
        "time_matrix": run_matrix(),
        "caps": {
            "max_solver_steps_total_if_all_four_runs_authorized": sum(BASE_STEPS * divisor for divisor in (1, 2, 4, 8)),
            "currently_requested_steps_h_h2_h4": sum(BASE_STEPS * divisor for divisor in (1, 2, 4)),
            "max_steps_per_run": BASE_STEPS * 8,
            "max_wall_seconds_per_run": WALL_CAP_PER_RUN_S,
            "max_campaign_wall_seconds": WALL_CAP_CAMPAIGN_S,
            "max_output_bytes_total": OUTPUT_CAP_BYTES,
            "one_attempt_per_resolution": True,
            "automatic_retry": False,
            "h_over_8_requires_separate_root_confirmation": True,
        },
        "output_forecast": forecast,
        "native_call_order": [
            "parseInput(frozen material.dat)",
            "read frozen MSH; construct SolidMechanicsModelCohesive",
            "initFull(_analysis_method=_explicit_lumped_mass, _is_extrinsic=True)",
            "assembleMassLumped() once at initial setup",
            "addPhysicalSurface('interface'); updateAutomaticInsertion(); require 8 eligible and 0 inserted",
            "apply zero-valued bottom z Dirichlet condition once",
            "setTimeStep(dt)",
            "per step A getter capture; set ramp time=step*dt; applyBC(ramp,'top'); B getter capture",
            "checkCohesiveStress(); capture return and native cohesive count delta",
            "if insertion occurred, refetch/copy all C arrays, then call the legacy native energy getters in their original position; cached v5.0.7 kinetic getter may internally reassemble lumped mass after the C array copy",
            "after C energy getters, copy post-getter nodal/topology arrays without an explicit assembly and compare them with C pre-energy arrays",
            "solveStep(); immediately refetch/copy D arrays, call the legacy energy getters, and copy post-getter nodal/topology arrays on insertion steps; capture every D through the common endpoint",
        ],
        "phase_clock": {
            "native_getTime": "unavailable in actual wheel; never called or presented as measured",
            "driver_time": "exact Decimal(step)*Decimal(dt), stored with binary64 dt and scheduled BC value",
            "phase_A": "before the current step's applyBC; native state is copied immediately at that phase",
            "phase_B": "after applyBC(step*dt), before checkCohesiveStress; bulk/facet stress are existing stress values",
            "phase_C": "immediately after checkCohesiveStress and any insertion, before solveStep; C_pre_energy arrays are copied before the native energy getters; the subsequent kinetic getter may be an internal mass-assembly barrier",
            "phase_D": "immediately after solveStep; D_pre_energy arrays precede the native energy getters; insertion-step post-getter arrays are copied without requesting another energy or force/mass assembly",
        },
        "capture_plan": {
            "every_A_and_B": "capture state to memory; retain current A and a rolling four previous no-insertion B states; persist dense A/B/C only on a reported insertion step",
            "every_C": "capture only after a positive insertion delta; C_pre_energy stores appended node/mass/force values before energy getters; C_post_energy_getter separately copies nodal/topology arrays after the legacy energy sequence, allowing observation of a possible kinetic-getter mass barrier without an explicit assembly",
            "every_D": "capture all nodal displacement, velocity, acceleration, mass, internal/external force and blocked DOFs; bulk material stress; all cohesive fields; energies; scalar trace; snapshot-local ancestry check; D arrays are pre-energy-getter and D post-getter nodal/topology arrays are saved on insertion steps",
            "after_insertion": "every native view is fetched again and copied after the operation; no pre-insertion AkaArray view is reused",
            "native_cohesive_fields": [name for name, _ in COHESIVE_FIELDS],
            "bulk_stress": "MaterialLinearElastic3D.getStress(_tetrahedron_4,_not_ghost), copied at each A/B/C/D sample",
            "native_energies": ENERGY_KEYS,
            "normals": "capture FacetsFEEngine native normals on pre-check candidate slots at B; capture cohesive internal `normal` after insertion; never infer cohesive orientation from MSH winding",
            "preinsert_candidate_stress": "capture only check-mask candidate rows from getStressOnFacets at B; preserve native slot and QP shape; no cross-insertion QP ordinal join",
            "connectivity_and_ancestry": "capture rest nodes/bulk T4/cohesive connectivity in event snapshots; hash/count and run generic exact source-cell-row/local-corner ancestry on every sample; declare unchanged bulk-row assumption",
            "candidate_facet_lineage": "public Mesh.getMeshFacets().getConnectivity(_triangle_3) rows are mapped through exact rest-node row identity to source Gmsh node IDs and stable source facet FIDs; the initial check mask must equal the eight frozen interface FIDs; live per-step candidates may only be subsets",
            "shape_maps": "copy read-only public FEEngine getShapes and integration-point counts; preserve all native facet slots and within-triangle QP order; no cross-insertion QP ordinal join",
            "coupon_specific_qp_measure": "use VERIFIED_COUPON_AFFINE_TRI3 only if exact native facet connectivity/FIDs and coordinates, three symmetric Tri3 shape points, and native unit normals aligned with the exact triangle geometry all pass; then each of the eight original affine coupon triangles uses its measured coordinate area divided by three per QP; otherwise retain raw arrays and report quadrature work UNKNOWN; do not generalize this rule",
            "assembly_policy": {
                "explicit_initial_assembleMassLumped_calls": 1,
                "extra_explicit_observation_assembly_calls": 0,
                "cached_v5_0_7_kinetic_getter": "getEnergy('kinetic') routes through assembleLumpedMatrix('M') to assembleMassLumped(); onNodesAdded marks mass for reassembly; treat C energy read as a possible MASS_ASSEMBLY_BARRIER after C_pre_energy arrays",
                "installed_5_0_7_post1_wheel_equivalence_to_cached_v5_0_7_source": "UNKNOWN; do not claim wheel implementation identity",
                "C_post_energy_getter_copy": "public getters only; no extra energy, force, or mass assembly call",
            },
        },
        "provisional_numerical_screens": {
            "matched_history_difference": {
                "limit": 0.02,
                "definition": "for each frozen scalar observable and each of h-vs-h/2, h/2-vs-h/4, and h-vs-h/4, take max over exact common times of abs(x_a-x_b)/scale; compare h step k to h/2 step 2k and h/4 step 4k for k=1..666; missing/nonfinite observables or failed identity mapping yield UNKNOWN; FAIL if any complete metric exceeds 0.02; PASS requires every mandatory metric complete and within threshold",
                "observable_scales": {
                    "bc_setpoint_m": "delta_c=2.1714285714285712e-7 m",
                    "internal_force_l2_N": "sigma_c*A=280 N",
                    "external_force_l2_N": "sigma_c*A=280 N",
                    "top_internal_reaction_estimate_N": "sigma_c*A=280 N; remains an estimate",
                    "cohesive_opening_l2_m": "delta_c=2.1714285714285712e-7 m",
                    "cohesive_traction_l2_Pa": "sigma_c=70 MPa",
                    "bulk_stress_max_abs_Pa": "sigma_c=70 MPa",
                    "mass_x_kg": "rho*V=4.0e-5 kg",
                    "mass_y_kg": "rho*V=4.0e-5 kg",
                    "mass_z_kg": "rho*V=4.0e-5 kg",
                    "kinetic_J": "Gc*A=3.04e-5 J",
                    "bulk_potential_J": "Gc*A=3.04e-5 J",
                    "cohesive_dissipated_J": "Gc*A=3.04e-5 J",
                    "cohesive_reversible_J": "Gc*A=3.04e-5 J",
                    "cohesive_contact_native_key_J": "Gc*A=3.04e-5 J",
                },
                "status": "numerical_screen_not_physical_calibration",
            },
            "energy_defect": {"formula": "max_t abs(K + U_bulk + D_cohesive + R_cohesive + E_contact - E_initial - W_boundary)", "limit_fraction_of_GcA": 0.01, "scale_J": 3.04e-5, "threshold_J": 3.04e-7, "status": "UNKNOWN unless actual actuator force, prescribed kinematic phase, and native energy ledger are independently verified"},
            "force_scale_N": 280.0,
            "opening_scale_m": 2.1714285714285712e-7,
            "cohesive_strength_scale_Pa": 70000000.0,
            "mass_scale_kg": 4.0e-5,
            "GcA_scale_J": 3.04e-5,
            "common_endpoint_s": str(DT_H * BASE_STEPS),
            "all_screens": "numerical screens only; no physical calibration, material qualification, or onset CFL certificate",
        },
        "candidate_facet_lineage_and_measure": {
            "stable_identity": "exact source node IDs from unchanged initial T4 row/local-corner order; FID from pgfrag.facet_ancestry source table anchored by frozen MSH SHA-256",
            "native_slot_source": "public Mesh.getMeshFacets().getConnectivity(_triangle_3) and the same triangle_3 ElementTypeMap ordering used by the insertion check mask, FacetsFEEngine normals/stress, and getShapes",
            "runtime_check": "initial native check mask must map exactly to all eight frozen interface FIDs; each later candidate mask must be a subset of that FID set; exact facet-mesh node coordinates must equal original native/source node rows",
            "coupon_qp_measure": "VERIFIED_COUPON_AFFINE_TRI3 conditionally; exact actual facet coordinates/connectivity determine each area, FEEngine exposes three symmetric shape points per affine Tri3, and public native normals must be unit and parallel to the exact connectivity cross-product with one signed orientation across QPs; then A/3 is recorded in observed native QP order",
            "scope_limit": "only the eight planar affine Tri3 facets in this exact frozen coupon; no claim for curved facets, other element types, or arbitrary native QP ordinals",
            "failure_semantics": "retain raw vectors and mark signed work/area-weighted quadrature diagnostics UNKNOWN; never invent a slot map, clamp, or take an absolute energy",
        },
        "ledger_and_native_limitations": {
            "prescribed_boundary_work": "capture top displacement history, blocked DOFs, native internal/external forces and their phases; the top internal-reaction sign estimate remains an estimate, not a verified actuator force",
            "energy_ledger_pass": "forbidden unless actual actuator force and prescribed kinematic phase are proven; otherwise report exact missing quantity and leave status UNKNOWN",
            "mass": "C_pre_energy arrays are captured before native energy getters and may be stale; cached source shows the kinetic getter can reassemble mass after insertion, so C_post_energy_getter arrays are recorded separately. D arrays are captured after solveStep and before energy getters; insertion-step D post-getter arrays test for changes. No energy getter is called a second time for these comparisons.",
            "energy_getter_side_effects": {
                "cached_native_source": "solid_mechanics_model.cc getKineticEnergy assembles the lumped mass; onNodesAdded marks it for reassembly. Cohesive energy getters integrate already-stored reversible/dissipated/contact fields without observed constitutive-history writes. Potential getter may update an energy cache from stored stress/grad_u but does not recompute mechanical stress in the inspected source.",
                "wheel_source_identity": "UNKNOWN; installed wheel is 5.0.7.post1 and cached source is tagged v5.0.7",
                "capture_interpretation": "native_state preserves old array-before-energy call order; C and insertion-step D post-getter public array copies are explicitly separate observations",
            },
            "weights": "no public per-QP integration-weight/Jacobian getter exists; the only permitted quadrature measure is the conditionally verified area/3 rule for these eight affine coupon Tri3 facets",
            "native_insert_stress": "capture native per-QP sigma_c_eff, delta_c_eff, insertion_stress, normal, opening, traction, contact opening/traction and damage; no sigma inferred from traction",
        },
        "comparison_outline": [
            "compare first insertion physical clock step*dt and stable facet identity set across h/h2/h4",
            "compute q_first only where native cohesive opening, actual native normal, and frozen beta=3 metric are populated; report q_first/h^2; retain missing evidence as UNKNOWN",
            "anchor any local onset ledger at the last verified pre-event D state inside the linear ramp; evaluate K_phys=0.5*sum_free(m*v^2)+0.5*sum_blocked(m*(g')^2), report native kinetic separately, and retain unresolved startup impulse; use R_b=m_b*g''-Fint_b-Fext_b only at a verified force phase (constant-velocity ramp with ext=0 reduces to -Fint_b); trapezoidal boundary work remains an estimate unless actuator force and phase are verified",
            "evaluate direction compatibility against native insertion-stress vector only after its field semantics are confirmed from the cached native law; do not create a circular reconstruction",
            "reconstruct signed QP work only with snapshot-local stable FIDs, validated native opening/traction semantics, and VERIFIED_COUPON_AFFINE_TRI3; use ΔW=sum((A/3)*0.5*(t_prev+t_cur)·(δ_cur-δ_prev)), R=sum((A/3)*0.5*t_cur·δ_cur), D=ΔW-R, preserving every signed increment",
            "compare force, opening, history, native mass, and energy at matched physical times; never compare row/QP ordinals across insertion transitions",
            "report the first negative signed-QP-work FID/QP and the complete captured vector path; treat under-resolution as one contributing mechanism, not a sole-cause conclusion",
            "report frozen post-insertion tangent/eigenvalue diagnostics separately; they are not onset CFL certification",
        ],
        "launch_gate": {
            "requires_reviewed_protocol_sha256": True,
            "required_flag": "--launch",
            "campaign_was_launched_by_prepare": False,
            "exact_commands": "README.md and prepare_validation.json carry the three currently planned h, h/2, h/4 commands and protocol SHA-256; a separate h/8 root-confirmation gate remains; no campaign was launched by preparation",
        },
    }


def prepare() -> int:
    DATA.mkdir(parents=True, exist_ok=True)
    protocol = build_protocol()
    write_json(PROTOCOL_PATH, protocol, exclusive=True)
    protocol_hash = sha256_file(PROTOCOL_PATH)
    if existing_bytes(DATA) > OUTPUT_CAP_BYTES:
        raise RuntimeError("prepared data folder already exceeds campaign output cap")
    commands = launch_commands(protocol_hash)
    summary = {
        "status": "PREPARED_NO_MODEL_NO_SOLVER_STEPS",
        "protocol_path": relative(PROTOCOL_PATH),
        "protocol_sha256": protocol_hash,
        "expected_campaign_solver_steps_h_h2_h4": sum(BASE_STEPS * divisor for divisor in (1, 2, 4)),
        "campaign_launched": False,
        "commands_for_root_review": commands,
        "output_forecast": protocol["output_forecast"],
    }
    write_json(DATA / "prepare_validation.json", summary, exclusive=True)
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


def launch_commands(protocol_hash: str) -> list[str]:
    linux_root = "/mnt/c/Users/xiemo/Desktop/ShardForge"
    python = linux_root + "/third_party/akantu-runtime/venv/bin/python"
    script = linux_root + "/pgfrag/experiments/continuation_20261016/extrinsic_onset_refinement/refinement_driver.py"
    common = f"wsl.exe --exec env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 {python} {script} --launch --reviewed-protocol-sha256 {protocol_hash}"
    return [
        common + " --resolution h",
        common + " --resolution h_over_2",
        common + " --resolution h_over_4",
    ]


def inspect_wheel_api() -> int:
    """Record module-level Python binding surface without constructing a model."""
    import akantu as aka

    expected_wheel = sha256_file(RUNTIME_WHEEL)
    if expected_wheel != LOCKED_HASHES["runtime_wheel"]:
        raise RuntimeError("installed API inventory requested against a wheel with wrong SHA-256")
    names = (
        "getShapes", "getShapesDerivatives", "getNbIntegrationPoints",
        "getNormalsOnIntegrationPoints", "getIntegrationPoints", "getWeights",
        "getIntegrationPointWeights", "getJacobians", "getJacobian",
        "integrate", "integrateOnIntegrationPoints", "getShapeFunctionsInterface",
        "getIntegratorInterface",
    )
    surface = {}
    for name in names:
        member = getattr(aka.FEEngine, name, None)
        doc = getattr(member, "__doc__", None) if member is not None else None
        if doc:
            doc = re.sub(r"array\(\[\[\d+\]\], dtype=int32\)", "array([[<process-local-default>]], dtype=int32)", doc)
        surface[name] = {"available": member is not None, "doc": doc}
    mesh_surface = {}
    for name in ("getMeshFacets", "getConnectivity", "getNodes", "getNbElement", "getNbNodes"):
        member = getattr(aka.Mesh, name, None)
        mesh_surface[name] = {
            "available": member is not None,
            "doc": getattr(member, "__doc__", None) if member is not None else None,
        }
    api = {
        "schema": "pgfrag.extrinsic_onset_refinement.wheel_api_inventory.v1",
        "status": "MODULE_LEVEL_API_ONLY_NO_MODEL",
        "wheel_path": relative(RUNTIME_WHEEL),
        "wheel_sha256": expected_wheel,
        "akantu_version": getattr(aka, "__version__", "UNKNOWN"),
        "python": sys.executable,
        "python_version": platform.python_version(),
        "model_initialized": False,
        "solver_steps": 0,
        "surface": surface,
        "mesh_surface": mesh_surface,
        "documented_shape_access": "getShapes(type, ghost_type=_not_ghost, id=0) -> AkaArray",
        "documented_integration_point_count_access": "getNbIntegrationPoints(type, ghost_type=_not_ghost) -> int",
        "native_quadrature_weights_available_in_public_python_surface": any(surface[key]["available"] for key in ("getWeights", "getIntegrationPointWeights")),
        "native_jacobian_getter_available_in_public_python_surface": any(surface[key]["available"] for key in ("getJacobians", "getJacobian")),
        "interpretation": "FEEngine shape maps and integration point count are public; no exposed per-QP integration-weight/Jacobian getter is claimed. integrate() is not used to fabricate per-point weights.",
    }
    DATA.mkdir(parents=True, exist_ok=True)
    write_json(API_EVIDENCE_PATH, api, exclusive=True)
    print(json.dumps(api, indent=2, sort_keys=True))
    return 0


def array_copy(value) -> np.ndarray:
    return np.asarray(value).copy()


def read_cohesive_field(model, aka, name: str, components: int) -> np.ndarray:
    field_map = model.flattenInternalReal(name, aka._ek_cohesive, aka._not_ghost)
    return array_copy(field_map(aka._cohesive_3d_6, aka._not_ghost)).reshape(-1, components)


def energy_getters(model) -> tuple[dict[str, float | None], dict[str, str]]:
    """Capture native model energies after arrays; kinetic may be a mass barrier.

    The cached v5.0.7 source routes the kinetic getter through lumped-mass
    assembly.  The installed wheel is 5.0.7.post1 with distinct provenance, so
    treat that side effect as possible and record post-getter arrays at events.
    """
    values: dict[str, float | None] = {}
    errors: dict[str, str] = {}
    for label, key in ENERGY_KEYS.items():
        try:
            values[label] = float(model.getEnergy(key))
        except Exception as error:
            values[label] = None
            errors[label] = repr(error)
    return values, errors


def get_material_handles(model):
    count = int(model.getNbMaterials())
    materials = [model.getMaterial(index) for index in range(count)]
    bulk = next((item for item in materials if type(item).__name__ == "MaterialLinearElastic3D"), None)
    cohesive = next((item for item in materials if "Cohesive" in type(item).__name__), None)
    if bulk is None or cohesive is None:
        raise RuntimeError(f"native material handles did not match expected bulk/cohesive objects: {[type(item).__name__ for item in materials]}")
    return materials, bulk, cohesive


def source_lineage(mesh_path: Path, initial_bulk: np.ndarray):
    sys.path.insert(0, str(ROOT / "pgfrag/src"))
    from pgfrag.cohesive_mesh import read_gmsh22
    from pgfrag.facet_ancestry import build_source_facet_table, map_native_cohesive_ancestry

    mesh = read_gmsh22(mesh_path)
    tets = np.asarray([item["nodes"] for item in mesh["tetrahedra"]], dtype=np.int64)
    if not np.array_equal(initial_bulk, tets):
        raise RuntimeError("initial native bulk T4 row/local-corner connectivity differs from frozen MSH; lineage stops")
    node_ids = [int(value) for value in mesh["node_tags"]]
    cell_ids = [int(item["element_id"]) for item in mesh["tetrahedra"]]
    interface = [tuple(map(int, triangle["node_tags"])) for triangle in mesh["triangles"] if triangle["physical_name"] == "interface"]
    table = build_source_facet_table(
        mesh["nodes"], tets, node_ids, cell_ids, sha256_file(mesh_path), eligible_internal_facets=interface
    )
    if len(table.eligible_internal_facets) != 8:
        raise RuntimeError("frozen source lineage did not produce exactly eight eligible interface FIDs")
    return table, map_native_cohesive_ancestry


def native_candidate_facet_map(mesh, aka, source_table, native_nodes: np.ndarray, check_mask: np.ndarray) -> dict[str, object]:
    """Bind the live triangle ElementTypeMap slots to exact source facet IDs.

    The row-to-source-node mapping is accepted only when the public facet mesh
    has the exact same rest-coordinate rows as the original mesh.  This avoids
    nearest-coordinate matching and preserves the source MSH identity anchor.
    """
    facets_mesh = mesh.getMeshFacets()
    facet_nodes = array_copy(facets_mesh.getNodes())
    if facet_nodes.shape != native_nodes.shape or not np.array_equal(facet_nodes, native_nodes):
        raise RuntimeError(
            "public facet-mesh node rows do not exactly match the native/source rest-node rows; "
            "stable facet lineage is unknown"
        )
    connectivity = array_copy(facets_mesh.getConnectivity(aka._triangle_3)).astype(np.int64, copy=False)
    if connectivity.ndim != 2 or connectivity.shape[1] != 3:
        raise RuntimeError(f"public facet connectivity has unexpected shape {connectivity.shape}")
    if check_mask.ndim != 1 or len(check_mask) != len(connectivity):
        raise RuntimeError(
            f"candidate check-mask/facet-connectivity slot mismatch: {check_mask.shape} vs {connectivity.shape}"
        )
    if connectivity.size and (int(connectivity.min()) < 0 or int(connectivity.max()) >= len(native_nodes)):
        raise RuntimeError("public facet connectivity references a node outside the exact source node rows")

    node_id_by_row = tuple(map(int, source_table.source_node_ids))
    slot_records: list[dict[str, object]] = []
    for slot, node_rows in enumerate(connectivity):
        source_ids = tuple(sorted(node_id_by_row[int(row)] for row in node_rows))
        source_facet = source_table.facets_by_node_ids.get(source_ids)
        if source_facet is None:
            raise RuntimeError(f"native triangle slot {slot} has no exact source T4 facet: {source_ids}")
        points = native_nodes[np.asarray(node_rows, dtype=np.int64)]
        cross = np.cross(points[1] - points[0], points[2] - points[0])
        twice_area = float(np.linalg.norm(cross))
        if not math.isfinite(twice_area) or twice_area <= 0.0:
            raise RuntimeError(f"native triangle slot {slot} is degenerate")
        slot_records.append({
            "native_triangle_slot": int(slot),
            "source_node_ids_sorted": list(source_ids),
            "source_facet_identity": source_facet.identity,
            "source_facet_classification": source_facet.classification,
            "area_m2_from_exact_native_facet_coordinates": 0.5 * twice_area,
            "connectivity_native_node_rows_ordered": list(map(int, node_rows)),
            "geometry_unit_normal_from_connectivity_order": (cross / twice_area).tolist(),
        })

    candidate_slots = np.flatnonzero(check_mask).astype(int).tolist()
    candidate_records = [slot_records[slot] for slot in candidate_slots]
    candidate_fids = [str(item["source_facet_identity"]) for item in candidate_records]
    expected_fids = sorted(source_table.eligible_internal_facets)
    if len(candidate_fids) != len(set(candidate_fids)) or sorted(candidate_fids) != expected_fids:
        raise RuntimeError(
            "native check mask does not select the exact eight eligible source FIDs: "
            f"selected={sorted(candidate_fids)}, expected={expected_fids}"
        )
    if any(item["source_facet_classification"] != "internal" for item in candidate_records):
        raise RuntimeError("native check mask selected a non-internal source facet")
    return {
        "status": "EXACT_NODE_ROW_AND_STABLE_FID_MAP",
        "method": "public Mesh.getMeshFacets().getConnectivity(_triangle_3); exact facet-mesh/rest-node rows; source Gmsh node IDs; source T4 facet table FIDs",
        "facet_element_type_slot_count": int(len(slot_records)),
        "candidate_native_slots_ascending": candidate_slots,
        "eligible_fids": expected_fids,
        "candidate_fids_match_frozen_interface_set": True,
        "candidate_records_in_native_slot_order": candidate_records,
        "all_slot_records": slot_records,
        "native_triangle_connectivity": connectivity.copy(),
        "candidate_check_mask": check_mask.copy(),
    }


def validate_symmetric_triangle3_shape_rule(shape_status: dict[str, object], shape_maps: dict[str, np.ndarray], facet_map: dict[str, object], shape_probe: dict[str, object]) -> dict[str, object]:
    """Normalize only the public getShapes layout observed by the zero-step probe."""
    label = "facet_triangle_3"
    key = label + "__shape_functions"
    status = shape_status.get(label, {})
    if not isinstance(status, dict) or status.get("status") != "ok" or key not in shape_maps:
        return {"status": "UNKNOWN", "reason": "native public Tri3 shape map unavailable"}
    raw = np.asarray(shape_maps[key], dtype=np.float64)
    nfacets = int(facet_map["facet_element_type_slot_count"])
    nqp = int(status.get("integration_points_per_element", -1))
    preflight_layout = shape_probe.get("public_api_layout", {}).get("FEEngine.getShapes", {}).get("layout", {})
    expected_layout = preflight_layout.get("layout")
    expected_shape = tuple(map(int, preflight_layout.get("raw_shape", [])))
    expected_candidate_slots = list(map(int, shape_probe.get("public_api_layout", {}).get("Mesh.getMeshFacets_and_getConnectivity", {}).get("candidate_slots_ascending", [])))
    actual_candidate_slots = list(map(int, facet_map.get("candidate_native_slots_ascending", [])))
    if expected_layout != "FLAT_ELEMENT_MAJOR_REPEATED_QP_NODE":
        return {"status": "UNKNOWN", "reason": f"unsupported preflight shape layout {expected_layout!r}; raw retained", "raw_shape": list(raw.shape)}
    if actual_candidate_slots != expected_candidate_slots:
        return {"status": "UNKNOWN", "reason": "candidate facet slots differ from the zero-step shape-layout evidence", "actual_slots": actual_candidate_slots, "preflight_slots": expected_candidate_slots, "raw_shape": list(raw.shape)}
    if nqp != QPS_PER_COHESIVE_ELEMENT or raw.shape != expected_shape or raw.shape != (nfacets * nqp, 3):
        return {
            "status": "UNKNOWN",
            "reason": "runtime shape-array dimensions differ from the zero-step preflight layout/count",
            "raw_shape": list(raw.shape),
            "preflight_raw_shape": list(expected_shape),
            "runtime_nqp": nqp,
        }
    grouped = raw.reshape(nfacets, nqp, 3)
    if not np.allclose(grouped, grouped[0:1], rtol=0.0, atol=1e-14):
        return {"status": "UNKNOWN", "reason": "runtime flattened shapes no longer repeat the verified element-major per-facet rule", "raw_shape": list(raw.shape)}
    shapes = grouped[actual_candidate_slots].copy()
    if status.get("integration_points_per_element") != QPS_PER_COHESIVE_ELEMENT:
        return {"status": "UNKNOWN", "reason": "native facet integration point count is not three", "raw_shape": list(raw.shape)}

    expected = np.asarray([
        [2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0],
        [1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0],
        [1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0],
    ])
    candidate_slots = actual_candidate_slots
    failures = []
    rows_record = []
    for index, slot in enumerate(candidate_slots):
        local = shapes[index]
        row_sums_ok = bool(np.allclose(local.sum(axis=1), 1.0, rtol=0.0, atol=1e-12))
        sorted_rows = np.sort(local, axis=1)
        expected_sorted = np.sort(expected, axis=1)
        values_ok = all(any(np.allclose(row, exp, rtol=0.0, atol=1e-12) for exp in expected_sorted) for row in sorted_rows)
        unique_expected_permutations = all(
            sum(np.allclose(row, expected_row, rtol=0.0, atol=1e-12) for row in local) == 1
            for expected_row in expected
        )
        if not (row_sums_ok and values_ok and unique_expected_permutations):
            failures.append(slot)
        rows_record.append({
            "native_triangle_slot": slot,
            "raw_qp_shape_values_in_native_order": local.tolist(),
            "partition_of_unity_ok": row_sums_ok,
            "three_symmetric_points_each_seen_once": bool(values_ok and unique_expected_permutations),
        })
    if failures:
        return {
            "status": "UNKNOWN",
            "reason": "one or more selected facet shape maps do not match the symmetric three-point Tri3 rule",
            "failed_native_triangle_slots": failures,
            "raw_shape": list(raw.shape),
            "candidate_shape_evidence": rows_record,
        }
    return {
        "status": "SYMMETRIC_THREE_POINT_TRI3_SHAPE_RULE_VERIFIED",
        "raw_layout": expected_layout,
        "raw_shape": list(raw.shape),
        "normalized_all_facet_shape": list(grouped.shape),
        "engine_candidate_shape": list(shapes.shape),
        "quadrature_points_per_element": 3,
        "candidate_shape_evidence": rows_record,
        "candidate_shape_functions_native_qp_order": shapes,
        "qp_order": "raw native within-facet QP order retained; not joined across insertion or resolutions",
        "reference_rule_measure": "equal 1/3 share per point for this symmetric three-point affine Tri3 rule",
    }


def validate_candidate_native_normals(normals_by_slot: np.ndarray, facet_map: dict[str, object]) -> dict[str, object]:
    normals = np.asarray(normals_by_slot, dtype=np.float64)
    candidate_slots = list(map(int, facet_map["candidate_native_slots_ascending"]))
    if normals.ndim != 3 or normals.shape[1:] != (3, 3):
        return {"status": "UNKNOWN", "reason": "native selected normals do not have (candidate facet, 3 QP, xyz) shape", "shape": list(normals.shape)}
    if normals.shape[0] != len(candidate_slots):
        return {"status": "UNKNOWN", "reason": "native selected-normal rows do not match the candidate slot list", "shape": list(normals.shape)}
    evidence = []
    failures = []
    records_by_slot = {int(item["native_triangle_slot"]): item for item in facet_map["candidate_records_in_native_slot_order"]}
    for row, slot in enumerate(candidate_slots):
        local = normals[row]
        geom = np.asarray(records_by_slot[slot]["geometry_unit_normal_from_connectivity_order"], dtype=np.float64)
        lengths = np.linalg.norm(local, axis=1)
        dots = local @ geom
        unit_ok = bool(np.allclose(lengths, 1.0, rtol=0.0, atol=1e-8))
        parallel_ok = bool(np.allclose(np.abs(dots), 1.0, rtol=0.0, atol=1e-8))
        one_orientation = bool(np.all(np.sign(dots) == np.sign(dots[0])))
        ok = bool(unit_ok and parallel_ok and one_orientation)
        if not ok:
            failures.append(slot)
        evidence.append({
            "native_triangle_slot": slot,
            "source_facet_identity": records_by_slot[slot]["source_facet_identity"],
            "native_normals_qp_order": local.tolist(),
            "native_normal_lengths": lengths.tolist(),
            "dot_with_native_connectivity_geometry_normal": dots.tolist(),
            "orientation_sign_against_native_connectivity": int(np.sign(dots[0])) if len(dots) else None,
            "unit_normal_ok": unit_ok,
            "parallel_to_exact_native_triangle_ok": parallel_ok,
            "one_orientation_for_all_three_qps_ok": one_orientation,
        })
    return {
        "status": "NATIVE_NORMALS_MATCH_EXACT_TRIANGLE_GEOMETRY" if not failures else "UNKNOWN",
        "method": "native FacetsFEEngine normals compared pointwise with normalized cross product of exact public facet connectivity coordinates; sign retained per FID",
        "failed_native_triangle_slots": failures,
        "candidate_evidence": evidence,
    }


def coupon_quadrature_measure(facet_map: dict[str, object], shape_rule: dict[str, object], normal_check: dict[str, object]) -> dict[str, object]:
    verified = (
        facet_map.get("status") == "EXACT_NODE_ROW_AND_STABLE_FID_MAP"
        and shape_rule.get("status") == "SYMMETRIC_THREE_POINT_TRI3_SHAPE_RULE_VERIFIED"
        and normal_check.get("status") == "NATIVE_NORMALS_MATCH_EXACT_TRIANGLE_GEOMETRY"
    )
    records = []
    for item in facet_map["candidate_records_in_native_slot_order"]:
        area = float(item["area_m2_from_exact_native_facet_coordinates"])
        records.append({
            "native_triangle_slot": int(item["native_triangle_slot"]),
            "source_facet_identity": item["source_facet_identity"],
            "source_node_ids_sorted": item["source_node_ids_sorted"],
            "area_m2_from_exact_native_facet_coordinates": area,
            "qp_weights_m2_in_native_order": [area / 3.0] * 3 if verified else None,
        })
    return {
        "status": "VERIFIED_COUPON_AFFINE_TRI3" if verified else "UNKNOWN_RAW_ONLY",
        "scope": "only these eight planar affine Tri3 facets in the frozen coupon; not generalizable to curved facets or arbitrary native QP rules",
        "derivation": "exact native facet coordinates/connectivity establish area; public shape map verifies the symmetric 3-point Tri3 rule; native FacetsFEEngine normals verify orientation and three-QP grouping; use area/3 per QP",
        "weights_source": "coupon-specific geometric measure; no public per-QP weight/Jacobian getter is claimed",
        "candidate_records": records,
        "shape_rule_status": shape_rule.get("status"),
        "normal_check_status": normal_check.get("status"),
    }


def native_state(mesh, model, aka, bulk_material, *, capture_candidate_facet: bool, inserter=None, candidate_facet_map=None, shape_rule=None) -> tuple[dict[str, np.ndarray], dict[str, object]]:
    """Copy public getters in legacy order; native energy getters may have side effects."""
    errors: dict[str, str] = {}
    arrays: dict[str, np.ndarray] = {}
    methods = {
        "nodes_rest": lambda: mesh.getNodes(),
        "bulk_tetra_connectivity": lambda: mesh.getConnectivity(aka._tetrahedron_4),
        "cohesive_connectivity": lambda: mesh.getConnectivity(aka._cohesive_3d_6),
        "nodal_displacement": model.getDisplacement,
        "nodal_velocity": model.getVelocity,
        "nodal_acceleration": model.getAcceleration,
        "nodal_mass": model.getMass,
        "nodal_internal_force": model.getInternalForce,
        "nodal_external_force": model.getExternalForce,
        "nodal_blocked_dofs": model.getBlockedDOFs,
    }
    for name, getter in methods.items():
        try:
            value = array_copy(getter())
            if name == "cohesive_connectivity" and value.size == 0:
                value = np.empty((0, 6), dtype=np.int64)
            arrays[name] = value
        except Exception as error:
            errors[name] = repr(error)

    try:
        arrays["bulk_stress_native"] = array_copy(bulk_material.getStress(aka._tetrahedron_4, aka._not_ghost))
    except Exception as error:
        errors["bulk_stress_native"] = repr(error)

    for name, components in COHESIVE_FIELDS:
        key = "cohesive_field__" + name.replace(" ", "_")
        try:
            arrays[key] = read_cohesive_field(model, aka, name, components)
        except Exception as error:
            errors[key] = repr(error)

    energy_values, energy_errors = energy_getters(model)
    errors.update({"energy__" + name: error for name, error in energy_errors.items()})

    counts = {
        "mesh_nodes": int(mesh.getNbNodes()),
        "bulk_tetrahedra": int(mesh.getNbElement(type=aka._tetrahedron_4)),
        "cohesive_elements": int(mesh.getNbElement(type=aka._cohesive_3d_6)),
    }
    meta: dict[str, object] = {
        "getter_errors": errors,
        "energy_J": energy_values,
        **counts,
        "energy_getter_order": "all native array getters copied before energy_getters; cached v5.0.7 kinetic getter may trigger lumped-mass assembly after those arrays",
        "actual_wheel_side_effect_equivalence": "UNKNOWN; installed wheel 5.0.7.post1 differs in provenance from cached v5.0.7 source",
    }

    if capture_candidate_facet:
        if inserter is None:
            meta["candidate_facet_status"] = "unavailable_no_inserter"
        elif candidate_facet_map is None:
            meta["candidate_facet_status"] = "UNKNOWN_no_exact_native_candidate_facet_map"
        else:
            try:
                candidate_mask = array_copy(inserter.getCheckFacets()(aka._triangle_3)).reshape(-1).astype(bool)
                arrays["candidate_triangle_check_mask"] = candidate_mask.copy()
                if len(candidate_mask) != int(candidate_facet_map["facet_element_type_slot_count"]):
                    raise RuntimeError(
                        f"check-mask/native-facet slot count changed: {len(candidate_mask)} vs "
                        f"{candidate_facet_map['facet_element_type_slot_count']}"
                    )
                candidate_slots = np.flatnonzero(candidate_mask).astype(int).tolist()
                all_records = candidate_facet_map["all_slot_records"]
                selected_records = [all_records[slot] for slot in candidate_slots]
                selected_fids = [str(item["source_facet_identity"]) for item in selected_records]
                eligible_fids = set(candidate_facet_map["eligible_fids"])
                if not set(selected_fids).issubset(eligible_fids):
                    raise RuntimeError("live native check mask selected a facet outside the frozen eligible FID set")
                meta["candidate_triangle_native_slots"] = candidate_slots
                meta["candidate_facet_slot_fid_map"] = selected_records
                meta["candidate_facet_map_status"] = "exact native Triangle_3 ElementTypeMap slot to stable source FID; live check mask may shrink after insertion"
                arrays["candidate_triangle_native_connectivity"] = np.asarray(
                    candidate_facet_map["native_triangle_connectivity"], dtype=np.int64
                )[candidate_mask].copy()
                try:
                    facet_engine = model.getFEEngine("FacetsFEEngine")
                    normals = array_copy(facet_engine.getNormalsOnIntegrationPoints(aka._triangle_3, aka._not_ghost))
                    if normals.ndim == 2 and normals.shape[1] == 3 and len(candidate_mask) and normals.shape[0] % len(candidate_mask) == 0:
                        qp_count = normals.shape[0] // len(candidate_mask)
                        normals_by_slot = normals.reshape(len(candidate_mask), qp_count, 3)
                        selected_normals = normals_by_slot[candidate_mask].copy()
                        arrays["candidate_facet_native_normals"] = selected_normals
                        meta["candidate_normals_status"] = "native FacetsFEEngine normals selected by same triangle ElementTypeMap check mask"
                        meta["candidate_facet_qp_count"] = qp_count
                        if qp_count == 3:
                            normal_check = validate_candidate_native_normals(selected_normals, {
                                "candidate_native_slots_ascending": candidate_slots,
                                "candidate_records_in_native_slot_order": selected_records,
                            })
                        else:
                            normal_check = {"status": "UNKNOWN", "reason": f"native candidate normals have {qp_count} QPs per triangle"}
                        meta["candidate_native_normal_orientation_check"] = normal_check
                        dynamic_facet_map = {
                            "status": candidate_facet_map["status"],
                            "candidate_records_in_native_slot_order": selected_records,
                        }
                        meta["coupon_specific_quadrature_measure"] = coupon_quadrature_measure(
                            dynamic_facet_map, shape_rule or {"status": "UNKNOWN"}, normal_check
                        )
                    else:
                        meta["candidate_normals_status"] = f"native shape/mask mismatch: normals={normals.shape}, mask={candidate_mask.shape}"
                except Exception as error:
                    meta["candidate_normals_status"] = "error:" + repr(error)
                try:
                    stress = array_copy(model.getStressOnFacets(aka._triangle_3, aka._not_ghost))
                    if stress.ndim == 2 and len(candidate_mask) and stress.shape[0] % len(candidate_mask) == 0:
                        qp_count = stress.shape[0] // len(candidate_mask)
                        arrays["candidate_facet_native_stress"] = stress.reshape(len(candidate_mask), qp_count, stress.shape[1])[candidate_mask].copy()
                        meta["candidate_facet_stress_status"] = "native getStressOnFacets rows selected by same triangle ElementTypeMap check mask; source stress not reassembled"
                        meta["candidate_stress_components"] = int(stress.shape[1])
                    else:
                        meta["candidate_facet_stress_status"] = f"native shape/mask mismatch: stress={stress.shape}, mask={candidate_mask.shape}"
                except Exception as error:
                    meta["candidate_facet_stress_status"] = "error:" + repr(error)
            except Exception as error:
                meta["candidate_facet_status"] = "error:" + repr(error)
    return arrays, meta


def post_energy_barrier_state(mesh, model, aka) -> tuple[dict[str, np.ndarray], dict[str, str]]:
    """Copy public arrays after energy getters without requesting explicit assembly."""
    getters = {
        "nodes_rest": mesh.getNodes,
        "bulk_tetra_connectivity": lambda: mesh.getConnectivity(aka._tetrahedron_4),
        "cohesive_connectivity": lambda: mesh.getConnectivity(aka._cohesive_3d_6),
        "nodal_displacement": model.getDisplacement,
        "nodal_velocity": model.getVelocity,
        "nodal_acceleration": model.getAcceleration,
        "nodal_mass": model.getMass,
        "nodal_internal_force": model.getInternalForce,
        "nodal_external_force": model.getExternalForce,
        "nodal_blocked_dofs": model.getBlockedDOFs,
    }
    arrays: dict[str, np.ndarray] = {}
    errors: dict[str, str] = {}
    for name, getter in getters.items():
        try:
            value = array_copy(getter())
            if name == "cohesive_connectivity" and value.size == 0:
                value = np.empty((0, 6), dtype=np.int64)
            arrays[name] = value
        except Exception as error:
            errors[name] = repr(error)
    return arrays, errors


def compare_state_arrays(pre_getter: dict[str, np.ndarray], post_getter: dict[str, np.ndarray]) -> dict[str, object]:
    comparison: dict[str, object] = {}
    for name in sorted(set(pre_getter) | set(post_getter)):
        before = pre_getter.get(name)
        after = post_getter.get(name)
        if before is None or after is None:
            comparison[name] = {"status": "UNKNOWN", "before_available": before is not None, "after_available": after is not None}
            continue
        if before.shape != after.shape or before.dtype.kind not in "biufc" or after.dtype.kind not in "biufc":
            comparison[name] = {"status": "UNKNOWN", "before_shape": list(before.shape), "after_shape": list(after.shape), "before_dtype": str(before.dtype), "after_dtype": str(after.dtype)}
            continue
        equal = bool(np.array_equal(before, after, equal_nan=True))
        if before.size:
            difference = np.asarray(after, dtype=np.float64) - np.asarray(before, dtype=np.float64)
            finite = difference[np.isfinite(difference)]
            max_abs = float(np.max(np.abs(finite))) if finite.size else None
            l2 = float(np.linalg.norm(finite)) if finite.size else None
        else:
            max_abs = 0.0
            l2 = 0.0
        comparison[name] = {
            "status": "UNCHANGED" if equal else "CHANGED",
            "exact_array_equal": equal,
            "shape": list(before.shape),
            "dtype_before": str(before.dtype),
            "dtype_after": str(after.dtype),
            "max_abs_difference": max_abs,
            "l2_difference": l2,
        }
    return comparison


def capture_fe_shape_maps(model, aka, *, include_cohesive: bool) -> tuple[dict[str, np.ndarray], dict[str, object]]:
    result: dict[str, np.ndarray] = {}
    status: dict[str, object] = {}
    plans = [
        ("bulk_tetrahedron_4", "", aka._tetrahedron_4),
        ("facet_triangle_3", "FacetsFEEngine", aka._triangle_3),
    ]
    if include_cohesive:
        plans.append(("cohesive_3d_6", "CohesiveFEEngine", aka._cohesive_3d_6))
    else:
        status["cohesive_3d_6"] = {"status": "not_requested_zero_initial_cohesive_elements"}
    for label, engine_name, element_type in plans:
        try:
            engine = model.getFEEngine(engine_name)
            nqp = int(engine.getNbIntegrationPoints(element_type, aka._not_ghost))
            shapes = array_copy(engine.getShapes(element_type, aka._not_ghost, 0))
            result[label + "__shape_functions"] = shapes
            status[label] = {"status": "ok", "engine_name": engine_name or "default FEEngine", "integration_points_per_element": nqp, "shape_functions_shape": list(shapes.shape), "shape_functions_dtype": str(shapes.dtype), "weights": "unavailable_no_public_getter"}
        except Exception as error:
            status[label] = {"status": "unavailable_or_error", "engine_name": engine_name or "default FEEngine", "error": repr(error), "weights": "unavailable_no_public_getter"}
    return result, status


def mass_sums(mass: np.ndarray) -> list[float | None]:
    if mass.ndim != 2:
        return [float(np.sum(mass))]
    return [float(np.sum(mass[:, axis])) for axis in range(mass.shape[1])]


def finite_norm(array: np.ndarray | None) -> float | None:
    if array is None:
        return None
    valid = np.asarray(array)[np.isfinite(array)]
    return float(np.linalg.norm(valid)) if valid.size else None


def ancestry_result(source_table, ancestry_function, state: dict[str, np.ndarray], initial_bulk: np.ndarray) -> dict[str, object]:
    try:
        report = ancestry_function(
            source_table,
            state["nodes_rest"],
            state["bulk_tetra_connectivity"],
            state["cohesive_connectivity"],
            native_initial_bulk_tets=initial_bulk,
        )
        return {"status": "ok", **report}
    except Exception as error:
        return {
            "status": "unknown_lineage_failure",
            "error": repr(error),
            "assumption": "source T4 row/local-corner identity must remain unchanged; exact coordinate consistency is required; no nearest-coordinate repair",
        }


def compressed_snapshot(path: Path, arrays: dict[str, np.ndarray]) -> int:
    with path.open("xb") as stream:
        np.savez_compressed(stream, **arrays)
    return path.stat().st_size


def save_event_snapshot(run_dir: Path, step: int, phase: str, arrays: dict[str, np.ndarray], record: dict[str, object]) -> str:
    path = run_dir / "event_states" / f"step_{step:05d}_{phase}.npz"
    path.parent.mkdir(parents=True, exist_ok=True)
    compressed_snapshot(path, arrays)
    record["native_array_snapshot"] = relative(path)
    return path.name


def state_arrays_for_chunk(state: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    selected: dict[str, np.ndarray] = {}
    if state.get("nodal_displacement", np.empty((0,))).shape[0] > MAX_NODES:
        raise RuntimeError(f"native node count exceeds conservative geometry cap {MAX_NODES}")
    for name in NODAL_ARRAY_NAMES:
        if name not in state:
            raise RuntimeError(f"mandatory nodal getter unavailable: {name}")
        selected[name] = state[name]
    if "bulk_stress_native" not in state:
        raise RuntimeError("mandatory native bulk stress getter unavailable")
    selected["bulk_stress_native"] = state["bulk_stress_native"]
    for name, _ in COHESIVE_FIELDS:
        key = "cohesive_field__" + name.replace(" ", "_")
        if key not in state:
            raise RuntimeError(f"mandatory native cohesive internal field unavailable: {key}")
        selected[key] = state[key]
    if state["cohesive_field__opening"].shape[0] > MAX_COHESIVE_QPS:
        raise RuntimeError(f"native cohesive QP rows exceed conservative cap {MAX_COHESIVE_QPS}")
    return selected


class DChunkWriter:
    def __init__(self, run_dir: Path, max_steps: int):
        self.run_dir = run_dir
        self.max_steps = max_steps
        self.rows: list[tuple[int, dict[str, np.ndarray]]] = []
        self.chunks: list[dict[str, object]] = []
        (run_dir / "d_chunks").mkdir(parents=True, exist_ok=True)

    def add(self, step: int, state: dict[str, np.ndarray]) -> None:
        self.rows.append((step, state_arrays_for_chunk(state)))
        if len(self.rows) >= D_CHUNK_STEPS:
            self.flush()

    def flush(self) -> None:
        if not self.rows:
            return
        start_step = self.rows[0][0]
        end_step = self.rows[-1][0]
        payload: dict[str, np.ndarray] = {}
        node_counts = np.asarray([row[1]["nodal_displacement"].shape[0] for row in self.rows], dtype=np.int16)
        qp_counts = np.asarray([row[1]["cohesive_field__opening"].shape[0] for row in self.rows], dtype=np.int16)
        payload["step"] = np.asarray([row[0] for row in self.rows], dtype=np.int32)
        payload["node_count"] = node_counts
        payload["cohesive_qp_count"] = qp_counts
        for name in CHUNK_ARRAY_NAMES:
            exemplar = self.rows[0][1][name]
            if name in NODAL_ARRAY_NAMES:
                shape = (len(self.rows), MAX_NODES) + exemplar.shape[1:]
                dtype = exemplar.dtype
                data = np.zeros(shape, dtype=dtype)
                for index, (_, row) in enumerate(self.rows):
                    value = row[name]
                    if value.shape[0] > MAX_NODES:
                        raise RuntimeError(f"{name} exceeds node padding cap")
                    data[index, :value.shape[0]] = value
            elif name.startswith("cohesive_field__"):
                components = exemplar.shape[1] if exemplar.ndim == 2 else 1
                data = np.zeros((len(self.rows), MAX_COHESIVE_QPS, components), dtype=exemplar.dtype)
                for index, (_, row) in enumerate(self.rows):
                    value = row[name]
                    if value.shape[0] > MAX_COHESIVE_QPS:
                        raise RuntimeError(f"{name} exceeds cohesive QP padding cap")
                    data[index, :value.shape[0]] = value.reshape(-1, components)
            else:
                data = np.stack([row[name] for _, row in self.rows])
            payload[name] = data
        path = self.run_dir / "d_chunks" / f"steps_{start_step:05d}_{end_step:05d}.npz"
        byte_count = compressed_snapshot(path, payload)
        self.chunks.append({"start_step": start_step, "end_step": end_step, "rows": len(self.rows), "path": relative(path), "bytes": byte_count})
        self.rows.clear()


def write_line(path: Path, value) -> None:
    with path.open("a", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(jsonable(value), sort_keys=True, separators=(",", ":")) + "\n")
        stream.flush()
        os.fsync(stream.fileno())


def geometry_and_top_nodes(rest_nodes: np.ndarray) -> np.ndarray:
    top = np.flatnonzero(np.abs(rest_nodes[:, 2] - 4.0e-3) < 1.0e-15)
    if len(top) != 9:
        raise RuntimeError(f"expected nine source top nodes, found {len(top)}")
    return top


def run_campaign(resolution: str, reviewed_hash: str) -> int:
    if not PROTOCOL_PATH.is_file():
        raise RuntimeError("frozen protocol missing; run --prepare first")
    actual_protocol_hash = sha256_file(PROTOCOL_PATH)
    if actual_protocol_hash.lower() != reviewed_hash.lower():
        raise RuntimeError("reviewed protocol SHA does not match frozen_protocol.json; solver not started")
    protocol = json.loads(PROTOCOL_PATH.read_text(encoding="utf-8"))
    if protocol.get("status") != "FROZEN_FOR_ROOT_REVIEW_NO_CAMPAIGN_LAUNCHED":
        raise RuntimeError("protocol status is not the frozen pre-review state")
    if protocol["driver"]["sha256"] != sha256_file(Path(__file__)):
        raise RuntimeError("driver has changed since protocol freeze; solver not started")
    if protocol["output_forecast"]["total_forecast_bytes"] > OUTPUT_CAP_BYTES:
        raise RuntimeError("frozen output forecast exceeds cap")
    entries = {item["resolution"]: item for item in protocol["time_matrix"]}
    if resolution not in entries:
        raise ValueError(f"unknown resolution {resolution!r}")
    entry = entries[resolution]
    if entry["status"] != "PLANNED_ROOT_REVIEW":
        raise RuntimeError("this resolution is reserved and requires separate root confirmation")
    run_dir = DATA / resolution
    if run_dir.exists():
        raise RuntimeError(f"single-attempt directory already exists; no overwrite/retry: {relative(run_dir)}")
    prior_elapsed = 0.0
    for other in DATA.iterdir():
        status_path = other / "run_status.json" if other.is_dir() else None
        if status_path and status_path.is_file():
            try:
                prior_elapsed += float(json.loads(status_path.read_text(encoding="utf-8")).get("elapsed_s", 0.0))
            except Exception:
                prior_elapsed += WALL_CAP_PER_RUN_S
    if prior_elapsed + WALL_CAP_PER_RUN_S > WALL_CAP_CAMPAIGN_S:
        raise RuntimeError("campaign cumulative wall budget would exceed 20 minutes")
    existing_output_bytes = existing_bytes(DATA)
    run_output_forecast = int(protocol["output_forecast"]["per_resolution"][resolution]["total_forecast_bytes"])
    if existing_output_bytes + run_output_forecast > OUTPUT_CAP_BYTES:
        raise RuntimeError("remaining forecast allowance would exceed total output cap")

    run_dir.mkdir(parents=False, exist_ok=False)
    started_utc = datetime.now(timezone.utc).isoformat()
    started = time.monotonic()
    dt = float(Decimal(entry["dt_decimal_s"]))
    steps = int(entry["steps_to_common_endpoint"])
    status: dict[str, object] = {
        "schema": "pgfrag.extrinsic_onset_refinement.run_status.v1",
        "run_id": resolution,
        "status": "SINGLE_NATIVE_ATTEMPT_RUNNING",
        "started_utc": started_utc,
        "argv": sys.argv,
        "cwd": os.getcwd(),
        "python": sys.executable,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "protocol_sha256": actual_protocol_hash,
        "driver_sha256": sha256_file(Path(__file__)),
        "resolution": resolution,
        "dt_decimal_s": entry["dt_decimal_s"],
        "dt_binary64_s": dt,
        "step_cap": steps,
        "solver_attempts": 1,
        "solver_retries": 0,
        "data_bytes_before_run": existing_output_bytes,
        "run_output_forecast_bytes": run_output_forecast,
        "total_output_cap_bytes": OUTPUT_CAP_BYTES,
        "thread_environment": {name: os.environ.get(name) for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
    }
    write_json(run_dir / "run_status.json", status, exclusive=True)
    chunks: DChunkWriter | None = None
    trace_stream = None
    try:
        import akantu as aka

        if getattr(aka, "__version__", None) != "5.0.7-post1":
            raise RuntimeError(f"unexpected runtime version: {getattr(aka, '__version__', None)}")
        if sha256_file(RUNTIME_WHEEL) != LOCKED_HASHES["runtime_wheel"]:
            raise RuntimeError("runtime wheel changed after protocol freeze")
        manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
        if manifest["material"]["beta"] != 3.0 or manifest["material"]["sigma_c_Pa"] != 70e6:
            raise RuntimeError("frozen material law parameters differ from accepted beta=3/sigma=70 MPa case")

        aka.parseInput(str(INPUT / "material.dat"))
        mesh = aka.Mesh(3)
        mesh.read(str(INPUT / "block_2x2x4mm.msh"))
        model = aka.SolidMechanicsModelCohesive(mesh)
        model.initFull(_analysis_method=aka._explicit_lumped_mass, _is_extrinsic=True)
        model.assembleMassLumped()
        inserter = model.getElementInserter()
        inserter.addPhysicalSurface("interface")
        model.updateAutomaticInsertion()
        check_mask0 = array_copy(inserter.getCheckFacets()(aka._triangle_3)).reshape(-1).astype(bool)
        inserted_mask0 = array_copy(inserter.getInsertionFacets()(aka._triangle_3)).reshape(-1).astype(bool)
        if int(check_mask0.sum()) != 8 or int(inserted_mask0.sum()) != 0:
            raise RuntimeError(f"initial eligibility mismatch: eligible={int(check_mask0.sum())}, inserted={int(inserted_mask0.sum())}")

        materials, bulk_material, cohesive_material = get_material_handles(model)
        initial_rest = array_copy(mesh.getNodes())
        initial_bulk = array_copy(mesh.getConnectivity(aka._tetrahedron_4))
        initial_coh = array_copy(mesh.getConnectivity(aka._cohesive_3d_6))
        if initial_coh.size != 0:
            raise RuntimeError("fresh extrinsic setup unexpectedly contains cohesive elements")
        source_table, map_ancestry = source_lineage(INPUT / "block_2x2x4mm.msh", initial_bulk)
        candidate_facet_map = native_candidate_facet_map(mesh, aka, source_table, initial_rest, check_mask0)
        top_nodes = geometry_and_top_nodes(initial_rest)
        model.applyBC(aka.FixedValue(0.0, aka._z), "bottom")

        class RampTop(aka.DirichletFunctor):
            def __init__(self, amplitude: float, duration: float):
                super().__init__(aka._z)
                self.axis = int(aka._z)
                self.amplitude = float(amplitude)
                self.duration = float(duration)
                self.time = 0.0

            def set_time(self, value: float) -> None:
                self.time = float(value)

            def __call__(self, node, flags, displacement, coordinates):
                flags[self.axis] = True
                displacement[self.axis] = self.amplitude * min(max(self.time / self.duration, 0.0), 1.0)

        ramp = RampTop(manifest["loading"]["top_displacement_m"], manifest["loading"]["ramp_duration_s"])
        model.setTimeStep(dt)
        shape_maps, shape_status = capture_fe_shape_maps(model, aka, include_cohesive=False)
        shape_probe = json.loads(SHAPE_PREFLIGHT_EVIDENCE.read_text(encoding="utf-8"))
        triangle_shape_rule = validate_symmetric_triangle3_shape_rule(shape_status, shape_maps, candidate_facet_map, shape_probe)
        if shape_maps:
            with (run_dir / "feengine_shape_maps.npz").open("xb") as stream:
                np.savez_compressed(stream, **shape_maps)
        write_json(run_dir / "feengine_shape_map_status.json", {
            "public_shape_maps": shape_status,
            "public_weight_or_jacobian_getter": "unavailable_in_wheel_python_binding; no fabricated weights",
            "coupon_specific_tri3_shape_rule": triangle_shape_rule,
            "nodal_force_reconstruction": "RAW_ONLY_UNTIL_COUPON_AFFINE_TRI3_MEASURE_AND_STABLE_MAPPING_ARE_VERIFIED",
        }, exclusive=True)
        candidate_records_by_slot = candidate_facet_map["all_slot_records"]
        selected_candidate_records = candidate_facet_map["candidate_records_in_native_slot_order"]
        initial_normal_check: dict[str, object]
        initial_normal_array: np.ndarray | None = None
        try:
            facet_engine = model.getFEEngine("FacetsFEEngine")
            all_initial_normals = array_copy(facet_engine.getNormalsOnIntegrationPoints(aka._triangle_3, aka._not_ghost))
            if (
                all_initial_normals.ndim == 2
                and all_initial_normals.shape[1] == 3
                and len(check_mask0)
                and all_initial_normals.shape[0] % len(check_mask0) == 0
            ):
                initial_qps = all_initial_normals.shape[0] // len(check_mask0)
                initial_normals_by_slot = all_initial_normals.reshape(len(check_mask0), initial_qps, 3)
                initial_normal_array = initial_normals_by_slot[check_mask0].copy()
                initial_normal_check = validate_candidate_native_normals(initial_normal_array, {
                    "candidate_native_slots_ascending": candidate_facet_map["candidate_native_slots_ascending"],
                    "candidate_records_in_native_slot_order": selected_candidate_records,
                }) if initial_qps == 3 else {
                    "status": "UNKNOWN",
                    "reason": f"native initial facet normals have {initial_qps} QPs per triangle",
                    "shape": list(initial_normal_array.shape),
                }
            else:
                initial_normal_check = {
                    "status": "UNKNOWN",
                    "reason": "native initial facet-normal array does not align with the exact triangle slot map",
                    "shape": list(all_initial_normals.shape),
                    "facet_slot_count": int(len(check_mask0)),
                }
        except Exception as error:
            initial_normal_check = {"status": "UNKNOWN", "reason": "public initial FacetsFEEngine normals getter failed", "error": repr(error)}
        initial_measure = coupon_quadrature_measure(candidate_facet_map, triangle_shape_rule, initial_normal_check)
        if initial_normal_array is not None:
            with (run_dir / "initial_facet_quadrature_evidence.npz").open("xb") as stream:
                np.savez_compressed(
                    stream,
                    candidate_triangle_native_connectivity=np.asarray(candidate_facet_map["native_triangle_connectivity"])[check_mask0].copy(),
                    candidate_triangle_check_mask=check_mask0.copy(),
                    candidate_facet_native_normals=initial_normal_array,
                    candidate_facet_native_shape_functions=np.asarray(
                        triangle_shape_rule.get("candidate_shape_functions_native_qp_order", np.empty((0, 3, 3)))
                    ).copy(),
                )
        write_json(run_dir / "candidate_facet_lineage_initial.json", {
            "status": candidate_facet_map["status"],
            "method": candidate_facet_map["method"],
            "facet_element_type_slot_count": candidate_facet_map["facet_element_type_slot_count"],
            "eligible_fids": candidate_facet_map["eligible_fids"],
            "candidate_native_slots_ascending": candidate_facet_map["candidate_native_slots_ascending"],
            "candidate_fids_match_frozen_interface_set": candidate_facet_map["candidate_fids_match_frozen_interface_set"],
            "candidate_records_in_native_slot_order": selected_candidate_records,
            "all_slot_source_fids": [item["source_facet_identity"] for item in candidate_records_by_slot],
            "initial_native_normal_check": initial_normal_check,
            "source_facet_identity_anchor": source_table.source_anchor_hash,
            "qp_order_policy": "native FEEngine slot and within-triangle QP order are retained; no order join is made across insertion/resolutions",
        }, exclusive=True)
        write_json(run_dir / "facet_quadrature_measure_initial.json", {
            **initial_measure,
            "initial_native_normal_check": initial_normal_check,
            "initial_triangle_shape_rule": triangle_shape_rule,
            "event_B_measure_is_recaptured": True,
            "solver_steps_at_capture": 0,
            "observation_assembly_called": False,
        }, exclusive=True)
        initial_state, initial_meta = native_state(mesh, model, aka, bulk_material, capture_candidate_facet=False)
        initial_ancestry = ancestry_result(source_table, map_ancestry, initial_state, initial_bulk)
        initial_state_path = run_dir / "initial_state.npz"
        compressed_snapshot(initial_state_path, initial_state)
        write_json(run_dir / "initial_state_record.json", {
            "phase": "INITIAL_BEFORE_STEP_1",
            "driver_time_s": "0",
            "arrays": relative(initial_state_path),
            "native_state": initial_meta,
            "ancestry": initial_ancestry,
            "mass_status": "initial mass after one required assembleMassLumped; no time step has run",
        }, exclusive=True)

        steps_name = f"steps_{steps}"
        trace_path = run_dir / "scalar_trace.csv"
        trace_fields = [
            "step", "dt_s", "driver_time_s", "bc_setpoint_m", "ramp_fraction",
            "cohesive_check_return", "insertions_this_step", "insertions_cumulative",
            "mesh_nodes", "bulk_tetrahedra", "cohesive_elements",
            "mass_x_kg", "mass_y_kg", "mass_z_kg", "mass_min_kg",
            "internal_force_l2_N", "external_force_l2_N", "top_internal_reaction_estimate_N",
            "blocked_dof_count", "cohesive_opening_l2_m", "cohesive_traction_l2_Pa",
            "bulk_stress_l2_Pa", "bulk_stress_max_abs_Pa",
            "kinetic_J", "bulk_potential_J", "cohesive_dissipated_J", "cohesive_reversible_J",
            "cohesive_contact_native_key_J", "contact_alias_diagnostic_J", "energy_status",
            "candidate_fid_set_json", "candidate_measure_status",
            "ancestry_status", "connectivity_sha256", "rest_nodes_sha256", "dense_d_chunk", "dense_d_row",
        ]
        trace_stream = trace_path.open("x", newline="", encoding="utf-8")
        trace = csv.DictWriter(trace_stream, fieldnames=trace_fields)
        trace.writeheader()
        chunks = DChunkWriter(run_dir, steps)
        events: list[dict[str, object]] = []
        previous_b: deque[tuple[int, dict[str, np.ndarray], dict[str, object], dict[str, object]]] = deque(maxlen=4)
        previous_active: list[str] = []
        cumulative_insertions = 0
        previous_d_state = initial_state
        previous_d_meta = initial_meta
        previous_d_ancestry = initial_ancestry
        current_step = 0

        for step in range(1, steps + 1):
            current_step = step
            scheduled = Decimal(step) * Decimal(entry["dt_decimal_s"])
            prescribed_time = float(scheduled)
            # A is copied before the step's boundary update so it can be persisted
            # if the following check inserts one or more eligible facets.
            a_state, a_meta = native_state(mesh, model, aka, bulk_material, capture_candidate_facet=False)
            a_ancestry = ancestry_result(source_table, map_ancestry, a_state, initial_bulk)
            ramp.set_time(prescribed_time)
            model.applyBC(ramp, "top")
            b_state, b_meta = native_state(
                mesh, model, aka, bulk_material, capture_candidate_facet=True,
                inserter=inserter, candidate_facet_map=candidate_facet_map, shape_rule=triangle_shape_rule,
            )
            b_ancestry = ancestry_result(source_table, map_ancestry, b_state, initial_bulk)

            api_return = int(model.checkCohesiveStress())
            cohesive_count = int(mesh.getNbElement(type=aka._cohesive_3d_6))
            insertion_delta = cohesive_count - int(b_meta["cohesive_elements"])
            if insertion_delta < 0 or api_return != insertion_delta:
                raise RuntimeError(f"native check return/count delta disagreement at step {step}: return={api_return}, delta={insertion_delta}")

            if insertion_delta:
                # Flush only the four preceding no-insertion B states, then save
                # the event-step A/B/C states with snapshot-local array ordering.
                for prior_step, prior_state, prior_meta, prior_ancestry in previous_b:
                    prior_record = {
                        "step": prior_step,
                        "phase": "B_PREINSERT_ROLLING",
                        "driver_time_s": str(Decimal(prior_step) * Decimal(entry["dt_decimal_s"])),
                        "native_state": prior_meta,
                        "ancestry": prior_ancestry,
                        "note": "copied after that step applyBC and before checkCohesiveStress",
                    }
                    prior_record["native_array_snapshot"] = save_event_snapshot(run_dir, prior_step, "B_rolling", prior_state, prior_record)
                    write_line(run_dir / "phase_events.jsonl", prior_record)
                previous_b.clear()

                event_a = {
                    "step": step,
                    "phase": "A",
                    "driver_time_before_bc_s": str(Decimal(step - 1) * Decimal(entry["dt_decimal_s"])),
                    "scheduled_boundary_time_s": str(scheduled),
                    "native_state": a_meta,
                    "ancestry": a_ancestry,
                    "note": "before applying the current step's top BC; exact step-local native getter snapshot",
                }
                event_a["native_array_snapshot"] = save_event_snapshot(run_dir, step, "A", a_state, event_a)
                event_b = {
                    "step": step,
                    "phase": "B",
                    "driver_time_s": str(scheduled),
                    "applied_top_displacement_m": manifest["loading"]["top_displacement_m"] * min(max(prescribed_time / manifest["loading"]["ramp_duration_s"], 0.0), 1.0),
                    "native_state": b_meta,
                    "ancestry": b_ancestry,
                    "native_array_snapshot": None,
                }
                event_b["native_array_snapshot"] = save_event_snapshot(run_dir, step, "B", b_state, event_b)
                # Preserve legacy C ordering: copy native arrays first, then
                # call energy getters. The cached kinetic getter can act as an
                # internal mass-assembly barrier after these pre-energy arrays.
                c_state, c_meta = native_state(mesh, model, aka, bulk_material, capture_candidate_facet=False)
                c_post_energy_state, c_post_energy_errors = post_energy_barrier_state(mesh, model, aka)
                c_post_energy_comparison = compare_state_arrays(
                    {name: c_state[name] for name in BARRIER_ARRAY_NAMES if name in c_state},
                    {name: c_post_energy_state[name] for name in BARRIER_ARRAY_NAMES if name in c_post_energy_state},
                )
                cohesive_shape_maps, cohesive_shape_status = capture_fe_shape_maps(model, aka, include_cohesive=True)
                for shape_name, shape_value in cohesive_shape_maps.items():
                    if shape_name.startswith("cohesive_3d_6__"):
                        c_state["feengine__" + shape_name] = shape_value
                event_c_shape_status = cohesive_shape_status.get("cohesive_3d_6")
                c_ancestry = ancestry_result(source_table, map_ancestry, c_state, initial_bulk)
                event_c = {
                    "step": step,
                    "phase": "C",
                    "driver_time_s": str(scheduled),
                    "native_state": c_meta,
                    "ancestry": c_ancestry,
                    "native_array_snapshot": None,
                    "cohesive_feengine_shape_map": event_c_shape_status,
                    "pre_energy_state_semantics": "C arrays were copied after insertion and before the legacy native energy getters; appended node mass/force values may be stale",
                    "energy_values_phase": "after C pre-energy arrays; getEnergy('kinetic') may internally reassemble lumped mass under cached v5.0.7 source; installed wheel side-effect equivalence is UNKNOWN",
                    "post_energy_getter_state_errors": c_post_energy_errors,
                    "post_energy_getter_state_comparison": c_post_energy_comparison,
                    "post_energy_state_semantics": "public nodal/topology getters copied after the legacy energy sequence; no additional energy/force/mass assembly was requested",
                    "mass_interpretation": "C_pre_energy mass is potentially stale; separately compare C_post_energy_getter mass after possible native kinetic-getter barrier",
                    "material_getter_crosschecks": {},
                }
                try:
                    event_c["material_getter_crosschecks"]["cohesive_getOpening"] = array_copy(cohesive_material.getOpening(aka._cohesive_3d_6, aka._not_ghost))
                    event_c["material_getter_crosschecks"]["cohesive_getTraction"] = array_copy(cohesive_material.getTraction(aka._cohesive_3d_6, aka._not_ghost))
                except Exception as error:
                    event_c["material_getter_crosschecks_error"] = repr(error)
                event_c["native_array_snapshot"] = save_event_snapshot(run_dir, step, "C_pre_energy", c_state, {})
                c_post_record: dict[str, object] = {}
                save_event_snapshot(run_dir, step, "C_post_energy_getter", c_post_energy_state, c_post_record)
                event_c["post_energy_getter_array_snapshot"] = c_post_record.get("native_array_snapshot")
                write_line(run_dir / "phase_events.jsonl", event_a)
                write_line(run_dir / "phase_events.jsonl", event_b)
                write_line(run_dir / "phase_events.jsonl", event_c)
                events.append({
                    "step": step,
                    "event_driver_time_s": str(scheduled),
                    "check_return": api_return,
                    "insertion_delta": insertion_delta,
                    "candidate_fids_before_check": [
                        item["source_facet_identity"] for item in b_meta.get("candidate_facet_slot_fid_map", [])
                    ],
                    "candidate_measure_status_before_check": b_meta.get("coupon_specific_quadrature_measure", {}).get("status", "UNKNOWN"),
                    "before_active_fids": b_ancestry.get("active_facet_identities", []),
                    "after_active_fids": c_ancestry.get("active_facet_identities", []),
                    "ancestry_status_C": c_ancestry.get("status"),
                    "C_pre_energy_mass_component_sums_kg": mass_sums(c_state["nodal_mass"]) if "nodal_mass" in c_state else None,
                    "C_post_energy_getter_mass_component_sums_kg": mass_sums(c_post_energy_state["nodal_mass"]) if "nodal_mass" in c_post_energy_state else None,
                    "C_post_energy_getter_mass_comparison": c_post_energy_comparison.get("nodal_mass"),
                })
                cumulative_insertions += insertion_delta
            else:
                previous_b.append((step, b_state, b_meta, b_ancestry))

            # The legacy driver calls solveStep immediately after insertion
            # check; no explicit assembly is called here to inspect C or D.
            model.solveStep()
            d_state, d_meta = native_state(mesh, model, aka, bulk_material, capture_candidate_facet=False)
            d_post_energy_state: dict[str, np.ndarray] = {}
            d_post_energy_errors: dict[str, str] = {}
            d_post_energy_comparison: dict[str, object] = {}
            if insertion_delta:
                d_post_energy_state, d_post_energy_errors = post_energy_barrier_state(mesh, model, aka)
                d_post_energy_comparison = compare_state_arrays(
                    {name: d_state[name] for name in BARRIER_ARRAY_NAMES if name in d_state},
                    {name: d_post_energy_state[name] for name in BARRIER_ARRAY_NAMES if name in d_post_energy_state},
                )
            d_ancestry = ancestry_result(source_table, map_ancestry, d_state, initial_bulk)
            if d_ancestry.get("status") != "ok":
                ancestry_label = "unknown_lineage_failure"
            else:
                active_now = d_ancestry["active_facet_identities"]
                ancestry_label = "ok"
                if active_now != previous_active:
                    write_line(run_dir / "ancestry_changes.jsonl", {
                        "step": step,
                        "driver_time_s": str(scheduled),
                        "active_facet_identities": active_now,
                        "missing_eligible_facet_identities": d_ancestry["missing_eligible_facet_identities"],
                        "native_node_to_source_node_id": d_ancestry["native_node_to_source_node_id"],
                        "lineage_assumptions": d_ancestry["lineage_assumptions"],
                        "qp_ordinal_stability_claimed": False,
                    })
                    previous_active = list(active_now)

            row_number = len(chunks.rows)
            chunks.add(step, d_state)
            chunk_name = f"steps_{(step - 1) // D_CHUNK_STEPS * D_CHUNK_STEPS + 1:05d}_{min(((step - 1) // D_CHUNK_STEPS + 1) * D_CHUNK_STEPS, steps):05d}.npz"
            mass = d_state.get("nodal_mass")
            mass_sum = mass_sums(mass) if mass is not None else [None, None, None]
            top_internal = d_state.get("nodal_internal_force")
            top_reaction_estimate = -float(top_internal[top_nodes, 2].sum()) if top_internal is not None and len(top_internal) >= len(top_nodes) else None
            opening = d_state.get("cohesive_field__opening")
            traction = d_state.get("cohesive_field__tractions")
            bulk_stress = d_state.get("bulk_stress_native")
            blocked = d_state.get("nodal_blocked_dofs")
            energy = d_meta["energy_J"]
            energy_error_labels = sorted(name for name in d_meta["getter_errors"] if name.startswith("energy__"))
            trace.writerow({
                "step": step,
                "dt_s": repr(dt),
                "driver_time_s": str(scheduled),
                "bc_setpoint_m": repr(manifest["loading"]["top_displacement_m"] * min(max(prescribed_time / manifest["loading"]["ramp_duration_s"], 0.0), 1.0)),
                "ramp_fraction": repr(min(max(prescribed_time / manifest["loading"]["ramp_duration_s"], 0.0), 1.0)),
                "cohesive_check_return": api_return,
                "insertions_this_step": insertion_delta,
                "insertions_cumulative": cumulative_insertions,
                "mesh_nodes": d_meta["mesh_nodes"],
                "bulk_tetrahedra": d_meta["bulk_tetrahedra"],
                "cohesive_elements": d_meta["cohesive_elements"],
                "mass_x_kg": mass_sum[0] if len(mass_sum) > 0 else None,
                "mass_y_kg": mass_sum[1] if len(mass_sum) > 1 else None,
                "mass_z_kg": mass_sum[2] if len(mass_sum) > 2 else None,
                "mass_min_kg": float(np.min(mass)) if mass is not None and mass.size else None,
                "internal_force_l2_N": finite_norm(top_internal),
                "external_force_l2_N": finite_norm(d_state.get("nodal_external_force")),
                "top_internal_reaction_estimate_N": top_reaction_estimate,
                "blocked_dof_count": int(np.count_nonzero(blocked)) if blocked is not None else None,
                "cohesive_opening_l2_m": finite_norm(opening),
                "cohesive_traction_l2_Pa": finite_norm(traction),
                "bulk_stress_l2_Pa": finite_norm(bulk_stress),
                "bulk_stress_max_abs_Pa": float(np.max(np.abs(bulk_stress))) if bulk_stress is not None and bulk_stress.size else None,
                **{name: energy.get(name) for name in ENERGY_KEYS},
                "energy_status": "OK" if not energy_error_labels else "UNKNOWN:" + ";".join(energy_error_labels),
                "candidate_fid_set_json": json.dumps(
                    [item["source_facet_identity"] for item in b_meta.get("candidate_facet_slot_fid_map", [])],
                    separators=(",", ":"),
                ),
                "candidate_measure_status": b_meta.get("coupon_specific_quadrature_measure", {}).get("status", "UNKNOWN"),
                "ancestry_status": ancestry_label,
                "connectivity_sha256": hashlib.sha256(d_state["bulk_tetra_connectivity"].tobytes() + d_state["cohesive_connectivity"].tobytes()).hexdigest(),
                "rest_nodes_sha256": hashlib.sha256(d_state["nodes_rest"].tobytes()).hexdigest(),
                "dense_d_chunk": chunk_name,
                "dense_d_row": row_number,
            })
            trace_stream.flush()
            if insertion_delta:
                event_d = {
                    "step": step,
                    "phase": "D",
                    "driver_time_s": str(scheduled),
                    "native_state": d_meta,
                    "ancestry_status": ancestry_label,
                    "ancestry": d_ancestry,
                    "pre_energy_state_semantics": "D arrays were copied immediately after solveStep and before the legacy native energy getters; these arrays are also represented in the dense D chunk",
                    "energy_values_phase": "after D pre-energy arrays; no second energy getter was called for side-effect comparison",
                    "post_energy_getter_state_errors": d_post_energy_errors,
                    "post_energy_getter_state_comparison": d_post_energy_comparison,
                    "post_energy_state_semantics": "public nodal/topology getters copied after the legacy energy sequence; no additional energy/force/mass assembly was requested",
                    "mass_interpretation": "D_pre_energy arrays are post-solve; compare against D_post_energy_getter arrays because the legacy kinetic getter may be an internal barrier if native mass is marked dirty",
                    "dense_d_chunk": chunk_name,
                    "dense_d_row": row_number,
                }
                d_post_record: dict[str, object] = {}
                save_event_snapshot(run_dir, step, "D_post_energy_getter", d_post_energy_state, d_post_record)
                event_d["post_energy_getter_array_snapshot"] = d_post_record.get("native_array_snapshot")
                write_line(run_dir / "phase_events.jsonl", event_d)
                events[-1]["after_solve_active_fids"] = d_ancestry.get("active_facet_identities", [])
                events[-1]["ancestry_status_D"] = ancestry_label
            if step == 1 or step % 25 == 0:
                write_line(run_dir / "phase_events.jsonl", {
                    "step": step,
                    "phase": "D",
                    "driver_time_s": str(scheduled),
                    "native_state": d_meta,
                    "ancestry_status": ancestry_label,
                    "native_node_to_source_node_id": d_ancestry.get("native_node_to_source_node_id"),
                    "active_facet_identities": d_ancestry.get("active_facet_identities"),
                    "mass_interpretation": "post-solve getter result; operative candidate after native solveStep, not a pass by itself",
                    "dense_d_chunk": chunk_name,
                    "dense_d_row": row_number,
                })
            previous_d_state, previous_d_meta, previous_d_ancestry = d_state, d_meta, d_ancestry
            if step == steps or time.monotonic() - started > WALL_CAP_PER_RUN_S:
                chunks.flush()
            if existing_bytes(DATA) > OUTPUT_CAP_BYTES:
                raise RuntimeError("actual data folder exceeded the frozen 256 MiB cap")
            if time.monotonic() - started > WALL_CAP_PER_RUN_S:
                raise TimeoutError(f"{WALL_CAP_PER_RUN_S}s single-run wall cap reached at step {step}; no retry")

        # Preserve the terminal rolling B history even if no later insertion
        # occurred to trigger its flush. This is pre-check evidence only.
        for prior_step, prior_state, prior_meta, prior_ancestry in previous_b:
            terminal_record = {
                "step": prior_step,
                "phase": "B_TERMINAL_ROLLING",
                "driver_time_s": str(Decimal(prior_step) * Decimal(entry["dt_decimal_s"])),
                "native_state": prior_meta,
                "ancestry": prior_ancestry,
                "note": "last four no-insertion B states retained at common endpoint",
            }
            terminal_record["native_array_snapshot"] = save_event_snapshot(run_dir, prior_step, "B_terminal", prior_state, terminal_record)
            write_line(run_dir / "phase_events.jsonl", terminal_record)
        chunks.flush()
        trace_stream.close()
        trace_stream = None
        write_json(run_dir / "events.json", {
            "events": events,
            "event_clock": "scheduled BC time step*dt at which checkCohesiveStress was called, not native model time",
            "event_count": len(events),
        }, exclusive=True)
        status.update({
            "status": "COMPLETED_COMMON_ENDPOINT",
            "ended_utc": datetime.now(timezone.utc).isoformat(),
            "elapsed_s": time.monotonic() - started,
            "steps_completed": steps,
            "endpoint_driver_time_s": str(Decimal(steps) * Decimal(entry["dt_decimal_s"])),
            "events": events,
            "d_chunks": chunks.chunks,
            "output_bytes_data_folder": existing_bytes(run_dir),
            "solver_retries": 0,
            "boundary_work_ledger": "UNKNOWN_UNVERIFIED_ACTUATOR_FORCE_AND_KINEMATIC_PHASE; top reaction is only an estimate",
        })
        write_json(run_dir / "run_status.json", status)
        return 0
    except BaseException as error:
        if chunks is not None:
            try:
                chunks.flush()
            except Exception as flush_error:
                status["flush_error"] = repr(flush_error)
        if trace_stream is not None:
            try:
                trace_stream.flush()
                trace_stream.close()
            except Exception:
                pass
        status.update({
            "status": "FAILED_PRESERVED_NO_RETRY",
            "ended_utc": datetime.now(timezone.utc).isoformat(),
            "elapsed_s": time.monotonic() - started,
            "last_step_started": current_step,
            "exception": repr(error),
            "traceback": traceback.format_exc(),
            "solver_retries": 0,
        })
        write_json(run_dir / "run_status.json", status)
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--inspect-wheel-api", action="store_true", help="module-level API inventory only; no model or solver")
    action.add_argument("--prepare", action="store_true", help="hash/geometry/budget validation only; no Akantu model")
    action.add_argument("--launch", action="store_true", help="run one root-reviewed native resolution")
    parser.add_argument("--resolution", choices=("h", "h_over_2", "h_over_4", "h_over_8_reserved"))
    parser.add_argument("--reviewed-protocol-sha256")
    args = parser.parse_args()
    if args.inspect_wheel_api:
        if args.resolution or args.reviewed_protocol_sha256:
            parser.error("API inspection takes no run parameters")
        return inspect_wheel_api()
    if args.prepare:
        if args.resolution or args.reviewed_protocol_sha256:
            parser.error("--prepare takes no run parameters")
        return prepare()
    if not args.resolution or not args.reviewed_protocol_sha256:
        parser.error("--launch requires --resolution and --reviewed-protocol-sha256")
    if not ROOT_REVIEW_REQUIRED:
        parser.error("launch gate disabled unexpectedly")
    return run_campaign(args.resolution, args.reviewed_protocol_sha256)


if __name__ == "__main__":
    raise SystemExit(main())
