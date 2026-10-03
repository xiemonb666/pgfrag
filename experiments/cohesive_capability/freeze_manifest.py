"""Freeze all inputs and the unvalidated tiny-case settings before simulation."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REVIEW = ROOT.parents[2] / "relicweave-synlab" / "tmp" / "review_20261004" / "implementation" / "cohesive_capability"
REVIEW.mkdir(parents=True, exist_ok=True)
MANIFEST = REVIEW / "manifest.json"
if MANIFEST.exists():
    raise SystemExit(f"refusing to overwrite frozen manifest: {MANIFEST}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


preflight = json.loads((ROOT / "preflight.json").read_text(encoding="utf-8"))
geometry = json.loads((ROOT / "input" / "generated_geometry.json").read_text(encoding="utf-8"))
stable_dt = float(preflight["stable_time_step_s"])
endtime = 1.0e-4
dt_limit_factor = 0.8
base_steps = math.ceil(endtime / (dt_limit_factor * stable_dt))
base_dt = endtime / base_steps
fine_dt = base_dt / 2.0
fine_steps = base_steps * 2
if base_steps > 10_000 or fine_steps > 10_000:
    raise SystemExit("frozen cases exceed the 10,000-step cap")

input_files = {
    "input/block_2x2x4mm.msh": sha256(ROOT / "input" / "block_2x2x4mm.msh"),
    "input/material.dat": sha256(ROOT / "input" / "material.dat"),
    "input/generated_geometry.json": sha256(ROOT / "input" / "generated_geometry.json"),
}
source_files = {
    "generate_inputs.py": sha256(ROOT / "generate_inputs.py"),
    "preflight.py": sha256(ROOT / "preflight.py"),
    "preflight.json": sha256(ROOT / "preflight.json"),
    "run_contract.py": sha256(ROOT / "run_contract.py"),
}
wheel_path = ROOT.parents[2] / "third_party" / "akantu-runtime" / "downloads" / (
    "akantu-5.0.7.post1-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl"
)

rho = float(geometry["material"]["rho_kg_m3"])
volume = float(geometry["bulk_volume_m3"])
interface_area = float(geometry["interface_area_m2"])
gc = float(geometry["material"]["Gc_J_m2"])
manifest = {
    "schema": "cohesive-capability-manifest-v1",
    "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
    "purpose": "Akantu extrinsic cohesive FEM tiny capability contract; not material calibration or free crack-path validation",
    "backend": {
        "distribution": "akantu",
        "version": preflight["akantu_version"],
        "wheel_filename": wheel_path.name,
        "wheel_sha256": sha256(wheel_path),
        "wheel_source": "https://pypi.org/pypi/akantu/json",
        "python_version": preflight["python_version"],
        "python_executable": preflight["python_executable"],
        "platform": preflight["platform"],
        "numpy_version": preflight["numpy_version"],
        "scipy_version": preflight["scipy_version"],
        "runtime_scope": "third_party/akantu-runtime isolated under WSL Ubuntu; no global/system install",
    },
    "geometry": {
        "units": "SI",
        "dimensions_m": geometry["dimensions_m"],
        "bulk_volume_m3": volume,
        "bulk_tetrahedra": geometry["bulk_tetrahedra"],
        "nodes_before_insertion": geometry["nodes_before_insertion"],
        "interface_plane_z_m": 2.0e-3,
        "interface_area_m2": interface_area,
        "interface_triangles": geometry["interface_facets"],
        "only_allowed_physical_insertion_surface": "interface",
        "central_interface_is_api_contract_test_only": True,
    },
    "expected_insertion_facets": int(geometry["interface_facets"]),
    "material": {
        **geometry["material"],
        "status": "UNVALIDATED_TOY",
        "mass_density_volume_expected_kg": rho * volume,
        "expected_energy_if_fully_debonded_J": gc * interface_area,
        "source_note": "illustrative values also shown in official Akantu examples; not ceramic properties",
    },
    "expected_interface_dissipation_J": gc * interface_area,
    "preflight": {
        "stable_time_step_s": stable_dt,
        "initial_cohesive_elements": preflight["cohesive_elements_before_insertion"],
        "initial_energies_J": preflight["initial_energies_J"],
        "mass_shape": preflight["mass_array_shape"],
        "mass_raw_sum_over_all_dof_components_kg": preflight["mass_raw_sum_over_all_dof_components_kg"],
        "mass_component_sums_kg": preflight["mass_component_sums_kg"],
        "expected_mass_rhoV_kg": preflight["expected_physical_body_mass_rhoV_kg"],
        "eligible_interface_facets_checked": 8,
        "all_triangle_facet_slots": 232,
    },
    "loading": {
        "analysis": "explicit lumped mass",
        "extrinsic_insertion": True,
        "zero_load_boundary_conditions": "bottom and top z fixed at 0 m; no imposed displacement",
        "mode": "mode-I by prescribed top z displacement; bottom z fixed",
        "top_displacement_m": 1.0e-5,
        "ramp_duration_s": endtime,
        "endtime_s": endtime,
        "dt_safety_factor_below_stable": dt_limit_factor,
    },
    "cases": {
        "zero_load": {"dt_s": base_dt, "steps": base_steps, "endtime_s": endtime},
        "modeI_dt": {"dt_s": base_dt, "steps": base_steps, "endtime_s": endtime},
        "modeI_dt2": {"dt_s": fine_dt, "steps": fine_steps, "endtime_s": endtime},
    },
    "limits": {"max_steps": 10_000, "max_wall_time_s_per_case": 300.0, "installation_and_start_budget_min": 15},
    "acceptance": {
        "dynamic_insertions_must_be_nonzero_in_modeI": True,
        "native_connectivity_must_show_bilateral_cohesive_facets": True,
        "parent_material_node_mapping": "only exact native IDs/connectivity; UNKNOWN on conflict or missing proof; never nearest-coordinate matching",
        "complete_debond_energy_relative_error_max": 0.01,
        "dt_refinement_relative_change_max": 0.01,
        "mass_component_sums_retained": True,
        "incomplete_debond_status": "NOT_REACHED",
    },
    "nonclaims": {
        "contact": "NOT_RUN",
        "roughness": "NOT_VALIDATED",
        "free_crack_path": "NOT_VALIDATED; only a central physical interface is enabled",
        "real_ceramic_qualification": "NOT_VALIDATED_TOY",
        "gpu": "UNKNOWN",
    },
    "input_sha256": input_files,
    "source_sha256": source_files,
}
serialized = json.dumps(manifest, indent=2) + "\n"
MANIFEST.write_text(serialized, encoding="utf-8")
(ROOT / "manifest.json").write_text(serialized, encoding="utf-8")
print(serialized, end="")
