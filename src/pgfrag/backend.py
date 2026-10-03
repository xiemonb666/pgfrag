from __future__ import annotations

import csv
import json
import math
import os
import re
import shutil
import subprocess
import time
import tomllib
from pathlib import Path
from typing import Any

from . import BACKEND_ID, __version__
from .assets import validate_asset
from .errors import PGFragError, UnsupportedCapability, ValidationError
from .hashing import file_record, sha256_file
from .jobs import load_job, safe_output_path


OUTPUT_DATA_DICTIONARY: dict[str, Any] = {
    "schema": "pgfrag.output_dictionary.v1",
    "frame": {"name": "right_handed_z_up", "origin": "center of the requested uniform-box bar fixture",
              "reference_current_relation": "x_current = x_reference + displacement in the same local frame"},
    "tables": {
        "points.csv": {
            "meaning": "one solver material point at the final requested time step",
            "columns": {
                "point_id": {"dtype": "int", "unit": "1", "role": "global point key"},
                "x_ref_m": {"dtype": "float64", "unit": "m", "role": "immutable reference/material x coordinate"},
                "y_ref_m": {"dtype": "float64", "unit": "m", "role": "immutable reference/material y coordinate"},
                "z_ref_m": {"dtype": "float64", "unit": "m", "role": "immutable reference/material z coordinate"},
                "x_current_m": {"dtype": "float64", "unit": "m", "role": "endpoint x coordinate in the same frame"},
                "y_current_m": {"dtype": "float64", "unit": "m", "role": "endpoint y coordinate in the same frame"},
                "z_current_m": {"dtype": "float64", "unit": "m", "role": "endpoint z coordinate in the same frame"},
                "vx_mps": {"dtype": "float64", "unit": "m/s", "role": "endpoint x velocity"},
                "vy_mps": {"dtype": "float64", "unit": "m/s", "role": "endpoint y velocity"},
                "vz_mps": {"dtype": "float64", "unit": "m/s", "role": "endpoint z velocity"},
                "damage": {"dtype": "float64", "unit": "1", "role": "upstream material-model point damage in [0,1]"},
                "n_active_bonds": {"dtype": "int", "unit": "1", "role": "upstream active neighbor count for this point"},
                "strain_energy_density_J_m3": {"dtype": "float64", "unit": "J/m^3", "role": "computed by upstream Peridynamics export_field"},
                "volume_m3": {"dtype": "float64", "unit": "m^3", "role": "solver point volume"},
                "mass_kg": {"dtype": "float64", "unit": "kg", "role": "rho times solver point volume"},
                "intact_bond_component_id": {"dtype": "int", "unit": "1", "role": "endpoint intact horizon-bond graph component; not a verified rigid fragment"},
            },
        },
        "bond_states.csv": {
            "meaning": "one canonical undirected pair for each initial native horizon bond, with its final activity state",
            "columns": {
                "point_i_global": {"dtype": "int", "unit": "1", "role": "lower global endpoint ID"},
                "point_j_global": {"dtype": "int", "unit": "1", "role": "higher global endpoint ID"},
                "active": {"dtype": "bool", "unit": "1", "role": "common state of the two reciprocal native directed bonds"},
                "reference_length_m": {"dtype": "float64", "unit": "m", "role": "initial native horizon-bond length"},
                "current_length_m": {"dtype": "float64", "unit": "m", "role": "endpoint pair length"},
            },
        },
        "broken_bonds.csv": {
            "meaning": "exact inactive-pair subset of bond_states.csv at the endpoint; no event time or surface is inferred",
            "columns": {
                "point_i_global": {"dtype": "int", "unit": "1", "role": "lower global endpoint ID"},
                "point_j_global": {"dtype": "int", "unit": "1", "role": "higher global endpoint ID"},
                "reference_length_m": {"dtype": "float64", "unit": "m", "role": "reference pair length"},
                "current_length_m": {"dtype": "float64", "unit": "m", "role": "endpoint pair length"},
            },
        },
        "components.csv": {
            "meaning": "per-point endpoint connected-component assignment in the intact native bond graph",
            "columns": {
                "point_id": {"dtype": "int", "unit": "1", "role": "global point key"},
                "component_id": {"dtype": "int", "unit": "1", "role": "endpoint intact-bond graph component label"},
                "point_mass_kg": {"dtype": "float64", "unit": "kg", "role": "point mass"},
                "damage": {"dtype": "float64", "unit": "1", "role": "upstream endpoint point damage"},
            },
        },
        "isolated_points.csv": {
            "meaning": "points of degree zero in the endpoint intact bond graph",
            "columns": {
                "point_id": {"dtype": "int", "unit": "1", "role": "global point key"},
                "x_ref_m": {"dtype": "float64", "unit": "m", "role": "immutable reference/material x coordinate"},
                "y_ref_m": {"dtype": "float64", "unit": "m", "role": "immutable reference/material y coordinate"},
                "z_ref_m": {"dtype": "float64", "unit": "m", "role": "immutable reference/material z coordinate"},
                "mass_kg": {"dtype": "float64", "unit": "kg", "role": "point mass"},
            },
        },
        "summary.toml": {
            "meaning": "solver provenance, actual requested/delivered discretization, endpoint metrics, and capability UNKNOWNs",
            "units": {"position": "m", "dx_m|horizon_m": "m", "point_center_span_mm|implied_full_point_cell_support_mm": "mm",
                      "point_volume_sum_mm3": "mm^3", "rho_kg_m3": "kg/m^3", "E_pa": "Pa", "Gc_j_m2": "J/m^2",
                      "dt_s|accepted_dt_s|native_stability_bound_dt_s|time_s": "s",
                      "kinetic_energy_j|strain_energy_j": "J", "wall_time_s": "s"},
        },
    },
    "absent_fields": {
        "crack_surface_mesh": "UNSUPPORTED: no geometric surface reconstruction",
        "paired_fracture_surface_lineage_or_historical_surface_adjacency": "UNKNOWN: endpoint horizon bonds do not provide positive-area paired surfaces",
        "per_step_fracture_event_timeline": "UNKNOWN: only endpoint bond state is emitted",
        "contact_state_and_fragment_pose": "UNKNOWN: no contact fixture was run",
        "GPU_state": "UNSUPPORTED/NOT_IMPLEMENTED",
    },
}

# Keep frame/null semantics at the column level so downstream readers need not
# infer coordinate meaning from CSV header spellings. Non-spatial scalars use
# frame=None; all emitted table cells are required and validated as finite.
for _table in OUTPUT_DATA_DICTIONARY["tables"].values():
    for _column_name, _column in _table.get("columns", {}).items():
        _column["frame"] = (
            "right_handed_z_up"
            if _column_name.startswith(("x_", "y_", "z_", "vx_", "vy_", "vz_"))
            else None
        )
        _column["missing"] = "not emitted as null; incomplete rows fail result validation"


def _atomic_json(path: Path, value: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def _project_root() -> Path:
    configured = os.environ.get("PGFRAG_RUNTIME_ROOT")
    if configured:
        root = Path(configured).expanduser().resolve(strict=True)
        if not all((root / name).is_file() for name in ("Project.toml", "Manifest.toml", "run_fixture.jl", "backend_provenance.json")):
            raise ValidationError("PGFRAG_RUNTIME_ROOT must be the PG-Frag Julia environment directory")
        return root
    package_resources = Path(__file__).resolve().parent / "resources" / "julia"
    if all((package_resources / name).is_file() for name in ("Project.toml", "Manifest.toml", "run_fixture.jl", "backend_provenance.json")):
        return package_resources
    checkout = Path(__file__).resolve().parents[2]
    if all((checkout / name).is_file() for name in ("Project.toml", "Manifest.toml", "run_fixture.jl", "backend_provenance.json")):
        return checkout
    raise ValidationError("PG-Frag Julia project resources are unavailable; install the wheel extras or set PGFRAG_RUNTIME_ROOT")


def _find_julia(project_root: Path) -> Path | None:
    configured = os.environ.get("PGFRAG_JULIA_EXE")
    if configured:
        candidate = Path(configured).expanduser()
        return candidate.resolve() if candidate.is_file() else None
    suffix = ".exe" if os.name == "nt" else ""
    search_roots: list[Path] = []
    checkout = Path(__file__).resolve().parents[3]
    search_roots.append(checkout / "third_party" / "pgfrag-runtime" / "extract" / "julia-1.10.12")
    for root in search_roots:
        candidate = root / "bin" / f"julia{suffix}"
        if candidate.is_file():
            return candidate.resolve()
    executable = shutil.which("julia")
    return Path(executable).resolve() if executable else None


def _isolated_depot(project_root: Path) -> Path:
    configured = os.environ.get("PGFRAG_JULIA_DEPOT_PATH")
    if configured:
        depot = Path(configured).expanduser().resolve()
    else:
        checkout_depot = Path(__file__).resolve().parents[3] / "third_party" / "pgfrag-runtime" / "depot"
        if checkout_depot.is_dir():
            depot = checkout_depot.resolve()
        else:
            if os.name == "nt":
                base = Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData" / "Local"))
            else:
                base = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))
            depot = base / "pgfrag" / "julia-depot"
    depot.mkdir(parents=True, exist_ok=True)
    return depot


def _verify_backend_provenance(project_root: Path) -> dict[str, Any]:
    provenance_path = project_root / "backend_provenance.json"
    try:
        provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
        project = tomllib.loads((project_root / "Project.toml").read_text(encoding="utf-8"))
        manifest = tomllib.loads((project_root / "Manifest.toml").read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, tomllib.TOMLDecodeError) as exc:
        raise ValidationError(f"invalid Julia backend provenance/project files: {exc}") from exc
    peridynamics_entry = manifest.get("deps", {}).get("Peridynamics")
    if isinstance(peridynamics_entry, list):
        peridynamics_entry = peridynamics_entry[0] if len(peridynamics_entry) == 1 else None
    if not isinstance(peridynamics_entry, dict) or peridynamics_entry.get("path") != "vendor/Peridynamics.jl":
        raise ValidationError("Manifest.toml must resolve Peridynamics from the relative vendor/Peridynamics.jl path")
    if peridynamics_entry.get("version") != provenance.get("backend_version"):
        raise ValidationError("Manifest Peridynamics version disagrees with backend provenance")
    changed = []
    verified_records = provenance.get("mechanics_source_files", []) + provenance.get("headless_changed_source_files", [])
    for record in verified_records:
        path = (project_root / record["path"]).resolve(strict=True)
        if project_root.resolve(strict=True) not in path.parents:
            raise ValidationError("backend source provenance path escapes the configured project")
        actual = sha256_file(path)
        expected = record.get("sha256", record.get("fork_sha256", ""))
        if actual.lower() != expected.lower():
            changed.append(record["path"])
    if changed:
        raise ValidationError(f"vendored mechanics/time integration sources differ from the pinned audit: {changed}")
    patch = project_root / provenance.get("headless_patch_path", "")
    if not patch.is_file() or sha256_file(patch).lower() != str(provenance.get("headless_patch_sha256", "")).lower():
        raise ValidationError("headless patch file is missing or its SHA-256 does not match backend provenance")
    license_path = project_root / provenance.get("vendored_license_path", "")
    if not license_path.is_file() or sha256_file(license_path).lower() != str(provenance.get("vendored_license_sha256", "")).lower():
        raise ValidationError("vendored upstream LICENSE file is missing or its SHA-256 does not match provenance")
    runner_record = provenance.get("adapter_runner", {})
    runner_path = (project_root / runner_record.get("path", "")).resolve(strict=True)
    if project_root.resolve(strict=True) not in runner_path.parents:
        raise ValidationError("fixture runner provenance path escapes the configured project")
    if sha256_file(runner_path).lower() != str(runner_record.get("sha256", "")).lower():
        raise ValidationError("fixture runner source differs from its recorded SHA-256")
    return {"upstream_tag": provenance.get("upstream_tag"), "backend_version": provenance.get("backend_version"),
            "upstream_source_archive_sha256": provenance.get("upstream_source_archive_sha256"),
            "mechanics_source_files_verified": len(provenance.get("mechanics_source_files", [])),
            "headless_io_source_files_verified": len(provenance.get("headless_changed_source_files", [])),
            "headless_patch_sha256": provenance.get("headless_patch_sha256"),
            "vendored_license_sha256": provenance.get("vendored_license_sha256"),
            "adapter_runner_sha256": runner_record.get("sha256"),
            "source_tree": "physics/** and time_solvers/** match the recorded upstream 0.5.4 source-file hashes"}


def _julia_version(executable: Path) -> str:
    try:
        completed = subprocess.run([str(executable), "--version"], capture_output=True, text=True, timeout=15,
                                   check=False, encoding="utf-8", errors="replace")
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise UnsupportedCapability(f"could not query the configured Julia runtime: {exc}") from exc
    version_text = (completed.stdout + completed.stderr).strip()
    match = re.search(r"julia version (\d+)\.(\d+)\.(\d+)", version_text, re.IGNORECASE)
    if completed.returncode != 0 or not match:
        raise UnsupportedCapability(f"configured Julia runtime did not report a semantic version: {version_text!r}")
    version = tuple(int(value) for value in match.groups())
    if version[0:2] != (1, 10):
        raise UnsupportedCapability(f"this locked Peridynamics environment requires Julia 1.10.x, got {version_text!r}")
    return f"{version[0]}.{version[1]}.{version[2]}"


def _under_root(root: Path, relative: str, *, must_exist: bool) -> Path:
    candidate = root.joinpath(*Path(relative.replace("\\", "/")).parts)
    resolved = candidate.resolve(strict=must_exist)
    root_resolved = root.resolve(strict=True)
    if resolved != root_resolved and root_resolved not in resolved.parents:
        raise ValidationError(f"package path escapes the job directory: {relative!r}")
    return resolved


def _load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _finite_float(row: dict[str, str], key: str, label: str) -> float:
    try:
        value = float(row[key])
    except (KeyError, ValueError) as exc:
        raise ValidationError(f"missing or invalid numeric field {label}.{key}") from exc
    if not math.isfinite(value):
        raise ValidationError(f"non-finite numeric field {label}.{key}")
    return value


def _validate_solver_outputs(output_dir: Path, summary: dict[str, Any], job: dict[str, Any]) -> dict[str, Any]:
    points = _load_csv(output_dir / "points.csv")
    bonds = _load_csv(output_dir / "bond_states.csv")
    broken = _load_csv(output_dir / "broken_bonds.csv")
    components = _load_csv(output_dir / "components.csv")
    isolated = _load_csv(output_dir / "isolated_points.csv")
    n = summary.get("n_points")
    steps = summary.get("n_steps")
    if isinstance(n, bool) or not isinstance(n, int) or n <= 0 or len(points) != n:
        raise ValidationError("point row count does not match summary n_points")
    if isinstance(steps, bool) or not isinstance(steps, int) or steps != job["solver"]["steps"]:
        raise ValidationError("reported step count does not match requested solver.steps")
    dt = float(summary["dt_s"])
    time_s = float(summary["time_s"])
    if not math.isfinite(dt) or dt <= 0 or not math.isfinite(time_s) or time_s <= 0 or not math.isclose(time_s, dt * steps, rel_tol=1e-11, abs_tol=1e-15):
        raise ValidationError("reported time step/duration is non-finite or inconsistent")
    if job["solver"]["dt_s"] is not None and not math.isclose(dt, job["solver"]["dt_s"], rel_tol=1e-10, abs_tol=1e-15):
        raise ValidationError("reported dt_s does not match requested solver.dt_s")
    for key in ("max_damage", "mean_damage", "strain_energy_j", "kinetic_energy_j", "point_mass_kg", "point_volume_sum_mm3", "wall_time_s",
                "native_stability_bound_dt_s", "accepted_dt_s"):
        value = summary.get(key)
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
            raise ValidationError(f"summary field {key} is not finite numeric data")
    if min(summary["strain_energy_j"], summary["kinetic_energy_j"], summary["wall_time_s"]) < 0:
        raise ValidationError("summary energy and wall time fields must be nonnegative")
    if summary["native_stability_bound_dt_s"] <= 0 or summary["accepted_dt_s"] <= 0:
        raise ValidationError("native time-step stability bound and accepted dt must be positive")
    requested_dt_mode = summary.get("requested_dt_mode")
    stability_check = summary.get("native_stability_check")
    expected_scope = "single-body fixture only; not a general OSB/contact stability guarantee"
    if summary.get("native_stability_scope") != expected_scope:
        raise ValidationError("native stability scope is missing or overclaims solver generality")
    tolerance = 128 * math.ulp(max(summary["native_stability_bound_dt_s"], summary["accepted_dt_s"]))
    if requested_dt_mode == "AUTO":
        if stability_check != "PASS_NATIVE_AUTO_DT" or not math.isclose(
            summary["accepted_dt_s"], summary["native_stability_bound_dt_s"], rel_tol=1e-12, abs_tol=tolerance
        ):
            raise ValidationError("automatic dt does not match the upstream native stable time-step bound")
    elif requested_dt_mode == "EXPLICIT":
        requested_dt = job["solver"].get("dt_s")
        if stability_check != "PASS_EXPLICIT_DT_WITHIN_NATIVE_BOUND" or requested_dt is None:
            raise ValidationError("explicit dt lacks a successful native bound preflight record")
        if not math.isclose(summary["accepted_dt_s"], requested_dt, rel_tol=1e-12, abs_tol=tolerance):
            raise ValidationError("accepted explicit dt disagrees with the requested job value")
        if summary["accepted_dt_s"] > summary["native_stability_bound_dt_s"] + tolerance:
            raise ValidationError("explicit dt exceeds the recorded upstream native stability bound")
    else:
        raise ValidationError("requested_dt_mode must be AUTO or EXPLICIT")

    point_ids: list[int] = []
    point_by_id: dict[int, dict[str, str]] = {}
    volume_sum = mass_sum = strain_energy = kinetic_energy = 0.0
    for row in points:
        try:
            point_id = int(row["point_id"])
            active_count = int(row["n_active_bonds"])
            component_id = int(row["intact_bond_component_id"])
        except (KeyError, ValueError) as exc:
            raise ValidationError("point state identifiers/counts are invalid") from exc
        if point_id <= 0 or active_count < 0 or component_id <= 0:
            raise ValidationError("point IDs, active bond counts and component IDs must be positive/nonnegative")
        point_ids.append(point_id)
        point_by_id[point_id] = row
        for field in ("x_ref_m", "y_ref_m", "z_ref_m", "x_current_m", "y_current_m", "z_current_m",
                      "vx_mps", "vy_mps", "vz_mps", "damage", "strain_energy_density_J_m3", "volume_m3", "mass_kg"):
            _finite_float(row, field, f"point[{point_id}]")
        damage = float(row["damage"])
        sed = float(row["strain_energy_density_J_m3"])
        volume = float(row["volume_m3"])
        mass = float(row["mass_kg"])
        if not 0.0 <= damage <= 1.0 or sed < 0 or volume <= 0 or mass <= 0:
            raise ValidationError(f"point[{point_id}] has damage/energy/volume/mass outside valid ranges")
        expected_mass = job["material"]["rho_kg_m3"] * volume
        if not math.isclose(mass, expected_mass, rel_tol=1e-10, abs_tol=1e-18):
            raise ValidationError(f"point[{point_id}] mass does not match rho times volume")
        volume_sum += volume
        mass_sum += mass
        strain_energy += sed * volume
        kinetic_energy += 0.5 * mass * sum(float(row[field]) ** 2 for field in ("vx_mps", "vy_mps", "vz_mps"))
    if len(set(point_ids)) != n or set(point_ids) != set(range(1, n + 1)):
        raise ValidationError("point IDs are duplicated, incomplete, or outside the global 1..n range")
    if not math.isclose(mass_sum, summary["point_mass_kg"], rel_tol=1e-10, abs_tol=1e-15):
        raise ValidationError("point mass rows do not match summary mass")
    if not math.isclose(volume_sum * 1e9, summary["point_volume_sum_mm3"], rel_tol=1e-10, abs_tol=1e-8):
        raise ValidationError("point volume rows do not match summary volume")
    if not math.isclose(strain_energy, summary["strain_energy_j"], rel_tol=1e-9, abs_tol=1e-18):
        raise ValidationError("serialized strain energy does not match point SED integration")
    if not math.isclose(kinetic_energy, summary["kinetic_energy_j"], rel_tol=1e-9, abs_tol=1e-18):
        raise ValidationError("serialized kinetic energy does not match point state integration")

    edge_states: dict[tuple[int, int], bool] = {}
    adjacency = {point_id: set() for point_id in point_ids}
    for row in bonds:
        try:
            i, j = int(row["point_i_global"]), int(row["point_j_global"])
        except (KeyError, ValueError) as exc:
            raise ValidationError("bond endpoint IDs are invalid") from exc
        if not (1 <= i <= n and 1 <= j <= n and i < j):
            raise ValidationError(f"bond pair {(i, j)} is noncanonical or outside the point set")
        key = (i, j)
        if key in edge_states:
            raise ValidationError(f"duplicate undirected bond pair: {key}")
        active_text = row.get("active", "").lower()
        if active_text not in {"true", "false"}:
            raise ValidationError(f"bond pair {key} has invalid active state")
        active = active_text == "true"
        edge_states[key] = active
        for field in ("reference_length_m", "current_length_m"):
            if _finite_float(row, field, f"bond{key}") <= 0:
                raise ValidationError(f"bond pair {key} has a nonpositive length")
        if active:
            adjacency[i].add(j)
            adjacency[j].add(i)
    if not edge_states:
        raise ValidationError("solver emitted no native initial bond pairs")
    broken_rows: dict[tuple[int, int], dict[str, str]] = {}
    for row in broken:
        key = (int(row["point_i_global"]), int(row["point_j_global"]))
        if key in broken_rows:
            raise ValidationError(f"duplicate broken bond pair {key}")
        broken_rows[key] = row
    expected_broken = {key for key, active in edge_states.items() if not active}
    if set(broken_rows) != expected_broken:
        raise ValidationError("broken_bonds.csv does not equal the inactive endpoint bond-state subset")

    parent = {point_id: point_id for point_id in point_ids}
    def root(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value
    for (i, j), active in edge_states.items():
        if active:
            ri, rj = root(i), root(j)
            if ri != rj:
                parent[max(ri, rj)] = min(ri, rj)
    root_to_component: dict[int, set[int]] = {}
    for point_id in point_ids:
        root_to_component.setdefault(root(point_id), set()).add(point_id)
    solver_groups: dict[int, set[int]] = {}
    component_mass_rows: dict[int, float] = {}
    component_damage_rows: dict[int, float] = {}
    if len(components) != n:
        raise ValidationError("components.csv must contain one row per point")
    for row in components:
        try:
            point_id = int(row["point_id"])
            component_id = int(row["component_id"])
        except (KeyError, ValueError) as exc:
            raise ValidationError("component row identifiers are invalid") from exc
        if point_id not in point_by_id or point_id in component_mass_rows or component_id <= 0:
            raise ValidationError("component rows have duplicate/unknown point IDs or invalid labels")
        mass = _finite_float(row, "point_mass_kg", f"component[{point_id}]")
        damage = _finite_float(row, "damage", f"component[{point_id}]")
        if not math.isclose(mass, float(point_by_id[point_id]["mass_kg"]), rel_tol=1e-10, abs_tol=1e-18):
            raise ValidationError(f"component row mass disagrees for point {point_id}")
        if not math.isclose(damage, float(point_by_id[point_id]["damage"]), rel_tol=1e-10, abs_tol=1e-12):
            raise ValidationError(f"component row damage disagrees for point {point_id}")
        if int(point_by_id[point_id]["intact_bond_component_id"]) != component_id:
            raise ValidationError(f"point/component CSV labels disagree for point {point_id}")
        solver_groups.setdefault(component_id, set()).add(point_id)
        component_mass_rows[point_id] = mass
        component_damage_rows[point_id] = damage
    graph_partition = {frozenset(group) for group in root_to_component.values()}
    solver_partition = {frozenset(group) for group in solver_groups.values()}
    if graph_partition != solver_partition:
        raise ValidationError("emitted component labels disagree with recomputed intact-bond graph connectivity")
    isolated_expected = {point_id for point_id, neighbors in adjacency.items() if not neighbors}
    isolated_ids = {int(row["point_id"]) for row in isolated}
    if isolated_ids != isolated_expected:
        raise ValidationError("isolated_points.csv does not match intact-bond graph degree zero")
    if summary.get("intact_bond_connected_component_count") != len(graph_partition):
        raise ValidationError("summary component count disagrees with graph and component map")
    return {
        "point_ids_complete_unique": True,
        "point_state_finite_ranges": True,
        "point_mass_volume_consistent": True,
        "energy_integrals_match_point_state": True,
        "native_edges_unique_valid": True,
        "broken_edges_equal_inactive_subset": True,
        "component_partition_matches_graph": True,
        "isolated_points_match_graph": True,
        "steps_and_duration_consistent": True,
        "native_time_step_bound_validated": True,
        "point_count": n,
        "bond_pair_count": len(edge_states),
        "broken_bond_pair_count": len(expected_broken),
        "component_count": len(graph_partition),
    }


def _write_result(output_dir: Path, result: dict[str, Any]) -> None:
    _atomic_json(output_dir / "result.json", result)


def _write_data_dictionary(output_dir: Path) -> None:
    _atomic_json(output_dir / "data_dictionary.json", OUTPUT_DATA_DICTIONARY)


def _run_unsupported_asset(job: dict[str, Any], job_path: Path, output_dir: Path,
                           job_digest: str) -> dict[str, Any]:
    initial = job["initial_condition"]
    package = _under_root(job_path.parent, initial["package_dir"], must_exist=True)
    audit = validate_asset(package)
    if audit["manifest_sha256"].lower() != initial["manifest_sha256"].lower():
        raise ValidationError("input asset manifest SHA-256 does not match the job")
    if audit["declared_geometry_validation"]["status"] != "PASS":
        return {
            "schema": "pgfrag.result.v1", "job_id": job["job_id"], "status": "REJECTED",
            "backend_id": BACKEND_ID, "core_version": __version__,
            "input": {"job_sha256": job_digest, "asset_manifest_sha256": audit["manifest_sha256"]},
            "validation": {"status": "FAIL", "reason": "asset declared geometry gate is not PASS", "asset_audit": audit},
            "missing_capabilities": [],
        }
    return {
        "schema": "pgfrag.result.v1", "job_id": job["job_id"], "status": "UNSUPPORTED",
        "backend_id": BACKEND_ID, "core_version": __version__,
        "input": {"job_sha256": job_digest, "asset_manifest_sha256": audit["manifest_sha256"]},
        "validation": {"status": "NOT_RUN", "asset_contract": "PASS", "asset_audit": audit},
        "missing_capabilities": [{
            "capability": "forge_asset_surface_to_volume_discretization",
            "status": "UNSUPPORTED",
            "reason": "surface-to-volume material point discretization and void/topology checks are not implemented; no solver was launched",
        }],
        "outputs": [],
    }


def run_job(job_file: str | Path) -> tuple[int, dict[str, Any]]:
    job, job_path = load_job(job_file)
    output_dir = safe_output_path(job_path, job["output_dir"])
    job_digest = sha256_file(job_path)
    if output_dir.exists() and any(output_dir.iterdir()):
        raise ValidationError(f"output_dir already exists and is non-empty: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    try:
        project_root = _project_root()
    except PGFragError as exc:
        result = {"schema": "pgfrag.result.v1", "job_id": job["job_id"], "status": "UNSUPPORTED",
                  "backend_id": BACKEND_ID, "core_version": __version__,
                  "input": {"job_sha256": job_digest},
                  "validation": {"status": "NOT_RUN", "reason": str(exc)},
                  "missing_capabilities": [{"capability": "locked_julia_project_resources", "status": "UNSUPPORTED", "reason": str(exc)}]}
        _write_result(output_dir, result)
        return 3, result
    runner = project_root / "run_fixture.jl"
    result_path = output_dir / "result.json"

    if job["initial_condition"]["kind"] == "forge_asset_surface":
        try:
            result = _run_unsupported_asset(job, job_path, output_dir, job_digest)
        except PGFragError as exc:
            result = {
                "schema": "pgfrag.result.v1", "job_id": job["job_id"], "status": "REJECTED",
                "backend_id": BACKEND_ID, "core_version": __version__,
                "input": {"job_sha256": job_digest},
                "validation": {"status": "FAIL", "reason": str(exc)}, "missing_capabilities": [],
            }
        _write_result(output_dir, result)
        code = 3 if result["status"] == "UNSUPPORTED" else (2 if result["status"] == "REJECTED" else 0)
        return code, result

    julia = _find_julia(project_root)
    if julia is None:
        result = {
            "schema": "pgfrag.result.v1", "job_id": job["job_id"], "status": "UNSUPPORTED",
            "backend_id": BACKEND_ID, "core_version": __version__,
            "input": {"job_sha256": job_digest},
            "validation": {"status": "NOT_RUN", "reason": "Julia 1.10+ runtime was not found; set PGFRAG_JULIA_EXE"},
            "missing_capabilities": [{"capability": "peridynamics_jl_runtime", "status": "UNSUPPORTED",
                                       "reason": "no solver process was launched"}],
        }
        _write_result(output_dir, result)
        return 3, result

    try:
        runtime_version = _julia_version(julia)
        source_provenance = _verify_backend_provenance(project_root)
        depot = _isolated_depot(project_root)
    except UnsupportedCapability as exc:
        result = {"schema": "pgfrag.result.v1", "job_id": job["job_id"], "status": "UNSUPPORTED",
                  "backend_id": BACKEND_ID, "core_version": __version__,
                  "input": {"job_sha256": job_digest},
                  "validation": {"status": "NOT_RUN", "reason": str(exc)},
                  "missing_capabilities": [{"capability": "compatible_julia_runtime", "status": "UNSUPPORTED", "reason": str(exc)}]}
        _write_result(output_dir, result)
        return 3, result
    except PGFragError as exc:
        result = {"schema": "pgfrag.result.v1", "job_id": job["job_id"], "status": "REJECTED",
                  "backend_id": BACKEND_ID, "core_version": __version__,
                  "input": {"job_sha256": job_digest},
                  "validation": {"status": "FAIL", "reason": str(exc)}, "missing_capabilities": []}
        _write_result(output_dir, result)
        return 2, result

    if not runner.is_file():
        raise ValidationError(f"Julia fixture runner missing: {runner}")
    initial = job["initial_condition"]
    material = job["material"]
    loading = job["loading"]
    solver = job["solver"]
    resources = job["resources"]
    dimensions = [str(v) for v in initial["dimensions_mm"]]
    args = [
        str(julia), "--startup-file=no", "--compiled-modules=no", f"--project={project_root}",
        f"--threads={resources['threads']}", str(runner), str(output_dir),
        *dimensions, str(initial["dx_mm"]), str(initial["horizon_mm"]),
        material["model"], str(material["rho_kg_m3"]), str(material["E_pa"]), str(material["nu"]),
        str(material["Gc_j_m2"]), loading["kind"], str(loading["speed_m_s"]),
        str(solver["steps"]), "AUTO" if solver["dt_s"] is None else str(solver["dt_s"]),
        str(solver["safety_factor"]), str(resources["max_points"]),
    ]
    env = os.environ.copy()
    env["JULIA_NUM_THREADS"] = str(resources["threads"])
    env["JULIA_DEPOT_PATH"] = str(depot)
    env.setdefault("JULIA_NUM_PRECOMPILE_TASKS", "1")
    start = time.monotonic()
    process = subprocess.Popen(args, cwd=project_root, env=env, stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE, text=True, encoding="utf-8", errors="replace")
    try:
        stdout, stderr = process.communicate(timeout=resources["timeout_seconds"])
        status = "COMPLETED" if process.returncode == 0 else "FAILED"
        reason = None if process.returncode == 0 else "solver process returned a nonzero exit code"
    except subprocess.TimeoutExpired:
        process.terminate()
        try:
            stdout, stderr = process.communicate(timeout=5)
        except subprocess.TimeoutExpired:
            process.kill()
            stdout, stderr = process.communicate()
        status, reason = "FAILED", "solver process exceeded timeout_seconds and was terminated"
    except KeyboardInterrupt:
        process.terminate()
        try:
            stdout, stderr = process.communicate(timeout=5)
        except subprocess.TimeoutExpired:
            process.kill()
            stdout, stderr = process.communicate()
        status, reason = "CANCELLED", "caller interrupted the solver process"
    elapsed = time.monotonic() - start
    (output_dir / "solver.stdout.log").write_text(stdout or "", encoding="utf-8")
    (output_dir / "solver.stderr.log").write_text(stderr or "", encoding="utf-8")
    summary_path = output_dir / "summary.toml"
    base_result: dict[str, Any] = {
        "schema": "pgfrag.result.v1", "job_id": job["job_id"], "status": status,
        "backend_id": BACKEND_ID, "backend_version": "0.5.4+pgfrag.1", "core_version": __version__,
        "input": {"job_sha256": job_digest, "job_path": str(job_path)},
        "process": {"returncode": process.returncode, "wall_time_s": elapsed,
                    "timeout_seconds": resources["timeout_seconds"], "cancel_supported": True,
                    "command_executable": str(julia), "julia_version": runtime_version,
                    "julia_depot_path": str(depot), "project_path": str(project_root)},
        "backend_provenance": source_provenance,
        "validation": {"status": "NOT_RUN" if status != "COMPLETED" else "PENDING",
                       "reason": reason},
        "capabilities": {"cpu": "USED", "gpu": "UNSUPPORTED/NOT_IMPLEMENTED",
                         "material_model": material["model"], "calibration_status": "UNVALIDATED_PROTOTYPE"},
        "missing_capabilities": [
            {"capability": "surface_geometry_and_crack_faces", "status": "UNKNOWN",
             "reason": "uniform_box bar fixture emits material points only; no geometry surface reconstruction"},
            {"capability": "forge_asset_discretization", "status": "UNSUPPORTED",
             "reason": "surface-to-volume discretization and void/topology checks are not implemented"},
            {"capability": "contact", "status": "UNKNOWN", "reason": "no contact fixture was run"},
            {"capability": "paired_fracture_surface_lineage_or_historical_surface_adjacency", "status": "UNKNOWN",
             "reason": "endpoint broken horizon bonds do not provide positive-area paired-surface mapping or historical adjacency; this neutral core does not derive upper-layer relation labels"},
            {"capability": "fracture_event_timeline", "status": "UNKNOWN",
             "reason": "this result stores an endpoint bond state, not per-step event chronology"},
            {"capability": "energy_conservation", "status": "UNKNOWN",
             "reason": "prescribed velocity boundaries do external work; no complete energy ledger is emitted"},
        ],
    }
    if status != "COMPLETED":
        base_result["validation"]["returncode"] = process.returncode
        if process.returncode is None:
            base_result["validation"]["reason"] = reason
        base_result["outputs"] = [file_record(path, output_dir) for path in sorted(output_dir.iterdir()) if path.is_file()]
        _write_result(output_dir, base_result)
        return (130 if status == "CANCELLED" else 4), base_result

    if not summary_path.is_file():
        base_result["status"] = "FAILED"
        base_result["validation"] = {"status": "FAIL", "reason": "solver exited 0 without summary.toml", "returncode": process.returncode}
        _write_result(output_dir, base_result)
        return 4, base_result

    try:
        summary = tomllib.loads(summary_path.read_text(encoding="utf-8"))
        checks = _validate_solver_outputs(output_dir, summary, job)
        bonds = _load_csv(output_dir / "bond_states.csv")
        broken = _load_csv(output_dir / "broken_bonds.csv")
    except (OSError, KeyError, ValueError, tomllib.TOMLDecodeError, ValidationError) as exc:
        base_result["status"] = "FAILED"
        base_result["validation"] = {"status": "FAIL", "reason": f"result consistency check failed: {exc}", "returncode": process.returncode}
        base_result["outputs"] = [file_record(path, output_dir) for path in sorted(output_dir.iterdir()) if path.is_file()]
        _write_result(output_dir, base_result)
        return 4, base_result

    base_result.update({
        "units": {"reference_coordinates": "m", "current_coordinates": "m", "velocity": "m/s",
                  "point_volume": "m^3", "point_mass": "kg", "time": "s", "energy": "J"},
        "frame": {"name": "right_handed_z_up", "origin": "center of requested uniform-box bar fixture",
                  "reference_and_current_coordinates_share_frame": True},
        "initial_condition": {"kind": "uniform_box_bar_fixture", "requested_dimensions_mm": initial["dimensions_mm"],
                              "spacing_mm": initial["dx_mm"], "horizon_mm": initial["horizon_mm"],
                              "solver_model": summary["solver_model"],
                              "discrete_center_span_mm": summary["point_center_span_mm"],
                              "implied_full_point_cell_support_mm": summary["implied_full_point_cell_support_mm"],
                              "point_volume_sum_mm3": summary["point_volume_sum_mm3"]},
        "material": {"rho_kg_m3": material["rho_kg_m3"], "E_pa": material["E_pa"], "nu": material["nu"],
                     "Gc_j_m2": material["Gc_j_m2"], "calibration_status": material["calibration_status"],
                     "provenance": material["provenance"]},
        "loading": loading,
        "solver": {"steps": summary["n_steps"], "dt_s": summary["dt_s"], "time_s": summary["time_s"],
                   "requested_dt_mode": summary["requested_dt_mode"],
                   "native_stability_bound_dt_s": summary["native_stability_bound_dt_s"],
                   "accepted_dt_s": summary["accepted_dt_s"],
                   "native_stability_check": summary["native_stability_check"],
                   "native_stability_scope": summary["native_stability_scope"],
                   "max_damage": summary["max_damage"], "mean_damage": summary["mean_damage"],
                   "strain_energy_j": summary["strain_energy_j"], "kinetic_energy_j": summary["kinetic_energy_j"],
                   "total_energy_conservation": "UNKNOWN"},
        "state": {"point_count": summary["n_points"], "canonical_initial_bond_pairs": len(bonds),
                  "broken_bond_pairs_at_endpoint": len(broken), "intact_bond_graph_component_count": summary["intact_bond_connected_component_count"],
                  "point_ids_complete_unique": checks["point_ids_complete_unique"],
                  "broken_state_consistent": checks["broken_edges_equal_inactive_subset"]},
        "validation": {"status": "PASS", "scope": "bounded endpoint serialization and solver process completion only",
                       "physical_material_calibration": "PENDING", "mesh_or_surface_validity": "NOT_APPLICABLE_FIXTURE",
                       "contact_or_fragment_surface_validity": "NOT_TESTED",
                       "time_step_stability": {"status": "PASS", "method": "upstream Peridynamics.calc_stable_timestep",
                                                "accepted_dt_s": summary["accepted_dt_s"],
                                                "native_stability_bound_dt_s": summary["native_stability_bound_dt_s"],
                                                "scope": summary["native_stability_scope"]},
                       "checks": checks},
        "result_path": str(result_path.resolve()),
        "output_hash_scope": "solver artifacts and logs, excluding result.json itself to avoid a recursive hash",
        "data_dictionary_path": "data_dictionary.json",
    })
    _write_data_dictionary(output_dir)
    base_result["outputs"] = [file_record(path, output_dir) for path in sorted(output_dir.iterdir()) if path.is_file()]
    _write_result(output_dir, base_result)
    return 0, base_result
