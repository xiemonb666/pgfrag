import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from pgfrag.bowl_dynamics import (
    BowlGeometry,
    _load_asset_spec,
    _load_raw_states,
    _verify_source_snapshot,
    bond_crosses_inner_void,
    generate_centroid_lattice,
    independent_pair_energy_joule,
    prepare_run,
    segment_minimum_radius_m,
)


def test_centroid_lattice_matches_frozen_bowl_counts_at_both_spacings():
    geometry = BowlGeometry(0.012, 0.008, (0.0, 0.0, 0.012))

    h1 = generate_centroid_lattice(geometry, 0.001)
    h08 = generate_centroid_lattice(geometry, 0.0008)

    assert h1.shape == (2516, 3)
    assert h08.shape == (5052, 3)
    radii = np.linalg.norm(h1 - np.asarray(geometry.center_m), axis=1)
    assert np.all(radii >= geometry.inner_radius_m)
    assert np.all(radii <= geometry.outer_radius_m)


def test_inner_void_filter_uses_closed_segment_closest_point():
    center = (0.0, 0.0, 0.0)
    geometry = BowlGeometry(2.0, 1.0, center)

    assert segment_minimum_radius_m((-1.5, 0.0, 0.0), (1.5, 0.0, 0.0), center) == 0.0
    assert bond_crosses_inner_void((-1.5, 0.0, 0.0), (1.5, 0.0, 0.0), geometry)
    assert segment_minimum_radius_m((-2.0, 1.0, 0.0), (2.0, 1.0, 0.0), center) == 1.0
    assert not bond_crosses_inner_void((-2.0, 1.0, 0.0), (2.0, 1.0, 0.0), geometry)


def test_independent_pair_energy_is_the_frozen_bond_potential():
    x = np.asarray([[0.0, 0.0, 0.0], [1.1, 0.0, 0.0]])
    i, j = np.asarray([0]), np.asarray([1])
    rest, volume = np.asarray([1.0]), np.asarray([0.2, 0.3])
    c = 12.0
    expected = 0.5 * c * volume[0] * volume[1] * rest[0] * 0.1**2

    assert independent_pair_energy_joule(x, i, j, rest, volume, c) == pytest.approx(expected)
    assert independent_pair_energy_joule(x, i, j, rest, volume, c, np.asarray([False])) == 0.0


def test_julia_column_major_bytes_preserve_state_point_and_coordinate_axes(tmp_path: Path):
    # Literal byte order from Julia vec(Array(2,2,3)); this is not a round trip
    # through the Python writer and catches a default-C reshape of native data.
    native_xyz = [1, 2, 10, 20, 3, 4, 30, 40, 5, 6, 50, 60]
    np.asarray(native_xyz, dtype="<f8").tofile(tmp_path / "position.f64le")
    np.asarray(native_xyz, dtype="<f8").tofile(tmp_path / "velocity.f64le")
    np.asarray([.1, .2, .3, .4], dtype="<f4").tofile(tmp_path / "damage.f32le")
    np.asarray([0, 1, 0, 1], dtype="<i4").tofile(tmp_path / "component.i32le")
    np.asarray([11, 12, 13, 14], dtype="<f4").tofile(tmp_path / "native_sed.f32le")
    (tmp_path / "checkpoints.csv").write_text("step,time_s\n", encoding="utf-8")
    states = _load_raw_states(tmp_path, n_steps=1, n_points=2)
    assert np.array_equal(states["position_m"][0], [[1, 3, 5], [10, 30, 50]])
    assert np.array_equal(states["position_m"][1], [[2, 4, 6], [20, 40, 60]])
    assert np.array_equal(states["component_label"], [[0, 0], [1, 1]])
    assert np.array_equal(states["native_sed_J_m3"], [[11, 13], [12, 14]])


def _write_minimal_asset(tmp_path: Path, mechanics: dict[str, float]) -> Path:
    tmp_path.mkdir(parents=True, exist_ok=True)
    payload = tmp_path / "payload.bin"
    payload.write_bytes(b"geometry-payload")
    manifest = {
        "schema": "forge.asset.v1",
        "units": {"length": "mm", "time": "s", "mass": "kg"},
        "frame": "right_handed_z_up",
        "validation": {"scope": "geometry_only", "status": "PASS"},
        "provenance": {
            "spec": {
                "shape": "concentric_lower_hemisphere_bowl",
                "outer_radius_mm": 12.0,
                "inner_radius_mm": 8.0,
                "center_xyz_mm": [0.0, 0.0, 12.0],
            }
        },
        "material": {"calibration_status": "UNVALIDATED_PROTOTYPE", "mechanics": mechanics},
        "files": {
            "payload.bin": {
                "bytes": payload.stat().st_size,
                "sha256": hashlib.sha256(payload.read_bytes()).hexdigest(),
            }
        },
    }
    manifest_path = tmp_path / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, separators=(",", ":")), encoding="utf-8")
    return manifest_path


def test_source_manifest_mechanics_are_exactly_bound_to_frozen_case(tmp_path: Path):
    expected = {"rho_kg_m3": 2300.0, "E_pa": 70.0e9, "nu": 0.25, "Gc_j_m2": 7.6}
    manifest_path = _write_minimal_asset(tmp_path, expected)

    path, _manifest, geometry, digest, mechanics = _load_asset_spec(manifest_path)

    assert path == manifest_path.resolve()
    assert digest == hashlib.sha256(manifest_path.read_bytes()).hexdigest()
    assert geometry.outer_radius_m == pytest.approx(0.012)
    assert mechanics == expected

    changed = dict(expected, E_pa=69.0e9)
    changed_path = _write_minimal_asset(tmp_path / "changed", changed)
    with pytest.raises(ValueError, match="frozen illustrative case"):
        _load_asset_spec(changed_path)


def test_prepared_attempt_pins_manifest_and_snapshot_runner_without_solving(tmp_path: Path):
    repo_root = Path(__file__).resolve().parents[2]
    asset = repo_root / "vesselforge" / "data" / "continuation_20261015" / "spherical_bowl" / "three_bands"
    run_dir = prepare_run(asset, runs_root=tmp_path / "runs", alpha_s_1=0.0)

    source_manifest = (asset / "manifest.json").read_bytes()
    snapshot_manifest = (run_dir / "source_snapshot" / "forge_manifest.json").read_bytes()
    prepared_argv = json.loads((run_dir / "argv.json").read_text(encoding="utf-8"))
    receipts = {entry["path"]: entry for entry in prepared_argv["source_receipts"]}

    assert snapshot_manifest == source_manifest
    assert "source_snapshot/forge_manifest.json" in receipts
    assert "source_snapshot/run_bowl_dynamics.jl" in receipts
    assert Path(prepared_argv["argv"][2]) == run_dir / "source_snapshot" / "run_bowl_dynamics.jl"
    assert _verify_source_snapshot(run_dir)
    assert not (run_dir / "native_status.json").exists()
    assert not (run_dir / "run_receipt.json").exists()
