# Independent static review: prepared point inputs

**Result:** the two prepared 125-point fixtures pass the independent arithmetic and metadata checks. They are internally consistent synthetic voxel quadrature fixtures. This review does not establish eligibility as a real artifact candidate or as a reconstructed specimen.

## Checked inputs and freeze

The audit used the stable generator, two CSV/TOML pairs, and common/case TOMLs listed in `prefreeze_sha256.txt`. All eight SHA-256 values matched the post-audit recheck in `postfreeze_verification.json`; the manifest's corresponding eight pins also matched at audit time.

| File | SHA-256 |
|---|---|
| `make_point_inputs.ps1` | `f218b55964634631fccdaf33ca24649c5c64bc0c2c89c38e7800942fa531b0a0` |
| `input/coupon_uniform_voxel_5cube_125.points.csv` | `9b9207727364c667ecc0b592ccff34405eaab75a0992e057a2d8a7a121b706da` |
| `input/coupon_uniform_voxel_5cube_125.metadata.toml` | `8f6235661f7ba1ff4a8ac03e08be84284b214ce089118f83381abeb5abda10d0` |
| `input/coupon_cell_jitter_5cube_seed17.points.csv` | `2adc7cf3f8a825a6b964cf6e5dc3cba099634054b9a10d60abee8b4590813ee5` |
| `input/coupon_cell_jitter_5cube_seed17.metadata.toml` | `47e4fe7d1d66063b2c234e48e65bae286fc4cdfda8d18b5913a3367da3a92e19` |
| `common_case.toml` | `ca2e4a64946d32f2ba6abdc1bd5b4d1768f5f7587812f3bcce490c0713a0f3bf` |
| `case_zero_load_control.toml` | `a2e6694d047e70ed445ee9037bbcbabbd99095dd060ebe1c78bed674199b2178` |
| `case_tinyload_event_diagnostic.toml` | `53f7aca2d222eb714c1640a27d389aa845669b7cbf26a449f831e4eb307f4a7c` |

The source and manifest were changing during the initial directory inspection, so that earlier observation is retained in `unstable_tree_observation.txt` and `unstable_full_tree_observation_sha256.txt`; it was not mixed into the arithmetic baseline. The single arithmetic audit read manifest SHA-256 `25d703e1ee43b8534ceeee7b86174d5f88b469457983c4387e2f0624e0e36f1b`. Its hash was unchanged immediately before and after execution, its status was `SOURCE_PREPARED_STATIC_REVIEW_ONLY_AWAITING_ROOT_FINAL_FREEZE`, and its fixture declarations and relevant input pins passed. A later recheck observed manifest SHA-256 `50b99e50325f8d44d2d45638ba61a1d5cba02a942d515556d87ac0d1e3fff854` (mtime `2026-10-04 16:49:16 +08:00`); the eight relevant source/input/config pins still match the frozen hashes. This version change is recorded in `manifest_version_recheck.json`; `post_manifest_input_recheck.json` confirms those eight current files still match their prefreeze hashes and that the final manifest pins still match. The arithmetic results remain tied to the unchanged frozen inputs, and do not claim a final launcher/runner freeze revision.

## Geometry, IDs, volume, and mass

Each CSV has exactly the schema `point_id,x_m,y_m,z_m,volume_m3`, 125 finite three-dimensional coordinates, and unique stable IDs `1..125` in row order. Both metadata files declare coordinates in metres and volumes in cubic metres; the generator records row order as the native index mapping and retains `point_id` as source identity (`make_point_inputs.ps1`, lines 25–46; metadata lines 5–9).

Every row has positive volume `1e-9 m^3`. The independently summed total is `1.2500000000000002e-7 m^3` in binary floating-point arithmetic, matching the declared `1.2500000000000002e-7 m^3` and the intended `125 × 1e-9 = 1.25e-7 m^3`. At the common declared density of `2100 kg/m^3`, total mass is `0.00026250000000000004 kg` (the expected `2.625e-4 kg`, subject to float representation); each point contributes `2.1e-6 kg`.

Both calculated volume-weighted COMs and unweighted centroids are `(0.0025, 0.0025, 0.0025) m`. The uniform coordinates span `0.0005` to `0.0045000000000000005 m` on each axis. The jittered coordinates span `0.00030000000000000003` to `0.0047 m` on each axis. All points remain inside their declared original cells and the 5 mm cube. Calculated geometry and totals agree with both metadata files and the manifest's 5 mm extent, 125-point count, total volume, and center declarations.

The axial load metadata uses axis 3 and the same stable IDs for both inputs: bottom `1..25` and top `101..125` (`make_point_inputs.ps1`, lines 61–64 and 104–107). These IDs map to the bottom and top source cell layers. Their z ranges are `0.0005 m` and `0.0045 m` for the uniform fixture, and `[0.0003, 0.0007] m` and `[0.0043, 0.0047] m` for the jitter fixture.

## Difference between the fixtures and jitter limits

The uniform fixture places one point at each 1 mm cell center. The second fixture applies the source's fixed seed-17 modular offset formula in each original cell, retaining the same cell volume (`make_point_inputs.ps1`, lines 31–46). All 125 jitter coordinates reproduce that formula; all 125 original cells contain exactly one point. 124 IDs move and one has a zero offset in all three components.

Each coordinate offset is one of `−0.2, −0.1, 0, +0.1, +0.2` cell widths, with 25 occurrences of each value on each axis. The largest computed absolute ratio was `0.2000000000000001` cell widths from floating-point subtraction; it is within `1e-12` of the formula's `0.2` bound. The mean offset is approximately `2.22e-18` cell widths per axis. The pattern is deterministic modular arithmetic, not independent random jitter. Balanced component histograms and near-zero mean do not establish that grid or orientation bias has been removed; this review did not compare rotated fixtures or independent jitter realizations.

The two point sets are geometrically distinct beyond rigid rotation/translation or point permutation: after sorting their 7,750 pairwise distances, 7,694 differ at `1e-12 m` resolution; the largest sorted-spectrum difference is `0.0004188989228913684 m`. That proves a change in geometry, not improved physical representativeness.

The metadata describes the weights as known voxel-cell quadrature for the uniform set and retained original voxel-cell quadrature after jitter (`make_point_inputs.ps1`, lines 66–75; metadata quadrature sections). Both explicitly set `voronoi_claim = false` and `real_shell_claim = false`. These weights are not exact arbitrary Voronoi volumes, do not infer a parent specimen's volume, and do not identify a measured shell.

## Shared conditions and comparison boundary

The two inputs can use the same shared declared setup: horizon `0.003 m`, OSB material values `rho=2100 kg/m^3`, `E=5e9 Pa`, `nu=0.24`, `Gc=45 J/m^2`, safety factor `0.8`, and `AUTO_NATIVE_STABLE_BOUND` time-step mode (`common_case.toml`, lines 1–16). The config labels the material values `UNVALIDATED_PROTOTYPE` and illustrative, without specimen calibration (lines 11–12).

The prepared case files specify different conditions: zero load for one step versus symmetric axial velocity `1e-4 m/s` on axis 3 for five steps (`case_zero_load_control.toml`, lines 4–10; `case_tinyload_event_diagnostic.toml`, lines 4–10). A controlled layout comparison must pass the same one of those case files to both pointsets. At the audit-time launcher revision, the launcher selected only the uniform CSV/metadata and zero-load case; it did not select the jitter input. Thus these files prepare comparable inputs and shared constants, but do not establish an executed same-case A/B comparison.

## Execution evidence and limits

The actual arithmetic audit used bundled Python `3.12.14`, with `audit_point_inputs.py`, and exited `0`; stdout is `audit_stdout.json`, stderr is empty, and the exact executable/argv are recorded in `audit_argv.json`. `audit_exit_code.txt` records the numeric exit code. A prior attempt through the WindowsApps `python.exe` alias exited `9009` before the audit script started; that separate attempt is preserved under `python_alias_attempt_*`.

No Julia, solver, native run, pytest, or launcher preflight was executed. The checks establish only static consistency of the frozen synthetic point data and declared metadata. A native runtime result, by itself, would not make these generated fixtures eligible as real artifact candidates.
