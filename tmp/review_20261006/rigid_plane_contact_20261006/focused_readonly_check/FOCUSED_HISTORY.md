# Focused read-only and rotation follow-up

This follow-up ran only `run_self_checks()` from `source_snapshot/run_focused_checks.py`; it did not call `simulate()` or the driver's `main()`. It completed with exit code 0 in 0.2021 s using the existing workspace Python 3.12.14 runtime. It produced `focused_check_receipt.json` and `focused_source_and_prior_artifact_pins.json`. The receipt and pins were moved from the runner's default nested output directory to this directory without changing their bytes.

The immutable plane point and normal are now arrays backed by immutable `bytes` storage. In the focused run, `.setflags(write=True)` raised `ValueError` for both fields and both arrays remained non-writeable. Six focused check groups passed, including the corrected general-orientation rotation covariance check.

For the z-axis rotation by `0.37 rad`, the unit normal changed from `[0, 1, 0]` to `[-0.361615431964962, 0.9323273456060345, 0]`. The signed gap remained `-0.02 m` (absolute difference `0`); the original force `[0, 2, 0] N` rotated to `[-0.723230863929924, 1.864654691212069, 0] N`, exactly matching the force recomputed from the rotated plane and point (reported covariance norm error `0 N`). The potential remained `0.02 J` (absolute difference `0`). Full values and source SHA-256 pins are in the JSON receipt and pin file.

## First toy-run history

At the original toy run, no immutable source snapshot or source SHA-256 was saved. Its old source hash is therefore unavailable; the present snapshot must not be treated as that historical source. The original `actual_trajectory_force_energy_impulse.csv`, `contact_trajectory_energy.png`, `diagnostic_summary.json`, and earlier `check_receipt.json` remain present and were SHA-256 pinned at this follow-up. They were not regenerated here.

`comparison_errors.csv` was edited after the toy run to make its rebound `measured` and `reference` columns match their names. The pre-edit CSV bytes and hash were not preserved. The pin file records the current post-edit hash and explicitly labels that history; it does not claim to recover or reproduce the earlier CSV. No source or data file is represented as an original-run snapshot unless it was actually saved at that time.
