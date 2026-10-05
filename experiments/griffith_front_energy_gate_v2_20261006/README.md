# Quasi-static energy gate v2

This isolated wrapper adds a strict quasi-static outer energy check to the one-step front controller. It reuses the v1 `propose_step` geometry candidate and explicit surface types. It deliberately does **not** call v1 `confirm_candidate`: the v1 `delta_potential_J` field has an unclear sign and meaning, so v2 declares `legacy_delta_potential_bridge = NOT_BRIDGED` and requires separate stored-elastic-energy and external-work inputs.

The only supported mode is `QUASISTATIC_ELASTIC_NO_INERTIA_NO_OTHER_DISSIPATION`. Dynamic and mixed-cohesive modes return `UNKNOWN` and leave the old surface active. In that supported mode, callers must declare the terms as `DELTA_STORED_ELASTIC_PLUS_EXTERNAL_WORK_SEPARATE`, provide:

```text
residual_J = external_work_increment_J
             - delta_stored_elastic_energy_J
             - Gc_delta_area_J
delta_total_potential_J = delta_stored_elastic_energy_J
                         - external_work_increment_J
allowed_abs_residual_J = abs_tol_J + relative_tol * reference_energy_J
```

`delta_total_potential_J` is derived for the audit record; it is never subtracted again in the residual. The gate checks finite nonnegative absolute and relative tolerances, a finite positive reference energy, units (`J`, `m^2`), area and `Gc × area` arithmetic, parent/candidate/state hashes, and that the absolute residual does not exceed the explicit limit. Zero tolerances are allowed when closure is exact. Missing parameters, bad units or bindings, and residuals outside tolerance return `NEEDS_RECOMPUTE`; the prior committed surface stays active and the candidate remains pending. Tolerance basis text and hash are recorded as provenance. They explain the caller's choice but do not validate its physical adequacy.

The evidence also records solver/source/energy-ledger hashes and caller provenance. These hashes and the `REANALYZED_CURRENT` string are declarations; this module does not verify an external solver, field source, or tolerance derivation. Caller provenance is metadata only. There is no `reviewed` flag that can substitute for numeric terms or units.

The external current-reanalysis path stays available and, after all bindings and the numeric quasi-static gate pass, activates an immutable surface version labelled `COMMITTED_EXPERIMENTAL_REANALYSIS`. It remains `EXPERIMENTAL_UNCALIBRATED` and the embedded surface contract remains `NOT_COUPLED`; the label is not physical validation. A separate `EXPLICIT_FIXTURE_ONLY` path can only produce `COMMITTED_FIXTURE_ONLY`. Targeted tests exercise the external API branch with synthetic hashes clearly marked `test_only`; they are not solver evidence.

## Required upstream data

From mechanical reanalysis of the exact proposed candidate, the adapter must provide the reanalyzed `ExplicitSurfaceState` with actual side displacement arrays; candidate ID and `REANALYZED_CURRENT`; solver-result/source hashes; parent-state, candidate-geometry and reanalyzed-state hashes; and an energy-ledger hash. The ledger must carry `delta_stored_elastic_energy_J`, signed `external_work_increment_J`, `Gc_delta_area_J`, area increment, and their declared units. It must separately provide finite `abs_tol_J`, `relative_tol`, positive `reference_energy_J`, a tolerance basis and its hash. The physics model must establish the supported no-inertia/no-other-dissipation assumption before selecting the supported mode. This gate checks inputs and arithmetic; it cannot establish that assumption itself.

`Gc_delta_area_J` counts the common new crack area once. The controller checks it against the proposed reference area and the `Gc` retained with that proposal. A total-potential change can be supplied only as `delta_total_potential_J` for explicit rejection; v2 never maps it to stored energy. The derived relation is `Δtotal_potential = Δstored_elastic - external_work` under the declared signs.

## Targeted fixtures

Run from the repository root:

```powershell
& 'third_party\vesselforge-runtime\venv\Scripts\python.exe' pgfrag\experiments\griffith_front_energy_gate_v2_20261006\run_fixture.py
```

The bounded suite covers balanced fixture and synthetic external-path commits, large positive and negative residuals, missing tolerances, total-potential mispass, dynamic/mixed-mode rejection, zero and invalid tolerances/reference energy, parent binding, and units. The runner writes an actual test receipt, state snapshots, and an SVG flow diagram under `pgfrag/tmp/review_20261006/griffith_front_energy_gate_v2_20261006/`. It does not call a solver/native runtime or write production outputs.
