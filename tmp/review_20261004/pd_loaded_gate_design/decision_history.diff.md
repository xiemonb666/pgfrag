# Decision revision diff record

This is a human-readable correction record, not a source-code patch. Previous drafts remain byte-for-byte in `decision_v1.md`, `decision_v2_pre_wavecorrection.md`, and `decision_v3_pre_a_run.md`; current text is `decision.md`.

## v1 → corrected plan

- Replaced the corner-only `0.915–0.987 μs` timestep estimate with the accepted native run01 `dt = 4.925124976599115e-7 s` and native graph counts `6,124 directed / 3,062 canonical` from `zero_load_control_run01/summary.toml`.
- Recomputed A at 8 steps: `T = 3.940099981279292e-6 s`, `εnom = 9.850249953198229e-6 = 0.0061088657 εc`; no-break dynamic response is the sole A objective.
- Removed the proposed 320-step velocity threshold crossing. At 320 steps and ±0.02 m/s, the corrected estimate is `εnom = 0.0015760399925117167 = 0.97741851 εc`, so it does not nominally reach threshold.
- Replaced B with three independent affine deformation evaluations `F=diag(1,1,1+η)`, `η={0, 0.5εc, 1.5εc}`, using native OSB failure/force evaluation without time integration or manual bond edits.
- Recorded the native weighted-volume update `m_i(active bonds)` as an independent post-failure force/energy caveat.

## v2 → v3 correction

- Corrected the wave-time sentence: A spans about `3.9401 μs`; with `cL≈1675 m/s` and the 4 mm center-plane gauge, one transit is about `2.39 μs`, hence A is about `1.65` gauge transits. It remains an abrupt, explicitly dynamic 8-step window with no hold/equilibrium check.
- Added source proof that native VelocityVerlet initializes both `velocity` and `velocity_half` to zero arrays; the loaded runner branch does not rely on the zero-load case's initial-condition setup.
- No frozen physics gate, old runner, source, input, or historical run result was changed.

## v3 → v4 actual A result

- Following the main reviewer’s explicit one-run authorization, executed only A with the pinned runner/case/common/launcher/runtime hashes. Wrapper and child exited 0; no timeout or retry. B remains NOT_RUN.
- Read-only audit of the nine recorded states passed its narrow observer criteria: 125 points/state, native `dt=4.925124976599115e-7 s`, end time `3.940099981279292e-6 s`, mass `0.00026249999999999966 kg` each state, coordinate/displacement maximum residual `1.608322006572579e-18 m`, reciprocal directed/canonical activity consistent, and zero directed/canonical transitions.
- Nonzero U/V/internal body-force density confirms a loaded response. The result remains `UNVALIDATED_PROTOTYPE`; no fracture, material, contact, full energy ledger, rough-surface, or R2 qualification follows.
- The actual receipt, stdout/stderr byte counts and hashes, raw-output hashes, and per-step extrema are recorded in `pd_loaded_diagnostics/reports/loaded_elastic_a8_readonly_audit.json`; `manifest.json` was updated to distinguish A completed from B not run. Frozen source/config bytes were not changed.
