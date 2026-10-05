# Native reversible elastic coupon — `dt/2` source-only comparison

**State: `SOURCE_ONLY_NOT_RUN`.** This is the sole time-step companion to the completed [baseline experiment](../free_evolution_coupon_20261005/). It reuses byte-identical driver and launcher sources. No Julia or solver was run for this companion, and its result files do not yet exist.

## Controlled comparison

The baseline case and the companion use the same analytic `uniform_box` geometry (`12 × 9 × 9 mm`, `h = 0.6 mm`, `δ = 1.8 mm`, `δ/h = 3`), 4,500-point grid, 1,134 full-horizon interior points, OSB `NoCorrection`, prototype `E = 5 GPa`, `ν = 0.24`, `ρ = 2,100 kg/m³`, `Gc = 45 J/m²`, failure-disabled intact response, native VelocityVerlet grip boundary conditions, continuous physical-time loading path, phase durations, held-grip dwell, single thread, and all point/step/output limits. All material values remain `UNCALIBRATED`; no surface, contact, fracture, or real-object qualification is implied.

The only numerical input change is `solver.dt_fraction_of_native_raw_bound`; descriptive `case_id` and `purpose` strings distinguish this companion. The baseline requested `0.25` and reached `Δt = 7.5e-8 s` because its shortest-phase resolution cap was also `7.5e-8 s`. The completed baseline recorded native raw stability bound `3.4242302121147243e-7 s`, 48 steps, and 49 states. This companion requests the Float64 ratio `0.10951366490292391`, calculated as `3.75e-8 / 3.4242302121147243e-7`; multiplying the recorded baseline bound by that ratio gives `3.75e-8 s`, below the unchanged phase cap. If the later run computes a different native bound or accepts a step other than exactly half the baseline `7.5e-8 s`, it is not to be reported as this matched `dt/2` comparison.

At `Δt = 3.75e-8 s`, the fixed physical duration `3.6e-6 s` resolves to 96 steps and 97 states. The 21-field Float64 tensor is expected to occupy 73,332,000 bytes (about 69.94 MiB). Those are design expectations; native bound, accepted step, step count, output size, and all response metrics remain unobserved for this case. The fixed case retains the same 5,000-point, 192-step, 190 MiB experiment-output, 200 MiB total-attempt, 180 s, and single-thread limits.

## Baseline evidence and comparison questions

The baseline's one authorized native run completed with exit code 0 at `pgfrag/experiments/free_evolution_coupon_20261005/run_attempt_01/`; its source freeze hash was `3d0a3cadf244eeb4612ec8999a9328b4be06cd087e16b08a731551250a4614b2`. Baseline integrity review found all 49 states finite, zero damage, and unchanged 409,556 directed active bonds. The baseline reported final OSB strain energy `2.4345788925163157e-4 J`, kinetic energy `2.5289129122172413e-4 J`, cumulative discrete BC kinetic work `4.606984761629457e-4 J`, actuator trapezoid surrogate `6.047855000282772e-4 J`, and diagnostic energy residual `3.565070431041004e-5 J`. The residual is not a qualified energy-closure or fracture-energy measure. Its stderr also contains the native message `stepsize specified! Please be sure that the CFD-condition holds!` from the VelocityVerlet timestep path.

After root review and a separate run authorization, compare the full displacement/velocity, reaction impulse, strain-energy, kinetic-energy, and diagnostic-residual histories at the same physical timestamps and loading path. Check that the observed native raw bound matches the baseline, the accepted `Δt` is exactly `3.75e-8 s`, the histories contain 97 states, the grip path returns to zero, all point data are finite, the failure-disabled damage and active-bond histories remain unchanged, and the free-DOF/momentum checks remain within the reported numerical tolerance. Record whether the energy residual or dynamic response changes; do not predict a scale factor or treat a smaller residual as energy qualification. This failure-disabled pair cannot test crack initiation, crack-path freedom, fracture dissipation, material calibration, or real-object physics.

## Source and reproduction boundaries

The driver and one-shot launcher are byte-identical copies of the completed baseline sources; the launcher resolves paths from its own experiment directory at the same nesting depth. The companion freeze carries forward all 78 shared environment/runtime/vendor pins and replaces only the four experiment-local pin paths and their hashes (fixed case, driver, launcher, and this review). No vendor source or baseline artifact is modified.

`run_once.ps1` must not be invoked until root has reviewed this freeze and separately authorized the one-shot run. It writes only `run_attempt_01/` and refuses an existing attempt. Any run that reports a changed stability bound, non-half step, time-step warning, nonzero exit, timeout, or cap event must be reported as observed; do not edit the case, driver, or launcher, and do not retry under this authorization.

```powershell
& 'C:\Users\xiemo\Desktop\ShardForge\pgfrag\experiments\free_evolution_coupon_dt_half_20261005\run_once.ps1'
```

This command was not run. The companion remains `SOURCE_ONLY_NOT_RUN` pending root's pre-run review.
