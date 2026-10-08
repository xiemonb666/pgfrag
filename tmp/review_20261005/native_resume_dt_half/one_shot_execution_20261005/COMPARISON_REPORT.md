# Native coupon `dt/2` comparison — read-only post-run report

The authorized companion run completed once with launcher/native exit code `0` in `38,077 ms`. It used 4,500 points, 96 steps, 97 states, `Δt = 3.75e-8 s`, and the exact baseline native raw bound `3.4242302121147243e-7 s`. Thus the accepted timestep is exactly half the baseline `7.5e-8 s`; this pair is `MATCHED_DT_HALF_DIAGNOSTIC_ONLY`. No timeout or output-cap event occurred. The native stderr timestep warning is preserved.

Before and after the run, all 82 unique source pins matched the freeze; the freeze SHA remained `684439db34e2375d5996eed0dd1d518a413565cbaf92cd96a42f4d3799c70aef`. Baseline summary/CSV/binary hashes were also identical before and after. The full pre/post inventories and fieldwise comparison are in [pre_run_hash_audit.json](pre_run_hash_audit.json), [post_run_hash_audit.json](post_run_hash_audit.json), and [baseline_vs_dt_half_comparison.json](baseline_vs_dt_half_comparison.json).

The fixed residual normalization for both runs is `E_ref = max_t(OSB strain energy + full-step kinetic energy)` over the baseline CSV states: `5.099516633924144e-4 J`. Each normalized residual below is `abs(energy_balance_residual_diagnostic_J) / E_ref`; the denominator is shared by both runs and held constant.

| Measure | Baseline, `Δt = 7.5e-8 s` | Companion, `Δt = 3.75e-8 s` | Companion minus baseline |
|---|---:|---:|---:|
| Steps / states | 48 / 49 | 96 / 97 | — |
| Final absolute energy residual | `3.5650704310e-5 J` (`6.990997% E_ref`) | `1.7510918924e-5 J` (`3.433839% E_ref`) | `−1.8139785386e-5 J` |
| Maximum sampled absolute energy residual | `3.5650704310e-5 J` (`6.990997% E_ref`) | `1.7510918924e-5 J` (`3.433839% E_ref`) | `−1.8139785386e-5 J` |
| Final strain + kinetic energy | `4.9634918047e-4 J` | `4.9635925141e-4 J` | `+1.0070940926e-8 J` |
| Maximum sampled strain + kinetic energy | `5.0995166339e-4 J` at `2.4 µs` | `5.1034608936e-4 J` at `2.4 µs` | `+3.9442596714e-7 J` |
| Cumulative discrete BC kinetic work | `4.6069847616e-4 J` | `4.7884833249e-4 J` | `+1.8149856327e-5 J` |
| Actuator trapezoid work surrogate | `6.0478550003e-4 J` | `5.6097628338e-4 J` | `−4.3809216646e-5 J` |
| Final momentum x component | `−3.9263050007e-5 kg·m/s` | `−3.9640979433e-5 kg·m/s` | `−3.7792942623e-7 kg·m/s` |
| Maximum momentum norm | `3.3043911149e-4 kg·m/s` at `1.575 µs` | `3.3509100770e-4 kg·m/s` at `1.6125 µs` | `+4.6511896210e-6 kg·m/s` |
| Maximum point displacement norm | `1.3771297958e-6 m` at `2.625 µs` | `1.3758768577e-6 m` at `2.625 µs` | `−1.2529381397e-9 m` |
| Maximum full-velocity norm | `1.9145615323 m/s` | `1.918541452962 m/s` | `+0.0039799206 m/s` |

At the 49 shared output times, the maximum pointwise displacement-vector difference was `3.1471619742e-9 m`; RMS difference over all points and shared times was `3.9052510238e-10 m`. Both runs had zero damage and exactly 409,556 directed active bonds (204,778 canonical pairs) at every state. Independent binary checks found all Float64 fields finite and every per-point active-bond count equal to its reference count. Both CSV histories were numeric and finite; maximum free-DOF impulse residual was zero. The independently recomputed maximum momentum-update residuals were `2.6402016966e-18 N·s` and `2.0549019300e-18 N·s`.

These measurements do not establish an energy-closure or fracture result. Failure was disabled, bond counts did not change, and the work/residual fields remain diagnostics for staggered VelocityVerlet states. Material inputs remain uncalibrated; contact was not run and real-object physics qualification remains unknown. The measured lower residual at `dt/2` is recorded without assigning a convergence order or qualification meaning.

## Artifact pointers

- New run: `pgfrag/experiments/free_evolution_coupon_dt_half_20261005/run_attempt_01/`
- New binary: `point_state_history.f64le`, 73,332,000 B, SHA-256 `674a8b00ba433bdca840384a6e6cd9437fe78fbbdfeabc0293151322170e367c`
- New energy/reaction CSV: 68,940 B, SHA-256 `05eae634e4ae8eb89a670241499947dd4b667399c310803b9dc11db2327d8a97`
- New summary: 19,310 B, SHA-256 `2b48fc147ca4d45ee72ef5addfc8708e31188c57d0ca1bec1675c7cd9e916dbe`
- Flat launcher receipt: 2,497 B, SHA-256 `12148d9ced4774da3495ea1df85f9fa3ce33b290e8f19c3e4962059fe5105100`
- Native stdout: 98 B, SHA-256 `f814a58a5ea66147cab9fef2c5c9dad3e39ea452043def93ce121312731dccec`
- Native stderr: 227 B, SHA-256 `b87722a6adc7800d3eea59db8e3e633b3e2cb9e160eaf980b5772bce56980937`
