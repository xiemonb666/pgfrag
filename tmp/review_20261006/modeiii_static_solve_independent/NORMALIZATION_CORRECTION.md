# Off-grid error normalization correction

The captured `static_solve_independent_audit.json` and `independent_run_receipt.json` remain unchanged. The audit JSON's `uz_reference_scaled_max_error` and `gradient_error_max_reference_scaled` fields divided by the largest sampled exact value. The frozen case instead defines reference scales as `A*sqrt(R)` for `uz` and `A/sqrt(R)` for `grad(uz)`.

Using the pinned fixed-case values, `A=(KIII/mu)*sqrt(2/pi)=0.00020744998580874498 sqrt(m)`, `A*sqrt(R)=2.07449985808745e-5 m`, and `A/sqrt(R)=0.0020744998580874497` (dimensionless). Re-scaling the saved absolute errors in the audit JSON gives:

| Volume level | `max |grad error| / (A/sqrt(R))` | `RMS |grad error| / (A/sqrt(R))` | `max |uz error| / (A*sqrt(R))` | `RMS |uz error| / (A*sqrt(R))` |
|---|---:|---:|---:|---:|
| `volume_gl4x8x4` (128 points) | 17.3612% | 6.77525% | 3.64003% | 1.53263% |
| `volume_gl4x16x4` (256 points) | 0.851322% | 0.352767% | 0.202954% | 0.0682155% |

These corrected ratios are post-processing of the saved audit result only; the primary diagnostic was not reread, and no solver or field values were recomputed.
