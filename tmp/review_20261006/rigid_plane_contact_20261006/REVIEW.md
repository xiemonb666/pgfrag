# Bounded rigid plane contact review

## Scope and implementation

`pgfrag/experiments/rigid_plane_contact_20261006/rigid_plane_contact.py` implements the assigned fixed infinite plane and point force/potential API. For signed distance `g=(x-x_plane)·n`, with `n` a unit normal pointing into free space, it uses `U(g)=0.5*K*min(g,0)^2` and `F=-grad_x(U)=-K*min(g,0)*n`. This follows by differentiating the stated piecewise potential; the force and potential are independent APIs. No production integration or complete collision engine is implied.

The toy's analytical reference follows from free flight until `g=0`, followed by `m*y_ddot + K*y=0` during penetration. With entry speed magnitude `v0`, the contact solution is `y(τ)=-(v0/ω)sin(ωτ)`, `ω=sqrt(K/m)`, so the first return to `y=0` occurs after `π/ω`, with ideal exit speed `+v0`. Its ideal contact impulse is `2*m*v0`. These are equations derived for this assigned model; no external calibration or damping/friction claim is used.

## One executed run

- Command: `C:\Users\xiemo\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe pgfrag\experiments\rigid_plane_contact_20261006\run_toy_diagnostic.py`
- Runtime: Python 3.12.14, NumPy 2.3.5, Pillow available; matplotlib was absent and was not needed.
- Exit code: 0; measured wall time: 0.229 s.
- Parameters: `m=1 kg`, `K=100 N/m`, initial gap `0.2 m`, initial normal velocity `-2 m/s`, end time `0.48 s`.
- Resolutions: `dt=0.002 s` (240 integration steps) and `dt=0.001 s` (480 steps), 720 total; 722 actual trajectory rows include each run's initial state.
- The focused API receipt covers free positive gap, linear penetration force/potential, rigid rotation behavior, invalid normals, and invalid SI scalar/nonfinite inputs. The actual files are `check_receipt.json`, `diagnostic_summary.json`, `actual_trajectory_force_energy_impulse.csv`, `comparison_errors.csv`, and `contact_trajectory_energy.png`.

## Measured results

| dt (s) | Exact rebound (s) | Detected rebound (s) | Time error (s) | Rebound speed error (m/s) | Analytic impulse error (N·s) | Max absolute energy error (J) |
|---:|---:|---:|---:|---:|---:|---:|
| 0.002 | 0.4141592654 | 0.4141540372 | 5.2281e-6 | 2.0214e-10 | 2.8435e-5 | 2.00003e-4 |
| 0.001 | 0.4141592654 | 0.4141579579 | 1.3075e-6 | 4.4229e-11 | 1.3301e-5 | 5.00012e-5 |

The maximum absolute total-energy residual over saved samples fell by a factor of `0.2500024` when dt was halved; the coarse and fine residuals are `2.00003e-4 J` and `5.00012e-5 J`. The final signed energy residuals were `+5.68710e-5 J` and `+2.66018e-5 J`. These are measured discrete Verlet errors, not a statement that energy conservation passed.

The contact-force trapezoid impulses were `4.0000284353 N·s` and `4.0000133009 N·s`; corresponding point momentum changes were `4.0000284353 N·s` and `4.0000133009 N·s`. Their balance residual was `8.88e-16 N·s` in each case. Relative to the ideal elastic reference `4.0 N·s`, the impulse errors were `2.84e-5 N·s` and `1.33e-5 N·s`.

## Execution receipt qualification

After the bounded run, visual/source review found the executed covariance assertion rotated about the plane normal, so it did not exercise a changed plane orientation. The driver source was strengthened to rotate the normal itself and later verified by a separate focused `run_self_checks()` call; see [FOCUSED_HISTORY.md](focused_readonly_check/FOCUSED_HISTORY.md) and its JSON receipt. No toy trajectory was rerun. I also corrected the two rebound rows in `comparison_errors.csv` after noticing their metric labels called for errors while the columns should contain the detected value and analytical reference; the corrected values were copied from the saved run summary and no simulation was rerun. The trajectory CSV and plot are direct outputs of the one executed run.

The API rejects unit-bearing objects and unit strings, but no runtime can infer the physical provenance of an ordinary bare float; callers must supply the SI units named by the API. This experiment does not implement damping, friction, moving planes, finite surfaces, multi-body/fragment contacts, self-contact, body rotation, or broad phase.
