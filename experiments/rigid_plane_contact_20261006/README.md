# Rigid plane / point-mass contact diagnostic

This isolated experiment provides a small force and potential API for a point at `position_m` and one fixed infinite plane. The plane's `unit_normal` points into allowed space. The signed gap is `g = (x - x_plane) dot n`; `g >= 0` is free, while `g < 0` gives `F = -K*g*n` and `U = 0.5*K*g^2`. `K` is positive and expressed in N/m. The derivative relation `F = -grad_x(U)` follows directly from this piecewise potential; there is no damping or friction term.

All dimensional API inputs are bare finite SI scalars or 3-vectors, identified by parameter names (`_m`, `_N_per_m`). Unit-bearing objects and strings are rejected so callers must convert before calling. A bare number carries no detectable provenance, so the API cannot determine whether a caller supplied a numerically plausible value in a mislabeled unit. Plane normals must already have unit length within `1e-10`; they are not silently normalized.

`RigidPlanePointContact` exposes signed gap, force, and potential methods. Top-level `signed_gap_m`, `contact_force_N`, and `contact_potential_J` convenience functions are also provided. These routines do not resolve rigid bodies, fragments, finite surfaces, moving planes, rotation, or multiple contacts. This component is not a complete collision engine and is not wired into production code.

Run the one bounded API-check and toy trajectory using the workspace Python runtime:

```powershell
& 'C:\Users\xiemo\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe' .\run_toy_diagnostic.py
```

The driver uses Velocity-Verlet for a mass initially `0.2 m` above a plane with normal velocity `-2 m/s`, `m=1 kg`, and `K=100 N/m`. It records two step sizes (`0.002 s`, `0.001 s`, 720 total steps) and compares the rebound time and speed to the piecewise analytical solution. It saves the actual samples, force/potential/kinetic energy, contact impulse, error CSV/JSON, and a PNG. The impulse balance is `m*(v_final-v_initial)` against the trapezoidal integral of the sampled contact force. Energy residuals and their halving-step trend are reported as measured discrete errors, with no conservation pass/fail claim.
