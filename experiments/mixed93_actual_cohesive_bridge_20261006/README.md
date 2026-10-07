# Mixed93 fixed Mode-I cohesive bridge

State: `SOURCE_READY_EXECUTION_DISABLED`. This package freezes the next bounded design after independent recovery accepted the 93-DOF cache. It performs no numerical assembly, MLS query, cohesive trial, or solve.

Run the metadata-only preparation check with:

```powershell
python -B pgfrag/experiments/mixed93_actual_cohesive_bridge_20261006/mixed93_actual_cohesive_bridge.py
```

The driver AST-parses Python sources, checks source and metadata SHA-256 pins, verifies the accepted recovery audit and receipt, the explicit 27 source coordinates, and all 93 semantic tokens in order. It does not open or hash the saved NPZ, import NumPy or any numerical module, instantiate the MLS/provider/trace primitive, or call a solver. The `--execute` option always refuses. The case keeps execution authorization false pending separate root review.

The future one-attempt chain is fixed to the finite sheet `y=0, x in [-R,0], z in [-L/2,L/2]`, `R=L=0.01 m`, normal `+y`, 16 mapped 4x4 Gauss points, and 16 new continuous BASE MLS queries. It reuses the frozen 27-site grid/configuration and exact upper/lower `theta=+/-pi` branch limits at identical coordinates. The old 162-DOF provider is not adapted as the 93 registry.

The recovered `pure_I_q_solution` at `KI=100 Pa sqrt(m)`, exterior load, and physical gauge target are scaled by 0.01 to `KI=1 Pa sqrt(m)`; `q_exact` remains a separate diagnostic. The elastic cohesive probe fixes `K0=1e8 Pa/m`, `tpeak=1000 Pa`, `Gc=0.01 J/m^2`. These are toy numerical values, not material data. The future path permits one 99-order KKT solve and at most three primitive assemblies: virgin baseline, synthetic fully-damaged zero-tangent control (discarded), and virgin coupled-final verification.

The plan requires actual six rigid trace checks; `KG` has only zero-based entry `[82,82]` and equals `2*K0*R*L=20,000 N/m`. It checks the positive internal potential gradient separately from physical face-traction signs, the unchanged six-zero/87-positive spectrum, equilibrium, and coupled energy. Shear, compression, mixed-mode response, front motion, real-material calibration, and crack growth remain out of scope.

The accepted array digest is recorded as metadata for the later authorized run. The static check confirms that the recovery audit and receipt report the same digest without reading array contents.
