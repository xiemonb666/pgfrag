# Scalar bilinear cohesive law (isolated review primitive)

This experiment implements a single opening-only Mode-I cohesive law. It does
not integrate with a solver or a crack front. API calls accept SI values only:
`k0_pa_per_m`, `tpeak_pa`, `gc_j_per_m2`, and `opening_m`. These unit-bearing
field names are part of the contract; strings, booleans, nonfinite values, and
unit-bearing external objects are rejected in place of silent conversion.

The envelope is `T = K0 δ` through `δ0 = tpeak/K0`, followed by linear softening
to zero at `δf = 2 Gc/tpeak`. Construction rejects parameters unless `δf > δ0`.
The state variable `kappa_m = max(opening_m)` controls irreversible damage. At
fixed history, traction and recoverable energy are `T = (1-d) K0 δ` and
`ψ = 0.5 (1-d) K0 δ²`. Cumulative dissipation follows the exact damage-work
integral and is clamped to `Gc` at complete failure. External work is integrated
along fixed-damage elastic unload/reload legs and the piecewise-linear envelope;
the ledger is `W = ψ + D`. Thus an elastic unload returns stored energy without
decreasing cumulative dissipation.

Compression and mixed-mode inputs raise `UnsupportedModeError`: those responses
are explicitly UNKNOWN here. Rate dependence, bulk material behavior, front
coupling, and a solver are outside the module.

The parameter fixture is illustrative software data, not a ceramic calibration.
Formula sources and the bounded run receipt are recorded in
`pgfrag/tmp/review_20261006/cohesive_scalar_law_20261006/run_receipt.json`.

## Reproduce

Run the focused tests and review-artifact writer with the existing VesselForge
Python environment (it supplies Matplotlib for the PNGs):

```powershell
$py = "third_party/vesselforge-runtime/venv/Scripts/python.exe"
& $py -m unittest discover -s pgfrag/experiments/cohesive_scalar_law_20261006 -p "test_*.py" -v
& $py pgfrag/experiments/cohesive_scalar_law_20261006/make_review_artifacts.py
```

The artifact path contains one CSV with 256 or fewer states, one JSON receipt,
and separate traction, history, and energy PNGs.
