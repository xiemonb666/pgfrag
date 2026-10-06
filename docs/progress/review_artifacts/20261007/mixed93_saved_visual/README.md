# Saved-array stress comparison (source prepared)

This renderer is deliberately gated. It only reads the failed producer's saved `arrays.npz` and a future independent recovery-audit JSON. It does not call the producer, sample new points, or interpolate a field. The NPZ is pinned to SHA-256 `3e1e6fd17751e4e06113010f20ea33c7b69914ca90431e26bbf76124951b5d41`; the producer returned rc 1 while writing `result.json`, so the NPZ is retained as recovered input, not a successful producer result.

`RENDER_PINS.json` records the independent audit at `pgfrag/tmp/review_20261006/mixed93_saved_independent/audit.json` by its exact SHA, byte count, schema, and status. The renderer checks the three exact `consistency_gates` keys, the failed NPZ path/size and pre/post hashes, producer rc 1 with no result JSON, `UNKNOWN_NOT_ASSESSED` fracture status, and the 93-DOF/99-KKT/27-source/4-branch layout before importing NumPy or Matplotlib. The audit says the recovered cache is independently consistent; its physics qualification remains unknown.

After authorization, it will create at most one PNG: four panels compare analytic versus solved full stress-tensor Frobenius norms at the same 128 saved off-grid points. For Voigt order `[sxx, syy, szz, sxy, syz, sxz]`, the norm is `sqrt(sxx^2+syy^2+szz^2+2*(sxy^2+syz^2+sxz^2))`; each panel reports the corresponding full-tensor relative L2 error and includes an identity line. Titles identify the manufactured-elasticity comparison, no crack growth, and producer rc=1. No surfaces or new coordinates are generated.

The source has not loaded NumPy, Matplotlib, or the NPZ and has not rendered a figure. Its normal no-argument path only prints a source-ready marker. The source and audit pins are for root review before any separate GO.
