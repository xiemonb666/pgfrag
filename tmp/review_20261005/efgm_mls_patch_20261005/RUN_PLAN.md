# Deferred one-shot diagnostic plan

This is a plan only. No run is performed by this source-freeze task.

1. Root reviews `SOURCE_REVIEW.md`, the fixed JSON case, adapter, runner, and
   source/dependency hashes in `source_freeze.json`.
2. If root authorizes the single diagnostic, invoke from the repository root:

   ```powershell
   python pgfrag/experiments/efgm_mls_patch_20261005/operator_diagnostic.py `
     --output pgfrag/experiments/efgm_mls_patch_20261005/results/operator_diagnostic.json
   ```

3. Preserve the JSON and its embedded source hashes. Review MLS status first.
   A blocked/UNKNOWN result is a valid fail-closed diagnostic outcome. Do not
   retry it, change the fixed support/material/quadrature, fill UNKNOWN rows,
   add a ridge/pseudoinverse, or force K positive.
4. For an `OK` operator result, review the complete reported 81-value
   spectrum, all six rigid residuals, any additional near-zero modes, MLS
   interpolation errors, and independent tensor-energy comparison. Treat the
   binary64 scales as descriptive references only, never gates.

The outer runner enforces a 30-second wall timeout and the worker checks the
32 MiB explicit numeric-array estimate before MLS/operator allocation. This
step invokes one MLS batch and one assembly only. It runs no boundary-value
solve, load, crack visibility, SIF, fracture law, benchmark sweep, or render.
