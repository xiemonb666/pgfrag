# Independent static-solver review — frozen v1

**Review status:** Frozen-v1 static review complete. Do not treat v1 as ready for an authorized run: the output needs the author's distinct v2 review snapshots, and the finite-diagnostic issue below needs to be resolved in that output revision. No solver, tests, numerical assembly, integration, render, or Git command was run.

## Freeze verification

The freeze sidecar SHA-256 for `source_freeze.json` is `24d1d07ec2d61036e77686b8cbecfb72260f30fdc31658385e09ccb4847445da` and matches the file. All four owned-source hashes in the manifest match, including `static_solve.py` SHA-256 `d9294d9a90d20505d6ecd01354a2d9e0897f56b5ff0bb98aad01a4fd31513787` and `static_solve_case.json` SHA-256 `93c51463d717ad2aa8e78bc163fb8915aaf26164f411c2b46a20d80492bda42d`. All seven declared read-only dependency hashes also match. The interim case-pin mismatch is fixed in this frozen source: `EXPECTED_SHA256` now matches the actual case hash (`static_solve.py:46-54`).

## Static mechanics review

I found no high-impact algebraic or assembly defect in the frozen v1 source:

- `_build_sources` and `_build_registry` produce the fixed 27 semantic sources and 54 UZ coefficients in the registry's canonical order (`static_solve.py:301-317`). `_manufactured_coefficients` sets BASE to zero and every unshifted-F1 coefficient to `A*sqrt(R)` (`706-718`); with the provider's `F1/sqrt(R)` row and partition of unity this has the intended manufactured-field scaling.
- The engineering-strain map is applied per provider contribution before `BᵀDB`, while the mean row uses `1/(pi R² L)` times the physical volume weights (`604-703`). The units are consistent: K is N/m, f is N, c is dimensionless, and `cᵀd` is a physical mean displacement in meters.
- The side force is assembled directly from analytic `sigma*n` values and provider values, independently of `K*d_star` (`469-594`). Cap `t_z` is checked before exclusion; full cap tractions and resultants are retained as diagnostics. This is the intended restricted UZ weak load, not a full-vector boundary solve.
- The translation vector uses one in each BASE slot, the gauge check verifies its physical mean is one, and the scaled saddle uses `S=diag(K_aa^-1/2)`, the symmetric part of `S K S`, and a normalized `S c` row. The coefficient and multiplier unscaling is consistent with that change of variables (`721-900`). The implementation calls `numpy.linalg.solve` directly and has no coefficient pin, pseudoinverse, or ridge.
- The spectrum reports all scaled eigenvalues and the positive quotient condition, and blocks a missing/additional near-zero mode before solving. The 128-to-256 comparison refines theta only, and the off-grid metrics use a separately built 256-point set. No threshold turns either fixed-level sensitivity or the manufactured energy into a convergence claim. The result classification and limitations stay diagnostic-only (`1146-1327`).

The load compatibility diagnostic is partial: the result reports the integrated side z-force and the saddle multiplier, but has no explicit compatibility residual or check tying the force vector's translation projection to that resultant/multiplier (`590-600`, `853-900`, `1046-1062`). Because the run is explicitly unqualified, I do not treat this alone as a false mechanics claim; the requested v2 raw force and mean-row snapshots will make the relation inspectable without reassembly.

## Blocking output issue before a run

The recursive JSON normalizer correctly converts nested NumPy arrays and `np.generic` values—including `np.bool_` through `.item()`—to JSON-native values (`static_solve.py:968-979`), and `json.dumps(..., allow_nan=False)` rejects any nonfinite value that survives normalization (`1331-1334`). However, `_json_safe` changes every nonfinite Python float to JSON `null` (`977-978`). The solution is marked `OK` after checking only the direct-solve vector for finiteness; subsequent weak residual and backward residual calculations (`853-900`), solution energy/work and their relative errors (`1046-1062`), and off-grid error reductions (`903-965`) are not explicitly checked for finiteness. A nonfinite diagnostic can therefore become `null` while the top-level status remains `SOLVED_DIAGNOSTIC_ONLY_UNQUALIFIED`. Add an explicit finite-diagnostic guard in the output revision so that this condition writes UNKNOWN instead of silently serializing a successful-looking receipt with null metrics. This is a source-level possibility, not an observed run failure.

V1 also records only a SHA-256 of raw K and omits the raw K, f, c, canonical coefficient IDs, and the off-grid pointwise coordinates/reference/solution fields (`982-1035`, `1276-1327`). Thus the stored v1 result cannot support independent matrix/load reconstruction or data visualization without reassembling. I agree with the root's requested v2 output-only snapshots before any run. The worker writes one JSON receipt to staging before returning; there is no optional render path (`1331-1397`).

## Evidence boundary

This is a static source review only. It does not establish numerical solve success, physical equilibrium, local enrichment convergence, material calibration, fracture behavior, SIF, or crack propagation. I did not modify solver inputs or implementation, and did not run any code.
