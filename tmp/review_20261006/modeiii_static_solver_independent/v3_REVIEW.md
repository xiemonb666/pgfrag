# Independent static review — staging-integrity v3

**Finding: no blocker for the requested criteria.** V3 rejects a worker failure, nonzero exit, absent/invalid/oversized staging JSON, or failed success-evidence contract before it can promote the staging file. The numerical mechanics are unchanged from v2 by direct function-body comparison. No solver, tests, render, numerical assembly, integration, or Git command was run; this is not execution authorization.

## Freeze verification

The SHA-256 sidecar matches `source_freeze.json`. All v3 owned-source pins match the on-disk files: solver `7a8ac7076a788a29d7f1bcae34bcede80b1481030910c37da1af673a2a370a9f`, case `13ca46e05b520624d903be5d71c55e58b4d23db621423b47da4b2a0c174a38a9`, README `f257ce8e173f1867badf358a4dcd696734b3e674bfa7868b6a981a0a7e923ae4`, and source review `5e2bff119031201056092fd565521d80cef22fda0495d4c7e77f85ca8a68340c`. Manifest SHA-256 is `1e57fe4e81fca80789550aec376934f3ac450b819ce7963744d30009c901c5a0`. All seven read-only dependency hashes match. The preserved v1 and v2 solver and case snapshots remain at their recorded hashes.

## Promotion path

The parent checks the worker return code first. A timeout or any nonzero exit writes an UNKNOWN receipt and preserves the staging file without promotion (`static_solve.py:1700-1732`). A zero exit with no staging file also writes UNKNOWN (`1734-1751`). Otherwise, the parent reads at most the 4 MiB cap, decodes strict UTF-8, rejects JSON constants such as NaN/Infinity, parses JSON, and validates the v3 schema, case ID, solved status, source hashes, finite required scalars, 54 canonical IDs, force vector, both full matrices and mean rows, both solution vectors, and the 256-point off-grid data (`1753-1802`). Any parse or contract failure writes UNKNOWN while preserving the staging file.

Before promotion, the parent checks that the destination is still absent, rehashes the staging bytes against the bytes already validated, then uses `os.replace`; promotion errors produce UNKNOWN and retain the staging file when present (`1804-1825`). Thus a worker that writes a partial JSON and exits nonzero cannot be promoted; a zero-exit partial or malformed file is rejected by JSON parsing or the required evidence contract. V2's existing finite-value serialization and 4 MiB result cap remain in force.

## Mechanics comparison

The fixed quadrature builders, boundary and volume assembly, manufactured coefficient map, scaled spectrum, direct saddle solver, off-grid calculation, finite-value output records, and worker result generation were compared with v2. The mechanics functions are byte-for-byte identical. `_run_static_solve`, `_write_result`, and `_worker` are also byte-for-byte identical after normalizing only the v2-to-v3 result schema strings. The v3 case changes the result schema/case identity and documents the promotion contract; frozen geometry, material, source count, 54 coefficients, gauge, 128/256 theta-only levels, side/cap quadrature, off-grid points, budget, and numerical thresholds remain unchanged.

No `pinv`, ridge, coefficient pin, or alternate solver was introduced; the source still calls `numpy.linalg.solve` for the 55×55 scaled saddle. The result remains an unqualified manufactured antiplane diagnostic. Root retains the execution decision.
