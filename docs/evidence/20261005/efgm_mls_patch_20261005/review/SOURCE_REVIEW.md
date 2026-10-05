# Limited source review: MLS-to-elastic adapter

Status: **source-only, awaiting root review**. This review does not represent
an execution, numerical result, test pass, operator acceptance, or stability
claim. No dependency source was changed.

## Scope and ownership

Owned implementation is limited to:

- `pgfrag/experiments/efgm_mls_patch_20261005/fixed_case.json`
- `pgfrag/experiments/efgm_mls_patch_20261005/mls_elastic_adapter.py`
- `pgfrag/experiments/efgm_mls_patch_20261005/operator_diagnostic.py`
- its adjacent README

The matching review/freeze records are in this directory. The read-only
dependencies are the phase-2 MLS shape API/configuration, the bounded elastic
operator, and the explicit-crack formulation contract. The selected paper's
published CPM/EFGM adoption is documented by the formulation contract; this
work makes no method-novelty claim.

## Static interface findings

The adapter uses the actual frozen `shape_functions_3d` API. It submits one
bounded batch with source coordinates `(27,3)` and query coordinates `(64,3)`.
The expected result shapes are weights `(64,27)` and gradients `(64,27,3)`;
the adapter maps source row `i` to canonical dense global coefficient ID `i`
and returns `support_ids = [0,...,26]`. The output gradient convention is
`gradients[q,i,j] = d(Phi_i)/dX_j` in `1/m`.

The adapter omits both visibility arguments for its explicitly no-crack,
identity-visibility fixture. It checks the operator's sequential query
against the cached batch row, propagates poor-support status, and preserves
the originating MLS UNKNOWN code in a diagnostic for every other UNKNOWN.
UNKNOWN shape arrays are passed with a non-OK status; the frozen operator
blocks before it validates/copies them and returns no partial stiffness.
Coordinate and ID snapshots are bytes-backed immutable arrays. Future crack
visibility is not supported by this adapter and cannot safely assert
`mask_locally_constant=True` without separate query-local evidence.

The fixed integration result is expected to have 81 displacement DOFs and
64 volume points. Coordinates are metres, volume weights m^3, material E is
Pa, K is N/m, and `0.5 d.T K d` is J. The operator uses engineering shear
`[eps_xx, eps_yy, eps_zz, gamma_xy, gamma_yz, gamma_xz]`, with
`gamma_ij = du_i/dX_j + du_j/dX_i` and shear constitutive diagonal `mu`. The
fixture material `(E=1e6 Pa, nu=0.3)` is a fixed engineering input, not a
ceramic material claim. Support radius `2.5 m` is a frozen numerical MLS
length, not a material value.

## Diagnostic bounds and interpretation

The runner's fixed-case guard covers cube bounds, grid, quadrature count and
volume, material, MLS reference values, and public operator caps. Its
preflight estimate is checked before coordinate, MLS, or operator arrays are
allocated; the estimate covers explicit numeric arrays up to 32 MiB and
plainly excludes Python object and vendor-BLAS workspace. The outer process
uses a 30-second subprocess timeout that includes the worker's numerical
module imports, so a stuck native call is terminated and produces a JSON
timeout record. The worker makes exactly one MLS batch and one operator
assembly call; the runner has no support tuning, retry, solver, load, crack
geometry, or matrix repair path.

The planned record covers MLS statuses/support/conditioning, affine
interpolation errors if all rows are known, K symmetry, six rigid-vector
residuals, the complete symmetric-part spectrum, negative and additional
near-zero mode counts, and independent tensor-vs-operator affine strain
energy. Binary64-derived scales are labels only, not pass/fail gates. Six
rigid modes being present does not imply that only six modes are null.
Eigenvalues are calculated from `(K + K.T)/2` for reporting only; the
assembled K is never symmetrized, regularized, clipped, or made positive
definite.

Unresolved until the single later diagnostic: actual MLS status/conditioning
on this fixture, numerical affine errors, K metrics and any additional null
modes. If an MLS row is UNKNOWN, the expected behavior is a blocked assembly
and an explicit per-point reason, not fabricated zeros. Even a fully known
result only characterizes this fixed no-crack operator diagnostic. It does
not prove EFGM stability, connect an explicit crack surface to mechanics, or
validate fracture behavior.

## Execution status

No pytest, unit test, import-time execution, assembly, solver, benchmark,
render, or crack operation was run for this integration package. The next
authorized step must be one parent-controlled invocation after review, under
the source hashes listed in `source_freeze.json` and `source_freeze.sha256`.
