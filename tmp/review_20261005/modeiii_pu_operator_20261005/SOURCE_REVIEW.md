# Mode III PU operator source review

Status: source prepared for root review; no tests, imports, worker, solver, or
numerical assembly have been run for this package.

## Owned source

- `pgfrag/experiments/modeiii_pu_operator_20261005/fixed_case.json`
- `pgfrag/experiments/modeiii_pu_operator_20261005/mls_modeiii_provider.py`
- `pgfrag/experiments/modeiii_pu_operator_20261005/modeiii_operator_diagnostic.py`
- `pgfrag/experiments/modeiii_pu_operator_20261005/README.md`

The experiment imports the pinned MLS phase-2, isotropic operator,
`DofRegistry`, and Mode III field-reference modules. Those dependencies remain
read-only and are recorded in `source_freeze.json`.

## Assembly contract

The provider implements the shared `AntiplaneShapeProvider.evaluate(Vector3)`
shape: input is an exact finite three-value tuple; each successful sample
returns canonical registry IDs and one `(value, gradient_xyz)` pair per UZ
coefficient. It precomputes a bounded, immutable batch, follows fixed query
order, and maps any MLS or field-reference UNKNOWN row to an empty UNKNOWN
sample with the original per-row status preserved in JSON diagnostics.

For each source `i`, the rows are `phi_i` and
`psi_i = phi_i * F1/sqrt(R)`. The enriched gradient uses the full product
rule. The factor is dimensionless, so every coefficient has displacement
units (m). Visibility and crack masks are absent. Source values of the front
function are not read.

The scalar strain column is
`[0, 0, dUz/dz, 0, dUz/dy, dUz/dx]` in the frozen engineering-shear ordering.
The 54-by-54 matrix is assembled directly from the frozen isotropic `D`; its
integrand includes `mu` for the x/y gradients and `lambda+2mu` for the z
gradient. No old three-DOF-per-node assembler is used. Unknown rows, malformed
shapes, nonfinite local matrices, or nonfinite accumulated matrices stop
assembly without returning a partial K.

## Fixed review inputs

- Cylinder radius and length are both `0.01 m`, centered at z=0.
- There are 27 sources on the fixed 3-by-3-by-3 grid and 54 scalar DOFs.
- MLS origin is `(0,0,0) m`, reference scale is `R=0.01 m`, and numerical
  support radius is `h=0.025 m`; `h` is not a material parameter.
- All 27 source support balls contain the straight front, so the fixed
  algebraic registry uses `T_F=all 27` and `J_H=empty`. This choice says
  nothing about local enrichment transitions or fracture qualification.
- Volume, lateral, and cap quadrature use 128, 64, and 32 points per cap,
  respectively, for 256 unique MLS queries. Radial nodes are positive, but
  the radial integration interval remains the full `[0,R]`; there is no hole
  or clipping. Open Gauss nodes avoid the branch cut and front.
- Material is `E=1e6 Pa`, `nu=0.3`, with `KIII=100 Pa sqrt(m)`. The fixed
  manufactured coefficient vector has zero BASE coefficients and
  `A sqrt(R) m` for each FRONT coefficient.

## Diagnostic boundary

The prepared runner would assemble one K and report actual MLS support and
conditioning, manufactured-field and gradient errors, tensor-integrated
energy versus K energy, traction work, symmetry, the single scalar UZ
translation residual, the complete eigenvalue spectrum, negative values, and
additional near-zero values. Its stated numerical scales are descriptive;
they are not acceptance thresholds. It does not solve, impose displacement
constraints, load crack faces, calculate reactions or a stress-intensity
factor, or claim full 3-D equilibrium. The UZ-only coefficient space cannot
represent the six full-vector rigid modes. There is no qualification of
EFGM stability, continuous-front fracture behavior, or material calibration.

The numerical outcome, including support rank, extra null modes, and any
negative eigenvalues, remains unknown until a separately authorized bounded
diagnostic run. No regularization, rank truncation, pseudoinverse, or
positive-definiteness forcing is present in the planned path.

## Review checklist

- [ ] Confirm the provider adapter follows the shared Vector3 tuple protocol.
- [ ] Confirm the scalar engineering-shear map and longitudinal term.
- [ ] Confirm all counts and array estimates are checked before numerical
  fixture allocation and the worker has a hard 30-second timeout.
- [ ] Confirm fixed dependency pins before any run.
- [ ] Only after review, perform at most one fixed-case diagnostic and retain
  its JSON, stdout/stderr, return code, timing, and pre/post source pins.
