# Source review record — local-front scalar operator

## Frozen question

Can the fixed local-front provider deliver all 256 valid MLS rows needed to assemble the bounded 140-DOF scalar antiplane Poisson stiffness matrix, and what do its full spectrum, affine patch, mean row, and symmetry diagnostics report? This package does not assume those outcomes. No provider, MLS, stiffness, solver, or pytest run has been performed in this source-preparation stage.

## Input and semantic scope

- Domain: `[-.01,.01] × [-.01,.01] × [-.008,.008] m`; finite crack `y=0, -.01<=x<0, -.008<=z<=.008 m`.
- Read-only provider: `pgfrag/experiments/local_front_provider_20261006/local_front_provider.py`, using its 100-source / 140-ID registry and `sample_at(query_index)` API. Every successful row must preserve `registry.canonical_ids` order.
- Geometric source roles: 60 BASE, 30 FRONT, 10 HEAVISIDE-eligible. Semantic DOFs: 100 BASE, 30 shifted-F1, 10 Heaviside.
- No visibility mask is passed. Volume points use no trace side because the split rule avoids the crack plane and front.
- Material fixture: `E=1e6 Pa`, `nu=.3`, so the scalar coefficient is `mu=E/(2(1+nu))`. This choice fixes the coefficient; it is not a calibration or a full elastic constitutive law.
- Bilinear form: `K_ab = sum_q mu * (grad N_a · grad N_b) * w_q`; shape gradients are 1/m, quadrature weights m³, K is N/m, coefficients m, and energy J.

## Integration and geometric coverage

The four x/y sign rectangles each have side lengths `0.01 m × 0.01 m`, and share the z interval of length `0.016 m`; each has volume `1.6e-6 m³`, so the full box is `6.4e-6 m³`. GL4 is mapped separately to every signed x/y interval and to the full z interval. Its nodes are strictly interior and its weights positive, so no query touches `y=0`, `x=0`, a box boundary, or a front endpoint. The code checks the 256-point count, positivity, non-duplication, open bounds, and volume sum before MLS evaluation.

The full source set is the Cartesian product of the frozen axes x,z=`[-.008,-.004,0,.004,.008] m` and y=`[-.008,-.004,.004,.008] m`. In each coordinate, the largest nearest-axis-site distance over the box is attained at a box endpoint or the midpoint of an adjacent source-site interval. Those maxima are `.002`, `.004`, and `.002 m`, so the tensor lattice nearest-site bound is `sqrt(.002²+.004²+.002²) = .004898979485566356 m`. The fixed `h=.007 m` leaves `.002101020514433644 m` of strict geometric support margin. The runtime source check verifies that the pinned records are exactly this complete lattice.

This geometric bound does not imply any particular MLS support count, polynomial rank, conditioning, or query status. Only runtime provider diagnostics can establish those values. A single split GL4 rule does not establish integration convergence: MLS products are not polynomials, and the shifted-front gradient is singular at the front. There is no radial hole, quadrature clipping, or epsilon displacement.

## Diagnostics and fail-closed behavior

- `K` is assembled only after all provider and MLS query rows report OK and carry finite values/gradients. A poor-support or malformed row blocks the entire operator; no partial matrix is returned.
- The BASE-only coefficient vector for `Uz=1 m` is the one expected scalar translation null mode. It does not assert six vector rigid modes.
- The BASE-only affine probe is `Uz=x`, with each BASE coefficient equal to its source x-coordinate and enriched coefficients zero. Its energy is compared to an independently evaluated scalar tensor rule and the closed-box expression `0.5*mu*6.4e-6`. It is diagnostic only and is not the local Mode III analytic coefficient vector.
- The full symmetric-part spectrum is retained. Significant negative eigenvalues or a near-zero count other than one yield `mechanical_eligibility=BLOCKED_INVALID_SPECTRUM`; `assembly_status=OK` still says only that K was formed. The worker returns rc=3, so the complete K/spectrum JSON stays in staging with a SHA-bound UNKNOWN receipt and cannot be promoted as the primary. The only promotable status requires exactly one near-zero translation mode, zero significant negative eigenvalues, normal eligibility, and `solve_allowed=true`; this means only that a later reviewed step may proceed. No solver is invoked here.
- The volume mean row and all probe vectors use the canonical registry order. No physical nodal Dirichlet clamping is inferred from MLS coefficients.
- The wrapper will only promote a complete, finite, schema/case/hash/dimension-validated JSON result after worker return code zero and within the hard deadline. Otherwise it preserves worker staging and writes a separate UNKNOWN receipt; existing output paths are never overwritten.

## Resource ceiling and remaining unknowns

The preflight bounds source count by 125, semantic scalar DOFs by 375, and volume queries by 256 before provider/query/operator allocations. For the frozen 100/140/256 fixture, the explicit-array estimate is 4,517,536 bytes, below 8,388,608 bytes. It includes 16 dense K-sized work arrays, and records that Python object overhead and vendor eigensolver workspace are excluded. The worker is externally timed at 60 seconds, with no retry.

Still unknown pending a separately authorized run: actual provider/MLS row statuses and original reasons, actual support counts/conditions, K entries/hash, eigenvalue signs and extra null modes, affine value/gradient error, affine energy mismatch, mean row, and symmetry residual. No solver, applied load, reaction, crack growth, or SIF is in scope.

## Pre-run review checklist

- [ ] Verify all local source and fixed-case hashes from `source_freeze.json`.
- [ ] Verify all provider, support, registry, MLS, and material pins match their read-only inputs.
- [ ] Verify a new empty run directory and output path; do not reuse historical artifacts.
- [ ] If authorized, run exactly once with Python `-B`; keep the JSON as the primary artifact, and retain argv, real return code, stdout/stderr, and elapsed time.
- [ ] Do not retry, change `h`, tune quadrature, regularize, clamp, or solve if support, spectrum, or patch diagnostics are poor.
