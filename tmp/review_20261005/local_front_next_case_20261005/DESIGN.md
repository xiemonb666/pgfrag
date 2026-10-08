# Local-front Mode III next-case design

Status: `SOURCE_ONLY / PROPOSED / NOT_COUPLED`. This is a formulation freeze for root review. It is not an implementation, solver run, convergence result, or authorization to promote the current global bridge.

## Decision

The next meaningful case should use continuous, unmasked MLS with **one local enrichment at most per source**: shifted F1 for supports containing the actual straight front, shifted Heaviside only for supports topologically split by the finite crack and separated from its front, and BASE for every source. The selected geometry below has 100 sources, 140 scalar UZ coefficients, and a front-core radius where the exact Mode III field is reproduced by a direct partition-of-unity identity.

This local construction does **not** promise exact Mode III reproduction over the whole domain. The collar where F1 support coverage blends into Heaviside support coverage must be measured and must converge. If root requires global exact inclusion as a condition of the next case, retain the existing all-site F1 bridge only as a separate control; that control cannot qualify local front enrichment.

## Frozen geometry and support roles

Use a rectangular domain with `s=0.004 m`:

```text
Ω = [-2.5s,2.5s] × [-2.5s,2.5s] × [-2s,2s]
Γ = {(x,0,z) ∈ Ω : x < 0}
F = {(0,0,z) : z ∈ [-2s,2s]}

x,z sources: {-2s,-s,0,s,2s}
y sources:   {-2s,-s,s,2s}
source count: 5 × 4 × 5 = 100
Wendland C2 support radius: h = 1.75s = 0.007 m
```

MLS remains linear, continuous, and unmasked, with origin `(0,0,0)` and reference scale `R0=2.5s=0.01 m`. These dimensions are a proposed bounded case, not evidence that every MLS row has full rank; source coverage and moment rank remain preflight gates.

For source `i`, define its closed physical support `S_i = Ω ∩ B_h(X_i)`. Determine roles from the **finite physical Γ and F**, never from an infinite-plane surrogate:

* `T_F`: `S_i` intersects `F`. Here that is `x∈{-s,0,s}`, `y∈{-s,s}`, for every z row: 30 sources.
* `J_H`: `S_i \ Γ` has exactly two nonempty connected components, one on each side of `y=0`; every point of `S_i` on `y=0` lies on Γ; and `dist(S_i,F)>0`. Also require the source itself is off Γ so shifted `H_i` is defined. Here that is `x=-2s`, `y∈{-s,s}`, for every z row: 10 sources. The support's maximum x is `-0.25s`, and its front clearance is at least `(√5−1.75)s ≈ 0.486s`.
* `BASE`: every source, including sources with no enrichment.
* Any support intersecting Γ that is neither in `T_F` nor certified for `J_H` is `UNKNOWN` and blocks the case. Tangencies and ambiguous topology do not default to H.

These source sets are disjoint. The frozen registry already has distinct `FRONT_F1_SHIFTED` and `HEAVISIDE` identities, so register 100 BASE keys, 30 shifted-F1 keys, and 10 H keys: **140 scalar DOFs**. No source receives both enrichments. Supports from different sources may overlap spatially; that is expected for MLS and must be checked for column dependence.

A cut-coverage bound is available from the proposed lattice. Along each crack-face trace, an enriched source center has x coordinate in `{-2s,-s,0}` and y coordinate `±s`; nearest x and z offsets are at most `s/2` and y offset is `s`. Thus an enriched support is within `√(1.5)s < h` of every point on Γ. The cut has no BASE-only gap.

## Basis and traces

For `r=√(x²+y²)`, use one global straight-front chart on `Ω\\Γ`:

```text
θ = atan2(y,x) ∈ (-π,π), with its cut only on the actual Γ
Fbar = √(r/R0) sin(θ/2)
H = sign(y) on each side of Γ
```

All source y coordinates are nonzero, so `Fbar(X_i)` and `H(X_i)` have unambiguous values. For `i∈T_F` and `j∈J_H`:

```text
BASE:       φ_i
F1:         φ_i (Fbar-Fbar(X_i))
Heaviside:  φ_j (H-H(X_j))
```

Use full product-rule gradients. Evaluate H and F1 traces with an explicit `PLUS` or `MINUS` side tag; never offset a query. On each open side, `∇H=0`, so the H row uses `(H-H_i)∇φ_i`; the crack jump is represented by separate traces in the broken-domain weak form, not by a delta term in the bulk gradient. F1 uses its analytic gradient, which is `O(r^-1/2)` and square-integrable in the transverse plane. At `r=0`, the value/gradient row is `UNKNOWN`; quadrature must avoid the exact front and resolve its singular integral by refinement.

The Heaviside support test means an H jump can occur only on Γ: every eligible support lies wholly behind x=0 at its y=0 cross-section, and Wendland weights vanish at support boundaries. F1 has its single branch cut on Γ. With continuous MLS and no visibility source dropping, neither term creates a shadow jump or an extension ahead of the front.

For the exact scalar Mode III reference, use `u_z=A√r sin(θ/2)`, `A=(KIII/μ)√(2/π)`, `u_x=u_y=0`, with the current bridge's isotropic values `E=1 MPa`, `ν=0.3`, and `KIII=100 Pa√m`; set `μ=E/[2(1+ν)]`. Apply the analytic Mode III traction on the rectangular outer boundary. The scalar traction on the z end caps is zero. On the physical crack faces:

```text
u_z(x,0+) = +A√(-x),  u_z(x,0-) = -A√(-x),  x<0
[u_z] = 2A√(-x), and t_z = μ ∂u_z/∂y = 0 on both faces.
```

The jump closes as `√r` at F. For x>0 the field is single-valued and continuous. End-cap scalar traction is zero; full-vector end-cap equilibrium remains outside this UZ-only case.

## Exact core versus blending collar

For this lattice, the nearest source not in `T_F` is `(0,±2s,z)` at distance `2s` from F. Therefore

```text
r_core = min_{i∉T_F} dist(X_i,F) - h = 0.25s = 0.001 m.
```

For any query with `dist(q,F)<r_core`, every source whose compact support contributes at q must belong to `T_F`. Hence, when MLS is OK, `Σ_{i∈T_F} φ_i(q)=1` there. Set the shifted-F1 coefficient to `c=A√R0` for each T_F source and its BASE coefficient to `c Fbar(X_i)`. Then

```text
Σ_T φ_i c Fbar(X_i) + Σ_T φ_i c(Fbar-Fbar(X_i)) = c Fbar.
```

This is a local exact-inclusion argument from the support geometry and MLS partition, not an inference from affine reproduction or a full-rank mask test. It establishes exact Mode III field reproduction only in the front core.

Outside `r_core`, BASE-only supports and fully split H supports can enter the row, so `Σ_T φ_i` need not equal one. The same coefficient vector is not guaranteed to reproduce global F1 there. This collar is the principal approximation risk. Report the exact field and energy errors separately in the front core, blending collar, crack region outside the front core, and uncracked x>0 region. Require independent quadrature refinement followed by support/source refinement; freeze numerical tolerances and runtime limits at root review. Do not hide a collar error with a ridge, pseudoinverse, or rank truncation.

## Alternatives and recommendation

**A — topology-selected shifted PU (recommended minimum case).** Use the exact sets and basis above. It reuses continuous MLS, the existing registry kinds, and the antiplane row contract. The local core is exactly included, H support cannot reach a ghost cut, and no source receives duplicate enrichment. The risk is a blending/conditioning error in the F1-to-H collar. Accept only if the collar, full-domain error, crack-face traces, spectrum, and quadrature all converge under separately recorded refinements. This is an approximation benchmark, not an exact-field pass.

**B — ramp-corrected front PU.** If Alternative A fails the collar gate, add a smooth front window and normalized front partition: `W_F=Σ_{i∈T_F}φ_i`, `φ̂_i=φ_i/W_F` where `W_F>0`, and `χ_F φ̂_i Fbar` with a C1 window `χ_F=1` in the front core and tapering to zero before the normalized denominator becomes small. Include all product and normalization derivatives. Keep H eligibility unchanged and blend its active region with a separately specified continuous window. Check window sums, derivative sums, traces, and augmented rank. This can remove the front-side partition loss in the chosen window, but adds geometry-dependent basis terms and risks large gradients/near-dependence at the taper; it is a separate formulation requiring review. It still does not imply exact global Mode III representation outside the window.

The existing global all-source F1 provider remains the simpler global-reproduction control. Do not combine global F1 and H as a production alternative: both create crack jumps over overlapping supports and invite redundant modes.

## Bounded work and acceptance gates

* Counts: 100 unique sources, 140 scalar coefficients, at most 256 quadrature queries per provider/assembly batch. Stream bounded batches. A 256×140×3 float64 gradient block is about 0.86 MB; the dense 140×140 stiffness matrix is about 0.16 MB.
* Keep volume integration broken at Γ. A first front-centered polar core batch can use GL4 in r, GL8 on open `θ∈(-π,π)`, and GL4 in z: 128 points with the cylindrical Jacobian r. Use positive radial nodes, then cut-aware exterior subcells; refine each batch independently without querying `r=0` or a fictitious offset. Crack-face trace checks use exact PLUS/MINUS side tags. Boundary and trace batches each remain under 256.
* Topology audit: record each selected support ID, source coordinate, front distance, front clearance, whether Γ intersects the support, connected-component evidence, and classification. Require exactly 30 `T_F` and 10 `J_H` for this geometry; any different count means geometry/implementation drift and blocks this frozen case.
* Basis audit: continuous MLS only; verify partition sum and gradient sum on OK rows; verify the local exact-core coefficient identity; confirm F1/H traces jump only on Γ, exact crack-face traction is zero, and no enrichment support reaches x>0 with H.
* Approximation/stability audit: report core, collar, far-crack, and x>0 errors; compare field, energy, equilibrium/work, and spectrum under quadrature then source/support refinement. Inspect enriched-column singular values or stiffness eigenvalues for additional near-zero modes. No regularization or unreported null mode. Affine patch and full-rank rows alone do not qualify this case.
* This restricted antiplane benchmark has one representable UZ translation and cannot establish six rigid modes or full 3-D equilibrium. It does not validate a curved front, penny crack, propagation, R2, or material calibration.

## Source interfaces to update after root accepts the formulation

1. Create a new fixed-case manifest for this geometry; leave the existing 27-site algebraic bridge unchanged.
2. Add an independent support-topology record/classifier that computes finite-Γ/F distances and the support-side connectivity evidence. It must return UNKNOWN on ambiguous/tangent cases; a bare eligibility boolean is not a certificate.
3. Generalize `ContinuousMLSFrontProvider`: replace `2*source_count` and “BASE plus unshifted F1 for every source” assumptions with BASE-all plus the frozen role manifest; support shifted F1, shifted H, product-rule gradients, and explicit PLUS/MINUS traces. Keep the dense row in canonical registry order and fail closed on any unsupported row.
4. Keep `DofRegistry` identity/caps, but make the provider validate exactly one enrichment key per eligible source and none for BASE-only sources. The existing source protocol already has `HEAVISIDE` and `FRONT_F1_SHIFTED` keys; it intentionally does not prove topology.
5. Parameterize the diagnostic runner that currently hard-codes 27 sources, 54 DOFs, all-site `T_F`, empty `J_H`, and 256 total queries. Add streaming ≤256-point cut/front-aware integration batches and separate trace checks; do not pass local coefficients into the legacy 3-DOF/node assembler.

Inputs reviewed: `astra_front_shadow_decision_20261005/DECISION.md`, `enriched_dof_contract_20261005/CONTRACT.md`, and `modeiii_pu_operator_20261005/{README.md,fixed_case.json,mls_modeiii_provider.py}`. No official source was needed for this formulation audit; no source code, solver, test, CAD, or numerical evaluation was run.
