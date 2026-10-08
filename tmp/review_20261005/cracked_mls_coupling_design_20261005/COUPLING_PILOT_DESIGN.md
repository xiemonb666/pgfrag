# Cracked MLS coupling pilot design — 2026-10-05

**Status: design only.** No source, test, solver, geometry, or Git state was changed or run for this review. This note defines a falsifiable assembly-level pilot and a front-shadow continuity gate. It does not define a calibrated crack law, physical fragment count, growth update, or measured/physical GT.

## Root decision and claim boundary

The proposed first case reuses the frozen 27-node, 64-volume-quadrature uncracked MLS/EFGM fixture and its exact MLS configuration. The facet is a finite internal triangle, with a three-dimensional front inside the body. The test changes only visibility and the split-DOF map, so it can reveal whether a surface enters the discrete operator.

That is not enough to call the result a correct traction-free crack. Sharp visibility creates source-to-front shadow boundaries away from the physical facet. A pointwise mask-stability certificate only justifies the derivative conditional on one fixed mask near that point; it does not show that the MLS trial field is continuous across the shadow boundary or belongs to the `H¹` space assumed by the ordinary volume-only Galerkin form. The front-shadow trace test below is therefore a **precondition** to interpreting `KΓ` as a traction-free crack operator. If the two one-sided shape rows differ outside the physical facet, the current Boolean-visibility/volume-only coupling is blocked pending a root decision on the mathematical method.

A passing pilot may support only this narrow statement: “the declared finite surface changes a bounded, discrete MLS elastic operator and has a separately addressable displacement jump under the tested masks.” It cannot establish fracture physics, growth, stress-intensity accuracy, contact, ceramic material calibration, vessel breakup, R2, or model R2. A K difference by itself is not evidence of a physically correct release of traction.

## Fixed candidate A — same 27 nodes and 64 quadrature points

Use the existing cube `Ω = [-1,1]³`, node coordinates `{-1,0,1}³` (27 coefficient sites), and the existing 64 Gauss points, each with volume weight `1/8`. The uncracked operator has 81 coefficient DOFs. Keep the baseline’s `MLS3DConfig` byte-for-byte: reference origin, reference scale, support radius, conditioning limits, and all tolerances come from the separately prepared uncracked certificate; do not change them for this pilot. Use the existing dimensionless material `E=1, ν=0.3`.

The visibility primitive documents coordinates and its absolute tolerance in meters. For this normalized fixture, record the explicit convention `1 reference-coordinate unit = 1 m` for geometric predicates and use one consistent conversion in visibility, MLS, and assembly. This is a numerical unit convention for a dimensionless pilot, not a physical specimen size or calibration.

The one crack surface is the finite triangle

```text
A = (0, -3/4, -1/2)
B = (0,  3/4, -1/2)
C = (0,  0,   3/4)
Γ = conv{A,B,C},   oriented normal n = +eₓ.
```

Its three edges are the explicit front. The complete triangle is inside the cube and does not touch the exterior. Among the 27 lattice sites, `(0,0,0)` lies strictly inside Γ; split that particle into separate `+` and `−` coefficient IDs. Do not split front sites. Use the accepted-article visibility-only offset `0.01 r₀` along `±n` when classifying the two copies; MLS weights and gradients still use their coincident original coordinates. Preserve the source-contract note that the thesis gives `0.001 r₀`; do not tune this offset. This yields 84 DOFs if the single added copy is supported by the current map.

The fixed 64-point rule has no point on `x=0`: its x coordinates are `±(1/2 ± 1/(2√3))`. It therefore samples both sides without putting a volume quadrature point on Γ or its front. Some projected points lie strictly inside the triangle and some outside it, including points whose source segments can pass beyond the finite front. Keep the quadrature points and weights unchanged. This exact coarse rule is an operator smoke fixture, not a convergence-qualified integration of a cut domain.

### Query support, split DOFs, and fail-closed mapping

For each quadrature query, classify each original source-to-query segment against the one finite triangle. Only `VISIBLE` and `BLOCKED` may be used. Any `UNKNOWN` visibility pair blocks that row. The separately prepared exact stability evidence is usable only if it explicitly covers this candidate’s finite-facet masks and the actual quadrature/source pairs; an uncracked support-boundary certificate alone is insufficient. Do not mark `mask_locally_constant=True` solely because the point is not visually on Γ. A mask switch, support-radius boundary, rank loss, condition failure, or duplicate active source coordinate returns `UNKNOWN` and stops the whole assembly with no partial K.

The current MLS API takes one source mask for a call, while the finite visibility result is query-major. The current EFGM provider is evaluated one material point at a time, so a provider can pass each query’s mask in a one-row MLS call. A multi-query adapter would have to group identical masks or evaluate each row separately; it must not broadcast one query’s mask across the 64 quadrature points.

The split-DOF adapter must present at most one active coefficient ID for each physical source coordinate in any one MLS fit. `shape_functions_3d` rejects duplicate active source positions. On the crack faces, plus/minus queries must select the corresponding split copy and give the two traces independent coefficients. At finite-front queries where both offset copies are `VISIBLE`, the current primitives provide no rule for collapsing or constraining the coincident copies while preserving front continuity. Do not pass both copies, silently choose one, average their displacements, or extend a half-space label past the finite front. Mark those rows `UNKNOWN_DUPLICATE_SOURCE_MAPPING` until the parent/root selects and justifies a mapping.

This coincident-copy issue is separate from moment-matrix rank. Even if each one-sided support has rank four and good conditioning, an unresolved copy map does not define a unique MLS field. Conversely, any `UNKNOWN_POOR_SUPPORT`, `UNKNOWN_RANK_DEFICIENT`, `UNKNOWN_ILL_CONDITIONED`, `UNKNOWN_SUPPORT_BOUNDARY`, `UNKNOWN_VISIBILITY_SWITCH`, `UNKNOWN_DUPLICATE_SOURCES`, or non-OK visibility result remains `UNKNOWN`; no pseudoinverse, regularization, all-visible fallback, or zero-filled failed row is allowed.

## Essential front-shadow continuity read

For a finite facet, the Boolean rule changes a source’s mask on the set of queries whose segment to that source grazes a facet edge/front. This shadow sheet extends from the physical front through the surrounding body. A local derivative on either side does not account for the jump in the MLS row at the sheet.

The mechanism is direct: in the MLS moment matrix, visibility multiplies each nonzero compact-kernel contribution. When source `s` switches at a query with `‖x-s‖ < r₀`, its contribution changes by `w(‖x-s‖/r₀) p(s)p(s)ᵀ` (up to sign); the right-hand-side/shape row changes too. Wendland C2 is positive for this strict interior-support case. There is no stated identity forcing the two MLS rows to match across that change, so a nonzero jump is the generic expectation, not an observed result. The local derivative formula differentiates `A` and `B` only while the mask is fixed; it contains no derivative of the Boolean switch.

Use this one fixed geometric path for the first read. Let the base front edge be `AB`, its midpoint `p=(0,0,-1/2)`, and the ordinary uncracked source `s=(1,0,0)`. Define

```text
q₀ = p + (1/16)(p-s) = (-1/16, 0, -17/32)
Π  = plane through s and AB: z = (x-1)/2
nΠ = (-1/2, 0, 1)
q(t) = q₀ + t nΠ.
```

`q₀` is strictly outside Γ and away from its front (its distance to `p` is `√5/32`). The segment `[s,q₀]` passes through the front midpoint `p`. For sufficiently small positive and negative `t`, the segment crosses opposite sides of the front edge, so the visibility of `s` switches while the query path remains outside the physical crack. The pair is active only if the frozen baseline support radius includes `s` at `q₀` (`‖s-q₀‖ = √1445/32`); if not, record `UNKNOWN_INACTIVE_SHADOW_SOURCE` and do not search over another source or path.

Obtain the two one-sided masks `M⁻, M⁺` from exact-sided geometric classification, with every other source’s status recorded. Do not evaluate the exact grazing query as an ordinary `OK` visibility row. Instead evaluate the one-sided limits of the fixed-mask MLS rows at `q₀`, using each side’s frozen mask and the same source/DOF map:

```text
Φ⁻ = lim(t→0⁻) Φ(q(t); M⁻)
Φ⁺ = lim(t→0⁺) Φ(q(t); M⁺)
JΦ = ‖Φ⁺ - Φ⁻‖₂.
```

Each frozen-mask continuation must independently pass MLS support/rank/conditioning and radial-support-boundary checks. The coefficient-row difference is the supremum displacement jump over unit coefficient vectors; also report the jump for fixed unit vectors on the switched source and split center copies. A nonzero `JΦ` outside Γ demonstrates a nonphysical shape discontinuity. If `JΦ=0` for this path, that only passes this falsification probe; it does not prove continuity over all shadow sheets.

If `JΦ ≠ 0`, the current ordinary `∫ BᵀDB dΩ` assembly integrates only the classical gradients on each fixed-mask region. It omits the distributional derivative/interface contribution of that shape jump and supplies no constraint or flux coupling on Π. The discrete trial/test space is then broken on a virtual surface outside Γ. A `K` change, hidden-weight-zero check, good local rank, symmetry, or six-mode rigid check cannot repair this mathematical defect or establish traction-free behavior on the intended crack. Stop the physical-coupling claim and report `BLOCKED_NONPHYSICAL_SHADOW_JUMP` to root. The reviewed formulation and primitives establish finite-facet visibility and split support, but do not establish a support/diffraction correction, enrichment, or continuous visibility ramp that removes or consistently treats these extra sheets; those prerequisites remain `UNKNOWN` for root.

If the path itself is `UNKNOWN` because of visibility ambiguity, coincident-copy mapping, an inactive source, or poor support, that is not a pass. The correct result is `UNKNOWN_SHADOW_CONTINUITY`, with the reason. Any method needed to remove or account for the shadow jump (for example a proven continuous visibility ramp/partition, or explicit treatment of every induced interface) is outside the supplied primitives and remains a root-level mathematical decision; this note does not select an algorithm.

## Operator evidence required after the mathematical gate

Only if all shape rows and the shadow-continuity precondition are resolved, assemble the same small-strain operator from the actual per-query MLS gradients:

```text
KΓ = Σq BΓ(q)ᵀ D BΓ(q) wq,       K0 = Σq B0(q)ᵀ D B0(q) wq.
```

There is no facet stiffness, cohesive penalty, contact, or added `Γ` term. `D` is unchanged. The finite facet is the only changed geometric input; no bond deletion, PD horizon, CriticalStretch, post-fit mesh cut, or output reconstruction participates.

Record all of the following before using “mechanically coupled” for the **discrete operator**:

1. **Hidden weights:** every blocked or out-of-support source has exactly zero shape value and gradient; no `UNKNOWN` is turned into zero. For every successful row, record support count, condition number, `ΣᵢΦᵢ=1`, `ΣᵢΦᵢXᵢ=x`, and their gradient identities. Check affine reproduction separately on both crack sides and on visible paths around the front.
2. **Independent traces:** at matched points with projection strictly inside Γ, evaluate the plus and minus limits and demonstrate a nonzero jump for a fixed coefficient vector that gives the two center copies different values. Show that shape gradients and assembly IDs really use the matching copy. This establishes that a constant Boolean mask alone is not being mislabeled as a displacement jump; the duplicate DOFs provide the independent traces.
3. **Release across the face:** for quadrature rows whose support segments strictly cross Γ’s relative interior, no local `B(q)` may contain both members of a split pair. Record per-quadrature support IDs and local stiffness contributions; report any cross-pair contribution explicitly. Do not infer zero global cross-blocks from a graph picture, since unsplit DOFs can carry legitimate finite-front coupling.
4. **Finite-front connection:** the off-facet shadow path must have equal one-sided field traces (see the gate above). The relative opening/sliding vector of the split pair must have finite, positive energy through the surrounding ligament; it must not become an additional rigid mode. These are discrete connectivity checks, not a stress-intensity or crack-tip accuracy claim.
5. **Null modes and matrix checks:** use the same predeclared eigenvalue/residual tolerance as the uncracked certificate. The finite internal triangle leaves one connected body, so `KΓ` should have the same six physical rigid modes as `K0`; any additional exact/near null mode is `UNKNOWN_EXTRA_NULL_MODE` unless independently explained and approved. Check finiteness, symmetry, and semidefinite spectrum. The existing six-rigid-vector diagnostic alone does not prove that the six modes are the only null modes.
6. **Surface-dependent stiffness:** publish `ΔK = PᵀKΓP-K0`, where `P` maps the two center-copy coefficients to the same intact parent coefficient. Also publish split-pair blocks and each `q` contribution. Do not require `ΔK=0`, and do not treat `ΔK≠0` as sufficient. Any claimed mechanical effect must be traceable to a changed finite-facet visibility row and survive the shadow continuity gate.

The current EFGM module assembles bounded dense K only; it has no load, constraint, solve, or traction-output path. This pilot therefore stops at the operator and trace evidence. A loaded static solution and a traction residual would be a separate next stage after root accepts the space and front treatment. In a conforming piecewise-`H¹` formulation, traction-free Γ is the natural boundary condition when no crack-face term is applied; it must not be asserted from a stiffness delta or missing cross-bonds alone.

## Candidate B — escalation only if A is support-limited

If candidate A blocks because only one interior particle is split, either side lacks rank, or the front test cannot be resolved, do not enlarge supports or alter tolerances. Ask root whether to replace the fixed 27/64 smoke case with a regular `5×5×5` layout at coordinates `{-1,-1/2,0,1/2,1}³`, retaining `Ω`, the `x=0` facet plane, material, and MLS kernel/configuration in physical units. Use this one fallback triangle, which contains seven lattice sites strictly in its interior:

```text
A₂ = (0, -9/10, -3/5),   B₂ = (0, 9/10, -3/5),   C₂ = (0, 0, 9/10).
```

This gives the front more than one adjacent interior site. Candidate B is a second fixed candidate only; there is no facet or support sweep.

That is a different discretization and not a same-`K` comparison. The matching two-Gauss-per-cell volume rule has 512 points; the split layout has 132 sites/396 DOFs if all seven interior facet particles are duplicated. Both exceed the present elastic operator caps (`256` quadrature points and `375` DOFs). Candidate B is therefore only a principled escalation proposal for root to decide; it requires a separately approved bounded assembly budget and a matching uncracked baseline. Do not silently lower quadrature, raise a cap, or run it under candidate A’s acceptance criteria.

## Fixed sequence and stop outcomes

1. Consume, but do not edit or rerun, the separate 27-node/64-point uncracked MLS/81-DOF and exact mask-stability artifacts. Freeze their config and tolerances.
2. Precompute the single stated facet, side-copy offsets, finite visibility classifications, and front-shadow one-sided masks. Any visibility or copy-map ambiguity is `UNKNOWN`.
3. Evaluate affine identities, split-side trace jump, and the off-facet shadow limits. A nonzero shadow jump blocks the ordinary volume-only Galerkin interpretation before K is assembled.
4. Only after the space gate passes, assemble `KΓ` on the exact same 27/64 points and produce the six-mode, extra-null, cross-face, front-connectivity, and `ΔK` records above.
5. Stop. No load solve, SIF, growth step, R2 label, full-vessel result, material calibration, physics claim, or model-GT promotion is included.

The source contract’s frozen inclined penny-crack/SIF target and unknown kernel/support provenance remain separate. This triangle is only a coupling falsification geometry; it is not the paper benchmark or a ground-truth crack.

## Exact source snapshots reviewed

The following SHA-256 values identify the source files read for this design. These values are provenance for this document, not claims that those sources implement the proposed coupling.

| File | SHA-256 |
|---|---|
| `pgfrag/experiments/explicit_crack_formulation_20261005/FORMULATION_CONTRACT.md` | `56e36e197985e4b829b95090f224e67c7b7a94bb6c977961b1dbb45f21b681eb` |
| `pgfrag/experiments/mls3d_kernel_20261005/mls3d_mls.py` | `404256d9b7546361a5704c4c3cec7c022069a9cf83a137c58a16e7a64f79ec63` |
| `pgfrag/experiments/mls3d_kernel_20261005/mls3d_shape.py` | `66e680616a988182875d17aeb84fadb541867aee13fc15e0ae098f7c98839ef5` |
| `pgfrag/experiments/explicit_visibility_20261005/visibility3d.py` | `9d5097ca33ce899cf788ad43ef0740a34296216a8fd7aa490868c27ffe1ebb03` |
| `pgfrag/experiments/efgm_elastic_operator_20261005/efgm_elastic_operator.py` | `41a5ac6ef6359ed1640ee681c94a3496c18652b29cc72c679b3e37c20dd9ac0d` |
| `pgfrag/tmp/review_20261005/astra_crack_surface_architecture/DECISION.md` | `284c444662960b2f15d8ae99d61aed2b25091048ff3c104dce72d2139f6d3b47` |
