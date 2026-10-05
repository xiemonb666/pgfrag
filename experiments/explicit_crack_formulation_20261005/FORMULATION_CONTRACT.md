# Explicit 3D crack formulation contract (research draft)

Status: source-first design only. No solver, test, CAD, or production code was run or changed. This is a new quasi-static LEFM/EFGM cracking-particle backend, not an extension of the existing OSB/CriticalStretch mechanics. OSB remains an event diagnostic.

The selected reference is Ai & Augarde (2018), DOI [10.1002/nme.5786](https://doi.org/10.1002/nme.5786), especially §3 and §4.2. The source evidence, page map, hashes, and unresolved source conflicts are in the corresponding review directory. The method explicitly carries a crack/front surface and changes EFGM supports through visibility. It is not the current peridynamic horizon and has no demonstrated cohesive/ceramic process-zone law.

## Benchmark contract

Use the §4.2 inclined penny-shaped crack as the first mechanics-facing target. Parameters confirmed in the author's official thesis, Ch. 6 §6.4.2 (print pp. 139–141; PDF pp. 164–166), which reproduces the same example:

| Item | Source value |
|---|---|
| Cubic domain | Side `2w = 200 mm`, centered in the figure |
| Initial crack | Circular disk, radius `a = w/5 = 20 mm`; represented by a circular crack cell for visibility |
| Inclination | `β = π/6`, from the x-axis in the xz plane; preserve the source convention and verify the figure orientation during implementation |
| Material | Isotropic, linear elastic; `E = 100 GPa`, `ν = 0.3` |
| Load | Uniaxial tension `σ = 100 MPa` on the upper and lower faces, other faces free (thesis §6.4.1 plus “same ... boundary conditions” in §6.4.2) |
| Initial point layout | Regular `15 × 15 × 15`; adaptive refinement follows in the published method |
| Volume quadrature | `14 × 14 × 14` background cells, `4 × 4 × 4` Gauss points per cell |
| Loaded-face quadrature | `14 × 14` cells, `8 × 8` Gauss points per cell |
| Rigid modes | Six integral constraints via Lagrange multipliers: three mean translations and three mean rotations, source Eqs. (6.7)–(6.8) |
| Adaptive target | Global relative error `η_t = 0.02` for the inclined case |
| Growth exponent | `m = 3.32` for the examples, from source Eq. (6.5) discussion |

The accepted-article web extraction renders the load unit as `N/m²`; the author's official thesis states `100 MPa` for the baseline penny crack and the inclined case inherits that boundary condition. Treat `100 MPa` as the author-thesis-supported benchmark value and retain the extraction discrepancy in provenance. A byte-verified accepted-article PDF and a visual check of that article page were unavailable in this pass.

The reference mixed-mode SIFs along the initial front are the infinite-body solutions cited by the paper (article Eqs. (20)–(22); thesis Eqs. (6.11a–c)):

```text
KI(θ)   =  2 σ cos²(β) √(a/π)
KII(θ)  = -4 σ sin(β)cos(β)/(2-ν) √(a/π) cos(θ)
KIII(θ) = -4(1-ν) σ sin(β)cos(β)/(2-ν) √(a/π) sin(θ)
```

`θ` is measured around the circular front from the x-axis in the xy plane. The cube is finite, so the thesis explicitly reports boundary-effect deviations near `θ = 0, π`; compare the full front profile and report that limitation rather than imposing a pointwise exact match. The paper provides curves, not a numeric SIF table or a predeclared acceptance tolerance.

`Δa_max` and the interaction-integral box dimensions `c, l_c` are not reported as numeric values in the primary source. They remain `UNKNOWN`; do not choose them and label the result a paper reproduction. The source also leaves the exact MLS kernel and support radius used in §4.2 unspecified. The paper’s extracted unit and the source's Eq. (8) mode-III coefficient require the explicit provenance notes below.

## Mechanics and MLS kernel

The backend solves the static small-strain isotropic elastic weak form:

```text
ε(u) = sym(∇u)
σ(u) = 2 μ ε(u) + λ tr(ε(u)) I
∫Ω ε(v) : C : ε(u) dΩ = ∫Ω v·b dΩ + ∫Γt v·t̄ dΓ
μ = E/[2(1+ν)],  λ = Eν/[(1+ν)(1-2ν)]
```

Use MLS EFGM shape functions for displacement and their analytic gradients for strain and assembly. The author's thesis gives the linear 3D polynomial basis (Ch. 2 §2.2, Eq. (2.1)):

```text
p(x) = [1, x, y, z]ᵀ
J(a;x) = Σᵢ wᵢ(x) [pᵀ(xᵢ)a(x) - uᵢ]²
A(x) = Σᵢ wᵢ(x) p(xᵢ)pᵀ(xᵢ)
B(x) = [w₁p(x₁), …, wₙp(xₙ)]
Φᵀ(x) = pᵀ(x) A⁻¹(x) B(x)
```

Keep all derivatives of the moment matrix and weight vector. From thesis Eq. (2.10), for coordinate `x_j`:

```text
∂ⱼΦᵀ = (∂ⱼpᵀ)A⁻¹B
       + pᵀ[-A⁻¹(∂ⱼA)A⁻¹B + A⁻¹(∂ⱼB)]
∂ⱼA = Σᵢ (∂ⱼwᵢ) p(xᵢ)pᵀ(xᵢ)
∂ⱼB = [ (∂ⱼw₁)p(x₁), …, (∂ⱼwₙ)p(xₙ) ]
```

The 2018 paper identifies MLS/EFGM and a linear 3D basis size (4×4 moment matrix), but does not state the §4.2 kernel or support scale. The author's thesis §2.2 provides this compact quartic spline as an EFGM example and says a spherical 3D support radius is typically 2–4 times average point spacing:

```text
rᵢ = ||x-xᵢ||/r₀
w(rᵢ) = 1 - 6rᵢ² + 8rᵢ³ - 3rᵢ⁴, 0 ≤ rᵢ ≤ 1; 0, rᵢ > 1
```

The quartic spline and radius range are author-source references, not confirmed §4.2 settings. The standalone [3D MLS primitive](../mls3d_kernel_20261005/README.md) selects Wendland C2 with explicit support configuration as the root-approved engineering kernel. Keep that implementation choice separate from the author's quartic example; neither establishes the kernel/support used in the paper benchmark. A separate author 2019 2D CPM paper uses `d_s = 2.2`; that value must not be silently transferred to this 3D example.

## Explicit crack state and visibility

The surface is authoritative mechanical state, stored in the reference/material frame:

```text
CrackSurface: surface_id, reference vertices/facets, front edges/nodes,
              oriented side normal, creation/update lineage, status
CrackedParticle: parent_particle_id, original x, side ∈ {−,+}, side DOF,
                 incident facet IDs, visibility-only offset metadata
```

For each MLS evaluation point and candidate support particle, perform a segment-versus-finite-crack-surface visibility query. If the connecting segment intersects the actual crack patch, exclude that particle from the opposite-side support; an intersection beyond the finite front is visible. Use the local triangulated segments, not an infinite plane. The source describes one circular crack cell for initial visibility, then triangular facets for the advancing surface. Particles inside a surface triangle use the circular discontinuity; edge and vertex particles use the incident nonplanar segments. Store at most eight incident facets per cracking particle, as in the paper. Front particles are not split; their influence reaches both sides.

Cracking particles away from the front have coincident `+`/`−` copies with separate displacement DOFs and side-specific visible supports. The displacement jump at `x` is:

```text
[[u(x)]] = Σᵢ∈S+ Φᵢ(xS+) uᵢ+  −  Σᵢ∈S− Φᵢ(xS−) uᵢ−
```

The selected accepted article specifies an offset of `0.01 × support size` solely to distinguish coincident copies during visibility/connectivity checks; shape functions use original coordinates. The author's thesis Ch. 6 gives `0.001 × support size`. Keep the accepted-article value for the selected source, record the thesis discrepancy, and do not perturb geometry or shape-function coordinates.

For a static crack the faces carry no cohesive law: discontinuity is introduced by split supports and visibility in the linear elastic weak form. Cohesive fracture, contact, branching, coalescence, and dynamics are outside this paper’s demonstrated formulation.

## SIFs and one propagation update

For each front particle, form local axes from its incident front segments and triangular facets (source Eq. (9)): area/length-weighted front extension `x̂`, area-weighted crack normal `ŷ`, orthogonalize `ŷ` against `x̂` and normalize, then `ẑ = x̂ × ŷ`. Evaluate the interaction integral over a front-centered Cartesian box with dimensions `2c × 2c × 2l_c`:

```text
q = (1-|x|/c)(1-|y|/c)(1-|z|/l_c)
I = (1/l_c) ∫V (σauxᵢⱼ uⱼ,1 + σᵢⱼ uauxⱼ,1
                - σauxⱼₖ εⱼₖ δ₁ᵢ) q,ᵢ dV
```

The auxiliary elastic fields and source conversion are in article Eqs. (7)–(8) / thesis Eqs. (6.1)–(6.2), with fields in thesis Appendix C. **Source conflict:** Eq. (6.2) visibly prints `+ μ KIII KIIIaux` (multiplication), although `μ` is defined as shear modulus and Eq. (6.4)’s equivalent-SIF coefficient implies inverse-`μ` energy scaling. Preserve that literal in provenance; do not claim the article used the corrected term. For the independent implementation, use the energy-consistent bilinear interaction under isotropic LEFM,

```text
G(K) = (1-ν²)/E (KI² + KII²) + KIII²/(2μ)
I(K,Kaux) = G(K+Kaux) - G(K) - G(Kaux)
         = 2(1-ν²)/E (KI KIaux + KII KIIaux)
           + KIII KIIIaux/μ
```

This is an explicit derivation for the new primitive, not a transcription of the source. Require an analytical pure-mode-III interaction test before connecting this coefficient to the solver. `c` and `l_c` are `UNKNOWN`.

For an accepted propagation experiment, the article uses maximum circumferential stress for `θ_c` (Eq. (10)), ignores mode III for direction, combines modes with Eq. (12), then sets `Δa_i = Δa_max (ΔK_eq,i/ΔK_eq,max)^m` (Eq. (13)), `m=3.32`, and advances in the local frame (Eq. (14)). Connect old and new front with triangles. A nonplanar quadrilateral strip is split into three triangles around its midpoint; that midpoint is not advanced in the published simple algorithm. Run one front update only after static-SIF checks; do not invent the unknown numeric `Δa_max`.

The 2018 paper is quasi-static LEFM. Its Paris-law-scaled step is not a measured ceramic toughness criterion; it does not provide a nucleation rule, critical `K`, `Gc` calibration, or a process-zone model.

## Implementation-facing task interface

Keep functions pure where practical and make the frame and kernel configuration explicit:

1. `build_crack_state(geometry, front, lineage) -> CrackSurface`: creates/validates a circular initial patch, front, orientations, and surface ID in material coordinates.
2. `visible_support(query_x, particle, CrackSurface, side) -> bool`: finite-facet visibility with deterministic intersection tolerances and a recorded visibility-only split offset.
3. `mls_shape(x, support_ids, kernel_config) -> (Φ, ∇Φ, conditioning)`: exact source-backed linear basis; returns analytic gradients and failure status for rank-deficient support.
4. `assemble_elastic_system(points, volumes, quadrature, material, CrackSurface, BC) -> K,f,constraints`: sparse symmetric small-strain EFGM weak form; side DOFs are selected from visible support; use no OSB bonds or CriticalStretch deletion.
5. `solve_static(K,f,constraints) -> displacement, reactions, solve_record`: no time integrator.
6. `interaction_sif(front_particle, displacement, stress, auxiliary_fields, c, l_c) -> SIFRecord`: returns all terms, chosen coefficient convention, local axes, and quadrature record; rejects unresolved coefficient configuration.
7. `propose_front_advance(front, SIFRecords, Δa_max, m) -> AdvanceProposal`: pure geometry proposal; no mutation until accepted; write old/new front and generated facet lineage.
8. `export_surface_pair(surface_id, side_mapping)`: duplicate the authoritative internal interface for output while preserving one common area and two opposed sides. Do not infer mates from OSB bonds.

Current status is `NOT_COUPLED`: this research contract and the standalone MLS primitive do not create a mechanics-grounded surface, and there is no model R2. After a future backend couples the explicit surface to mechanics and passes surface-consistency and material-partition gates, that declared model output may qualify as `SIMULATED_GEOMETRIC_GT`, with `MATERIAL_UNCALIBRATED` retained until material validation. Surface/interface identity, front, fragment identity, event labels, and relations stay on the target/provenance side. Candidate ranking receives only observation-derived fragment features in the existing RelicWeave path; no GT, surface IDs, or relation fields enter features.

## Minimum testable slice and budget proposal

First implement only the 3D linear MLS value/gradient, finite-surface visibility, split-support displacement jump, and sparse static linear-elastic assembly. Begin with a no-crack affine patch/rigid-mode check; then a frozen inclined circular crack under dimensionless `2w=2`, `a=0.2`, `β=π/6`, `E=1`, `ν=0.3`, `σ=1e-3`. The ratio `σ/E=1e-3` matches the paper-supported `100 MPa / 100 GPa`, keeping the engineering fixture in a small-strain regime; its normalized SIF/reference scaling is unchanged. This fixture validates mathematical primitives only, not a paper reproduction, calibrated material fracture, or GT. Check affine reproduction away from the crack, a nonzero side jump with no cross-crack support contribution, force balance, symmetry/conditioning, and the initial-front `KI/KII/KIII` profile against Eqs. (6.11) while recording the finite-cube boundary deviations. Do not run propagation until the static SIF coefficient and convergence gate are resolved.

Suggested first execution cap, to be approved separately before any solver run: at most 4,000 original particles and 15,000 total displacement DOFs after splitting; at most 50,000 volume quadrature points and 13,000 points per loaded face; one static load solve; 10 minutes wall time and 4 GiB peak memory, stop on cap without retries. This is a development cap, not the paper’s error target. Full §4.2 reproduction still requires the unknown exact kernel/support, `c,l_c`, `Δa_max`, and a predeclared tolerance digitized from the published SIF figure; the `η_t=0.02` adaptive target alone does not define those settings.
