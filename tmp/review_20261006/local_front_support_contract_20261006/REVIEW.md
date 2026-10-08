# Support-topology certificate review

Status: `SOURCE_ONLY / GEOMETRY_CERTIFIED / NOT_COUPLED`. This review covers only finite support topology and the fixed source-role manifest. It does not implement or qualify an MLS provider, operator, solver, crack growth, or the local finite-element approximation.

## Decision

The independent certificate assigns exactly 100 source IDs to the accepted Alternative A geometry: 30 `FRONT_F1_SHIFTED`, 10 `HEAVISIDE`, and 60 `BASE` only. Every source has BASE and at most one enrichment, for 140 scalar DOFs. All source coordinates use metres, all source y coordinates are nonzero, and IDs bind to lattice indices and coordinates. The frozen role manifest is pinned by the fixed-case bytes, classifier-source bytes, record-set digest, and full-manifest digest.

The certificate is finite-geometry based. For each source it evaluates distance to the finite front segment F and finite crack rectangle Γ, then reasons about the convex clipped support `C_i = Ω ∩ closed_ball(X_i,h)`. It never treats crossing the infinite plane `y=0` as sufficient for H. H is certified only when the full box-clipped section `C_i ∩ {y=0}` lies inside the finite Γ, both open y-side pieces are nonempty, and C_i has positive clearance from F. Convexity then proves exactly two connected components after removing Γ. Front/cut tangencies, finite-end overruns, partial crack intersections, and source-side values undefined for shifted bases return `UNKNOWN`.

For the ten H records, the source-to-front distance is `√5 s`; the ball front clearance is `(√5−1.75)s = 1.9442719 mm`. Their y=0 plane sections have maximum x `−0.5638593s`, or 2.2554374 mm before the front; supports at z endpoints are clipped by Ω but their entire plane sections remain in Γ. The two-side margin is `h−|y|=3 mm`. The front-core bound from the nearest non-front source is `0.25s=1 mm`. Every crack-face point is within `√1.5 s = 4.8989795 mm` of an enriched support center, leaving 2.1010205 mm support margin.

The crack-face values `u_z^±=±A√(-x)` and zero scalar traction recorded in the prior design are continuum reference targets. This role manifest makes no claim that a finite local approximation reproduces those traces exactly outside the proven front core.

## Ramp alternative B is not accepted

The normalized ramp alternative in the prior design is **not accepted**. It divides by `W_F=Σ_{i∈T_F} φ_i`, but MLS shape values can be negative and support coverage alone does not prove that W_F is positive or bounded away from zero. Before any use, provide a certified positive lower bound for W_F on the entire proposed window, including signed MLS weights, and show its derivatives remain bounded through the taper. A support-center coverage argument is insufficient. The original `DESIGN.md` is unchanged; this review records the later rejection of alternative B.

## Artifacts

* `pgfrag/experiments/local_front_support_contract_20261006/support_topology.py` contains the independent geometry certificate and manifest builder. It imports only Python standard-library modules; it imports no MLS, provider, or solver code.
* `pgfrag/experiments/local_front_support_contract_20261006/fixed_case.json` freezes Ω, finite Γ/F, h, topology tolerance, lattice, counts, and bounds.
* `pgfrag/experiments/local_front_support_contract_20261006/role_manifest.json` stores each source ID, lattice index, coordinate, basis role, distance/clearance values, plane-section evidence, and split-component proof.

The module's first CLI use materialized the manifest. No unit or regression test, MLS query, provider, assembly, or solve was run. Before running one primitive focused test, root should read `support_topology.py`, `fixed_case.json`, and `role_manifest.json`.

## Targeted regression definitions (not executed)

1. **Manifest binding and counts.** Verify fixed-case and classifier file pins, source-record digest, manifest digest, IDs-to-lattice-coordinate binding, unique IDs/coordinates, one-enrichment maximum, y≠0 on every shifted source, and counts `30/10/60`, `140`.
2. **Finite crack versus infinite plane.** At `(+2s,+s,0)`, the support intersects the infinite plane but misses finite Γ: expect `BASE`, never `HEAVISIDE`.
3. **Front and crack tangencies.** A source at `(h,0,0)` and a source at `(-2s,h,0)` are respectively tangent to finite F and finite Γ: expect `UNKNOWN`. Perturb either distance into the tolerance band: still `UNKNOWN`.
4. **Finite endpoint overrun.** Use a test geometry whose crack z interval is strictly shorter than Ω and a non-front support whose plane section runs beyond that interval: expect `UNKNOWN`, not H.
5. **Box-clipped endpoints.** For manifest H sources at `z=±2s`, verify the recorded plane section is clipped to Ω, remains inside finite Γ, has two nonempty sides, and certifies exactly two components.
6. **Shifted-source values.** Every manifest F1/H source has nonzero y. A front/cut source with an otherwise eligible support must preserve geometric evidence but return a basis-role `UNKNOWN` for the shifted value; no side guess or epsilon offset.
7. **Permutation and tampering.** Reordering source input rows while preserving IDs/index/coordinates produces the same canonical record digest. Swapping an ID's coordinate or lattice index fails. Editing a role/coordinate without updating pins fails binding.
8. **Finite-front priority.** A point ahead of the tip whose support crosses the infinite plane but remains wholly at x>0 must not gain an H role. A support that reaches the actual front takes the F1 role before H eligibility.

**First root check:** after reading the sources, run one primitive source-only smoke check that loads the frozen case and manifest, verifies file/record bindings, and confirms one front record, one endpoint-clipped H record, and the finite-ahead BASE record. Do not run an operator or solver as part of that check.

## Source pins

Hashes are SHA-256 over exact file bytes. The prior design and DOF contract hashes are included to bind the review context; the original design remains unchanged.

| File | SHA-256 |
|---|---|
| `support_topology.py` | `EA20E51D2B774C7E6B8B48403359B56BAB1BA6F5405DBA41F34813B6AEC4578E` |
| `fixed_case.json` | `0BDBB20032E105851D2909FF54850EBC54D3F2041C25AB5E42E9102244D90CA5` |
| `role_manifest.json` | `E5E318F80B5DFEDBC8BA122F82DF3CCA1585E604272399D21BC35907D1D28F01` |
| Prior `DESIGN.md` | `6DF2E9384D91886E8DC5AA5165082C78081E4184090901B491433360CFB18532` |
| Frozen DOF `CONTRACT.md` | `3C4110DD6D2A6804C45FA1FE938592D785FD6CF3B84B83343C14492E320AFAB3` |
