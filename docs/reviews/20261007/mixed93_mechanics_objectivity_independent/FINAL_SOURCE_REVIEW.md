# Final source review: Mixed93 objectivity adapter

**Disposition:** source-accepted for a coordinate-covariance diagnostic after one narrow record-construction fix. Numeric execution remains blocked by the case’s explicit authorization fields (`false`, with no root-go token). This review performed no numerical run.

**Identity:** reviewed the original 321-line adapter SHA-256 `f9b35bd271e1ac9e86bedb813435f3a0b9b045bc84f4b1e1e56c1f10244c346f`, case SHA-256 `a189527e639190192312dd8e9685d493f2b241aaed7680c2979e42a63390574e`, and matching `source_pins.json` SHA-256 `2815ceb641c56568057793a8f345d7751b43f43b4ecadda7530f3d422c320f69`. All 17 declared source-pin entries match. The cached NPZ was byte-hashed only; its archive/payload was not opened or parsed. The prior `REVIEW.md` and `review.json` remain unchanged.

## Confirmed coordinate-covariance contract

The adapter correctly uses Q=Rz(0.37)Rx(−0.23) and t=(0.003,−0.007,0.004) m. It co-transforms source, volume, and boundary points; normals; body/front/crack directions; crack axes, origin, and tip; MLS reference origin; and gauge center. It keeps scalar MLS scales, R, L, E, ν, and local KI/KII/KIII fixed. The provider receives the actual transformed world points and source cloud, with unchanged scalar quadrature measures.

With H storing world crack axes in columns, `(x−tip) @ H` gives the row form of crack-local coordinates. The adapter holds those local coordinates fixed across the rigid motion, evaluates the four branch functions there, and maps local branch gradients back with H′ᵀ in row form. The analytic field uses a transformed `CrackFrame`; its returned stress tensor is in world coordinates and produces traction as σ′n′. This fixes the earlier identity-frame assumptions.

Its engineering-strain conversion halves the three engineering shear entries to build ε; its stress conversion copies the three shear stresses directly. It rotates the resulting tensors with QTQᵀ. The coefficient, stiffness, constraint, load, and physical-multiplier comparisons use A=I₃₁⊗Q and B=diag(Q,Q): q′=Aq, K′=AKAᵀ, C′=BCAᵀ, f′=Af, and λ′=Bλ. The moment helper receives x−c, with c′=Qc+t, and the direct moments use the same relative coordinates. The wrench helper checks both centered M′c=QMc and global-origin M′₀=QM₀+t×F′, with force and moment scales based on the absolute traction integral. These formulas are correct for the pinned isotropic fixture.

## Confirmed runtime blocker

The adapter constructs `FrameState` with `mls_config` but no `mls_origin_world_m` field. `make_frame_state` computes the transformed MLS origin and stores it only inside `mls_config`; later, record construction accesses `x0["frame"].mls_origin_world_m` and `x1["frame"].mls_origin_world_m`. That attribute does not exist. A call that reaches record construction will raise `AttributeError` after the assemblies, spectrum gates, KKT solves, and comparisons, before returning the result. Add the origin to `FrameState` or serialize `mls_config.reference_origin`. This is a source-confirmed blocker; it was not exercised.

## Helper contracts and gates

The adapter’s private helper calls match the pinned source contracts: `_physical_moment_matrix(np, values, relative_xyz, weights, radius)` returns `(C, volume, zero_center)`; `_direct_physical_moments(np, displacement, relative_xyz, weights, radius, volume)` returns six moments; `_scaled_spectrum_gate(np, K)` returns `(Kscaled, scale_vector, eigenvalues, metrics)`; `_solve_case(np, K, Kscaled, C, f, scale_vector, g)` returns `(q, physical_lambda, solve_metrics, arrays)`; and `_reconstruct_field(np, values, gradients, q, D)` returns displacement `(n,3)`, engineering strain `(n,6)`, stress `(n,6)`, and batch sizes. The adapter passes the expected argument counts and uses the private solver’s restored physical λ.

Each fresh spectrum must pass with exactly six zero modes, 87 positive modes, and no significant negative mode before the KKT solve. The private solve then enforces scaled saddle backward error, physical weak residual, gauge residual, and energy identity at 1e−9; failure raises instead of returning a pass. The adapter performs two fresh raw K/C/f assemblies, two spectrum checks, and two 99×99 solves. Its workload arithmetic is consistent: 1,280 MLS rows/state in five 256-row chunks, 2,560 rows total, and 1,536 analytic boundary-field evaluations. No spectrum is reused.

One scope limit remains: the fixture’s gauge target is zero in both states, so `g′=Bg` is satisfied trivially here. Nonzero target covariance is not exercised by this diagnostic. The adapter’s capture function also delegates source-pin and dependency-origin verification to its caller, as its docstring states; the outer runner must verify those pins before calling it.

This is still only a proposed coordinate-covariance diagnostic. Relative source-cloud orientation, quadrature-rule orientation, tetrahedral mesh direction, crack-front/path bias, and material or crack-growth qualification remain UNKNOWN. Co-rotating the cloud and quadrature does not resolve those questions. No provider calls, field evaluations, assemblies, solves, tests, or NPZ payload access occurred in this review.
