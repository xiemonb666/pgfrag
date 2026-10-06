# Enriched-face Mode-I cohesive trial primitive

**State: `SOURCE_ONLY_NOT_COUPLED`.** This is a bounded experimental trace/traction primitive. It is not a solver, crack-growth rule, vessel coupling, or ceramic calibration.

`EnrichedFacePoint` carries a frozen ordered tuple of semantic coefficient IDs and paired dimensionless trace operators `H+`, `H-` with shape `3 × ndof`. The coefficient amplitudes are in metres, so the common-point jump is `J q = (H+ - H-) q` in metres. Both operators must declare the same world reference point. A point has exactly one positive `one_side_area_weight_m2`; the API has no second face area to add accidentally. Duplicate `(interface_id, quadrature_id)` keys and history-key mismatches fail closed.

For a unit normal `n`, opening is `δ = n · Jq`. A point contributes

`f = A Jᵀ n t` and `K = A (Jᵀ n) k_alg (Jᵀ n)ᵀ`.

With metre-valued generalized coefficients, force is N and tangent is N/m. Traction is Pa, area is m², and stored energy, dissipation, and work are J. All state updates and tractions come from the pinned scalar law. The algorithmic derivative is reused from the private `_algorithmic_tangent` helper in the pinned nodal experiment; this shares the already reviewed branch convention without copying a second damage law. Softening may produce a negative tangent, so the result makes no SPD claim. Trials read the caller's immutable committed state and never commit it. Kinks report nondifferentiability and the selected one-sided branch.

The caller must provide `RigidTraceCandidates` in frozen `Tx, Ty, Tz, Rx, Ry, Rz` order and, at each QP, six expected point displacements. Both face operators must reproduce each caller-declared displacement and `J` must annihilate it. This deliberately uses no node coordinates as proxies for enriched coefficients. Passing verifies only the six vectors the caller supplied; it does **not** establish that they are the physical rigid modes of the actual 93-coefficient front spectrum. The mandatory evidence is required for both `BASE` and `GLOBALFRONT` tags; those tags describe the caller's trace family, not a generated basis.

The module accepts at most 93 unique semantic coefficient IDs and 16 QPs, checking these limits before allocating the dense tangent. It supports only positive Mode-I opening. Compression, tangential slip, and mixed mode are `UNKNOWN` and rejected. The illustrative `case_globalfront_toy.json` is only a future fixture manifest. That 24-coefficient future fixture is evaluated at a non-origin point, where its six base translation/rotation candidates have nonzero rotation traces; three enriched columns form `J = n ⊗ a` for `a = (0.5, -0.25, 1)`. The unrun micro-checker compares three forward finite-difference columns with the assembled tangent and independent `A t a` / `A k_alg a aᵀ` expectations, while preserving a 24-call ceiling.

Source dependencies are pinned to scalar law `78433638608ad4216979607412fe0d4966f1423706c47f734f50dd2a8121fa17` and nodal derivative source `d2bf0dd79507e5ab26efe69b5489c3284b7c9118656ab49fffe942303ba11854`. Neither dependency, nor this primitive or the future micro-checker, is executed by the static contract check.
