# Local-front provider review

This is a source-only provider check for the single frozen Mode III local-front case. It adds a dense 140-ID antiplane row over the existing continuous MLS primitive and frozen support-role manifest. No operator, stiffness matrix, solver, CAD model, or crack-growth logic was run or added.

The provider refuses changed geometry, classifier, manifest, registry, MLS, and manufactured-field files by raw SHA-256. It validates manifest bindings and uses the manifest sequence for source coordinates. The registered rows are 100 BASE, 30 `FRONT_F1_SHIFTED`, and 10 `HEAVISIDE` keys. Each source has BASE and at most one enrichment. The MLS input has all 100 sources and no visibility mask. Registry IDs are identical when keys are registered in reverse order.

For each source, the provider maps `phi_i` to BASE; `phi_i*(H(q)-H(X_i))` to HEAVISIDE; and `phi_i*(Fbar(q)-Fbar(X_i))`, with `Fbar=F1/sqrt(R0)`, to shifted F1. The front gradient uses the full product rule. Heaviside gradients are the one-sided open-domain derivatives. PLUS/MINUS crack traces are evaluated at exact `y=0` with no coordinate offset. The provider does not divide by a subset sum of signed MLS weights. Geometry preflight returns empty UNKNOWN rows at the front, missing/ambiguous traces, finite crack ends, and out-of-domain points. MLS support-boundary status also remains UNKNOWN.

The 242-point batch used 236 MLS queries; six geometrically invalid rows were rejected before MLS, and one additional row was UNKNOWN at an MLS support switch. The result has 235 OK rows and 7 empty UNKNOWN rows. It completed in 0.396 seconds. Budgets were 100 sources, 140 scalar DOFs, and 242 query records, within 125/375/256 caps.

The focused checks measured:

- BASE affine reproduction over 235 OK rows: maximum value error `9.71e-15`; maximum gradient-component error `3.49e-11`.
- Shifted-F1 core reproduction at four diagonal points, each `r=0.0006364 m`: maximum field error `1.69e-21 m`; maximum gradient error `2.17e-18`. The observed front BASE partition sums were within roundoff of one, but no partition normalization was applied and this sample does not prove exactness elsewhere.
- At three interior crack-face positions, the base rows have no PLUS/MINUS jump; the shifted-F1 and Heaviside per-key jumps match their direct trace identities. The maximum Heaviside identity residual was zero. With unit coefficients on all ten H rows, the observed row-sum jumps were `2.0`, `1.0`, and `0.08576` at `x=-0.008`, `-0.006`, and `-0.004 m`.
- At `(0.004,0,0)`, ahead of the finite front, no explicit side was requested; the provider returned OK and every H value and gradient was exactly zero. The point was not classified as an infinite-plane crack trace.
- The whole-box `z=0` section used 224 direct samples, all OK. Its reported local-front reconstruction errors were: mean absolute `9.46e-6 m`, RMS `1.22e-5 m`, 95th percentile `2.43e-5 m`, and maximum `2.53e-5 m`. These are observed collar/outside-core errors, not acceptance targets for exactness.

The PNG shows the exact Mode III section, the local shifted-F1 PU reconstruction, and their absolute error. The finite crack and front are marked; the figure is based on the saved direct query samples. Its SHA-256 is `46c02dbd9344d975e46f5d26fd35ef85615d18c8a07ba7558bbde3593555713c`.

The source-topology controls include front and crack tangencies, a near-tangent front support inside the tolerance band, and an outside-domain source. They returned UNKNOWN. A manifest H source whose support-plane section is clipped at the closed box end was rechecked and remained certified; the recorded clipped `z` bound is exactly `0.008 m`. Query rows exactly at finite crack ends still return UNKNOWN. These are checks of this exact fixture and pinned classifier only.

Root's independent topology review identified two tolerance-band weaknesses in the classifier when used as a generic arbitrary-geometry API: a barely overhanging finite-end section can be certified, and a near-front source just outside the box can be classified as FRONT without clipping. This provider does not claim general geometry support; it binds to the exact fixed case/classifier/manifest hashes, whose frozen source clearances are well above those bands. The end-source check is deliberately limited to the pinned box-clipped endpoint record.

Numerical results and all raw dependency pins are in [`provider_diagnostic.json`](../../../experiments/local_front_provider_20261006/provider_diagnostic.json); the reviewed 2-D section is [`local_front_section.png`](../../../experiments/local_front_provider_20261006/local_front_section.png). Main implementation and runner are [`local_front_provider.py`](../../../experiments/local_front_provider_20261006/local_front_provider.py) and [`run_provider_diagnostic.py`](../../../experiments/local_front_provider_20261006/run_provider_diagnostic.py).
