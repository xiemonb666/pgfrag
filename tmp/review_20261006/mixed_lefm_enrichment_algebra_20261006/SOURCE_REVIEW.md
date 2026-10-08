# Source review: mixed LEFM enrichment algebra

Status: source-only helper prepared for main-agent review; no numerical test was run.

## Contract

- Inputs are one stable `source_ids` tuple, row-aligned base values and world-frame base gradients, either world point plus tip or exact local coordinates relative to the tip, a proper crack frame `Q` whose columns are local axes in world coordinates, a positive reference length `R`, and an optional explicit branch side. The local route preserves exact cut traces without snapping or tolerances.
- Output source rows are sorted lexically by the real source ID. Each semantic coefficient ID is JSON for `[source_id, component, basis_kind, basis_index]`; no dense slot number is used as an identity.
- One source contributes five scalar functions (`BASE` plus four branches) times three vector displacement components, or 15 coefficients. `mixed_dof_count(27)` is 405. The helper has no 375-DOF cap so it cannot silently drop the five-function mixed basis.
- The six B-column entries follow the existing full 3-D engineering convention `[eps_xx, eps_yy, eps_zz, gamma_xy, gamma_yz, gamma_xz]`. For example, an `ux` coefficient with gradient `(gx, gy, gz)` maps to `(gx, 0, 0, gy, 0, gz)`.

## Algebra checked by inspection

For `F_a = q h_a(theta)`, `q = sqrt(r/R)`, the four angular functions and their product-rule angular derivatives match Moës, Dolbow, and Belytschko (1999), Eq. (16), pp. 135–136. With `q/r` scaling, the local gradient is

```text
dF/dx = (q/r) [0.5 cos(theta) h - sin(theta) h']
dF/dy = (q/r) [0.5 sin(theta) h + cos(theta) h']
```

so branch values and gradients use the same normalization. World scalar gradients are `Q grad_local`; enriched rows use `F grad(N) + N grad(F)` in world coordinates. The B columns directly encode the corresponding component-wise symmetric gradient with engineering shears.

At the front the function returns `UNKNOWN_FRONT` with empty coefficient arrays. At exact local `x < 0, y = 0`, the no-side path returns `UNKNOWN_BRANCH_FACE`; `upper/lower` explicitly select `+pi/-pi`. There is no tolerance or point shift. Improper frames, invalid selectors, malformed base rows, duplicate IDs, and nonfinite inputs are rejected.

## Source and integration boundary

The branch set is from Moës, Dolbow, and Belytschko, “A finite element method for crack growth without remeshing,” *International Journal for Numerical Methods in Engineering* 46 (1999), Eq. (16): [publisher DOI](https://doi.org/10.1002/%28SICI%291097-0207%2819990910%2946%3A1%3C131%3A%3AAID-NME726%3E3.0.CO%3B2-J), [author-uploaded paper](https://www.researchgate.net/publication/51992357_A_Finite_Element_Method_for_Crack_Growth_without_Remeshing). The basis normalization by `sqrt(R)` is documented in the helper README.

This remains an algebraic row builder. It does not call the existing solver or shape provider, build MLS values, determine support/visibility, or model fracture. The current elastic operator limit is 375 DOFs, while 27 sources with this basis require 405; operator-budget changes and any mechanical coupling need separate review.

See `FUTURE_TINY_VALIDATION_BUDGET.json` for the capped follow-up checks. Its status is planned, not executed.
