# Narrow source fix acceptance

**Disposition:** accepted for the source-only coordinate-covariance diagnostic. Numeric execution remains unauthorized; mesh, source-cloud, quadrature, and crack-bias findings remain UNKNOWN.

**Frozen identity:** original 321-line adapter snapshot SHA-256 `f9b35bd271e1ac9e86bedb813435f3a0b9b045bc84f4b1e1e56c1f10244c346f`; current 323-line adapter `6ca66a86dcf801464afc1104d69b55b3421aadda6860482bc898170ce5c253ce`; unchanged case `a189527e639190192312dd8e9685d493f2b241aaed7680c2979e42a63390574e`; refreshed source pins `367ca3c1dfd694a760941b42bb2b54b5f00a64560978cd11f5f1f6413d874178`.

`fc.exe /N` confirms exactly two inserted lines and no other source differences: `mls_origin_world_m: object` in `FrameState`, and the already computed `mls_origin` passed into its constructor. This supplies the attribute used by record serialization; no transform law, arrays, math, gate, or control flow changed. The original snapshot remains preserved.

All 20 current pin entries match. The CASE hash is unchanged. The cached NPZ was checked by raw-byte hash only; its payload was not opened or parsed. The prior `FINAL_SOURCE_REVIEW.md/json` remain unchanged.

The helper contracts and fresh spectrum/KKT gates remain as documented in the prior final review. Case counts remain 2,560 MLS rows in five 256-row chunks per state, 1,536 analytic boundary evaluations total, two fresh K/C/f assemblies, two spectrum checks, and two 99×99 solves. No run or tests occurred. Root accepted this source scope after the record fix; this is not a numerical PASS or fracture qualification.
