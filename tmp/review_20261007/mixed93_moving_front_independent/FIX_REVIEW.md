# Narrow review of the moving-front source fix

**Disposition:** the two reported source issues are corrected in this revision. No new blocker was found in the reviewed delta. This is not numerical authorization; all mechanical outcomes remain unknown.

Current pins verified: entry `a124f68a4e605a244123c30f9256838c343816531d7cda78445d9a755f39442f` (49,479 bytes/617 lines), case `1c967f9b7a7ed8d786e93ecb4c90e7ef57e8e8a10eaefa80ff4c49265e7b6662`, and `SOURCE_PINS.json` `604d484b534ab16d7b0c151816b04ea35f2a17e35d0dac436897193465dac311`. Entry/case hashes match the pin. All 9 source references and 14 case dependencies match their current pinned SHA-256 values. Execution authorization is still false.

Batching now iterates over the coarse and fine ranges separately. For 448 coarse rows, chunks are `[256,192]`; for 896 fine rows, chunks are `[256,256,256,128]`. The concatenated query order remains coarse then fine, giving 1,344 unique rows and the exact asserted six-batch list without crossing the 448-row boundary. The status/support/condition checks and max-256 limit remain intact.

The energy metric keeps its existing field key, but now explicitly defines `external_work_J` as the fixed-dead-load `q·f` force-displacement pairing, not path-integrated work. The case has matching wording, and the potential remains `Π=U−qᵀf`; `G_fd` remains the finite difference divided by `LΔa`. These additive metadata fields remain compatible with the v1 result schema and the unchanged source accesses.

The narrow diff introduces no incompatible case, MLS, analytic-field, semantic-registry, or elastic-assembler access. The output schema identifier still matches the case contract. The earlier report’s nonblocking note about relying on the pinned cache metadata contract for cached MLS/source-ID provenance remains applicable.

Before any numerical attempt, root must make a separate authorization decision and update the case/pin authorization consistently. Until a run is separately authorized, MLS status, wrench balance, volume measure, raw spectrum, KKT residuals, and quadrature sensitivity remain **UNKNOWN**. No NPZ payload was read; no numerical modules, provider, field, assembly, solver, or tests were run.
