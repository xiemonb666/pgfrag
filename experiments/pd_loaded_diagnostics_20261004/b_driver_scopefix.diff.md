# B driver static-review correction record

The early draft prefix `c1855a8c…` was seen by root but its exact bytes were not retained and it was never executed; it is not reconstructed. The actual pre-fix local bytes are preserved as `run_affine_critical_stretch_probe_pre_scopefix_d7c8.jl` (12,438 B, SHA256 `d7c8c11fe9d24ff3bd395240942d625bae0e17e6cbf70619b2d87aade09f03aa`). The pre-fix design manifest is preserved as `manifest_pre_scopefix_2c246c38.json` (10,170 B, SHA256 `2c246c38f5f298bc33d069307a8ff908a1f6a4157fc533bf152ab6632ec2a01b`). Both describe an unexecuted draft.

The final driver makes only the two requested corrections:

```diff
- require(volume[i] > 0 && all(isfinite, position[:, i]), ...)
+ require(isfinite(volume[i]) && volume[i] > 0 && all(isfinite, position[:, i]), ...)
```

```diff
+ directed_bond_counts = [state.directed_state_rows for state in state_summaries]
+ all_graph_counts_match = length(directed_bond_counts) == 3 &&
+     all(count -> count == Int(expected["directed_bonds"]), directed_bond_counts)
+ require(all_graph_counts_match, ...)
- "native_directed_bonds" => length(system.bonds),
+ "native_directed_bonds" => state_summaries[1].directed_state_rows,
```

This avoids reading `system` after the loop's local scope and asserts all three fresh native handlers contain the expected 6,124 directed bonds before reporting. No Julia/B driver execution or solver run was performed; the final manifest hashes the corrected draft.
