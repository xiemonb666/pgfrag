# Independent review — output-evidence v2

**Status:** V2's output-evidence changes pass static review for a root decision on whether to authorize one diagnostic run. No solver, tests, numerical assembly, integration, render, or Git operation was run. This review does not authorize execution.

## Freeze and delta verification

The v2 `source_freeze.sha256` sidecar matches its `source_freeze.json`; every owned-source hash and all seven read-only dependency hashes match the manifest. V2 solver SHA-256 is `5369c9ff661959c93721ff2a4986cbeffeb9419b78d9dca0a018e64b31fecd2e`, case SHA-256 is `d4246341d4e5c3c8cb38ff6caf681010393cb73bd589981eeb07cf6a55fbe3b1`, and freeze-manifest SHA-256 is `d8de9597d82e5db9a983f6be83f5c8407ec9cd56e9a8148ea22caf86d642b990`. The manifest's preserved v1 solver, case, README, review, and freeze hashes still match their earlier snapshots.

Static text comparison confirms the volume and boundary quadrature builders, boundary load assembly, volume stiffness assembly, manufactured coefficient map, spectrum diagnostics, and direct scaled saddle solver are byte-for-byte identical to the reviewed v1 function bodies. V2 adds output evidence and serialization controls without changing the frozen mechanics or numerical rules.

## Output evidence and finite-value handling

The solved-result path now writes each 54×54 K matrix with units N/m, its 54-entry physical-mean row, the 54-entry load vector in N, and the ordered canonical IDs that index all vectors and both matrix axes. It also retains the 256 off-grid coordinates, exact displacement and gradient, and each quadrature level's solved displacement and gradient. For this UZ field, `grad(u_z)` is dimensionless; the provider's basis gradients are 1/m, multiplied by displacement coefficients in m. The corrected v2 labels use dimensionless units consistently.

The extra array estimate is consistent with the saved shapes: two 54×54 matrices, three 54-vectors, and 256×(3 coordinates + 4 reference fields + 2×4 solved fields) float64 values, totaling 78,672 bytes. The result has a 4 MiB cap and an oversize candidate becomes an UNKNOWN receipt. The added arrays and vector ordering are therefore reviewable without reassembly or runtime output.

The recursive JSON conversion handles NumPy arrays and scalar types, including nested `np.bool_`. Unlike v1, it raises on nonfinite Python floats. V2 checks the retained matrices/vectors and weak defect, solution coefficients, energy/work/residual scalars, and off-grid reference/solution/error arrays. Serialization failure is converted to an UNKNOWN result. There is no render path, so no optional visualization can fail after the primary JSON save.

## Remaining receipt and diagnostic caveats

The staging-file wrapper still promotes any existing worker staging file without checking the worker return code or parsing/validating the JSON. The worker uses `Path.write_bytes` directly on that staging path. If an I/O failure or abrupt non-timeout worker termination leaves a partial file, the parent can promote malformed JSON. This is a low-probability filesystem failure path, not an observed failure; if atomic receipt integrity must cover that case, the writer/parent needs a write-to-temp-and-rename or post-write validation step.

The scalar translation compatibility identity is still not reported as a single residual. V2 exposes the raw force and mean vectors, side resultant, and Lagrange multiplier, so root can check the compatibility relation from the saved data. No numerical tolerance or new solve gate is inferred here.

The result remains explicitly `SOLVED_DIAGNOSTIC_ONLY_UNQUALIFIED`; fixed theta sensitivity is not convergence evidence, exact manufactured energy is not solve qualification, and the source does not claim full-vector equilibrium, finite crack-end accuracy, local enrichment convergence, calibration, SIF, fracture, or propagation.
