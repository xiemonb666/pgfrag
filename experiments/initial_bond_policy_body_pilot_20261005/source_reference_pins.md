# Source and artifact pins

SHA values are lowercase SHA-256. New pilot files are source-only. The existing fork overlay, approved patch, root project/manifest, and runtime package source are read-only references.

## New pilot source and fixture files

| File | SHA-256 |
|---|---|
| `README.md` | `7f4468532784369efe981854f38f30790f8742f288a830929a75b5e90153e6ca` |
| `planned_assertions.md` | `f1a3e3b96213d05f3b318e48f96f0eed55c6e6cc7b21cb631703de4ef45c9d2f` |
| `fixedcase_points.csv` | `0f635481d102e98c4ace79dc880e3cdfc5cbe80fe5a1ee73eac463c113c5fab7` |
| `fixedcase_components.csv` | `1f55f30c2803618fa689b541c79757834d7dd86945f5c614a4c4b3802d541940` |
| `fixedcase_source_marker.txt` | `2091ec6772be7411e42efa37038017460050c26cbd977fc1cab9edbdcbb830fe` |
| `generate_fixture_policy.py` | `570dba4a4360b1d22b9aea336c4ee4539364c8c7acb9282900c5258bb1488a5d` |
| `native_body_pilot.jl` | `038f289af08b91f5d4cbee599aaab0e1e37061d5f1615316c62d15343faae4bf` |
| `../../tmp/review_20261005/initial_bond_policy_body_pilot/static_review.md` | `37e0f6577f7baa24fc1558948a13b93a8631648928ea92eb1cad1fb5e08dd016` |

## Existing isolated bridge/fork input pins (read-only)

| Source | SHA-256 |
|---|---|
| `initial_bond_policy_native_20261005/source_pin_manifest.md` (23-entry source ledger) | `fd9f1278577ee8266d1a9767c6f4cedcc17dff716c2debbb69ed7921c844466f` |
| Approved isolated overlay patch | `68748b1cebe261b76f4adc8d1056e7319ee8c69176700bdd34f8008830dcc557` |
| Overlay module | `87cd708d01b2d766a8da1bce28d19b22262a2ab1b9b001ccbdc96674b0bc30ca` |
| Overlay `Body` | `e1b79560030e2675c40c7127d1e237edda64fcbe562e443ed5f2d21cd3ca8ed8` |
| Overlay `BondSystem` / filter lifecycle | `544de92a6611f5ba9ca906c54145cd5e720715a8fb9b8d7fa4a07e4b33cabfda` |
| Overlay strict initial policy loader/scope/filter | `7057585eb58b6ef3d005d0988586b659e54a28b01911187b3cbafc05afc2ec17` |
| Original root `pgfrag/Project.toml` | `eebff26de1f99a56b29a43ad767d2877bba5b6a0e4a03de1f909eab00f59f4ae` |
| Original root `pgfrag/Manifest.toml` | `eb7d8902b84cfc64faaa0d45b0e6a3ae8534a252d06399f3db074da5a738e594` |
| Frozen bridge TOML compiler | `6779b3242ef87c13d5626a1e6d89414e1cdd7d6ed8f7e264fb01ed83edef7ea9` |

## Horizon, graph, material, and correction source pins

`PointNeighbors` is version `0.6.7`, manifest git-tree SHA-1 `e1273088c75fe6f79e34dfdd6cb67a1cd565619b`, installed at the pinned depot package tree `STFG7`. `Peridynamics` is the isolated `0.5.4+pgfrag.1` source path from the root manifest, with its approved overlay supplied by the run03 disposable project.

| Read-only source | SHA-256 | Reviewed location |
|---|---|---|
| PointNeighbors squared-radius filter | `9a516dff1675cca5e12fe5f6ce0274239ab689a10d22bef021e400bfc42cf24d` | `third_party/pgfrag-runtime/depot/packages/PointNeighbors/STFG7/src/nhs_grid.jl:549` (`distance2 <= search_radius^2`) |
| PointNeighbors API prose and neighbor call | `040525485b7b4d4eceaf57f7e55f8e14c67ced570668989843c3baa4616b813f` | `.../PointNeighbors/STFG7/src/neighborhood_search.jl:23-25, 237-263` |
| Native OSB weighted-volume formula | `65ffd614270c8f409a28f83e5def1ea1abc5f94ddcb29eac7b2cd7320e720846` | `.../Peridynamics/9JARS/src/physics/ordinary_state_based.jl:207-219` |
| One-chunk decomposition | `d766e8c582355a00725abd1ed17343fd0b7a2fac2de9a0947920bd8761031e22` | `.../Peridynamics/9JARS/src/discretization/decomposition.jl:15-63` |
| BodyChunk construction and initialization seam | `e93ec9f1c56e880bdeac9dccd6c7b673d4ab26d7e649a727678e004d56309da4` | `.../Peridynamics/9JARS/src/discretization/body_chunk.jl:41-51` |
| Material horizon parameter assignment | `4a032a539214e3d1f1291d6e63b19ff6254533f849950c10d5a8e4b98fde96ab` | `.../Peridynamics/9JARS/src/physics/material_parameters.jl:64-92, 100-106` |
| Default correction factor | `63c3a3c458f84dcaa0ad3345c141f67fad2ae42f165b1df3c32ecc6e9c508186` | `.../Peridynamics/9JARS/src/discretization/bond_system_corrections.jl:17-23` (`NoCorrection` gives `1`) |
| Default linear kernel | `37a68ed8ecdcc3ad77ccff08ba35a43201b7d99df6002ee732fc49a1d77261fb` | `.../Peridynamics/9JARS/src/discretization/kernels.jl:13-22` |
| Chunk point/localizer ordering | `ee37278c8966dedc3248d561e75a5897584eec9695827c521d5c7a26dc12f28d` | `.../Peridynamics/9JARS/src/discretization/chunk_handler.jl:19-37, 74-80` |
| Threads handler path not constructed by this pilot | `f341855e3eac7fcf4ac7ffa139150c9f9c171c49ba6ad1bd146285422355d5f5` | `.../Peridynamics/9JARS/src/core/threads_body_data_handler.jl:31-59` |

## Execution state

No Julia, Python fixture compiler, pytest, CAD, solver, or test command was executed for this pilot. The root's previously accepted run03 result (27 focused assertions) concerns the isolated helper/runtime probes only; it is not native graph readback evidence. The new pilot remains `NOT_RUN` until separately accepted and run in the disposable fork.
