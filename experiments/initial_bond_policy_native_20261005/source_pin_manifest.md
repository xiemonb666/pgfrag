# Source pin manifest

SHA-256 values are lowercase hex. Base files are read-only inputs under the pinned vendored package. Overlay files are staged copies or additions under this experiment directory.

| Role | Exact path | SHA-256 |
|---|---|---|
| Base module | `pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/Peridynamics.jl` | `195127ad9eaaafea064c9b7087c83d7b80b3330b3a46f62b38e3838b00b1f4cb` |
| Overlay module | `overlay/src/Peridynamics.jl` | `87cd708d01b2d766a8da1bce28d19b22262a2ab1b9b001ccbdc96674b0bc30ca` |
| Base Body | `pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/discretization/body.jl` | `60762a4f52c9b6c0b8f23d295b266c15bc697b1e2a356ddd759a561236bd3421` |
| Overlay Body | `overlay/src/discretization/body.jl` | `e1b79560030e2675c40c7127d1e237edda64fcbe562e443ed5f2d21cd3ca8ed8` |
| Base bond system | `pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/discretization/bond_system.jl` | `c98bed960c6b6789f015f02e400d6b3a99e09df22c529f386c384a1564bb69aa` |
| Overlay bond system | `overlay/src/discretization/bond_system.jl` | `544de92a6611f5ba9ca906c54145cd5e720715a8fb9b8d7fa4a07e4b33cabfda` |
| Base Threads multibody handler | `pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/core/threads_multibody_data_handler.jl` | `245eb4e3d8ab28c3b89f3f97df09996bee51e44520587cde14c669c7ff165bee` |
| Overlay Threads multibody handler | `overlay/src/core/threads_multibody_data_handler.jl` | `e2bc105098eb766d0eb7683d56e0848ff39666fb1aa297c5baa91d4700e177a7` |
| Base MPI multibody handler | `pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/core/mpi_multibody_data_handler.jl` | `e32b0c3c83dab0657fcf9b6219f9b10a16b1aa52d84232533f710efb56082272` |
| Overlay MPI multibody handler | `overlay/src/core/mpi_multibody_data_handler.jl` | `1f480688a354869aaf8527d0883e44cc436c87140c28f6bdb03e75e01910f627` |
| Base MPI single-body handler | `pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/core/mpi_body_data_handler.jl` | `84b8303ddeaa18902d40904b2e540591678a81566df38c72303e0585eefb9150` |
| Overlay MPI single-body handler | `overlay/src/core/mpi_body_data_handler.jl` | `60397cd45c1da41c4829d0d9d8be0ba85ed085d5a216f43d5ce91c89d663cf83` |
| Base LICENSE | `pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/LICENSE` | `0dc42bc9160a936d5fa5ec6f622e94d2baa1e02d13a675e9080795c156da483e` |
| Overlay LICENSE | `overlay/LICENSE` | `0dc42bc9160a936d5fa5ec6f622e94d2baa1e02d13a675e9080795c156da483e` |
| Base fork Project input (not copied) | `pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/Project.toml` | `b888cd311dad495d5c1283ce2972a184f04ed42b8b180c907be2bb12e20cc986` |
| Base construction-order reference (not copied) | `pgfrag/src/pgfrag/resources/julia/vendor/Peridynamics.jl/src/discretization/body_chunk.jl` | `e076d40b5e39d2154ab341342987ff7addb4d75bd375109477acd6de7a804cf4` |
| Source-only new helper | `overlay/src/discretization/initial_bond_policy.jl` | `7057585eb58b6ef3d005d0988586b659e54a28b01911187b3cbafc05afc2ec17` |
| Source-only Julia regression file (run03 probe-shape correction) | `tests/test_initial_bond_policy.jl` | `5c40945c0a31460fad662487fa513b3a724e2b4275912906ce5e2adb79d4ef74` |
| Applicable isolated patch | `initial_bond_policy_overlay.patch` | `68748b1cebe261b76f4adc8d1056e7319ee8c69176700bdd34f8008830dcc557` |
| Isolated fork compile/test plan (not executed) | `isolated_fork_compile_test_plan.md` | `11d7fc4f2181c994a72c56a5dca7a1358f888438eb919b15467524fbb4954ed0` |
| Architecture summary including MPI body guard | `architecture_report.md` | `65f6114cb5f92a2b06f13fb479039d2534ba6454c7229021ce4b9869cba90c83` |
| Static review note with follow-up delta | `pgfrag/tmp/review_20261005/bridge_independent_review/initial_bond_policy_native_static_review.md` | `3624d6d2299090d0bf9929fc6f03a7200f9c83e6c31a62f82fe2a4765410de00` |
| TOML transport compiler source | `pgfrag/experiments/material_visibility_native_bridge_20261005/native_mask_transport.py` | `6779b3242ef87c13d5626a1e6d89414e1cdd7d6ed8f7e264fb01ed83edef7ea9` |
