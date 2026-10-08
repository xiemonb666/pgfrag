# Isolated fork compile/test execution plan

Status: plan only. No Julia process, package manager, source patch application, runtime copy, or test has been run by this task. Root authorization is required before execution.

## Inputs and isolation

1. Use the exact base package and source files recorded in `source_pin_manifest.md`. Before copying, verify the source file hashes against the manifest, including the original `Project.toml`, `Manifest.toml` (if present), and each changed Julia source file. Copy the complete pinned headless vendor package into a disposable task-owned path such as `pgfrag/tmp/initial_bond_policy_native_20261005/fork/Peridynamics.jl`; verify the copy's source pins before continuing.
2. Apply `initial_bond_policy_overlay.patch` only to that disposable copy. The patch may change only the six Julia source files and add the policy helper listed in the manifest, plus add the `SHA` and `TOML` stdlib UUIDs to the copied `Project.toml`. It must not alter or regenerate the copied `Manifest.toml`. Do not write to the vendored source, production tree, repository `Project.toml`/`Manifest.toml`, Git state, or shared Julia depot.
3. Use the existing Julia depot as a read-only secondary depot and a fresh task-owned scratch depot first in `JULIA_DEPOT_PATH`. Set package resolution offline. Do not run `Pkg.instantiate`, `Pkg.resolve`, `Pkg.develop`, `Pkg.test`, or any installation step. If required dependencies are absent, stop and report the blocker without fetching or modifying dependencies.

## Bounded source test

After root authorizes execution, launch the Julia binary with a hard 180-second process timeout, `--project=<disposable fork>`, `--startup-file=no`, and `--history-file=no`. The only test payload is `tests/test_initial_bond_policy.jl` from this experiment directory. Load `Peridynamics`, assert `pathof(Peridynamics)` resolves to the disposable fork's `src/Peridynamics.jl`, then include that test source. This compiles the package source and exercises the focused pure helper/guard tests without constructing a solver or native body.

Treat a timeout, unexpected package path, precompile/import failure, missing offline dependency, or source test failure as a stop condition. Preserve the disposable copy and its logs for root review; do not broaden into the package suite, run a solver, submit a production job, run CAD/classifier code, or claim native graph readback. Julia cache writes, if any, must remain in the task-owned scratch depot.
