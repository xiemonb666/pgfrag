# Isolated import/precompile validation plan

Status: plan only; do not execute as part of this candidate task. The parent agent may decide later whether to perform exactly one import/precompile run after reviewing the candidate.

## Isolation layout

Create a new scratch directory under this review directory, for example isolation_20261004/, by copying the complete pinned environment from pgfrag/src/pgfrag/resources/julia/. This copy contains the current pinned Project.toml and Manifest.toml plus vendor/Peridynamics.jl. Make no edits to the original resource environment or package tree.

In the copied package only:
1. Replace src/auxiliary/precompile_workload.jl with the candidate source in pgfrag/experiments/headless_precompile_workload_fix_20261004/precompile_workload.jl.
2. Change the copied Peridynamics Project.toml version from 0.5.4+pgfrag.1 to the explicit candidate marker 0.5.4+pgfrag.1.headlessprecompilefix.20261004.
3. Change only the copied [[deps.Peridynamics]] Manifest version to that same marker. Keep the UUID, local path vendor/Peridynamics.jl, dependency list, and every other pinned package entry unchanged.

Capture SHA-256 for the source copy, candidate workload, copied package Project, and copied environment Project/Manifest before the Julia invocation. Verify the baseline environment Project and Manifest hashes still match candidate_freeze_manifest.json.

## Preventing an old compiled-cache hit

Use a new empty depot below the isolation directory and set JULIA_DEPOT_PATH to that single path for the one Julia process. Do not append the existing user depot: Julia searches the depot stack for compiled package images, so including the old depot could expose its Peridynamics cache. The scratch depot's compiled directory must be absent before launch.

Materialize dependencies in that isolated depot from the copied pinned Manifest (for example, instantiate the copied environment with automatic precompilation disabled inside the same process before importing Peridynamics). The copied local vendor path resolves the candidate source. Retaining the package UUID keeps the project dependency valid; the distinct version marker and an isolated depot with no previous compiled cache make the candidate identity and cache provenance explicit. Record the exact Julia executable/version and the effective depot path in the run receipt.

The Julia manual documents that import can create/use package cache files and that they are stored under the compiled subfolder of the first depot; the Pkg depot documentation describes the depot path as a search stack. See https://docs.julialang.org/en/v1/manual/modules/ and https://pkgdocs.julialang.org/v1/depots/.

## Single planned invocation and preserved evidence

Run exactly one Julia process with the isolated project, startup/history disabled, isolated depot, and one using Peridynamics. If dependency materialization is needed, run Pkg.instantiate(; allow_autoprecomp=false) in that same process before the import so it does not create a separate precompile pass. Print the loaded Peridynamics version and pathof(Peridynamics) after import to bind the result to the candidate copy. Package precompilation will execute the package's existing tiny workload as designed; do not separately call submit or run another solver invocation.

Capture stdout, stderr, and process exit code to distinct raw files without filtering. Preserve those files even on failure, plus the exact invocation, copied source/package/environment hashes, and resulting compiled-cache file listing. A stderr precompilation errored message or nonzero exit remains a failure; do not suppress it, edit the log, or translate a caught error into success.

This one run checks only whether the headless package import/precompile workload completes with the candidate. It is not evidence of nonzero mechanics, constitutive response, crack propagation, or real fracture, and it does not authorize a production backend update.