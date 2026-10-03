# Headless Peridynamics.jl source pin

The adapter backend identifier is `peridynamics_jl_headless_0_5_4`; it is not
the unmodified upstream package. The upstream source is Peridynamics.jl tag
`v0.5.4` from [kaipartmann/Peridynamics.jl](https://github.com/kaipartmann/Peridynamics.jl).
The archived source tarball is
`third_party/pgfrag-runtime/Peridynamics.jl-v0.5.4-source.tar.gz`, SHA-256
`1b089cd3ab945b34f31458c0533f2d9d3b90eb69fd51b5d546746c55f30c68df`.
The upstream LICENSE is retained at
`vendor/Peridynamics.jl/LICENSE`, SHA-256
`0dc42bc9160a936d5fa5ec6f622e94d2baa1e02d13a675e9080795c156da483e`.

The checked-in patch is
[`vendor/patches/peridynamics-0.5.4-headless.patch`](../vendor/patches/peridynamics-0.5.4-headless.patch).
It disables VTK/LightXML reader and writer paths and removes VTK-only cell
metadata. The output API is intentionally unavailable; the adapter writes its
own CSV/JSON artifacts. A replay check copied the pristine archived `src/`
tree to an isolated temporary directory, ran `git apply --check` and `git
apply`, then compared all 63 Julia source files against the vendored source;
normalized source text matched in all 63.

All 17 `src/physics/*.jl` and all three `src/time_solvers/*.jl` files match
their recorded upstream file hashes exactly (20/20). The five changed files
are `src/Peridynamics.jl`, `src/auxiliary/io.jl`, `src/core/job.jl`,
`src/core/threads_body_data_handler.jl`, and `src/discretization/body_chunk.jl`.
They only remove/guard VTK output, the VTK reader export/include, and
VTK-specific cell objects. `backend_provenance.json` contains upstream,
fork, patch and license hashes; the Python launcher verifies these before
starting a solver process.

`Project.toml` and `Manifest.toml` resolve the package by the relative path
`vendor/Peridynamics.jl`. No package depot path is embedded in those files.
The locked project uses Julia 1.10.x and the launcher checks the executable's
reported version. It sets `JULIA_DEPOT_PATH` to `PGFRAG_JULIA_DEPOT_PATH`, the
checked-out isolated runtime depot, or a PG-Frag-specific cache; it never
inherits Julia's default global depot implicitly.

This verifies only source identity, import/run viability and I/O patch replay.
It is not upstream regression-suite parity, solver verification, material
calibration, or proof that unsupported VTK APIs work.
