# PG-Frag

PG-Frag is a research SDK and evidence workflow for peridynamic and cohesive
fracture experiments. It is separate from VesselForge and RelicWeave/SynLab:
it consumes versioned geometry and source records, and can produce candidate
fracture evidence, but a solver run or a closed mesh does not establish
material calibration, physical validity, or training eligibility.

**Project status (2026-10-16): development is paused.** The next-step order,
verified experiment state, and evidence paths are in
[`docs/handover/20261016/DEVELOPMENT_STATUS.md`](docs/handover/20261016/DEVELOPMENT_STATUS.md).

## What is in this repository

The installable Python package is `pgfrag-core` 0.1.0. Its base dependencies
are Python 3.11 or newer and NumPy 1.24 or newer. The `pgfrag` command exposes
asset and job contract checks plus a bounded CPU solver adapter. The adapter's
supported first fixture is a homogeneous bar using the pinned, headless
Peridynamics.jl 0.5.4 backend. Surface-asset mesh-to-volume conversion and
production vessel fracture are outside this CLI capability; unsupported
requests return a structured status.

Research code in `experiments/` includes a separate native cohesive-vessel
prototype. These scripts are not installed as the `pgfrag` CLI, and use of an
experimental solver path does not promote a result to a physical pass.

## Install the Python package

From the `pgfrag` repository directory, create a virtual environment and
install the local package in editable mode:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e .
pgfrag --help
```

These commands follow `pyproject.toml`: the package requires `Python >=3.11`
and declares `numpy>=1.24`. To build the wheel instead, install the package
from the resulting wheel with `python -m pip install <wheel-path>`.

Example contract checks:

```powershell
pgfrag validate-asset C:\path\to\forge-asset
pgfrag validate-job C:\path\to\job.json
pgfrag run C:\path\to\job.json
```

`validate-asset` and `validate-job` perform contract checks. `run` can start a
Julia process only for a supported job/backend combination; a valid asset is
not by itself evidence that its geometry can be solved.

## Optional solver runtimes

The Python package does not bundle or install the Julia executable. The
Peridynamics.jl project, manifest, runner, provenance record, and vendored
solver source are included with the source checkout and packaged resources.
The manifest was generated with Julia 1.10.12; the launcher accepts Julia
1.10.x. Install Julia 1.10.12 separately, then set its executable and a
PG-Frag-specific depot. From a source checkout, instantiate the locked Julia
environment once in that depot:

```powershell
$env:PGFRAG_JULIA_EXE = 'C:\tools\julia-1.10.12\bin\julia.exe'
$env:PGFRAG_JULIA_DEPOT_PATH = 'C:\path\to\pgfrag-julia-depot'
& $env:PGFRAG_JULIA_EXE --project=. -e 'using Pkg; Pkg.instantiate()'
```

Run the instantiation command from the `pgfrag` repository directory. If
`PGFRAG_JULIA_EXE` is unset, the launcher looks for Julia on `PATH` and in the
workspace's `third_party/pgfrag-runtime/extract/julia-1.10.12` location. A
custom `PGFRAG_RUNTIME_ROOT` may point at a directory containing
`Project.toml`, `Manifest.toml`, `run_fixture.jl`,
`backend_provenance.json`, and `vendor/`. The backend uses an isolated depot;
it does not silently reuse Julia's global default depot.

The separate native cohesive experiments use Akantu 5.0.7.post1 under WSL
Ubuntu with the workspace's Python 3.12.15 runtime and NumPy 2.5.3. This is an
optional research layer: it is not a dependency in `pyproject.toml`, is not
installed by `pip install -e .`, and is not needed for the Python contract
commands or Julia bar fixture. The prepared workspace runtime and wheel live
outside this repository under `../third_party/akantu-runtime`; experiment
READMEs contain their recorded WSL invocation. A portable fresh Akantu setup
is not currently part of the PG-Frag package installation workflow.

The vendored Peridynamics.jl build has a narrow headless I/O patch. Its
mechanics-source and patch provenance are described in
[`docs/headless_backend.md`](docs/headless_backend.md).

## Current research boundary

The repository has evidence of native evolution on a prescribed weak-seam
coupon and export of two closed native components. That case uses an
uncalibrated, preinserted seam; it is not a calibrated ceramic fracture result.
A separate 46-step, zero-load reference-vessel run passed its narrow numerical
control and cost audit. A newer v6 mesh has geometry evidence only, and its
one-shot native initialization preflight has not run. The planned h/h2/h4
extrinsic-onset campaign is prepared but has not launched. These boundaries
and the exact workspace evidence paths are recorded in the handover document.

PG-Frag is one stage in a sibling-repository workflow. VesselForge supplies
source geometry and mesh artifacts. PG-Frag records solver states, events,
and lineage candidates. RelicWeave/SynLab owns downstream fragment retrieval
and training research. Artifacts cross these boundaries through explicit
hash-bound files; the repositories do not import each other as runtime
dependencies, and PG-Frag does not grant downstream training eligibility.

## Pause handover and remaining work

The [cross-repository handover](https://github.com/xiemonb666/relicweave-synlab/blob/main/docs/handover/20261016/START_HERE.md)
contains installation boundaries, evidence, source binding, priorities and the updated V1.48 plan.
The pause delivery passed 55 focused PG checks (VF 66 and A12 90, 211 total); these are offline
software checks, not a full-suite result or a physical qualification. Complete qualified physical
packages remain **0**. Resume priorities are root review of the v6 mother binding, the frozen
native initialization preflight, signed-work/onset diagnosis, material calibration, then validated
fracture and downstream observations. Prepared solver campaigns remain unstarted while paused.

Clone the three repositories as siblings. A12's `workspace_support/` archive preserves shared
workspace scripts, small evidence and dependency setup metadata; its restore utility checks
hashes and conflicts before copying. Large raw histories, solver executables and runtime caches
are not uploaded; the handover inventory records what remains local. Peridynamics.jl and other
upstream resources retain their own licenses and provenance. This repository currently has no
project-level `LICENSE` file; upstream licenses do not grant rights to this project's own code.
