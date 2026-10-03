# PG-Frag Core prototype

This is an independent Python SDK/CLI around the pinned Peridynamics.jl 0.5.4
CPU solver. It does not import VesselForge or RelicWeave/SynLab.

The first executable slice is intentionally a bounded homogeneous bar fixture.
It supports Peridynamics.jl `BBMaterial` (`nu=0.25` in this 3D implementation)
and native `OSBMaterial` for valid `0 < nu < 0.5`. All values must be marked
`UNVALIDATED_PROTOTYPE`; running the solver is not material calibration or a
physical ceramic fracture pass. The bar reports material-point state and
endpoint intact/broken horizon bonds. It does not emit a geometric crack
surface, dynamic fragment shapes, contact, historical A12 `R2`, or a per-step
fracture event chronology.

The `validate-asset` command checks the current `forge.asset.v1` package
contract, hashes, required arrays, material binding and geometry-report
status. That is a contract check only. `run` explicitly returns `UNSUPPORTED`
for a surface asset because mesh-to-volume discretization and void/topology
classification are not implemented. Code `255` remains an asset region code;
PG-Frag never converts it into a dynamic fracture label.

## Run

Install the Python package from a wheel or source tree. Provide Julia 1.10.x with
`PGFRAG_JULIA_EXE`, or make it available as `julia` on `PATH`. A checkout may
use the local runtime under `third_party/pgfrag-runtime`.

```powershell
$env:PYTHONPATH = "C:\path\to\pgfrag\src"
$env:PGFRAG_JULIA_EXE = "C:\path\to\julia-1.10.12\bin\julia.exe"
$env:PGFRAG_JULIA_DEPOT_PATH = "C:\path\to\isolated\pgfrag-depot"
python -m pgfrag.cli validate-job C:\path\to\job.json
python -m pgfrag.cli run C:\path\to\job.json
```

`PGFRAG_RUNTIME_ROOT` may point to a project directory containing
`Project.toml`, `Manifest.toml`, `run_fixture.jl`, `backend_provenance.json`,
and `vendor/`; the installed wheel carries these files. `JULIA_DEPOT_PATH` is
always set explicitly to the configured isolated depot or a PG-Frag-specific
user cache. `local_smoke` limits a run to 1000 points, 500 steps, one Julia
thread/chunk, and 120 seconds. Multi-chunk export is not implemented. Ctrl+C requests child-process cancellation and the
CLI records the resulting status and logs.

`forge.asset.v1` files and result units are versioned independently. Surface
asset vertices enter in mm; the solver works in SI. Numerical bar result
coordinates and point volumes are serialized in m and m³. Unsupported result
fields are marked `UNKNOWN` or `UNSUPPORTED` with a reason; they are never
filled with synthetic zero values. Fragment-count targets are forbidden.

The vendored headless backend has its own identifier,
`peridynamics_jl_headless_0_5_4`. See
[`docs/headless_backend.md`](docs/headless_backend.md) for the source pin,
license, I/O-only patch and mechanics-source hash evidence.
