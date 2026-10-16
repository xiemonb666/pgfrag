# PG-Frag handover — 2026-10-16

## Scope and pause state

Development is paused. This handover and the package README are documentation
work only. No solver run was started. During this handover, the main agent ran
the offline test set recorded below; resumption and technical acceptance
remain with the main agent.

Evidence paths below refer to the current local workspace unless marked as
tracked documentation. At the time of this handover, `data/` and
`experiments/continuation_20261016/` are untracked in the PG-Frag Git checkout.
The in-progress `src/pgfrag/mother_binding.py` and
`tests/test_mother_binding.py` are also untracked. A documentation-only GitHub
change will not publish those payloads or source files. Keep the paths as
provenance pointers; curate and add evidence or code separately if it is to be
part of a public release.

## Installation and runtime map

The install command follows the tracked `pyproject.toml`. From the PG-Frag
repository directory on Windows PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e .
```

The package metadata declares Python `>=3.11` and only one runtime dependency,
`numpy>=1.24`. The source package includes the Julia project resources and
headless Peridynamics.jl source. The locked `Manifest.toml` records Julia
1.10.12; `src/pgfrag/backend.py` accepts Julia 1.10.x. Install the Julia
executable separately. For a fresh PG-Frag Julia depot, set the paths and
instantiate the locked environment once from the repository directory:

```powershell
$env:PGFRAG_JULIA_EXE = 'C:\tools\julia-1.10.12\bin\julia.exe'
$env:PGFRAG_JULIA_DEPOT_PATH = 'C:\path\to\pgfrag-julia-depot'
& $env:PGFRAG_JULIA_EXE --project=. -e 'using Pkg; Pkg.instantiate()'
```

`PGFRAG_RUNTIME_ROOT` may instead name a Julia project directory containing
`Project.toml`, `Manifest.toml`, `run_fixture.jl`,
`backend_provenance.json`, and `vendor/`. The backend searches Julia on `PATH`
and the workspace extraction path before returning an unsupported-runtime
result. It selects an isolated PG-Frag depot rather than Julia's implicit
global depot.

Akantu 5.0.7.post1 is a separate, optional native-cohesive research runtime.
The recorded experiments use WSL Ubuntu, Python 3.12.15, NumPy 2.5.3, and the
local Linux wheel at `../third_party/akantu-runtime/downloads/`; their runner
invocations use `../third_party/akantu-runtime/venv/bin/python`. This runtime
is outside the PG-Frag repository and is not installed by the package command
above. The experiment README files are the source for exact WSL invocations.
The Python SDK, Julia bar fixture, and Akantu cohesive experiments are
separate capability layers; no single pip command installs all three.

Command provenance: `pyproject.toml` declares Python and NumPy requirements;
`src/pgfrag/backend.py` defines Julia discovery, version checking, and depot
selection; `src/pgfrag/resources/julia/Manifest.toml` pins the Julia package
environment; `docs/headless_backend.md` records the Peridynamics.jl source
and patch provenance. The Akantu runtime version and hash are recorded in
`data/continuation_20261016/extrinsic_capture_preflight/REPORT.md` and
`data/continuation_20261016/reference_v6_preflight/run_control.json`.

## Verified progress and limits

| Workstream | Evidence-backed state | Limit that remains |
| --- | --- | --- |
| Python SDK/CLI | Package metadata, job/asset validation, provenance checks, and a bounded homogeneous-bar adapter exist. The CLI is a narrow prototype. | `validate-asset` is a contract check. Surface-to-volume discretization is unsupported by the CLI, and a valid asset does not imply a solver run or material pass. |
| Prescribed weak-seam coupon | `data/continuation_20261015/cohesive_weak_seam_v4/result.json` records a completed 3,244-step native Akantu evolution. `data/continuation_20261015/native_cohesive_fragments_v5/manifest.json` records two components and 158 paired failed facets; its separation record reports the saved post-debond states separated. `data/continuation_20261015/native_cohesive_fragments_v6/control_comparison.json` reports a narrow transient-compliance control pass. | The seam was preinserted and its material was uncalibrated. Mesh sensitivity, continuous-time separation, full physical qualification, typed R2, and PRETRAIN eligibility are not established. A narrow control pass does not change that status. |
| Reference-vessel zero-load control | `experiments/continuation_20261016/reference_vessel_benchmark/native/run_attempt_01/completion.json` records 46 completed steps. `data/continuation_20261016/reference_vessel_benchmark/native_evolution_audit.json` reports `ZERO_LOAD_NATIVE_EVOLUTION_NUMERICAL_PASS`; zero-load states, unchanged unfailed QPs, and the independent force/momentum checks passed the stated tolerances. | This is a cost and zero-load control on the 56,139-T4 v5 mesh with placeholder values. It is not fracture, material, or physical qualification. |
| Reference-vessel v6 geometry | `data/continuation_20261016/reference_v6_preflight/README.md` records the geometry candidate: 86,624 T4, 16 source face IDs, and 0.240308% meshed CAD-volume error. The MSH alias parser-equivalence record passes its byte/tag/connectivity checks. | Geometry evidence only. It uses a candidate profile, local foot-fillet area errors remain about −9.7%, material and interface values are placeholders, and a prescribed weak seam is present. |
| Reference-vessel v6 native preflight | The current `run_control.json` is `READY_FOR_SINGLE_NATIVE_INITIALIZATION_PREFLIGHT`; `preflight_budget.json` says `WHOLE_OUTPUT_FORECAST_PASS_RUN_NOT_STARTED` and `native_solver_execution_started: false`. `preparation_failure_01.json` and `preparation_failure_02.json` report failures before native execution. There is no native `started.json`, result, exit-code record, or initialization audit in this attempt. | **The v6 native preflight was prepared but has not been executed.** The two preparation failures are not solver failures or passes. The one-shot initialization attempt remains unused in the recorded control. |
| Extrinsic onset / signed work | `data/continuation_20261016/extrinsic_event_audit/REPORT.md` reconstructs the historical 96-T4 step-652 negative dissipation from signed QP work; it does not clamp or take magnitudes. `data/continuation_20261016/extrinsic_capture_preflight/REPORT.md` records completed zero-step getter and native-shape probes: `getShapes()` returned 696 rows (232 Tri3 facets × 3 QPs) for the frozen coupon's affine shape. | The signed-work accounting explains the negative getter. It leaves the initial directional-stiffness/onset behavior unresolved; the report describes a possible near-zero-opening mechanism, not a sole established cause. The `getEnergy("kinetic")` getter may rebuild the mass matrix according to cached v5.0.7 source; whether the installed 5.0.7.post1 wheel has the same side effect is unknown. Root/Astra's current decision is to preserve the original C call order and capture state before and after getters, including mass-matrix state. |
| Extrinsic h/h2/h4 refinement | `data/continuation_20261016/extrinsic_onset_refinement/shape_only_preflight/preflight_status_final.json` records a completed zero-step shape-only preflight. `prepare_validation.json` says `PREPARED_NO_MODEL_NO_SOLVER_STEPS` and `campaign_launched: false`. | **The h/h2/h4 native campaign was only prepared; it has not started and is not root-approved.** The shape probe supports only the frozen coupon's affine Tri3 layout. Keep it behind main-agent review of the frozen protocol and current Akantu wheel/source provenance. |
| Facet ancestry | The focused 25-case check was reported PASS for its declared contract. `experiments/continuation_20261016/facet_ancestry/README.md` states the native assumption; `saved_96t4_run_receipt.json` holds the historical C/D connectivity mapping. | Mapping depends on the native backend preserving initial T4 row and local-corner order. This is connectivity ancestry, not physical fracture, material state, or stable quadrature-point order. |
| Mother binding | `src/pgfrag/mother_binding.py` contains a geometry-only source/mesh binding implementation; `tests/test_mother_binding.py` defines nine focused cases for reordering, role/FID conflicts, hash binding, claim rejection, and bounded output. Those tests were included in the main agent's offline 55-test run below. | This remains work in progress, not an accepted integration result. The source and tests are currently untracked. The main agent has not completed full source review or actual run05 binding acceptance. The unit-test result does not establish an integrated mother binding or a solver/physical pass. |

## Offline verification reported by the main agent

The main agent ran the following tests with `PYTHONPATH=src` and the A12
virtual environment's Python; no native solver was started:

```text
pytest tests/test_mother_binding.py tests/test_facet_ancestry.py tests/test_cohesive_mesh.py tests/test_cohesive_stability.py tests/test_cohesive_fragments.py tests/test_cohesive_vessel.py -q
55 passed in 1.07s
```

This is offline test evidence for the listed unit and contract checks. The nine
mother-binding tests are included in the 55; their pass does not replace source
review or actual run05 binding validation. No saved test receipt path was
provided with this result.

## Remaining qualification work

The current evidence does not establish a calibrated material model or a
freely stress-driven fracture of the reference vessel. The demonstrated
fragmentation uses an imposed weak plane. Real-profile fidelity, native
onset/directional response, independent mesh and time-step behavior, global
energy closure with verified loading phases, and contact/separation behavior
still need evidence. No PG-Frag result is currently qualified as a physical
R2 label or RelicWeave/SynLab training data.

## Sibling-repository collaboration

- **VesselForge** supplies source geometry, lineage, and meshing artifacts.
  PG-Frag consumes frozen files and hashes; it does not import VesselForge as
  a runtime dependency.
- **PG-Frag** records numerical states, fracture events, fragment topology,
  and source-facet ancestry as separate, reviewable evidence. Geometry and
  native-run status remain distinct qualification stages.
- **RelicWeave/SynLab** owns downstream fragment observation, retrieval, and
  training research. PG-Frag can supply candidate evidence only after the
  geometry, numeric, material, fracture-surface, and data gates pass. Solver
  output, a closed mesh, or lineage code cannot make training eligibility
  automatic.

The tracked `docs/research/INNOVATION_LEDGER.md` describes the research
boundary and the required downstream retrieval evidence. Licensing and
redistribution review for external runtimes, source assets, and generated
data must be completed before promising a public data bundle.

## Resume order

1. Re-read this handover and inspect the current workspace records before
   doing technical work. Keep the development pause in force until the main
   agent explicitly resumes it.
2. Review `reference_v6_preflight/run_control.json` and both preparation
   failure receipts. If the prepared one-shot remains in scope after resume,
   confirm the frozen inputs and budget first; record any single native
   initialization attempt and its outcome. Do not reinterpret preparation
   records as a native run.
3. Have the main agent review the extrinsic capture interpretation and the
   frozen refinement protocol. Do not launch h/h2/h4 based on
   `prepare_validation.json` or the zero-step shape probe alone.
4. Review `mother_binding.py` and its nine focused test cases, then validate
   actual run05 binding against its frozen inputs. The main agent owns source
   review and acceptance; preserve any failure evidence. The offline unit-test
   pass above is not integration acceptance.
5. Resume physical work only with an explicit material/geometry qualification
   plan: real source-profile checks, calibrated parameters, a freely
   stress-driven case, verified energy/loading phases, and independent mesh
   and time-step comparisons. Keep downstream RelicWeave/SynLab data use
   closed until those gates are evidenced.
