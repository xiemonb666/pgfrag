# Mixed93 opening-history saved visual

This is a source-only renderer for the planned JSON result `pgfrag/tmp/review_20261007/mixed93_opening_control_20261007/run01/result.json`. Its SHA-256 remains `null`; the path was absent at initial preflight and became present during concurrent preparation, but its contents and hash were not inspected. Execution and root GO remain disabled. No opening-history result was read, no numeric or plotting package was imported, and no figure was rendered.

The renderer follows the frozen producer contract in `run_mixed93_opening_control.py` SHA-256 `ba00de0d635d2daecfd113ab28c9518d196655c89406a09472fbcf0af5fcd291`. It requires schema `mixed93_opening_control_result.v1`, status `DIAGNOSTIC_ONLY`, the fixed-face toy classification, matching producer source hashes, the exact 12-step Qstar path, and sixteen stable saved QPs at every step. It plots only saved JSON fields; it does not solve, query a provider, open an NPZ, or reconstruct a cohesive response.

Panel (a) selects three real saved quadrature points by minimum, upper-median rank, and maximum analytic opening factor `alpha`. Each point shows all 12 saved traction-opening states, connected in history order; marker color uses the shared zero-based step scale. Ties are resolved by saved quadrature ID. Panel (b) keeps the cohesive generalized internal-force component `gamma[82]` separate from the total positive-opening control force, both against prescribed `Q` in metres. Their signs and meanings remain distinct.

Panel (c) shows only cohesive stored energy, cohesive dissipated energy, and cumulative independent cohesive-law path work in microjoules. The cumulative curve is the sum of the saved `cohesive_work_increment_independent_law_j` entries; bulk work and total control work are excluded. Panel (d) maps the sixteen saved `x/z` QPs at zero-based step 3, with a fixed damage range `[0,1]` and equal physical aspect. The footer identifies full damage at steps 8, 10, and 11, including step 10 at `Q=0`.

The figure header states `FIXED_FACE`, `TOY`, `NOT A FRACTURE OR FRAGMENTATION MODEL`, `NO_GROWTH`, and `NO_MATERIAL_QUALIFICATION`. This is an opening-control toy history display, not physical fracture evidence.

`run_saved_opening_history_visual.py` uses only the standard library until an explicitly authorized child render starts. Its default source check does not read a pending result or import Matplotlib. A future render requires root to pin the result SHA and authorize one run with an exact token. The fixed limits are 15 seconds, 2 MiB including capture and raw child logs, and one attempt. The runner refuses to overwrite `run01` and records input hashes before and after execution.

Planned future render invocation after the actual result is reviewed and root enables GO:

```powershell
C:/Users/xiemo/Desktop/ShardForge/third_party/vesselforge-runtime/venv/Scripts/python.exe -B pgfrag/tmp/review_20261007/mixed93_opening_history_saved_visual/run_saved_opening_history_visual.py --render --root-go <ROOT_GO_TOKEN>
```
