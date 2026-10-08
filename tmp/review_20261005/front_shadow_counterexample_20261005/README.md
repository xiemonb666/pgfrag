# Front-shadow counterexample (2026-10-05)

This folder contains a frozen five-source SI engineering fixture and an independent algebra note for a source-front visibility shadow switch. The query path stays at least `sqrt(2) m` from the finite crack triangle at `q0`, while the segment from the special source crosses the triangle interior for `t < 0` and misses it for `t > 0`.

The one-sided `q0` MLS values in `math.json` are labeled `ONE_SIDED_ALGEBRAIC_TRACE_ONLY`. Exact visibility at `q0` is an edge contact and is `UNKNOWN`; no `q0` gradient is claimed. The fixture is not a physical-fracture claim, not mechanics-coupled, and does not establish solver or R2 acceptance.

Root reviewed the independent algebra and approved one bounded source-only helper run. `run_front_shadow_counterexample.py` reads only the frozen case and the three SHA-256-pinned helper sources, evaluates the six fixed nonzero queries, and writes one transparent PNG with `--plot`. It refuses to overwrite prior results. Limits are five sources, six queries, 30 seconds, and 2 MiB total output.
