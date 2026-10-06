# Source and focused-run review record

This review area records a static implementation snapshot for the independent continuous-polyline Hausdorff-bounds metric. The experiment source, planned checks, and qualitative counterexample remain in `pgfrag/experiments/continuous_polyline_bias_bounds_20261006/`.

The initial SOURCE ONLY pass performed Python AST parsing, line counts, SHA-256 pins, and JSON structure parsing without calling the metric or tests. The ROOT-approved v1 stdlib unittest run is recorded in `focused_run01/`; the separate ROOT-approved v2 run is recorded in `focused_run02/`. Both retain raw output, return codes, elapsed times, pre/post source pins, and the metric return records from their own test suites. No pytest, rendering, provider, or solver was used.

The snapshot uses the repository evidence location `pgfrag/tmp/review_20261006/...`.

The result bounds apply only to the supplied piecewise-linear polylines in a shared frame. The module hard-codes `physical_admission=NOT_ADMITTED`; it makes no claim about an underlying continuous crack surface, solver lineage, or absence of orientation bias.

See `pins.json` for exact file hashes and line counts. The pre-fix source and pin receipt are retained as `source_snapshot_before_narrow_fix.py` and `pins_before_narrow_fix.json`; their hashes remain under `historical_record`. The current file list deliberately excludes the manifest itself because a file cannot carry its own stable SHA-256 digest.
