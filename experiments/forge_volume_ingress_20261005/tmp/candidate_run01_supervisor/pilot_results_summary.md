# Bounded VF-to-PG ingress pilot: retained-run summary

Status: **RECOVERED_DIAGNOSTIC_FROM_FAILED_REPORT_RUN**. This is a read-only summary of the one authorized pilot; it does not promote the blocked stage or authorize a rerun.

The worker used the pinned VF Python and exited 2 after 4.616230 s. It processed the fixed 96-cell ROI and wrote its ledgers and pointset, then failed report finalization with `could not stabilize report-inclusive output byte count`. The retained stage is 325,188 bytes including all files, below the 20 MiB cap. Inputs and source-freeze hashes match before and after.

## Local geometry result

- ROI `[7.5,10.5] × [-1,1] × [4,6] mm`, h = 0.5 mm; 96 / 96 cells processed.
- Cell outcomes: 80 positive, 16 zero, 0 blocked.
- Positive connected components: 80; centroid classes: {'IN': 80}; component topology: {'PASS': 80}.
- ROI common volume: 8.01676271854943 mm³. Sum over 96 cell commons: 8.01676271854445 mm³; difference -4.9809e-12 mm³, tolerance 8.01676e-08 mm³. Retained raw ledger result: `PASS_LOCAL_ROI_BREP_VOLUME_RECONCILIATION`.
- Candidate points: 80; point-volume sum 8.01676271854445 mm³; mass sum at declared 2100 kg/m³ is 1.68352017089433e-05 kg.
- Scope remains local B-rep geometry only. Surface-parent mapping, native PG integration, pair classification, fragment extraction, global wall, full-mother conservation, and physics qualification remain unestablished or unknown.

## Failure and preserved evidence

The worker raised `IngressBlocked`: `could not stabilize report-inclusive output byte count`. `_write_final_report` resets `output.total_bytes` to zero on every iteration, so it does not feed the newly computed total into the next fixed-point step. The worker wrote `worker_failure.json`; the parent wrote `parent_execution_receipt.json` and retained the staging directory. The stage has no `ingress_report.json`, and the pointset file is only an intermediate raw candidate in the blocked stage.

Supervisor receipt: `tmp/candidate_run01_supervisor/final_receipt.json`. Stdout: `tmp/candidate_run01_supervisor/stdout.log`; stderr: `tmp/candidate_run01_supervisor/stderr.log`.
