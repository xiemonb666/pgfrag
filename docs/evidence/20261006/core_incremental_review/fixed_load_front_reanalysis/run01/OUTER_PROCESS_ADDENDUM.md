# Run01 outer-process receipt addendum

This is a post-run reconstruction from the existing Codex tool response and the saved result, not a contemporaneous raw process log. No process was rerun to create it. The structured record is [outer_process_receipt.json](outer_process_receipt.json).

The sole invocation used `C:/Users/xiemo/Desktop/ShardForge/third_party/vesselforge-runtime/venv/Scripts/python.exe` with `-B`, script `pgfrag/experiments/fixed_load_front_reanalysis_20261006/fixed_load_front_reanalysis.py`, `--execute`, and output `pgfrag/tmp/review_20261006/fixed_load_front_reanalysis_20261006/run01/result.json`, from workspace root `C:/Users/xiemo/Desktop/ShardForge`. The wrapper reported child exit code 0 and elapsed time 0.5507252 s under its 60 s timeout.

The textual tool response is transcribed below as returned. It preserves the displayed text only; stdout and stderr raw bytes were **NOT_PERSISTED**. The response included no `STDERR=` line, so this record makes no claim about the exact stderr byte stream.

```text
EXEC_SECONDS=0.5507252
EXEC_RC=0
STDOUT={"case_id": "modeiii_pu_full3d_fixed_load_front_reanalysis_20261006", "status": "FINITE_CHANGE_RELEASE_DIAGNOSTIC_ONLY_NOT_COMMITTED", "result_path": "C:\\Users\\xiemo\\Desktop\\ShardForge\\pgfrag\\tmp\\review_20261006\\fixed_load_front_reanalysis_20261006\\run01\\result.json", "operator_sha256": "3b31e7d40c1d8e4b5d829f0b093fee4edadc966bde4086a4b3f4ba48fdf89cd2", "baseline_static_solve_sha256": "ef66783e2961dd77f34785a99e1b397210ccc0beb6e8fa2979678aa4992a0ce3"}
exit_code=0
```

Preflight stdlib-only checks recorded core SHA-256 `9707fa3c47403e0dada8081e3fba9273249e9a05fc2055b0e3a20a69b62c952c`, case SHA-256 `8f919a55127b6f310e6f262b4a550f8daa9cd214496d824c3e1848ed0b1b46df`, frozen `SOURCE_PINS.json` SHA-256 `16355c781749fb4940a187afb3d985d9f1362fd3a33c1e6cc425cffdf9ed0ac7`, operator-result SHA-256 `3b31e7d40c1d8e4b5d829f0b093fee4edadc966bde4086a4b3f4ba48fdf89cd2`, and baseline-result SHA-256 `ef66783e2961dd77f34785a99e1b397210ccc0beb6e8fa2979678aa4992a0ce3`. The child response and saved result repeat the operator and baseline pins. The child did not echo a post-run `SOURCE_PINS.json` digest; its source-validation path completed without blocking, and the preflight pin remains the separately recorded observation.

The result is 539,566 bytes with SHA-256 `a60e51c85c26510c358ee57cafa4bfdda6ed89ea463d99bcc86551f64220732c`. Its residual criterion is a numerical solve-trust check only, not a physical fracture or growth PASS. The main review flagged an approximately `0.013` finite-G discrepancy and old-tip quadrature risk; an independent critic is recalculating. This receipt does not resolve that review, alter the result or threshold, compare against `Gc`, or claim physical validity. A future immutable per-function evidence record should be planned and persisted before execution; this post-hoc receipt is not a substitute.
