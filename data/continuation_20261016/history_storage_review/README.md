# Native vessel history storage review

This read-only audit covers the current history writer, independent audit, fragment exporters, and the actual zero-step full-size preflight. No solver evolution or storage change to PG source was made. Material values remain uncalibrated placeholders; this is a cost/layout review, not specimen qualification.

The 50 μs forecast contains 22,505 frames, 1,623 interface quadrature points, 3,708 patch nodes, and 18,634 native nodes after seam insertion. All primary histories are float64. Each frame writes 211,760 bytes of QP/work/ledger/patch values; each saved nodal frame adds 1,788,864 bytes.

| Node stride | Saved node frames | Raw numeric histories | Current wrapper history estimate | Current full-run allocation with margin |
|---:|---:|---:|---:|---:|
| 50 | 452 | 5,574,225,328 B | 5,574,233,662 B | 7 GiB |
| 1 | 22,505 | 45,024,043,120 B | 45,024,227,877 B | 51 GiB |

The raw chunked-NPY alternative preserves every numeric value and adds exact per-chunk NPY headers plus the declared manifest and saved-step index. It improves fail-closed prefix handling and chunked reads, but does not materially shrink the forecast. Compressed NPZ chunks are lossless and the synthetic test passed bitwise equality, shape/dtype/hash checks, and incomplete-temp exclusion; the actual vessel compression ratio is unknown because there are no evolved histories.

The independent auditor expects individual NPY files and performs several full-history NumPy operations. Both v5/v6 fragment exporters directly memory-map named NPY files. A chunk format therefore needs a versioned audit/export reader and streaming reductions before it can replace that contract. All QP frames remain necessary for failure/work checks; node state remains exactly at the declared stride (stride 1 for strict full-step nodal identities, stride 50 for sampled reconstruction), including the final frame.

See `history_storage_review.json` for field-by-field shapes, dtype/lifetime/consumer, exact chunk/header forecasts, compatibility constraints, source hashes, and migration gates. The prototype is under `pgfrag/experiments/continuation_20261016/history_storage_review/`.
