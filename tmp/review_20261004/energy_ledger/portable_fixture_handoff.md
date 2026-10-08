# Energy ledger portable fixture handoff

- Diff: `pgfrag/tmp/review_20261004/energy_ledger/portable_fixture.diff`.
- New test freeze: `pgfrag/tmp/review_20261004/energy_ledger/freeze_v2_portable/`.
- Fixture identity: hand-authored energy-ledger regression fixture; copied byte-for-byte from the original review input.
- Original and portable fixture: 3707 bytes each, SHA-256 `559ddce00f53217e0c0992248d4b9b4078566c0557cc0bf09018e991c12ab81e`; exact byte equality verified.
- Updated test SHA-256: `e647e819e79e6fedf4879314410c4327a88552dbc6aeb786b0d05c6d685f24cd`. The only test edit changes `INPUT_PATH` to the checked-in test fixture.
- Energy source SHA-256 remains `6d7cea82c0e1bd105e143e3a3f5324595d2594c99201d9de7a0003d28496ad0e`.
- No tests were run in this handoff. The prior v2 receipt remains intact at `pgfrag/tmp/review_20261004/energy_ledger/v2_run/`; it covers the pre-portability test hash and reports 14 passed. The fixture-path-only edit adds no distinct tests.
