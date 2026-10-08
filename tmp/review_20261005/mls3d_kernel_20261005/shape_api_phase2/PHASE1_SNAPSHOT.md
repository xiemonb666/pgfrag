# Phase 1 source preservation

Phase 1's implementation, acceptance source, and README are copied verbatim to
`phase1_snapshot/` before adding the direct shape-function API. Their SHA-256
values match the independently recorded and reviewed phase1 pre/post hashes.
The phase1 focused pytest record remains in
`../focused_test/RUN_RECORD.md`; it records one Syn-venv invocation, exit code
0, and six JUnit testcase elements. New phase2 sources do not alter those
phase1 bytes or retroactively change that run's claim.

The exact phase1 snapshot digests are listed in
`PHASE1_SNAPSHOT_SHA256.txt`. Phase1 passed files are only historical inputs
for this phase; the new phase2 acceptance source has not been run.
