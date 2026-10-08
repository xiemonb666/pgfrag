# Frozen verification receipt

The only pytest process that executed tests used the frozen `source/energy_ledger.py`, `tests/test_energy_ledger.py`, and `input/audit_cases.json` whose copies are under `freeze/`.

Final test command was launched through uv because the system `python.exe` name is a Microsoft Store placeholder and the bundled Python runtime has no pytest. The complete launcher argv, nested pytest argv, working directory, and relevant environment are in `argv.json`. Pytest `8.4.2` was supplied through a temporary uv cache under this review directory. It ran `tests/test_energy_ledger.py` with cache provider disabled and bytecode writes disabled.

`stdout.txt` records `7 passed in 0.02s`; `stderr.txt` is empty; `exit.json` records exit code `0`. `hashes_before.json` and `hashes_after.json` give byte counts and SHA-256 for the live source, test, and fixture input. All three hashes match before and after. The frozen copies under `freeze/` are the exact pre-run snapshots.

Two earlier launches did not execute pytest: `attempt_01_store_alias/` records the Microsoft Store Python placeholder exit `9009`, and `attempt_02_uv_module_as_program/` records uv treating `-m` as the program and exiting `2`. The final invocation corrected the launcher by running `python -m pytest` inside uv. No source, test, or input edits followed the successful run.
