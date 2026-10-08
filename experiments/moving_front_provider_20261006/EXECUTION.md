# Moving-front provider one-shot execution

**Outcome: BLOCKED before producing a diagnostic result.** The one authorized driver invocation returned process exit code `1` after `0.237516 s`. The run-once marker remains present; no result JSON was produced. This is a pin-guard failure, not a diagnostic PASS.

## Frozen inputs and command

The exact pending spec bytes are preserved in [moving_front_check_spec.source_only.snapshot.json](moving_front_check_spec.source_only.snapshot.json), SHA-256 `114e45e2046e696fad655a5d1976364abc50345475d253a2e43879fe56ed4109`. The sole edit to the live spec was `status: SOURCE_ONLY_AWAITING_ROOT_READY` → `status: READY`; READY spec SHA-256 is `5a3b480d36ab1830904d6f76c9b4d27005327546aae6f3e15dd16f8a7e677669`.

VF interpreter: `C:\Users\xiemo\Desktop\ShardForge\third_party\vesselforge-runtime\venv\Scripts\python.exe`, SHA-256 `2b83f3f84ea2f51132fa1aae73815be54126a799521de85ff6ea70627cde6da7` (matches the prior VF run receipt). The one invocation was:

```text
C:\Users\xiemo\Desktop\ShardForge\third_party\vesselforge-runtime\venv\Scripts\python.exe -B C:\Users\xiemo\Desktop\ShardForge\pgfrag\experiments\moving_front_provider_20261006\moving_front_check_driver.py
```

stdout/stderr were redirected to [moving_front_check.stdout.txt](moving_front_check.stdout.txt) and [moving_front_check.stderr.txt](moving_front_check.stderr.txt). The stdout file is empty (0 bytes; SHA-256 `e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855`). Stderr is preserved byte-for-byte (1,279 bytes; SHA-256 `e5c3b3a9aa040770807f014f58f01c8c7396e2c3f5f034c25147aca6cda7c9db`):

```text
Traceback (most recent call last):
  File "C:\Users\xiemo\Desktop\ShardForge\pgfrag\experiments\moving_front_provider_20261006\moving_front_check_driver.py", line 1014, in <module>
    raise SystemExit(main())
                     ^^^^^^
  File "C:\Users\xiemo\Desktop\ShardForge\pgfrag\experiments\moving_front_provider_20261006\moving_front_check_driver.py", line 1002, in main
    result = run_check(frozen)
             ^^^^^^^^^^^^^^^^^
  File "C:\Users\xiemo\Desktop\ShardForge\pgfrag\experiments\moving_front_provider_20261006\moving_front_check_driver.py", line 874, in run_check
    bundle = mf.load_pinned_inputs()
             ^^^^^^^^^^^^^^^^^^^^^^^
  File "C:\Users\xiemo\Desktop\ShardForge\pgfrag\experiments\moving_front_provider_20261006\moving_front_provider.py", line 251, in load_pinned_inputs
    "moving_front_cases.json": _check_pin(
                               ^^^^^^^^^^^
  File "C:\Users\xiemo\Desktop\ShardForge\pgfrag\experiments\moving_front_provider_20261006\moving_front_provider.py", line 189, in _check_pin
    raise MovingFrontBlockedError(f"pinned source changed: {relative_path}")
moving_front_provider.MovingFrontBlockedError: pinned source changed: pgfrag/experiments/moving_front_provider_20261006/moving_front_cases.json
```

The guard stopped in `load_pinned_inputs`, before provider queries or the two planned MLS batches; no check records or result JSON exist. The 10 spec-listed files matched their hashes both before and after the invocation. Thus the failure is the provider's own pinned-source check for `moving_front_cases.json`, despite that file matching the READY spec pin; no hashes were changed to bypass it.

## Source pins before and after

All listed hashes matched the READY spec before the run and remained unchanged after it:

| Source | SHA-256 |
|---|---|
| `moving_front_check_driver.py` | `02314dbc2c4160f6687f27cc67c25b6c7f2ddd3dcde0b9766dbad427512d539b` |
| `moving_front_provider.py` | `b3e88622220fae0ed63961b943ba710e5689b470d7df00361c297934c42dfac6` |
| `moving_front_cases.json` | `2766245b25bf93eefcf6836befe81536917888936b4543666c86ebe9d6b72100` |
| safe-v2 `support_topology.py` | `919cb60a6d65ae6b55a0bd12b569ec6627066b2a72dfcca52145209336588664` |
| `fixed_case.json` (source contract) | `0bdbb20032e105851d2909ff54850ebc54d3f2041c25ab5e42e9102244d90ca5` |
| `role_manifest.json` | `e5e318f80b5dfedbc8ba122f82df3cca1585e604272399d21bc35907d1d28f01` |
| `enriched_dof_contract.py` | `e3ddb3b067465bb5b269576a1b8a5100e2bb92c9c924bc09c949702c19efa9d8` |
| `mls3d_mls.py` | `404256d9b7546361a5704c4c3cec7c022069a9cf83a137c58a16e7a64f79ec63` |
| `mls3d_shape.py` | `66e680616a988182875d17aeb84fadb541867aee13fc15e0ae098f7c98839ef5` |
| material `fixed_case.json` | `dce0a1f3009246a52198bc0ed833ed67998144549a653e7bce9bcc17babd71b9` |

## Run-once artifacts

The driver-created [moving_front_check.run-once](moving_front_check.run-once) marker is preserved (241 bytes; SHA-256 `7c6ee6012eaef6cc1c06d433358ea2c71c1c0e9c5c2d0f08c9d11bf67f111066`). It binds the READY spec SHA and driver SHA and states that failures do not clear the marker. `moving_front_check_result.json` is absent. Captured artifacts are well below the 2 MiB output cap. No retry, solver/operator run, coefficient transfer, or physical crack-propagation claim was made.