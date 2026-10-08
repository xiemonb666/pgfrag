# Moving-front provider source review

**Disposition requested: root review and READY decision before one actual check.** These artifacts define two source-only front states and a one-shot bounded provider check. The provider and driver have not been imported or executed; no source classification, MLS query batch, provider sample, operator, or solve was run. The check spec remains `SOURCE_ONLY_AWAITING_ROOT_READY`; the driver rejects execution until root changes only that status to `READY`.

## Files in scope

- `pgfrag/experiments/moving_front_provider_20261006/moving_front_cases.json` — immutable input contract for fronts at `x=0` and `x=0.0005 m`, shared box and 100 material-rest source coordinates, and 25-query-per-state / 254-record plan budget.
- `pgfrag/experiments/moving_front_provider_20261006/moving_front_provider.py` — source-only provider proposal. It reads only source IDs/lattice indices/positions from the old coordinate manifest; it reclassifies all 100 sources per front through hash-pinned safe-v2, loaded under a private hash-specific module name. Each front receives a fresh BASE/H/shifted-F1 registry and current-front shifted source values. H and F1 gradients use product rules. Query traces require exact PLUS/MINUS tags on strict finite-Γ interior; front, endpoints, ambiguous cut-band, and outside queries fail closed. It does not implement coefficient mapping or solve behavior.
- `pgfrag/experiments/moving_front_provider_20261006/moving_front_check_driver.py` — source-only one-shot driver. It counts safe-v2 source classifications, topology-only controls, and calls to the MLS shape API; it accepts no CLI inputs, writes an exclusive one-run marker only after READY and hash validation, and refuses a repeat run.
- `pgfrag/experiments/moving_front_provider_20261006/moving_front_check_spec.json` — frozen driver/provider/case/dependency hashes, exact tolerances, query layout, core-radius rule, memory/time/query caps, and pending READY status.
- `pgfrag/experiments/moving_front_provider_20261006/FUTURE_CHECK_PLAN.md` — the same bounded check in human-readable form.

The future workload is 50 provider point/trace rows + 200 current-role classifications + 4 synthetic tangency probes = exactly 254 records, with a hard cap of 256. Each fresh role state is built once; the reversed-order semantic-ID check rebuilds only the registry and adds no topology classification or MLS batch. The four exact-side rows are checked against one-sided H/F1 basis formulas. Core radii use `min(0.25 mm, 0.25*r_core)` and `min(0.50 mm, 0.50*r_core)`; a nonpositive core blocks without searching for replacement radii. Dynamic roles, DOFs, core radii, and errors remain unknown until root approves the check.

Frozen tolerances are affine displacement `1e-10 m`, affine gradient `1e-9` dimensionless, core displacement `1e-12 m`, core gradient `1e-9` dimensionless, trace basis value `1e-12`, and trace basis gradient `1e-9 m^-1`. The driver caps time at 30 seconds and diagnostic-owned arrays + a 512 KiB transient reserve + serialized result at 2 MiB; interpreter/import baseline is excluded. Result JSON is limited to 256 KiB. Collar errors are recorded without a threshold. The ahead-of-front ±`1e-8 m` rows are a finite-sample diagnostic, not a global continuity proof.

## Source pins recorded

| Input | SHA-256 |
|---|---|
| safe-v2 `support_topology.py` | `919cb60a6d65ae6b55a0bd12b569ec6627066b2a72dfcca52145209336588664` |
| fixed rest-source case | `0bdbb20032e105851d2909ff54850ebc54d3f2041c25ab5e42e9102244d90ca5` |
| coordinate-manifest file | `e5e318f80b5dfedbc8ba122f82df3cca1585e604272399d21bc35907d1d28f01` |
| coordinate records | `a8687a23353b6465aa0cc3fcda49a57c240297c38c141a017a167f3c526b23bb` |
| registry contract | `e3ddb3b067465bb5b269576a1b8a5100e2bb92c9c924bc09c949702c19efa9d8` |
| MLS core | `404256d9b7546361a5704c4c3cec7c022069a9cf83a137c58a16e7a64f79ec63` |
| MLS shape | `66e680616a988182875d17aeb84fadb541867aee13fc15e0ae098f7c98839ef5` |
| manufactured material fixture | `dce0a1f3009246a52198bc0ed833ed67998144549a653e7bce9bcc17babd71b9` |

The provider also verifies imported registry/MLS module origins against the pinned files, so a same-name cached module from another path blocks use. Its source-position matrix is backed by immutable `bytes`, so the NumPy writeable flag cannot be re-enabled. The case-pair JSON pin embedded in the source is `2766245b25bf93eebcf6836befe81536917888936b4543666c86ebe9d6b72100`.

## Limits for root review

This is a two-geometry source path, not a generic moving-front API. It depends on the fixed closed box, fixed material-rest coordinates, straight finite front, `h=0.007 m`, `R0=0.01 m`, and safe-v2's current classifier contract. Hashes and coefficient provenance force each geometry/basis generation to be treated independently; cross-front old coefficients are rejected. No coefficient reprojection is supplied. The manufactured field is centered on each current front for mathematical verification only; it is not a new solution for the original fixed loading.

The check plan treats exact crack-face samples as one-sided basis checks and tests local manufactured-field agreement only strictly inside the measured current core. Cross-front provenance rejection uses a synthetic provenance tag only; no coefficient values are loaded, transferred, or solved. Collar errors are observations, not a global exactness criterion. No growth, production coupling, or general crack topology is qualified here.

## Artifact hashes

| Artifact | SHA-256 |
|---|---|
| `moving_front_cases.json` | `2766245b25bf93eebcf6836befe81536917888936b4543666c86ebe9d6b72100` |
| `moving_front_provider.py` | `b3e88622220fae0ed63961b943ba710e5689b470d7df00361c297934c42dfac6` |
| `moving_front_check_driver.py` | `02314dbc2c4160f6687f27cc67c25b6c7f2ddd3dcde0b9766dbad427512d539b` |
| `moving_front_check_spec.json` (pending status) | `114e45e2046e696fad655a5d1976364abc50345475d253a2e43879fe56ed4109` |
| `FUTURE_CHECK_PLAN.md` | `9a2c6ed958b6781bc5ed1b81f0424873c8fe8b0a0f3eab6e596eedb992c1c494` |
