# PG energy-ledger audit core

This artifact adds a reusable **history auditor**, not a solver adapter. Its input is a finite, strictly ordered sample clock, per-channel signed raw values with units/definitions/source keys/source/clock IDs, explicit tolerances, and a case qualification. Missing expected channels are emitted as `UNKNOWN`; the auditor does not synthesize zero values or fuzzy-match getter names.

## Source-grounded boundaries

- Akantu v5.0.7 `MaterialCohesive::getEnergy` defines `reversible`, `dissipated`, and the exact `cohesive contact` selector. Its fallback returns zero for any unrecognized name. Therefore a trace obtained with `getEnergy("contact")` cannot establish a zero contact energy. The exact source is `third_party/akantu-runtime/source-v5.0.7/src/model/solid_mechanics/solid_mechanics_model_cohesive/materials/material_cohesive.cc`; the existing capability report says wheel-to-source correspondence remains UNKNOWN and also records the contact selector issue in `relicweave-synlab/docs/a12/COHESIVE_CAPABILITY_20261004.md`.
- Akantu's reviewed `potential` is the bulk material energy. Cohesive reversible, dissipated, and contact terms come from the cohesive law. The model-level `external work` getter in the inspected source is an increment per call; the retained extrinsic trace used its own cumulative reaction/displacement estimate. The negative cohesive dissipation values remain raw: `-5.572968137645039e-5 J` and `-9.726359554870569e-6 J`; the report already classifies both as FAIL. The auditor checks those signed endpoints against the explicit `1e-12 J` roundoff floor without applying `abs()` or clamping.
- The intrinsic 201-point Mode-I run's retained summary reports native keys `kinetic`, `potential`, `dissipated`, `reversible`, and `cohesive contact`; final `D=3.0399999999999994e-5 J` matches `Gc*A=3.0399999999999997e-5 J` to relative error `1.1145170358609216e-16`. All XYZ DOFs are prescribed, so its free-energy balance is `NOT_APPLICABLE`. This is a toy constitutive comparison only; material calibration remains `UNKNOWN`.
- Peridynamics.jl OSB computes a point strain-energy density from dilatational and deviatoric terms, and its exporter returns that field. `pgfrag/src/pgfrag/resources/julia/run_fixture.jl` multiplies it by point volumes, computes kinetic energy from `0.5*m*v²`, and records total-energy conservation as UNKNOWN because prescribed boundary velocities do external work and the fixture lacks a complete ledger. `short_range_force_contact.jl` assembles contact force density but does not define/export a contact-energy history. The auditor therefore requires any PD contact-energy channel to be separately measured and qualified before including it.

The exact retained numbers and source paths used by the tests are frozen in `input/audit_cases.json`.

## Interface

Import `Clock`, `EnergyChannel`, `EnergyHistory`, `EnergyTolerance`, `EvidenceRef`, `BoundaryWorkDefinition`, `Qualification`, and `audit_energy_history` from `pgfrag.energy_ledger`.

Channels declare an exact upstream `source_key` and one role (`kinetic`, `stored`, `dissipated`, `external_work`, or `diagnostic`). `Qualification.required_balance_channels` names the complete channel set for that case. Every required key must have values, matching clock metadata, and its own structured `EvidenceRef(source, method)`. A status string such as `PASS` does not count as evidence. For a free-dynamics calculation, qualification must additionally provide a physical time clock, independent-balance evidence, cumulative dissipation/work, and an explicit boundary-work definition with sign convention `positive_into_system`.

Residuals are evaluated against the first supplied row, which is reported as the baseline; the same baseline is subtracted from each cumulative channel. The residual is `Δ(kinetic + stored) + Δdissipated - Δexternal_work`. At each row the allowance is `max(abs_tol_J, rel_tol * sum(abs(each baseline delta term)))`. The result describes arithmetic closure for caller-supplied channels; it cannot authenticate their provenance, establish solver correctness, or calibrate material parameters.

The all-DOF-prescribed mode always reports free energy balance `NOT_APPLICABLE`. A separate final `D` versus `Gc*A` comparison can be reported as `MATCHED_REFERENCE`; it does not change calibration or free-fracture status.

## Bounded verification

One short pytest invocation covers the two retained negative-dissipation magnitudes, raw tiny negative roundoff retention, missing contact/clock, incompatible energy units, all-DOF `NOT_APPLICABLE` with `Gc*A` arithmetic, a synthetic independent analytic conservation control, and rejection of status-only evidence. No native solver is invoked, no existing experiment or report is modified, and no long history is loaded.
