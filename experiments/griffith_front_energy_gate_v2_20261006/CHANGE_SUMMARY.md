# V2 change summary relative to v1

V2 keeps the v1 geometry proposal flow and replaces its confirm policy with a separately typed quasi-static energy gate. The ambiguous v1 `delta_potential_J` is not consumed, remapped, or sent through v1 `confirm_candidate`; the v2 snapshot declares `legacy_delta_potential_bridge = NOT_BRIDGED`.

The v2 energy receipt names stored elastic energy and external work separately. It derives `delta_total_potential = delta_stored_elastic - external_work` only for reporting and computes the closure residual exactly once as `external_work - delta_stored_elastic - Gc_delta_area`. Commit requires finite numeric energy terms, units, positive reference energy, nonnegative absolute/relative tolerances, a tolerance basis, exact candidate/parent/state bindings, geometric area consistency, and `abs(residual) <= abs_tol + relative_tol*reference_energy`.

The gate supports only `QUASISTATIC_ELASTIC_NO_INERTIA_NO_OTHER_DISSIPATION`; dynamic or mixed-cohesive inputs return `UNKNOWN`. Missing inputs or over-tolerance residuals return `NEEDS_RECOMPUTE` and preserve the active parent surface. A fixture path yields `COMMITTED_FIXTURE_ONLY`; the external current-reanalysis path remains available as `COMMITTED_EXPERIMENTAL_REANALYSIS`, with the parent surface contract still `NOT_COUPLED` and no physical validation claim.

All v2 code, tests, and artifacts are in `griffith_front_energy_gate_v2_20261006` directories. The v1 source, test, and receipt were not edited. Receipt hashes freeze the two read-only source dependencies used by the runner.
