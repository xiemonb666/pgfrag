# Griffith front step: one-step experimental controller

This source-only experiment implements one bounded proposal/confirm/reject cycle for a straight finite Mode III front. The accepted frame is a finite line parallel to `z`, with one declared propagation direction, `+x`. It does not model crack initiation, turning, branching, front splitting, crack speed, or fragment counts. Every generated surface version stays `NOT_COUPLED` and `EXPERIMENTAL_UNCALIBRATED`.

The existing [`explicit_surface_state.py`](../explicit_surface_state_20261005/explicit_surface_state.py) contract requires numeric side-displacement arrays when it constructs a `CrackPatch`. A proposal therefore stores an immutable geometry-only strip candidate and leaves the active committed state untouched. It does not fill a future displacement field with zeros. Confirmation accepts a new `ExplicitSurfaceState` only after it contains the candidate strip, unchanged old geometry and assignments, and a supplied hash for that exact reanalyzed state.

## State transitions

`propose_step` needs a current driving assessment bound to the exact active surface hash. It requires finite `G` and positive `Gc` in `J/m^2`, solver-result and source hashes, an explicit `COUPLED_CURRENT` declaration, and reviewed contour/path sensitivity receipts. `REVIEWED` means that a named review record is present; it is not an assertion that a numerical threshold passed or that J is path independent. The controller computes the strict comparison itself: `G < Gc` records arrest, `G == Gc` holds at equality, and `G > Gc` creates one geometry candidate. Invalid/missing inputs produce `UNKNOWN`, `NEEDS_RECOMPUTE`, or a named rejection status.

A `PROPOSED_WAITING_REANALYSIS` candidate is never the active surface. `confirm_candidate` requires a current reanalysis receipt bound to the same parent state and candidate geometry, the hash of the supplied reanalyzed surface state, solver/source hashes, an energy-ledger hash, and finite `delta_potential_J`, `external_work_J`, and `Gc_delta_area_J` plus an energy-review receipt. An external receipt must state `REANALYZED_CURRENT`; the test path uses the separate `EXPLICIT_FIXTURE_ONLY` status. The supplied new state must include the candidate patch with its newly solved side displacements while preserving old patch geometry and ownership. Missing or mismatched evidence leaves the prior immutable version active and marks the candidate `NEEDS_RECOMPUTE`; `reject_candidate` records rejection and also retains the prior active version.

The required sign convention is `delta_potential = potential(candidate) - potential(parent)`, positive `external_work` means work into the modeled system, and positive `Gc_delta_area` counts the common new crack area once. The energy ledger records

```text
external_work_J - delta_potential_J - Gc_delta_area_J
```

as a signed residual. No tolerance or `physics_passed` flag is accepted. `COMMITTED` means the required evidence fields were supplied and their hashes/geometry bindings agreed; the controller does not verify solver receipts, assess the physical adequacy of the energy method, or claim energy balance. Production promotion remains blocked by the surface-state contract.

## Inputs required from real components

The current Mode III J integrator consumes an already sampled field and returns a status-bearing J value; it does not run the mechanics solve or independently verify field provenance. Before a real proposal can be made, a J/mechanics integration layer must provide:

- `G` from a current J result with status `OK` in `J/m^2`, `J_result_id`, the exact `ExplicitSurfaceState` hash on which the field was solved, solver-result SHA-256, and source SHA-256 for the field producer/J component;
- contour sensitivity status and a hash-addressed review record, plus a separate path sensitivity status and review record, all tied to that same field/surface revision; the existing integrator does not supply these reviews;
- `Gc` in `J/m^2`, a separate material/source SHA-256, material provenance, and the explicit `EXPERIMENTAL_UNCALIBRATED` qualification;
- the declared finite front coordinates, `+x` direction, allowed step cap, and reference-domain bounds.

After a proposal, the mechanical reanalysis must solve that exact candidate geometry and return the new `ExplicitSurfaceState` with actual side-displacement arrays; its solver-result/source hashes; the parent-state and candidate-geometry hashes; and signed `delta_potential_J`, `external_work_J`, `Gc_delta_area_J`, and area-increment evidence with units, an energy-ledger hash, and a hash-addressed energy review. The controller recomputes the candidate strip area and `Gc × area` arithmetic for consistency, but it leaves the physical energy interpretation and review method to that upstream record. Any missing drive, absent coupling, stale surface/field, or unchecked energy record cannot activate a new surface.

The hashes and status strings are receipt declarations. This isolated module has no verifier that proves those declarations came from a genuine solver. The upstream adapter must provide that trust boundary before these transitions are used with real results.

## Checks and artifacts

Run from the repository root with the bundled Python runtime:

```powershell
& 'third_party\vesselforge-runtime\venv\Scripts\python.exe' pgfrag\experiments\griffith_front_step_20261006\run_fixture.py
```

The 15 short tests use explicit in-memory fixtures only. Their successful transition is labelled `COMMITTED_FIXTURE_ONLY` / `EXPLICIT_FIXTURE_NOT_SOLVER`; it validates state-machine wiring and is not a solver success. The test runner writes an actual execution receipt, four serialized state snapshots, and a proposed-versus-committed diagram under `pgfrag/tmp/review_20261006/griffith_front_step_20261006/`.
