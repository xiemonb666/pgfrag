# PG-Frag independent-discretization bias protocol

This review module defines pgfrag.bias-study.v1. It checks a pre-registered
comparison design and the provenance attached to reported metrics. It does not
generate a mesh or point cloud and does not invoke a solver.

All valid SHA-256 inputs are compared as lowercase hexadecimal values. The
manifest itself is not modified; reports normalize comparison keys separately.

## Manifest contract

Each study has a registration record, a fixed physical_case, and a list of
planned or completed runs.

physical_case fixes the source physical object, boundary/load definition,
material model and parameter hash, rest-coordinate material-field hash and
seed/hash, and solver contract hash. These values must agree across studies
before their runs can be compared. Missing identity evidence leaves physical
comparability UNKNOWN; a material or physical-case mismatch returns
NOT_COMPARABLE.

Every run declares its direction, discretization kind, generator and source
hash, generation seed, explicit spacing (h for tetrahedral meshes or dx for
PD points), solver dt_s and end_time_s, and construction lineage. For PD
contrasts, an explicit, valid horizon_m must be present and equal across the
pair. A missing horizon makes that contrast UNKNOWN_NOT_READY; a changed
horizon makes it CONFOUNDED_NOT_COMPARABLE.
Discretization signatures include an exact artifact SHA-256. A
rotation_invariant_discretization declaration also carries its digest,
algorithm, source artifact digest, and calculation-evidence digest. The
current checker validates their shape and linkage, then reports the
rotation-invariant claim as UNVERIFIED: it does not recompute signatures
from mesh or point payloads.

Lineage may be independent_generation, rigid_rotation, refinement, or
unknown. An exact artifact hash match is a duplicate even if seeds differ.
A rigid-rotation parent hash that points to another run is classified as the
same source rotation. Equal rotation-invariant digests are also excluded from
independent-remesh evidence, while their caller-provided digest remains
unverified. Different hashes plus distinct seed/direction and generated
lineage only make an UNVERIFIED_INDEPENDENT_CANDIDATE; they do not prove
independent construction.

PD point runs carry explicit dx. A native horizon_m and process-zone length,
definition hash, and source hash are separate fields. If both runs omit the
process-zone length, they may still form a diagnostic design contrast, while
peridynamic_resolution_qualification stays UNKNOWN. If one side is known and
the other is unknown, the contrast is UNKNOWN_NOT_READY. Two known
process-zones must match in length, definition hash, and source hash; any
difference makes the contrast CONFOUNDED_NOT_COMPARABLE. The checker never
substitutes horizon for process-zone length or infers an exact wall-layer
count from thickness, t, or a maximum mesh edge.

## Readiness rules

comparison_readiness reports READY_NOT_RUN when the physical case is
comparable, a pre-registration record is present, the manifest includes the
required contrasts, and no run evidence is recorded:

1. distinct generated seeds and directions at matching spacing and dt;
   direction comparisons also hold end_time_s fixed;
2. at least two explicit spacings with matching direction, dt, and end_time_s;
   and
3. at least two dt values on the same exact discretization, explicit horizon,
   process-zone identity, and planned end_time_s.

When completed dt runs are present, their actual_end_time_s values must also
match. A changed actual endpoint makes the dt contrast
CONFOUNDED_NOT_COMPARABLE; missing or incomplete endpoint evidence makes it
UNKNOWN_NOT_READY.

A shared rigid rotation cannot meet the independent seed/direction criterion.
These checks establish planned contrasts only. They do not claim that
discretization independence was verified or directional bias was removed.
If any run evidence is recorded, the status becomes
RUN_EVIDENCE_PRESENT_ANALYSIS_PENDING until a separate result analysis is done.
A separate objectivity check is marked DESIGNED_NOT_RUN when a run's
rigid-rotation lineage points to another recorded artifact. Otherwise it is
NOT_TESTED. Neither status carries an objectivity result before actual outputs
are supplied.

## Metric evidence

The report recognizes energy, response, mass_distribution, crack_direction,
and fragment_distribution. A metric value is surfaced only when its
provenance identifies an actual completed run, one of that run's output
hashes, and a qualified surface whose source receipt matches the run receipt.
Otherwise it is reported NOT_RUN with a null value, including geometry-only
data. A complete chain is PROVENANCE_LINKED_UNVERIFIED because this checker
does not rehash the source output or surface bytes. No fragment count is
required or synthesized.

The pgfrag point-output contract currently provides material coordinates,
point mass, energy density, bond state, and endpoint bond components; it does
not provide a reconstructed crack surface or verified rigid-fragment
geometry. The current mesh/point signature and lineage inputs are declarations
until a future protocol reader can recompute them from the frozen artifacts.
