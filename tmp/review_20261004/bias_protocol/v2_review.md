# Bias protocol v2 review package

V2 addresses two comparison blockers identified during review:

1. All valid SHA-256 equality checks normalize hexadecimal case in comparison
   records. The input manifest remains unchanged.
2. Peridynamic direction and spatial contrasts hold explicit horizon fixed.
   A missing horizon yields UNKNOWN_NOT_READY; a changed horizon yields
   CONFOUNDED_NOT_COMPARABLE. Process-zone length, definition hash, and source
   hash remain separate from horizon. Two omitted process zones can support a
   diagnostic preregistration while resolution qualification stays UNKNOWN.
   One known/one unknown is UNKNOWN_NOT_READY; differing known process-zone
   identity yields CONFOUNDED_NOT_COMPARABLE.

Direction pairs also match explicit spacing, dt, and planned end time.
Spatial pairs match dt and planned end time. Time-step pairs use the same
artifact, spacing, horizon, process-zone identity, and planned end time. When
both completed run records are present, their actual_end_time_s must match;
missing or incomplete endpoint evidence is UNKNOWN_NOT_READY.

New cases cover mixed-case SHA values and unchanged caller data, upper/lower
spellings of the same artifact digest, a horizon change, missing horizon,
both process zones unknown, different known process-zone definitions, and
different actual dt-run end times. V1's six-pass output, argv receipt, source
hashes, and original source/test files are preserved in v1_freeze/.

The V2 Python sources passed AST parsing while the review diff was prepared.
V2 pytest has not been run and is waiting for root review. No mechanics solver
was invoked. Declared discretization signatures remain UNVERIFIED, and no
result claims that independence or bias removal has been established.
