# Mode III scalar traction solve, staging-integrity v3

V3 is a separate source snapshot based on the frozen v2 scalar UZ solver. Its only behavioral change is in the parent process: a worker staging result is promoted only after exit code zero, strict UTF-8/JSON parsing, exact v3 schema/case/success status, and validation of required finite replayable evidence shapes. Nonzero exit, timeout, missing or invalid staging JSON, failed result-contract checks, or promotion errors write a separate UNKNOWN receipt and preserve the staging file where it exists.

The coefficient registry, continuous MLS/front-PU provider, stiffness/load/mean assembly, direct scaled 55×55 Lagrange solve, quadrature, numerical gates, finite-value policy, 30-second wall budget, 32 MiB numeric-array budget, 4 MiB result limit, and saved K/c/f/ID/off-grid evidence are unchanged from v2. This source remains unrun pending root review and authorization.
