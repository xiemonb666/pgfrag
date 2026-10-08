# Independent free-fracture pilot review

Independent audit status: **FAIL**. This status covers source/result byte binding and numerical reconstruction only. It does not establish physical fragments, a crack surface, convergence, calibrated material response, surface validity, grid/horizon independence, or fracture dissipation.

Native arithmetic validation status: **NOT_COMPLETED** (0/0 required independent numeric checks completed). This is reported separately from the execution-receipt bookkeeping status.


Receipt bookkeeping discrepancy (exact, zero-byte tolerance): the receipt records 88924215 attempt bytes, while the independently counted completed attempt tree has 88924268 bytes (observed − recorded = 53 bytes). The original receipt remains unchanged; numeric review status is shown separately.

The audit uses the completed one-shot receipt as an execution/provenance gate, verifies frozen source pins and output hashes, then reconstructs activity directly from the little-endian state tensor, canonical bond roster, point coordinates, and event threshold rule. Producer summary numbers are hash-bound for provenance but are not used as numerical oracles.

The state and graph tensors are independently decoded in Julia column-major order and checked against the sidecar field names, units, counts, shapes, byte totals, and SHA-256 digests. At every state the reviewer recomputes each active pair's stretch from the recorded deformed point coordinates; first breaks must be permitted and strictly beyond the CriticalStretch threshold, except that events within a 128-ULP comparison band are explicitly counted as roundoff-ambiguous. It then rebuilds per-point active degree and damage, irreversible transitions, union-find components, minimum-point-ID graph labels, component mass, event totals, and graph CSV metrics.

The native grip path is checked at all 97 state times against the fixed physical-time cubic ramp/hold/unload/dwell profile. The right grip's stored half-step velocity is checked against the secant of that path; free displacement drift, free half-step velocity, stored full velocity, acceleration/force-density relation, discrete grip impulses and work, momentum update, and the energy/reaction CSV are independently recomputed. Raw full-step kinetic energy is split into free and grip-masked DOFs. Since the native boundary condition constrains vhalf and the stored full-step velocity includes the post-force acceleration half-step, physical kinetic energy at constrained grip DOFs remains **UNKNOWN**.

Graph labels describe only connectivity in the active canonical-bond graph. A component count, even one isolated from an end grip, is not evidence of a physical fragment or an internal crack. Event locations and component bounding boxes are shown only as actual geometry and actual sampled first-break bond segments; the reviewer fits no crack plane and invents no surface. Material calibration, physical loading validity, NoCorrection free-surface response, grid/horizon bias, and fracture dissipation remain **UNKNOWN** or **NOT COMPUTED**.

The geometry figure is pending and uses equal world-unit millimetres with no displacement exaggeration. The timeline/graph/raw-energy figure is pending.
