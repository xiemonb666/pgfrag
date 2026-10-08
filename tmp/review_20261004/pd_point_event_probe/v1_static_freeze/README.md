# Native PD point-event probe

This isolated experiment keeps the Peridynamics.jl 0.5.4 native `OSBMaterial` and `VelocityVerlet` path. It does not change the frozen PD baseline, the SDK runner, or vendored source. The observer is attached through experiment-local, more-specific methods for the native `export_reference_results` and `export_results` generics. Native `submit(job; quiet=true)` still builds and initializes the data handler, chooses the accepted time step, applies boundary conditions, computes force density, and advances every step.

The current fixture is a solid 5 mm × 5 mm × 5 mm coupon represented by 125 1 mm³ material cells. Its primary input uses cell-center uniform voxel quadrature. A separate deterministic candidate jitters each point within its original cell while retaining that cell's positive 1 mm³ quadrature volume. The jittered candidate is explicitly a declared quadrature experiment; its volumes are not exact Voronoi volumes and it is not a reconstructed shell or measured specimen.

## Frozen point-input contract

Each input uses a UTF-8 CSV with the exact header `point_id,x_m,y_m,z_m,volume_m3`. Row order defines the native Body index; `point_id` is the preserved source identity and must be a unique positive integer. Coordinates are rest-space SI metres and every volume is a positive SI cubic metre value. The required sibling TOML declares the schema, point-source kind/description/reference/hash, coordinate frame and units, point-ID semantics, quadrature method and volume semantics, point count, volume sum, volume-weighted center of mass, and explicit top/bottom source point IDs for a future axial case. CSV bytes are pinned in the experiment manifest and rechecked by the launcher before execution.

An external frozen point set can replace the selected CSV and metadata pair after its source, units, IDs, quadrature meaning, file digest, volume sum, center of mass, and any driven point sets are recorded in the same contract. The probe maps the input row order to native indices and carries the supplied point IDs into all state and bond outputs. It does not infer ancestry from coordinates.

## Output and event semantics

The step-zero reference hook records state after native initialization. Each subsequent result hook records a completed native VelocityVerlet step at `time_s = step × accepted dt`. The observer stores point reference/current positions, displacement, velocity, acceleration, native internal/external force densities (`b_int`, `b_ext` in N/m³), damage as a descriptive field, positive point volume, and mass. It writes every directed native bond state and the corresponding canonical undirected pair state at every sample, then writes directed and canonical state transitions from the immediately previous native sample. Reciprocal directed entries must exist and agree at every sample or the run stops.

Files named `canonical_bond_*` describe the native intact-bond graph only. They do not certify physical fragments, crack surfaces, contact, energy conservation, verified lineage, R2, or training eligibility. Missing energy quantities remain `UNKNOWN`; no missing energy term is filled with zero. The parameter values in `common_case.toml` are `UNVALIDATED_PROTOTYPE` values.

## Runtime gates

The frozen launcher targets only `case_zero_load_control.toml`, one native step, one CPU thread, and at most 512 points. Its working directory is `pgfrag/tmp/review_20261004/pd_point_event_probe/OUT`, which must be empty. Snapshots and process receipts are written to sibling directories named `snapshots` and `receipts`. The launcher records the actual numeric Julia process exit code and refuses to overwrite any prior run directory or receipt.

The zero-load runtime gate is prepared but still awaits final freeze review. It must show exactly step 0 and step 1, one native chunk, complete stable point-ID mapping, positive mass closure, correct accepted native dt, complete reciprocal directed/canonical graph state at both samples, and no activity transitions. Do not launch Julia, the solver, or a preflight until the root agent explicitly releases that freeze.

`case_tinyload_event_diagnostic.toml` is prepared as a later diagnostic configuration only. It is not reachable through the frozen launcher and requires a separate runtime authorization. Neither case assumes a target fragment count.
