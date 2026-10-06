# Mixed global-front elastic operator (source only)

This experiment assembles raw three-dimensional elastic stiffness and surface-traction force for a semantic 93-DOF mixed basis. With 27 continuous MLS sources it has 27 scalar `SOURCE_BASE` functions times three displacement components (81 coefficients), plus four independent global tip-branch functions times three components (12 coefficients). The total is 31 scalar functions and 93 vector coefficients. `GLOBAL_FRONT` columns are global amplitudes and are never represented as per-source enrichments.

`adapt_mixed_shape_row` consumes the reviewed four-branch helper row. It extracts each source's base value/gradient once and reads the four branch values/gradients once, then emits stable, scope-tagged IDs. Source order is canonical lexical ID order; component order is `ux, uy, uz`; global branches are indexed 1 through 4. IDs encode `SOURCE_BASE` versus `GLOBAL_FRONT`, so a global branch cannot be mistaken for one of 27 source-owned columns.

For transparency, the adapter also supports the separately named `per_source_front_diagnostic` layout. This reproduces the old 27 × (BASE + four branches) × three-component 405-column row and may assemble its raw matrix. It does not drop columns to fit another cap or certify that layout as mechanically eligible. Under linear reproduction the four global functions satisfy

\[
xF_3-yF_4+yF_1=0,\qquad xF_4+yF_3-yF_2=0.
\]

In a per-source partition of unity, these yield two zero-field coefficient vectors per displacement component (six additional null candidates) beyond the six physical rigid modes. The operator exposes those candidates for diagnostics and applies no truncation, ridge, pseudoinverse, load projection, KKT solve, or eigenvalue filtering. The global 93 layout removes these source-wise coefficient dependencies, but it is not yet physically qualified: that requires an authorized spectrum check demonstrating exactly six physical zeros for the actual 93-DOF matrix.

The constitutive matrix is full isotropic 3-D `C(E, nu)` in `[eps_xx, eps_yy, eps_zz, gamma_xy, gamma_yz, gamma_xz]` engineering-strain order. With `B` in 1/m, `C` in Pa, and positive volume weights in m³, raw `K = sum(w B.T C B)` has units N/m. Surface force uses `sum(w_m2 Nvec.T traction_pa)` and has units N. Each input chunk has at most 256 rows; the bounded API accepts at most 256 chunks. Matrix and vector payloads are immutable little-endian float64 `bytes`; K is returned as accumulated, without forced symmetrization.

Two existing 375-DOF checks apply to the frozen full-3D route. The registry/problem-budget guard is `MAX_TOTAL_DOF` in `pgfrag/experiments/enriched_dof_contract_20261005/enriched_dof_contract.py::validate_problem_budget`; the reused array assembler independently defines `MAX_DOF` in `pgfrag/experiments/efgm_elastic_operator_20261005/efgm_elastic_operator.py` and rejects an oversized system before allocating K. The full-3D wrapper calls both paths. They are explicit source-level bounds, not a limit implied by the B/C/K equations. Neither check nor any frozen operator is changed here. The 405 layout is confined to this separate explicitly diagnostic implementation and does not alter either existing cap.

This source-only package has no MLS provider, physical field generator, stiffness run, force run, solve, crack-growth model, or fracture calibration. The small review checker is prepared but not executed; root authorization and a separate 93-DOF spectrum gate remain prerequisites for any mechanical benchmark claim.
