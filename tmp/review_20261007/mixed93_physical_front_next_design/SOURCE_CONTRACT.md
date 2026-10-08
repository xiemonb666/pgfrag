# Mixed93 finite-front diagnostic: exact source and load contract

Contract ID: mixed93_old_tip_pureI_fixed_traction_zero_g.v1. Owner: source-only preparation; numerical source and run receipts belong to the full3d operator agent. Root decision is to run an uncommitted pure-LEFM diagnostic, then use cohesive-history dissipation for later real growth.

## Frozen identity and geometry

- Input case: pgfrag/experiments/mixed93_cached_lefm_solve_20261006/case.json; SHA-256 67ff4df51a0269e88b8d107a2a95929c3990a7cd9605585a31a2df583434cc1c.
- Consumer case fcc4: pgfrag/experiments/mixed93_actual_cohesive_bridge_20261006/case.json; SHA-256 fcc4b46ed101b0f6067a9f3ec10816909c9fa931edd5a8cdde4c7f46cb3f3a71. It binds the accepted array digest 3e1e6fd17751e4e06113010f20ea33c7b69914ca90431e26bbf76124951b5d41; this preparation did not open the NPZ.
- Geometry: cylinder R=L=0.01 m, z=[−0.005,0.005] m; y=0 crack plane; old front x=0 along z; candidate front x=0.0005 m. Source sites are the ordered 3×3×3 product grid with axes −0.005, 0, 0.005 m.
- Basis: 93 semantic DOFs = 81 BASE + 12 GLOBAL_FRONT. Rebuild shifted branch basis at candidate tip; keep material source set, physical domain, and external boundary fixed. No 162-DOF Mode-III operator/cache may be substituted.

## Physical load is frozen; coefficient load is not

- Source rows are the accepted mixed93 pure-I boundary traction array in Pa, generated for old-tip KI=100 Pa√m. Apply the fcc4 scale 0.01 exactly once, giving the synthetic KI=1 Pa√m boundary load. The field is a benchmark input, not calibrated material data.
- Keep the same physical boundary positions, outward unit normals, area weights, and scaled traction samples in both geometries. Row order/slices are side [0,256), bottom cap [256,512), top cap [512,768); total 768 rows, 256 per region.
- Frozen components are positions_xyz_m, boundary_normals, boundary_area_weights_m2, and pure_I_boundary_traction_Pa_scaled. Crack faces remain traction-free; they are not part of the outer-load sample array.
- For each basis, re-integrate its own coefficient vector: f(a)_j = sum_q area_weight_q * N_j(a,x_q) · traction_q. Therefore f_old and f_candidate may differ while the physical load is identical. Never copy old f into the shifted basis, recompute t from the shifted solution, retune KI, project the load, or drop components.

## Exact quadrature and query budget

- Frozen a0/theta32 selection: 512 volume rows, GL4 radial parameter × GL32 theta × GL4 z. At a=0, r=R*s² and dV=2*R²*s³ ds dtheta dz. The physical domain is the full disk times z; open nodes avoid the front without cutting a hole. Side has 256 rows from GL64(theta) × GL4(z), normal=(cos(theta),sin(theta),0), dA=R dtheta dz. Each cap uses GL4 on inner/outer radial segments split at c=0.5 mm and GL16 on each theta semicircle; both caps total 512. Inner r=c sin²(pi*s/2), outer r=c+(R-c)s²; dA=r*(dr/ds)ds dtheta. Cap normals are -ez at bottom and +ez at top. Boundary row order is side, bottom cap, top cap; 768 external-boundary rows total. There are 128 off-grid diagnostic rows. Raw selected-row uses=1408; 1280 rows enter bulk or external-boundary integration.
- Old and candidate geometries use the same integration xyz, volume weights, boundary normals, and area weights for this finite comparison. The candidate branch factor/gradient is recomputed analytically at all required rows; BASE MLS values/gradients reuse the exact frozen row IDs and 27-site support.
- Added MLS/provider queries = 0 only when cache digest, source order, MLS configuration, row IDs, coordinates, weights, and normals match exactly. Otherwise block. A changed quadrature or new query is outside this budget and requires a new reviewed contract.
- Two independent zero-gauge re-equilibrations are required: a=0 baseline and a=0.5 mm candidate. Budget two bulk K assemblies and two 99×99 KKT solves; one physical load integration per basis; zero cohesive-primitive calls; no retry.

## Gauge, energy, terminal status

- For each basis, assemble C for the same functional V^-1 ∫[ux,uy,uz,((x-center)×u)/R]dV, center=(0,0,0), R=0.01m, order tx,ty,tz,rx,ry,rz. Set g=(0,0,0,0,0,0) in both solves. Do not use cached manufactured nonzero-g q_solution as baseline.
- Verify per-state U=0.5 qᵀKq, W=qᵀf and the independent boundary integral against the frozen traction, Π=U−W, equilibrium/gauge/KKT residuals, rigid-mode count, and finite 93-field reconstruction.
- Report ΔΠ and optionally −ΔΠ/(Lδa) only as a finite-potential diagnostic. Do not compare Gc, add GcΔA, combine virgin CZM stored energy with a Griffith term, create/commit a new ExplicitSurfaceState, or call this autonomous growth.
- This phase is FIXED_LOAD_FINITE_PERTURBATION_DIAGNOSTIC_ONLY. The later calibrated cohesive-growth stage uses irreversible cohesive-history ΔDcoh; do not also add GcΔA for that same damage. Nonconverged, unstable/snapback, stale, or unbound inputs are UNKNOWN.

## Nonclaims

No material calibration, J/Gc acceptance, convergence proof from this single quadrature rule, physical growth, cohesive damage, front turning/branching, arbitrary multi-crack geometry, crack initiation, crack speed, or fragment count is delivered.


