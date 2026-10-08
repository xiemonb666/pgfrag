# Moving-front 93-DOF cache mapping

Status: named cache rows and the independent 27-source BASE MLS table are verified. The cached cap rule is not bound to the pending dual-tip cap contract, and no mechanical solve or moving-front result is claimed.

## Audited inputs

The read-only inspection verified these hashes:

- tip_aligned_quadrature_study_20261006/run01/arrays.npz: a9de2e3eaf88ba0eb130fca6f160d5d327a963bb00203266b7390cadbe6693ea
- Its result.json: 2510eabd3e3406f243d9fe900ac5af3cfdf17071cdbe50383093de0a13adff4b
- mixed93_cached_lefm_solve_20261006/run01/arrays.npz: 3e1e6fd17751e4e06113010f20ea33c7b69914ca90431e26bbf76124951b5d41
- mixed93_actual_cohesive_bridge_20261006/case.json: fcc4b46ed101b0f6067a9f3ec10816909c9fa931edd5a8cdde4c7f46cb3f3a71

The 93-array file came from a producer that failed during result serialization. Its independent cache audit is RECOVERED_CACHE_INDEPENDENTLY_CONSISTENT_PHYSICS_UNQUALIFIED; it is metadata/cache evidence, not a successful solve receipt.

## Source and BASE MLS binding

Both NPZ files contain the same exact 27×3 source-position array in the same order. Its float64 C-order SHA-256 is 997d598f8672b307437f5dfeb68fa47d0d53859e873b077ab7cc90ffab93e9f3. The canonical IDs are node-000 through node-026, matching the frozen result and fcc4 case. They enumerate the Cartesian product of x, y, z coordinates [-0.005, 0, 0.005] m, with z varying fastest.

The tip-study NPZ has 4,224 unique query positions and raw BASE MLS arrays:

- base_mls_values: float64 [4224, 27], dimensionless.
- base_mls_gradients_per_m: float64 [4224, 27, 3], 1/m.

Only their NPY headers were read; the payloads were not materialized. Each saved cloud has its own unique-row IDs. For a selected cloud, retrieve raw BASE weights and gradients by indexing those row IDs into the two raw arrays, then evaluate the current 93-DOF GLOBAL_FRONT branches at the cloud coordinates. Do not split 54-column enriched tables. The 93-DOF archive stores an old-front 93-column basis representation, not a substitute for this independent BASE table.

The saved tip-study provider record reports 4,224 successful rows, support count 27 for every row, and finite condition numbers from 6.686204652 to 12.623216411. The frozen provider inputs in the study driver are reference origin (0,0,0) m, reference scale R=0.01 m, support radius 0.025 m, max condition 1e12, rank tolerance 1e-12, switch tolerance 1e-10, and duplicate tolerance 1e-12. These values apply only when source positions and query coordinates match the frozen arrays exactly.

## Verified row map

All XYZ rows were compared exactly with unique_query_positions_xyz_m[unique_row_ids]. Units come from the frozen result: coordinates in m, volume weights in m³, boundary area weights in m².

| Cloud | XYZ / row-ID / weight members | Rows | Unique-row-ID range | Exact row lookup |
|---|---|---:|---:|---|
| Old-front body, a=0, theta32 | cloud_a0_volume32_* | 512 | 0–511 | Pass |
| Candidate body, a=0.0005 m, theta32 | cloud_a0_0005_volume32_* | 512 | 2304–2815 | Pass |
| Physical side | cloud_side_* | 256 | 3840–4095 | Pass |
| Old-front cap16 | cloud_a0_cap16_* | 512 | 1792–2303 | Pass |
| Candidate-tagged cap16 | cloud_a0_0005_cap16_* | 512 | 1792–2303 | Pass |

The failed 93 cache’s 512-row old-front volume XYZ, row IDs, and m³ weights match the tip-study a=0 volume exactly. Its boundary’s first 256 rows match the side XYZ, row IDs, and m² weights exactly. The source IDs and source positions also match exactly. The 93 cache has 93 semantic IDs, but it contains no successful moving-front solve.

Therefore the saved 27-column BASE table can supply MLS values and gradients at the listed candidate body rows without new MLS queries, provided the geometry and provider pins above remain identical. The table does not establish that theta32 is sufficiently converged. Other selected coordinates require a separately bounded query decision; row IDs cannot be transferred to a different point set.

## Load and cap limits

The fcc4 bridge contract binds a fixed old-tip pure-I load with KI=1, KII=KIII=0, scaled by 0.01 from the saved KI=100 cache. It uses the 27 sources and 93-entry registry (81 BASE plus 12 GLOBAL_FRONT). For the old and candidate geometry, keep the physical surface traction field fixed and re-integrate its generalized force in each basis. Equality of physical traction does not imply equality of the two coefficient-force vectors. The ModeIII tip-study tractions and forces are not the fcc4 ModeI load and must not be reused.

The inspected cap16 point cloud is shared by both front labels, with the same 512 coordinates and row IDs. Its frozen construction is centered at the old physical tip and splits old-tip radius at c=0.0005 m; the shifted basis tip is at the split radius and angular seam. This is not the operator’s pending two-tip Voronoi/union-chart cap rule and no cap sensitivity pass was inspected. Keep cap binding pending the root decision or use a newly specified dual-tip rule. The audit also did not load cap/side normals or traction values.

The 512-row candidate body cloud is a reusable coordinate/weight set, not evidence of energy release, front acceptance, material calibration, or a committed geometry change. The 162-DOF study is synthetic ModeIII with E=1 MPa; only its source coordinates, raw BASE MLS table, and exactly indexed geometric rows are candidates for reuse. Recompute the 93-DOF ModeI branches and reanalysis under the fixed old-tip physical load. This stage remains a finite-front diagnostic with no Gc, crack growth, or physical qualification.

## Inspection receipt

The successful named-array/header inspection and preserved failed preflight are documented in receipt.json; captured argv, stdout, stderr, and exit status are under rawlogs/.
