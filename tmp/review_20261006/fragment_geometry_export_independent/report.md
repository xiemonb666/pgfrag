# Independent fragment geometry review — 2026-10-06

## Scope and method

Read the fixture inputs, producer source, saved JSON/OBJ/ledger/checks, and the current fracture-artifact contract source/schema/README/tests. No solver or experiment/contract test was run. The only execution was a bounded independent Node.js readback against the saved small fixture; it verifies output geometry and hashes its source inputs. Its receipt is [independent_readback.json](/C:/Users/xiemo/Desktop/ShardForge/pgfrag/tmp/review_20261006/fragment_geometry_export_independent/independent_readback.json), with the reproducible audit script beside it.

The readback script initially stopped before writing a receipt because my audit code looked for a computed area on the raw fixture input. I corrected that local audit-script reference to the generated output, then the readback passed in under one second. No producer artifacts or source were changed.

## Geometry result

For this explicit box split, I found no geometry/export defect:

- The supplied cut plane strictly crosses the 8-vertex box. Independently reconstructed edge intersections equal the four declared cut vertices, and all four cut-boundary edges lie on the cube exterior.
- Both exported pieces are closed, consistently wound, connected genus-zero triangle meshes: each has V=8, E=18, F=12, Euler characteristic 2. Independently recomputed signed volumes are 4.44 m³ and 3.56 m³, both positive, totaling the parent 8 m³. The computed centroids match the saved values.
- At the supplied density of 2400 kg/m³, fragment masses sum to the parent 19,200 kg. Original-face lineage covers the parent exterior area of 24 m²; each side’s clipped per-face areas match its source-face lineage.
- Both sides contain the same two input cut triangles in mother coordinates. Per-triangle area is 2.1774296773948865 m² on both sides; paired geometric normals have dot products -1 and -0.9999999999999999. Total paired interface area equals the input patch area, 4.354859354789773 m².
- OBJ vertices/faces exactly match the JSON arrays, and OBJ groups agree with original-face versus explicit-fixture-cut lineage. Triangle indices in the pair ledger point to those exported faces.

These results verify the small deterministic geometry fixture only. The source says the plane is an explicit input and tags its cut faces as explicit_fracture_face_fixture; it also emits FIXTURE_GEOMETRY_ONLY, NOT_COUPLED, and NOT_CLAIMED for physical fracture. See [fixture README](/C:/Users/xiemo/Desktop/ShardForge/pgfrag/experiments/explicit_fragment_fixture_20261006/README.md:5), [clippedHalf and createPairLedger](/C:/Users/xiemo/Desktop/ShardForge/pgfrag/experiments/explicit_fragment_fixture_20261006/explicit_fragment_fixture.js:330), and the saved independent values above.

## Contract review boundary

The current contract is appropriately strict about declared metadata: graph components cannot define fragments; a fragment must reference a declared closed mesh; the observation projection returns only the observation channel while clean labels require the explicit audit path; and fixture/simulation classes cannot qualify for pretraining. The latest contract also includes per-side surface triangle rows, reciprocal reversed triangle connectivity, and agreement between computed triangle area and reference_area_m2. See [validator](/C:/Users/xiemo/Desktop/ShardForge/pgfrag/experiments/fracture_artifact_contract_20261006/fracture_artifact_contract.py:371), [pair validation](/C:/Users/xiemo/Desktop/ShardForge/pgfrag/experiments/fracture_artifact_contract_20261006/fracture_artifact_contract.py:474), and [A12 read projections](/C:/Users/xiemo/Desktop/ShardForge/pgfrag/experiments/fracture_artifact_contract_20261006/fracture_artifact_contract.py:737).

The remaining boundary is material for integration: closed_mesh is currently a hash, format, frame, watertight boolean, and closure-method string. The validator does not load that mesh or prove its hash matches the supplied bytes. Pair triangles are checked against their paired surface-coordinate arrays, but the artifact has no references from those triangles to face IDs in either fragment mesh. It therefore cannot establish that the pair is actually a boundary of the two declared fragments, that the fragment partition covers the mother, or that fragment mass/volume agrees with mesh geometry. The README calls this out correctly: watertight=true is source-provided closure evidence, and FIXTURE_ASSERTION_ONLY is not geometry QA. Keep that distinction explicit in downstream qualification.

Also, the validator checks per-vertex normals for unit length and opposing pair values, but does not compare them with the oriented geometric normals computed from the triangle winding. Triangle rows are reciprocal and their areas are measured; the current checks still do not require each surface patch to be a manifold, connected, or free of zero-area triangles. These are concrete surface-integrity gaps, separate from physical fracture truth.

## Evidence boundary and next decision

Neither this fixture nor the artifact contract recovers a crack from native PG-Frag bond state. The existing PG-Frag README says the bar runner emits point/bond state but no geometric crack surface or fragment shapes, and that forge_asset_surface is unsupported. The plane fixture has zero information about material-generated roughness or lattice/discretization directional bias; its cut is flat because it is prescribed to be flat. Do not map its two pieces or cut labels into physical fracture GT.

For A12 integration, choose whether the artifact is a metadata envelope or a geometry-verified partition. If the goal is actual paired fracture surfaces attached to ranked candidates, the next contract revision should tie the supplied fragment mesh bytes to a geometry validator receipt and reference each paired surface triangle to the corresponding face(s) in each mesh. That validator should compute, rather than trust, closure/orientation/manifold status, paired-coordinate and normal agreement, interface-area equality, and parent partition mass/volume conservation. Only after that interface exists can a surface reconstructed from solver state be evaluated; uncertainty about unresolved sub-point roughness and orientation bias must remain explicit. Until then, accept the current schema as an export/anti-leak boundary, not as evidence that surfaces or fragments are physically identified.

The separate actual-run contract validation was left to its author; this review inspected the current contract source and fixture definitions but did not execute its tests.

## RelicWeave/VesselForge handoff

The current A12 observation projection is not yet a consumable inputfragment-to-ranked-candidates adapter. Observation rows contain opaque observed-mesh and visibility-mask hashes plus sample/visible-area metadata; they do not embed or resolve those mesh bytes, and they do not link observation IDs to producer fragment IDs or provide candidate rows/ranks. That is a sound leakage boundary for inference, but a downstream adapter still needs an authorized way to load the observation geometry and produce ranked candidates/evidence without joining clean GT. No such runtime handoff was exercised in this review.
