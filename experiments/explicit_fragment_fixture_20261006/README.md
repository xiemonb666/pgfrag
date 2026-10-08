# Explicit planar surface to paired fragment fixture

This bounded experiment consumes a closed box and an explicit triangulated planar cut patch, clips the box on the two sides of the patch, and exports two independently indexed closed triangle meshes. It exists to demonstrate the export and annotation mechanics only.

Every generated artifact is `FIXTURE_GEOMETRY_ONLY`. The count of two pieces is prescribed by this one input plane. It says nothing about physical crack growth, real fragment counts, model R2, or production qualification. Mechanics are `NOT_COUPLED`; no production event, relation, or observed-feature label is emitted. The lineage label `explicit_fracture_face_fixture` means “the supplied fixture cut face,” not a measured or simulated physical fracture.

## Input and interface contract

`input_fixture.json` is the runnable explicit input. It carries a closed, outward-oriented polygon mesh for the original box and a `cut_surface` containing its stable `surface_id`, plane, mother-coordinate vertices, and oriented triangles. The input triangle winding follows the plane normal. This implementation accepts one complete planar section through the box: the patch boundary must exactly follow the full plane/box intersection. It rejects inset patches, planes that miss or only touch the box, ambiguous cuts through box vertices, nonplanar triangles, holes, nonmanifold boundaries, and partial triangle coverage.

The current adapter is `buildFixture(input)` in `explicit_fragment_fixture.js`. Its bounded planar geometry stage is `validateCutSurface(...)` followed by `clippedHalf(...)`. The stable output contract is two closed oriented meshes, per-face source lineage and labels, volume/centroid/mass, and an exact reciprocal vertex/triangle interface ledger.

To accept a supplied curved surface later, retain the same input identity and output contract but replace the planar stage with a separately qualified adapter of the form:

```text
splitClosedMeshByExplicitSurface({ originalClosedMesh, cutTriangleMesh, sourceLineage })
  -> { negativeClosedMesh, positiveClosedMesh, interfacePairLedger, metrics }
```

That adapter must perform robust surface/edge intersections, conforming subdivision of crossed original faces, side classification, and paired cap construction from the *same* authoritative input surface. It must reject unsupported self-intersections, branches, gaps, and incomplete partitions. The current fixture does not implement those curved-surface operations and does not approximate them by fitting a plane or by repairing a broken-bond graph.

## Reproduction

From the repository root in PowerShell:

```powershell
node pgfrag/experiments/explicit_fragment_fixture_20261006/explicit_fragment_fixture.js
node pgfrag/experiments/explicit_fragment_fixture_20261006/test_explicit_fragment_fixture.js
```

The implementation and the one bounded verification use only Node.js built-ins. The fixture output is written to `pgfrag/tmp/review_20261006/explicit_fragment_fixture_20261006/`.

## Checks and deliverables

The splitter verifies oriented closed manifold topology, connectedness, Euler characteristic 2, positive piece volumes and centroid-side signs, volume and mass conservation, non-overlapping half-space interiors, identical interface mother coordinates, reciprocal piece-local vertex maps, and opposite outward interface normals. It enforces 64 total declared vertices and 128 total declared faces across the box, input cut patch, and both outputs. `checks.json` records the checks and residuals.

The review directory contains the full mesh/lineage JSON, separate OBJ meshes, a standalone strict interface ledger, the check report, and `exploded_mesh_preview.png`. The PNG is rendered from the exported mesh arrays; its small visual gap is only an illustration transform and does not alter mesh coordinates or volume.
