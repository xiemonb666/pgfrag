# Native cohesive fragment extraction

`pgfrag.cohesive_fragments.extract_cohesive_fragments` reconstructs bulk
mechanical components and their boundary triangles directly from native T4 and
COH3D6 connectivity. It does not use volume labels or infer crack faces from a
plane or from coordinates.

## Inputs and failure rule

The required arrays are `rest_nodes (N, 3)`, `native_tets (T, 4)`,
`native_cohesive (F, 6)`, `facet_damage_qp (F, Q)`, and
`tractions_qp (F, Q, 3)`. The caller must declare `damage_tolerance` and
`traction_tolerance`. A facet is removed from the mechanical graph only when
every quadrature point satisfies both inclusive checks:

```text
damage >= 1 - damage_tolerance
norm(traction) <= traction_tolerance
```

The result keeps the per-QP damage, traction norms, pass bits, and thresholds
so the split decision can be audited. An ordinary native face with two T4
parents joins those elements. Each COH3D6 side must match one exact native T4
face with exactly one parent; the two sides must have distinct node IDs and
distinct parents. An unfailed cohesive facet joins its two parent tetrahedra.
All T4 rows are partitioned into components and their volumes are summed and
checked against the full native bulk volume.

## Surface and lineage rules

For every boundary triangle the extractor preserves exact native node IDs,
parent T4 row, local face ordinal, outward orientation, rest area, and current
area when `current_nodes` is supplied. A failed facet produces a paired record
for both side faces, with exact native node triples, both parent incidences,
the optional parent facet/source labels, and per-side plus mean area.

Unfailed cohesive side triangles are internal and omitted from the component
surface. Their corresponding nodes are joined only for the boundary edge audit,
which closes the shell across the duplicated native IDs. With
`native_node_to_source`, the extractor verifies that each side has the same
three source nodes and pairs nodes by source identity. Without that lineage, it
uses the COH3D6 ordinal convention (side node 1 pairs with side node 4, and so
on). That convention should only be used when the native element ordering is
known to follow the COH3D6 pairing convention. Either pairing route must also
prove exact rest-coordinate equality for every paired node; the paired rest
side areas must then match exactly. An ancestry match alone is insufficient
when the native rest geometry contains a gap or offset.

Any other singly incident T4 face is labeled `original_exterior`; it is never
silently relabeled as a crack. Such faces carry
`source_exterior_ancestry_verified: false` unless the caller supplies
`source_exterior_faces`. For verification, provide a mapping of source face ID
to its three source node IDs, plus `native_node_to_source`. The source exterior
list must match the non-cohesive native boundary faces exactly, with no missing
or extra ancestry.

The native edge audit checks that every actual boundary edge appears twice
with opposite orientation using exact native node IDs. This audit drives
`closed_solid` and qualification. A second diagnostic,
`intact_cohesive_welded_edge_audit`, shows the closure obtained by virtually
pairing nodes across unfailed cohesive facets; it does not turn a spring-linked
surface into a geometrically closed native T4 solid. Components with active
intact cohesive facets receive the explicit reason
`active_cohesive_interface_not_geometrically_closed`.

`qualified_closed_mechanical_fragment` additionally requires verified source
exterior face ancestry and no exact native node shared with another component.
If source exterior ancestry is omitted, faces remain labeled
`original_exterior`, but qualification is false with
`source_exterior_ancestry_unverified`. Components that share a node without
sharing a face remain separate in the graph and are marked ineligible. Euler
characteristic is reported from the native boundary audit; it is descriptive
and does not impose a shape or genus assumption.

T4 connectivity is required to have positive signed volume in `rest_nodes`;
when `current_nodes` is supplied, each deformed T4 must also have positive
signed volume. Inverted or zero-volume cells raise `ValueError`. `rho`, when
provided, is a nonnegative bulk density and yields rest bulk mass. Current
coordinates and per-component current bulk volume are retained when supplied.

The qualification is limited to native connectivity and topological closure.
It does not establish material calibration or an R2 result.

## Example

```python
from pgfrag.cohesive_fragments import extract_cohesive_fragments

result = extract_cohesive_fragments(
    rest_nodes,
    native_tets,
    native_cohesive,
    facet_damage_qp,
    tractions_qp,
    damage_tolerance=1e-8,
    traction_tolerance=1e-6,
    current_nodes=current_nodes,
    native_node_to_source=native_node_to_source,
    parent_facet_ids=parent_facet_ids,
    parent_source_ids=parent_source_ids,
    source_exterior_faces=source_exterior_faces,
    rho=rho,
)
```

Each entry in `result["components"]` includes `tet_rows`, `native_node_ids`,
`face_indices`, `boundary_faces`, `oriented_native_triangles`, the edge/Euler
audit, rest/current bulk volumes, optional mass, and explicit qualification
status/reasons. `result["failed_interface_face_pairs"]` links each failed
cohesive facet to the two generated boundary face indices.

## Focused validation

From the repository root in PowerShell:

```powershell
$env:PYTHONPATH = (Resolve-Path 'pgfrag/src').Path
& 'third_party/vesselforge-runtime/venv/Scripts/python.exe' -m pytest pgfrag/tests/test_cohesive_fragments.py -q
```

The focused suite covers intact and fully failed interfaces, one nonpassing QP,
shared native nodes across components, complete rest/current volume and mass
accounting, missing exterior ancestry, incorrect cohesive pairing, exact paired
rest coordinates and areas, and nonpositive rest/current tetrahedron volume.
