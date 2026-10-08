# PG-Frag static bond visibility classifier: source and fixture review

**Review state: implementation drafted; execution frozen pending root review.** This work adds an isolated static classifier and its tests. It does not call the native solver or change an existing module. No pytest, CAD construction, STEP export, Julia, solver, or Git command was run for this review package.

## Source basis

- `../cavity_bond_safety/design.md` SHA-256: `6035dfa7c3fdeade43a61ff653b03f55f002326ddebe95b5228db0d0b245e0b`.
- `vesselforge/src/vesselforge/wall_reference.py` SHA-256: `a4f06671fc138fcf6912d9d075841fa958185b44a1d9c7e93b0b62e2fe29e8c8` (read-only; no changes).
- The frozen design selects canonical STEP/B-rep all-hit intersections plus exact open-interval solid occupancy before any native `BondSystem` initialization. It requires reciprocal undirected decisions, explicit source hashes/units/frame/tolerance, and separate static topology records for OUT/VOID and UNKNOWN.
- `wall_reference.py` binds its stations to a canonical STEP hash, requires declared mm coordinates and `right_handed_z_up`, checks the source semantics SHA, and explicitly leaves mesh binding and physics unassessed. This classifier follows the same source-hash/unit/frame discipline and makes the mm-to-metre conversion explicit; it does not import or mutate that frozen module.

## Added files and behavior

- `pgfrag/src/pgfrag/bond_visibility.py` defines `load_canonical_step`, `classify_candidate_pairs`, source records, and the static ledger. STEP bytes are hashed before import; mismatches raise. CadQuery is asked to normalize imported STEP geometry to mm. The source must attest mm, a frame, and a source accuracy basis. One closed solid may use its own domain ID; multiple solids require unique explicit domain IDs. Open, invalid, nonmanifold, mixed leaf, or precision-unknown sources remain UNKNOWN.
- Candidate pairs are canonicalized by global point IDs. For every pair, the classifier preserves all OCC line/face hits, their face/solid references and parameters, clusters potential ambiguities without discarding raw rows, and classifies the open interval between every distinct event. Any void interval gives OUT_OR_VOID/AIR_PATH. Tangency, ON, ambiguous hits, overlap, unresolved faces, or intervals no wider than twice the combined uncertainty band give UNKNOWN.
- Pair records include both directed orientations with the same decision. The initial topology ledger counts IN_MATERIAL, OUT_OR_VOID, and UNKNOWN separately and states that no native graph, neighbor count, physics, or whole-asset qualification was produced. Only IN_MATERIAL is marked as a static eligible candidate.
- `pgfrag/tests/test_bond_visibility.py` contains unrun in-memory fixture tests and a hash-mismatch check that rejects before attempting STEP import. The fixtures are code only; they have not been constructed or exported.

## Fixture plan and expected controls

| Fixture | Candidate path | Expected decision |
|---|---|---|
| Hollow box shell: outer 20 mm cube minus inner 12 mm cube | Points `(-8,0,0)` to `(8,0,0)`, both inside opposite walls | Preserve two inner-wall hits; intervals IN / OUT / IN; `OUT_OR_VOID`, `AIR_PATH` |
| Same hollow shell, full section | `(-11,0,0)` to `(11,0,0)` | Preserve all four outer/inner B-rep crossings at about 1, 5, 17, and 21 mm; interval occupancy OUT / IN / OUT / IN / OUT |
| Same wall and solid coupon | Two points whose exact segment remains inside one validated solid | `IN_MATERIAL` |
| Handle-gap control | Two separated analytic box arms, each assigned its own explicit solid/domain | `OUT_OR_VOID`; preserve void interval and distinct-domain reason |
| Boundary controls | Sphere tangent path; shell points on the inner face | `UNKNOWN` for tangent/non-crossing or ON classification |
| Thin-material control | 0.01 mm slab with a declared accuracy band larger than the slab | `UNKNOWN`, retaining the short interval and its reason |
| Source qualification controls | Non-solid wire; otherwise valid solid with missing source accuracy | `UNKNOWN`; no whole-asset physics status other than `NOT_ASSESSED` |
| Hash control | Temporary non-STEP bytes with a deliberately wrong expected SHA | Reject on actual byte hash before CadQuery/STEP parsing |

The candidate coordinates are fixed undeformed reference coordinates. The fixture dimensions do not select or tune `δ`, point spacing, process-zone length, or any material parameter.

## Freeze boundary and next allowed work

Before root review, the only permitted work is source/test/design text review and edits within the three owned paths: `bond_visibility.py`, `test_bond_visibility.py`, and `pgfrag/tmp/review_20261004/bond_visibility/`. The next execution checkpoint should authorize at most one minimal CAD/pytest pass against the hollow-shell and solid-coupon controls. Do not export a STEP fixture; construct the small analytic shapes in memory in that later test pass. Expand to the handle/tangent/thin controls only after reviewing the first results.

If the later CAD pass finds an API or topology issue, preserve the fixture reports and stop before native integration. No post hoc geometry-mask result may be described as a fracture event, damage, fragment, R2, energy, or physical PASS. `δ`, `dx`, process-zone length, and constitutive parameters remain independent inputs.
