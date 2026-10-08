# PG-Frag static bond visibility classifier: source and fixture review

**Review state: v2 source/test draft; execution frozen pending root review.** This work adds an isolated static classifier and its tests. It does not call the native solver or change an existing module. The v1 source, tests, and review files were archived byte-for-byte before v2 edits. No pytest, CAD construction, STEP export, Julia, solver, or Git command has been run for v2.

## Source basis

- `../cavity_bond_safety/design.md` SHA-256 (recomputed from the current read-only file): `6035fa7c3fdead24a43a61ff653b03f55f002326ddebe95b5228db0d0b245e0b`.
- `vesselforge/src/vesselforge/wall_reference.py` SHA-256: `a4f06671fc138fcf6912d9d075841fa958185b44a1d9c7e93b0b62e2fe29e8c8` (read-only; no changes).
- The frozen design selects canonical STEP/B-rep all-hit intersections plus exact open-interval solid occupancy before any native `BondSystem` initialization. It requires reciprocal undirected decisions, explicit source hashes/units/frame/tolerance, and separate static topology records for OUT/VOID and UNKNOWN.
- `wall_reference.py` binds its stations to a canonical STEP hash, requires declared mm coordinates and `right_handed_z_up`, checks the source semantics SHA, and explicitly leaves mesh binding and physics unassessed. This classifier follows the same source-hash/unit/frame discipline and makes the mm-to-metre conversion explicit; it does not import or mutate that frozen module.

## Added files and behavior

- `pgfrag/src/pgfrag/bond_visibility.py` defines `load_canonical_step`, `classify_candidate_pairs`, source records, and the static ledger. STEP bytes are hashed before import; mismatches raise. CadQuery is asked to normalize imported STEP geometry to mm. The source must attest mm, a frame, and a source accuracy basis. One closed solid may use its own domain ID; multiple solids require an explicit non-empty material-domain ID for each solid. Open, invalid, nonmanifold, mixed leaf, or precision-unknown sources remain UNKNOWN.
- Candidate pairs are canonicalized by global point IDs. For every pair, the classifier preserves all OCC line/face hits, their face/solid references and parameters, clusters potential ambiguities without discarding raw rows, and classifies the open interval between every distinct event. Any void interval gives OUT_OR_VOID/AIR_PATH. Tangency, ON, ambiguous hits, overlap, unresolved faces, or intervals no wider than twice the combined uncertainty band give UNKNOWN.
- V2 closes a fail-open path in single-link hit clustering: any cluster with more than one distinct hit parameter is unresolved, regardless of its total span or face identity. Its full finite span, raw hit rows, and reason remain in the ledger; the representative midpoint is diagnostic-only and can never make a pair eligible. Exact-parameter duplicates remain subject to face-identity ambiguity checks.
- V2 maps arithmetic outside the finite representable range to explicit UNKNOWN reasons. Segment vectors are checked before `math.hypot`; exact coincident endpoints alone are `ZERO_LENGTH_SEGMENT`. The combined accuracy band and its doubled interval threshold must both be finite. Interval/hit midpoints and coordinates use bounded arithmetic, while finite raw endpoint and hit parameters are retained. Tests serialize these negative controls with JSON `allow_nan=False`.
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
| Pure-math cluster chain | Same-face parameters `[0, 0.75 tol, 1.5 tol]` | Keep all three raw hits; unresolved diagnostic cluster; UNKNOWN, never eligible |
| Near-distinct short-cavity control | Same-face parameters `10` and `10 + 0.25 tol` | Preserve the sub-tolerance span and raw hits; force UNKNOWN even if mocked diagnostic interval midpoints classify IN |
| Numeric-range controls | Opposite extreme finite endpoints; overflowing accuracy+tolerance or doubled band; large finite hit parameters | Structured UNKNOWN; no fake zero length, no nonfinite JSON, finite raw values retained |

The candidate coordinates are fixed undeformed reference coordinates. The fixture dimensions do not select or tune `δ`, point spacing, process-zone length, or any material parameter.

## Freeze boundary and next allowed work

Before root review, the only permitted work is source/test/design text review and edits within the three owned paths: `bond_visibility.py`, `test_bond_visibility.py`, and `pgfrag/tmp/review_20261004/bond_visibility/`. Do not run the check below until root explicitly reviews and authorizes it. The proposed next checkpoint is one bounded pytest invocation covering the pure-math class plus only the hollow-shell all-hits/air-path and solid-coupon controls. From PowerShell, with the repository root as the working directory:

```powershell
$env:PYTHONPATH = "pgfrag/src"
& "third_party/vesselforge-runtime/venv/Scripts/python.exe" -m pytest `
  pgfrag/tests/test_bond_visibility.py::BondVisibilityPureMathRegressionTests `
  pgfrag/tests/test_bond_visibility.py::BondVisibilityGeometryTests::test_hollow_shell_endpoints_inside_but_segment_crosses_void `
  pgfrag/tests/test_bond_visibility.py::BondVisibilityGeometryTests::test_hollow_shell_preserves_all_outer_and_inner_hits `
  pgfrag/tests/test_bond_visibility.py::BondVisibilityGeometryTests::test_same_wall_and_solid_coupon_remain_inside_material
```

This command is a proposed plan only; no test was run. It uses in-memory analytic shapes and exports no STEP fixture. Expand to the handle/tangent/thin controls only after reviewing those results.

If the later CAD pass finds an API or topology issue, preserve the fixture reports and stop before native integration. No post hoc geometry-mask result may be described as a fracture event, damage, fragment, R2, energy, or physical PASS. `δ`, `dx`, process-zone length, and constitutive parameters remain independent inputs.
