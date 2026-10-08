"use strict";

const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const { buildFixture } = require("./explicit_fragment_fixture.js");

const inputPath = path.join(__dirname, "input_fixture.json");
const input = JSON.parse(fs.readFileSync(inputPath, "utf8"));
const result = buildFixture(input);

assert.equal(result.qualification, "FIXTURE_GEOMETRY_ONLY");
assert.equal(result.mechanics_qualification, "NOT_COUPLED");
assert.equal(result.physical_fracture_claim, "NOT_CLAIMED");
assert.equal(result.relation_labels, "NONE");
assert.equal(result.pieces.length, 2);
assert.equal(result.checks.closed_outward_manifold_meshes, true);
assert.equal(result.checks.piece_euler_characteristic_2, true);
assert.equal(result.checks.volume_sum_matches_original, true);
assert.equal(result.checks.mass_sum_matches_original, true);
assert.equal(result.checks.pair_vertex_map_reciprocal, true);
assert.equal(result.checks.interiors_disjoint_with_positive_thickness, true);
assert.ok(result.checks.bounded_geometry_caps_respected.total_declared_vertices <= 64);
assert.ok(result.checks.bounded_geometry_caps_respected.total_declared_faces <= 128);

const [negative, positive] = result.pieces;
assert.ok(negative.metrics.volume_m3 > 0 && positive.metrics.volume_m3 > 0);
assert.ok(Math.abs(negative.metrics.volume_m3 + positive.metrics.volume_m3 - result.original_body.volume_m3) < 1e-12);
assert.equal(negative.vertex_ids.some((id) => positive.vertex_ids.includes(id)), false);
for (let i = 0; i < result.interface_pair_ledger.vertex_pairs.length; i += 1) {
  const pair = result.interface_pair_ledger.vertex_pairs[i];
  const reverse = result.interface_pair_ledger.reciprocal_vertex_pairs[i];
  assert.deepEqual(reverse, {
    negative_vertex_id: pair.positive_vertex_id,
    positive_vertex_id: pair.negative_vertex_id,
    mother_coordinate_m: pair.mother_coordinate_m,
  });
  assert.deepEqual(negative.vertices_m[pair.negative_local_vertex_index], positive.vertices_m[pair.positive_local_vertex_index]);
}
for (const pair of result.interface_pair_ledger.triangle_pairs) {
  assert.ok(pair.outward_normal_dot < -0.999999);
  assert.equal(pair.mother_coordinates_m.length, 3);
}
for (const mesh of result.pieces) {
  assert.equal(mesh.triangles.length, mesh.face_lineage.length);
  assert.ok(mesh.face_lineage.every((label) => label.fixture_geometry_only === true));
  assert.ok(mesh.face_lineage.some((label) => label.kind === "original_surface_face"));
  assert.ok(mesh.face_lineage.some((label) => label.kind === "explicit_fracture_face_fixture"));
}

const insetSurface = JSON.parse(JSON.stringify(input));
const center = insetSurface.cut_surface.vertices_m.reduce((sum, p) => sum.map((v, axis) => v + p[axis] / insetSurface.cut_surface.vertices_m.length), [0, 0, 0]);
insetSurface.cut_surface.vertices_m = insetSurface.cut_surface.vertices_m.map((p) => center.map((v, axis) => v + 0.2 * (p[axis] - v)));
assert.throws(() => buildFixture(insetSurface), /does not match the full through-section|inset or its boundary/);

process.stdout.write("explicit fragment fixture checks passed (including closure, pairing, conservation, caps, and inset-surface rejection)\n");
