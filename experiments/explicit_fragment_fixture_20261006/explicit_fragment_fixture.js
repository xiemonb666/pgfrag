"use strict";

// Bounded, dependency-free planar split fixture. This is deliberately not a
// fracture solver or a production truth-label producer.
const fs = require("node:fs");
const path = require("node:path");
const zlib = require("node:zlib");

const QUALIFICATION = "FIXTURE_GEOMETRY_ONLY";
const LABEL_NAMESPACE = "fixture_geometry_label.explicit_fragment.v1";
const OUTPUT_VERTEX_CAP = 64;
const OUTPUT_FACE_CAP = 128;
const EPS = 1e-9;

class FixtureGeometryError extends Error {
  constructor(message) {
    super(message);
    this.name = "FixtureGeometryError";
  }
}

function fail(message) {
  throw new FixtureGeometryError(message);
}

function add(a, b) { return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]; }
function sub(a, b) { return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]; }
function mul(a, k) { return [a[0] * k, a[1] * k, a[2] * k]; }
function dot(a, b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
function cross(a, b) {
  return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
}
function norm(a) { return Math.sqrt(dot(a, a)); }
function unit(a, name) {
  const length = norm(a);
  if (!Number.isFinite(length) || length <= EPS) fail(`${name} must have a finite nonzero length`);
  return mul(a, 1 / length);
}
function mean(points) {
  return mul(points.reduce((sum, p) => add(sum, p), [0, 0, 0]), 1 / points.length);
}
function close(a, b, tolerance = 1e-8) {
  return norm(sub(a, b)) <= tolerance;
}
function coordinateKey(p) { return p.map((v) => Math.round(v * 1e10)).join(","); }
function edgeKey(a, b) { return a < b ? `${a},${b}` : `${b},${a}`; }

function boxInput() {
  const vertices_m = [
    [-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
    [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1],
  ];
  const faces = [
    { id: "original:z_min", vertices: [0, 3, 2, 1] },
    { id: "original:z_max", vertices: [4, 5, 6, 7] },
    { id: "original:y_min", vertices: [0, 1, 5, 4] },
    { id: "original:y_max", vertices: [3, 7, 6, 2] },
    { id: "original:x_min", vertices: [0, 4, 7, 3] },
    { id: "original:x_max", vertices: [1, 2, 6, 5] },
  ];
  return {
    body_id: "fixture_box_01",
    density_kg_m3: 2400,
    vertices_m,
    faces,
  };
}

function boxUniqueEdges(body) {
  const edges = new Map();
  for (const face of body.faces) {
    for (let i = 0; i < face.vertices.length; i += 1) {
      const a = face.vertices[i];
      const b = face.vertices[(i + 1) % face.vertices.length];
      const key = edgeKey(a, b);
      if (!edges.has(key)) edges.set(key, [a, b]);
    }
  }
  return [...edges.values()];
}

function orderedPlaneSection(body, normal, offset) {
  const signed = body.vertices_m.map((p) => dot(normal, p) - offset);
  if (signed.some((distance) => Math.abs(distance) <= 1e-7)) {
    fail("cut plane through an original box vertex/edge is unsupported in this bounded fixture");
  }
  const hits = [];
  for (const [ia, ib] of boxUniqueEdges(body)) {
    const da = signed[ia];
    const db = signed[ib];
    if ((da < 0 && db > 0) || (da > 0 && db < 0)) {
      const t = da / (da - db);
      hits.push(add(body.vertices_m[ia], mul(sub(body.vertices_m[ib], body.vertices_m[ia]), t)));
    }
  }
  const unique = new Map();
  for (const p of hits) unique.set(coordinateKey(p), p);
  const points = [...unique.values()];
  if (points.length < 3 || points.length > 8) fail("plane does not produce one bounded polygonal section of the box");
  const center = mean(points);
  const reference = Math.abs(normal[2]) < 0.9 ? [0, 0, 1] : [0, 1, 0];
  const u = unit(cross(reference, normal), "plane basis u");
  const v = cross(normal, u);
  points.sort((a, b) => {
    const da = sub(a, center);
    const db = sub(b, center);
    return Math.atan2(dot(da, v), dot(da, u)) - Math.atan2(dot(db, v), dot(db, u));
  });
  return points;
}

function makeDefaultInput() {
  const original = boxInput();
  const rawNormal = [1, 0.37, 0.22];
  const length = norm(rawNormal);
  const normal = mul(rawNormal, 1 / length);
  const offset = 0.11 / length;
  const cutVertices = orderedPlaneSection(original, normal, offset);
  const triangles = [];
  for (let i = 1; i < cutVertices.length - 1; i += 1) triangles.push([0, i, i + 1]);
  return {
    schema_version: "pgfrag.explicit-fragment-fixture-input.v1",
    qualification: QUALIFICATION,
    fixture_note: "One declared planar cut through a small closed box. Its two output pieces are a geometry mechanism fixture only.",
    original_body: original,
    cut_surface: {
      surface_id: "explicit_cut_surface_01",
      geometry_kind: "explicit_planar_triangle_patch",
      plane: { normal, offset_m: offset },
      vertices_m: cutVertices,
      triangles,
      orientation: "input triangles point along plane normal; negative-side output keeps winding",
    },
  };
}

function validateOriginalBody(body) {
  if (!body || typeof body.body_id !== "string" || !body.body_id.trim()) fail("original body needs a body_id");
  if (!Array.isArray(body.vertices_m) || body.vertices_m.length < 4 || body.vertices_m.length > 64) fail("original body vertex count is outside the fixture cap");
  if (!Array.isArray(body.faces) || body.faces.length < 4 || body.faces.length > 128) fail("original body face count is outside the fixture cap");
  if (!(Number.isFinite(body.density_kg_m3) && body.density_kg_m3 > 0)) fail("density_kg_m3 must be finite and positive");
  for (const p of body.vertices_m) {
    if (!Array.isArray(p) || p.length !== 3 || !p.every(Number.isFinite)) fail("original body vertices must be finite 3D coordinates");
  }
  const edges = new Map();
  const ids = new Set();
  for (const face of body.faces) {
    if (!face || typeof face.id !== "string" || !face.id.trim() || ids.has(face.id)) fail("original face IDs must be non-empty and unique");
    ids.add(face.id);
    if (!Array.isArray(face.vertices) || face.vertices.length < 3) fail(`original face ${face.id} must be a polygon`);
    for (const index of face.vertices) {
      if (!Number.isInteger(index) || index < 0 || index >= body.vertices_m.length) fail(`original face ${face.id} has an invalid vertex index`);
    }
    const polygon = face.vertices.map((index) => body.vertices_m[index]);
    const areaVector = polygon.reduce((sum, p, i) => add(sum, cross(p, polygon[(i + 1) % polygon.length])), [0, 0, 0]);
    if (norm(areaVector) <= EPS) fail(`original face ${face.id} has zero area`);
    for (let i = 0; i < face.vertices.length; i += 1) {
      const a = face.vertices[i];
      const b = face.vertices[(i + 1) % face.vertices.length];
      const key = edgeKey(a, b);
      if (!edges.has(key)) edges.set(key, []);
      edges.get(key).push([a, b]);
    }
  }
  for (const [key, uses] of edges) {
    if (uses.length !== 2 || uses[0][0] !== uses[1][1] || uses[0][1] !== uses[1][0]) {
      fail(`original body is not a closed consistently oriented polygon mesh at edge ${key}`);
    }
  }
  return true;
}

function triangleArea(a, b, c) { return 0.5 * norm(cross(sub(b, a), sub(c, a))); }
function triangleUnitNormal(a, b, c) { return unit(cross(sub(b, a), sub(c, a)), "triangle normal"); }

function boundaryCycleFromTriangles(triangles, vertexCount) {
  const edges = new Map();
  for (const tri of triangles) {
    for (let i = 0; i < 3; i += 1) {
      const a = tri[i];
      const b = tri[(i + 1) % 3];
      const key = edgeKey(a, b);
      if (!edges.has(key)) edges.set(key, []);
      edges.get(key).push([a, b]);
    }
  }
  const boundary = [];
  for (const [key, uses] of edges) {
    if (uses.length > 2) fail(`cut surface is nonmanifold at edge ${key}`);
    if (uses.length === 2 && (uses[0][0] !== uses[1][1] || uses[0][1] !== uses[1][0])) fail(`cut surface winding disagrees at interior edge ${key}`);
    if (uses.length === 1) boundary.push(uses[0]);
  }
  if (!boundary.length) fail("cut surface must have one boundary loop that reaches the original exterior");
  const outgoing = new Map();
  const incoming = new Map();
  for (const [a, b] of boundary) {
    outgoing.set(a, (outgoing.get(a) || 0) + 1);
    incoming.set(b, (incoming.get(b) || 0) + 1);
  }
  const boundaryVertices = new Set([...incoming.keys(), ...outgoing.keys()]);
  if ([...boundaryVertices].some((vertex) => incoming.get(vertex) !== 1 || outgoing.get(vertex) !== 1)) fail("cut surface boundary must be one simple directed cycle");
  if (boundaryVertices.size !== boundary.length || boundaryVertices.size > vertexCount) fail("cut surface boundary topology is inconsistent");
  return boundary.map(([a, b]) => [a, b]);
}

function validateCutSurface(body, surface) {
  if (!surface || typeof surface.surface_id !== "string" || !surface.surface_id.trim()) fail("cut surface needs a non-empty surface_id");
  if (surface.geometry_kind !== "explicit_planar_triangle_patch") fail("this fixture accepts only an explicit planar triangle patch");
  if (!surface.plane || !Array.isArray(surface.plane.normal) || surface.plane.normal.length !== 3 || !surface.plane.normal.every(Number.isFinite) || !Number.isFinite(surface.plane.offset_m)) fail("cut plane must have a finite normal and offset");
  const originalNormalLength = norm(surface.plane.normal);
  const normal = unit(surface.plane.normal, "cut plane normal");
  const offset = surface.plane.offset_m / originalNormalLength;
  if (!Array.isArray(surface.vertices_m) || surface.vertices_m.length < 3 || surface.vertices_m.length > 16) fail("cut surface vertex count is outside the fixture cap");
  if (!Array.isArray(surface.triangles) || surface.triangles.length < 1 || surface.triangles.length > 64) fail("cut surface triangle count is outside the fixture cap");
  const vertices = surface.vertices_m;
  const triangles = surface.triangles;
  const scale = Math.max(...body.vertices_m.flatMap((p) => p.map(Math.abs)), 1);
  const tolerance = 1e-8 * scale;
  for (const p of vertices) {
    if (!Array.isArray(p) || p.length !== 3 || !p.every(Number.isFinite)) fail("cut surface vertices must be finite 3D coordinates");
    if (Math.abs(dot(normal, p) - offset) > tolerance) fail("cut surface vertices are not coplanar with the declared plane");
    if (p.some((value, axis) => value < Math.min(...body.vertices_m.map((q) => q[axis])) - tolerance || value > Math.max(...body.vertices_m.map((q) => q[axis])) + tolerance)) fail("cut surface lies outside the original body bounds");
  }
  for (let i = 0; i < vertices.length; i += 1) {
    for (let j = i + 1; j < vertices.length; j += 1) {
      if (close(vertices[i], vertices[j], tolerance)) fail("cut surface must not contain duplicate coordinates");
    }
  }
  const area = triangles.reduce((sum, tri, index) => {
    if (!Array.isArray(tri) || tri.length !== 3 || tri.some((v) => !Number.isInteger(v) || v < 0 || v >= vertices.length)) fail(`cut triangle ${index} has an invalid vertex index`);
    if (new Set(tri).size !== 3) fail(`cut triangle ${index} repeats a vertex`);
    const a = vertices[tri[0]];
    const b = vertices[tri[1]];
    const c = vertices[tri[2]];
    const faceArea = triangleArea(a, b, c);
    if (faceArea <= tolerance * tolerance) fail(`cut triangle ${index} is degenerate`);
    if (dot(triangleUnitNormal(a, b, c), normal) < 1 - 1e-8) fail(`cut triangle ${index} must point along the declared plane normal`);
    return sum + faceArea;
  }, 0);
  const boundary = boundaryCycleFromTriangles(triangles, vertices.length);
  if (boundary.length !== vertices.length) fail("the supported cut patch has boundary vertices only; interior patch vertices are not supported by this fixture");
  const expected = orderedPlaneSection(body, normal, offset);
  if (expected.length !== boundary.length) fail("explicit cut surface does not cover the full box plane section");
  const expectedKeys = new Set(expected.map(coordinateKey));
  const givenKeys = new Set(vertices.map(coordinateKey));
  if (expectedKeys.size !== givenKeys.size || [...expectedKeys].some((key) => !givenKeys.has(key))) fail("explicit cut surface boundary does not match the full through-section of the original body");
  const expectedEdges = new Set(expected.map((_, i) => edgeKey(coordinateKey(expected[i]), coordinateKey(expected[(i + 1) % expected.length]))));
  const boundaryEdges = new Set(boundary.map(([a, b]) => edgeKey(coordinateKey(vertices[a]), coordinateKey(vertices[b]))));
  if (expectedEdges.size !== boundaryEdges.size || [...expectedEdges].some((edge) => !boundaryEdges.has(edge))) fail("explicit cut surface is inset or its boundary does not follow the complete box intersection");
  let expectedArea = 0;
  for (let i = 1; i < expected.length - 1; i += 1) expectedArea += triangleArea(expected[0], expected[i], expected[i + 1]);
  if (Math.abs(area - expectedArea) > 1e-8 * Math.max(expectedArea, 1)) fail("cut triangles do not tile the complete box cross-section exactly once");
  return { normal, offset, area_m2: area, boundary, expected_section: expected };
}

function cleanPolygon(points) {
  const result = [];
  for (const p of points) {
    if (!result.length || !close(result[result.length - 1], p, 1e-8)) result.push(p);
  }
  if (result.length > 1 && close(result[0], result[result.length - 1], 1e-8)) result.pop();
  return result;
}

function clipPolygonToPlane(polygon, normal, offset, keep) {
  const signed = (p) => dot(normal, p) - offset;
  const inside = (p) => keep === "negative" ? signed(p) <= EPS : signed(p) >= -EPS;
  let previous = polygon[polygon.length - 1];
  let previousInside = inside(previous);
  const output = [];
  for (const current of polygon) {
    const currentInside = inside(current);
    if (currentInside !== previousInside) {
      const da = signed(previous);
      const db = signed(current);
      const t = da / (da - db);
      output.push(add(previous, mul(sub(current, previous), t)));
    }
    if (currentInside) output.push(current);
    previous = current;
    previousInside = currentInside;
  }
  return cleanPolygon(output);
}

function makeMeshBuilder(side, input) {
  const vertices = [];
  const lookup = new Map();
  const faces = [];
  const face_lineage = [];
  const surfaceVertices = input.cut_surface.vertices_m;
  const originalVertices = input.original_body.vertices_m;
  function canonical(p) {
    for (const source of surfaceVertices) if (close(source, p, 2e-8)) return source;
    for (const source of originalVertices) if (close(source, p, 2e-8)) return source;
    return p;
  }
  function addVertex(raw) {
    const p = canonical(raw);
    const key = coordinateKey(p);
    if (lookup.has(key)) return lookup.get(key);
    const index = vertices.length;
    vertices.push([...p]);
    lookup.set(key, index);
    return index;
  }
  function addTriangle(points, lineage) {
    const indices = points.map(addVertex);
    if (new Set(indices).size !== 3) fail("triangulation produced a repeated local vertex");
    const a = vertices[indices[0]];
    const b = vertices[indices[1]];
    const c = vertices[indices[2]];
    if (triangleArea(a, b, c) <= 1e-12) fail("triangulation produced a zero-area face");
    faces.push(indices);
    face_lineage.push(lineage);
  }
  function finish() {
    return {
      piece_id: side === "negative" ? "fixture_piece_negative" : "fixture_piece_positive",
      side,
      vertex_ids: vertices.map((_, index) => `${side}:v${index}`),
      vertices_m: vertices,
      triangles: faces,
      face_lineage,
    };
  }
  return { vertices, faces, face_lineage, addTriangle, finish };
}

function clippedHalf(input, plane, side) {
  const builder = makeMeshBuilder(side, input);
  for (const sourceFace of input.original_body.faces) {
    const polygon = sourceFace.vertices.map((index) => input.original_body.vertices_m[index]);
    const clipped = clipPolygonToPlane(polygon, plane.normal, plane.offset, side);
    if (clipped.length < 3) continue;
    for (let i = 1; i < clipped.length - 1; i += 1) {
      builder.addTriangle([clipped[0], clipped[i], clipped[i + 1]], {
        kind: "original_surface_face",
        label_namespace: LABEL_NAMESPACE,
        original_body_id: input.original_body.body_id,
        original_surface_face_id: sourceFace.id,
        fixture_geometry_only: true,
      });
    }
  }
  for (let sourceIndex = 0; sourceIndex < input.cut_surface.triangles.length; sourceIndex += 1) {
    const tri = input.cut_surface.triangles[sourceIndex];
    let points = tri.map((index) => input.cut_surface.vertices_m[index]);
    if (side === "positive") points = [points[0], points[2], points[1]];
    builder.addTriangle(points, {
      kind: "explicit_fracture_face_fixture",
      label_namespace: LABEL_NAMESPACE,
      surface_id: input.cut_surface.surface_id,
      input_surface_triangle_index: sourceIndex,
      pair_id: `${input.cut_surface.surface_id}:triangle:${sourceIndex}`,
      fixture_geometry_only: true,
    });
  }
  return builder.finish();
}

function meshTopology(mesh) {
  const edgeUses = new Map();
  const adjacency = Array.from({ length: mesh.vertices_m.length }, () => new Set());
  for (let faceIndex = 0; faceIndex < mesh.triangles.length; faceIndex += 1) {
    const tri = mesh.triangles[faceIndex];
    if (!Array.isArray(tri) || tri.length !== 3 || new Set(tri).size !== 3) fail(`${mesh.piece_id} has a malformed triangle`);
    for (const v of tri) if (!Number.isInteger(v) || v < 0 || v >= mesh.vertices_m.length) fail(`${mesh.piece_id} triangle index is invalid`);
    for (let i = 0; i < 3; i += 1) {
      const a = tri[i];
      const b = tri[(i + 1) % 3];
      const key = edgeKey(a, b);
      if (!edgeUses.has(key)) edgeUses.set(key, []);
      edgeUses.get(key).push([a, b]);
      adjacency[a].add(b);
      adjacency[b].add(a);
    }
  }
  for (const [edge, uses] of edgeUses) {
    if (uses.length !== 2) fail(`${mesh.piece_id} is not closed/manifold at edge ${edge}`);
    if (uses[0][0] !== uses[1][1] || uses[0][1] !== uses[1][0]) fail(`${mesh.piece_id} has inconsistent outward winding at edge ${edge}`);
  }
  const seen = new Set([0]);
  const stack = [0];
  while (stack.length) for (const next of adjacency[stack.pop()]) if (!seen.has(next)) { seen.add(next); stack.push(next); }
  if (seen.size !== mesh.vertices_m.length) fail(`${mesh.piece_id} is disconnected`);
  const euler = mesh.vertices_m.length - edgeUses.size + mesh.triangles.length;
  if (euler !== 2) fail(`${mesh.piece_id} Euler characteristic is ${euler}, expected 2`);
  return { vertex_count: mesh.vertices_m.length, edge_count: edgeUses.size, face_count: mesh.triangles.length, euler_characteristic: euler, closed: true, orientable_with_consistent_winding: true, connected_components: 1 };
}

function volumeCentroid(mesh) {
  let signedVolume = 0;
  let weighted = [0, 0, 0];
  for (const tri of mesh.triangles) {
    const a = mesh.vertices_m[tri[0]];
    const b = mesh.vertices_m[tri[1]];
    const c = mesh.vertices_m[tri[2]];
    const tetraVolume = dot(a, cross(b, c)) / 6;
    signedVolume += tetraVolume;
    weighted = add(weighted, mul(add(add(a, b), c), tetraVolume / 4));
  }
  if (!Number.isFinite(signedVolume) || signedVolume <= 0) fail(`${mesh.piece_id} has nonpositive signed volume; faces are not outward oriented`);
  return { volume_m3: signedVolume, centroid_m: mul(weighted, 1 / signedVolume) };
}

function faceUnitNormal(mesh, faceIndex) {
  const tri = mesh.triangles[faceIndex];
  return triangleUnitNormal(mesh.vertices_m[tri[0]], mesh.vertices_m[tri[1]], mesh.vertices_m[tri[2]]);
}

function findVertexByCoordinate(mesh, point) {
  const key = coordinateKey(point);
  const index = mesh.vertices_m.findIndex((candidate) => coordinateKey(candidate) === key && close(candidate, point, 1e-9));
  if (index < 0) fail(`${mesh.piece_id} is missing an explicit interface mother-coordinate`);
  return index;
}

function createPairLedger(input, negative, positive, normal) {
  const negativeByInput = new Map();
  const positiveByInput = new Map();
  for (let i = 0; i < negative.face_lineage.length; i += 1) {
    const label = negative.face_lineage[i];
    if (label.kind === "explicit_fracture_face_fixture") negativeByInput.set(label.input_surface_triangle_index, i);
  }
  for (let i = 0; i < positive.face_lineage.length; i += 1) {
    const label = positive.face_lineage[i];
    if (label.kind === "explicit_fracture_face_fixture") positiveByInput.set(label.input_surface_triangle_index, i);
  }
  const vertex_pairs = input.cut_surface.vertices_m.map((mother, inputVertexIndex) => ({
    input_surface_vertex_index: inputVertexIndex,
    mother_coordinate_m: [...mother],
    negative_vertex_id: negative.vertex_ids[findVertexByCoordinate(negative, mother)],
    positive_vertex_id: positive.vertex_ids[findVertexByCoordinate(positive, mother)],
    negative_local_vertex_index: findVertexByCoordinate(negative, mother),
    positive_local_vertex_index: findVertexByCoordinate(positive, mother),
  }));
  const pairs = input.cut_surface.triangles.map((tri, inputTriangleIndex) => {
    const negativeFace = negativeByInput.get(inputTriangleIndex);
    const positiveFace = positiveByInput.get(inputTriangleIndex);
    if (negativeFace === undefined || positiveFace === undefined) fail("an input interface triangle is not represented on both pieces");
    const negativeTri = negative.triangles[negativeFace];
    const positiveTri = positive.triangles[positiveFace];
    const negativePoints = negativeTri.map((index) => negative.vertices_m[index]);
    const positivePoints = positiveTri.map((index) => positive.vertices_m[index]);
    const sortedKey = (points) => points.map(coordinateKey).sort().join("|");
    if (sortedKey(negativePoints) !== sortedKey(positivePoints)) fail("paired interface triangles do not share identical mother coordinates");
    const nNegative = faceUnitNormal(negative, negativeFace);
    const nPositive = faceUnitNormal(positive, positiveFace);
    if (dot(nNegative, normal) < 1 - 1e-8 || dot(nPositive, normal) > -1 + 1e-8) fail("paired interface faces do not have opposite outward normals");
    return {
      pair_id: `${input.cut_surface.surface_id}:triangle:${inputTriangleIndex}`,
      input_surface_triangle_index: inputTriangleIndex,
      negative_face_index: negativeFace,
      positive_face_index: positiveFace,
      negative_face_id: `${negative.piece_id}:f${negativeFace}`,
      positive_face_id: `${positive.piece_id}:f${positiveFace}`,
      mother_coordinates_m: negativePoints,
      negative_outward_normal: nNegative,
      positive_outward_normal: nPositive,
      outward_normal_dot: dot(nNegative, nPositive),
    };
  });
  const reversePairs = vertex_pairs.map((pair) => ({
    negative_vertex_id: pair.negative_vertex_id,
    positive_vertex_id: pair.positive_vertex_id,
    mother_coordinate_m: [...pair.mother_coordinate_m],
  }));
  return {
    schema_version: "pgfrag.explicit-fragment-interface-pair-ledger.v1",
    qualification: QUALIFICATION,
    surface_id: input.cut_surface.surface_id,
    pairing_rule: "exact input-surface triangle identity and exact canonical mother coordinates; piece-local vertex IDs are distinct",
    vertex_pairs,
    reciprocal_vertex_pairs: reversePairs.map((pair) => ({ ...pair, negative_vertex_id: pair.positive_vertex_id, positive_vertex_id: pair.negative_vertex_id })),
    triangle_pairs: pairs,
  };
}

function originalMetrics(body) {
  const pseudo = { vertices_m: body.vertices_m, triangles: [] };
  for (const face of body.faces) for (let i = 1; i < face.vertices.length - 1; i += 1) pseudo.triangles.push([face.vertices[0], face.vertices[i], face.vertices[i + 1]]);
  return volumeCentroid(pseudo);
}

function buildFixture(input) {
  if (!input || input.schema_version !== "pgfrag.explicit-fragment-fixture-input.v1") fail("input schema_version is not supported");
  if (input.qualification !== QUALIFICATION) fail("input qualification must be explicitly FIXTURE_GEOMETRY_ONLY");
  validateOriginalBody(input.original_body);
  const plane = validateCutSurface(input.original_body, input.cut_surface);
  const negative = clippedHalf(input, plane, "negative");
  const positive = clippedHalf(input, plane, "positive");
  const negativeTopology = meshTopology(negative);
  const positiveTopology = meshTopology(positive);
  const negativeMetrics = volumeCentroid(negative);
  const positiveMetrics = volumeCentroid(positive);
  const original = originalMetrics(input.original_body);
  const totalVertexCount = input.original_body.vertices_m.length + input.cut_surface.vertices_m.length + negative.vertices_m.length + positive.vertices_m.length;
  const totalFaceCount = input.cut_surface.triangles.length + negative.triangles.length + positive.triangles.length;
  if (totalVertexCount > OUTPUT_VERTEX_CAP) fail(`combined declared geometry has ${totalVertexCount} vertices, above cap ${OUTPUT_VERTEX_CAP}`);
  if (totalFaceCount > OUTPUT_FACE_CAP) fail(`combined declared geometry has ${totalFaceCount} faces, above cap ${OUTPUT_FACE_CAP}`);
  const ledger = createPairLedger(input, negative, positive, plane.normal);
  const density = input.original_body.density_kg_m3;
  negativeMetrics.mass_kg = negativeMetrics.volume_m3 * density;
  positiveMetrics.mass_kg = positiveMetrics.volume_m3 * density;
  original.mass_kg = original.volume_m3 * density;
  const volumeError = Math.abs(negativeMetrics.volume_m3 + positiveMetrics.volume_m3 - original.volume_m3);
  const massError = Math.abs(negativeMetrics.mass_kg + positiveMetrics.mass_kg - original.mass_kg);
  const sideNegative = dot(plane.normal, negativeMetrics.centroid_m) - plane.offset;
  const sidePositive = dot(plane.normal, positiveMetrics.centroid_m) - plane.offset;
  if (!(negativeMetrics.volume_m3 > 0 && positiveMetrics.volume_m3 > 0 && sideNegative < -1e-8 && sidePositive > 1e-8)) fail("cut does not produce two pieces with positive interior thickness");
  if (volumeError > 1e-9 * Math.max(original.volume_m3, 1)) fail("piece volumes do not sum to the original body volume");
  if (massError > 1e-9 * Math.max(original.mass_kg, 1)) fail("piece masses do not sum to the original body mass");
  const allPieceSignsValid = [
    ...negative.vertices_m.map((p) => dot(plane.normal, p) - plane.offset <= 1e-8),
    ...positive.vertices_m.map((p) => dot(plane.normal, p) - plane.offset >= -1e-8),
  ].every(Boolean);
  if (!allPieceSignsValid) fail("piece vertex lies across the opposite half-space; interior overlap check failed");
  const interfaceAreaNegative = negative.face_lineage.reduce((sum, label, i) => label.kind === "explicit_fracture_face_fixture" ? sum + triangleArea(...negative.triangles[i].map((v) => negative.vertices_m[v])) : sum, 0);
  const interfaceAreaPositive = positive.face_lineage.reduce((sum, label, i) => label.kind === "explicit_fracture_face_fixture" ? sum + triangleArea(...positive.triangles[i].map((v) => positive.vertices_m[v])) : sum, 0);
  if (Math.abs(interfaceAreaNegative - interfaceAreaPositive) > 1e-10) fail("paired interface areas disagree");
  const result = {
    schema_version: "pgfrag.explicit-fragment-fixture-output.v1",
    qualification: QUALIFICATION,
    mechanics_qualification: "NOT_COUPLED",
    production_promotion: "BLOCKED_FIXTURE_ONLY",
    physical_fracture_claim: "NOT_CLAIMED",
    relation_labels: "NONE",
    fixture_note: "The requested plane is the declared generative input. A split count of two is a fixture consequence and does not model physical growth or fragmentation.",
    original_body: { body_id: input.original_body.body_id, density_kg_m3: density, ...original },
    cut_surface: {
      surface_id: input.cut_surface.surface_id,
      geometry_kind: input.cut_surface.geometry_kind,
      plane: { normal: plane.normal, offset_m: plane.offset },
      vertices_m: input.cut_surface.vertices_m,
      triangles: input.cut_surface.triangles,
      area_m2: plane.area_m2,
      lineage: { parent_body_id: input.original_body.body_id, source_kind: "explicit_fixture_input_surface", production_fracture_event_id: null },
    },
    pieces: [
      { ...negative, metrics: negativeMetrics, topology: negativeTopology },
      { ...positive, metrics: positiveMetrics, topology: positiveTopology },
    ],
    interface_pair_ledger: ledger,
    checks: {
      closed_outward_manifold_meshes: true,
      piece_euler_characteristic_2: true,
      volume_sum_matches_original: true,
      mass_sum_matches_original: true,
      paired_coordinates_identical: true,
      pair_vertex_map_reciprocal: true,
      interface_outward_normals_opposed: true,
      interface_area_matches: true,
      interiors_disjoint_with_positive_thickness: true,
      all_output_faces_have_source_labels: true,
      bounded_geometry_caps_respected: { total_declared_vertices: totalVertexCount, maximum_vertices: OUTPUT_VERTEX_CAP, total_declared_faces: totalFaceCount, maximum_faces: OUTPUT_FACE_CAP },
      numerical_errors: { volume_error_m3: volumeError, mass_error_kg: massError, interface_area_difference_m2: Math.abs(interfaceAreaNegative - interfaceAreaPositive) },
    },
  };
  return result;
}

function writeObj(mesh, filePath) {
  const lines = ["# FIXTURE_GEOMETRY_ONLY; not a physical fracture surface or production label", `o ${mesh.piece_id}`];
  for (const p of mesh.vertices_m) lines.push(`v ${p[0].toPrecision(17)} ${p[1].toPrecision(17)} ${p[2].toPrecision(17)}`);
  let activeGroup = "";
  for (let i = 0; i < mesh.triangles.length; i += 1) {
    const lineage = mesh.face_lineage[i];
    const group = lineage.kind === "original_surface_face" ? lineage.original_surface_face_id : `${lineage.surface_id}_fixture_cut`;
    if (group !== activeGroup) { lines.push(`g ${group.replace(/[^A-Za-z0-9_.-]/g, "_")}`); activeGroup = group; }
    lines.push(`f ${mesh.triangles[i].map((index) => index + 1).join(" ")}`);
  }
  fs.writeFileSync(filePath, `${lines.join("\n")}\n`, "utf8");
}

const FONT = {
  A: ["01110", "10001", "10001", "11111", "10001", "10001", "10001"], B: ["11110", "10001", "10001", "11110", "10001", "10001", "11110"],
  C: ["01111", "10000", "10000", "10000", "10000", "10000", "01111"], E: ["11111", "10000", "10000", "11110", "10000", "10000", "11111"],
  F: ["11111", "10000", "10000", "11110", "10000", "10000", "10000"], G: ["01111", "10000", "10000", "10111", "10001", "10001", "01111"],
  H: ["10001", "10001", "10001", "11111", "10001", "10001", "10001"], I: ["11111", "00100", "00100", "00100", "00100", "00100", "11111"],
  L: ["10000", "10000", "10000", "10000", "10000", "10000", "11111"], M: ["10001", "11011", "10101", "10101", "10001", "10001", "10001"],
  N: ["10001", "11001", "10101", "10011", "10001", "10001", "10001"], O: ["01110", "10001", "10001", "10001", "10001", "10001", "01110"],
  P: ["11110", "10001", "10001", "11110", "10000", "10000", "10000"], R: ["11110", "10001", "10001", "11110", "10100", "10010", "10001"],
  S: ["01111", "10000", "10000", "01110", "00001", "00001", "11110"], T: ["11111", "00100", "00100", "00100", "00100", "00100", "00100"],
  U: ["10001", "10001", "10001", "10001", "10001", "10001", "01110"], X: ["10001", "10001", "01010", "00100", "01010", "10001", "10001"], Y: ["10001", "10001", "01010", "00100", "00100", "00100", "00100"],
  0: ["01110", "10001", "10011", "10101", "11001", "10001", "01110"], 1: ["00100", "01100", "00100", "00100", "00100", "00100", "01110"],
  2: ["01110", "10001", "00001", "00010", "00100", "01000", "11111"], 3: ["11110", "00001", "00001", "01110", "00001", "00001", "11110"],
  "-": ["00000", "00000", "00000", "11111", "00000", "00000", "00000"], ":": ["00000", "00100", "00100", "00000", "00100", "00100", "00000"],
};

function pngCrc32(buffer) {
  let crc = 0xffffffff;
  for (const byte of buffer) {
    crc ^= byte;
    for (let bit = 0; bit < 8; bit += 1) crc = (crc >>> 1) ^ ((crc & 1) ? 0xedb88320 : 0);
  }
  return (crc ^ 0xffffffff) >>> 0;
}

function pngChunk(type, data) {
  const typeBytes = Buffer.from(type, "ascii");
  const body = Buffer.concat([typeBytes, data]);
  const length = Buffer.alloc(4);
  length.writeUInt32BE(data.length, 0);
  const crc = Buffer.alloc(4);
  crc.writeUInt32BE(pngCrc32(body), 0);
  return Buffer.concat([length, body, crc]);
}

function writeRgbaPng(width, height, rgba, filePath) {
  const rows = Buffer.alloc((width * 4 + 1) * height);
  for (let y = 0; y < height; y += 1) {
    const offset = y * (width * 4 + 1);
    rows[offset] = 0;
    Buffer.from(rgba.buffer, rgba.byteOffset + y * width * 4, width * 4).copy(rows, offset + 1);
  }
  const header = Buffer.alloc(13);
  header.writeUInt32BE(width, 0);
  header.writeUInt32BE(height, 4);
  header[8] = 8;
  header[9] = 6;
  header[10] = 0;
  header[11] = 0;
  header[12] = 0;
  const png = Buffer.concat([
    Buffer.from([137, 80, 78, 71, 13, 10, 26, 10]),
    pngChunk("IHDR", header),
    pngChunk("IDAT", zlib.deflateSync(rows, { level: 7 })),
    pngChunk("IEND", Buffer.alloc(0)),
  ]);
  fs.writeFileSync(filePath, png);
}

function color(hex) { return [parseInt(hex.slice(1, 3), 16), parseInt(hex.slice(3, 5), 16), parseInt(hex.slice(5, 7), 16)]; }
function drawImageText(rgba, width, height, text, x, y, scale, rgb) {
  const ink = text.toUpperCase();
  let cursor = x;
  for (const character of ink) {
    if (character === " ") { cursor += 4 * scale; continue; }
    const glyph = FONT[character];
    if (!glyph) { cursor += 6 * scale; continue; }
    for (let gy = 0; gy < 7; gy += 1) for (let gx = 0; gx < 5; gx += 1) {
      if (glyph[gy][gx] !== "1") continue;
      for (let dy = 0; dy < scale; dy += 1) for (let dx = 0; dx < scale; dx += 1) {
        const px = cursor + gx * scale + dx;
        const py = y + gy * scale + dy;
        if (px >= 0 && px < width && py >= 0 && py < height) {
          const i = (py * width + px) * 4;
          rgba[i] = rgb[0]; rgba[i + 1] = rgb[1]; rgba[i + 2] = rgb[2]; rgba[i + 3] = 255;
        }
      }
    }
    cursor += 6 * scale;
  }
}

function renderExplodedPreview(result, filePath) {
  const width = 1280;
  const height = 800;
  const rgba = new Uint8Array(width * height * 4);
  for (let y = 0; y < height; y += 1) for (let x = 0; x < width; x += 1) {
    const i = (y * width + x) * 4;
    const t = y / height;
    rgba[i] = Math.round(249 - t * 8); rgba[i + 1] = Math.round(251 - t * 8); rgba[i + 2] = Math.round(253 - t * 8); rgba[i + 3] = 255;
  }
  const view = unit([0.55, -0.78, 0.35], "preview camera");
  const worldUp = [0, 0, 1];
  const right = unit(cross(view, worldUp), "preview right axis");
  const up = unit(cross(right, view), "preview up axis");
  const normal = result.cut_surface.plane.normal;
  const gap = 0.23;
  const panels = [
    { mesh: result.pieces[0], shift: mul(normal, -gap / 2), fill: color("#7198bd"), title: "PIECE A" },
    { mesh: result.pieces[1], shift: mul(normal, gap / 2), fill: color("#50aa9c"), title: "PIECE B" },
  ];
  const projected = [];
  for (const panel of panels) for (let fi = 0; fi < panel.mesh.triangles.length; fi += 1) {
    const tri = panel.mesh.triangles[fi];
    const world = tri.map((index) => add(panel.mesh.vertices_m[index], panel.shift));
    const lineage = panel.mesh.face_lineage[fi];
    const normal3 = triangleUnitNormal(world[0], world[1], world[2]);
    const center3 = mean(world);
    projected.push({
      points: world.map((p) => [dot(p, right), dot(p, up), dot(p, view)]),
      world,
      normal3,
      center3,
      panel,
      cut: lineage.kind === "explicit_fracture_face_fixture",
      mesh_face_index: fi,
    });
  }
  let minX = Infinity; let maxX = -Infinity; let minY = Infinity; let maxY = -Infinity;
  for (const face of projected) for (const p of face.points) { minX = Math.min(minX, p[0]); maxX = Math.max(maxX, p[0]); minY = Math.min(minY, p[1]); maxY = Math.max(maxY, p[1]); }
  const availableWidth = 1040;
  const availableHeight = 590;
  const scale = Math.min(availableWidth / (maxX - minX), availableHeight / (maxY - minY));
  const midX = (minX + maxX) / 2;
  const midY = (minY + maxY) / 2;
  const project = (p) => [Math.round(width / 2 + (p[0] - midX) * scale), Math.round(height / 2 + 15 - (p[1] - midY) * scale), p[2]];
  const depth = new Float32Array(width * height);
  depth.fill(-Infinity);
  function putPixel(x, y, z, rgb, depthBias = 0) {
    if (x < 0 || x >= width || y < 0 || y >= height) return;
    const pixel = y * width + x;
    if (z + depthBias < depth[pixel] - 1e-5) return;
    depth[pixel] = Math.max(depth[pixel], z + depthBias);
    const i = pixel * 4;
    rgba[i] = rgb[0]; rgba[i + 1] = rgb[1]; rgba[i + 2] = rgb[2]; rgba[i + 3] = 255;
  }
  for (const face of projected) {
    const pts = face.points.map(project);
    const minPX = Math.max(0, Math.floor(Math.min(...pts.map((p) => p[0]))));
    const maxPX = Math.min(width - 1, Math.ceil(Math.max(...pts.map((p) => p[0]))));
    const minPY = Math.max(0, Math.floor(Math.min(...pts.map((p) => p[1]))));
    const maxPY = Math.min(height - 1, Math.ceil(Math.max(...pts.map((p) => p[1]))));
    const [a, b, c] = pts;
    const denom = (b[1] - c[1]) * (a[0] - c[0]) + (c[0] - b[0]) * (a[1] - c[1]);
    if (Math.abs(denom) < 1e-8) continue;
    const light = unit([0.45, -0.3, 0.84], "preview light");
    const shade = 0.57 + 0.43 * Math.max(0, Math.abs(dot(face.normal3, light)));
    const base = face.cut ? color("#ed785b") : face.panel.fill;
    const shaded = base.map((channel) => Math.round(channel * shade));
    for (let y = minPY; y <= maxPY; y += 1) for (let x = minPX; x <= maxPX; x += 1) {
      const wa = ((b[1] - c[1]) * (x - c[0]) + (c[0] - b[0]) * (y - c[1])) / denom;
      const wb = ((c[1] - a[1]) * (x - c[0]) + (a[0] - c[0]) * (y - c[1])) / denom;
      const wc = 1 - wa - wb;
      if (wa >= -1e-6 && wb >= -1e-6 && wc >= -1e-6) putPixel(x, y, wa * a[2] + wb * b[2] + wc * c[2], shaded);
    }
  }
  function drawLine(a, b, rgb, bias = 1e-4) {
    const pa = project(a); const pb = project(b);
    const steps = Math.max(Math.abs(pb[0] - pa[0]), Math.abs(pb[1] - pa[1]), 1);
    for (let i = 0; i <= steps; i += 1) {
      const t = i / steps;
      const x = Math.round(pa[0] + (pb[0] - pa[0]) * t);
      const y = Math.round(pa[1] + (pb[1] - pa[1]) * t);
      const z = pa[2] + (pb[2] - pa[2]) * t;
      for (let oy = -1; oy <= 1; oy += 1) for (let ox = -1; ox <= 1; ox += 1) putPixel(x + ox, y + oy, z, rgb, bias);
    }
  }
  for (const panel of panels) {
    const edgeFaces = new Map();
    for (let fi = 0; fi < panel.mesh.triangles.length; fi += 1) {
      const tri = panel.mesh.triangles[fi];
      for (let i = 0; i < 3; i += 1) {
        const a = tri[i]; const b = tri[(i + 1) % 3];
        const key = edgeKey(a, b);
        if (!edgeFaces.has(key)) edgeFaces.set(key, []);
        edgeFaces.get(key).push({ faceIndex: fi, a, b });
      }
    }
    for (const uses of edgeFaces.values()) {
      if (uses.length !== 2) continue;
      const labels = uses.map((use) => panel.mesh.face_lineage[use.faceIndex]);
      const lineageKey = (label) => label.kind === "original_surface_face"
        ? `${label.kind}:${label.original_surface_face_id}`
        : `${label.kind}:${label.surface_id}`;
      if (lineageKey(labels[0]) === lineageKey(labels[1])) continue;
      const isCutBoundary = labels.some((label) => label.kind === "explicit_fracture_face_fixture");
      const shift = panel.shift;
      const a = add(panel.mesh.vertices_m[uses[0].a], shift);
      const b = add(panel.mesh.vertices_m[uses[0].b], shift);
      drawLine(a, b, color(isCutBoundary ? "#9e432f" : "#244054"), 2e-4);
    }
  }
  const ink = color("#263a4b");
  drawImageText(rgba, width, height, "FIXTURE GEOMETRY ONLY", 42, 36, 4, ink);
  drawImageText(rgba, width, height, panels[0].title, 66, 735, 3, panels[0].fill);
  drawImageText(rgba, width, height, panels[1].title, 980, 735, 3, panels[1].fill);
  drawImageText(rgba, width, height, "ORANGE: PAIRED INPUT CUT FACES", 340, 754, 2, color("#9e432f"));
  writeRgbaPng(width, height, rgba, filePath);
}

function writeReviewOutputs(result, outputDir) {
  fs.mkdirSync(outputDir, { recursive: true });
  fs.writeFileSync(path.join(outputDir, "explicit_fragment_fixture_output.json"), `${JSON.stringify(result, null, 2)}\n`, "utf8");
  writeObj(result.pieces[0], path.join(outputDir, "piece_negative.obj"));
  writeObj(result.pieces[1], path.join(outputDir, "piece_positive.obj"));
  fs.writeFileSync(path.join(outputDir, "interface_pair_ledger.json"), `${JSON.stringify(result.interface_pair_ledger, null, 2)}\n`, "utf8");
  fs.writeFileSync(path.join(outputDir, "checks.json"), `${JSON.stringify(result.checks, null, 2)}\n`, "utf8");
  renderExplodedPreview(result, path.join(outputDir, "exploded_mesh_preview.png"));
}

function main(argv) {
  const here = __dirname;
  const inputPath = path.join(here, "input_fixture.json");
  const defaultInput = makeDefaultInput();
  if (argv.includes("--write-default-input")) {
    fs.writeFileSync(inputPath, `${JSON.stringify(defaultInput, null, 2)}\n`, "utf8");
    process.stdout.write(`${inputPath}\n`);
    return;
  }
  if (!fs.existsSync(inputPath)) fail("input_fixture.json is missing; run this module once with --write-default-input");
  const input = JSON.parse(fs.readFileSync(inputPath, "utf8"));
  const result = buildFixture(input);
  const outputDir = path.resolve(here, "../../tmp/review_20261006/explicit_fragment_fixture_20261006");
  writeReviewOutputs(result, outputDir);
  process.stdout.write(`${JSON.stringify({ outputDir, checks: result.checks }, null, 2)}\n`);
}

if (require.main === module) {
  try { main(process.argv.slice(2)); }
  catch (error) {
    process.stderr.write(`${error.name || "Error"}: ${error.message}\n`);
    process.exitCode = 1;
  }
}

module.exports = { FixtureGeometryError, buildFixture, makeDefaultInput, validateCutSurface, meshTopology, volumeCentroid };
