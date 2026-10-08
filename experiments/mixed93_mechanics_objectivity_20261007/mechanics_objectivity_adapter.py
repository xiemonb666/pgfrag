"""Frame-aware two-state mixed93 adapter. Importing this file performs no work.

capture_arrays is deliberately locked until exact ROOT authorization. The caller
must verify source pins and provide frozen quadrature arrays; this module never
opens the cached NPZ, and never repairs K or reuses a spectrum.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass

MAX_CHUNK = 256


class AdapterBlocked(RuntimeError):
    """Fail-closed support, geometry, or numerical-state result."""


@dataclass(frozen=True)
class FrameState:
    q: object
    t: object
    axes: object
    origin: object
    tip: object
    center: object
    source_xyz: object
    volume_xyz: object
    boundary_xyz: object
    normals: object
    mls_config: object
    mls_origin_world_m: object
    field_case: object
    body_axis_world: object
    front_direction_world: object
    crack_plane_normal_world: object


def require_root_go(case: dict, token: str | None) -> None:
    auth = case.get("execution_authorization", {})
    if (case.get("numeric_execution_authorized") is not True
            or auth.get("authorized") is not True
            or not isinstance(auth.get("root_go_token"), str)
            or not auth["root_go_token"] or token != auth["root_go_token"]):
        raise PermissionError("blocked until ROOTFULLREAD and exact root-go token")


def rotation_matrix(np, record: dict):
    az, ax = float(record["angle_z_rad"]), float(record["angle_x_rad"])
    rz = np.array([[math.cos(az), -math.sin(az), 0.0],
                   [math.sin(az), math.cos(az), 0.0], [0.0, 0.0, 1.0]])
    rx = np.array([[1.0, 0.0, 0.0],
                   [0.0, math.cos(ax), -math.sin(ax)],
                   [0.0, math.sin(ax), math.cos(ax)]])
    q = rz @ rx
    if not np.allclose(q.T @ q, np.eye(3), rtol=0.0, atol=1e-12) or abs(float(np.linalg.det(q))-1.0) > 1e-12:
        raise AdapterBlocked("rotation must be proper orthogonal")
    return q


def _vector(np, value, name):
    v = np.asarray(value, dtype=np.float64)
    if v.shape != (3,) or not bool(np.isfinite(v).all()):
        raise AdapterBlocked(f"{name} must be finite xyz")
    return v


def _points(xyz, q, t):
    return xyz @ q.T + t


def _strain_tensor(np, v):
    v = np.asarray(v, dtype=np.float64)
    t = np.zeros(v.shape[:-1] + (3, 3), dtype=np.float64)
    t[..., 0, 0], t[..., 1, 1], t[..., 2, 2] = v[..., 0], v[..., 1], v[..., 2]
    t[..., 0, 1] = t[..., 1, 0] = 0.5 * v[..., 3]
    t[..., 1, 2] = t[..., 2, 1] = 0.5 * v[..., 4]
    t[..., 0, 2] = t[..., 2, 0] = 0.5 * v[..., 5]
    return t


def _stress_tensor(np, v):
    v = np.asarray(v, dtype=np.float64)
    t = np.zeros(v.shape[:-1] + (3, 3), dtype=np.float64)
    t[..., 0, 0], t[..., 1, 1], t[..., 2, 2] = v[..., 0], v[..., 1], v[..., 2]
    t[..., 0, 1] = t[..., 1, 0] = v[..., 3]
    t[..., 1, 2] = t[..., 2, 1] = v[..., 4]
    t[..., 0, 2] = t[..., 2, 0] = v[..., 5]
    return t


def _field_case(field, case, axes, origin):
    fr = field.CrackFrame(tuple(origin), tuple(tuple(float(v) for v in row) for row in axes))
    g, m, load = case["geometry"], case["material"], case["load"]
    return field.MixedLEFMCase(float(g["radius_m"]), float(g["length_m"]),
        *map(float, load["K_Pa_sqrt_m"]), float(m["young_modulus_Pa"]),
        float(m["poisson_ratio"]), frame=fr)


def make_frame_state(np, deps, case, inputs, q, t, base_frame):
    mls, field = deps["mls"], deps["field"]
    axes0 = np.asarray(base_frame["axes_world_from_local"], dtype=np.float64)
    origin0 = _vector(np, base_frame["origin_world_m"], "crack origin")
    tip0 = _vector(np, base_frame["tip_world_m"], "crack tip")
    if not np.array_equal(origin0, tip0):
        raise AdapterBlocked("fixture requires branch tip and CrackFrame origin to coincide")
    center0 = _vector(np, case["gauge"]["center_xyz_m"], "moment center")
    axes, origin, tip = q @ axes0, q @ origin0 + t, q @ tip0 + t
    center = q @ center0 + t
    xyz0 = inputs["volume_xyz_m"]
    bxyz0 = inputs["boundary_xyz_m"]
    source = _points(inputs["source_positions_xyz_m"], q, t)
    volume, boundary = _points(xyz0, q, t), _points(bxyz0, q, t)
    normals = inputs["boundary_normals"] @ q.T
    body_axis = q @ _vector(np, case["geometry"]["body_axis_world"], "body axis")
    front_direction = q @ _vector(np, case["geometry"]["front_direction_world"], "front direction")
    crack_normal = q @ _vector(np, case["geometry"]["crack_plane_normal_world"], "crack plane normal")
    cfg = case["mls"]
    mls_origin = q @ _vector(np, cfg["reference_origin_m"], "MLS origin") + t
    config = mls.MLS3DConfig(tuple(mls_origin), float(cfg["reference_scale_m"]),
        float(cfg["support_radius_m"]), float(cfg["max_condition_number"]),
        float(cfg["rank_relative_tolerance"]), float(cfg["switch_tolerance"]),
        float(cfg["duplicate_tolerance"]))
    return FrameState(q, t, axes, origin, tip, center, source, volume, boundary,
                      normals, config, mls_origin,
                      _field_case(field, case, axes, origin),
                      body_axis, front_direction, crack_normal)


def _b(component, g):
    x, y, z = map(float, g)
    if component == "ux": return (x, 0., 0., y, 0., z)
    if component == "uy": return (0., y, 0., x, z, 0.)
    return (0., 0., z, 0., y, x)


def _row(np, op, field, source_ids, base_values, base_gradients, local, axes, radius):
    if local[0] < 0.0 and local[1] == 0.0:
        raise AdapterBlocked("crack-cut row needs an explicitly pinned trace side")
    br = field.branch_functions_local(float(local[0]), float(local[1]), reference_length_m=radius)
    vals = np.concatenate((base_values, np.asarray(br.values, dtype=np.float64)))
    bgrad = np.asarray(br.gradients_per_m, dtype=np.float64) @ axes.T
    grads = np.concatenate((base_gradients, bgrad), axis=0)
    slots = {s: i for i, s in enumerate(source_ids)}
    ids = op.semantic_coefficient_ids(source_ids)
    items = []
    for sid in ids:
        scope, owner, comp, index = json.loads(sid)
        if scope == "SOURCE_BASE": i, value, grad = slots[owner], vals[slots[owner]], grads[slots[owner]]
        elif scope == "GLOBAL_FRONT": i = 27 + index - 1; value, grad = vals[i], grads[i]
        else: raise AdapterBlocked(f"unknown semantic scope {scope}")
        g = tuple(float(v) for v in grad)
        items.append(op.SemanticContribution(sid, scope, owner, comp, index,
                                             float(value), g, _b(comp, g)))
    return op.SemanticShapeRow(op.GLOBAL_FRONT_LAYOUT, source_ids, ids, tuple(items), "OK"), vals, grads


def _group(np, deps, case, source_ids, state, xyz, local, weights, normals, surface):
    shape, field, op = deps["shape"], deps["field"], deps["operator"]
    vals, grads = np.empty((len(xyz), 31)), np.empty((len(xyz), 31, 3))
    points, tractions, cmax = [], [], 0.0
    for start in range(0, len(xyz), MAX_CHUNK):
        end = min(start + MAX_CHUNK, len(xyz))
        r = shape.shape_functions_3d(state.source_xyz,
                                     np.asarray(xyz[start:end], dtype=np.float64), state.mls_config)
        if any(getattr(s, "value", s) != "OK" for s in r.statuses) or not bool(np.all(r.support_counts == 27)) or not bool(np.isfinite(r.condition_numbers).all()):
            raise AdapterBlocked("fresh MLS support/status/condition gate failed")
        cmax = max(cmax, float(np.max(r.condition_numbers)))
        chunk = []
        for j, ix in enumerate(range(start, end)):
            row, vals[ix], grads[ix] = _row(np, op, field, source_ids,
                r.weights[j], r.gradients[j], local[ix], state.axes,
                float(case["geometry"]["radius_m"]))
            if surface:
                exact = field.evaluate_local_point(tuple(local[ix]), case=state.field_case)
                if getattr(exact.status, "value", exact.status) != "OK" or exact.stress_tensor_world_Pa is None:
                    raise AdapterBlocked("analytic boundary stress status is not OK")
                tr = np.asarray(exact.stress_tensor_world_Pa) @ normals[ix]
                tractions.append(tr)
                chunk.append(op.SurfaceTractionPoint(float(weights[ix]), tuple(map(float, tr)), row))
            else:
                chunk.append(op.VolumePoint(float(weights[ix]), row))
        points.append(tuple(chunk))
    return tuple(points), vals, grads, np.asarray(tractions) if surface else None, cmax


def _rel(np, x, y):
    return float(np.linalg.norm(x-y) / max(float(np.linalg.norm(y)), np.finfo(float).tiny))


def _weighted(np, x, y, weights):
    w = weights.reshape((-1,) + (1,) * (x.ndim-1))
    return math.sqrt(float(np.sum(w*(x-y)**2)) / max(float(np.sum(w*y*y)), np.finfo(float).tiny))


def _wrench(np, xyz, weights, traction, center):
    f = np.einsum("q,qi->i", weights, traction)
    mc = np.einsum("q,qi->i", weights, np.cross(xyz-center, traction))
    return f, mc, mc + np.cross(center, f), float(np.dot(weights, np.linalg.norm(traction, axis=1)))


def _state(np, deps, case, inputs, source_ids, q, t, base_frame):
    op, worker = deps["operator"], deps["worker"]
    s = make_frame_state(np, deps, case, inputs, q, t, base_frame)
    # Local coordinates come from the frozen baseline, so crack-face branch selection is identical.
    vlocal = (inputs["volume_xyz_m"] - np.asarray(base_frame["tip_world_m"], dtype=np.float64)) @ np.asarray(base_frame["axes_world_from_local"], dtype=np.float64)
    blocal = (inputs["boundary_xyz_m"] - np.asarray(base_frame["tip_world_m"], dtype=np.float64)) @ np.asarray(base_frame["axes_world_from_local"], dtype=np.float64)
    vchunks, vv, vg, _, cv = _group(np, deps, case, source_ids, s, s.volume_xyz,
        vlocal, inputs["volume_weights_m3"], None, False)
    bchunks, bv, bg, traction, cb = _group(np, deps, case, source_ids, s, s.boundary_xyz,
        blocal, inputs["boundary_weights_m2"], s.normals, True)
    d = np.asarray(op.isotropic_constitutive_matrix_pa(
        case["material"]["young_modulus_Pa"], case["material"]["poisson_ratio"]), dtype=np.float64)
    kr = op.assemble_raw_stiffness(source_ids, vchunks, tuple(map(tuple, d)))
    fr = op.assemble_raw_force(source_ids, bchunks)
    k = np.frombuffer(kr.values_le_f64, dtype="<f8").copy().reshape(kr.matrix_shape)
    f = np.frombuffer(fr.values_le_f64, dtype="<f8").copy().reshape(fr.vector_shape)
    relxyz = s.volume_xyz - s.center
    c, vmeasure, _ = worker._physical_moment_matrix(np, vv, relxyz, inputs["volume_weights_m3"], case["geometry"]["radius_m"])
    ks, scale, _, spectrum = worker._scaled_spectrum_gate(np, k)
    if not spectrum["passed"] or (spectrum["zero_count"], spectrum["positive_count"], spectrum["significant_negative_count"]) != (6,87,0):
        raise AdapterBlocked("fresh spectrum requires 6 zero, 87 positive, 0 negative modes")
    g = np.zeros(6)
    u, lam, solve, _ = worker._solve_case(np, k, ks, c, f, scale, g)
    disp, eps, sig, _ = worker._reconstruct_field(np, vv, vg, u, d)
    moments = worker._direct_physical_moments(np, disp, relxyz, inputs["volume_weights_m3"], case["geometry"]["radius_m"], vmeasure)
    wrench = _wrench(np, s.boundary_xyz, inputs["boundary_weights_m2"], traction, s.center)
    return {"frame":s, "K":k, "C":c, "f":f, "g":g, "q":u, "lambda":lam,
        "spectrum":spectrum, "solve":solve, "values_v":vv, "grad_v":vg,
        "values_b":bv, "grad_b":bg, "disp":disp, "eps":eps, "sig":sig,
        "moments":moments, "traction":traction, "wrench":wrench,
        "condition_max":max(cv,cb), "volume_m3":float(vmeasure),
        "numeric_array_bytes":sum(a.nbytes for a in (k,c,f,g,u,lam,vv,vg,bv,bg,disp,eps,sig,traction))}


def capture_arrays(np, deps, case, inputs, root_token=None):
    """Assemble fresh baseline/rotated states from caller-verified frozen arrays."""
    require_root_go(case, root_token)
    inputs = dict(inputs)
    ids = tuple(inputs["source_ids"])
    if len(ids) != 27 or ids != tuple(sorted(ids)):
        raise AdapterBlocked("source IDs must be canonical 27-site order")
    required = {"source_positions_xyz_m":(27,3), "volume_xyz_m":(512,3),
        "volume_weights_m3":(512,), "boundary_xyz_m":(768,3),
        "boundary_weights_m2":(768,), "boundary_normals":(768,3)}
    for name, shape in required.items():
        value = np.asarray(inputs[name], dtype=np.float64)
        if value.shape != shape or not bool(np.isfinite(value).all()):
            raise AdapterBlocked(f"invalid frozen input {name}")
        inputs[name] = value
    if bool(np.any(inputs["volume_weights_m3"] <= 0)) or bool(np.any(inputs["boundary_weights_m2"] <= 0)):
        raise AdapterBlocked("quadrature measures must be positive")
    if float(np.max(np.abs(np.linalg.norm(inputs["boundary_normals"], axis=1)-1.0))) > 1e-10:
        raise AdapterBlocked("boundary normals must remain unit vectors")
    baseline = case["baseline_crack_frame"]
    q = rotation_matrix(np, case["rigid_transform"])
    t = _vector(np, case["rigid_transform"]["translation_m"], "translation")
    z = np.zeros(3)
    x0 = _state(np,deps,case,inputs,ids,np.eye(3),z,baseline)
    x1 = _state(np,deps,case,inputs,ids,q,t,baseline)
    a = np.kron(np.eye(31),q); b = np.zeros((6,6)); b[:3,:3]=q; b[3:,3:]=q
    expected_k, expected_c = a@x0["K"]@a.T, b@x0["C"]@a.T
    expected_f, expected_u, expected_lam = a@x0["f"], a@x0["q"], b@x0["lambda"]
    g0 = np.concatenate((x0["grad_v"],x0["grad_b"]))
    g1 = np.concatenate((x1["grad_v"],x1["grad_b"]))
    ge = np.einsum("ij,nkj->nki",q,g0)
    strain0, strain1 = _strain_tensor(np,x0["eps"]), _strain_tensor(np,x1["eps"])
    stress0, stress1 = _stress_tensor(np,x0["sig"]), _stress_tensor(np,x1["sig"])
    strain_e = np.einsum("ij,njk,lk->nil",q,strain0,q)
    stress_e = np.einsum("ij,njk,lk->nil",q,stress0,q)
    u_e = x0["disp"]@q.T
    force_e=q@x0["wrench"][0]; mc_e=q@x0["wrench"][1]
    m0_e=q@x0["wrench"][2]+np.cross(t,force_e)
    fscale=max(x0["wrench"][3],np.finfo(float).tiny)
    mscale=float(case["geometry"]["radius_m"])*fscale
    errors={"K":_rel(np,x1["K"],expected_k),"C":_rel(np,x1["C"],expected_c),
        "f":_rel(np,x1["f"],expected_f),"q":_rel(np,x1["q"],expected_u),
        "lambda_scaled_by_force":float(np.linalg.norm(x1["lambda"]-expected_lam)/max(float(np.linalg.norm(x0["f"])),np.finfo(float).tiny)),
        "basis_values_max_abs":float(np.max(np.abs(np.concatenate((x1["values_v"],x1["values_b"]))-np.concatenate((x0["values_v"],x0["values_b"]))))),
        "basis_gradient_scale_max_abs":float(np.max(np.abs(g1-ge))*case["mls"]["reference_scale_m"]),
        "volume_displacement":_weighted(np,x1["disp"],u_e,inputs["volume_weights_m3"]),
        "volume_strain_tensor":_weighted(np,strain1,strain_e,inputs["volume_weights_m3"]),
        "volume_stress_tensor":_weighted(np,stress1,stress_e,inputs["volume_weights_m3"]),
        "force_wrench":float(np.linalg.norm(x1["wrench"][0]-force_e)/fscale),
        "centered_moment_wrench":float(np.linalg.norm(x1["wrench"][1]-mc_e)/mscale),
        "origin_moment_translation_law":float(np.linalg.norm(x1["wrench"][2]-m0_e)/mscale),
        "direct_physical_moment_scaled_by_R":float(np.linalg.norm(x1["moments"]-x1["C"]@x1["q"])/case["geometry"]["radius_m"]),
        "raw_K_asymmetry":max(_rel(np,x0["K"],x0["K"].T),_rel(np,x1["K"],x1["K"].T)),
        "twice_energy":_rel(np,np.asarray([x1["solve"]["twice_strain_energy_J"]]),np.asarray([x0["solve"]["twice_strain_energy_J"]])),
        "load_work":_rel(np,np.asarray([x1["solve"]["load_work_J"]]),np.asarray([x0["solve"]["load_work_J"]])),
        "energy_identity":max(x0["solve"]["energy_identity_relative_error"],x1["solve"]["energy_identity_relative_error"])}
    tol=case["tolerances"]
    passed=(errors["K"]<=tol["K_C_relative"] and errors["C"]<=tol["K_C_relative"]
        and errors["f"]<=tol["f_relative_L2"] and errors["q"]<=tol["q_relative_L2"]
        and errors["lambda_scaled_by_force"]<=tol["lambda_scaled_by_force"]
        and errors["basis_values_max_abs"]<=tol["provider_weights"]
        and errors["basis_gradient_scale_max_abs"]<=tol["provider_gradient_scale"]
        and max(errors["volume_displacement"],errors["volume_strain_tensor"],errors["volume_stress_tensor"])<=tol["volume_fields_relative"]
        and max(errors["force_wrench"],errors["centered_moment_wrench"],errors["origin_moment_translation_law"])<=tol["wrench_scaled"]
        and errors["direct_physical_moment_scaled_by_R"]<=tol["gauge_residual"]
        and errors["raw_K_asymmetry"]<=tol["raw_K_relative_asymmetry"]
        and max(errors["twice_energy"],errors["load_work"])<=tol["energy_invariance"]
        and errors["energy_identity"]<=tol["energy_identity"])
    record = {"schema":"mixed93_frame_aware_objectivity_capture.v1",
        "classification":"ACTUAL_MECHANICS_COORDINATE_COVARIANCE_ONLY_NO_RELATIVE_CLOUD_QUADRATURE_TETRA_OR_CRACK_BIAS",
        "status":"PASS_COORDINATE_COVARIANCE_DIAGNOSTIC_ONLY" if passed else "FAIL_COORDINATE_COVARIANCE_DIAGNOSTIC",
        "provider_queries":2560,"analytic_boundary_field_evaluations":1536,
        "raw_assemblies":2,"fresh_spectrum_checks":2,"KKT_solves_99x99":2,
        "Q":q.tolist(),"translation_m":t.tolist(),"errors":errors,
        "baseline":{"K_N_per_m":x0["K"].tolist(),"C":x0["C"].tolist(),"f_N":x0["f"].tolist(),"q_m":x0["q"].tolist(),"lambda_N":x0["lambda"].tolist(),"g_m":x0["g"].tolist(),"spectrum":x0["spectrum"],"solve_metrics":x0["solve"],"volume_moments_m":x0["moments"].tolist(),"volume_displacement_m":x0["disp"].tolist(),"volume_strain_engineering":x0["eps"].tolist(),"volume_stress_Pa":x0["sig"].tolist(),"boundary_traction_Pa":x0["traction"].tolist(),"volume_measure_m3":x0["volume_m3"],"condition_number_max":x0["condition_max"]},
        "rotated":{"K_N_per_m":x1["K"].tolist(),"C":x1["C"].tolist(),"f_N":x1["f"].tolist(),"q_m":x1["q"].tolist(),"lambda_N":x1["lambda"].tolist(),"g_m":x1["g"].tolist(),"spectrum":x1["spectrum"],"solve_metrics":x1["solve"],"volume_moments_m":x1["moments"].tolist(),"volume_displacement_m":x1["disp"].tolist(),"volume_strain_engineering":x1["eps"].tolist(),"volume_stress_Pa":x1["sig"].tolist(),"boundary_traction_Pa":x1["traction"].tolist(),"volume_measure_m3":x1["volume_m3"],"condition_number_max":x1["condition_max"]},
        "frame_vectors":{"baseline":{"body_axis":x0["frame"].body_axis_world.tolist(),"front_direction":x0["frame"].front_direction_world.tolist(),"crack_plane_normal":x0["frame"].crack_plane_normal_world.tolist(),"crack_axes":x0["frame"].axes.tolist(),"crack_tip_m":x0["frame"].tip.tolist(),"MLS_origin_m":x0["frame"].mls_origin_world_m.tolist(),"moment_center_m":x0["frame"].center.tolist()},"rotated":{"body_axis":x1["frame"].body_axis_world.tolist(),"front_direction":x1["frame"].front_direction_world.tolist(),"crack_plane_normal":x1["frame"].crack_plane_normal_world.tolist(),"crack_axes":x1["frame"].axes.tolist(),"crack_tip_m":x1["frame"].tip.tolist(),"MLS_origin_m":x1["frame"].mls_origin_world_m.tolist(),"moment_center_m":x1["frame"].center.tolist()}},
        "wrenches":{"baseline":{"force_N":x0["wrench"][0].tolist(),"centered_moment_Nm":x0["wrench"][1].tolist(),"origin_moment_Nm":x0["wrench"][2].tolist(),"absolute_traction_integral_N":x0["wrench"][3]},"rotated":{"force_N":x1["wrench"][0].tolist(),"centered_moment_Nm":x1["wrench"][1].tolist(),"origin_moment_Nm":x1["wrench"][2].tolist(),"absolute_traction_integral_N":x1["wrench"][3]},"moment_center_transform":"c'=Qc+t; M_c'=Q M_c; M_0'=Q M_0+t cross F'","force_scale":"integral |traction| dA","moment_scale":"R times integral |traction| dA"},
        "limits":["toy manufactured isotropic mechanics only","co-rotated cloud and quadrature: no relative orientation bias","no material or crack-growth qualification"]}
    record["explicit_numeric_array_bytes_both_states"] = x0["numeric_array_bytes"] + x1["numeric_array_bytes"]
    encoded_size = len(json.dumps(record, allow_nan=False, separators=(",", ":")).encode("utf-8"))
    record["encoded_result_bytes"] = encoded_size
    if encoded_size > int(case["resource_caps_provisional"]["result_bytes"]):
        raise AdapterBlocked("serialized record exceeds provisional output cap")
    if record["explicit_numeric_array_bytes_both_states"] > int(case["resource_caps_provisional"]["explicit_dense_array_bytes"]):
        raise AdapterBlocked("explicit numeric arrays exceed provisional memory cap")
    return record
