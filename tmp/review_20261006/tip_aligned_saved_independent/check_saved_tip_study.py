#!/usr/bin/env python3
"""Independent saved-array audit; no production imports, queries, or solves."""
import argparse, ast, hashlib, json, math, struct, time, zipfile
from pathlib import Path
np = None
CASE_SHA = "9d1b2153a1a3db36387d838f584714dae8595c722f542155e83bfe1f037e1c5d"
RESULT_SHA = "2510eabd3e3406f243d9fe900ac5af3cfdf17071cdbe50383093de0a13adff4b"
NPZ_SHA = "a9de2e3eaf88ba0eb130fca6f160d5d327a963bb00203266b7390cadbe6693ea"
R, L, E, NU, KIII, A_SHIFT = 0.01, 0.01, 1.0e6, 0.3, 100.0, 0.0005
CASES = [("a0_theta32", 0.0, 32), ("a0_theta64", 0.0, 64),
         ("a0.0005_theta32", A_SHIFT, 32), ("a0.0005_theta64", A_SHIFT, 64)]
ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[3]
DRIVER = REPO / "pgfrag/experiments/tip_aligned_quadrature_study_20261006/tip_aligned_quadrature_study.py"
CASE = DRIVER.parent / "study_case.json"
DEFAULT_RESULT = REPO / "pgfrag/tmp/review_20261006/tip_aligned_quadrature_study_20261006/run01/result.json"

def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()
def npy_shape(z, key):
    with z.open(key + ".npy") as f:
        if f.read(6) != b"\x93NUMPY": raise ValueError("invalid NPY member: " + key)
        major, minor = f.read(2); n = struct.unpack("<H" if major == 1 else "<I", f.read(2 if major == 1 else 4))[0]
        if major not in (1,2,3): raise ValueError("unsupported NPY version: " + repr((major,minor)))
        header = ast.literal_eval(f.read(n).decode("latin1").strip())
        if header["fortran_order"] or header["descr"].endswith("O"):
            raise ValueError("fortran/object array is not allowed: " + key)
        return tuple(header["shape"])
def preflight_npz(path, result):
    n = int(result["query_budget"]["unique_xyz_rows"]); shapes = {
        "source_positions_xyz_m":(27,3), "unique_query_positions_xyz_m":(n,3),
        "base_mls_values":(n,27), "base_mls_gradients_per_m":(n,27,3), "exact_stress_Pa":(n,3,3),
        "rigid_mode_coefficient_vectors":(162,6), "gauge_target_m":(6,),
        "registry_ordered_keys_json":(162,), "basis_row_ids_json":(54,),
        "cloud_side_xyz_m":(256,3), "cloud_side_weights":(256,), "cloud_side_normals":(256,3),
        "cloud_side_unique_row_ids":(256,), "cloud_offgrid_xyz_m":(128,3), "cloud_offgrid_unique_row_ids":(128,)}
    for a in (0.0,A_SHIFT):
        tag=f"a{a:g}".replace(".","_")
        for suffix,count in (("cap8",256),("cap16",512)):
            k=f"cloud_{tag}_{suffix}_"; shapes.update({k+"xyz_m":(count,3),k+"weights":(count,),k+"normals":(count,3),k+"unique_row_ids":(count,)})
        for theta in (32,64):
            k=f"cloud_{tag}_volume{theta}_"; count=512 if theta==32 else 1024
            shapes.update({k+"xyz_m":(count,3),k+"weights":(count,),k+"unique_row_ids":(count,)})
    for name,a,theta in CASES:
        k=name+"__"; nv=512 if theta==32 else 1024
        shapes.update({k+"volume_xyz_m":(nv,3),k+"volume_weights_m3":(nv,),k+"volume_basis_values":(nv,54),k+"volume_basis_gradients_per_m":(nv,54,3),k+"K_N_per_m":(162,162),k+"boundary_load_N":(162,),k+"constraint_matrix":(6,162),k+"rigid_mode_coefficients":(162,6),k+"coefficients_m":(162,),k+"multipliers":(6,),k+"offgrid_xyz_m":(128,3),k+"offgrid_basis_values":(128,54),k+"offgrid_basis_gradients_per_m":(128,54,3),k+"solved_offgrid_gradient_m_per_m":(128,3,3),k+"solved_offgrid_stress_Pa":(128,3,3)})
        for suffix,count in (("side",256),("cap8",256),("cap16",512)):
            shapes.update({k+suffix+"_basis_values":(count,54),k+suffix+"_basis_gradients_per_m":(count,54,3),k+suffix+"_traction_Pa":(count,3)})
        shapes.update({k+"side_load_N":(162,),k+"cap16_load_N":(162,)})
    with zipfile.ZipFile(path) as z:
        entries={x[:-4] for x in z.namelist() if x.endswith(".npy")}; missing=set(shapes)-entries
        if missing: raise ValueError("NPZ members missing runtime keys: "+repr(sorted(missing)))
        bad={key:(shape, got) for key,shape in shapes.items() if (got:=npy_shape(z,key))!=shape}
        if bad: raise ValueError("NPZ member shape mismatch: "+repr(bad))
    if set(result["cases"]) != {name for name,_,_ in CASES} or any(not {"solve"} <= set(result["cases"][name]) or not {"strain_energy_J","external_work_J"} <= set(result["cases"][name]["solve"]) for name,_,_ in CASES) or any(not {"Pi_a0_J","Pi_a0.5mm_J","G_delta_a_J_per_m2"} <= set(result["cross_case_finite_delta_a_energy"][f"theta{t}"]) for t in (32,64)):
        raise ValueError("JSON case/cross-case keys differ from the expected result schema")
    return {"npy_member_count":len(entries),"runtime_key_count":len(shapes),"unique_xyz_rows":n}
def rel(a, b): return float(np.linalg.norm(a - b) / max(float(np.linalg.norm(b)), np.finfo(float).tiny))
def make_D():
    mu = E / (2.0 * (1.0 + NU)); lam = E * NU / ((1.0 + NU) * (1.0 - 2.0 * NU))
    D = np.zeros((6, 6))
    D[:3, :3] = lam; np.fill_diagonal(D[:3, :3], lam + 2.0 * mu); D[3:, 3:] = np.eye(3) * mu
    return D
def shifted_rows(base, base_g, xyz, a):
    dx, y = xyz[:, 0] - a, xyz[:, 1]
    rho = np.hypot(dx, y)
    if np.any(rho <= 0.0) or np.any((y == 0.0) & (dx < 0.0)):
        raise ValueError("cached F1 query is singular or on its undeclared branch cut")
    theta, root = np.arctan2(y, dx), np.sqrt(rho)
    f = root * np.sin(theta / 2.0) / math.sqrt(R)
    gf = np.column_stack((-np.sin(theta / 2.0) / (2.0 * root * math.sqrt(R)),
                          np.cos(theta / 2.0) / (2.0 * root * math.sqrt(R)),
                          np.zeros(len(xyz))))
    val, grad = np.empty((len(base), 54)), np.empty((len(base), 54, 3))
    val[:, 0::2], grad[:, 0::2] = base, base_g
    val[:, 1::2] = base * f[:, None]
    grad[:, 1::2] = base_g * f[:, None, None] + base[:, :, None] * gf[:, None, :]
    return val, grad
def dof_maps(arr, result):
    keys = [json.loads(str(x)) for x in arr["registry_ordered_keys_json"]]; basis = [tuple(json.loads(str(x))) for x in arr["basis_row_ids_json"]]
    if len(keys) != 162 or len(set(tuple(x) for x in keys)) != 162 or len(basis) != 54:
        raise ValueError("semantic registry must contain 162 unique DOFs and 54 basis rows")
    bmap = {key: i for i, key in enumerate(basis)}
    if len(bmap) != 54:
        raise ValueError("basis row IDs are not unique")
    comp = {"ux": 0, "uy": 1, "uz": 2}; slots = []
    for sid, field, kind, _basis_index in keys:
        if field not in comp or (sid, kind) not in bmap:
            raise ValueError("semantic key does not map to a known field/basis row")
        slots.append((bmap[(sid, kind)], comp[field]))
    registry = result["registry"]
    expected = [(sid, kind) for sid in registry["source_ids"] for kind in ("base", "front_f1_unshifted")]
    if basis != expected:
        raise ValueError("basis rows must interleave each source as BASE then FRONT_F1_UNSHIFTED")
    if keys != registry["ordered_keys"] or basis != [tuple(x) for x in registry["basis_row_ids"]]:
        raise ValueError("cached semantic registry differs from result metadata")
    return keys, slots
def independent_rigid(source_xyz, source_ids, keys):
    center = source_xyz.mean(axis=0); axes = np.eye(3)
    index = {sid: i for i, sid in enumerate(source_ids)}; out = np.zeros((162, 6))
    for j, (sid, field, basis, _bi) in enumerate(keys):
        if basis != "base":
            continue
        c, p = {"ux": 0, "uy": 1, "uz": 2}[field], source_xyz[index[sid]] - center
        out[j, c] = 1.0
        for k in range(3):
            out[j, 3 + k] = np.cross(axes[k], p)[c]
    return out
def rebuild_K(grads, weights, slots, D):
    K = np.zeros((162, 162))
    for lo in range(0, len(weights), 128):
        hi = min(lo + 128, len(weights)); g, w = grads[lo:hi], weights[lo:hi]
        B = np.zeros((len(w), 6, 162))
        for j, (row, c) in enumerate(slots):
            gx, gy, gz = g[:, row, 0], g[:, row, 1], g[:, row, 2]
            if c == 0:
                B[:, 0, j], B[:, 3, j], B[:, 5, j] = gx, gy, gz
            elif c == 1:
                B[:, 3, j], B[:, 1, j], B[:, 4, j] = gx, gy, gz
            else:
                B[:, 5, j], B[:, 4, j], B[:, 2, j] = gx, gy, gz
        DB = np.einsum("ij,qja->qia", D, B, optimize=True); K += np.einsum("qia,qib,q->ab", B, DB, w, optimize=True)
    return K
def rebuild_C(values, xyz, weights, slots):
    volume = float(weights.sum()); centered = xyz - np.einsum("q,qc->c", weights, xyz) / volume
    modes = np.zeros((len(xyz), 6, 3)); C = np.zeros((6, 162))
    modes[:, 0, 0] = modes[:, 1, 1] = modes[:, 2, 2] = 1.0
    for k, axis in enumerate(np.eye(3)):
        modes[:, 3 + k, :] = np.cross(axis, centered) / R
    for j, (row, c) in enumerate(slots):
        C[:, j] = np.einsum("q,q,qm->m", weights, values[:, row], modes[:, :, c]) / volume
    return C
def load_vector(values, traction, weights, slots):
    return np.array([np.dot(weights * traction[:, c], values[:, row]) for row, c in slots])

def audit(args):
    global np
    t0 = time.monotonic(); rp = Path(args.result).resolve(); result_hash_pre = sha(rp)
    if result_hash_pre.lower() != RESULT_SHA: raise ValueError("primary JSON differs from the pinned run01 result")
    driver_hash_pre, case_hash_pre = sha(DRIVER), sha(CASE); result = json.loads(rp.read_text(encoding="utf-8"))
    if result.get("schema") != "tip_aligned_full3d_quadrature_study_result.v1":
        raise ValueError("unexpected result schema")
    if case_hash_pre.lower() != CASE_SHA or result.get("study_case_sha256", "").lower() != case_hash_pre.lower():
        raise ValueError("result was not produced from the reviewed study case")
    if result.get("source_sha256", {}).get("driver", "").lower() != driver_hash_pre.lower():
        raise ValueError("result driver pin differs from reviewed source")
    if result.get("status") != "OK_DIAGNOSTIC_ONLY_UNQUALIFIED":
        raise ValueError("primary result status is not the declared diagnostic status")
    npz_path = rp.parent / result["arrays_npz"]["path"]; npz_hash_pre = sha(npz_path)
    if npz_hash_pre.lower() != NPZ_SHA or npz_hash_pre.lower() != result["arrays_npz"]["sha256"].lower():
        raise ValueError("cached NPZ SHA-256 mismatch")
    if npz_path.stat().st_size != result["arrays_npz"]["bytes"]:
        raise ValueError("cached NPZ byte count mismatch")
    if not (ROOT == Path(args.out).resolve() or ROOT in Path(args.out).resolve().parents):
        raise ValueError("output must remain inside this independent-review directory")
    outdir = Path(args.out).resolve()
    planned=["audit.json","a0_offgrid_stress_scatter.png","finite_G_theta_ledger.png"]
    if any((outdir/n).exists() for n in planned): raise FileExistsError("reviewer output exists; refusing overwrite")
    outdir.mkdir(parents=True, exist_ok=True)
    schema = preflight_npz(npz_path,result); import numpy as np
    with np.load(npz_path, allow_pickle=False) as z:
        a = {k: z[k] for k in z.files}
    for name, x in a.items():
        if x.dtype.kind in "fiu" and not np.isfinite(x).all():
            raise ValueError("non-finite cached numeric array: " + name)
    keys, slots = dof_maps(a, result); D = make_D()
    unique = a["unique_query_positions_xyz_m"]; exactS = a["exact_stress_Pa"]; source = a["source_positions_xyz_m"]
    source_ids = result["registry"]["source_ids"]
    if source.shape != (27, 3) or len(source_ids) != 27 or exactS.shape != (len(unique), 3, 3):
        raise ValueError("source or reference array dimensions violate the frozen semantic schema")
    rigid = independent_rigid(source, source_ids, keys); rigid_cache = a["rigid_mode_coefficient_vectors"]
    rigid_rel = rel(rigid, rigid_cache); basis_checks, algebra, offgrid = {}, {}, {}
    def point(tag):
        p = "cloud_" + tag + "_"
        return (a[p + "xyz_m"], a.get(p + "weights"), a.get(p + "normals"), a[p + "unique_row_ids"])
    for name, trial_a, theta in CASES:
        p = name + "__"
        vxyz, vw = a[p + "volume_xyz_m"], a[p + "volume_weights_m3"]; vg, vv = a[p + "volume_basis_gradients_per_m"], a[p + "volume_basis_values"]
        K = rebuild_K(vg, vw, slots, D); C = rebuild_C(vv, vxyz, vw, slots)
        tagvol = f"a{trial_a:g}_volume{theta}".replace(".", "_")
        checks = {}
        for label, tag, saved_v, saved_g in (
            ("volume", tagvol, vv, vg), ("offgrid", "offgrid", a[p+"offgrid_basis_values"], a[p+"offgrid_basis_gradients_per_m"]),
            ("side", "side", a[p+"side_basis_values"], a[p+"side_basis_gradients_per_m"]),
            ("cap8", f"a{trial_a:g}_cap8".replace(".", "_"), a[p+"cap8_basis_values"], a[p+"cap8_basis_gradients_per_m"]),
            ("cap16", f"a{trial_a:g}_cap16".replace(".", "_"), a[p+"cap16_basis_values"], a[p+"cap16_basis_gradients_per_m"])):
            xyz, _, _, ids = point(tag)
            if label == "volume" and not np.array_equal(vxyz, xyz):
                raise ValueError(f"{name} {label} coordinates differ from shared cloud")
            ev, eg = shifted_rows(a["base_mls_values"][ids], a["base_mls_gradients_per_m"][ids], xyz, trial_a)
            checks[label] = {"value_rel": rel(saved_v, ev), "gradient_rel": rel(saved_g, eg),
                             "base_value_max_abs": float(np.max(np.abs(saved_v[:, 0::2]-a["base_mls_values"][ids]))),
                             "base_gradient_max_abs": float(np.max(np.abs(saved_g[:, 0::2]-a["base_mls_gradients_per_m"][ids])))}
        basis_checks[name] = checks
        # Independently reconstruct old-tip prescribed traction and generalized f.
        _side_xyz, sw, sn, si = point("side")
        cap_tag = f"a{trial_a:g}_cap16".replace(".", "_")
        _cap_xyz, cw, cn, ci = point(cap_tag)
        st = np.einsum("qij,qj->qi", exactS[si], sn); ct = np.einsum("qij,qj->qi", exactS[ci], cn)
        if rel(st, a[p+"side_traction_Pa"]) > 1e-10 or rel(ct, a[p+"cap16_traction_Pa"]) > 1e-10:
            raise ValueError(f"{name} saved boundary traction differs from old-tip stress dot normal")
        fs = load_vector(a[p+"side_basis_values"], st, sw, slots); fc = load_vector(a[p+"cap16_basis_values"], ct, cw, slots)
        f = fs + fc
        K0, C0, f0 = a[p+"K_N_per_m"], a[p+"constraint_matrix"], a[p+"boundary_load_N"]
        d, m, g = a[p+"coefficients_m"], a[p+"multipliers"], a["gauge_target_m"]
        diagonal = np.diag(K)
        scale = 1.0 / np.sqrt(diagonal); Ks = scale[:, None] * K * scale[None, :]
        eig = np.linalg.eigvalsh(0.5 * (Ks + Ks.T)); tol = 64 * np.finfo(float).eps * 162 * max(float(np.max(np.abs(eig))), np.finfo(float).tiny)
        rigid_scaled = rigid / scale[:, None]
        rigid_res = np.linalg.norm(Ks @ rigid_scaled, axis=0) / np.maximum(
            max(float(np.max(np.abs(eig))), np.finfo(float).tiny) * np.linalg.norm(rigid_scaled, axis=0), np.finfo(float).tiny)
        eq = K @ d + C.T @ m - f; cr = C @ d - g
        saddle = np.block([[K, C.T], [C, np.zeros((6, 6))]])
        rhs = np.r_[f, g]; sol = np.r_[d, m]
        U, W = 0.5 * float(d @ K @ d), float(d @ f); ident = 2.0 * U - W + float(g @ m)
        og = a[p+"offgrid_basis_gradients_per_m"]; grad = np.zeros((len(og),3,3))
        for j,(row,c) in enumerate(slots): grad[:,c,:] += d[j]*og[:,row,:]
        strain = np.column_stack((grad[:,0,0],grad[:,1,1],grad[:,2,2],grad[:,0,1]+grad[:,1,0],grad[:,1,2]+grad[:,2,1],grad[:,0,2]+grad[:,2,0]))
        sv = strain @ D.T; so = a[p+"solved_offgrid_stress_Pa"]
        ss = np.column_stack((so[:,0,0],so[:,1,1],so[:,2,2],so[:,0,1],so[:,1,2],so[:,0,2]))
        raw_solve = result["cases"][name]["solve"]
        algebra[name] = {
            "K_rebuild_rel": rel(K, K0), "C_rebuild_rel": rel(C, C0),
            "f_rebuild_rel": rel(f, f0), "side_f_rel": rel(fs, a[p+"side_load_N"]),
            "cap16_f_rel": rel(fc, a[p+"cap16_load_N"]), "rigid_coeff_rel": rigid_rel,
            "scaled_eigen_counts": {"near_zero": int(np.count_nonzero(np.abs(eig) <= tol)),
                                     "positive": int(np.count_nonzero(eig > tol)),
                                     "negative": int(np.count_nonzero(eig < -tol)), "zero_tolerance": float(tol)},
            "six_scaled_rigid_residuals": rigid_res.tolist(),
            "KKT168_backward_residual": float(np.linalg.norm(saddle @ sol-rhs)/max(np.linalg.norm(rhs),np.finfo(float).tiny)),
            "equilibrium_rel": float(np.linalg.norm(eq)/max(np.linalg.norm(K@d)+np.linalg.norm(C.T@m)+np.linalg.norm(f),np.finfo(float).tiny)),
            "constraint_rel": float(np.linalg.norm(cr)/max(np.linalg.norm(g),np.linalg.norm(C@d),np.finfo(float).tiny)),
            "U_J": U, "W_J": W, "Pi_J": U-W, "identity_defect_J": ident,
            "identity_rel": abs(ident)/max(abs(2*U),abs(W),abs(float(g@m)),np.finfo(float).tiny),
            "saved_U_rel_delta": abs(U-raw_solve["strain_energy_J"])/max(abs(raw_solve["strain_energy_J"]),np.finfo(float).tiny),
            "saved_W_rel_delta": abs(W-raw_solve["external_work_J"])/max(abs(raw_solve["external_work_J"]),np.finfo(float).tiny),
            "coeff_rebuilt_offgrid_gradient_rel": rel(grad,a[p+"solved_offgrid_gradient_m_per_m"]),
            "coeff_rebuilt_offgrid_stress_rel": rel(sv,ss),
        }
        oi = point("offgrid")[3]; ox = a[p+"offgrid_xyz_m"]
        if not np.array_equal(ox, unique[oi]):
            raise ValueError("offgrid points do not map to their shared physical coordinates")
        if trial_a == 0.0:
            ex = exactS[oi]
            so = a[p+"solved_offgrid_stress_Pa"]
            normal = (so-ex)[:, (0,1,2), (0,1,2)]
            offgrid[name] = {"shared_rows": int(len(oi)), "full_stress_rel": rel(so, ex),
                             "normal_stress_error_over_exact_full_norm": float(np.linalg.norm(normal)/max(np.linalg.norm(ex),np.finfo(float).tiny)),
                             "max_abs_sigma_yy_Pa": float(np.max(np.abs(so[:,1,1]))),
                             "max_abs_normal_component_error_Pa": float(np.max(np.abs(normal)))}
    # Confirm fixed-load point clouds and prescribed tractions across a and theta.
    fixed_load = {}
    for label, tags in (("side", ["side"]*4), ("cap8", ["a0_cap8","a0_0005_cap8"]*2), ("cap16", ["a0_cap16","a0_0005_cap16"]*2)):
        refs = [point(tag) for tag in tags]
        fixed_load[label] = {"cloud_xyz_weights_normals_identical": all(
            all(np.array_equal(refs[0][j], r[j]) for j in range(3)) for r in refs[1:]),
            "stored_traction_max_delta_Pa": max(float(np.max(np.abs(a[c+"__"+label+"_traction_Pa"]-a[CASES[0][0]+"__"+label+"_traction_Pa"])))
                for c, _, _ in CASES)}
    # Cross-case physical offgrid field deltas and finite delta-a G from cached solved states.
    ledger = {}
    for theta in (32, 64):
        n0, n5 = f"a0_theta{theta}", f"a0.0005_theta{theta}"; p0, p5 = n0+"__", n5+"__"
        pi0, pi5 = algebra[n0]["Pi_J"], algebra[n5]["Pi_J"]
        ledger[str(theta)] = {"Pi_a0_J": pi0, "Pi_a0.5mm_J": pi5,
                              "G_delta_a_J_per_m2": -(pi5-pi0)/(L*A_SHIFT),
                              "theta_fixed_offgrid_stress_change_rel": rel(a[p5+"solved_offgrid_stress_Pa"],a[p0+"solved_offgrid_stress_Pa"])}
    for x in (0.0, A_SHIFT):
        n32, n64 = f"a{x:g}_theta32", f"a{x:g}_theta64"
        ledger[f"theta_change_a{x:g}"] = {"Pi64_minus_Pi32_J": algebra[n64]["Pi_J"]-algebra[n32]["Pi_J"],
            "raw_K_theta_rel": rel(a[n64+"__K_N_per_m"], a[n32+"__K_N_per_m"]),
            "fixed_f_theta_rel": rel(a[n64+"__boundary_load_N"],a[n32+"__boundary_load_N"])}
    primary_g = result["cross_case_finite_delta_a_energy"]
    for theta in (32,64):
        q=primary_g[f"theta{theta}"]; v=ledger[str(theta)]
        v["primary_Pi0_rel_delta"]=abs(v["Pi_a0_J"]-q["Pi_a0_J"])/max(abs(q["Pi_a0_J"]),np.finfo(float).tiny); v["primary_PiShift_rel_delta"]=abs(v["Pi_a0.5mm_J"]-q["Pi_a0.5mm_J"])/max(abs(q["Pi_a0.5mm_J"]),np.finfo(float).tiny); v["primary_G_rel_delta"]=abs(v["G_delta_a_J_per_m2"]-q["G_delta_a_J_per_m2"])/max(abs(q["G_delta_a_J_per_m2"]),np.finfo(float).tiny)
    # Consistency flags cover stored-array math only; physical error direction is never gated.
    algebra_flags = {n: (v["K_rebuild_rel"] < 1e-8 and v["C_rebuild_rel"] < 1e-10 and v["rigid_coeff_rel"] < 1e-12 and max(v["six_scaled_rigid_residuals"]) < 128*np.finfo(float).eps*162 and
                         v["f_rebuild_rel"] < 1e-10 and v["KKT168_backward_residual"] < 1e-8 and v["equilibrium_rel"] < 1e-8 and v["constraint_rel"] < 1e-8 and v["saved_U_rel_delta"] < 1e-8 and v["saved_W_rel_delta"] < 1e-8 and
                         v["identity_rel"] < 1e-8 and v["coeff_rebuilt_offgrid_gradient_rel"] < 1e-10 and v["coeff_rebuilt_offgrid_stress_rel"] < 1e-10 and v["scaled_eigen_counts"] ==
                         {"near_zero": 6, "positive": 156, "negative": 0, "zero_tolerance": v["scaled_eigen_counts"]["zero_tolerance"]})
                     for n,v in algebra.items()}
    basis_ok = all(max(q["value_rel"],q["gradient_rel"]) < 1e-10 for c in basis_checks.values() for q in c.values())
    fixed_load_ok = all(v["cloud_xyz_weights_normals_identical"] and v["stored_traction_max_delta_Pa"] < 1e-10
                        for v in fixed_load.values())
    ledger_ok = all(max(ledger[str(t)]["primary_Pi0_rel_delta"],ledger[str(t)]["primary_PiShift_rel_delta"]) < 1e-8 and ledger[str(t)]["primary_G_rel_delta"] < 1e-6 for t in (32,64))
    if not basis_ok or not all(algebra_flags.values()) or not fixed_load_ok or not ledger_ok:
        audit_status = "BLOCKED_CACHED_NUMERIC_OR_STRUCTURAL_MISMATCH"
    else:
        audit_status = "CACHED_ARRAYS_INDEPENDENTLY_CONSISTENT_PHYSICS_UNQUALIFIED"
    # Figures use only persisted points, stresses, and energies.
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.7))
    off_id = point("offgrid")[3]
    exact_off = exactS[off_id]
    for ax, (i,j,title) in zip(axes, ((0,2,r"$\sigma_{xz}$ (Pa)"),(1,2,r"$\sigma_{yz}$ (Pa)"),(1,1,r"$\sigma_{yy}$ (Pa)"))):
        for name, color, marker in (("a0_theta32","tab:blue","o"),("a0_theta64","tab:orange","x")):
            pred = a[name+"__solved_offgrid_stress_Pa"][:,i,j]
            ax.scatter(exact_off[:,i,j],pred,s=15,alpha=.7,label=name.split("_")[-1],c=color,marker=marker)
        ax.set_title(title); ax.set_xlabel("old-tip manufactured exact"); ax.set_ylabel("saved solved"); ax.grid(alpha=.2)
    axes[0].legend(title="a=0, theta"); fig.tight_layout(); fig.savefig(outdir/"a0_offgrid_stress_scatter.png",dpi=140); plt.close(fig)
    fig, axes = plt.subplots(1,3,figsize=(13,3.7))
    axes[0].plot([32,64],[ledger[str(t)]["G_delta_a_J_per_m2"] for t in (32,64)],"o-")
    axes[0].set(xlabel="volume theta order",ylabel="finite-delta-a G (J/m²)",title="fixed-load finite offset")
    axes[1].plot([32,64],[ledger[str(t)]["Pi_a0_J"] for t in (32,64)],"o-",label="Pi(a=0)")
    axes[1].plot([32,64],[ledger[str(t)]["Pi_a0.5mm_J"] for t in (32,64)],"x-",label="Pi(a=0.5mm)")
    axes[1].legend(); axes[1].set(xlabel="volume theta order",ylabel="Pi (J)",title="fixed-theta potentials")
    names=["a0","a0.5mm"]
    vals=[ledger["theta_change_a0"]["Pi64_minus_Pi32_J"],ledger["theta_change_a0.0005"]["Pi64_minus_Pi32_J"]]
    axes[2].bar(names,vals); axes[2].set(ylabel="Pi64 − Pi32 (J)",title="theta-change ledger")
    fig.tight_layout(); fig.savefig(outdir/"finite_G_theta_ledger.png",dpi=140); plt.close(fig)
    result_hash_post, npz_hash_post = sha(rp), sha(npz_path); driver_hash_post, case_hash_post = sha(DRIVER), sha(CASE)
    if result_hash_pre != result_hash_post or npz_hash_pre != npz_hash_post or driver_hash_pre != driver_hash_post or case_hash_pre != case_hash_post:
        raise ValueError("primary JSON or NPZ changed during cached-array audit")
    doc = {"schema":"tip_aligned_saved_independent_audit.v1","status":audit_status,
        "physical_fracture_status":"UNKNOWN_NOT_ASSESSED","interpretation":"Synthetic E=1 MPa, nu=0.3, KIII=100 benchmark only; algebraic/array consistency is not fracture validation or material calibration.",
        "pins":{"result_json":str(rp),"result_sha256_pre":result_hash_pre,"result_sha256_post":result_hash_post,
                "npz":str(npz_path),"npz_sha256_pre":npz_hash_pre,"npz_sha256_post":npz_hash_post,
                "study_driver_sha256_pre":driver_hash_pre,"study_driver_sha256_post":driver_hash_post,
                "study_case_sha256_pre":case_hash_pre,"study_case_sha256_post":case_hash_post,"reviewer_sha256":sha(Path(__file__))},
        "npz_schema_preflight":schema,"review_scope_flags":{"basis_product_gradient":basis_ok,"fixed_load_clouds_and_tractions":fixed_load_ok,"primary_finite_G_Pi_match":ledger_ok},"array_consistency_flags":algebra_flags,"basis_product_gradient_checks":basis_checks,
        "independent_raw_K_C_f_KKT_energy":algebra,"fixed_load_cloud_and_traction_checks":fixed_load,
        "a0_offgrid_old_tip_exact_comparison":offgrid,"finite_delta_a_and_theta_ledger":ledger,
        "plots":["a0_offgrid_stress_scatter.png","finite_G_theta_ledger.png"],"elapsed_seconds":time.monotonic()-t0}
    audit_path=outdir/"audit.json"; payload=json.dumps(doc,indent=2,allow_nan=False)+"\n"
    if len(payload.encode("utf-8"))+sum((outdir/n).stat().st_size for n in doc["plots"]) > 2*1024*1024: raise ValueError("reviewer output exceeds 2 MiB")
    audit_path.write_text(payload,encoding="utf-8")
    print(json.dumps({"status":audit_status,"audit_json":str(outdir/"audit.json"),"elapsed_seconds":doc["elapsed_seconds"]}))
    return 2 if audit_status.startswith("BLOCKED") else 0

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--result",type=Path,default=DEFAULT_RESULT)
    p.add_argument("--out",type=Path,default=ROOT)
    return audit(p.parse_args())
if __name__ == "__main__":
    raise SystemExit(main())
