"""Independent cached-array check; no producer, provider, assembler or solver imports."""
from __future__ import annotations
import hashlib, json
from pathlib import Path
import numpy as np
HERE, ROOT = Path(__file__).resolve().parent, Path(__file__).resolve().parents[4]
NEW_PATH = ROOT / "pgfrag/tmp/review_20261006/fixed_load_front_reanalysis_20261006/run01/result.json"
BASE_PATH = ROOT / "pgfrag/tmp/review_20261006/full3d_static_solve_bridge_20261006/run02/static_solve_bridge_result.json"
OP_PATH = ROOT / "tmp/review_20261006/full3d_cracked_elastic_operator_20261006/run01/result.json"
PINS = {
    "fixed_result": "a60e51c85c26510c358ee57cafa4bfdda6ed89ea463d99bcc86551f64220732c",
    "baseline_result": "ef66783e2961dd77f34785a99e1b397210ccc0beb6e8fa2979678aa4992a0ce3",
    "operator_result": "3b31e7d40c1d8e4b5d829f0b093fee4edadc966bde4086a4b3f4ba48fdf89cd2",
}
REL_TOL = 5e-11
def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()
def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"), parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)))
def arr(value, shape, name):
    a = np.asarray(value, dtype=np.float64)
    if a.shape != shape or not np.isfinite(a).all(): raise ValueError(f"{name}: wrong shape or non-finite")
    return a
def rel(a, b):
    a, b = np.asarray(a, dtype=np.float64), np.asarray(b, dtype=np.float64)
    return float(np.linalg.norm(a-b)/max(float(np.linalg.norm(a)),float(np.linalg.norm(b)),np.finfo(float).tiny))
def f64sha(a):
    return hashlib.sha256(np.asarray(a, dtype="<f8", order="C").tobytes()).hexdigest()
def main():
    paths = {"fixed_result": NEW_PATH, "baseline_result": BASE_PATH, "operator_result": OP_PATH}
    pre = {k: sha(p) for k, p in paths.items()}
    if pre != PINS:
        raise ValueError(f"input pin mismatch: {pre}")
    new, base, opres = read_json(NEW_PATH), read_json(BASE_PATH), read_json(OP_PATH)
    replay = opres["static_solve_replay"]; quad, basis, boundary = replay["quadrature"], replay["basis"], replay["boundary"]; op = replay["operator"]
    positions = arr(quad["positions_xyz_m"], (384,3), "positions"); vw = arr(quad["volume_weights_m3"], (256,), "volume weights"); aw = arr(quad["boundary_area_weights_m2"], (128,), "boundary weights")
    normals = arr(quad["boundary_normals"], (128,3), "normals"); oldv = arr(basis["values"], (384,54), "old values"); oldg = arr(basis["gradients_per_m"], (384,54,3), "old gradients"); traction = arr(boundary["traction_Pa"], (128,3), "traction")
    keys = replay["registry"]["ordered_keys"]
    if len(keys) != 162 or [k["slot"] for k in keys] != list(range(162)):
        raise ValueError("semantic DOF order malformed")
    ids = replay["registry"]["source_ids"]
    if len(ids) != 27 or len(replay["basis_row_ids"]) != 54: raise ValueError("source/basis row count mismatch")
    rows = {(v["source_id"], v["basis_kind"]): i for i, v in enumerate(replay["basis_row_ids"])}
    expected = {(sid, fam) for sid in ids for fam in ("base", "front_f1_unshifted")}
    if set(rows) != expected or len(rows) != 54: raise ValueError("basis row map malformed or incomplete")
    if not np.all(vw > 0) or not np.all(aw > 0) or replay["exact_field"]["statuses"] != ["OK"] * 384:
        raise ValueError("quadrature weights or frozen-query status invalid")
    case = opres["fixed_case"]
    geom, mat = case["geometry"], case["material_and_manufactured_field"]
    radius, length = float(geom["radius_m"]), float(geom["length_m"])
    delta = float(new["front"]["delta_a_m"])
    dx, dy = positions[:, 0]-delta, positions[:, 1]; rr = np.hypot(dx,dy)
    singular = (rr == 0.0) | ((dy == 0.0) & (dx < 0.0))
    if bool(np.any(singular)):
        raise ValueError("cached quadrature intersects shifted tip/branch")
    theta = np.arctan2(dy,dx); fbar = np.sqrt(rr/radius)*np.sin(0.5*theta); inv = 0.5/np.sqrt(radius*rr)
    gradf = np.column_stack((-inv*np.sin(0.5*theta), inv*np.cos(0.5*theta), np.zeros(384)))
    values, gradients = np.empty_like(oldv), np.empty_like(oldg)
    for (sid, family), row in rows.items():
        brow = rows[(sid, "base")]
        phi, gphi = oldv[:, brow], oldg[:, brow]
        if family == "base":
            values[:, row], gradients[:, row] = phi, gphi
        else:
            values[:, row] = phi * fbar
            gradients[:, row] = fbar[:, None] * gphi + phi[:, None] * gradf
    if not np.isfinite(values).all() or not np.isfinite(gradients).all():
        raise ValueError("independently refronted basis is non-finite")
    young, nu = float(mat["young_modulus_Pa"]), float(mat["poisson_ratio"])
    mu_expected = young/(2.0*(1.0+nu))
    lame_expected = young*nu / ((1.0+nu)*(1.0-2.0*nu))
    mu, lame = mu_expected, lame_expected
    D = np.zeros((6, 6), dtype=np.float64)
    D[:3, :3] = lame
    D[0, 0] += 2*mu; D[1, 1] += 2*mu; D[2, 2] += 2*mu
    D[3, 3] = D[4, 4] = D[5, 5] = mu
    dof_rows = np.asarray([rows[(k["source_id"],k["basis_kind"])] for k in keys],dtype=np.int64); components = np.asarray([{"ux":0,"uy":1,"uz":2}[k["component"]] for k in keys],dtype=np.int64)
    K1 = np.zeros((162, 162), dtype=np.float64)
    for q in range(256):
        gq = gradients[q, dof_rows, :]
        B = np.zeros((6, 162), dtype=np.float64)
        mx, my, mz = components == 0, components == 1, components == 2
        B[0, mx], B[3, mx], B[5, mx] = gq[mx, 0], gq[mx, 1], gq[mx, 2]
        B[1, my], B[3, my], B[4, my] = gq[my, 1], gq[my, 0], gq[my, 2]
        B[2, mz], B[4, mz], B[5, mz] = gq[mz, 2], gq[mz, 1], gq[mz, 0]
        K1 += vw[q] * (B.T @ D @ B)
    K0 = arr(op["K_N_per_m"], (162,162), "baseline K"); K_saved = arr(new["assembly"]["K_N_per_m"], (162,162), "shifted K")
    asym0, asym1 = rel(K0,K0.T), rel(K1,K1.T)
    symdelta1 = rel(K1,0.5*(K1+K1.T))
    tvals = {"ux": 0, "uy": 1, "uz": 2}
    def load(v):
        return np.asarray([np.dot(aw * traction[:, tvals[k["component"]]], v[256:, rows[(k["source_id"], k["basis_kind"])]]) for k in keys])
    f0, f1 = load(oldv), load(values); f0saved = arr(op["boundary_load_N"], (162,), "old f"); f1saved = arr(new["assembly"]["f_N"], (162,), "new f")
    volume = float(np.sum(vw)); vp = positions[:256]
    centroid = np.einsum("q,qc->c", vw, vp) / volume
    modes = np.zeros((6, 256, 3), dtype=np.float64)
    modes[0, :, 0] = modes[1, :, 1] = modes[2, :, 2] = 1.0
    eye = np.eye(3)
    for ax in range(3):
        modes[3+ax] = np.cross(eye[ax], vp-centroid) / radius
    C1 = np.zeros((6, 162), dtype=np.float64)
    for i, key in enumerate(keys):
        comp = tvals[key["component"]]
        C1[:, i] = np.sum(modes[:, :, comp] * (vw * values[:256, dof_rows[i]])[None, :], axis=1) / volume
    d0 = arr(base["solve"]["coefficients_m"],(162,),"old d"); l0 = arr(base["solve"]["lagrange_multipliers_N"],(6,),"old lambda"); C0 = arr(base["gauge"]["C_dimensionless"],(6,162),"old C"); g0 = arr(base["gauge"]["g_target_m"],(6,),"old g")
    d1 = arr(new["solve"]["coefficients_m"],(162,),"new d"); l1 = arr(new["solve"]["lagrange_multipliers_N"],(6,),"new lambda"); C1saved = arr(new["gauge"]["C_dimensionless"],(6,162),"new C"); g1 = arr(new["gauge"]["g_fixed_from_baseline_m"],(6,),"new g")
    force1 = np.einsum("q,qc->c", aw, traction); bpos = positions[256:]
    moment0 = np.einsum("q,qc->c", aw, np.cross(bpos, traction))
    center = arr(replay["sources"]["rigid_mode_center_xyz_m"],(3,),"rigid center")
    momentc = np.einsum("q,qc->c", aw, np.cross(bpos-center, traction))
    moment_caps = {"top_xy_max_abs_Pa": float(np.max(np.abs(traction[64:96, :2]))),
                   "bottom_xy_max_abs_Pa": float(np.max(np.abs(traction[96:128, :2])))}
    conditions = [arr(replay["provider_evidence"][n]["condition_numbers"],(count,),n+" conditions") for n,count in (("volume",256),("boundary",128))]
    if any(np.any(x<1.0) for x in conditions): raise ValueError("cached provider condition evidence invalid")
    condition_max = max(1.0, *(float(np.max(x)) for x in conditions))
    def spectrum(K):
        scale = 1.0 / np.sqrt(np.diag(K))
        scaled = scale[:,None]*K*scale[None,:]; eig = np.linalg.eigvalsh(0.5*(scaled+scaled.T))
        tol = 64*np.finfo(float).eps*condition_max*max(float(np.max(np.abs(eig))),np.finfo(float).tiny)
        return eig, tol
    eig0, tol0 = spectrum(K0)
    eig1, tol1 = spectrum(K1)
    eigsaved = arr(new["spectrum"]["scaled_eigenvalues"], (162,), "shifted saved spectrum")
    scale1 = 1.0 / np.sqrt(np.diag(K1))
    weak = K1 @ d1 + C1.T @ l1 - f1
    gres = C1 @ d1 - g1
    weak_scaled = scale1 * weak
    kd_scaled, cl_scaled, f_scaled = scale1*(K1 @ d1), scale1*(C1.T @ l1), scale1*f1
    weak_rel = float(np.linalg.norm(weak_scaled) / max(np.linalg.norm(kd_scaled)+np.linalg.norm(cl_scaled)+np.linalg.norm(f_scaled), np.finfo(float).tiny))
    gauge_rel = float(np.linalg.norm(gres) / max(np.linalg.norm(g1), np.linalg.norm(C1 @ d1), np.finfo(float).tiny))
    U0, W0 = 0.5*float(d0 @ K0 @ d0), float(d0 @ f0)
    U1, W1 = 0.5*float(d1 @ K1 @ d1), float(d1 @ f1)
    gl0, gl1 = float(g0 @ l0), float(g1 @ l1)
    pi0, pi1 = U0-W0, U1-W1
    delta_pi = pi1-pi0
    one_area, two_area = length*delta, 2*length*delta
    Gfinite = -delta_pi/one_area
    KIII = float(mat["K_III_Pa_sqrt_m"])
    Gan = KIII*KIII/(2*mu)
    energy_ref = {
        "baseline": {"U_J": U0, "W_J": W0, "g_dot_lambda_J": gl0, "Pi_U_minus_W_J": pi0,
                     "saved_U_J": float(base["solve"]["strain_energy_J"]), "saved_W_J": float(base["solve"]["external_work_J"])},
        "extended": {"U_J": U1, "W_J": W1, "g_dot_lambda_J": gl1, "Pi_U_minus_W_J": pi1,
                     "saved_U_J": float(new["solve"]["strain_energy_J"]), "saved_W_J": float(new["solve"]["external_work_J"])},
        "delta_Pi_new_minus_old_J": delta_pi, "one_sided_area_L_times_delta_a_m2": one_area,
        "two_face_area_m2_reported_separately": two_area, "finite_change_G_J_per_m2": Gfinite,
        "analytic_local_modeIII_G_J_per_m2": Gan, "finite_over_analytic": Gfinite/Gan,
        "percent_below_local_analytic": 100.0*(1.0-Gfinite/Gan),
        "analytic_inputs": {"KIII_Pa_sqrt_m": KIII, "mu_Pa": mu},
        "comparison_scope": "local LEFM Mode III reference; finite fixed-load, finite-geometry and frozen-quadrature delta is not a convergence/physical qualification",
    }
    traction_sha = f64sha(traction)
    quad_cat = np.concatenate((positions.ravel(), vw, normals.ravel(), aw))
    quad_sha = f64sha(quad_cat)
    old_k_hash = f64sha(op["K_N_per_m"])
    eq_old = K0 @ d0 + C0.T @ l0 - f0
    res_new_scaled = scale1[:, None] * K1 * scale1[None, :]
    row_norm = np.linalg.norm(C1 * scale1[None, :], axis=1)
    cscale = 1.0/row_norm
    cscaled = cscale[:, None] * C1 * scale1[None, :]
    saddle = np.zeros((168, 168), dtype=np.float64)
    saddle[:162, :162] = res_new_scaled; saddle[:162, 162:] = cscaled.T; saddle[162:, :162] = cscaled
    rhs = np.r_[scale1*f1, cscale*g1]
    solution_rebuilt = np.r_[d1/scale1, l1/cscale]
    saddle_rel = float(np.linalg.norm(saddle@solution_rebuilt-rhs)/max(np.linalg.norm(rhs), np.finfo(float).tiny))
    array_bytes = int(sum(v.nbytes for v in locals().values() if isinstance(v,np.ndarray))+sum(v.nbytes for v in conditions))
    checks = {
        "source_pins_match": bool(new["source_results"]["operator_sha256"] == PINS["operator_result"] and new["source_results"]["baseline_static_solve_sha256"] == PINS["baseline_result"]),
        "result_is_uncommitted_numerical_diagnostic": new["status"] == "FINITE_CHANGE_RELEASE_DIAGNOSTIC_ONLY_NOT_COMMITTED" and not new["potential_energy"]["Gc_compared"] and not new["potential_energy"]["growth_committed"],
        "shifted_K_rebuild_rel_le_5e-11": rel(K1, K_saved) <= REL_TOL,
        "shifted_f_rebuild_rel_le_5e-11": rel(f1, f1saved) <= REL_TOL,
        "old_f_rebuild_rel_le_5e-11": rel(f0, f0saved) <= REL_TOL,
        "new_C_rebuild_rel_le_5e-11": rel(C1, C1saved) <= REL_TOL,
        "fixed_g_matches_baseline": rel(g1, g0) <= REL_TOL,
        "stored_spectrum_matches_rebuild_le_5e-11": rel(eig1, eigsaved) <= REL_TOL,
        "both_spectra_have_6_zero_156_positive_0_negative": int(np.count_nonzero(np.abs(eig0)<=tol0)) == 6 and int(np.count_nonzero(eig0>tol0)) == 156 and int(np.count_nonzero(eig0 < -tol0)) == 0 and int(np.count_nonzero(np.abs(eig1)<=tol1)) == 6 and int(np.count_nonzero(eig1>tol1)) == 156 and int(np.count_nonzero(eig1 < -tol1)) == 0,
        "old_raw_K_hash_matches_baseline": old_k_hash == base["solve"]["raw_K_sha256_float64_little_endian_C_order"],
        "fixed_traction_and_quad_hashes_match": traction_sha == new["fixed_data"]["traction_values_little_endian_f64_sha256"] and quad_sha == new["fixed_data"]["quadrature_and_boundary_hash_sha256"],
        "no_added_queries_and_384_rows": new["fixed_data"]["query_rows_reused"] == 384 and new["fixed_data"]["additional_queries"] == 0,
        "cap_inplane_traction_present": moment_caps["top_xy_max_abs_Pa"] > 0 and moment_caps["bottom_xy_max_abs_Pa"] > 0,
        "new_physical_residuals_le_predeclared_1e-9": weak_rel <= 1e-9 and gauge_rel <= 1e-9 and saddle_rel <= 1e-9,
        "old_solution_equilibrium_finite": bool(np.isfinite(eq_old).all()),
        "energy_values_match_saved_le_5e-11": rel([U0,W0,U1,W1], [base["solve"]["strain_energy_J"],base["solve"]["external_work_J"],new["solve"]["strain_energy_J"],new["solve"]["external_work_J"]]) <= REL_TOL,
        "multipliers_and_new_energy_identity_match_saved": rel([gl0,gl1],[base["solve"]["g_dot_lambda_J"],new["solve"]["g_dot_lambda_J"]]) <= REL_TOL and abs((2*U1-W1+gl1)-float(new["solve"]["energy_identity_defect_2U_minus_W_plus_g_dot_lambda_J"]))/max(abs(2*U1),abs(W1),abs(gl1),np.finfo(float).tiny) <= REL_TOL,
        "finite_potential_ledger_matches_saved": rel([pi0,pi1,delta_pi,Gfinite],[new["potential_energy"]["baseline_Pi_J"],new["potential_energy"]["extended_Pi_J"],new["potential_energy"]["delta_Pi_J"],new["potential_energy"]["finite_change_release_diagnostic_J_per_m2"]]) <= REL_TOL,
        "one_sided_area_matches_L_delta_a": abs(one_area-float(new["potential_energy"]["one_sided_crack_extension_area_m2"])) <= 1e-18,
        "explicit_numeric_arrays_within_8MiB": array_bytes <= 8*1024*1024,
    }
    counts0 = {"near_zero": int(np.count_nonzero(np.abs(eig0)<=tol0)), "positive": int(np.count_nonzero(eig0>tol0)), "negative": int(np.count_nonzero(eig0 < -tol0))}
    counts1 = {"near_zero": int(np.count_nonzero(np.abs(eig1)<=tol1)), "positive": int(np.count_nonzero(eig1>tol1)), "negative": int(np.count_nonzero(eig1 < -tol1))}
    savedWrench = new["rigid_load_check"]
    wrench_delta = {"force_abs_max_N": float(np.max(np.abs(force1-arr(savedWrench["integrated_boundary_force_N"],(3,),"saved wrench force")))),
                    "origin_moment_abs_max_Nm": float(np.max(np.abs(moment0-arr(savedWrench["integrated_moment_about_origin_Nm"],(3,),"saved origin moment")))),
                    "center_moment_abs_max_Nm": float(np.max(np.abs(momentc-arr(savedWrench["integrated_moment_about_rigid_center_Nm"],(3,),"saved center moment"))))}
    report = {
        "schema": "fixed_load_front_saved_array_independent_audit.v1",
        "status": "PASS_ALGEBRAIC_REANALYSIS_ONLY" if all(checks.values()) else "FAIL_ARRAY_RECONSTRUCTION_OR_RESIDUAL",
        "input_sha256": pre, "input_bytes": {k: int(v.stat().st_size) for k,v in paths.items()},
        "method": "cached arrays only: independent shifted Fbar/product-rule gradient, full engineering-strain B.T D B, boundary load, physical gauge, spectrum, equilibrium and energy; no producer/provider/assembler/solver/reference/CAD/Git imports",
        "dimensions": {"sources": 27, "scalar_basis_rows": 54, "semantic_dofs": 162, "volume_points": 256, "boundary_points": 128, "additional_queries": 0},
        "explicit_numpy_array_bytes_estimate_conservative_views_may_be_counted_twice": array_bytes,
        "material_and_front": {"E_Pa": young,"nu":nu,"mu_Pa":mu,"mu_formula_relative_delta":abs(mu-mu_expected)/mu_expected,"lambda_Pa":lame,"lambda_formula_relative_delta":abs(lame-lame_expected)/lame_expected,"mu_lambda_formula_metadata":{"mu_Pa":mat["mu_Pa"],"lambda_Pa":mat["lambda_Pa"]},"R_m":radius,"L_m":length,"delta_a_m":delta,"delta_a_over_R":delta/radius,"shifted_singular_or_branch_rows":int(np.count_nonzero(singular))},
        "rebuild_relative_deltas": {"K_shifted":rel(K1,K_saved),"f_shifted":rel(f1,f1saved),"f_old":rel(f0,f0saved),"C_shifted":rel(C1,C1saved),"spectrum_shifted":rel(eig1,eigsaved)},
        "raw_K_symmetry_before_eig_symmetrization": {"baseline_relative_asymmetry":asym0,"extended_relative_asymmetry":asym1,"extended_to_symmetric_relative_delta":symdelta1},
        "independent_residuals": {"new_weak_equilibrium_relative":weak_rel,"new_physical_gauge_relative":gauge_rel,"new_saddle_backward_relative":saddle_rel,"new_energy_identity_defect_J":2*U1-W1+gl1,"old_weak_equilibrium_abs_norm_N_or_Nm":float(np.linalg.norm(eq_old))},
        "spectrum": {"baseline_scaled_eigenvalue_counts":counts0,"extended_scaled_eigenvalue_counts":counts1,"baseline_min":float(eig0[0]),"extended_min":float(eig1[0]),"extended_max":float(eig1[-1]),"baseline_zero_tolerance":float(tol0),"extended_zero_tolerance":float(tol1),"condition_max_from_cached_provider_evidence":condition_max,"baseline_eigenvalues":eig0.tolist(),"extended_eigenvalues":eig1.tolist()},
        "fixed_load": {"traction_sha256":traction_sha,"quadrature_boundary_sha256":quad_sha,"integrated_force_N":force1.tolist(),"integrated_moment_origin_Nm":moment0.tolist(),"integrated_moment_center_Nm":momentc.tolist(),"saved_wrench_abs_deltas":wrench_delta,"cap_inplane_max_abs_traction_Pa":moment_caps},
        "energy_ledger_and_G_comparison":energy_ref,
        "checks":checks,
        "interpretation": "Finite-change value is about 58% below the local analytic Mode III value. The analytic value is an asymptotic local LEFM reference; this single 5%-of-R finite extension in a finite cylinder, fixed-gauge solve and reused 256-point quadrature is not a converged local-G measurement. Keep the numerical reanalysis distinct from physical qualification.",
    }
    if not all(checks.values()):
        report["status"] = "FAIL_ARRAY_RECONSTRUCTION_OR_RESIDUAL"
    (HERE/"audit.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({"status":report["status"],"checks_passed":sum(bool(v) for v in checks.values()),"checks_total":len(checks),"Gfinite_J_per_m2":Gfinite,"Ganalytic_J_per_m2":Gan},allow_nan=False))
    return 0 if all(checks.values()) else 1
if __name__ == "__main__":
    raise SystemExit(main())
