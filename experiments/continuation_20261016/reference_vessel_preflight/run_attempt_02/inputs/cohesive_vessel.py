"""Bounded native Akantu experiments on conformal, preinserted vessel seams.

All bulk DOFs are free. Opposite smooth nodal patch loads supply the energy;
fracture is decided by the native cohesive law, never by the geometry tags.
Run a frozen copy with sibling ``cohesive_mesh`` and ``cohesive_stability``.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time
import traceback

import numpy as np

try:
    from .cohesive_mesh import (read_gmsh22, build_native_lineage, sha256,
                                input_interface_contract, quadrature_signed_work)
    from .cohesive_stability import mass_normalized_stiffness_bound
except ImportError:
    from cohesive_mesh import (read_gmsh22, build_native_lineage, sha256,
                               input_interface_contract, quadrature_signed_work)
    from cohesive_stability import mass_normalized_stiffness_bound

ENERGY_KEYS = ("kinetic", "potential", "dissipated", "reversible", "cohesive contact")


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def positive_centroid_weights(points_yz, base_weights, target_yz, *, tolerance=1e-11):
    """Relative-entropy projection to a declared force line, with positive weights.

    w_i proportional to b_i exp(lambda dot (p_i-target)). Newton's Hessian is
    the weighted covariance. This avoids negative forces from a linear fixup.
    Failure to reach an interior target rejects the load; it does not move it.
    """
    p = np.asarray(points_yz, dtype=float)
    b = np.asarray(base_weights, dtype=float)
    target = np.asarray(target_yz, dtype=float)
    if p.ndim != 2 or p.shape[1] != 2 or b.shape != (len(p),) or target.shape != (2,):
        raise ValueError("invalid load weight dimensions")
    if not all(np.isfinite(x).all() for x in (p, b, target)) or np.any(b <= 0):
        raise ValueError("patch coordinates/base weights must be finite and weights positive")
    scale = max(float(np.ptp(p, axis=0).max()), 1e-12)
    q = (p-target)/scale
    lam = np.zeros(2)
    logbase = np.log(b/b.sum())

    def state(l):
        logits = logbase + q @ l
        logits -= logits.max()
        w = np.exp(logits); w /= w.sum()
        g = w @ q
        h = (q-g).T @ (w[:, None]*(q-g))
        return w, g, h

    for _ in range(100):
        w, g, h = state(lam)
        if np.linalg.norm(g, ord=np.inf)*scale <= tolerance:
            if np.any(w <= 0) or not np.isclose(w.sum(), 1., atol=1e-13):
                raise ValueError("load projection lost positive normalization")
            return w
        try:
            direction = np.linalg.solve(h, g)
        except np.linalg.LinAlgError as exc:
            raise ValueError("force-line target is infeasible or patch is degenerate") from exc
        error = float(g @ g)
        fraction = 1.
        for _ in range(35):
            candidate = lam-fraction*direction
            _, newg, _ = state(candidate)
            if float(newg @ newg) < error:
                lam = candidate
                break
            fraction *= .5
        else:
            raise ValueError("positive force-line projection did not converge")
    raise ValueError("force-line target not reached")


def make_balanced_patch_load(source, native_nodes, native_to_source, mass, threshold_m):
    """Area-based outward x loads, with a common reference COM force line."""
    x = np.asarray(native_nodes)
    source_x = source["nodes"]
    ancestry_raw = np.asarray(native_to_source)
    mass = np.asarray(mass, dtype=float)
    if (x.ndim != 2 or x.shape[1] != 3 or not np.isfinite(x).all()
        or ancestry_raw.shape != (len(x),) or ancestry_raw.dtype.kind not in "iu"
        or np.any(ancestry_raw < 0) or np.any(ancestry_raw >= len(source_x))
        or mass.shape != (len(x),) or not np.isfinite(mass).all() or np.any(mass <= 0)
        or not np.isfinite(threshold_m) or threshold_m <= 0):
        raise ValueError("invalid native/source ancestry, positive point mass, or patch threshold")
    target = np.sum(mass[:, None]*x, axis=0)/mass.sum()
    source_weights = np.zeros(len(source_x))
    exterior = [e for e in source["triangles"] if e["physical_name"] == "exterior"]
    if not exterior:
        raise ValueError("declared exterior ancestry is required for loading")
    for e in exterior:
        ids = e["nodes"]; p = source_x[ids]
        area = .5*np.linalg.norm(np.cross(p[1]-p[0], p[2]-p[0]))
        source_weights[ids] += area/3
    ancestry = np.asarray(native_to_source, dtype=np.int64)
    load = np.zeros_like(x)
    patches = []
    for sign in (-1, 1):
        ids = np.flatnonzero((sign*x[:, 0] > threshold_m) & (source_weights[ancestry] > 0))
        if len(ids) < 4:
            raise ValueError("finite exterior load patch has too few nodes")
        if len(np.unique(ancestry[ids])) != len(ids):
            raise ValueError("loading patch contains duplicated seam nodes")
        base = source_weights[ancestry[ids]]
        weights = positive_centroid_weights(x[ids, 1:], base, target[1:])
        load[ids, 0] = sign*weights
        patches.append({"sign": sign, "native_node_ids": ids.tolist(), "weights": weights.tolist(),
                        "base_lumped_exterior_area_m2": float(base.sum()),
                        "equivalent_uniform_force_area_m2": float(1./np.sum(weights**2/base)),
                        "maximum_force_weight": float(weights.max()),
                        "reference_force_line_yz_m": (weights @ x[ids, 1:]).tolist()})
    if np.linalg.norm(load.sum(axis=0)) > 1e-12 or np.linalg.norm(np.cross(x, load).sum(axis=0)) > 1e-11:
        raise ValueError("reference force/moment balance failed")
    return load, {"patches": patches, "target_reference_com_m": target.tolist(),
                  "reference_unit_force_sum": load.sum(axis=0).tolist(),
                  "reference_unit_torque_sum_m": np.cross(x, load).sum(axis=0).tolist(),
                  "current_configuration_torque_policy": "READ_ACTUAL_TORQUE_EACH_STEP; NO_FOLLOWER_REBALANCING"}


def load_factor(time_s, duration_s, peak_N):
    if not 0 < time_s < duration_s:
        return 0.
    return float(peak_N*np.sin(np.pi*time_s/duration_s)**2)


def deformation_metrics(rest, current, tets):
    """Small strain, exact affine stretch and polar angle for each linear T4."""
    x, y = np.asarray(rest)[tets], np.asarray(current)[tets]
    j0 = np.swapaxes(x[:, 1:]-x[:, :1], 1, 2)
    j1 = np.swapaxes(y[:, 1:]-y[:, :1], 1, 2)
    f = j1 @ np.linalg.inv(j0)
    strain = .5*(f+np.swapaxes(f, 1, 2))-np.eye(3)
    small = float(np.abs(np.linalg.eigvalsh(strain)).max())
    u, singular, vt = np.linalg.svd(f)
    rot = u @ vt
    cosine = np.clip((np.trace(rot, axis1=1, axis2=2)-1)/2, -1, 1)
    return {"max_abs_principal_small_strain": small,
            "min_det_F": float(np.linalg.det(f).min()),
            "max_polar_rotation_rad": float(np.arccos(cosine).max()),
            "min_principal_stretch": float(singular.min()),
            "max_principal_stretch": float(singular.max())}


def contact_energy_check(contact_opening, contact_traction, qp_area_m2, native_energy_J):
    """Compare the native contact energy with its quadrature point formula."""
    opening = np.asarray(contact_opening, dtype=np.float64)
    traction = np.asarray(contact_traction, dtype=np.float64)
    area = np.asarray(qp_area_m2, dtype=np.float64)
    native = float(native_energy_J)
    if opening.ndim != 2 or opening.shape[1] != 3 or traction.shape != opening.shape:
        raise ValueError("contact opening and traction must share shape (n_qp, 3)")
    if area.shape != (len(opening),) or np.any(area <= 0.0):
        raise ValueError("contact quadrature areas must be positive and match the number of points")
    if not (np.all(np.isfinite(opening)) and np.all(np.isfinite(traction))
            and np.all(np.isfinite(area)) and np.isfinite(native)):
        raise ValueError("contact energy inputs must be finite")
    independent = 0.5 * float(np.sum(area * np.einsum("qi,qi->q", traction, opening)))
    return np.asarray([independent, native, independent - native], dtype=np.float64)


def validate_contact_penalty(penalty_Pa_m, normal_stiffness_Pa_m):
    """Validate the fixed-normal contact tangent against the Kn envelope."""
    try:
        penalty = float(penalty_Pa_m)
        normal_stiffness = float(normal_stiffness_Pa_m)
    except (TypeError, ValueError, OverflowError) as error:
        raise ValueError("contact penalty and Kn must be finite scalars") from error
    if (not np.isfinite(penalty) or not np.isfinite(normal_stiffness)
            or normal_stiffness <= 0.0 or not 0.0 <= penalty <= normal_stiffness):
        raise ValueError("contact penalty_Pa_m must be finite and satisfy 0 <= penalty <= Kn")
    return penalty


def interface_work_increment(
    opening_previous,
    traction_previous,
    contact_opening_previous,
    contact_traction_previous,
    opening_current,
    traction_current,
    contact_opening_current,
    contact_traction_current,
    qp_area_m2,
):
    """Return full, cohesive, contact, and cross interface work increments.

    Each increment uses the trapezoid rule over native quadrature points. The
    cross term is retained because changing normals or branches can make it
    nonzero.
    """
    opening_0 = np.asarray(opening_previous, dtype=np.float64)
    traction_0 = np.asarray(traction_previous, dtype=np.float64)
    contact_opening_0 = np.asarray(contact_opening_previous, dtype=np.float64)
    contact_traction_0 = np.asarray(contact_traction_previous, dtype=np.float64)
    opening_1 = np.asarray(opening_current, dtype=np.float64)
    traction_1 = np.asarray(traction_current, dtype=np.float64)
    contact_opening_1 = np.asarray(contact_opening_current, dtype=np.float64)
    contact_traction_1 = np.asarray(contact_traction_current, dtype=np.float64)
    if not (opening_0.shape == traction_0.shape == contact_opening_0.shape == contact_traction_0.shape
            == opening_1.shape == traction_1.shape == contact_opening_1.shape == contact_traction_1.shape):
        raise ValueError("cohesive/contact opening and traction frames must share shape (n_qp, 3)")
    area = np.asarray(qp_area_m2, dtype=np.float64)
    if area.shape != (len(opening_0),):
        raise ValueError("quadrature areas must match the number of points")
    cohesive = float(np.sum(quadrature_signed_work(opening_0, traction_0, opening_1, traction_1, area)))
    contact = float(np.sum(quadrature_signed_work(
        contact_opening_0, contact_traction_0, contact_opening_1, contact_traction_1, area
    )))
    full = float(np.sum(quadrature_signed_work(
        opening_0 + contact_opening_0,
        traction_0 + contact_traction_0,
        opening_1 + contact_opening_1,
        traction_1 + contact_traction_1,
        area,
    )))
    return np.asarray([full, cohesive, contact, full - cohesive - contact], dtype=np.float64)


def fields(model, aka, facets):
    result = {}
    for key, name, comp in (("opening", "opening", 3), ("traction", "tractions", 3),
                            ("damage", "damage", 1), ("delta_max", "delta max", 1),
                            ("contact_traction", "contact_tractions", 3),
                            ("contact_opening", "contact_opening", 3)):
        a = np.asarray(model.flattenInternalReal(name, aka._ek_cohesive, aka._not_ghost)
                       (aka._cohesive_3d_6, aka._not_ghost)).copy()
        if a.ndim != 2 or a.shape[1] != comp or len(a) % facets:
            raise ValueError(f"unexpected native {name} shape {a.shape}")
        result[key] = a
    if len(result["opening"]) != facets*3:
        raise ValueError("this Tri3 contract requires measured three quadrature rows per facet")
    return result


def file_inventory(root):
    return {str(p.relative_to(root)).replace("\\", "/"):
            {"bytes": p.stat().st_size, "sha256": sha256(p)}
            for p in sorted(root.rglob("*")) if p.is_file() and p.name != "receipt.json"}


def run(case_path, out, preflight_only=False):
    import akantu as aka
    started = time.monotonic()
    root = case_path.parent
    case = json.loads(case_path.read_text(encoding="utf-8"))
    if case["schema"] != "pgfrag.native_cohesive_vessel_case.v1":
        raise ValueError("case schema mismatch")
    for name, expected in case["input_hashes"].items():
        if sha256(root/name) != expected:
            raise ValueError(f"frozen input hash mismatch: {name}")
    source = read_gmsh22(root/"mesh.msh")
    contract = input_interface_contract(source)
    aka.parseInput(str(root/"material.dat"))
    mesh = aka.Mesh(3); mesh.read(str(root/"mesh.msh"))
    before_nodes = np.asarray(mesh.getNodes()).copy()
    before_tets = np.asarray(mesh.getConnectivity(aka._tetrahedron_4)).copy()
    if not np.array_equal(before_nodes, source["nodes"]) or not np.array_equal(before_tets, [e["nodes"] for e in source["tetrahedra"]]):
        raise ValueError("native initial row/local ordinal differs from frozen MSH")
    model = aka.SolidMechanicsModelCohesive(mesh)
    model.initFull(_analysis_method=aka._explicit_lumped_mass, _is_extrinsic=False)
    model.assembleMassLumped()
    rest = np.asarray(mesh.getNodes()).copy()
    tets = np.asarray(mesh.getConnectivity(aka._tetrahedron_4)).copy()
    coh = np.asarray(mesh.getConnectivity(aka._cohesive_3d_6)).copy()
    mass3 = np.asarray(model.getMass()).copy()
    blocked = np.asarray(model.getBlockedDOFs()); blocked[:] = False
    if not (len(coh) == contract["input_face_count"] and len(tets) == len(before_tets)):
        raise ValueError("native insertion counts differ from declared interface")
    lineage = build_native_lineage(source, rest, tets, coh, native_initial_bulk_tets=before_tets,
                    expected_area_m2=contract["rest_area_sum_m2"], plane_axis=0, plane_coordinate_m=0.)
    ancestry = np.asarray([lineage["native_node_to_input_node_zero_based"][str(i)] for i in range(len(rest))])
    pairs = []
    for face in coh:
        side_b = {ancestry[n]: n for n in face[3:]}
        pairs.append([face[:3], [side_b[ancestry[n]] for n in face[:3]]])
    pairs = np.asarray(pairs, dtype=np.int64)
    material = case["bulk_material"]
    kn = case["interface"]["sigma_c_Pa"]/case["interface"]["delta0_m"]
    penalty = validate_contact_penalty(case["interface"].get("penalty_Pa_m", 0.0), kn)
    # beta=kappa=1 gives disjoint fixed-normal cohesive normal/shear branches.
    # The compression-contact tangent is bounded by the same Kn I envelope
    # when penalty <= Kn. This only extends the fixed-normal envelope; it is
    # not a nonlinear or full-trajectory stability guarantee.
    if case["interface"]["beta"] != 1 or case["interface"]["kappa"] != 1:
        raise ValueError("this reviewed vessel bound requires beta=kappa=1")
    bound = mass_normalized_stiffness_bound(rest, tets, pairs, material["rho_kg_m3"],
                  material["E_pa"], material["nu"], kn, kn, blocked, mass3)
    bulk_dt = float(model.getStableTimeStep())
    maximum_dt = min(case["time"]["initial_bound_fraction"]*bound.critical_time_step_s,
                     .4*bulk_dt)
    total_time = case["time"]["total_s"]
    steps = math.ceil(total_time/maximum_dt)*case["time"].get("refinement_factor", 1)
    if steps > case["limits"]["max_steps"]:
        raise ValueError(f"required steps {steps} exceed declared cap")
    dt = total_time/steps
    unit_load, load_report = make_balanced_patch_load(source, rest, ancestry, mass3[:, 0], case["loading"]["patch_x_threshold_m"])
    initial = fields(model, aka, len(coh))
    if any(np.max(np.abs(initial[k])) > 1e-15 for k in (
        "opening", "traction", "damage", "contact_opening", "contact_traction"
    )):
        raise ValueError("nonzero native initial field")
    areas = np.repeat(np.asarray([p["positive_rest_area_m2"] for p in lineage["parent_facets"]])/3, 3)
    bound_report = {"lambda_upper_bound_per_s2": bound.lambda_upper_bound_per_s2,
                    "critical_time_step_s": bound.critical_time_step_s,
                    "bulk_only_native_dt_s": bulk_dt, "dt_s": dt, "steps": steps,
                    "mass_kg": bound.bulk_total_lumped_mass_kg,
                    "interface_area_m2": bound.cohesive_total_area_m2,
                    "scope": (bound.stability_scope + "; FIXED_NORMAL_COMPRESSION_BRANCH_PENALTY_LE_KN_COVERED_BY_KN_ENVELOPE; "
                              "NO_NONLINEAR_OR_FULL_TRAJECTORY_GUARANTEE"),
                    "contact_penalty_Pa_per_m": penalty,
                    "all_dofs_free": bool(not blocked.any())}
    np.savez_compressed(out/"initial_state.npz", rest_nodes_m=rest, tets=tets, cohesive=coh,
                        mass_kg=mass3, blocked=blocked, native_to_source=ancestry,
                        source_node_tags=np.asarray(source["node_tags"]),
                        source_tet_ids=np.asarray([e["element_id"] for e in source["tetrahedra"]]),
                        source_exterior_faces=np.asarray([e["nodes"] for e in source["triangles"] if e["physical_name"] == "exterior"]),
                        cohesive_pairs=pairs, qp_area_m2=areas, unit_external_load=unit_load)
    write_json(out/"lineage.json", lineage)
    write_json(out/"initial_bound.json", bound_report)
    write_json(out/"loading.json", load_report)
    if preflight_only:
        return {"status": "ZERO_STEP_PREFLIGHT_PASS", "solver_steps": 0, "counts": [len(rest), len(tets), len(coh)], **bound_report}

    stride = case["output"]["node_stride"]
    saved_steps = np.unique(np.r_[np.arange(0, steps+1, stride), steps]).astype(np.int64)
    specs = {"opening_m": ((steps+1, len(areas), 3), np.float64),
             "traction_Pa": ((steps+1, len(areas), 3), np.float64),
             "damage": ((steps+1, len(areas)), np.float64),
             "delta_max_m": ((steps+1, len(areas)), np.float64),
             "contact_traction_Pa": ((steps+1, len(areas), 3), np.float64),
             "contact_opening_m": ((steps+1, len(areas), 3), np.float64),
             "contact_energy_check": ((steps+1, 3), np.float64),
             "interface_work_check": ((steps+1, 4), np.float64),
             "energy_J": ((steps+1, 5), np.float64),
             "ledger": ((steps+1, 28), np.float64),
             "patch_displacement_x_m": ((steps+1, int(np.count_nonzero(unit_load[:, 0]))), np.float64),
             "node_state": ((len(saved_steps), len(rest), 4, 3), np.float64)}
    forecast = sum(np.prod(shape)*np.dtype(dtype).itemsize+256 for shape, dtype in specs.values())
    if forecast + 2_000_000 > case["limits"]["output_cap_bytes"]:
        raise ValueError(f"preallocated history forecast {forecast} exceeds cap")
    hist = out/"history"; hist.mkdir()
    arrays = {name: np.lib.format.open_memmap(hist/f"{name}.npy", mode="w+", dtype=dtype, shape=shape)
              for name, (shape, dtype) in specs.items()}
    np.save(hist/"node_saved_steps.npy", saved_steps)
    ledger_names = ["time_s", "load_factor_N", "interface_work_J", "external_work_J", "algorithmic_kinetic_J",
                    "momentum_x", "momentum_y", "momentum_z", "angular_momentum_x", "angular_momentum_y", "angular_momentum_z",
                    "external_torque_x_Nm", "external_torque_y_Nm", "external_torque_z_Nm",
                    "ma_minus_f_internal_external_max_N", "position_recurrence_max_m", "velocity_recurrence_max_mps",
                    "external_work_minus_total_alg_J", "max_damage", "min_damage",
                    "internal_work_cumulative_J", "kinetic_alg_increment_identity_residual_J",
                    "internal_work_plus_material_energy_J", "native_minus_nodal_kinetic_J",
                    "internal_current_torque_x_Nm", "internal_current_torque_y_Nm", "internal_current_torque_z_Nm",
                    "max_contact_traction_Pa"]
    write_json(hist/"layout.json", {"ledger_columns": ledger_names, "energy_columns": ENERGY_KEYS,
                "node_state_columns": ["displacement_m", "velocity_mps", "acceleration_mps2", "internal_force_N"],
                "contact_energy_check_columns": ["independent_quadrature_contact_energy_J", "native_contact_energy_J", "independent_minus_native_J"],
                "interface_work_check_columns": ["full_interface_work_J", "cohesive_work_J", "contact_work_J", "cross_work_J"],
                "qp_order": "native element-major, three quadrature points per cohesive_3d_6; no averaging",
                "node_stride": stride, "all_qp_states_retained": True})
    model.setTimeStep(dt)
    work = 0.; contact_work = 0.; full_interface_work = 0.
    extwork = 0.; internalwork = 0.; previous = initial
    prev_u = np.asarray(model.getDisplacement()).copy()
    prev_v = np.asarray(model.getVelocity()).copy()
    prev_a = np.asarray(model.getAcceleration()).copy()
    prev_f = np.zeros_like(rest)
    prev_i = np.asarray(model.getInternalForce()).copy()
    prev_kalg = 0.
    first_failure = np.full(len(areas), -1, dtype=np.int64)
    selected_patch = np.flatnonzero(unit_load[:, 0])
    snapshot = 0
    maxima = {"max_abs_principal_small_strain": 0., "max_polar_rotation_rad": 0., "min_det_F": 1.}
    last_metrics = None
    completed = 0
    try:
        for step in range(steps+1):
            factor = load_factor(step*dt, case["loading"]["pulse_duration_s"], case["loading"]["peak_force_N"])
            applied = factor*unit_load
            if step:
                np.asarray(model.getExternalForce())[:] = applied
                model.solveStep()
                completed = step
            current = fields(model, aka, len(coh)) if step else initial
            u = np.asarray(model.getDisplacement()).copy()
            v = np.asarray(model.getVelocity()).copy()
            a = np.asarray(model.getAcceleration()).copy()
            internal = np.asarray(model.getInternalForce()).copy()
            external = np.asarray(model.getExternalForce()).copy()
            energies = np.asarray([float(model.getEnergy(key)) for key in ENERGY_KEYS])
            if step == 0 and np.max(np.abs(energies)) > 1e-15:
                raise ValueError("this stress-free zero-energy experiment requires zero initial native energies")
            if not all(np.isfinite(z).all() for z in (*current.values(), u, v, a, internal, external, energies)):
                raise ValueError(f"nonfinite native state at {step}")
            if blocked.any() or not np.array_equal(external, applied):
                raise ValueError("unexpected boundary constraint or external-force alteration")
            if step:
                increments = interface_work_increment(
                    previous["opening"], previous["traction"],
                    previous["contact_opening"], previous["contact_traction"],
                    current["opening"], current["traction"],
                    current["contact_opening"], current["contact_traction"], areas,
                )
                full_interface_work += float(increments[0])
                work += float(increments[1])
                contact_work += float(increments[2])
                extwork += float(np.sum(.5*(prev_f+external)*(u-prev_u)))
                internalwork += float(np.sum(.5*(prev_i+internal)*(u-prev_u)))
            damage = current["damage"][:, 0]
            first_failure[(first_failure < 0) & (damage >= 1-case["gates"]["damage_tolerance"])] = step
            momentum = (mass3[:, :1]*v).sum(axis=0)
            angular = np.cross(rest+u, mass3[:, :1]*v).sum(axis=0)
            torque = np.cross(rest+u, external).sum(axis=0)
            kalg = float(energies[0] - dt**2*np.sum(mass3*a*a)/8)
            force_error = float(np.abs(mass3*a-internal-external).max())
            position_error = float(np.abs(u-prev_u-dt*prev_v-.5*dt**2*prev_a).max()) if step else 0.
            velocity_error = float(np.abs(v-prev_v-.5*dt*(prev_a+a)).max()) if step else 0.
            arrays["opening_m"][step] = current["opening"]
            arrays["traction_Pa"][step] = current["traction"]
            arrays["damage"][step] = damage
            arrays["delta_max_m"][step] = current["delta_max"][:, 0]
            arrays["contact_traction_Pa"][step] = current["contact_traction"]
            arrays["contact_opening_m"][step] = current["contact_opening"]
            arrays["contact_energy_check"][step] = contact_energy_check(
                current["contact_opening"], current["contact_traction"], areas, energies[4]
            )
            arrays["interface_work_check"][step] = [
                full_interface_work, work, contact_work, full_interface_work - work - contact_work
            ]
            arrays["energy_J"][step] = energies
            arrays["patch_displacement_x_m"][step] = u[selected_patch, 0]
            arrays["ledger"][step] = np.r_[step*dt, factor, work, extwork, kalg, momentum, angular, torque,
                force_error, position_error, velocity_error, extwork-kalg-energies[1:].sum(), damage.max(), damage.min(),
                internalwork, kalg-prev_kalg-float(np.sum(.5*(prev_f+external+prev_i+internal)*(u-prev_u))) if step else 0.,
                internalwork+energies[1:].sum(), float(energies[0]-.5*np.sum(mass3*v*v)),
                np.cross(rest+u, internal).sum(axis=0), float(np.abs(current["contact_traction"]).max())]
            if snapshot < len(saved_steps) and step == saved_steps[snapshot]:
                arrays["node_state"][snapshot] = np.stack((u, v, a, internal), axis=1)
                snapshot += 1
                last_metrics = deformation_metrics(rest, rest+u, tets)
                for key in maxima:
                    maxima[key] = min(maxima[key], last_metrics[key]) if key == "min_det_F" else max(maxima[key], last_metrics[key])
                if (last_metrics["max_abs_principal_small_strain"] > case["gates"]["max_abs_principal_small_strain"]
                    or last_metrics["max_polar_rotation_rad"] > case["gates"]["max_polar_rotation_rad"]
                    or last_metrics["min_det_F"] < case["gates"]["min_det_F"]):
                    raise ValueError(f"small-strain model applicability exceeded at {step}: {last_metrics}")
                if not np.array_equal(np.asarray(mesh.getNodes()), rest) or not np.array_equal(np.asarray(mesh.getConnectivity(aka._tetrahedron_4)), tets) or not np.array_equal(np.asarray(mesh.getConnectivity(aka._cohesive_3d_6)), coh):
                    raise ValueError("native topology/rest coordinates changed")
            if step % case["output"]["flush_stride"] == 0 or step == steps:
                for value in arrays.values():
                    value.flush()
                write_json(out/"checkpoint.json", {"completed_steps": step, "qp_frames_committed": step+1,
                             "node_frames_committed": snapshot, "elapsed_s": time.monotonic()-started})
                print(json.dumps({"step": step, "of": steps, "damage_max": float(damage.max()),
                                  "elapsed_s": time.monotonic()-started}), flush=True)
            if time.monotonic()-started > case["limits"]["wall_time_s"]:
                raise TimeoutError("declared wall-time cap reached")
            previous = current; prev_u = u; prev_v = v; prev_a = a; prev_f = external; prev_i = internal; prev_kalg = kalg
    finally:
        for value in arrays.values():
            value.flush()
        np.save(out/"qp_first_full_damage_step.npy", first_failure)
        write_json(out/"completion.json", {"steps_completed": completed, "qp_frames_written": completed+1,
                                          "node_frames_written": snapshot, "deformation_maxima_sampled": maxima,
                                          "last_metrics": last_metrics})
    return {"status": "NATIVE_RUN_COMPLETE_DIAGNOSTIC", "solver_steps": completed,
            "counts": [len(rest), len(tets), len(coh)], "dt_s": dt, "elapsed_s": time.monotonic()-started,
            "all_qp_full_damage": bool((first_failure >= 0).all()),
            "raw_history_bytes_forecast": int(forecast), "deformation_maxima_sampled": maxima,
            "qualification": "PREINSERTED_WEAK_SEAM_UNCALIBRATED; NO_PHYSICAL_OR_TRAINING_QUALIFICATION"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--preflight-only", action="store_true")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    os.chdir(args.out)
    started = time.time()
    snapshot_hashes = {p.name: sha256(p) for p in args.case.parent.iterdir() if p.is_file()}
    write_json(args.out/"started.json", {"argv": sys.argv, "started_unix_s": started,
                "source_input_inventory": snapshot_hashes, "python": sys.version})
    result = {}
    try:
        result = run(args.case.resolve(), args.out.resolve(), args.preflight_only)
        result["source_inputs_unchanged"] = all(sha256(args.case.parent/name) == h for name, h in snapshot_hashes.items())
        if not result["source_inputs_unchanged"]:
            raise ValueError("executed frozen sources/inputs changed during run")
    except Exception as exc:
        result = {"status": "FAILED_PRESERVED_NO_RETRY", "error": repr(exc), "traceback": traceback.format_exc()}
    write_json(args.out/"result.json", result)
    receipt = {"argv": sys.argv, "started_unix_s": started, "ended_unix_s": time.time(),
               "status": result["status"], "files": file_inventory(args.out)}
    receipt["persistent_bytes_before_receipt"] = sum(x["bytes"] for x in receipt["files"].values())
    write_json(args.out/"receipt.json", receipt)
    print(json.dumps(result, allow_nan=False), flush=True)
    return 2 if result["status"] == "FAILED_PRESERVED_NO_RETRY" else 0


if __name__ == "__main__":
    raise SystemExit(main())
