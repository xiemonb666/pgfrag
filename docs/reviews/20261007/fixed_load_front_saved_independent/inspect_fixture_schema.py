"""Standard-library-only schema inspection of pinned cached JSON inputs."""
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=Path(__file__).resolve().parents[4]
paths={
    "fixed":ROOT/"pgfrag/tmp/review_20261006/fixed_load_front_reanalysis_20261006/run01/result.json",
    "baseline":ROOT/"pgfrag/tmp/review_20261006/full3d_static_solve_bridge_20261006/run02/static_solve_bridge_result.json",
    "operator":ROOT/"tmp/review_20261006/full3d_cracked_elastic_operator_20261006/run01/result.json",
}
data={k:json.loads(v.read_text(encoding="utf-8")) for k,v in paths.items()}
def at(root,path):
    value=root
    for key in path.split("."): value=value[int(key)] if isinstance(value,list) else value[key]
    return value
def shape(value):
    if isinstance(value,list): return [len(value)]+(shape(value[0]) if value else [])
    return []
def desc(value):
    d={"type":type(value).__name__}
    if isinstance(value,(list,dict)): d["shape_or_size"]=shape(value) if isinstance(value,list) else len(value)
    elif isinstance(value,(str,int,float,bool)) or value is None: d["value"]=value
    return d
want={
"fixed":["schema","status","source_results.operator_sha256","source_results.baseline_static_solve_sha256","front.delta_a_m","fixed_data.query_rows_reused","fixed_data.additional_queries","fixed_data.traction_values_little_endian_f64_sha256","fixed_data.quadrature_and_boundary_hash_sha256","assembly.K_N_per_m","assembly.f_N","spectrum.scaled_eigenvalues","gauge.C_dimensionless","gauge.g_fixed_from_baseline_m","solve.coefficients_m","solve.lagrange_multipliers_N","solve.strain_energy_J","solve.external_work_J","solve.g_dot_lambda_J","solve.energy_identity_defect_2U_minus_W_plus_g_dot_lambda_J","solve.post_solve_residual_gate_passed","potential_energy.baseline_Pi_J","potential_energy.extended_Pi_J","potential_energy.delta_Pi_J","potential_energy.one_sided_crack_extension_area_m2","potential_energy.finite_change_release_diagnostic_J_per_m2","potential_energy.Gc_compared","potential_energy.growth_committed","rigid_load_check.integrated_boundary_force_N","rigid_load_check.integrated_moment_about_origin_Nm","rigid_load_check.integrated_moment_about_rigid_center_Nm"],
"baseline":["schema","status","source_result.sha256","gauge.C_dimensionless","gauge.g_target_m","solve.coefficients_m","solve.lagrange_multipliers_N","solve.strain_energy_J","solve.external_work_J","solve.g_dot_lambda_J","solve.raw_K_sha256_float64_little_endian_C_order"],
"operator":["schema","status","static_solve_replay.exact_field.statuses.0","static_solve_replay.exact_field.statuses.383","static_solve_replay.quadrature.positions_xyz_m","static_solve_replay.quadrature.volume_weights_m3","static_solve_replay.quadrature.boundary_area_weights_m2","static_solve_replay.quadrature.boundary_normals","static_solve_replay.basis.values","static_solve_replay.basis.gradients_per_m","static_solve_replay.basis_row_ids","static_solve_replay.basis_row_ids.0.source_id","static_solve_replay.basis_row_ids.0.basis_kind","static_solve_replay.basis_row_ids.53.source_id","static_solve_replay.basis_row_ids.53.basis_kind","static_solve_replay.registry.ordered_keys","static_solve_replay.registry.ordered_keys.0.slot","static_solve_replay.registry.ordered_keys.0.source_id","static_solve_replay.registry.ordered_keys.0.component","static_solve_replay.registry.ordered_keys.0.basis_kind","static_solve_replay.registry.ordered_keys.161.slot","static_solve_replay.registry.ordered_keys.161.source_id","static_solve_replay.registry.ordered_keys.161.component","static_solve_replay.registry.ordered_keys.161.basis_kind","static_solve_replay.registry.source_ids","static_solve_replay.boundary.traction_Pa","static_solve_replay.operator.K_N_per_m","static_solve_replay.operator.boundary_load_N","static_solve_replay.provider_evidence.volume.condition_numbers","static_solve_replay.provider_evidence.boundary.condition_numbers","static_solve_replay.sources.rigid_mode_center_xyz_m","fixed_case.geometry.radius_m","fixed_case.geometry.length_m","fixed_case.material_and_manufactured_field.young_modulus_Pa","fixed_case.material_and_manufactured_field.poisson_ratio","fixed_case.material_and_manufactured_field.mu_Pa","fixed_case.material_and_manufactured_field.lambda_Pa","fixed_case.material_and_manufactured_field.K_III_Pa_sqrt_m"]}
out={"schema":"fixed_load_front_stdlib_fixture_schema.v1","numeric_libraries_imported":False,"fields":{name:{p:desc(at(data[name],p)) for p in plist} for name,plist in want.items()}}
(HERE/"fixture_schema.json").write_text(json.dumps(out,indent=2,allow_nan=False)+"\n",encoding="utf-8")
print(json.dumps({"status":"SCHEMA_INSPECTED_STDLIB_ONLY","material":{k:desc(data["operator"]["fixed_case"]["material_and_manufactured_field"][k]) for k in ("young_modulus_Pa","poisson_ratio","mu_Pa","lambda_Pa","K_III_Pa_sqrt_m")}},allow_nan=False))
