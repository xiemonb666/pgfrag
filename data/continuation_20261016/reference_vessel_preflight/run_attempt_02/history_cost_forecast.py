from __future__ import annotations
import io
import json
import math
from pathlib import Path
import numpy as np

root = Path('/mnt/c/Users/xiemo/Desktop/ShardForge/pgfrag/experiments/continuation_20261016/reference_vessel_preflight/run_attempt_02')
case_dir = root / 'inputs'
run_dir = root / 'native_preflight/run01'
case = json.loads((case_dir / 'case.json').read_text(encoding='utf-8'))
result = json.loads((run_dir / 'result.json').read_text(encoding='utf-8'))
bound = json.loads((run_dir / 'initial_bound.json').read_text(encoding='utf-8'))
loading = json.loads((run_dir / 'loading.json').read_text(encoding='utf-8'))
assert result['status'] == 'ZERO_STEP_PREFLIGHT_PASS' and result['solver_steps'] == 0
with np.load(run_dir / 'initial_state.npz', allow_pickle=False) as state:
    rest = np.asarray(state['rest_nodes_m'])
    tets = np.asarray(state['tets'])
    cohesive = np.asarray(state['cohesive'])
    mass = np.asarray(state['mass_kg'])
    qp_area = np.asarray(state['qp_area_m2'])
    unit_load = np.asarray(state['unit_external_load'])
    ancestry = np.asarray(state['native_to_source'])
    assert len(rest) == result['counts'][0] and len(tets) == result['counts'][1] and len(cohesive) == result['counts'][2]
    assert mass.shape == (len(rest), 3) and qp_area.shape == (len(cohesive) * 3,)
    assert np.allclose(mass[:, 0], mass[:, 1], rtol=1e-8, atol=0.0)
    node_count, tet_count, facet_count = len(rest), len(tets), len(cohesive)
    qp_count = len(qp_area)
    patch_node_count = int(np.count_nonzero(unit_load[:, 0]))
    patch_negative_count = int(np.count_nonzero(unit_load[:, 0] < 0))
    patch_positive_count = int(np.count_nonzero(unit_load[:, 0] > 0))
    mass_kg = float(mass[:, 0].sum())
    unit_force_balance = np.sum(unit_load, axis=0)
    unit_torque_balance = np.cross(rest, unit_load).sum(axis=0)
    min_positive_unit_weight = float(np.min(np.abs(unit_load[unit_load[:, 0] != 0, 0])))
    qp_area_sum = float(qp_area.sum())
    native_source_nodes = int(len(np.unique(ancestry)))

steps = int(bound['steps'])
frames = steps + 1
peak = float(case['loading']['peak_force_N'])
assert peak == 0.0, 'This remains a zero-force cost preflight.'
assert patch_node_count == sum(len(p['native_node_ids']) for p in loading['patches'])
assert all(float(w) > 0 for p in loading['patches'] for w in p['weights'])

ledger_names = ['time_s', 'load_factor_N', 'interface_work_J', 'external_work_J', 'algorithmic_kinetic_J',
                'momentum_x', 'momentum_y', 'momentum_z', 'angular_momentum_x', 'angular_momentum_y', 'angular_momentum_z',
                'external_torque_x_Nm', 'external_torque_y_Nm', 'external_torque_z_Nm',
                'ma_minus_f_internal_external_max_N', 'position_recurrence_max_m', 'velocity_recurrence_max_mps',
                'external_work_minus_total_alg_J', 'max_damage', 'min_damage',
                'internal_work_cumulative_J', 'kinetic_alg_increment_identity_residual_J',
                'internal_work_plus_material_energy_J', 'native_minus_nodal_kinetic_J',
                'internal_current_torque_x_Nm', 'internal_current_torque_y_Nm', 'internal_current_torque_z_Nm',
                'max_contact_traction_Pa']
layout_base = {
    'ledger_columns': ledger_names,
    'energy_columns': ['kinetic', 'potential', 'dissipated', 'reversible', 'cohesive contact'],
    'node_state_columns': ['displacement_m', 'velocity_mps', 'acceleration_mps2', 'internal_force_N'],
    'contact_energy_check_columns': ['independent_quadrature_contact_energy_J', 'native_contact_energy_J', 'independent_minus_native_J'],
    'interface_work_check_columns': ['full_interface_work_J', 'cohesive_work_J', 'contact_work_J', 'cross_work_J'],
    'qp_order': 'native element-major, three quadrature points per cohesive_3d_6; no averaging',
    'all_qp_states_retained': True,
}
strides = {}
for stride in (50, 1):
    saved_steps = np.unique(np.r_[np.arange(0, steps + 1, stride), steps]).astype(np.int64)
    specs = {
        'opening_m': ((frames, qp_count, 3), np.float64),
        'traction_Pa': ((frames, qp_count, 3), np.float64),
        'damage': ((frames, qp_count), np.float64),
        'delta_max_m': ((frames, qp_count), np.float64),
        'contact_traction_Pa': ((frames, qp_count, 3), np.float64),
        'contact_opening_m': ((frames, qp_count, 3), np.float64),
        'contact_energy_check': ((frames, 3), np.float64),
        'interface_work_check': ((frames, 4), np.float64),
        'energy_J': ((frames, 5), np.float64),
        'ledger': ((frames, 28), np.float64),
        'patch_displacement_x_m': ((frames, patch_node_count), np.float64),
        'node_state': ((len(saved_steps), node_count, 4, 3), np.float64),
    }
    raw_arrays = sum(int(np.prod(shape, dtype=np.int64)) * np.dtype(dtype).itemsize for shape, dtype in specs.values())
    wrapper_forecast = sum(int(np.prod(shape, dtype=np.int64)) * np.dtype(dtype).itemsize + 256 for shape, dtype in specs.values())
    node_steps_buffer = io.BytesIO()
    np.save(node_steps_buffer, saved_steps)
    node_saved_steps_bytes = len(node_steps_buffer.getvalue())
    layout = dict(layout_base)
    layout['node_stride'] = stride
    layout_json_bytes = len((json.dumps(layout, indent=2, allow_nan=False) + '\n').encode('utf-8'))
    full_history_folder_estimate = wrapper_forecast + node_saved_steps_bytes + layout_json_bytes
    static_preflight_bytes = int(json.loads((run_dir / 'receipt.json').read_text(encoding='utf-8'))['persistent_bytes_before_receipt'])
    estimated_run_directory_bytes = static_preflight_bytes + full_history_folder_estimate
    margin_bytes = 2_000_000
    cap = int(case['limits']['output_cap_bytes'])
    recommended_bytes_20pct = math.ceil((estimated_run_directory_bytes + margin_bytes) * 1.2)
    recommended_allocation_gib_rounded = math.ceil(recommended_bytes_20pct / (1024 ** 3))
    strides[str(stride)] = {
        'node_stride': stride,
        'steps': steps,
        'time_frames': frames,
        'node_state_frames': int(len(saved_steps)),
        'qp_count': qp_count,
        'patch_node_count': patch_node_count,
        'raw_numeric_array_bytes': raw_arrays,
        'wrapper_raw_history_bytes_forecast': wrapper_forecast,
        'node_saved_steps_npy_bytes': node_saved_steps_bytes,
        'layout_json_bytes': layout_json_bytes,
        'full_history_folder_estimate_bytes': full_history_folder_estimate,
        'actual_zero_step_preflight_output_bytes_excluding_receipt': static_preflight_bytes,
        'estimated_run_directory_bytes_before_safety_margin': estimated_run_directory_bytes,
        'wrapper_2mb_safety_margin_bytes': margin_bytes,
        'configured_output_cap_bytes': cap,
        'wrapper_history_cap_check_passes': wrapper_forecast + margin_bytes <= cap,
        'bytes_over_configured_history_cap': max(0, wrapper_forecast + margin_bytes - cap),
        'recommended_minimum_bytes_with_20pct_margin': recommended_bytes_20pct,
        'recommended_output_allocation_gib_rounded_up': recommended_allocation_gib_rounded,
    }

report = {
    'schema': 'pgfrag.native_cohesive_vessel_history_cost_forecast.v1',
    'status': 'COST_FORECAST_FROM_ZERO_STEP_NATIVE_PREFLIGHT',
    'classification': 'PLACEHOLDER_COST_ESTIMATE_UNCALIBRATED; NO SPECIMEN PHYSICAL DATA',
    'solver_evolution_performed': False,
    'native_preflight_status': result['status'],
    'native_counts_nodes_tets_cohesive_facets': result['counts'],
    'native_quadrature_points': qp_count,
    'native_source_nodes_after_seam_insertion': node_count,
    'unique_source_ancestor_nodes': native_source_nodes,
    'mass_kg': mass_kg,
    'interface_area_m2': float(bound['interface_area_m2']),
    'qp_area_sum_m2': qp_area_sum,
    'balanced_patch': {
        'status': 'FEASIBLE_POSITIVE_RELATIVE_ENTROPY_PROJECTION',
        'force_peak_N': peak,
        'unit_load_is_used_for_patch_feasibility_only': True,
        'patch_node_count_total': patch_node_count,
        'negative_x_patch_nodes': patch_negative_count,
        'positive_x_patch_nodes': patch_positive_count,
        'min_positive_unit_weight': min_positive_unit_weight,
        'target_reference_com_yz_m': loading['target_reference_com_m'][1:],
        'negative_reference_force_line_yz_m': loading['patches'][0]['reference_force_line_yz_m'],
        'positive_reference_force_line_yz_m': loading['patches'][1]['reference_force_line_yz_m'],
        'reference_unit_force_sum_N_per_N': unit_force_balance.tolist(),
        'reference_unit_torque_sum_Nm_per_N': unit_torque_balance.tolist(),
        'current_configuration_policy': loading['current_configuration_torque_policy'],
    },
    'stability_and_time': {
        'lambda_upper_bound_per_s2': bound['lambda_upper_bound_per_s2'],
        'critical_time_step_s': bound['critical_time_step_s'],
        'bulk_only_native_dt_s': bound['bulk_only_native_dt_s'],
        'selected_dt_s': bound['dt_s'],
        'selected_steps_for_50us': steps,
        'limiting_bound': '0.4 * bulk-only native dt',
        'scope': bound['scope'],
        'all_dofs_free': bound['all_dofs_free'],
    },
    'per_frame_history_cost': {
        'native_qp_count': qp_count,
        'qp_field_doubles_per_point': 14,
        'qp_field_bytes_per_time_frame': qp_count * 14 * 8,
        'energy_work_ledger_doubles_per_time_frame': 3 + 4 + 5 + 28,
        'energy_work_ledger_bytes_per_time_frame': (3 + 4 + 5 + 28) * 8,
        'patch_displacement_bytes_per_time_frame': patch_node_count * 8,
        'node_state_bytes_per_saved_frame': node_count * 4 * 3 * 8,
        'wrapper_array_header_allowance_bytes': 12 * 256,
    },
    'node_stride_forecasts': strides,
    'resource_recommendation': {
        'storage_only': True,
        'node_stride_50_allocation_gib': strides['50']['recommended_output_allocation_gib_rounded_up'],
        'node_stride_1_allocation_gib': strides['1']['recommended_output_allocation_gib_rounded_up'],
        'configured_64mib_cap_supports_either_forecast': False,
        'recommendation': 'For the next solver-capacity review, reserve the rounded allocation shown for the chosen nodal stride plus the existing preflight artifacts; these are storage forecasts only and provide no wall-time prediction.',
    },
    'source_geometry_limitations': {
        'candidate': case['source_geometry']['candidate'],
        'exterior_source_face_ancestry': case['source_geometry']['source_face_ancestry']['status'],
        'unknown_exterior_triangles': case['source_geometry']['source_face_ancestry']['exterior_triangle_counts']['unknown'],
        'qualified_fragments_or_surface_paint': False,
    },
}
out = Path('/mnt/c/Users/xiemo/Desktop/ShardForge/pgfrag/data/continuation_20261016/reference_vessel_preflight/run_attempt_02/history_cost_forecast.json')
out.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n', encoding='utf-8')
print(json.dumps({
    'status': report['status'],
    'native_counts': result['counts'],
    'mass_kg': mass_kg,
    'steps': steps,
    'dt_s': bound['dt_s'],
    'critical_dt_s': bound['critical_time_step_s'],
    'bulk_dt_s': bound['bulk_only_native_dt_s'],
    'balanced_patch_nodes': patch_node_count,
    'positive_negative_counts': [patch_negative_count, patch_positive_count],
    'unit_force_balance': unit_force_balance.tolist(),
    'unit_torque_balance': unit_torque_balance.tolist(),
    'qp_count': qp_count,
    'qp_frame_bytes': qp_count * 14 * 8,
    'stride50_forecast': strides['50'],
    'stride1_forecast': strides['1'],
    'report_path': str(out),
}, separators=(',', ':')))