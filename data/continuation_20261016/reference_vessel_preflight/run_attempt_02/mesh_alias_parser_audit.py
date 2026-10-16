from __future__ import annotations
import json
import sys
from pathlib import Path
import numpy as np

root = Path('/mnt/c/Users/xiemo/Desktop/ShardForge/pgfrag/experiments/continuation_20261016/reference_vessel_preflight/run_attempt_02')
inputs = root / 'inputs'
archive = root / 'source_archive'
sys.path.insert(0, str(inputs))
from cohesive_mesh import read_gmsh22

source = read_gmsh22(archive / 'mesh.msh.original')
executed = read_gmsh22(inputs / 'mesh.msh')
assert np.array_equal(source['nodes'], executed['nodes'])
assert source['node_tags'] == executed['node_tags']

def signature(elements):
    return [(int(e['element_id']), int(e['element_type']), int(e['physical_tag']), tuple(e['node_tags']), tuple(e['nodes'])) for e in elements]
assert signature(source['elements']) == signature(executed['elements'])
source_ids = [(int(d), int(t)) for (d, t) in source['physical_names']]
executed_ids = [(int(d), int(t)) for (d, t) in executed['physical_names']]
assert source_ids == executed_ids
assert source['physical_names'] == {(2, 2): 'cut_interface', (2, 3): 'source_exterior', (3, 1): 'partitioned_body'}
assert executed['physical_names'] == {(2, 2): 'interface', (2, 3): 'exterior', (3, 1): 'partitioned_body'}
source_cut = [e for e in source['triangles'] if e['physical_name'] == 'cut_interface']
executed_cut = [e for e in executed['triangles'] if e['physical_name'] == 'interface']
source_ext = [e for e in source['triangles'] if e['physical_name'] == 'source_exterior']
executed_ext = [e for e in executed['triangles'] if e['physical_name'] == 'exterior']
assert len(source_cut) == len(executed_cut) == 541
assert len(source_ext) == len(executed_ext) == 35784
report = {
    'status': 'PASS_RAW_AND_PARSER_EQUIVALENCE',
    'parser': 'frozen cohesive_mesh.read_gmsh22',
    'nodes_count': len(source['nodes']),
    'nodes_coordinates_and_order_equal': True,
    'node_tags_order_equal': True,
    'tetrahedra_count': len(source['tetrahedra']),
    'tetrahedra_connectivity_element_ids_physical_tags_and_order_equal': signature(source['tetrahedra']) == signature(executed['tetrahedra']),
    'triangles_count': len(source['triangles']),
    'triangles_connectivity_element_ids_physical_tags_and_order_equal': signature(source['triangles']) == signature(executed['triangles']),
    'interface_triangles_original': len(source_cut),
    'interface_triangles_derived': len(executed_cut),
    'exterior_triangles_original': len(source_ext),
    'exterior_triangles_derived': len(executed_ext),
    'physical_names_original': [{'dimension': d, 'physical_tag': t, 'name': n} for (d, t), n in source['physical_names'].items()],
    'physical_names_derived': [{'dimension': d, 'physical_tag': t, 'name': n} for (d, t), n in executed['physical_names'].items()],
    'physical_dimension_tag_order_equal': source_ids == executed_ids,
}
out = Path('/mnt/c/Users/xiemo/Desktop/ShardForge/pgfrag/data/continuation_20261016/reference_vessel_preflight/run_attempt_02/parser_equivalence.json')
out.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
print(json.dumps(report, separators=(',', ':')))