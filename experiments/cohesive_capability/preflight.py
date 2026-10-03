"""Runtime import, mesh-reader, and cohesive-model initialization preflight."""

from __future__ import annotations

from pathlib import Path
import json
import importlib.metadata as metadata
import platform
import sys

import akantu as aka
import numpy as np

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "input"
REVIEW = ROOT.parents[2] / "relicweave-synlab" / "tmp" / "review_20261004" / "implementation" / "cohesive_capability"
REVIEW.mkdir(parents=True, exist_ok=True)

aka.parseInput(str(INPUT / "material.dat"))
mesh = aka.Mesh(3)
mesh.read(str(INPUT / "block_2x2x4mm.msh"))
bulk_count = mesh.getNbElement(type=aka._tetrahedron_4)
node_count = mesh.getNbNodes()
model = aka.SolidMechanicsModelCohesive(mesh)
model.initFull(
    _analysis_method=aka._explicit_lumped_mass,
    _is_extrinsic=True,
)
model.assembleMassLumped()
mass = np.asarray(model.getMass()).copy()
initial_energies = {
    name: float(model.getEnergy(name))
    for name in ("kinetic", "potential", "dissipated", "reversible", "contact")
}
inserter = model.getElementInserter()
inserter.addPhysicalSurface("interface")
model.updateAutomaticInsertion()
insertion_facets = inserter.getInsertionFacets()
facet_types = insertion_facets.elementTypes()
check_facets = inserter.getCheckFacets()


def describe_array(value):
    return {
        "type": type(value).__name__,
        "shape": getattr(value, "shape", None),
        "size": getattr(value, "size", None),
        "dtype": str(getattr(value, "dtype", None)),
        "true_count": int(np.count_nonzero(value)),
    }


facet_arrays = []
for element_type in facet_types:
    try:
        facet_arrays.append(
            {"element_type": str(element_type), **describe_array(insertion_facets(element_type))}
        )
    except Exception as error:  # include binding diagnostics in the preflight record
        facet_arrays.append({"element_type": str(element_type), "error": repr(error)})

check_facet_types = check_facets.elementTypes()
check_facet_arrays = []
for element_type in check_facet_types:
    try:
        check_facet_arrays.append(
            {"element_type": str(element_type), **describe_array(check_facets(element_type))}
        )
    except Exception as error:
        check_facet_arrays.append({"element_type": str(element_type), "error": repr(error)})

result = {
    "akantu_version": metadata.version("akantu"),
    "numpy_version": metadata.version("numpy"),
    "scipy_version": metadata.version("scipy"),
    "module_path": aka.__file__,
    "python_executable": sys.executable,
    "python_version": platform.python_version(),
    "platform": platform.platform(),
    "spatial_dimension": 3,
    "bulk_tetrahedra": bulk_count,
    "nodes_before_insertion": node_count,
    "cohesive_elements_before_insertion": mesh.getNbElement(type=aka._cohesive_3d_6),
    "mass_array_shape": list(mass.shape),
    "mass_array_dtype": str(mass.dtype),
    "mass_raw_sum_over_all_dof_components_kg": float(mass.sum()),
    "mass_component_sums_kg": [float(mass[:, component].sum()) for component in range(mass.shape[1])],
    "expected_physical_body_mass_rhoV_kg": 2500.0 * 16.0e-9,
    "initial_energies_J": initial_energies,
    "stable_time_step_s": model.getStableTimeStep(),
    "cohesive_model_api": [
        name
        for name in dir(model)
        if any(key in name.lower() for key in ("cohes", "energy", "mesh", "material"))
    ],
    "element_inserter_type": type(inserter).__name__,
    "element_inserter_api": dir(inserter),
    "eligible_insertion_facet_types": [str(element_type) for element_type in facet_types],
    "eligible_insertion_facet_arrays": facet_arrays,
    "check_facet_types": [str(element_type) for element_type in check_facet_types],
    "check_facet_arrays": check_facet_arrays,
    "physical_element_groups": [
        {"name": group.getName(), "element_count": len(group)}
        for group in mesh.elementGroups()
    ],
    "element_group_api": dir(aka.ElementGroup),
    "element_type_names": aka.getElementTypes(),
}
serialized = json.dumps(result, indent=2, default=str) + "\n"
(ROOT / "preflight.json").write_text(serialized, encoding="utf-8")
(REVIEW / "preflight.json").write_text(serialized, encoding="utf-8")
print(serialized, end="")
