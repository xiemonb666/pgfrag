#!/usr/bin/env python3
"""Zero-step native API capture for the frozen 96-T4 coupon.

The extrinsic mode reproduces the old driver's model initialization and facet
selection but does not solve, set boundary conditions, test insertion stress,
assemble internal forces, or evolve the model.  The optional intrinsic mode is
a separate fresh initialization used only to inspect allocated cohesive fields.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import inspect
import json
import math
import os
from pathlib import Path
import platform
import sys
import time
import traceback

import numpy as np


ROOT = Path(__file__).resolve().parents[4]
WORK = ROOT / "pgfrag/experiments/continuation_20261016/extrinsic_capture_preflight"
DATA = ROOT / "pgfrag/data/continuation_20261016/extrinsic_capture_preflight"
RUN = ROOT / "relicweave-synlab/tmp/review_20261004/implementation/cohesive_event_probe_20261004/coarse654_20261004_01"
INPUT = RUN / "snapshot/input"
MANIFEST = RUN / "snapshot/manifest.json"
OLD_DRIVER = ROOT / "pgfrag/experiments/cohesive_event_probe_20261004/probe_driver.py"
WHEEL = ROOT / "third_party/akantu-runtime/downloads/akantu-5.0.7.post1-cp312-cp312-manylinux_2_17_x86_64.manylinux2014_x86_64.whl"
AKANTU_SOURCE = ROOT / "third_party/akantu-runtime/source-v5.0.7"
OUTPUT_LIMIT_BYTES = 16 * 1024 * 1024

SOURCE_FILES = (
    "src/model/solid_mechanics/solid_mechanics_model.cc",
    "src/model/solid_mechanics/solid_mechanics_model_cohesive/solid_mechanics_model_cohesive.cc",
    "src/model/solid_mechanics/solid_mechanics_model_cohesive/materials/material_cohesive.cc",
    "src/model/solid_mechanics/solid_mechanics_model_cohesive/materials/constitutive_laws/material_cohesive_linear.cc",
    "src/model/solid_mechanics/solid_mechanics_model_cohesive/materials/constitutive_laws/material_cohesive_linear.hh",
)

FIELD_SPECS = (
    ("sigma_c_eff", "cohesive", "sigma_c_eff", 1,
     "per-cohesive-quadrature-point effective critical stress; registered by MaterialCohesiveLinear"),
    ("delta_c_eff", "cohesive", "delta_c_eff", 1,
     "per-cohesive-quadrature-point effective critical displacement; populated by insertion checks"),
    ("insertion_stress", "cohesive", "insertion_stress", 3,
     "stress vector recorded at insertion; populated by insertion checks"),
    ("normal", "cohesive", "normal", 3,
     "cohesive integration-point geometric normal, as registered by MaterialCohesive"),
    ("delta_max", "cohesive", "delta max", 1,
     "maximum effective opening history; native internal field name includes a space"),
    ("damage", "cohesive", "damage", 1,
     "cohesive scalar damage internal variable"),
    ("opening", "cohesive", "opening", 3,
     "current cohesive opening/separation vector, updated during cohesive traction evaluation"),
    ("tractions", "cohesive", "tractions", 3,
     "current cohesive traction vector, updated during cohesive traction evaluation"),
    ("contact_tractions", "cohesive", "contact_tractions", 3,
     "penalty-contact traction vector, separate from cohesive traction"),
    ("contact_opening", "cohesive", "contact_opening", 3,
     "penalty-contact opening vector, separate from cohesive opening"),
    ("reversible_energy", "cohesive", "reversible_energy", 1,
     "native cohesive reversible-energy internal field"),
    ("total_energy", "cohesive", "total_energy", 1,
     "native cohesive total-energy internal field"),
    ("bulk_stress", "bulk", "stress", 9,
     "bulk material stress tensor at integration points; read before any force assembly"),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_hash(path: Path) -> dict:
    return {
        "path": path.resolve().relative_to(ROOT.resolve()).as_posix(),
        "exists": path.is_file(),
        "bytes": path.stat().st_size if path.is_file() else None,
        "sha256": sha256(path) if path.is_file() else None,
    }


def jsonable(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        value = float(value)
    if isinstance(value, float) and not math.isfinite(value):
        return repr(value)
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    return repr(value)


def write_json(path: Path, value) -> None:
    path.write_text(json.dumps(jsonable(value), indent=2, sort_keys=True) + "\n", encoding="utf-8")


def describe_callable(obj, name: str) -> dict:
    member = getattr(obj, name, None)
    if not callable(member):
        return {"available": False, "attribute_type": type(member).__name__ if member is not None else None}
    try:
        signature = str(inspect.signature(member))
    except Exception as error:
        signature = None
        signature_error = repr(error)
    else:
        signature_error = None
    doc = inspect.getdoc(member)
    return {
        "available": True,
        "signature": signature,
        "signature_error": signature_error,
        "doc": doc[:1600] if doc else None,
    }


def array_summary(array: np.ndarray) -> dict:
    if np.issubdtype(array.dtype, np.number):
        finite = np.isfinite(array)
    elif np.issubdtype(array.dtype, np.bool_):
        finite = np.ones(array.shape, dtype=bool)
    else:
        finite = np.zeros(array.shape, bool)
    finite_values = array[finite]
    return {
        "shape": list(array.shape),
        "dtype": str(array.dtype),
        "size": int(array.size),
        "finite_count": int(finite.sum()),
        "finite_min": float(finite_values.min()) if finite_values.size else None,
        "finite_max": float(finite_values.max()) if finite_values.size else None,
        "finite_max_abs": float(np.max(np.abs(finite_values))) if finite_values.size else None,
        "finite_l2_norm": float(np.linalg.norm(finite_values)) if finite_values.size else None,
    }


def copy_native_array(arrays: dict, logical_name: str, getter) -> dict:
    try:
        value = np.asarray(getter()).copy()
        arrays[logical_name] = value
        return {"status": "ok", **array_summary(value)}
    except Exception as error:
        return {"status": "error", "error": repr(error)}


def mesh_facts(mesh, aka) -> dict:
    return {
        "nodes": int(mesh.getNbNodes()),
        "tetrahedron_4": int(mesh.getNbElement(type=aka._tetrahedron_4)),
        "cohesive_3d_6": int(mesh.getNbElement(type=aka._cohesive_3d_6)),
    }


def physical_interface_geometry(mesh_path: Path) -> dict:
    """Recover directed normals/areas from the saved Gmsh 2.2 surface elements."""
    lines = mesh_path.read_text(encoding="utf-8").splitlines()
    physical_names: dict[tuple[int, int], str] = {}
    nodes: dict[int, np.ndarray] = {}
    triangles: list[tuple[int, tuple[int, int, int]]] = []
    i = 0
    while i < len(lines):
        if lines[i] == "$PhysicalNames":
            count = int(lines[i + 1]); i += 2
            for line in lines[i:i + count]:
                dim, tag, name = line.split(maxsplit=2)
                physical_names[(int(dim), int(tag))] = name.strip().strip('"')
            i += count
        elif lines[i] == "$Nodes":
            count = int(lines[i + 1]); i += 2
            for line in lines[i:i + count]:
                values = line.split()
                nodes[int(values[0])] = np.asarray([float(v) for v in values[1:4]], dtype=np.float64)
            i += count
        elif lines[i] == "$Elements":
            count = int(lines[i + 1]); i += 2
            for line in lines[i:i + count]:
                values = [int(v) for v in line.split()]
                element_type, tag_count = values[1], values[2]
                tags = values[3:3 + tag_count]
                conn = values[3 + tag_count:]
                if element_type == 2 and len(conn) == 3 and tags:
                    physical = physical_names.get((2, tags[0]))
                    if physical == "interface":
                        triangles.append((tags[0], tuple(conn)))
            i += count
        i += 1

    records = []
    for physical_tag, conn in triangles:
        xyz = np.stack([nodes[node_id] for node_id in conn])
        vector = np.cross(xyz[1] - xyz[0], xyz[2] - xyz[0])
        norm = float(np.linalg.norm(vector))
        records.append({
            "physical_tag": physical_tag,
            "node_ids_in_msh_order": list(conn),
            "area_m2": 0.5 * norm,
            "directed_unit_normal_from_msh_triangle_order": (vector / norm).tolist() if norm else None,
        })
    return {
        "source": "saved Gmsh 2.2 interface triangles and node coordinates, directed by file connectivity order",
        "triangle_count": len(records),
        "total_area_m2": float(sum(item["area_m2"] for item in records)),
        "triangles": records,
        "normal_caveat": "This directed mesh-surface normal is a geometric fallback; only a native cohesive `normal` field can establish the post-insertion material orientation.",
    }


def model_and_material_api(model, materials) -> dict:
    relevant = (
        "flattenInternalReal", "getAcceleration", "getExternalForce", "getInternalForce",
        "getMass", "getDisplacement", "getVelocity", "getBlockedDOFs", "getBlockedDofs",
        "getTime", "getStress", "getNormals", "getNormalsOnIntegrationPoints",
        "getFEEngine", "getElementInserter", "getNbMaterials", "getMaterial", "getMaterials",
        "assembleMassLumped", "assembleInternalForces",
    )
    model_api = {name: describe_callable(model, name) for name in relevant}
    model_api["other_public_stress_normal_methods"] = {
        name: describe_callable(model, name)
        for name in sorted(n for n in dir(model) if not n.startswith("_") and
                           ("stress" in n.lower() or "normal" in n.lower()))
        if name not in model_api
    }
    model_api["other_public_time_methods"] = {
        name: describe_callable(model, name)
        for name in sorted(n for n in dir(model) if not n.startswith("_") and "time" in n.lower())
        if name not in model_api
    }
    material_api = []
    for index, material in enumerate(materials):
        names = [name for name in dir(material) if not name.startswith("_") and
                 any(key in name.lower() for key in ("stress", "sigma", "delta", "traction", "opening", "damage", "normal", "energy"))]
        material_api.append({
            "index": index,
            "type": type(material).__name__,
            "public_field_related_methods": {name: describe_callable(material, name) for name in sorted(names)},
            "requested_cohesive_getter_availability": {
                name: describe_callable(material, name)
                for name in ("getSigmaCEff", "getDeltaCEff", "getInsertionStress", "getNormal",
                             "getDeltaMax", "getDamage", "getOpening", "getTraction",
                             "getContactTraction", "getContactOpening", "getStress")
            },
        })
    return {"model": model_api, "materials": material_api}


def get_materials(model) -> tuple[list, dict]:
    access = {"status": "not_exposed", "count": None, "getter": None, "errors": []}
    count_fn = getattr(model, "getNbMaterials", None)
    if not callable(count_fn):
        count_fn = getattr(model, "getNbMaterial", None)
    count = None
    if callable(count_fn):
        try:
            count = int(count_fn())
            access["count"] = count
        except Exception as error:
            access["errors"].append({"method": getattr(count_fn, "__name__", "getNbMaterials"), "error": repr(error)})
    getter = getattr(model, "getMaterial", None)
    if callable(getter):
        access["getter"] = "getMaterial"
        handles = []
        for index in range(count if count is not None else 4):
            try:
                handles.append(getter(index))
            except Exception as error:
                access["errors"].append({"method": "getMaterial", "index": index, "error": repr(error)})
                break
        if handles:
            access["status"] = "observed"
            return handles, access
    get_all = getattr(model, "getMaterials", None)
    if callable(get_all):
        access["getter"] = "getMaterials"
        try:
            result = get_all()
            if isinstance(result, (list, tuple)):
                access["status"] = "observed"
                return list(result), access
            access["errors"].append({"method": "getMaterials", "return_type": type(result).__name__})
        except Exception as error:
            access["errors"].append({"method": "getMaterials", "error": repr(error)})
    return [], access


def attempt_flatten(model, aka, category: str, native_name: str, elem_type, expected_components: int,
                    semantics: str, arrays: dict) -> dict:
    result = {
        "native_field_name": native_name,
        "expected_components_from_cached_source": expected_components,
        "semantics_from_cached_source": semantics,
        "element_type": str(elem_type),
        "status": "not_attempted",
    }
    flatten = getattr(model, "flattenInternalReal", None)
    if not callable(flatten):
        result.update(status="unavailable", reason="public model.flattenInternalReal not exposed")
        return result
    try:
        element_kind = aka._ek_cohesive if category == "cohesive" else aka._ek_regular
        selected = flatten(native_name, element_kind, aka._not_ghost)
        value = np.asarray(selected(elem_type, aka._not_ghost)).copy()
        key = "field__" + ("bulk_stress" if category == "bulk" else native_name.replace(" ", "_"))
        arrays[key] = value
        result.update(status="ok", array_key=key, **array_summary(value))
        if value.size == 0:
            result["population_status"] = "getter_returned_empty_array; no values are available for constitutive interpretation"
        elif value.ndim == 2 and value.shape[1] != expected_components:
            result["population_status"] = "nonempty_shape_component_count_differs_from_cached-source_expectation"
        else:
            result["population_status"] = "array_available; values are the zero-step initialized state only"
    except Exception as error:
        result.update(status="error", error=repr(error))
    return result


def read_field_arrays(model, aka, arrays: dict) -> dict:
    result = {}
    for logical, category, native_name, components, semantics in FIELD_SPECS:
        elem_type = aka._cohesive_3d_6 if category == "cohesive" else aka._tetrahedron_4
        result[logical] = attempt_flatten(model, aka, category, native_name, elem_type, components, semantics, arrays)
    return result


def read_public_material_getters(materials, aka, arrays: dict) -> dict:
    """Call only documented field getters; retain their signatures in public_api."""
    result = {}
    if materials:
        bulk_get_stress = getattr(materials[0], "getStress", None)
        if callable(bulk_get_stress):
            try:
                value = np.asarray(bulk_get_stress(aka._tetrahedron_4, aka._not_ghost)).copy()
                arrays["material__bulk_getStress"] = value
                result["bulk_material_getStress"] = {
                    "status": "ok", "method": "MaterialLinearElastic3D.getStress(_tetrahedron_4, _not_ghost)",
                    **array_summary(value),
                }
            except Exception as error:
                result["bulk_material_getStress"] = {"status": "error", "error": repr(error)}
    if len(materials) > 1:
        material = materials[1]
        for method_name in ("getOpening", "getTraction"):
            method = getattr(material, method_name, None)
            if not callable(method):
                result["cohesive_material_" + method_name] = {"status": "unavailable"}
                continue
            try:
                value = np.asarray(method(aka._cohesive_3d_6, aka._not_ghost)).copy()
                arrays["material__" + method_name] = value
                result["cohesive_material_" + method_name] = {
                    "status": "ok", "method": f"{type(material).__name__}.{method_name}(_cohesive_3d_6, _not_ghost)",
                    **array_summary(value),
                }
            except Exception as error:
                result["cohesive_material_" + method_name] = {"status": "error", "error": repr(error)}
    return result


def run(mode: str) -> int:
    import akantu as aka

    DATA.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    mesh_path = INPUT / "block_2x2x4mm.msh"
    material_path = INPUT / "material.dat"
    aka.parseInput(str(material_path))
    mesh = aka.Mesh(3)
    mesh.read(str(mesh_path))
    model = aka.SolidMechanicsModelCohesive(mesh)
    if mode == "intrinsic-capability":
        # Separate fresh run: register just the saved physical interface before
        # initFull, allowing native intrinsic initialization to create its fields.
        inserter = model.getElementInserter()
        inserter.addPhysicalSurface("interface")
        model.initFull(_analysis_method=aka._explicit_lumped_mass, _is_extrinsic=False)
    else:
        # Match old probe_driver.build_model ordering exactly.
        model.initFull(_analysis_method=aka._explicit_lumped_mass, _is_extrinsic=True)
    model.assembleMassLumped()

    eligible = None
    eligible_total = None
    already_inserted = None
    eligible_mask = None
    if mode == "extrinsic":
        inserter = model.getElementInserter()
        inserter.addPhysicalSurface("interface")
        model.updateAutomaticInsertion()
        checks = np.asarray(inserter.getCheckFacets()(aka._triangle_3)).copy()
        eligible_mask = checks.reshape(-1).astype(bool, copy=True)
        inserted = np.asarray(inserter.getInsertionFacets()(aka._triangle_3)).copy()
        eligible = int(np.count_nonzero(checks))
        eligible_total = int(checks.size)
        already_inserted = int(np.count_nonzero(inserted))

    arrays: dict[str, np.ndarray] = {}
    native_getter_names = (
        "getAcceleration", "getExternalForce", "getInternalForce", "getMass",
        "getDisplacement", "getVelocity", "getBlockedDOFs", "getBlockedDofs",
    )
    native_getters = {}
    for method_name in native_getter_names:
        if method_name in native_getters:
            continue
        method = getattr(model, method_name, None)
        if not callable(method):
            continue
        captured = copy_native_array(arrays, "native__" + method_name, method)
        captured["method"] = method_name
        native_getters[method_name] = captured
    if "getBlockedDOFs" not in native_getters and "getBlockedDofs" not in native_getters:
        native_getters["blocked_dofs"] = {"status": "unavailable", "reason": "neither documented capitalization is exposed"}
    arrays["native__nodes"] = np.asarray(mesh.getNodes()).copy()

    materials, material_access = get_materials(model)
    api = model_and_material_api(model, materials)
    fields = read_field_arrays(model, aka, arrays)
    if mode == "intrinsic-capability":
        fields["sigma_c_eff"]["value_interpretation"] = "zero-initialized here; no extrinsic checkInsertion call occurred"
        fields["delta_c_eff"]["value_interpretation"] = "equals the explicit nominal delta_c default in this initialization; not an insertion-derived effective value"
        fields["insertion_stress"]["value_interpretation"] = "zero-initialized here; this separate intrinsic setup did not capture an extrinsic insertion stress"
        fields["normal"]["value_interpretation"] = "zero cohesive internal field before traction evaluation; use the separately captured FacetsFEEngine geometric normals for orientation"
        fields["opening"]["value_interpretation"] = "zero opening in the unloaded initial intrinsic setup; getter capability only"
        fields["tractions"]["value_interpretation"] = "zero traction in the unloaded initial intrinsic setup; getter capability only"
    material_getter_captures = read_public_material_getters(materials, aka, arrays)

    # Try the documented FEEngine normal accessor only if exposed.  This is a
    # geometry read; do not call any force or constitutive assembly routine.
    facet_normal_api = {"status": "unavailable"}
    get_stress_on_facets_result = {"status": "unavailable"}
    get_fe = getattr(model, "getFEEngine", None)
    if callable(get_fe):
        try:
            facet_engine = get_fe("FacetsFEEngine")
            facet_normal_api = describe_callable(facet_engine, "getNormalsOnIntegrationPoints")
            normal_getter = getattr(facet_engine, "getNormalsOnIntegrationPoints", None)
            if callable(normal_getter):
                try:
                    value = np.asarray(normal_getter(aka._triangle_3, aka._not_ghost)).copy()
                    arrays["native__facet_normals_all_triangle_3"] = value
                    facet_normal_api.update(status="ok", native_mesh_scope="all non-ghost triangle facets in FacetsFEEngine",
                                            **array_summary(value))
                    facet_normal_api["interface_filter_mapping"] = "the full array was not joined to node IDs; the eligible subset below uses the native check-facet mask and same triangle ElementTypeMap ordering"
                    if eligible_mask is not None and value.ndim == 2 and value.shape[1] == 3 and eligible_mask.size:
                        if value.shape[0] % eligible_mask.size == 0:
                            qps_per_facet = value.shape[0] // eligible_mask.size
                            grouped = value.reshape(eligible_mask.size, qps_per_facet, 3)
                            selected = grouped[eligible_mask].copy()
                            arrays["native__facet_normals_eligible_interface"] = selected.reshape(-1, 3).copy()
                            facet_normal_api["eligible_mask_join"] = {
                                "status": "shape_aligned_by_same_triangle_3 ElementTypeMap ordering; no exact node-ID join was attempted",
                                "eligible_facets": int(eligible_mask.sum()),
                                "qps_per_facet": int(qps_per_facet),
                                "selected_shape": list(selected.shape),
                                "selected_normals_by_facet_mean": selected.mean(axis=1).tolist(),
                            }
                        else:
                            facet_normal_api["eligible_mask_join"] = {
                                "status": "unavailable_shape_mismatch",
                                "eligible_facets": int(eligible_mask.sum()),
                                "mask_slots": int(eligible_mask.size),
                            }
                except Exception as error:
                    facet_normal_api.update(status="error", error=repr(error))
        except Exception as error:
            facet_normal_api = {"status": "error", "error": repr(error)}

    get_stress_on_facets = getattr(model, "getStressOnFacets", None)
    if callable(get_stress_on_facets):
        try:
            value = np.asarray(get_stress_on_facets(aka._triangle_3)).copy()
            arrays["native__stress_on_facets"] = value
            get_stress_on_facets_result = {
                "status": "ok",
                "method": "SolidMechanicsModelCohesive.getStressOnFacets(_triangle_3)",
                "semantics": "public cohesive-model getter for stress values on triangular facets; this zero-step capture does not call checkCohesiveStress or recompute/assemble stress",
                **array_summary(value),
            }
        except Exception as error:
            get_stress_on_facets_result = {"status": "error", "error": repr(error)}

    # Read initial native clock if documented; no call changes it.
    native_time = None
    get_time = getattr(model, "getTime", None)
    if callable(get_time):
        try:
            native_time = {"status": "ok", "value": float(get_time())}
        except Exception as error:
            native_time = {"status": "error", "error": repr(error)}
    else:
        native_time = {"status": "unavailable", "reason": "actual wheel model exposes no public getTime() method"}

    mode_key = "extrinsic" if mode == "extrinsic" else "intrinsic_capability"
    array_path = DATA / f"{mode_key}_native_arrays.npz"
    np.savez_compressed(array_path, **arrays)
    source_files = [file_hash(AKANTU_SOURCE / rel) for rel in SOURCE_FILES]
    input_files = [file_hash(mesh_path), file_hash(material_path), file_hash(INPUT / "generated_geometry.json"),
                   file_hash(MANIFEST), file_hash(OLD_DRIVER)]
    runtime = {
        "python_executable": sys.executable,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "akantu_version": getattr(aka, "__version__", None),
        "akantu_distribution_version": importlib.metadata.version("akantu"),
        "numpy_version": np.__version__,
        "omp_num_threads": os.environ.get("OMP_NUM_THREADS"),
        "openblas_num_threads": os.environ.get("OPENBLAS_NUM_THREADS"),
        "mkl_num_threads": os.environ.get("MKL_NUM_THREADS"),
        "runtime_wheel": file_hash(WHEEL),
    }
    no_solver_steps = {
        "solver_steps": 0,
        "boundary_conditions_applied": 0,
        "stress_checks": 0,
        "automatic_insertions_during_capture": 0,
        "intrinsic_cohesive_elements_created_during_init": int(mesh.getNbElement(type=aka._cohesive_3d_6)) if mode == "intrinsic-capability" else 0,
        "internal_force_assemblies_after_initialization": 0,
        "mass_assembly_calls_after_initFull": 1,
        "post_step_mass_reassemblies": 0,
        "evolution_calls": 0,
        "note": "assembleMassLumped was called once after initFull to match the old probe setup. In the separate intrinsic capability run, initFull itself created 8 cohesive elements before this initial mass assembly. No solve, load, BC, stress test, explicit step, or force assembly was called.",
    }
    result = {
        "schema": "pgfrag.extrinsic_capture_preflight.v1",
        "mode": mode,
        "status": "CAPTURED_ZERO_STEP_API_PREFLIGHT",
        "runtime": runtime,
        "frozen_inputs": input_files,
        "cached_native_sources": source_files,
        "source_distinction": "The installed runtime is the 5.0.7.post1 wheel. The cached inspected source is tagged v5.0.7 and is recorded separately; it is not represented as the wheel's exact source.",
        "missing_cached_source": {
            "path": "src/model/solid_mechanics/solid_mechanics_model_cohesive/materials/constitutive_laws/material_cohesive_linear_inline_impl.hh",
            "status": "not present in local source cache; no source-hash claim",
        },
        "native_steps_and_mutations": no_solver_steps,
        "native_time_at_capture": native_time,
        "native_mesh_counts": mesh_facts(mesh, aka),
        "interface_geometry": physical_interface_geometry(mesh_path),
        "extrinsic_facet_selection": {
            "eligible_interface_facets": eligible,
            "eligible_face_slots": eligible_total,
            "already_inserted_facets": already_inserted,
            "expected_eligible_from_frozen_manifest": manifest.get("expected_insertion_facets"),
            "interpretation": "Only populated for mode=extrinsic. Intrinsic mode is a separate initial-condition capability probe, not an extrinsic insertion/event record.",
        },
        "material_access": material_access,
        "public_api": api,
        "field_capture": fields,
        "public_material_getter_captures": material_getter_captures,
        "normal_access_fallback": {
            "native_cohesive_field": fields["normal"],
            "facets_feengine_getNormalsOnIntegrationPoints": facet_normal_api,
            "cohesive_model_getStressOnFacets": get_stress_on_facets_result,
            "geometry_recovery": "Saved Gmsh 2.2 physical interface triangles: area and directed unit normal are computed from actual node coordinates and triangle connectivity in file order.",
        },
        "native_array_copies": native_getters,
        "array_archive": {
            "path": array_path.resolve().relative_to(ROOT.resolve()).as_posix(),
            "bytes": array_path.stat().st_size,
            "sha256": sha256(array_path),
            "keys": {name: array_summary(value) for name, value in arrays.items()},
        },
        "capture_plan": {
            "first_preflight_result": "The current run establishes public API availability and initial shapes only; its zero stress/force values are not insertion evidence.",
            "future_native_step_capture": [
                "The frozen driver order is A before applyBC; B immediately after applyBC and before checkCohesiveStress; C immediately after checkCohesiveStress/insertion and before solveStep; D after solveStep.",
                "Before each step, record the driver's prescribed time step*dt and the applied BC value/vector. The 5.0.7.post1 wheel exposes no public model.getTime(), so do not claim a native model clock; retain the driver's step/dt clock and any separate time_before field as unavailable.",
                "At B, capture displacement, velocity, acceleration, mass, blocked DOFs, external/internal force, bulk stress, facet stress, and cohesive fields before insertion. applyBC occurs first; it changes the prescribed displacement, while the existing bulk stress is not recomputed by that setter.",
                "checkCohesiveStress interpolates the already available bulk stress to candidate facets; capture both the pre-check bulk/facet stress and the B displacement. Do not call assembleInternalForces merely to observe state.",
                "At C, immediately after insertion but before solveStep/reassembly, capture duplicated node IDs and public parent mapping if exposed, displacement, velocity, acceleration, blocked DOFs, and the appended zero-mass rows.",
                "At D, after solveStep and its operative mass update, capture the same arrays and native cohesive fields again. Attach phase label, step, driver's prescribed time, applied BC, and observed node/cohesive counts.",
            ],
            "insertion_field_interpretation": "An intrinsic initialized field array proves getter capability only. In this run delta_c_eff equals the nominal explicit delta_c default, while sigma_c_eff, insertion_stress, and cohesive normal are zero before any check/traction evaluation. It does not prove the old extrinsic checkCohesiveStress insertion values, their update phase, or the wheel-specific 651 law.",
        },
        "elapsed_seconds": time.monotonic() - started,
        "limits": {"wall_seconds": 60, "total_output_bytes": OUTPUT_LIMIT_BYTES},
    }
    result["script"] = file_hash(Path(__file__).resolve())
    output_path = DATA / f"{mode_key}_capture.json"
    write_json(output_path, result)

    # Verify our deliverables remain well below the task's total output cap.
    produced = [p for p in DATA.glob(f"{mode_key}_*") if p.is_file()]
    output_bytes = sum(p.stat().st_size for p in produced)
    folder_output_bytes = sum(p.stat().st_size for p in DATA.iterdir() if p.is_file())
    if folder_output_bytes > OUTPUT_LIMIT_BYTES:
        raise RuntimeError(f"output cap exceeded: {folder_output_bytes} > {OUTPUT_LIMIT_BYTES}")
    print(json.dumps({
        "status": result["status"], "mode": mode, "solver_steps": 0,
        "native_mesh_counts": result["native_mesh_counts"],
        "eligible_interface_facets": eligible,
        "already_inserted_facets": already_inserted,
        "field_status": {k: v["status"] for k, v in fields.items()},
        "array_archive_bytes": array_path.stat().st_size,
        "output_bytes": output_bytes,
        "data_folder_bytes": folder_output_bytes,
        "elapsed_seconds": result["elapsed_seconds"],
    }, sort_keys=True))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("extrinsic", "intrinsic-capability"), required=True)
    args = parser.parse_args()
    try:
        return run(args.mode)
    except Exception:
        traceback.print_exc()
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
