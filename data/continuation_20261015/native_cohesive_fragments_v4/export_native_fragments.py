"""Export closed native tet boundaries; no arbitrary cut or artificial caps."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"pgfrag/src"))
sys.path.insert(0, str(ROOT/"vesselforge/src"))
from pgfrag.cohesive_fragments import extract_cohesive_fragments
from vesselforge.spherical_bowl import evaluate_bowl_paint_at_xyz


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path, x):
    path.write_text(json.dumps(x, indent=2, allow_nan=False)+"\n", encoding="utf-8")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--frozen", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    run, frozen, out = a.run, a.frozen, a.out
    if out.exists():
        raise FileExistsError(out)
    receipt = json.loads((run/"receipt.json").read_text())
    for name, row in receipt["files"].items():
        path = run/name
        if path.stat().st_size != row["bytes"] or sha(path) != row["sha256"]:
            raise ValueError(f"native receipt mismatch: {name}")
    if json.loads((run/"result.json").read_text())["status"] != "NATIVE_RUN_COMPLETE_DIAGNOSTIC":
        raise ValueError("incomplete native source")
    case = json.loads((frozen/"case.json").read_text())
    manifest = json.loads((frozen/"forge_manifest.json").read_text())
    initial = np.load(run/"initial_state.npz", allow_pickle=False)
    node = np.load(run/"history/node_state.npy", mmap_mode="r", allow_pickle=False)
    damage = np.load(run/"history/damage.npy", mmap_mode="r", allow_pickle=False)
    traction = np.load(run/"history/traction_Pa.npy", mmap_mode="r", allow_pickle=False)
    saved = np.load(run/"history/node_saved_steps.npy", allow_pickle=False)
    lineage = json.loads((run/"lineage.json").read_text())
    rest = initial["rest_nodes_m"]; current = rest+node[-1, :, 0, :]
    facets = len(initial["cohesive"])
    result = extract_cohesive_fragments(rest, initial["tets"], initial["cohesive"],
              damage[-1].reshape(facets, 3), traction[-1].reshape(facets, 3, 3),
              damage_tolerance=case["gates"]["damage_tolerance"],
              traction_tolerance=case["gates"]["traction_zero_tolerance_Pa"],
              current_nodes=current, native_node_to_source=initial["native_to_source"],
              parent_facet_ids=[row["input_face_id"] for row in lineage["parent_facets"]],
              source_exterior_faces=initial["source_exterior_faces"], rho=case["bulk_material"]["rho_kg_m3"])
    if not result["all_components_qualified_closed_mechanical_fragments"]:
        raise ValueError("native topological fragment eligibility failed")
    result["schema"] = "pgfrag.native_cohesive_fragment_topology.v1"
    faces = result["boundary_faces"]
    tri = np.asarray([f["native_node_ids"] for f in faces], dtype=np.int32)
    component = np.asarray([f["component_index"] for f in faces], dtype=np.int32)
    kind = np.asarray([1 if f["surface_kind"] == "failed_cohesive_interface" else 0 for f in faces], dtype=np.uint8)
    parent = np.asarray([f["facet_row"] if f["facet_row"] is not None else -1 for f in faces], dtype=np.int32)
    region = np.zeros(len(tri), dtype=np.uint8)
    x = rest[tri]; center = np.asarray([0., 0., .012])
    radius = np.linalg.norm(x-center, axis=2)
    region[np.all(np.abs(radius-.012) < 1e-8, axis=1)] = 1
    region[np.all(np.abs(radius-.008) < 1e-8, axis=1)] = 2
    region[np.all(np.abs(x[:, :, 2]-.012) < 1e-8, axis=1)] = 3
    region[kind == 1] = 255
    if np.any(region == 0):
        raise ValueError("a native boundary face lacks an analytic source surface class")
    color = np.tile(manifest["material"]["appearance"]["base_albedo_linear_rgb"], (len(tri), 1))
    color[region == 1] = evaluate_bowl_paint_at_xyz(x[region == 1].mean(axis=1)*1000, manifest["ornament"]["recipe"])
    # Bare fracture albedo is a placeholder, independent of damage/IDs.
    color[region == 255] = [.58, .39, .24]
    out.mkdir(parents=True, exist_ok=False)
    np.savez_compressed(out/"state.npz", rest_nodes_m=rest, current_nodes_m=current,
             velocity_mps=node[-1, :, 1, :], tets=initial["tets"], triangles_native=tri,
             face_component=component, face_kind=kind, face_parent_id=parent,
             face_region_code=region, face_color_linear_rgb=color,
             native_to_source=initial["native_to_source"])
    dump(out/"fragments.json", result)
    shutil.copy2(frozen/"forge_manifest.json", out/"forge_manifest.json")
    shutil.copy2(frozen/"case.json", out/"case.json")
    shutil.copy2(run/"receipt.json", out/"receipt.json")
    first_full = int(np.flatnonzero(np.all(damage >= 1-case["gates"]["damage_tolerance"], axis=1))[0])
    separation = {"method": "global x-plane witness over every saved nodal state after all QPs fail; no volume-tag fragment labels",
                  "node_history_stride": case["output"]["node_stride"], "first_all_qp_damage_step": first_full,
                  "full_qp_damage_is_not_alone_a_traction_gate": True,
                  "self_intersection_continuous_time_certificate": False}
    if len(result["components"]) == 2:
        comps = sorted(result["components"], key=lambda c: rest[c["native_node_ids"], 0].mean())
        left, right = [np.asarray(c["native_node_ids"]) for c in comps]
        ids = np.flatnonzero(saved >= first_full)
        min_gap = min(float((rest[right, 0]+node[i, right, 0, 0]).min()-(rest[left, 0]+node[i, left, 0, 0]).max()) for i in ids)
        separation.update({"minimum_saved_global_x_gap_m": min_gap,
                           "all_saved_postdebond_components_x_separated": bool(min_gap > 0),
                           "final_global_x_gap_m": float(current[right, 0].min()-current[left, 0].max())})
    dump(out/"separation.json", separation)
    for c in result["components"]:
        ids = c["native_node_ids"]; index = {n: i+1 for i, n in enumerate(ids)}
        text = "# Native deformed T4 boundary in metres; uncalibrated weak seam\n"
        text += "".join(f"v {current[n,0]:.17g} {current[n,1]:.17g} {current[n,2]:.17g}\n" for n in ids)
        text += "".join("f "+" ".join(str(index[n]) for n in t)+"\n" for t in c["oriented_native_triangles"])
        (out/f"fragment_{c['component_index']:02d}.obj").write_text(text, encoding="utf-8")
    for source in (Path(__file__), ROOT/"pgfrag/src/pgfrag/cohesive_fragments.py"):
        shutil.copy2(source, out/source.name)
    files = {p.name: {"sha256": sha(p), "bytes": p.stat().st_size} for p in out.iterdir() if p.is_file()}
    package = {"schema": "pgfrag.native_cohesive_fragments.v1", "units": "SI", "files": files,
               "source": {"run_receipt_sha256": sha(run/"receipt.json"), "frozen_case_sha256": sha(frozen/"case.json"),
                          "forge_manifest_sha256": sha(frozen/"forge_manifest.json")},
               "qualification": {"status": "PREINSERTED_WEAK_SEAM_UNCALIBRATED", "preinserted_weak_seam": True,
                   "material_calibration_status": "UNCALIBRATED", "physical_pass": False,
                   "closed_native_topology": True, "mesh_sensitivity": "NOT_ESTABLISHED",
                   "typed_R2_qualification": "NOT_ESTABLISHED", "pretrain_eligible": False},
               "counts": {"nodes": len(rest), "tetrahedra": len(initial["tets"]), "fragments": len(result["components"]),
                          "paired_failed_facets": len(result["failed_interface_face_pairs"])},
               "paint": "continuous rest-coordinate recipe; per-face colours are preview only",
               "limitations": ["prescribed weak plane, uncalibrated material", "linear geometry approximates CAD",
                  "no continuous-time self-intersection proof or contact dynamics", "no global physical/PRETRAIN/R2 pass"]}
    dump(out/"manifest.json", package)
    print(json.dumps({"out": str(out), "counts": package["counts"], "separation": separation,
                      "mass_kg": result["total_rest_bulk_mass"], "manifest_sha256": sha(out/"manifest.json")}))


if __name__ == "__main__":
    main()
