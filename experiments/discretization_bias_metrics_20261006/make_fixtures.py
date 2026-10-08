"""Generate evaluator-only fixture descriptors. These are not solver results."""
from __future__ import annotations
import copy, json
from pathlib import Path

ROOT=Path(__file__).parent/"fixtures"


def base(run_id):
    return {
      "schema":"pgfrag.discretization_bias.descriptor.v1",
      "run_id":run_id,
      "evidence":{"kind":"FIXTURE","scope":"FIXTURE","physical_result":False},
      "physical_frame":{
        "common_frame_id":"fixture-VF-rest-frame-v1","coordinate_units":"mm",
        "local_to_common":[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]},
      "geometry":{"reference":"fixture://vf-mother/cup-wall-notch","sha256":"fixture-geometry-hash-0001","units":"mm"},
      "load":{"sha256":"fixture-load-hash-0001","units":"m/s","description":"FIXTURE grip waveform"},
      "material":{"parameter_sha256":"fixture-material-parameters-0001",
        "field_sha256":"fixture-rest-field-hash-0001","random_field_seed":17,
        "units":{"E":"Pa","Gc":"J/m^2","density":"kg/m^3"}},
      "solver":{"build_sha256":"fixture-build-hash-0001","failure_law_sha256":"fixture-law-hash-0001",
        "surface_correction":"FIXTURE_NoCorrection","run_status":"COMPLETE"},
      "discretization":{"mesh_sha256":"fixture-mesh-"+run_id,"h":0.5,"h_units":"mm",
        "delta":1.5,"delta_units":"mm","dt":1e-8,"dt_units":"s",
        "orientation_id":"orientation-0","translation_id":"translation-0"},
      "fracture":{"fracture_occurred":True,"history_complete":True,
        "initiation_time":0.002,"initiation_time_units":"s",
        "surface_extractor_sha256":"fixture-surface-extractor-0001",
        "path_sampling_policy_sha256":"fixture-fixed-spacing-samples-0001",
        "path_units":"mm","path_samples_local":[[0,0,0],[1,0,0],[2,0,0]],
        "surface":{"coordinate_units":"mm",
          "vertices":[{"id":"v0","xyz":[0,0,0]},{"id":"v1","xyz":[1,0,0]},
                      {"id":"v2","xyz":[1,1,0]},{"id":"v3","xyz":[0,1,0]}],
          "triangles":[["v0","v1","v2"],["v0","v2","v3"]]}},
      "fragment_evidence":{"kind":"closed_surface_geometric_components","parent_mass":1.0,
        "parent_mass_units":"kg","validation":{"closed_surfaces":True,"pairwise_nonoverlap":True,
        "parent_coverage_error":0.01},
        "components":[{"id":"frag-a","mass":0.6,"mass_units":"kg","surface_sha256":"fixture-frag-a-surface"},
                      {"id":"frag-b","mass":0.4,"mass_units":"kg","surface_sha256":"fixture-frag-b-surface"}]},
      "energy_terms":{"external_work":{"value":1.0,"units":"J"},
        "kinetic":{"value":0.1,"units":"J"},"elastic_strain":{"value":0.8,"units":"J"},
        "fracture_dissipation":{"value":0.08,"units":"J"},"residual":{"value":0.02,"units":"J"}}
    }


def put(name,value):
    (ROOT/name).write_text(json.dumps(value,indent=2)+"\n",encoding="utf-8")


def main():
    ROOT.mkdir(parents=True,exist_ok=True)
    b=base("fixture_baseline"); put("baseline.json",b)
    # A 90-degree local frame. Local points are inverse-rotated so their
    # transported coordinates reproduce the baseline specimen-frame geometry.
    r=copy.deepcopy(b); r["run_id"]="fixture_rotated"
    r["solver"]["build_sha256"]="fixture-build-hash-0001"
    r["discretization"].update(mesh_sha256="fixture-mesh-rotated-0002",orientation_id="orientation-r90")
    r["physical_frame"]["local_to_common"]=[[0,-1,0,0],[1,0,0,0],[0,0,1,0],[0,0,0,1]]
    r["fracture"]["path_samples_local"]=[[0,0,0],[0,-1,0],[0,-2,0]]
    r["fracture"]["surface"]["vertices"]=[
      {"id":"r0","xyz":[0,0,0]},{"id":"r1","xyz":[0,-1,0]},
      {"id":"r2","xyz":[1,-1,0]},{"id":"r3","xyz":[1,0,0]}]
    r["fracture"]["surface"]["triangles"]=[["r0","r1","r2"],["r0","r2","r3"]]
    put("rotated.json",r)
    # Physical path/area/mass values vary slightly, while only translation
    # placement changes in the mesh factor metadata.
    t=copy.deepcopy(b); t["run_id"]="fixture_translated"
    t["discretization"].update(mesh_sha256="fixture-mesh-translated-0003",translation_id="translation-half-cell")
    t["physical_frame"]["local_to_common"]=[[1,0,0,10],[0,1,0,20],[0,0,1,5],[0,0,0,1]]
    t["fracture"]["path_samples_local"]=[[-10,-20,-5],[-8.98,-20,-5],[-7.96,-20,-5]]
    t["fracture"]["surface"]["vertices"]=[
      {"id":"t0","xyz":[-10,-20,-5]},{"id":"t1","xyz":[-8.98,-20,-5]},
      {"id":"t2","xyz":[-8.98,-19,-5]},{"id":"t3","xyz":[-10,-19,-5]}]
    t["fracture"]["surface"]["triangles"]=[["t0","t1","t2"],["t0","t2","t3"]]
    t["fracture"]["initiation_time"]=0.0021
    t["fragment_evidence"]["validation"]["parent_coverage_error"]=0.02
    t["fragment_evidence"]["components"][0]["mass"]=0.65
    t["fragment_evidence"]["components"][1]["mass"]=0.35
    t["energy_terms"]["external_work"]["value"]=1.03
    t["energy_terms"]["kinetic"]["value"]=0.11
    t["energy_terms"]["elastic_strain"]["value"]=0.81
    t["energy_terms"]["fracture_dissipation"]["value"]=0.09
    t["energy_terms"]["residual"]["value"]=0.02
    put("translated.json",t)
    h=copy.deepcopy(b); h["run_id"]="fixture_h_seed_changed"
    h["discretization"].update(mesh_sha256="fixture-mesh-h0400-0004",h=0.4)
    h["material"].update(random_field_seed=18,field_sha256="fixture-rest-field-hash-0002")
    put("h_seed_changed.json",h)
    n=copy.deepcopy(b); n["run_id"]="fixture_no_fracture"
    n["discretization"].update(mesh_sha256="fixture-mesh-no-fracture-0005",orientation_id="orientation-r45")
    n["fracture"]={"fracture_occurred":False,"history_complete":True,
      "surface_extractor_sha256":"fixture-surface-extractor-0001",
      "path_sampling_policy_sha256":"fixture-fixed-spacing-samples-0001"}
    put("no_fracture.json",n)
    g=copy.deepcopy(r); g["run_id"]="fixture_graph_components"
    g["discretization"]["mesh_sha256"]="fixture-mesh-graph-components-0006"
    g["discretization"]["orientation_id"]="orientation-r90-graph"
    g["fragment_evidence"]={"kind":"bond_graph_components","components":[
        {"id":"component-1","mass":0.6,"mass_units":"kg"},
        {"id":"component-2","mass":0.4,"mass_units":"kg"}]}
    put("graph_only_fragments.json",g)
    tol={}
    for key in ("path_symmetric_p95_mm","path_symmetric_hausdorff_mm","surface_area_abs_delta_mm2",
      "initiation_time_abs_delta_s","fragment_mass_w1_fraction","fragment_mass_closure_abs_delta",
      "reference_fragment_mass_closure_error","candidate_fragment_mass_closure_error",
      "reference_parent_coverage_error","candidate_parent_coverage_error",
      "energy.external_work_abs_delta_J","energy.kinetic_abs_delta_J",
      "energy.elastic_strain_abs_delta_J","energy.fracture_dissipation_abs_delta_J","energy.residual_abs_delta_J"):
        tol[key]={"max_abs":1.0,"basis":"FIXTURE-only plumbing tolerance; not a scientific threshold"}
    pairs=[]
    for pid,cid,factor in (("rotate","fixture_rotated","orientation"),
       ("translate","fixture_translated","translation"),("h-seed-mismatch","fixture_h_seed_changed","h"),
       ("no-fracture","fixture_no_fracture","orientation"),("graph-only","fixture_graph_components","orientation")):
        pairs.append({"pair_id":pid,"reference_run_id":"fixture_baseline",
          "candidate_run_id":cid,"factor":factor,"tolerances":tol})
    study={"schema":"pgfrag.discretization_bias.study.v1","study_id":"fixture-evaluator-plumbing",
      "descriptors":["baseline.json","rotated.json","translated.json","h_seed_changed.json",
                     "no_fracture.json","graph_only_fragments.json"],"pairs":pairs}
    put("fixture_study.json",study)


if __name__=="__main__": main()

