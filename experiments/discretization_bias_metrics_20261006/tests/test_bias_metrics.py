from __future__ import annotations
import copy, json, sys, unittest
from pathlib import Path

HERE=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(HERE))
import bias_metrics as bm

FIX=HERE/"fixtures"


def descriptor(name):
    return json.loads((FIX/name).read_text(encoding="utf-8"))


def pair(candidate,factor):
    study=json.loads((FIX/"fixture_study.json").read_text(encoding="utf-8"))
    return next(p for p in study["pairs"] if p["candidate_run_id"]==candidate and p["factor"]==factor)


class BiasMetricsFixtureTests(unittest.TestCase):
    def test_rigid_rotation_is_transported_to_common_physical_frame(self):
        result=bm.compare({"fixture_baseline":descriptor("baseline.json"),
                           "fixture_rotated":descriptor("rotated.json")},
                          pair("fixture_rotated","orientation"))
        self.assertEqual(result["pairing_status"],"COMPARABLE")
        self.assertEqual(result["assessment_status"],"FIXTURE_ONLY")
        self.assertLess(result["metrics"]["path_symmetric_p95_mm"],1e-12)
        self.assertLess(result["metrics"]["surface_area_abs_delta_mm2"],1e-12)
        self.assertEqual(result["metrics"]["initiation_time_abs_delta_s"],0.0)
        self.assertEqual(result["metrics"]["energy.external_work_abs_delta_J"],0.0)

    def test_dt_pair_holds_mesh_and_physical_inputs_fixed(self):
        a=descriptor("baseline.json")
        b=copy.deepcopy(a)
        b["run_id"]="fixture_dt_half"
        b["discretization"]["dt"]=5e-9
        result=bm.compare({a["run_id"]:a,b["run_id"]:b},
                          {"pair_id":"dt-half","reference_run_id":a["run_id"],
                           "candidate_run_id":b["run_id"],"factor":"dt"})
        self.assertEqual(result["pairing_status"],"COMPARABLE")
        self.assertEqual(result["assessment_status"],"FIXTURE_ONLY")

    def test_native_candidate_without_qualification_stays_unknown(self):
        a=descriptor("baseline.json"); b=descriptor("rotated.json")
        for d in (a,b):
            d["evidence"]={"kind":"NATIVE_PHYSICAL_RESULT","scope":"UNQUALIFIED_TEST_ONLY"}
        result=bm.compare({a["run_id"]:a,b["run_id"]:b},
                          {"pair_id":"unqualified","reference_run_id":a["run_id"],
                           "candidate_run_id":b["run_id"],"factor":"orientation",
                           "tolerances":pair("fixture_rotated","orientation")["tolerances"]})
        self.assertEqual(result["pairing_status"],"COMPARABLE")
        self.assertEqual(result["assessment_status"],"UNKNOWN")
        self.assertTrue(any("qualification attestations missing" in x for x in result["reasons"]))

    def test_tolerances_are_not_defaulted(self):
        status,reasons=bm.tolerance_status({"metric_x":0.0},None)
        self.assertEqual(status,"UNKNOWN")
        self.assertTrue(any("no explicit caller tolerances" in x for x in reasons))

    def test_vertex_id_and_path_sample_order_do_not_change_metrics(self):
        a=descriptor("baseline.json")
        b=copy.deepcopy(a)
        b["run_id"]="fixture_reordered"
        b["fracture"]["path_samples_local"].reverse()
        surface=b["fracture"]["surface"]
        rename={"v0":"node-D","v1":"node-C","v2":"node-B","v3":"node-A"}
        for vertex in surface["vertices"]:
            vertex["id"]=rename[vertex["id"]]
        surface["vertices"].reverse()
        surface["triangles"]=[[rename[v] for v in face] for face in reversed(surface["triangles"])]
        metrics,unavailable=bm.metrics(a,b)
        self.assertEqual(unavailable,[])
        self.assertEqual(metrics["path_symmetric_p95_mm"],0.0)
        self.assertEqual(metrics["surface_area_abs_delta_mm2"],0.0)
        self.assertEqual(metrics["fragment_mass_w1_fraction"],0.0)

    def test_h_pair_rejects_changed_material_field_and_seed(self):
        result=bm.compare({"fixture_baseline":descriptor("baseline.json"),
                           "fixture_h_seed_changed":descriptor("h_seed_changed.json")},
                          pair("fixture_h_seed_changed","h"))
        self.assertEqual(result["pairing_status"],"INCOMPARABLE")
        self.assertIn("material.field_sha256", " ".join(result["reasons"]))
        self.assertNotEqual(result["assessment_status"],"WITHIN_TOLERANCES_NOT_BIAS_ELIMINATED")

    def test_no_fracture_is_unknown_not_a_zero_metric_pass(self):
        result=bm.compare({"fixture_baseline":descriptor("baseline.json"),
                           "fixture_no_fracture":descriptor("no_fracture.json")},
                          pair("fixture_no_fracture","orientation"))
        self.assertEqual(result["pairing_status"],"COMPARABLE")
        self.assertEqual(result["assessment_status"],"UNKNOWN")
        self.assertIn("contains no fracture event"," ".join(result["reasons"]))
        self.assertEqual(result["metrics"],{})

    def test_graph_components_are_not_geometric_fragments(self):
        result=bm.compare({"fixture_baseline":descriptor("baseline.json"),
                           "fixture_graph_components":descriptor("graph_only_fragments.json")},
                          pair("fixture_graph_components","orientation"))
        self.assertEqual(result["assessment_status"],"FIXTURE_ONLY")
        self.assertNotIn("fragment_mass_w1_fraction",result["metrics"])
        self.assertTrue(any("not geometric fragment evidence" in x for x in result["unavailable_metrics"]))

    def test_fixture_study_can_never_admit_or_claim_bias_removed(self):
        report=bm.run_study(FIX/"fixture_study.json")
        self.assertEqual(report["overall_status"],"FIXTURE_ONLY")
        self.assertEqual(report["admission"],"NOT_ADMITTED")
        self.assertFalse(report["bias_eliminated"])
        self.assertNotIn("PASS",{r["assessment_status"] for r in report["pair_results"]})


if __name__=="__main__":
    unittest.main()

