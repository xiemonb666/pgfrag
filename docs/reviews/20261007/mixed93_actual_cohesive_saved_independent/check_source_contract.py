"""Stdlib-only source check: read saved JSON contracts, never open NPZ or run math."""
import ast, hashlib, json
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1<<20),b""): h.update(b)
    return h.hexdigest()
def need(ok,msg):
    if not ok: raise SystemExit("SOURCE_FAIL: "+msg)

case=json.loads((HERE/"case.json").read_text(encoding="utf-8"))
pins=json.loads((HERE/"pins.json").read_text(encoding="utf-8"))
saved_path=HERE/"check_saved_result.py"
saved_text=saved_path.read_text(encoding="utf-8")
tree=ast.parse(saved_text,filename=str(saved_path))
need(sha(saved_path)==pins["saved_check_sha256"] and sha(HERE/"case.json")==pins["case_sha256"] and sha(HERE/"check_source_contract.py")==pins["source_check_sha256"],"own source pins mismatch")
need(len(saved_text.splitlines())<=280,"saved checker exceeds 280 source lines")
need(pins["execution_authorized"] is False and pins["root_go"]["authorized"] is False and pins["root_go"]["token"] is None,"source phase must remain disabled")
need(pins["result"]["sha256"] and pins["receipt"]["sha256"],"worker artifacts are not frozen")
need(case["classification"]=="DIAGNOSTIC_ONLY" and case["status"]=="SOURCE_READY_AWAITING_INDEPENDENT_ROOT_GO","case scope/status")
need(case["geometry"]["qp_count"]==16 and case["geometry"]["one_side_area_m2"]==case["geometry"]["R_m"]*case["geometry"]["L_m"],"single-side area contract")
need(case["material_fixture"]["qualification"]=="toy fixture; not ceramic calibration","fixture qualification")
need(case["expected_interface"]["analytic_KG_opening_N_per_m"]==20000.0 and case["geometry"]["opening_index_zero_based"]==82,"analytic interface contract")
need(case["tolerances"]["absolute_energy_J"]<=1e-24 and case["tolerances"]["absolute_energy_density_J_per_m2"]<=1e-20 and case["tolerances"]["kkt_backward_error"]==1e-9,"energy/KKT tolerances must remain discriminating")
need(len(case["required_gates"])==15,"expected producer gate inventory")
need(len(case["source_pins"])>=18,"frozen source inventory incomplete")
for name,ref in case["source_pins"].items():
    path=ROOT/ref["path"]
    need(path.is_file(),f"missing source pin {name}")
    need(sha(path)==ref["sha256"],f"source hash mismatch {name}")
worker_pins_path=ROOT/pins["author_numeric_source_pins"]["path"]
need(worker_pins_path.is_file() and sha(worker_pins_path)==pins["author_numeric_source_pins"]["snapshot_sha256_at_review"],"author pins review snapshot hash")
worker_pins=json.loads(worker_pins_path.read_text(encoding="utf-8"))
need(case["author_numeric_source_pins"]["snapshot_sha256_at_review"]==pins["author_numeric_source_pins"]["snapshot_sha256_at_review"],"author pins snapshot cross-reference")
need(worker_pins["execution_authorized"] is True and worker_pins["root_go"]["authorized"] is True and worker_pins["root_go"]["token"]==pins["author_numeric_source_pins"]["producer_root_go_token"] and pins["author_numeric_source_pins"]["was_root_authorized_for_producer_run"] is True,"author run authorization record")
need(worker_pins["references"]==case["source_pins"],"producer source references differ from frozen independent case")
need(worker_pins["accepted_cache_binding"]["arrays_sha256"]==case["expected_hashes"]["producer_npz"],"accepted NPZ hash metadata")
need(worker_pins["accepted_cache_binding"]["audit_sha256"]==case["expected_hashes"]["recovery_audit"],"recovery audit hash metadata")
need(worker_pins["accepted_cache_binding"]["receipt_sha256"]==case["expected_hashes"]["recovery_receipt"],"recovery receipt hash metadata")
need(worker_pins["accepted_cache_binding"]["actual_array_keys"]==case["cache_array_keys"],"accepted NPZ key metadata")
npz=ROOT/case["producer_npz"]
need(npz.is_file(),"accepted NPZ path missing")
result_path=ROOT/pins["result"]["path"]; receipt_path=ROOT/pins["receipt"]["path"]
need(sha(result_path)==pins["result"]["sha256"] and sha(receipt_path)==pins["receipt"]["sha256"],"frozen result/receipt hashes")
result=json.loads(result_path.read_text(encoding="utf-8")); receipt=json.loads(receipt_path.read_text(encoding="utf-8"))
need(result.get("schema")=="mixed93_actual_cohesive_bridge_numeric_result.v1" and result.get("status")=="ACCEPTED_DIAGNOSTIC_ONLY" and result.get("classification")=="FIXED_MODE_I_ELASTIC_COHESIVE_BRIDGE_ONLY_NOT_MATERIAL_OR_GROWTH_QUALIFICATION","saved worker result contract")
need(receipt.get("outer_return_code")==0 and receipt.get("timed_out") is False and receipt.get("source_pins_equal") is True and receipt.get("budget_passed") is True,"saved worker receipt contract")
need(receipt.get("files",{}).get("result.json",{}).get("sha256")==pins["result"]["sha256"] and receipt.get("files",{}).get("result.json",{}).get("bytes")==result_path.stat().st_size,"receipt result byte binding")
need(receipt.get("source_pre_pins",{}).get("references")==receipt.get("source_post_pins",{}).get("references"),"worker source pre/post pins differ")
for n in ast.walk(tree):
    if isinstance(n,ast.Import): need(not any(a.name.split(".")[0] in {"scipy","importlib","subprocess"} for a in n.names),"forbidden module import")
    if isinstance(n,ast.ImportFrom): need((n.module or "").split(".")[0] not in {"scipy","importlib","subprocess"},"forbidden module import")
    if isinstance(n,ast.Call) and isinstance(n.func,ast.Attribute):
        need(n.func.attr not in {"solve","eig","eigh","eigvals","eigvalsh","load_registered"},"solver/eigensolver/dynamic-loader call found")
numpy_imports=[n for n in ast.walk(tree) if isinstance(n,ast.Import) and any(a.name=="numpy" for a in n.names)]
run_fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=="run")
parents={child:parent for parent in ast.walk(tree) for child in ast.iter_child_nodes(parent)}
need(len(numpy_imports)==1 and parents.get(numpy_imports[0]) is run_fn,"NumPy import must be inside the post-GO audit path")
need("npz_opened\":False" in saved_text and "result_read\":False" in saved_text,"default path must remain read-disabled")
for name in ("result.json","receipt.json"):
    need(name in (pins["result"]["path"] if name=="result.json" else pins["receipt"]["path"]),f"expected frozen {name} path")
print(json.dumps({"status":"SOURCE_READY","npz_opened":False,"npz_hashed":False,"numeric_result_math_checked":False,"saved_result_schema_read":True,"saved_receipt_schema_read":True,"saved_checker_lines":len(saved_text.splitlines()),"source_pin_count":len(case["source_pins"]),"independent_authorization_enabled":False},sort_keys=True))
