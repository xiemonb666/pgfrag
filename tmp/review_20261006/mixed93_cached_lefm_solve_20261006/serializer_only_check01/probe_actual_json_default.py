import ast, json, pathlib, numpy as np
source_path = pathlib.Path(r"C:\Users\xiemo\Desktop\ShardForge\pgfrag\experiments\mixed93_cached_lefm_solve_20261006\mixed93_cached_lefm_solve.py")
tree = ast.parse(source_path.read_text(encoding="utf-8"), filename=str(source_path))
writer = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "_write_result")
default_node = next(node for node in writer.body if isinstance(node, ast.FunctionDef) and node.name == "json_default")
namespace = {"np": np}
exec(compile(ast.Module(body=[default_node], type_ignores=[]), str(source_path), "exec"), namespace)
default = namespace["json_default"]
payload = {"bool": np.bool_(False), "int": np.int64(3), "float": np.float64(0.25)}
encoded = json.dumps(payload, default=default, allow_nan=False, sort_keys=True)
decoded = json.loads(encoded)
assert type(decoded["bool"]) is bool and decoded["bool"] is False
assert type(decoded["int"]) is int and decoded["int"] == 3
assert type(decoded["float"]) is float and decoded["float"] == 0.25
try:
    json.dumps({"nan": np.float64(float("nan"))}, default=default, allow_nan=False)
except ValueError:
    nan_rejected = True
else:
    raise AssertionError("NaN should still fail with allow_nan=False")
try:
    json.dumps({"unsupported": object()}, default=default, allow_nan=False)
except TypeError:
    unsupported_rejected = True
else:
    raise AssertionError("unsupported object should still raise TypeError")
print(json.dumps({"status":"PASS","decoded":decoded,"decoded_types":{k:type(v).__name__ for k,v in decoded.items()},
                  "nan_rejected":nan_rejected,"unsupported_rejected":unsupported_rejected},sort_keys=True))
