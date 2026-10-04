from __future__ import annotations

import hashlib
import inspect
import json
import sys
import traceback
from pathlib import Path

import akantu as aka

OUT = Path("/mnt/c/Users/xiemo/Desktop/ShardForge/relicweave-synlab/tmp/review_20261004/implementation/cohesive_bilinear_patch_20261004/api_symbol_audit_v1")
TARGETS = (
    "initFull", "assembleMassLumped", "setTimeStep", "addDumpFieldToDumper",
    "getEnergy", "getMass", "dump", "solveStep", "getInternalForce", "getTime",
    "getDisplacement", "getBlockedDOFs", "getBlockedDofs", "displacement", "blocked_dofs",
)


def describe(cls, name):
    present = hasattr(cls, name)
    result = {"name": name, "present_on_class": bool(present)}
    if not present:
        return result
    try:
        member = getattr(cls, name)
        result["descriptor_type"] = type(member).__name__
        result["callable"] = callable(member)
        result["docstring"] = inspect.getdoc(member)
    except Exception as error:
        result["inspection_error"] = repr(error)
    return result


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    target = OUT / "native_model_symbol_audit_v1.json"
    if target.exists():
        raise SystemExit("audit output exists; refusing repeat")
    cls = aka.SolidMechanicsModelCohesive
    names = dir(cls)
    result = {
        "status": "READ_ONLY_SYMBOL_AUDIT_COMPLETE",
        "argv": sys.argv,
        "python": sys.version.split()[0],
        "akantu_module": str(Path(aka.__file__).resolve()),
        "akantu_version": getattr(aka, "__version__", "UNKNOWN"),
        "model_class": f"{cls.__module__}.{cls.__qualname__}",
        "mro": [f"{c.__module__}.{c.__qualname__}" for c in cls.__mro__],
        "target_members": [describe(cls, name) for name in TARGETS],
        "class_dir_blocked_names": [n for n in names if "blocked" in n.lower()],
        "class_dir_displacement_names": [n for n in names if "displacement" in n.lower()],
        "all_model_public_names_used_by_driver": {
            n: describe(cls, n) for n in
            ("initFull", "assembleMassLumped", "setTimeStep", "addDumpFieldToDumper",
             "getEnergy", "getMass", "dump", "solveStep", "getInternalForce", "getTime")
        },
        "model_or_mesh_constructed": False,
        "solver_calls": 0,
        "limitations": "Class-level symbols/docstrings only; no instance view writability or mutation semantics inferred.",
        "source_to_tag_correspondence": "UNKNOWN",
    }
    target.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps(result, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())