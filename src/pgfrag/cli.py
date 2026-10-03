from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Sequence

from .assets import validate_asset
from .backend import run_job
from .errors import PGFragError
from .jobs import load_job, safe_output_path


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="pgfrag", description="Independent PG-Frag prototype CLI")
    sub = parser.add_subparsers(dest="command", required=True)
    asset = sub.add_parser("validate-asset", help="validate forge.asset.v1 package layout and hashes")
    asset.add_argument("package_dir")
    job = sub.add_parser("validate-job", help="validate a pgfrag.job.v1 JSON job")
    job.add_argument("job_file")
    run = sub.add_parser("run", help="run a bounded CPU solver job or write a structured unsupported result")
    run.add_argument("job_file")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    try:
        if args.command == "validate-asset":
            result = validate_asset(args.package_dir)
            print(json.dumps(result, indent=2, sort_keys=True))
            return 0
        if args.command == "validate-job":
            job, job_path = load_job(args.job_file)
            output = safe_output_path(job_path, job["output_dir"])
            print(json.dumps({"status": "PASS", "schema": job["schema"], "job_id": job["job_id"],
                              "backend_id": job["backend_id"], "estimated_points": job.get("initial_condition", {}).get("estimated_points"),
                              "output_dir": str(output)}, indent=2, sort_keys=True))
            return 0
        if args.command == "run":
            job, job_path = load_job(args.job_file)
            expected_result_path = safe_output_path(job_path, job["output_dir"]) / "result.json"
            code, result = run_job(args.job_file)
            print(json.dumps({"status": result.get("status"), "job_id": result.get("job_id"),
                              "result_path": result.get("result_path", str(expected_result_path.resolve(strict=False)))},
                             indent=2, sort_keys=True, allow_nan=False))
            return code
    except PGFragError as exc:
        print(f"pgfrag: {exc}", file=sys.stderr)
        return 3 if exc.__class__.__name__ == "UnsupportedCapability" else 2
    except (OSError, ValueError) as exc:
        print(f"pgfrag: {exc}", file=sys.stderr)
        return 2
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
