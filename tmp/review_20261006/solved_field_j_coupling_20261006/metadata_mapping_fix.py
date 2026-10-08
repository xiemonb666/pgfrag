"""Source-only mapper for the frozen solve's actual quadrature metadata.

This helper reads and structurally validates the frozen solve JSON. It does
not import or invoke the provider, J integrator, or static solver, and writes
no artifacts. It is a candidate for a future adapter patch; it does not change
the frozen adapter or its execution receipt.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[4]
FROZEN_RESULT = (
    PROJECT_ROOT
    / "pgfrag"
    / "tmp"
    / "review_20261006"
    / "modeiii_pu_static_solve_v3_20261006"
    / "runs"
    / "run01"
    / "static_solve_diagnostic.json"
)
FROZEN_RESULT_SHA256 = "39cb2dabe767998baf5e84b038f39510d6ef7db10bcd6c77298b480045cb24d3"
FROZEN_RESULT_SCHEMA = "modeiii_pu_antiplane_static_solve_result.v3"
FROZEN_CASE_ID = "modeiii_pu_uz_static_traction_solve_20261006_v3"
FROZEN_RESULT_STATUS = "SOLVED_DIAGNOSTIC_ONLY_UNQUALIFIED"
EXPECTED_LEVELS = (("volume_gl4x8x4", 8, 128), ("volume_gl4x16x4", 16, 256))


class MetadataMappingError(ValueError):
    """Raised when the frozen source does not match the reviewed structure."""


def map_volume_level(level: dict[str, Any]) -> dict[str, Any]:
    """Derive the descriptor from source fields theta_order and volume_point_count."""
    if type(level) is not dict or level.get("status") != "OK":
        raise MetadataMappingError("volume_level_not_ok_object")
    solution = level.get("solution")
    if type(solution) is not dict or solution.get("status") != "OK":
        raise MetadataMappingError("volume_level_solution_not_ok")
    theta = level.get("theta_order")
    count = level.get("volume_point_count")
    if type(theta) is not int or type(count) is not int:
        raise MetadataMappingError("theta_order_and_volume_point_count_must_be_integers")
    if (level.get("name"), theta, count) not in EXPECTED_LEVELS or count != 4 * theta * 4:
        raise MetadataMappingError("volume_level_metadata_does_not_match_frozen_gl4_theta_gl4_rule")
    return {
        "name": level["name"],
        "quadrature_rule": f"GL4x{theta}x4",
        "theta_order": theta,
        "volume_point_count": count,
        "descriptor_source_fields": ["theta_order", "volume_point_count"],
    }


def map_frozen_result(path: Path = FROZEN_RESULT) -> list[dict[str, Any]]:
    """Read the pinned result and map only its two reviewed volume levels."""
    payload = path.read_bytes()
    if hashlib.sha256(payload).hexdigest() != FROZEN_RESULT_SHA256:
        raise MetadataMappingError("frozen_solve_result_hash_mismatch")
    try:
        result = json.loads(
            payload.decode("utf-8"),
            parse_constant=lambda value: (_ for _ in ()).throw(ValueError(f"nonfinite:{value}")),
        )
    except (UnicodeDecodeError, ValueError) as exc:
        raise MetadataMappingError(f"frozen_solve_json_invalid:{exc}") from exc
    if (
        type(result) is not dict
        or result.get("schema") != FROZEN_RESULT_SCHEMA
        or result.get("case_id") != FROZEN_CASE_ID
        or result.get("status") != FROZEN_RESULT_STATUS
    ):
        raise MetadataMappingError("frozen_solve_identity_mismatch")
    levels = result.get("volume_levels")
    if type(levels) is not list or len(levels) != len(EXPECTED_LEVELS):
        raise MetadataMappingError("frozen_volume_level_count_mismatch")
    by_name = {level.get("name"): level for level in levels if type(level) is dict}
    if len(by_name) != len(EXPECTED_LEVELS) or set(by_name) != {item[0] for item in EXPECTED_LEVELS}:
        raise MetadataMappingError("frozen_volume_level_names_mismatch")
    return [map_volume_level(by_name[name]) for name, _, _ in EXPECTED_LEVELS]
