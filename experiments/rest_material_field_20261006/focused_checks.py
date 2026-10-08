"""Small, deterministic contract checks for the rest-frame material field."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any

import numpy as np

from rest_material_field import (
    MaterialFieldParameters,
    MaterialFieldValidationError,
    RestFrameMaterialField,
    RigidFrameTransform,
)


@dataclass
class QueryLedger:
    """Counts requested point evaluations, including rejected coordinate probes."""

    count: int = 0


class CountedField:
    def __init__(self, field: RestFrameMaterialField, ledger: QueryLedger) -> None:
        self.field = field
        self.ledger = ledger

    def evaluate(self, point_m: Any) -> float:
        self.ledger.count += 1
        return self.field.evaluate(point_m)

    def evaluate_many(self, points_m: Any) -> tuple[float, ...]:
        # Count one attempted evaluation per row, even if field validation fails.
        try:
            query_count = len(points_m)
        except TypeError:
            query_count = 1
        self.ledger.count += query_count
        return self.field.evaluate_many(points_m)

    def evaluate_in_frame(self, point_frame_m: Any, transform: RigidFrameTransform) -> float:
        self.ledger.count += 1
        return self.field.evaluate_in_frame(point_frame_m, transform)

    def evaluate_many_in_frame(self, points_frame_m: Any, transform: RigidFrameTransform) -> tuple[float, ...]:
        try:
            query_count = len(points_frame_m)
        except TypeError:
            query_count = 1
        self.ledger.count += query_count
        return self.field.evaluate_many_in_frame(points_frame_m, transform)

    def evaluate_positive_parameter(self, nominal_value: float, point_m: Any) -> float:
        self.ledger.count += 1
        return self.field.evaluate_positive_parameter(nominal_value, point_m)


def _expect_validation_error(callable_: Any, label: str) -> None:
    try:
        callable_()
    except MaterialFieldValidationError:
        return
    raise AssertionError(f"{label}: expected MaterialFieldValidationError")


def run_focused_checks() -> dict[str, Any]:
    """Run the focused invariance/validation checks and return evidence.

    The tests make 31 point-evaluation attempts total.  Review artifact
    generation adds 196 actual sample points; with the recorded 26-query
    failed implementation pass, cumulative use remains below 256.
    """

    ledger = QueryLedger()
    field = RestFrameMaterialField(MaterialFieldParameters(seed=20261006, correlation_length_m=0.012, amplitude=0.30))
    counted = CountedField(field, ledger)
    completed: list[str] = []

    # Reordering node identities or sample order cannot change the value at a
    # physical coordinate.  A repeated query is identical as well.
    coordinates = np.asarray(
        [[-0.012, 0.0, 0.0], [-0.004, 0.003, 0.0], [0.002, -0.005, 0.001], [0.011, 0.002, -0.003]],
        dtype=np.float64,
    )
    order = np.asarray([2, 0, 3, 1])
    values = counted.evaluate_many(coordinates)
    reordered = counted.evaluate_many(coordinates[order])
    assert np.array_equal(np.asarray(reordered), np.asarray(values)[order])
    repeated = counted.evaluate(coordinates[1])
    assert repeated == values[1]
    completed.append("node_order_and_repeated_point_are_invariant")

    # Two different samplings contain a shared material point.  Only its
    # physical rest coordinate is used, so changed node counts/positions do
    # not resample the material field.
    mesh_a = np.asarray([[-0.01, 0.0, 0.0], [0.0, 0.0, 0.0], [0.01, 0.0, 0.0]])
    mesh_b = np.asarray([[-0.007, 0.0, 0.0], [0.0, 0.0, 0.0], [0.013, 0.0, 0.0], [0.026, 0.0, 0.0]])
    mesh_a_values = counted.evaluate_many(mesh_a)
    mesh_b_values = counted.evaluate_many(mesh_b)
    assert mesh_a_values[1] == mesh_b_values[1]
    completed.append("resampling_at_same_physical_point_is_invariant")

    # Current-frame coordinates are explicitly inverted to the immutable rest
    # frame before field evaluation.
    theta = 0.53
    rotation = np.asarray(
        [[math.cos(theta), -math.sin(theta), 0.0], [math.sin(theta), math.cos(theta), 0.0], [0.0, 0.0, 1.0]],
        dtype=np.float64,
    )
    translation = np.asarray([0.12, -0.04, 0.08])
    transform = RigidFrameTransform(
        rotation_rest_to_frame=rotation,
        translation_rest_to_frame_m=translation,
        frame_id="rigidly-moved-current-frame",
    )
    rest_points = np.asarray([[0.0, 0.0, 0.0], [0.003, -0.002, 0.001], [-0.008, 0.004, 0.002]])
    frame_points = rest_points @ rotation.T + translation
    rest_values = counted.evaluate_many(rest_points)
    frame_values = counted.evaluate_many_in_frame(frame_points, transform)
    assert np.allclose(frame_values, rest_values, rtol=0.0, atol=2e-14)
    completed.append("rigid_frame_transform_preserves_rest_binding")

    # Zero amplitude is exactly homogeneous, including for the positive map.
    homogeneous = RestFrameMaterialField(
        MaterialFieldParameters(seed=3, correlation_length_m=0.012, amplitude=0.0)
    )
    homogeneous_counted = CountedField(homogeneous, ledger)
    homogeneous_values = homogeneous_counted.evaluate_many(np.asarray([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, -2.0, 0.4], [0.2, 0.3, 0.4]]))
    assert homogeneous_values == (0.0, 0.0, 0.0, 0.0)
    assert homogeneous.positive_parameter_bounds(125.0) == (125.0, 125.0)
    assert homogeneous_counted.evaluate_positive_parameter(125.0, [0.0, 0.0, 0.0]) == 125.0
    # This final positive-map point is deliberately evaluated with the main
    # field so the upper and lower mapped bounds are checked on real values.
    nominal = 125.0
    bounds = field.positive_parameter_bounds(nominal)
    for point in ([0.001, 0.002, 0.003], [-0.006, 0.004, 0.002], [0.009, -0.001, 0.005]):
        mapped = counted.evaluate_positive_parameter(nominal, point)
        assert math.isfinite(mapped) and bounds[0] <= mapped <= bounds[1] and mapped > 0.0
    completed.append("zero_amplitude_and_bounded_positive_map")

    _expect_validation_error(lambda: MaterialFieldParameters(seed=1, correlation_length_m=0.0), "zero correlation length")
    _expect_validation_error(lambda: MaterialFieldParameters(seed=1, correlation_length_m=-0.01), "negative correlation length")
    _expect_validation_error(lambda: MaterialFieldParameters(seed=1, correlation_length_m=math.nan), "nonfinite correlation length")
    _expect_validation_error(
        lambda: MaterialFieldParameters(seed=1, correlation_length_m=0.01, coordinate_units="mm"),
        "non-metre coordinate units",
    )
    _expect_validation_error(
        lambda: MaterialFieldParameters(seed=1, correlation_length_m=0.01, correlation_length_units="mm"),
        "non-metre correlation length units",
    )
    _expect_validation_error(
        lambda: MaterialFieldParameters(seed=1, correlation_length_m=0.01, amplitude=math.inf),
        "nonfinite amplitude",
    )
    # Count the rejected query too, so the receipt documents all evaluation
    # attempts rather than only successful points.
    _expect_validation_error(lambda: counted.evaluate([math.nan, 0.0, 0.0]), "nonfinite query coordinate")
    completed.append("invalid_units_lengths_amplitude_and_coordinates_are_rejected")

    assert ledger.count == 31, f"focused check query budget changed: {ledger.count} != 31"
    return {
        "status": "PASS",
        "checks": completed,
        "assertion_groups": len(completed),
        "point_evaluation_attempts": ledger.count,
        "invalid_query_attempts_included": True,
    }


if __name__ == "__main__":
    import json

    print(json.dumps(run_focused_checks(), indent=2))
