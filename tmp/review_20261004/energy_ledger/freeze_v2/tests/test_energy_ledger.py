from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path
import unittest

from pgfrag.energy_ledger import (
    BoundaryWorkDefinition,
    Clock,
    EnergyChannel,
    EnergyHistory,
    EnergyTolerance,
    EvidenceRef,
    Qualification,
    audit_energy_history,
)


REPO = Path(__file__).resolve().parents[1]
INPUT_PATH = REPO / "tmp" / "review_20261004" / "energy_ledger" / "input" / "audit_cases.json"
INPUTS = json.loads(INPUT_PATH.read_text(encoding="utf-8"))
TOLERANCE_DATA = INPUTS["tolerances"]


def tolerance(*, constitutive: bool = False) -> EnergyTolerance:
    return EnergyTolerance(
        negative_dissipation_abs_j=TOLERANCE_DATA["negative_dissipation_abs_j"],
        negative_dissipation_definition=TOLERANCE_DATA["negative_dissipation_definition"],
        negative_dissipation_source=TOLERANCE_DATA["negative_dissipation_source"],
        balance_abs_j=TOLERANCE_DATA["balance_abs_j"],
        balance_rel=TOLERANCE_DATA["balance_rel"],
        balance_definition=TOLERANCE_DATA["balance_definition"],
        balance_source=TOLERANCE_DATA["balance_source"],
        constitutive_relative_max=(TOLERANCE_DATA["constitutive_relative_max"] if constitutive else None),
        constitutive_tolerance_definition=(TOLERANCE_DATA["constitutive_tolerance_definition"] if constitutive else None),
        constitutive_tolerance_source=(TOLERANCE_DATA["constitutive_tolerance_source"] if constitutive else None),
    )


def evidence(name: str, method: str | None = None) -> EvidenceRef:
    return EvidenceRef(
        source=f"fixture:{name}",
        method=method or f"explicit {name} fixture construction and source-key review",
    )


def channel(
    values,
    *,
    role: str,
    name: str,
    source_key: str | None = None,
    unit: str = "J",
    clock_id: str | None = "case-clock",
    cumulative: bool = False,
    definition: str | None = None,
    source: str | None = None,
) -> EnergyChannel:
    return EnergyChannel(
        values=values,
        unit=unit,
        definition=definition or f"Signed {role} energy history for {name}.",
        source=source or f"fixture:{name}:raw series",
        source_key=source_key or name,
        clock_id=clock_id,
        role=role,
        cumulative=cumulative,
    )


def qualify(
    regime: str,
    required: list[str],
    *,
    verified: list[str] | None = None,
    clock_id: str = "case-clock",
    include_work: bool = True,
    evidence_method: str | None = None,
    independent_basis: EvidenceRef | None = None,
    target_j: float | None = None,
    target_definition: str | None = None,
    target_source: str | None = None,
) -> Qualification:
    verified = required if verified is None else verified
    return Qualification(
        regime=regime,
        source=f"fixture:{regime}:scope",
        method=evidence_method or f"explicit test fixture classification as {regime}",
        required_balance_channels=required,
        verified_channels={name: evidence(name) for name in verified},
        independent_balance_basis=(
            independent_basis
            if independent_basis is not None
            else evidence("independent balance basis", "closed-form analytic contributions assembled separately from the energy channels")
        ),
        boundary_work=(
            BoundaryWorkDefinition(
                boundary="synthetic driven boundary",
                force_or_reaction="independently supplied analytic boundary force history",
                displacement_or_velocity="independently supplied prescribed displacement history",
                integration_rule="cumulative exact values in analytic positive-control fixture",
                sign_convention="positive_into_system",
                source="fixture:analytic-boundary-work",
                clock_id=clock_id,
            )
            if include_work
            else None
        ),
        constitutive_target_j=target_j,
        constitutive_target_definition=target_definition,
        constitutive_target_source=target_source,
    )


def source_endpoint_history(values: list[float], source_key: str = "dissipated") -> EnergyHistory:
    return EnergyHistory(
        clock=Clock(
            values=[0.0, 1.0],
            unit="s",
            definition="Two endpoint samples from an existing loaded run; not represented as a complete trace.",
            source="relicweave-synlab/docs/a12/COHESIVE_CAPABILITY_20261004.md",
            identity="case-clock",
        ),
        channels={
            "cohesive_dissipated": channel(
                values,
                role="dissipated",
                name="cohesive_dissipated",
                source_key=source_key,
                cumulative=True,
                source="relicweave-synlab/docs/a12/COHESIVE_CAPABILITY_20261004.md:loaded-case native signed endpoint",
            )
        },
        tolerance=tolerance(),
    )


def analytic_free_history(
    series: dict[str, tuple[list[float], str, str, bool]],
    *,
    clock_values: list[float] | None = None,
    balance_rel: float | None = None,
) -> tuple[EnergyHistory, Qualification]:
    """Build an explicitly synthetic, independently specified ledger fixture."""

    clock_id = "synthetic-analytic-clock"
    channels = {
        name: channel(
            values,
            role=role,
            name=name,
            unit=unit,
            cumulative=cumulative,
            clock_id=clock_id,
            definition=f"Synthetic analytic {role} series for the {name} test channel.",
        )
        for name, (values, role, unit, cumulative) in series.items()
    }
    selected_tolerance = tolerance()
    if balance_rel is not None:
        selected_tolerance = replace(selected_tolerance, balance_rel=balance_rel)
    history = EnergyHistory(
        clock=Clock(
            values=[0.0, 1.0] if clock_values is None else clock_values,
            unit="s",
            definition="Synthetic analytic test clock.",
            source="fixture:synthetic-analytic-clock",
            identity=clock_id,
        ),
        channels=channels,
        tolerance=selected_tolerance,
    )
    qualification = qualify("free_dynamics", list(series), clock_id=clock_id)
    return history, qualification


class EnergyLedgerTests(unittest.TestCase):
    def test_both_retained_extrinsic_negative_dissipation_cases_fail_signed(self) -> None:
        cases = (
            INPUTS["source_grounded_values"]["extrinsic_dt_dissipated_endpoint_J"],
            INPUTS["source_grounded_values"]["extrinsic_dt2_dissipated_endpoint_J"],
        )
        for values in cases:
            with self.subTest(raw_final=values[-1]):
                result = audit_energy_history(
                    source_endpoint_history(values),
                    qualify("prescribed_boundary_dynamics", ["cohesive_dissipated"]),
                )
                self.assertEqual(result["dissipation"]["status"], "FAIL")
                self.assertEqual(result["overall_status"], "FAIL")
                self.assertEqual(result["channels"]["cohesive_dissipated"]["values_raw"], values)
                self.assertEqual(
                    result["channels"]["cohesive_dissipated"]["values_raw"][-1],
                    values[-1],
                )
                self.assertLess(result["dissipation"]["channels"][0]["minimum_J"], -1e-12)
                self.assertFalse(result["dissipation"]["raw_values_are_clamped_or_abs_transformed"])
                self.assertEqual(result["free_energy_balance"]["status"], "UNKNOWN")

    def test_tiny_negative_roundoff_is_preserved_and_not_failed(self) -> None:
        raw = INPUTS["source_grounded_values"]["run03_roundoff_minimum_J"]
        history = EnergyHistory(
            clock=Clock(
                values=[4.0],
                unit="sample_index",
                definition="Native sample index for the retained roundoff-minimum observation.",
                source=INPUTS["evidence_sources"]["intrinsic_run03_summary"],
                identity="run03-point-index",
            ),
            channels={
                "cohesive_dissipated": channel(
                    [raw],
                    role="dissipated",
                    name="cohesive_dissipated",
                    source_key="dissipated",
                    clock_id="run03-point-index",
                    cumulative=True,
                    source=INPUTS["evidence_sources"]["intrinsic_run03_summary"],
                )
            },
            tolerance=tolerance(),
        )
        result = audit_energy_history(
            history,
            qualify("all_dof_prescribed_constitutive", ["cohesive_dissipated"], verified=["cohesive_dissipated"]),
        )
        self.assertEqual(result["dissipation"]["status"], "PASS_WITH_ROUNDOFF")
        self.assertEqual(result["channels"]["cohesive_dissipated"]["values_raw"], [raw])
        self.assertEqual(result["dissipation"]["channels"][0]["minimum_raw_value"], raw)
        self.assertEqual(result["free_energy_balance"]["status"], "NOT_APPLICABLE")

    def test_missing_contact_and_clock_keep_free_balance_unknown(self) -> None:
        clock = Clock(
            values=[0.0, 1.0],
            unit="s",
            definition=None,
            source=None,
            identity=None,
        )
        channels = {
            "kinetic": channel([0.0, 0.1], role="kinetic", name="kinetic", clock_id=None),
            "bulk_potential": channel([1.0, 0.9], role="stored", name="bulk_potential", clock_id=None),
            "dissipated": channel([0.0, 0.01], role="dissipated", name="dissipated", clock_id=None, cumulative=True),
            "external_work": channel([0.0, 0.11], role="external_work", name="external_work", clock_id=None, cumulative=True),
            # This is an explicit legacy alias diagnostic, not the required native cohesive-contact channel.
            "contact_alias_diagnostic": channel(
                [0.0, 0.0], role="diagnostic", name="contact_alias_diagnostic",
                source_key="contact", clock_id=None,
            ),
        }
        history = EnergyHistory(clock=clock, channels=channels, tolerance=tolerance())
        required = ["kinetic", "bulk_potential", "dissipated", "cohesive_contact", "external_work"]
        result = audit_energy_history(
            history,
            qualify("free_dynamics", required, include_work=False),
        )
        self.assertEqual(result["free_energy_balance"]["status"], "UNKNOWN")
        self.assertIn("cohesive_contact", result["missing_required_channels"])
        self.assertEqual(result["channels"]["cohesive_contact"]["status"], "UNKNOWN")
        self.assertIsNone(result["channels"]["cohesive_contact"]["values_raw"])
        self.assertEqual(result["channels"]["contact_alias_diagnostic"]["values_raw"], [0.0, 0.0])
        self.assertEqual(result["channels"]["contact_alias_diagnostic"]["role"], "diagnostic")
        self.assertEqual(result["clock"]["status"], "UNKNOWN")
        self.assertEqual(result["material_calibration"]["status"], "UNKNOWN")

    def test_incompatible_units_are_rejected(self) -> None:
        history = EnergyHistory(
            clock=Clock([0.0], "s", "single observation time", "fixture", "clock"),
            channels={"kinetic": channel([2.0], role="kinetic", name="kinetic", unit="N", clock_id="clock")},
            tolerance=tolerance(),
        )
        with self.assertRaisesRegex(ValueError, "not a supported energy unit"):
            audit_energy_history(history, qualify("prescribed_boundary_dynamics", ["kinetic"]))

    def test_intrinsic_201_point_toy_is_constitutive_only_and_matches_gc_area(self) -> None:
        values = INPUTS["source_grounded_values"]
        final_d = values["run03_dissipated_initial_and_final_J"][1]
        target = values["run03_gc_times_area_J"]
        history = EnergyHistory(
            clock=Clock(
                values=[200.0],
                unit="sample_index",
                definition="Endpoint at native point index 200 from a complete 201-point run.",
                source=INPUTS["evidence_sources"]["intrinsic_run03_summary"],
                identity="run03-point-index",
            ),
            channels={
                "cohesive_dissipated": channel(
                    [final_d],
                    role="dissipated",
                    name="cohesive_dissipated",
                    source_key="dissipated",
                    clock_id="run03-point-index",
                    cumulative=True,
                    source=INPUTS["evidence_sources"]["intrinsic_run03_summary"],
                )
            },
            tolerance=tolerance(constitutive=True),
        )
        qualification = qualify(
            "all_dof_prescribed_constitutive",
            ["cohesive_dissipated"],
            verified=["cohesive_dissipated"],
            target_j=target,
            target_definition="Gc times the declared interface area.",
            target_source=INPUTS["evidence_sources"]["intrinsic_run03_summary"],
        )
        result = audit_energy_history(history, qualification)
        self.assertEqual(result["free_energy_balance"]["status"], "NOT_APPLICABLE")
        self.assertEqual(result["constitutive_reference_comparison"]["status"], "MATCHED_REFERENCE")
        self.assertAlmostEqual(
            result["constitutive_reference_comparison"]["relative_error"],
            values["run03_reported_relative_error"],
            delta=1e-14,
        )
        self.assertEqual(result["constitutive_reference_comparison"]["final_raw_value"], final_d)
        self.assertEqual(result["material_calibration"]["status"], "UNKNOWN")
        self.assertEqual(result["overall_status"], "CONSTITUTIVE_ONLY")

    def test_one_sample_free_dynamics_has_no_evaluable_interval(self) -> None:
        one_sample_clock = INPUTS["v2_numeric_regressions"]["single_sample_clock_s"]
        history, qualification = analytic_free_history(
            {
                "kinetic": ([0.0], "kinetic", "J", False),
                "bulk_potential": ([1.0], "stored", "J", False),
                "dissipated": ([0.0], "dissipated", "J", True),
                "external_work": ([0.0], "external_work", "J", True),
            },
            clock_values=one_sample_clock,
        )
        result = audit_energy_history(history, qualification)
        self.assertEqual(
            result["free_energy_balance"]["status"],
            "NOT_EVALUATED_INSUFFICIENT_HISTORY",
        )
        self.assertEqual(result["free_energy_balance"]["sample_count"], 1)
        self.assertEqual(result["overall_status"], "UNKNOWN")

    def test_finite_overflow_in_balance_arithmetic_cannot_pass(self) -> None:
        overflow_input = INPUTS["v2_numeric_regressions"]["finite_balance_overflow"]
        controls = {
            # The specified regression: the residual's storage-plus-D sum overflows.
            "residual_overflow": {
                "kinetic": (overflow_input["kinetic_J"], "kinetic", "J", False),
                "bulk_potential": (overflow_input["stored_J"], "stored", "J", False),
                "dissipated": (overflow_input["dissipated_J"], "dissipated", "J", True),
                "external_work": (overflow_input["external_work_J"], "external_work", "J", True),
            },
            # The residual stays finite, but the absolute scale sum overflows.
            "scale_overflow": {
                "kinetic": ([0.0, 1.0e308], "kinetic", "J", False),
                "bulk_potential": ([0.0, -1.0e308], "stored", "J", False),
                "dissipated": ([0.0, 1.0e308], "dissipated", "J", True),
                "external_work": ([0.0, 1.0e308], "external_work", "J", True),
            },
        }
        for name, series in controls.items():
            with self.subTest(control=name):
                history, qualification = analytic_free_history(
                    series,
                    clock_values=overflow_input["clock_s"] if name == "residual_overflow" else None,
                )
                result = audit_energy_history(history, qualification)
                self.assertEqual(result["free_energy_balance"]["status"], "UNKNOWN_NUMERIC_RANGE")
                self.assertEqual(result["overall_status"], "UNKNOWN_NUMERIC_RANGE")
                self.assertTrue(result["free_energy_balance"]["raw_values_are_preserved"])
                for channel_name, (values, _role, _unit, _cumulative) in series.items():
                    self.assertEqual(result["channels"][channel_name]["values_raw"], values)

    def test_aggregation_baseline_difference_and_cumulative_increment_overflow_are_unknown(self) -> None:
        controls = {
            "aggregation": {
                "kinetic_a": ([0.0, 1.0e308], "kinetic", "J", False),
                "kinetic_b": ([0.0, 1.0e308], "kinetic", "J", False),
                "bulk_potential": ([0.0, 0.0], "stored", "J", False),
                "dissipated": ([0.0, 0.0], "dissipated", "J", True),
                "external_work": ([0.0, 0.0], "external_work", "J", True),
            },
            "baseline_difference": {
                "kinetic": ([0.0, 0.0], "kinetic", "J", False),
                "bulk_potential": ([1.0e308, -1.0e308], "stored", "J", False),
                "dissipated": ([0.0, 0.0], "dissipated", "J", True),
                "external_work": ([0.0, 0.0], "external_work", "J", True),
            },
            "dissipation_increment": {
                "kinetic": ([0.0, 0.0], "kinetic", "J", False),
                "bulk_potential": ([0.0, 0.0], "stored", "J", False),
                "dissipated": ([1.0e308, -1.0e308], "dissipated", "J", True),
                "external_work": ([0.0, 0.0], "external_work", "J", True),
            },
        }
        for name, series in controls.items():
            with self.subTest(control=name):
                history, qualification = analytic_free_history(series)
                result = audit_energy_history(history, qualification)
                if name == "dissipation_increment":
                    self.assertEqual(result["dissipation"]["status"], "FAIL")
                    self.assertEqual(
                        result["dissipation"]["numeric_range_status"],
                        "UNKNOWN_NUMERIC_RANGE",
                    )
                    self.assertEqual(result["overall_status"], "FAIL")
                else:
                    self.assertEqual(result["overall_status"], "UNKNOWN_NUMERIC_RANGE")
                self.assertIn(
                    result["free_energy_balance"]["status"],
                    {"UNKNOWN_NUMERIC_RANGE", "UNKNOWN"},
                )

    def test_overflowing_relative_limit_and_constitutive_ratio_are_unknown(self) -> None:
        history, qualification = analytic_free_history(
            {
                "kinetic": ([0.0, 0.0], "kinetic", "J", False),
                "bulk_potential": ([0.0, 2.0], "stored", "J", False),
                "dissipated": ([0.0, 1.0], "dissipated", "J", True),
                "external_work": ([0.0, 0.0], "external_work", "J", True),
            },
            balance_rel=1.0e308,
        )
        result = audit_energy_history(history, qualification)
        self.assertEqual(result["free_energy_balance"]["status"], "UNKNOWN_NUMERIC_RANGE")
        self.assertEqual(result["overall_status"], "UNKNOWN_NUMERIC_RANGE")

        all_dof_history = EnergyHistory(
            clock=Clock([0.0], "sample_index", "one constitutive endpoint", "fixture", "c"),
            channels={
                "dissipated": channel([1.0e308], role="dissipated", name="dissipated", clock_id="c", cumulative=True)
            },
            tolerance=replace(tolerance(constitutive=True), constitutive_relative_max=1.0),
        )
        all_dof_qualification = qualify(
            "all_dof_prescribed_constitutive",
            ["dissipated"],
            verified=["dissipated"],
            target_j=1.0e-308,
            target_definition="synthetic tiny positive comparison target",
            target_source="fixture:constitutive-overflow-control",
        )
        constitutive = audit_energy_history(all_dof_history, all_dof_qualification)
        self.assertEqual(
            constitutive["constitutive_reference_comparison"]["status"],
            "UNKNOWN_NUMERIC_RANGE",
        )
        self.assertEqual(constitutive["overall_status"], "UNKNOWN_NUMERIC_RANGE")

    def test_finite_inputs_that_overflow_unit_conversion_are_rejected(self) -> None:
        history, qualification = analytic_free_history(
            {
                "kinetic": ([0.0, 0.0], "kinetic", "J", False),
                "bulk_potential": (
                    INPUTS["v2_numeric_regressions"]["finite_stored_kwh_conversion_overflow_raw"],
                    "stored",
                    "kWh",
                    False,
                ),
                "dissipated": ([0.0, 0.0], "dissipated", "J", True),
                "external_work": ([0.0, 0.0], "external_work", "J", True),
            }
        )
        with self.assertRaisesRegex(ValueError, "overflows while converting"):
            audit_energy_history(history, qualification)

    def test_bool_numeric_inputs_and_non_bool_cumulative_flags_are_rejected(self) -> None:
        valid_series = {
            "kinetic": ([0.0, 0.0], "kinetic", "J", False),
            "bulk_potential": ([0.0, 0.0], "stored", "J", False),
            "dissipated": ([0.0, 0.0], "dissipated", "J", True),
            "external_work": ([0.0, 0.0], "external_work", "J", True),
        }
        history, qualification = analytic_free_history(valid_series)
        bool_energy_series = dict(valid_series)
        bool_energy_series["kinetic"] = ([0.0, True], "kinetic", "J", False)
        bool_energy_history, bool_energy_qualification = analytic_free_history(bool_energy_series)
        with self.assertRaisesRegex(ValueError, "non-boolean real number"):
            audit_energy_history(bool_energy_history, bool_energy_qualification)

        bool_clock_history = replace(
            history,
            clock=Clock([False, True], "s", "bool clock control", "fixture", "synthetic-analytic-clock"),
        )
        with self.assertRaisesRegex(ValueError, "non-boolean real number"):
            audit_energy_history(bool_clock_history, qualification)

        for field_name in (
            "negative_dissipation_abs_j",
            "balance_abs_j",
            "balance_rel",
        ):
            with self.subTest(tolerance_field=field_name):
                bool_tolerance_history = replace(
                    history,
                    tolerance=replace(history.tolerance, **{field_name: True}),
                )
                with self.assertRaisesRegex(ValueError, "non-boolean real number"):
                    audit_energy_history(bool_tolerance_history, qualification)
        bool_constitutive_tolerance = replace(
            tolerance(constitutive=True),
            constitutive_relative_max=True,
        )
        with self.assertRaisesRegex(ValueError, "non-boolean real number"):
            audit_energy_history(
                replace(history, tolerance=bool_constitutive_tolerance),
                qualification,
            )

        bool_target = replace(
            qualify(
                "all_dof_prescribed_constitutive",
                ["dissipated"],
                target_j=1.0,
                target_definition="fixture target",
                target_source="fixture target source",
            ),
            constitutive_target_j=True,
        )
        with self.assertRaisesRegex(ValueError, "non-boolean real number"):
            audit_energy_history(history, bool_target)

        string_cumulative_series = dict(valid_series)
        string_cumulative_series["dissipated"] = ([0.0, 0.0], "dissipated", "J", "false")
        string_cumulative_history, string_cumulative_qualification = analytic_free_history(string_cumulative_series)
        with self.assertRaisesRegex(ValueError, "cumulative must be a bool"):
            audit_energy_history(string_cumulative_history, string_cumulative_qualification)

    def test_required_balance_channels_must_be_nonempty_unique_legal_keys(self) -> None:
        history, qualification = analytic_free_history(
            {
                "kinetic": ([0.0, 0.0], "kinetic", "J", False),
                "bulk_potential": ([0.0, 0.0], "stored", "J", False),
                "dissipated": ([0.0, 0.0], "dissipated", "J", True),
                "external_work": ([0.0, 0.0], "external_work", "J", True),
            }
        )
        for required in ([], ["dissipated", "dissipated"], ["bad key"]):
            with self.subTest(required=required):
                invalid = replace(qualification, required_balance_channels=required)
                with self.assertRaises(ValueError):
                    audit_energy_history(history, invalid)

    def test_analytic_independent_free_dynamics_positive_control_closes(self) -> None:
        data = INPUTS["positive_control"]
        names_and_roles = {
            "kinetic": ("kinetic", data["kinetic_J"], False),
            "bulk_potential": ("stored", data["bulk_stored_J"], False),
            "cohesive_reversible": ("stored", data["cohesive_reversible_J"], False),
            "cohesive_contact": ("stored", data["cohesive_contact_J"], False),
            "dissipated": ("dissipated", data["dissipated_J"], True),
            "external_work": ("external_work", data["external_work_J"], True),
        }
        channels = {
            name: channel(
                values,
                role=role,
                name=name,
                source_key=("cohesive contact" if name == "cohesive_contact" else name),
                cumulative=cumulative,
                clock_id="analytic-clock",
                definition=f"Closed-form analytical {role} energy used only for the positive control: {name}.",
            )
            for name, (role, values, cumulative) in names_and_roles.items()
        }
        history = EnergyHistory(
            clock=Clock(
                values=data["clock_s"],
                unit="s",
                definition="Exact analytic positive-control sample clock.",
                source="fixture:analytic-clock",
                identity="analytic-clock",
            ),
            channels=channels,
            tolerance=tolerance(),
        )
        required = list(names_and_roles)
        qualification = qualify("free_dynamics", required, clock_id="analytic-clock")
        result = audit_energy_history(history, qualification)
        self.assertEqual(result["free_energy_balance"]["status"], "PASS")
        self.assertEqual(result["free_energy_balance"]["residuals_J"], data["expected_residuals_J"])
        self.assertEqual(result["free_energy_balance"]["baseline"]["sample_index"], 0)
        self.assertEqual(result["free_energy_balance"]["sample_count"], 3)
        self.assertEqual(result["dissipation"]["status"], "PASS")
        self.assertEqual(result["overall_status"], "PASS_ENERGY_AUDIT_ONLY")
        self.assertEqual(result["material_calibration"]["status"], "UNKNOWN")

    def test_status_text_alone_is_not_channel_evidence(self) -> None:
        history = source_endpoint_history([0.0, 1e-6])
        qualification = qualify(
            "free_dynamics",
            ["cohesive_dissipated"],
            verified=[],
            evidence_method="PASS",
            independent_basis=EvidenceRef(source="fixture", method="PASS"),
        )
        result = audit_energy_history(history, qualification)
        self.assertEqual(result["channels"]["cohesive_dissipated"]["semantic_evidence_status"], "UNKNOWN")
        self.assertEqual(result["free_energy_balance"]["status"], "UNKNOWN")


if __name__ == "__main__":
    unittest.main()
