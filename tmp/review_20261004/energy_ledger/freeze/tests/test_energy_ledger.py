from __future__ import annotations

import json
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
                    qualify("prescribed_boundary_dynamics", []),
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
            qualify("all_dof_prescribed_constitutive", [], verified=["cohesive_dissipated"]),
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
            audit_energy_history(history, qualify("prescribed_boundary_dynamics", []))

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
            [],
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
