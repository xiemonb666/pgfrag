"""Focused standard-library fake-field checks against the frozen contract source."""

from __future__ import annotations

import importlib.util
import os
from pathlib import Path
import sys
import unittest
import xml.etree.ElementTree as ET


REVIEW_DIR = Path(__file__).resolve().parents[2] / "tmp" / "review_20261005" / "enriched_dof_contract_20261005"
SNAPSHOT_PATH = Path(
    os.environ.get("SHARDFORGE_DOF_SNAPSHOT_PATH", REVIEW_DIR / "enriched_dof_contract_frozen.py")
)
RESULT_DIR = Path(os.environ.get("SHARDFORGE_DOF_RESULT_DIR", REVIEW_DIR / "fake_field_run_20261005"))


def _load_snapshot():
    spec = importlib.util.spec_from_file_location("enriched_dof_contract_frozen", SNAPSHOT_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load frozen snapshot at {SNAPSHOT_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


contract = _load_snapshot()


def base_key(source_id: str, component=None):
    if component is None:
        component = contract.FieldComponent.UZ
    return contract.DofKey(source_id, component, contract.BasisKind.BASE, 0)


def front_key(source_id: str, kind=None):
    if kind is None:
        kind = contract.BasisKind.FRONT_F1_SHIFTED
    return contract.DofKey(source_id, contract.FieldComponent.UZ, kind, 1)


def frozen_registry(*keys):
    registry = contract.DofRegistry()
    for key in keys:
        registry.register(key)
    registry.freeze()
    return registry


class FakeFieldContractTests(unittest.TestCase):
    def test_insertion_order_does_not_change_canonical_ids_or_slots(self):
        keys = (
            base_key("source-z"),
            front_key("source-z", contract.BasisKind.FRONT_F1_SHIFTED),
            front_key("source-a", contract.BasisKind.FRONT_F1_UNSHIFTED),
            front_key("source-a", contract.BasisKind.HEAVISIDE),
        )
        forward = frozen_registry(*keys)
        reverse = frozen_registry(*reversed(keys))
        self.assertEqual(forward.canonical_ids, reverse.canonical_ids)
        for key in keys:
            self.assertEqual(forward.slot_of(key), reverse.slot_of(key))

    def test_base_front_and_heaviside_identities_share_one_source(self):
        keys = (
            base_key("same-source"),
            front_key("same-source", contract.BasisKind.FRONT_F1_SHIFTED),
            front_key("same-source", contract.BasisKind.FRONT_F1_UNSHIFTED),
            front_key("same-source", contract.BasisKind.HEAVISIDE),
        )
        registry = frozen_registry(*keys)
        ids = tuple(key.canonical_id for key in keys)
        self.assertEqual(registry.source_count, 1)
        self.assertEqual(registry.dof_count, 4)
        self.assertEqual(len(set(ids)), 4)
        self.assertEqual(len(set(registry.canonical_ids)), 4)

    def test_fake_field_row_matches_registry_and_keeps_per_coefficient_columns(self):
        registry = frozen_registry(base_key("s"), front_key("s"))
        sample = contract.AntiplaneShapeSample(
            coefficient_ids=registry.canonical_ids,
            contributions=(
                contract.ScalarShapeContribution(2.0, (1.0, 2.0, 3.0)),
                contract.ScalarShapeContribution(4.0, (4.0, 5.0, 6.0)),
            ),
            status=contract.ProviderStatus.OK,
        )
        self.assertIs(contract.validate_antiplane_sample(sample, registry), sample)
        self.assertEqual(len(sample.contributions), registry.dof_count)
        self.assertEqual(sample.contributions[0].value, 2.0)
        self.assertEqual(sample.contributions[1].value, 4.0)

    def test_malformed_or_mismatched_ids_fail_closed(self):
        registry = frozen_registry(base_key("s"), front_key("s"))
        contributions = (
            contract.ScalarShapeContribution(1.0, (0.0, 0.0, 0.0)),
            contract.ScalarShapeContribution(2.0, (0.0, 0.0, 0.0)),
        )
        reversed_sample = contract.AntiplaneShapeSample(
            tuple(reversed(registry.canonical_ids)), contributions, contract.ProviderStatus.OK
        )
        with self.assertRaises(contract.ProviderContractError):
            contract.validate_antiplane_sample(reversed_sample, registry)

        list_ids_sample = contract.AntiplaneShapeSample(
            list(registry.canonical_ids), contributions, contract.ProviderStatus.OK
        )
        with self.assertRaises(contract.ProviderContractError):
            contract.validate_antiplane_sample(list_ids_sample, registry)

        wrong_length_sample = contract.AntiplaneShapeSample(
            registry.canonical_ids,
            (contract.ScalarShapeContribution(1.0, (0.0, 0.0)),),
            contract.ProviderStatus.OK,
        )
        with self.assertRaises(contract.ProviderContractError):
            contract.validate_antiplane_sample(wrong_length_sample, registry)

        ux_registry = frozen_registry(base_key("ux-source", contract.FieldComponent.UX))
        ux_sample = contract.AntiplaneShapeSample(
            ux_registry.canonical_ids,
            (contract.ScalarShapeContribution(1.0, (0.0, 0.0, 0.0)),),
            contract.ProviderStatus.OK,
        )
        with self.assertRaises(contract.ProviderContractError):
            contract.validate_antiplane_sample(ux_sample, ux_registry)

    def test_non_ok_rows_refuse_partial_data(self):
        registry = frozen_registry(base_key("s"), front_key("s"))
        unknown = contract.AntiplaneShapeSample((), (), contract.ProviderStatus.UNKNOWN, "front singularity")
        blocked = contract.AntiplaneShapeSample((), (), contract.ProviderStatus.BLOCKED, "unsupported point")
        self.assertIs(contract.validate_antiplane_sample(unknown, registry), unknown)
        self.assertIs(contract.validate_antiplane_sample(blocked, registry), blocked)

        partial = contract.AntiplaneShapeSample(
            (),
            (contract.ScalarShapeContribution(0.0, (0.0, 0.0, 0.0)),),
            contract.ProviderStatus.UNKNOWN,
            "partial provider result",
        )
        with self.assertRaises(contract.ProviderContractError):
            contract.validate_antiplane_sample(partial, registry)

        missing_reason = contract.AntiplaneShapeSample((), (), contract.ProviderStatus.UNKNOWN)
        with self.assertRaises(contract.ProviderContractError):
            contract.validate_antiplane_sample(missing_reason, registry)

    def test_engineering_strain_order_matches_frozen_elastic_convention(self):
        strain = contract.antiplane_engineering_strain((1.0, 2.0, 3.0))
        self.assertEqual(strain, (0.0, 0.0, 3.0, 0.0, 2.0, 1.0))

    def test_dof_and_source_cap_failures_leave_registry_unchanged(self):
        dof_limited = contract.DofRegistry()
        components = (
            contract.FieldComponent.UX,
            contract.FieldComponent.UY,
            contract.FieldComponent.UZ,
        )
        for source_index in range(contract.MAX_SOURCE_NODES):
            for component in components:
                dof_limited.register(base_key(f"source-{source_index:03d}", component))
        self.assertEqual(dof_limited.source_count, 125)
        self.assertEqual(dof_limited.dof_count, 375)
        before = (dof_limited.source_count, dof_limited.dof_count, dof_limited.state)
        with self.assertRaisesRegex(ValueError, "dof_count exceeds 375"):
            dof_limited.register(front_key("source-000"))
        self.assertEqual((dof_limited.source_count, dof_limited.dof_count, dof_limited.state), before)

        source_limited = contract.DofRegistry()
        for source_index in range(contract.MAX_SOURCE_NODES):
            source_limited.register(base_key(f"node-{source_index:03d}"))
        before = (source_limited.source_count, source_limited.dof_count, source_limited.state)
        with self.assertRaisesRegex(ValueError, "source_count exceeds 125"):
            source_limited.register(base_key("node-125"))
        self.assertEqual((source_limited.source_count, source_limited.dof_count, source_limited.state), before)

    def test_quadrature_budget_and_dof_budget_caps(self):
        contract.validate_problem_budget(125, 375, 256)
        for counts, reason in (
            ((125, 376, 1), "dof_count exceeds 375"),
            ((125, 375, 257), "quadrature_count exceeds 256"),
            ((126, 126, 1), "source_count exceeds 125"),
        ):
            with self.subTest(counts=counts):
                with self.assertRaisesRegex(ValueError, reason):
                    contract.validate_problem_budget(*counts)


def _flatten(suite):
    for item in suite:
        if isinstance(item, unittest.TestSuite):
            yield from _flatten(item)
        else:
            yield item


def _write_junit(result, tests):
    RESULT_DIR.mkdir(parents=True, exist_ok=True)
    failures = {id(test): trace for test, trace in result.failures}
    errors = {id(test): trace for test, trace in result.errors}
    skipped = {id(test): reason for test, reason in result.skipped}
    root = ET.Element(
        "testsuite",
        {
            "name": "enriched_dof_contract_fake_field",
            "tests": str(len(tests)),
            "failures": str(len(result.failures)),
            "errors": str(len(result.errors)),
            "skipped": str(len(result.skipped)),
            "time": "0",
        },
    )
    for test in tests:
        case = ET.SubElement(
            root,
            "testcase",
            {"classname": type(test).__name__, "name": test._testMethodName, "time": "0"},
        )
        if id(test) in failures:
            ET.SubElement(case, "failure", {"message": "assertion failure"}).text = failures[id(test)]
        elif id(test) in errors:
            ET.SubElement(case, "error", {"message": "test error"}).text = errors[id(test)]
        elif id(test) in skipped:
            ET.SubElement(case, "skipped", {"message": skipped[id(test)]})
    ET.ElementTree(root).write(RESULT_DIR / "junit.xml", encoding="utf-8", xml_declaration=True)


def main() -> int:
    if not SNAPSHOT_PATH.is_file():
        print(f"frozen source snapshot missing: {SNAPSHOT_PATH}", file=sys.stderr)
        return 2
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(FakeFieldContractTests)
    tests = list(_flatten(suite))
    result = unittest.TextTestRunner(stream=sys.stdout, verbosity=2).run(suite)
    _write_junit(result, tests)
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    raise SystemExit(main())

