"""Run the v2 targeted tests and emit explicit state-machine review artifacts."""

from __future__ import annotations

import contextlib
import hashlib
import io
import json
import sys
import time
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import griffith_front_energy_gate_v2 as gate  # noqa: E402
import griffith_front_step as front  # noqa: E402
from test_griffith_front_energy_gate_v2 import (  # noqa: E402
    GriffithFrontEnergyGateV2Tests,
    base_evidence,
    setup_proposal,
)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _snapshot_case(name, snapshot, *, extra=None):
    row = {
        "name": name,
        "origin": "EXPLICIT_FIXTURE_NOT_SOLVER",
        "active_surface_version_id": snapshot.active_surface.version_id,
        "pending_candidate_id": snapshot.pending_candidate_id,
        "last_controller_status": snapshot.controller.event_ledger[-1].status if snapshot.controller.event_ledger else "INITIAL",
        "last_energy_gate_status": snapshot.energy_gate_ledger[-1].status if snapshot.energy_gate_ledger else None,
        "snapshot": snapshot.to_payload(include_surfaces=True),
    }
    if extra:
        row.update(extra)
    return row


def make_states():
    initial, proposed, candidate, state = setup_proposal()
    states = [_snapshot_case("proposed_geometry_waiting_for_energy_gate", proposed, extra={
        "candidate_geometry_only": True,
        "active_surface_still_parent": proposed.active_surface.version_id == initial.active_surface.version_id,
    })]

    balanced = gate.confirm_candidate(
        proposed, candidate.candidate_id,
        reanalyzed_surface=state,
        evidence=base_evidence(proposed, candidate, state),
    )
    states.append(_snapshot_case("balanced_fixture_committed_fixture_only", balanced))

    external_evidence = base_evidence(
        proposed, candidate, state,
        evidence_origin="EXTERNAL_RECEIPT",
        mechanical_status="REANALYZED_CURRENT",
        caller_provenance=(("test_only", "synthetic hashes; no solver invoked"),),
    )
    external_path = gate.confirm_candidate(
        proposed, candidate.candidate_id,
        reanalyzed_surface=state, evidence=external_evidence,
    )
    states.append(_snapshot_case("external_api_branch_synthetic_test_only", external_path, extra={
        "not_solver_output": True,
        "caller_provenance_test_only": True,
    }))

    positive = gate.confirm_candidate(
        proposed, candidate.candidate_id,
        reanalyzed_surface=state,
        evidence=base_evidence(
            proposed, candidate, state,
            external_work_increment_J=1.0e6, abs_tol_J=1.0e-6,
        ),
    )
    states.append(_snapshot_case("large_positive_residual_rejected", positive))

    negative = gate.confirm_candidate(
        proposed, candidate.candidate_id,
        reanalyzed_surface=state,
        evidence=base_evidence(
            proposed, candidate, state,
            delta_stored_elastic_energy_J=1.0e6, abs_tol_J=1.0e-6,
        ),
    )
    states.append(_snapshot_case("large_negative_residual_rejected", negative))

    no_tolerance = gate.confirm_candidate(
        proposed, candidate.candidate_id,
        reanalyzed_surface=state,
        evidence=base_evidence(proposed, candidate, state, abs_tol_J=None),
    )
    states.append(_snapshot_case("missing_tolerance_rejected", no_tolerance))

    total_potential = gate.confirm_candidate(
        proposed, candidate.candidate_id,
        reanalyzed_surface=state,
        evidence=base_evidence(
            proposed, candidate, state,
            delta_stored_elastic_energy_J=None,
            delta_total_potential_J=-0.01,
        ),
    )
    states.append(_snapshot_case("legacy_total_potential_not_bridged", total_potential))

    dynamic = gate.confirm_candidate(
        proposed, candidate.candidate_id,
        reanalyzed_surface=state,
        evidence=base_evidence(proposed, candidate, state, physics_mode="DYNAMIC"),
    )
    states.append(_snapshot_case("dynamic_mode_unknown", dynamic))

    wrong_parent = gate.confirm_candidate(
        proposed, candidate.candidate_id,
        reanalyzed_surface=state,
        evidence=base_evidence(
            proposed, candidate, state,
            bound_parent_surface_state_sha256="0" * 64,
        ),
    )
    states.append(_snapshot_case("wrong_parent_binding_rejected", wrong_parent))
    return states


def make_svg() -> str:
    return """<svg xmlns="http://www.w3.org/2000/svg" width="1320" height="455" viewBox="0 0 1320 455">
<rect width="1320" height="455" fill="#fbfcfe"/>
<style>text{font-family:Segoe UI,Arial,sans-serif;fill:#1e293b}.title{font-size:23px;font-weight:700}.head{font-size:17px;font-weight:700}.sub{font-size:13px}.foot{font-size:13px;fill:#7c2d12}.formula{font-size:15px;font-family:Consolas,monospace}</style>
<text x="32" y="38" class="title">Quasi-static outer energy gate before surface activation</text>
<text x="32" y="63" class="sub">The candidate strip remains pending until a current reanalysis and an explicit numeric tolerance gate both agree.</text>
<rect x="25" y="91" width="300" height="286" rx="12" fill="white" stroke="#d3dce6"/>
<text x="43" y="120" class="head">1  PROPOSED</text><text x="43" y="142" class="sub">old surface stays active</text>
<rect x="90" y="181" width="145" height="112" fill="#dceaf3" stroke="#28628a" stroke-width="3"/>
<rect x="235" y="181" width="34" height="112" fill="#f6d89e" fill-opacity="0.6" stroke="#d59027" stroke-width="3" stroke-dasharray="8 6"/>
<line x1="235" y1="168" x2="235" y2="305" stroke="#16202a" stroke-width="2"/>
<text x="238" y="317" class="sub">proposed +x strip</text>
<text x="43" y="350" class="sub">No displacement field is invented.</text>
<path d="M337 232 L365 232" stroke="#64748b" stroke-width="3"/><path d="M365 232 L354 225 M365 232 L354 239" stroke="#64748b" stroke-width="3"/>
<rect x="380" y="91" width="410" height="286" rx="12" fill="white" stroke="#d3dce6"/>
<text x="398" y="120" class="head">2  SUPPORTED QUASI-STATIC BALANCE</text>
<text x="405" y="177" class="formula">R = W_ext − ΔU_stored − Gc·ΔA</text>
<text x="405" y="215" class="formula">|R| ≤ abs_tol + relative_tol·E_ref</text>
<text x="405" y="262" class="sub">ΔΠ = ΔU_stored − W_ext</text>
<text x="405" y="284" class="sub">ΔΠ is derived for reporting; it is not</text>
<text x="405" y="304" class="sub">subtracted again in R.</text>
<text x="398" y="350" class="sub">All terms in J; E_ref &gt; 0; tolerances ≥ 0.</text>
<path d="M802 232 L830 232" stroke="#64748b" stroke-width="3"/><path d="M830 232 L819 225 M830 232 L819 239" stroke="#64748b" stroke-width="3"/>
<rect x="845" y="91" width="215" height="286" rx="12" fill="white" stroke="#d3dce6"/>
<text x="862" y="120" class="head">3  FAIL / UNKNOWN</text>
<text x="862" y="151" class="sub">Missing terms/tolerance,</text><text x="862" y="172" class="sub">bad parent hash, or |R|</text><text x="862" y="193" class="sub">above tolerance:</text>
<text x="862" y="238" class="head" fill="#a32929">NEEDS_RECOMPUTE</text>
<text x="862" y="267" class="sub">Dynamic / mixed cohesive:</text>
<text x="862" y="290" class="head">UNKNOWN</text>
<text x="862" y="350" class="sub">Old committed surface retained.</text>
<path d="M1072 232 L1100 232" stroke="#64748b" stroke-width="3"/><path d="M1100 232 L1089 225 M1100 232 L1089 239" stroke="#64748b" stroke-width="3"/>
<rect x="1115" y="91" width="180" height="286" rx="12" fill="white" stroke="#d3dce6"/>
<text x="1132" y="120" class="head">4  PASS</text>
<text x="1132" y="151" class="sub">Energy and binding gate</text><text x="1132" y="172" class="sub">pass: activate new state.</text>
<text x="1132" y="218" class="head" fill="#28734b">COMMITTED_*</text>
<text x="1132" y="251" class="sub">Fixture: fixture-only.</text>
<text x="1132" y="273" class="sub">External: experimental</text><text x="1132" y="294" class="sub">reanalysis only.</text>
<text x="1132" y="350" class="sub">Never physical validation.</text>
<text x="32" y="423" class="foot">Diagram flow is exercised with explicit in-memory fixtures only. No solver/native run. Surface contract remains NOT_COUPLED; v1 delta_potential_J is NOT_BRIDGED.</text>
</svg>"""


def main() -> int:
    out_dir = Path(__file__).resolve().parents[2] / "tmp" / "review_20261006" / "griffith_front_energy_gate_v2_20261006"
    out_dir.mkdir(parents=True, exist_ok=True)
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(GriffithFrontEnergyGateV2Tests)
    stream = io.StringIO()
    start = time.perf_counter()
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    elapsed = time.perf_counter() - start
    states = make_states()
    state_payload = {
        "schema_version": "pgfrag.griffith-front-energy-gate-v2.fixture-states.v1",
        "provenance": "EXPLICIT_FIXTURE_NOT_SOLVER",
        "meaning": "state-machine tests only; external-receipt branch uses synthetic test hashes",
        "legacy_delta_potential_bridge": "NOT_BRIDGED",
        "state_count": len(states),
        "states": states,
    }
    state_path = out_dir / "fixture_states.json"
    state_path.write_text(json.dumps(state_payload, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    svg_path = out_dir / "quasistatic_energy_gate_flow.svg"
    svg_path.write_text(make_svg() + "\n", encoding="utf-8")
    v1_path = _V1_DIR / "griffith_front_step.py"
    surface_path = _V1_DIR.parent / "explicit_surface_state_20261005" / "explicit_surface_state.py"
    receipt = {
        "schema_version": "pgfrag.griffith-front-energy-gate-v2.actual-test-receipt.v1",
        "execution": "targeted Python unittest fixtures",
        "python": sys.version.split()[0],
        "test_count": result.testsRun,
        "failure_count": len(result.failures),
        "error_count": len(result.errors),
        "elapsed_seconds": elapsed,
        "limit_seconds": 10.0,
        "fixture_state_count": len(states),
        "limit_fixture_states": 64,
        "solver_invoked": False,
        "native_runtime_invoked": False,
        "production_writes": False,
        "external_commit_test_uses_synthetic_hashes": True,
        "external_test_metadata": "caller_provenance marks test_only; this is API-branch coverage, not solver evidence",
        "physical_validation_claim": False,
        "source_sha256": {
            "griffith_front_energy_gate_v2.py": _sha(HERE / "griffith_front_energy_gate_v2.py"),
            "test_griffith_front_energy_gate_v2.py": _sha(HERE / "test_griffith_front_energy_gate_v2.py"),
            "run_fixture.py": _sha(Path(__file__).resolve()),
        },
        "read_only_dependency_sha256": {
            "griffith_front_step.py": _sha(v1_path),
            "explicit_surface_state.py": _sha(surface_path),
        },
        "fixture_states_sha256": _sha(state_path),
        "diagram_sha256": _sha(svg_path),
        "test_output": stream.getvalue(),
    }
    receipt_path = out_dir / "ACTUAL_TEST_RECEIPT.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(f"tests={result.testsRun} failures={len(result.failures)} errors={len(result.errors)} elapsed_s={elapsed:.6f}")
    print(f"fixture_states={state_path}")
    print(f"diagram={svg_path}")
    print(f"receipt={receipt_path}")
    return 0 if result.wasSuccessful() and elapsed <= 10.0 and len(states) <= 64 else 1


_V1_DIR = Path(__file__).resolve().parents[1] / "griffith_front_step_20261006"


if __name__ == "__main__":
    raise SystemExit(main())

