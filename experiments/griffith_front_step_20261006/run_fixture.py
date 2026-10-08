"""Run the bounded source fixture suite and emit inspectable state receipts."""

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
import griffith_front_step as front  # noqa: E402
from test_griffith_front_step import GriffithFrontStepTests, fixture_setup  # noqa: E402


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _json_sha(payload: object) -> str:
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")
    return hashlib.sha256(raw).hexdigest()


def _fixture_transition_snapshots() -> list[dict[str, object]]:
    initial, _, frame, bounds = fixture_setup(G=15.0, Gc=10.0)
    states: list[dict[str, object]] = []

    arrest = front.propose_step(
        initial,
        front.make_fixture_driving_evidence(initial.active_surface, G=7.0, Gc=10.0),
        frame=frame, requested_advance_m=0.01, step_cap_m=0.02, bounds=bounds,
    )
    states.append({
        "name": "arrest_fixture",
        "origin": "EXPLICIT_FIXTURE_NOT_SOLVER",
        "final_event_status": arrest.event_ledger[-1].status,
        "active_version_id": arrest.active_surface_version_id,
        "pending_candidate_id": arrest.pending_candidate_id,
        "snapshot": arrest.to_payload(include_surfaces=True),
    })

    evidence = front.make_fixture_driving_evidence(initial.active_surface, G=15.0, Gc=10.0)
    proposed = front.propose_step(
        initial, evidence, frame=frame,
        requested_advance_m=0.01, step_cap_m=0.02, bounds=bounds,
    )
    candidate = proposed.candidates[-1]
    states.append({
        "name": "proposed_geometry_not_committed",
        "origin": "EXPLICIT_FIXTURE_NOT_SOLVER",
        "final_event_status": proposed.event_ledger[-1].status,
        "active_version_id": proposed.active_surface_version_id,
        "pending_candidate_id": proposed.pending_candidate_id,
        "candidate_has_displacements": False,
        "snapshot": proposed.to_payload(include_surfaces=True),
    })

    incomplete_surface = front.make_fixture_reanalyzed_state(proposed.active_surface, candidate)
    incomplete = front.make_fixture_reanalysis_evidence(
        candidate, incomplete_surface, proposed.active_surface.state_sha256, energy_complete=False
    )
    needs_recompute = front.confirm_candidate(
        proposed, candidate.candidate_id,
        reanalyzed_surface=incomplete_surface, evidence=incomplete,
    )
    states.append({
        "name": "missing_energy_keeps_old_commit",
        "origin": "EXPLICIT_FIXTURE_NOT_SOLVER",
        "final_event_status": needs_recompute.event_ledger[-1].status,
        "active_version_id": needs_recompute.active_surface_version_id,
        "pending_candidate_id": needs_recompute.pending_candidate_id,
        "snapshot": needs_recompute.to_payload(include_surfaces=True),
    })

    fixture_surface = front.make_fixture_reanalyzed_state(
        proposed.active_surface, candidate, displacement_z_m=2.0e-6
    )
    complete = front.make_fixture_reanalysis_evidence(
        candidate, fixture_surface, proposed.active_surface.state_sha256
    )
    committed = front.confirm_candidate(
        proposed, candidate.candidate_id,
        reanalyzed_surface=fixture_surface, evidence=complete,
    )
    states.append({
        "name": "fixture_only_commit_branch",
        "origin": "EXPLICIT_FIXTURE_NOT_SOLVER",
        "final_event_status": committed.event_ledger[-1].status,
        "active_version_id": committed.active_surface_version_id,
        "pending_candidate_id": committed.pending_candidate_id,
        "snapshot": committed.to_payload(include_surfaces=True),
    })
    return states


def _svg(states: list[dict[str, object]]) -> str:
    proposed = next(state for state in states if state["name"] == "proposed_geometry_not_committed")
    pending = next(state for state in states if state["name"] == "missing_energy_keeps_old_commit")
    committed = next(state for state in states if state["name"] == "fixture_only_commit_branch")
    old_id = str(states[0]["active_version_id"])
    candidate_id = str(proposed["pending_candidate_id"])
    panels = [
        ("1  COMMITTED START", "Initial immutable surface", "#28628a", False, "active version: " + old_id[:29]),
        ("2  PROPOSED", "Geometry candidate only", "#d59027", True, "active stays: " + old_id[:25]),
        ("3  NEEDS RECOMPUTE", "Missing energy evidence", "#d59027", True, "pending: " + str(pending["pending_candidate_id"])[:28]),
        ("4  COMMITTED · FIXTURE", "Synthetic transition exercise", "#34885b", False, "candidate: " + candidate_id[:27]),
    ]
    parts = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1320" height="480" viewBox="0 0 1320 480">',
        '<rect width="1320" height="480" fill="#fbfcfe"/>',
        '<style>text{font-family:Segoe UI,Arial,sans-serif;fill:#1e293b}.title{font-size:23px;font-weight:700}.head{font-size:17px;font-weight:700}.sub{font-size:13px}.id{font-size:11px;fill:#526273}.foot{font-size:13px;fill:#7c2d12}</style>',
        '<text x="36" y="37" class="title">Finite Mode III front: proposed geometry is not yet the active surface</text>',
        '<text x="36" y="61" class="sub">Each pane is a state-machine snapshot. The amber strip is proposed; it turns green only after a bound reanalysis and complete energy evidence.</text>',
    ]
    panel_w = 300
    for i, (heading, subtitle, color, is_proposed, identity) in enumerate(panels):
        x = 24 + i * 324
        y = 92
        parts.append(f'<rect x="{x}" y="{y}" width="300" height="310" rx="12" fill="#ffffff" stroke="#d3dce6"/>')
        parts.append(f'<text x="{x+16}" y="{y+30}" class="head">{heading}</text>')
        parts.append(f'<text x="{x+16}" y="{y+51}" class="sub">{subtitle}</text>')
        old_x, old_y, old_w, old_h = x + 70, y + 112, 142, 122
        # The base crack region is solid blue in every panel.
        parts.append(f'<rect x="{old_x}" y="{old_y}" width="{old_w}" height="{old_h}" fill="#dceaf3" stroke="#28628a" stroke-width="3"/>')
        # Its current finite front is the right edge; growth is along +x.
        ext_x, ext_w = old_x + old_w, 34
        if i in (1, 2):
            parts.append(f'<rect x="{ext_x}" y="{old_y}" width="{ext_w}" height="{old_h}" fill="#f6d89e" fill-opacity="0.55" stroke="#d59027" stroke-width="3" stroke-dasharray="8 6"/>')
            parts.append(f'<text x="{old_x+old_w+4}" y="{old_y+old_h+19}" class="id">candidate strip</text>')
        elif i == 3:
            parts.append(f'<rect x="{ext_x}" y="{old_y}" width="{ext_w}" height="{old_h}" fill="#cfe9d8" stroke="#34885b" stroke-width="3"/>')
            parts.append(f'<text x="{old_x+old_w+1}" y="{old_y+old_h+19}" class="id">committed fixture strip</text>')
        if i == 3:
            parts.append(f'<line x1="{ext_x}" y1="{old_y-7}" x2="{ext_x}" y2="{old_y+old_h+7}" stroke="#8492a0" stroke-width="1.5" stroke-dasharray="4 4"/>')
            parts.append(f'<line x1="{ext_x+ext_w}" y1="{old_y-13}" x2="{ext_x+ext_w}" y2="{old_y+old_h+7}" stroke="#16202a" stroke-width="2"/>')
            front_x = ext_x + ext_w
        else:
            parts.append(f'<line x1="{ext_x}" y1="{old_y-13}" x2="{ext_x}" y2="{old_y+old_h+7}" stroke="#16202a" stroke-width="2"/>')
            front_x = ext_x
        parts.append(f'<text x="{front_x-20}" y="{old_y-20}" class="id">front</text>')
        parts.append(f'<text x="{front_x+8}" y="{old_y+old_h/2+5}" class="head">+x</text>')
        parts.append(f'<text x="{x+16}" y="{y+275}" class="id">{identity}</text>')
        parts.append(f'<text x="{x+16}" y="{y+294}" class="id">active version: {str(committed["active_version_id"] if i == 3 else old_id)[:28]}</text>')
    parts.append('<text x="36" y="435" class="foot">All values and the green commit path are explicit test fixtures, not a solver run or calibrated fracture result. Surface contract stays NOT_COUPLED.</text>')
    parts.append('</svg>')
    return "\n".join(parts)


def main() -> int:
    out_dir = Path(__file__).resolve().parents[2] / "tmp" / "review_20261006" / "griffith_front_step_20261006"
    out_dir.mkdir(parents=True, exist_ok=True)
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(GriffithFrontStepTests)
    output = io.StringIO()
    start = time.perf_counter()
    with contextlib.redirect_stdout(output):
        result = unittest.TextTestRunner(stream=output, verbosity=2).run(suite)
    elapsed = time.perf_counter() - start
    states = _fixture_transition_snapshots()
    states_payload = {
        "schema_version": "pgfrag.griffith-front-step.fixture-states.v1",
        "provenance": "EXPLICIT_FIXTURE_NOT_SOLVER",
        "meaning": "illustrative state-machine snapshots only; do not use as mechanics or calibrated evidence",
        "state_count": len(states),
        "states": states,
    }
    states_path = out_dir / "fixture_states.json"
    states_path.write_text(json.dumps(states_payload, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    svg_path = out_dir / "proposed_vs_committed.svg"
    svg_path.write_text(_svg(states) + "\n", encoding="utf-8")
    receipt = {
        "schema_version": "pgfrag.griffith-front-step.actual-test-receipt.v1",
        "execution": "bounded Python unittest fixture suite",
        "python": sys.version.split()[0],
        "test_count": result.testsRun,
        "failure_count": len(result.failures),
        "error_count": len(result.errors),
        "elapsed_seconds": elapsed,
        "fixture_state_count": len(states),
        "limit_fixture_states": 64,
        "solver_invoked": False,
        "native_runtime_invoked": False,
        "production_writes": False,
        "all_fixture_states_labelled": all(s["origin"] == "EXPLICIT_FIXTURE_NOT_SOLVER" for s in states),
        "source_sha256": {
            "griffith_front_step.py": _sha(HERE / "griffith_front_step.py"),
            "test_griffith_front_step.py": _sha(HERE / "test_griffith_front_step.py"),
            "run_fixture.py": _sha(Path(__file__).resolve()),
        },
        "fixture_states_sha256": _sha(states_path),
        "diagram_sha256": _sha(svg_path),
        "test_output": output.getvalue(),
    }
    receipt_path = out_dir / "ACTUAL_TEST_RECEIPT.json"
    receipt_path.write_text(json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(f"tests={result.testsRun} failures={len(result.failures)} errors={len(result.errors)} elapsed_s={elapsed:.6f}")
    print(f"fixture_states={states_path}")
    print(f"diagram={svg_path}")
    print(f"receipt={receipt_path}")
    return 0 if result.wasSuccessful() and len(states) <= 64 and elapsed < 30.0 else 1


if __name__ == "__main__":
    raise SystemExit(main())

