from __future__ import annotations

import argparse
import csv
import hashlib
import math
import sys
import tomllib
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
DEFAULT_OUTPUT = REPO_ROOT.parent / "relicweave-synlab" / "tmp" / "review_20261003" / "implementation" / "pg" / "native_contact_minimal"
CASE_NAMES = ("contact_baseline", "contact_halfdt", "nocontact_control")
RHO = 1000.0


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_case(root: Path, name: str) -> tuple[dict, list[dict[str, str]]]:
    case_dir = root / name
    with (case_dir / "summary.toml").open("rb") as stream:
        summary = tomllib.load(stream)
    with (case_dir / "points.csv").open("r", newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return summary, rows


def vec(row: dict[str, str], prefix: str) -> tuple[float, float, float]:
    return tuple(float(row[f"{prefix}_{axis}_m"]) for axis in "xyz")


def norm(a: tuple[float, float, float]) -> float:
    return math.sqrt(sum(x * x for x in a))


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def scale(a, value):
    return tuple(x * value for x in a)


def body_sums(rows, vector_fn):
    totals = {"a": (0.0, 0.0, 0.0), "b": (0.0, 0.0, 0.0)}
    for row in rows:
        body = row["body_name"]
        totals[body] = add(totals[body], scale(vector_fn(row), float(row["mass_kg"])))
    return totals


def center_z(rows: list[dict[str, str]], body: str, prefix: str) -> float:
    selected = [r for r in rows if r["body_name"] == body]
    mass = sum(float(r["mass_kg"]) for r in selected)
    return sum(float(r["mass_kg"]) * float(r[f"{prefix}_z_m"]) for r in selected) / mass


def contact_oracle(rows, summary):
    """Re-evaluate the vendored endpoint contact equation; this does not alter native state."""
    a_rows = [r for r in rows if r["body_name"] == "a"]
    b_rows = [r for r in rows if r["body_name"] == "b"]
    radius = float(summary["contact_radius_m"])
    penalty = float(summary["contact_penalty_factor"])
    horizon = float(summary["horizon_m"])
    volume = float(summary["point_spacing_m"]) ** 3
    c0 = 9.0 * penalty / math.pi
    force_a = (0.0, 0.0, 0.0)
    force_b = (0.0, 0.0, 0.0)
    pair_distances = []
    if summary["contact_enabled"]:
        for ra in a_rows:
            xa = vec(ra, "final")
            for rb in b_rows:
                delta = sub(vec(rb, "final"), xa)
                distance = norm(delta)
                if distance < radius:
                    pair_distances.append(distance)
                    # Source: b_int[a] += C*(r-L)/L * V[b] * (x[b]-x[a]).
                    factor = c0 / horizon**4 * (radius - distance) / distance
                    pair_force = scale(delta, factor * volume * volume)
                    force_a = add(force_a, pair_force)
                    force_b = sub(force_b, pair_force)
    return {
        "active_pairs": len(pair_distances),
        "minimum_pair_distance_m": min(pair_distances) if pair_distances else None,
        "force_a_N": force_a, "force_b_N": force_b,
        "action_reaction_residual_N": norm(add(force_a, force_b)),
    }


def analyze_case(summary, rows):
    initial = body_sums(rows, lambda r: tuple(float(r[f"initial_v{x}_m_per_s"]) for x in "xyz"))
    final = body_sums(rows, lambda r: tuple(float(r[f"final_v{x}_m_per_s"]) for x in "xyz"))
    delta_a, delta_b = sub(final["a"], initial["a"]), sub(final["b"], initial["b"])
    total_before, total_after = add(initial["a"], initial["b"]), add(final["a"], final["b"])
    force = body_sums(rows, lambda r: tuple(float(r[f"b_int_{x}_N_per_m3"]) / RHO for x in "xyz"))
    free_error = 0.0
    for row in rows:
        for axis in "xyz":
            expected = float(row[f"rest_{axis}_m"]) + float(row[f"initial_v{axis}_m_per_s"]) * float(summary["end_time_s"])
            free_error = max(free_error, abs(float(row[f"final_{axis}_m"]) - expected))
    return {
        "summary": summary, "rows": rows, "oracle": contact_oracle(rows, summary),
        "delta_momentum_a": delta_a, "delta_momentum_b": delta_b,
        "total_momentum_drift": norm(sub(total_after, total_before)),
        "delta_momentum_residual": norm(add(delta_a, delta_b)),
        "native_force_a": force["a"], "native_force_b": force["b"],
        "native_force_residual": norm(add(force["a"], force["b"])),
        "free_motion_error": free_error,
    }


def draw_font(size: int):
    for path in (r"C:\Windows\Fonts\segoeui.ttf", r"C:\Windows\Fonts\arial.ttf"):
        if Path(path).is_file():
            return ImageFont.truetype(path, size)
    return ImageFont.load_default()


def write_xz_png(rows, out: Path):
    image = Image.new("RGB", (1040, 690), "white")
    draw = ImageDraw.Draw(image)
    font, small = draw_font(21), draw_font(14)
    left, top, width, height = 105, 110, 850, 490
    color = {"a": (35, 105, 166), "b": (211, 104, 44)}
    xs = [float(r[f"{s}_x_m"]) * 1000 for r in rows for s in ("rest", "final")]
    zs = [float(r[f"{s}_z_m"]) * 1000 for r in rows for s in ("rest", "final")]
    xmin, xmax = min(xs)-1.0, max(xs)+1.0
    zmin, zmax = min(zs)-1.0, max(zs)+1.0
    factor = min(width/(xmax-xmin), height/(zmax-zmin))
    xmid, zmid = (xmin+xmax)/2, (zmin+zmax)/2
    cx, cy = left+width/2, top+height/2

    def xy(xmm, zmm):
        return round(cx+(xmm-xmid)*factor), round(cy-(zmm-zmid)*factor)

    draw.text((left, 28), "Native point contact: X-Z projection", fill="black", font=font)
    draw.text((left, 61), "All discrete points shown; y is projected. Baseline reference to endpoint; axes in mm.", fill=(70,70,70), font=small)
    for tick in range(math.floor(xmin), math.ceil(xmax)+1):
        x, _ = xy(tick, zmid)
        draw.line((x, top, x, top+height), fill=(233,236,239), width=1)
        draw.text((x-10, top+height+8), str(tick), fill=(60,60,60), font=small)
    for tick in range(math.floor(zmin), math.ceil(zmax)+1):
        _, y = xy(xmid, tick)
        draw.line((left, y, left+width, y), fill=(233,236,239), width=1)
        draw.text((left-37, y-8), str(tick), fill=(60,60,60), font=small)
    draw.rectangle((left,top,left+width,top+height), outline=(60,60,60), width=2)
    draw.text((left+width//2-28, top+height+37), "x (mm)", fill="black", font=small)
    draw.text((22, top+height//2), "z (mm)", fill="black", font=small)
    for row in rows:
        body = row["body_name"]
        p0 = xy(float(row["rest_x_m"])*1000, float(row["rest_z_m"])*1000)
        p1 = xy(float(row["final_x_m"])*1000, float(row["final_z_m"])*1000)
        faded = tuple(int(c*.68+255*.32) for c in color[body])
        draw.line((*p0,*p1), fill=faded, width=2)
        draw.ellipse((p0[0]-6,p0[1]-6,p0[0]+6,p0[1]+6), outline=(95,95,95), width=2)
        draw.ellipse((p1[0]-6,p1[1]-6,p1[0]+6,p1[1]+6), fill=color[body], outline="white", width=1)
    draw.ellipse((left,645,left+12,657), outline=(95,95,95), width=2)
    draw.text((left+18,641), "reference", fill="black", font=small)
    draw.ellipse((left+125,645,left+137,657), fill=color["a"])
    draw.text((left+144,641), "endpoint a", fill="black", font=small)
    draw.ellipse((left+260,645,left+272,657), fill=color["b"])
    draw.text((left+279,641), "endpoint b", fill="black", font=small)
    image.save(out)


def write_3d_png(rows, out: Path):
    image = Image.new("RGB", (1040, 690), "white")
    draw = ImageDraw.Draw(image)
    font, small = draw_font(21), draw_font(13)
    left, top, width, height = 90, 115, 860, 500
    color = {"a": (35,105,166), "b": (211,104,44)}

    def project(point):
        x, y, z = (value*1000 for value in point)
        return (x-y)/math.sqrt(2), (x+y-2*z)/math.sqrt(6)

    projected = [project(vec(row,state)) for row in rows for state in ("rest","final")]
    umin, umax = min(p[0] for p in projected)-1, max(p[0] for p in projected)+1
    vmin, vmax = min(p[1] for p in projected)-1, max(p[1] for p in projected)+1
    factor = min(width/(umax-umin), height/(vmax-vmin))
    umid, vmid = (umin+umax)/2, (vmin+vmax)/2
    cx, cy = left+width/2, top+height/2

    def screen(q):
        return round(cx+(q[0]-umid)*factor), round(cy-(q[1]-vmid)*factor)

    draw.text((left,30), "Native point contact: 3D point cloud", fill="black", font=font)
    draw.text((left,64), "Orthographic isometric projection; x, y, z and movement are in mm.", fill=(70,70,70), font=small)
    origin = screen(project((0.0,0.0,0.0)))
    for axis, endpoint in (("x",(0.002,0,0)),("y",(0,0.002,0)),("z",(0,0,0.002))):
        end = screen(project(endpoint))
        draw.line((*origin,*end), fill=(90,90,90), width=2)
        draw.text((end[0]+4,end[1]-12), axis+" (mm)", fill=(55,55,55), font=small)
    for row in rows:
        body = row["body_name"]
        p0, p1 = screen(project(vec(row,"rest"))), screen(project(vec(row,"final")))
        faded = tuple(int(c*.72+255*.28) for c in color[body])
        draw.line((*p0,*p1), fill=faded, width=2)
        draw.ellipse((p0[0]-5,p0[1]-5,p0[0]+5,p0[1]+5), outline=(90,90,90), width=2)
        draw.ellipse((p1[0]-6,p1[1]-6,p1[0]+6,p1[1]+6), fill=color[body], outline="white", width=1)
    draw.ellipse((left,645,left+12,657), outline=(90,90,90), width=2)
    draw.text((left+18,641),"reference points",fill="black",font=small)
    draw.ellipse((left+145,645,left+157,657),fill=color["a"])
    draw.text((left+164,641),"endpoint a",fill="black",font=small)
    draw.ellipse((left+270,645,left+282,657),fill=color["b"])
    draw.text((left+289,641),"endpoint b",fill="black",font=small)
    image.save(out)


def toml_quote(value: str) -> str:
    return '"' + value.replace("\\", "\\\\").replace('"', '\\"').replace("\n", "\\n") + '"'


def verify_vendor(repo: Path, output: Path) -> None:
    project = repo / "src" / "pgfrag" / "resources" / "julia"
    sys.path.insert(0, str(repo / "src"))
    from pgfrag.backend import _verify_backend_provenance

    record = _verify_backend_provenance(project)
    lines = ["status = PASS", f"upstream_tag = {record['upstream_tag']}",
             f"backend_version = {record['backend_version']}",
             f"mechanics_source_files_verified = {record['mechanics_source_files_verified']}",
             f"headless_io_source_files_verified = {record['headless_io_source_files_verified']}",
             f"headless_patch_sha256 = {record['headless_patch_sha256']}",
             f"vendored_license_sha256 = {record['vendored_license_sha256']}",
             f"adapter_runner_sha256 = {record['adapter_runner_sha256']}",
             "verification = existing pgfrag.backend._verify_backend_provenance read-only helper",
             "contact_native_physics_sha256 = " + sha256(project/"vendor"/"Peridynamics.jl"/"src"/"physics"/"short_range_force_contact.jl")]
    output.mkdir(parents=True, exist_ok=True)
    (output/"vendor_integrity.txt").write_text("\n".join(lines)+"\n", encoding="utf-8")
    print("VENDOR_INTEGRITY_PASS", " ".join(lines[:8]))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("--verify-vendor","--gate-preflight","--summarize"))
    parser.add_argument("--repo-root", type=Path, default=REPO_ROOT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    repo, output = args.repo_root.resolve(), args.output_root.resolve()
    if args.mode == "--verify-vendor":
        verify_vendor(repo, output)
        return 0
    if args.mode == "--gate-preflight":
        summary, rows = read_case(output, "preflight")
        valid = (summary.get("submit_completed") is True and summary.get("case") == "preflight"
                 and summary.get("step_count") == 1 and summary.get("contact_enabled") is True
                 and len(rows) == 16 and {r["point_label"] for r in rows} == {f"{b}{i}" for b in "ab" for i in range(1,9)})
        if not valid:
            raise SystemExit("STOP: native one-step preflight did not return 16 uniquely identified point rows")
        print("PREFLIGHT_GATE_PASS: one native submit step and 16 points; this is not contact-success evidence")
        return 0

    cases = {}
    for name in CASE_NAMES:
        summary, rows = read_case(output, name)
        if summary.get("submit_completed") is not True or len(rows) != 16:
            raise SystemExit(f"STOP: malformed or incomplete native result for {name}")
        cases[name] = analyze_case(summary, rows)
    base, half, control = (cases[name] for name in CASE_NAMES)
    control_by_id = {r["point_label"]: r for r in control["rows"]}
    contact_control_delta = max(norm(sub(vec(r,"final"),vec(control_by_id[r["point_label"]],"final"))) for r in base["rows"])
    half_position_delta = max(norm(sub(vec(a,"final"),vec(b,"final"))) for a,b in zip(base["rows"],half["rows"]))
    half_velocity_delta = max(norm(tuple(float(a[f"final_v{x}_m_per_s"])-float(b[f"final_v{x}_m_per_s"]) for x in "xyz")) for a,b in zip(base["rows"],half["rows"]))
    center_separation = {name:center_z(case["rows"],"a","final")-center_z(case["rows"],"b","final") for name,case in cases.items()}
    pairs = {name:case["oracle"]["active_pairs"] for name,case in cases.items()}
    momentum_scale = sum(float(r["mass_kg"]) for r in base["rows"]) * 0.2
    momentum_relative = base["total_momentum_drift"] / max(momentum_scale,1e-30)
    force_scale = norm(base["native_force_a"])+norm(base["native_force_b"])
    force_relative = base["native_force_residual"] / max(force_scale,1e-30)
    classification = "POINT_CONTACT_NUMERICAL_PROTOTYPE" if pairs["contact_baseline"] > 0 and contact_control_delta > 1e-12 and center_separation["contact_baseline"] > 0 else "INSUFFICIENT_ENDPOINT_EVIDENCE"
    lines = [
        "Native contact minimal endpoint results",
        f"Classification: {classification}",
        "Three fixed jobs were submitted once each after the genuine native one-step preflight.",
        "Preflight proves native setup, include/compile, submit, and terminal extraction only; it is not contact-success evidence.",
        f"Endpoint active cross-body pairs: baseline={pairs['contact_baseline']}, half-dt={pairs['contact_halfdt']}, control={pairs['nocontact_control']}.",
        f"Baseline minimum active pair distance (m): {base['oracle']['minimum_pair_distance_m']}.",
        f"Baseline/control center separation (m): {center_separation['contact_baseline']}, {center_separation['nocontact_control']}.",
        f"Baseline max point-position difference from control (m): {contact_control_delta}.",
        f"dt/2 maximum endpoint position and velocity differences: {half_position_delta} m; {half_velocity_delta} m/s (reported without a preset tolerance).",
        f"Control free-motion maximum position error (m): {control['free_motion_error']}.",
        f"Baseline total momentum drift / (total mass x 0.2 m/s): {momentum_relative}.",
        f"Baseline native terminal body-force action-reaction residual (N): {base['native_force_residual']}; relative residual: {force_relative}.",
        f"Baseline opposite-body Δmomentum residual (kg m/s): {base['delta_momentum_residual']}.",
        "The endpoint pair-force oracle independently evaluates the vendored equation from terminal positions and point volumes. It is not a native force history or contact impulse.",
        "Contact path history, impulse from saved native force samples, and contact potential energy are UNKNOWN.",
        "Boundary: two 8-point BBMaterial cubes and discrete cross-body point contact only; no ceramic calibration, resolved surface contact, self-contact, or whole-vessel fracture.",
    ]
    (output/"summary.txt").write_text("\n".join(lines)+"\n",encoding="utf-8")
    with (output/"case_comparison.csv").open("w",newline="",encoding="utf-8") as stream:
        writer=csv.writer(stream)
        writer.writerow(("case","steps","dt_s","end_time_s","endpoint_active_pairs","center_separation_m","momentum_drift_kg_m_per_s","native_body_force_residual_N","oracle_force_a_N","oracle_force_b_N"))
        for name,case in cases.items():
            s,o=case["summary"],case["oracle"]
            writer.writerow((name,s["step_count"],s["solver_dt_s"],s["end_time_s"],o["active_pairs"],center_separation[name],case["total_momentum_drift"],case["native_force_residual"],norm(o["force_a_N"]),norm(o["force_b_N"])))
    write_xz_png(base["rows"],output/"native_points_xz.png")
    write_3d_png(base["rows"],output/"native_points_3d.png")

    project=repo/"src"/"pgfrag"/"resources"/"julia"
    hashes={
        "julia_script":SCRIPT_DIR/"run_native_contact_minimal.jl",
        "powershell_runner":SCRIPT_DIR/"run_native_contact_minimal.ps1",
        "python_summary":Path(__file__).resolve(),
        "julia_project":project/"Project.toml","julia_manifest":project/"Manifest.toml",
        "backend_provenance":repo/"backend_provenance.json",
        "contact_physics":project/"vendor"/"Peridynamics.jl"/"src"/"physics"/"short_range_force_contact.jl",
        "velocity_verlet":project/"vendor"/"Peridynamics.jl"/"src"/"time_solvers"/"velocity_verlet.jl",
        "multibody_data_handler":project/"vendor"/"Peridynamics.jl"/"src"/"core"/"threads_multibody_data_handler.jl",
        "submit":project/"vendor"/"Peridynamics.jl"/"src"/"core"/"submit.jl",
        "upstream_contact_unit_test":repo/"vendor"/"Peridynamics.jl"/"test"/"unit"/"physics"/"test_short_range_force_contact.jl",
        "upstream_license":repo/"vendor"/"Peridynamics.jl"/"LICENSE",
        "process_log":output/"process_log.csv","vendor_integrity_record":output/"vendor_integrity.txt",
    }
    for name in CASE_NAMES+("preflight",):
        hashes[name+"_summary"]=output/name/"summary.toml"
        hashes[name+"_points"]=output/name/"points.csv"
    prov=["[scope]","verification = \"existing pgfrag backend provenance helper plus explicit source SHA-256 records\"",
          "case_jobs = \"three prescribed jobs submitted once each after one native preflight\"",
          "preflight_contact_success = \"not claimed\"",""]
    for name,path in hashes.items():
        prov.extend((f"[sha256.{name}]",f"path = {toml_quote(str(path))}",f"value = {toml_quote(sha256(path))}",""))
    (output/"provenance.toml").write_text("\n".join(prov),encoding="utf-8")
    print("SUMMARY_WRITTEN",classification)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
