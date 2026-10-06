"""Render one bounded figure from the pinned failed-run NPZ after audit approval."""
from __future__ import annotations
import argparse, hashlib, json, math
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
PINS_PATH = HERE / "RENDER_PINS.json"
NPZ_REL = "pgfrag/tmp/review_20261006/mixed93_cached_lefm_solve_20261006/run01/arrays.npz"
NPZ_SHA = "3e1e6fd17751e4e06113010f20ea33c7b69914ca90431e26bbf76124951b5d41"
NPZ_BYTES = 1569232
CASES = ("pure_I", "pure_II", "pure_III", "mixed_70_40_100")
OUTPUT = HERE / "mixed93_saved_stress_norms_2x2.png"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""): digest.update(block)
    return digest.hexdigest()


def verify_audit(spec: dict) -> dict:
    if any(spec.get(key) is None for key in ("path", "sha256", "schema", "accepted_status")):
        raise RuntimeError("independent recovery audit is not pinned yet")
    path = (ROOT / spec["path"]).resolve()
    path.relative_to(ROOT)
    if sha256(path) != spec["sha256"]: raise RuntimeError("independent recovery audit SHA mismatch")
    audit = json.loads(path.read_text(encoding="utf-8"))
    if (audit.get("schema"), audit.get("status"), audit.get("passed")) != (spec["schema"], spec["accepted_status"], True):
        raise RuntimeError("independent recovery audit schema/status/pass gate failed")
    if audit.get("input_npz_path") != NPZ_REL or audit.get("input_npz_sha256") != NPZ_SHA:
        raise RuntimeError("independent recovery audit is not bound to the frozen NPZ")
    return audit


def render(npz_path: Path, audit: dict) -> dict:
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    with np.load(npz_path, allow_pickle=False) as data:
        xyz = np.asarray(data["offgrid_xyz_m"], dtype=np.float64)
        fields = {case: (np.asarray(data[f"{case}_offgrid_exact_stress"], dtype=np.float64),
                         np.asarray(data[f"{case}_offgrid_solution_stress_Pa"], dtype=np.float64)) for case in CASES}
    if xyz.shape != (128, 3) or not np.isfinite(xyz).all(): raise RuntimeError("saved common off-grid coordinates must be finite 128x3")
    stats = {}
    for case, (exact, solved) in fields.items():
        if exact.shape != (128, 6) or solved.shape != (128, 6) or not np.isfinite(exact).all() or not np.isfinite(solved).all():
            raise RuntimeError(f"saved stress shape/finiteness gate failed for {case}")
        x = np.sqrt(np.sum(exact[:, :3] ** 2, axis=1) + 2.0 * np.sum(exact[:, 3:] ** 2, axis=1))
        y = np.sqrt(np.sum(solved[:, :3] ** 2, axis=1) + 2.0 * np.sum(solved[:, 3:] ** 2, axis=1))
        diff = solved - exact
        rel = math.sqrt(float(np.sum(diff[:, :3] ** 2) + 2.0 * np.sum(diff[:, 3:] ** 2))) / max(
            math.sqrt(float(np.sum(exact[:, :3] ** 2) + 2.0 * np.sum(exact[:, 3:] ** 2))), np.finfo(float).tiny)
        stats[case] = {"relative_tensor_l2": rel, "analytic_norm_min_Pa": float(x.min()),
                       "analytic_norm_max_Pa": float(x.max()), "solved_norm_min_Pa": float(y.min()),
                       "solved_norm_max_Pa": float(y.max())}
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))
    for ax, case in zip(axes.flat, CASES):
        exact, solved = fields[case]
        x = np.sqrt(np.sum(exact[:, :3] ** 2, axis=1) + 2.0 * np.sum(exact[:, 3:] ** 2, axis=1))
        y = np.sqrt(np.sum(solved[:, :3] ** 2, axis=1) + 2.0 * np.sum(solved[:, 3:] ** 2, axis=1))
        span = max(float(x.max()), float(y.max()), 1.0); lo = min(float(x.min()), float(y.min()))
        hi = max(float(x.max()), float(y.max())); pad = 0.04 * max(hi - lo, span * 1.0e-6)
        ax.scatter(x, y, s=20, alpha=0.8)
        ax.plot((lo - pad, hi + pad), (lo - pad, hi + pad), "k--", linewidth=1)
        ax.set(xlim=(lo - pad, hi + pad), ylim=(lo - pad, hi + pad), xlabel="analytic stress Frobenius norm (Pa)",
               ylabel="solved stress Frobenius norm (Pa)", title=case.replace("_", " "))
        ax.text(0.04, 0.96, f"tensor relative L2 = {stats[case]['relative_tensor_l2']:.3%}\n128 shared saved points",
                transform=ax.transAxes, va="top", fontsize=9)
        ax.set_aspect("equal", adjustable="box")
    fig.suptitle("Manufactured elasticity · no crack growth · producer rc=1, recovered arrays\nIndependent audit: " + audit["status"])
    fig.tight_layout(rect=(0, 0, 1, 0.93)); fig.savefig(OUTPUT, dpi=150); plt.close(fig)
    return {"output": OUTPUT.name, "bytes": OUTPUT.stat().st_size, "case_metrics": stats,
            "offgrid_rows": int(len(xyz)), "audit_status": audit["status"]}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--render", action="store_true")
    args = parser.parse_args(argv)
    if not args.render: print("SOURCE_READY_RENDER_DISABLED"); return 0
    pins = json.loads(PINS_PATH.read_text(encoding="utf-8"))
    if pins.get("execution_authorized") is not True: raise RuntimeError("saved renderer is not authorized")
    if pins.get("source_sha256") != sha256(Path(__file__).resolve()): raise RuntimeError("renderer source SHA mismatch")
    npz_path = ROOT / NPZ_REL
    if sha256(npz_path) != NPZ_SHA or npz_path.stat().st_size != NPZ_BYTES: raise RuntimeError("failed-run NPZ SHA/size mismatch")
    audit = verify_audit(pins.get("independent_audit", {}))
    if OUTPUT.exists(): raise RuntimeError("refusing to overwrite an existing PNG")
    result = render(npz_path, audit)
    if result["bytes"] > 2 * 1024 * 1024: raise RuntimeError("PNG exceeds 2 MiB output cap")
    print(json.dumps(result, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__": raise SystemExit(main())
