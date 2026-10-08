from __future__ import annotations

import ast
from difflib import unified_diff
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
REVIEW = Path(__file__).resolve().parent
V1 = REVIEW / "v1_freeze"
PAIRS = (
    (
        "bias_protocol.py",
        ROOT / "pgfrag" / "src" / "pgfrag" / "bias_protocol.py",
        "pgfrag/src/pgfrag/bias_protocol.py",
    ),
    (
        "test_bias_protocol.py",
        ROOT / "pgfrag" / "tests" / "test_bias_protocol.py",
        "pgfrag/tests/test_bias_protocol.py",
    ),
    (
        "protocol.md",
        REVIEW / "protocol.md",
        "pgfrag/tmp/review_20261004/bias_protocol/protocol.md",
    ),
)


def main() -> None:
    chunks: list[str] = []
    for snapshot_name, current_path, display_path in PAIRS:
        old_text = (V1 / snapshot_name).read_text(encoding="utf-8")
        new_text = current_path.read_text(encoding="utf-8")
        if current_path.suffix == ".py":
            ast.parse(new_text, filename=display_path)
        chunks.extend(
            unified_diff(
                old_text.splitlines(keepends=True),
                new_text.splitlines(keepends=True),
                fromfile=f"v1_freeze/{snapshot_name}",
                tofile=display_path,
                n=3,
            )
        )
    (REVIEW / "v2.diff").write_text("".join(chunks), encoding="utf-8", newline="\n")
    print(f"syntax_parse=PASS; diff={REVIEW / 'v2.diff'}; changed_files={len(PAIRS)}")


if __name__ == "__main__":
    main()
