#!/usr/bin/env python3
"""Reconcile local two-row replays with the preserved source reports."""

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RESEARCH = ROOT / "papers" / "dismantlable-ample" / "research"
NAMES = (
    "two_row_relative_extension",
    "two_row_both_nonrelative",
    "two_row_parent_minor_filter",
)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    rows = []
    for name in NAMES:
        replay_path = RESEARCH / f"{name}.json"
        source_path = RESEARCH / f"{name}_source.json"
        replay = json.loads(replay_path.read_text())
        source = json.loads(source_path.read_text())
        replay.pop("seconds", None)
        source.pop("seconds", None)
        assert replay == source
        rows.append(
            {
                "name": name,
                "replay_sha256": sha(replay_path),
                "source_sha256": sha(source_path),
                "comparison": "all fields agree after removing measured seconds",
            }
        )
    report = {
        "status": "PASS",
        "reports": rows,
        "scope": "All mathematical fields, explicit families, orders, minor certificates, source digests, and exclusion records agree.",
    }
    (ROOT / "replayed_report_reconciliation.json").write_text(
        json.dumps(report, indent=2) + "\n"
    )
    print(f"PASS: reconciled {len(rows)} source reports")


if __name__ == "__main__":
    main()
