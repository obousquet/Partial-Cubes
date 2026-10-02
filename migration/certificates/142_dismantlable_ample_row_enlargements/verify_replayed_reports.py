#!/usr/bin/env python3
"""Compare local row-enlargement replays with the preserved source reports."""

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
RESEARCH = ROOT / "papers" / "dismantlable-ample" / "research"
NAMES = (
    "one_row_enlargement",
    "one_row_enlargement_verification",
    "simultaneous_row_enlargement",
    "two_row_fiber_shrink",
    "two_row_verification",
    "two_row_fiber_positivity",
    "two_row_corner_compatibility",
)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mathematical_payload(name: str, value: dict) -> dict:
    value = dict(value)
    value.pop("seconds", None)
    if name == "one_row_enlargement_verification":
        # The verifier hashes its generated input report; that report differs
        # from the source copy only in its measured run time.
        value.pop("source_sha256", None)
    return value


def main() -> None:
    rows = []
    for name in NAMES:
        replay_path = RESEARCH / f"{name}.json"
        source_path = RESEARCH / f"{name}_source.json"
        replay = json.loads(replay_path.read_text())
        source = json.loads(source_path.read_text())
        assert mathematical_payload(name, replay) == mathematical_payload(name, source)
        rows.append(
            {
                "name": name,
                "replay_sha256": sha(replay_path),
                "source_sha256": sha(source_path),
                "comparison": "equal after removing measured seconds"
                + (
                    " and the derived input-report digest"
                    if name == "one_row_enlargement_verification"
                    else ""
                ),
            }
        )
    report = {
        "status": "PASS",
        "reports": rows,
        "scope": "All mathematical fields, explicit families, orders, certificates, counts, and source-script digests agree with the preserved source reports.",
    }
    target = ROOT / "replayed_report_reconciliation.json"
    target.write_text(json.dumps(report, indent=2) + "\n")
    print(f"PASS: reconciled {len(rows)} source reports")


if __name__ == "__main__":
    main()
