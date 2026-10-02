#!/usr/bin/env python3
"""Audit the LaTeX-owned Q7 unique-corner DRAT certificate archives.

The database preserves the manifests and generator sources in supplemental
snapshot 2026-10-02-b003.  The large CNF, DRAT, and checker-log files remain
owned by the PartialCubes LaTeX repository; this audit checks every archived
file against the digest recorded in its preserved manifest.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
SNAPSHOT_ROOT = (
    REPO_ROOT
    / "migration/snapshots/2026-10-02-b003/PartialCubes/papers/ample-corners/evidence"
)
DEFAULT_OUTPUT = Path(__file__).with_name("unique_corner_q7_certificate_audit.json")
EXPECTED = {
    4: {
        "directory": "unique_corner_q7_degree4_certificates",
        "manifest_sha256": "6924119d4eb46ea27dceaefa21a3cdfea70f1184d47ccd766190b9f7a78b8b62",
        "case_count": 18,
        "sufficient_case_count": 17,
        "variables": 10943,
        "clauses": 65378,
    },
    3: {
        "directory": "unique_corner_q7_degree3_certificates",
        "manifest_sha256": "f6a88849d29952704d81b5fb2f33c51e8950e8cc7abb05ec2e41309be71bdb78",
        "case_count": 185,
        "sufficient_case_count": 185,
        "variables": 10943,
        "clauses": 65379,
    },
}


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def audit_degree(degree: int, latex_root: Path) -> dict:
    expected = EXPECTED[degree]
    directory = expected["directory"]
    compact = SNAPSHOT_ROOT / directory
    archive = latex_root / "papers/ample-corners/evidence" / directory
    manifest_path = compact / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    assert digest(manifest_path) == expected["manifest_sha256"]
    assert len(manifest["cases"]) == expected["case_count"]
    assert all(case["status"] == "UNSAT-DRAT-checked" for case in manifest["cases"])
    assert all(case["checker_returncode"] == 0 for case in manifest["cases"])
    assert all(case["variables"] == expected["variables"] for case in manifest["cases"])
    assert all(case["clauses"] == expected["clauses"] for case in manifest["cases"])

    compact_source_hashes = {}
    for name, expected_hash in manifest["dependencies"].items():
        source = compact / "sources" / name
        actual = digest(source)
        assert actual == expected_hash, (source, expected_hash, actual)
        compact_source_hashes[name] = actual

    archive_bytes = 0
    checked_files = 0
    for case in manifest["cases"]:
        representative = case["representative"]
        checks = [
            (archive / f"{representative}.cnf.gz", case["cnf_archive_sha256"]),
            (archive / f"{representative}.drat.gz", case["drat_archive_sha256"]),
            (archive / f"{representative}.check.log", case["check_log_sha256"]),
        ]
        for path, expected_hash in checks:
            actual = digest(path)
            assert actual == expected_hash, (path, expected_hash, actual)
            archive_bytes += path.stat().st_size
            checked_files += 1

    return {
        "degree": degree,
        "case_count": len(manifest["cases"]),
        "sufficient_case_count": expected["sufficient_case_count"],
        "representatives": [case["representative"] for case in manifest["cases"]],
        "status": "all UNSAT-DRAT-checked",
        "variables_per_formula": expected["variables"],
        "clauses_per_formula": expected["clauses"],
        "checker_sha256": manifest["checker_sha256"],
        "generator_sha256": manifest["script_sha256"],
        "dependency_sha256": compact_source_hashes,
        "manifest_sha256": expected["manifest_sha256"],
        "archive_file_count": checked_files,
        "archive_bytes": archive_bytes,
        "all_archive_hashes_match": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--latex-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    reports = [audit_degree(degree, args.latex_root) for degree in (4, 3)]
    output = {
        "schema": "partial-cubes-q7-unique-corner-certificate-audit-v1",
        "scope": "Full digest audit of the LaTeX-owned seven-coordinate unique-corner DRAT archives against database-preserved manifests and generator sources.",
        "degrees": reports,
        "totals": {
            "case_count": sum(report["case_count"] for report in reports),
            "sufficient_case_count": sum(
                report["sufficient_case_count"] for report in reports
            ),
            "archive_file_count": sum(report["archive_file_count"] for report in reports),
            "archive_bytes": sum(report["archive_bytes"] for report in reports),
        },
        "ownership": "CNF, DRAT, and checker-log bytes remain LaTeX-owned; the database owns this audit, exact manifests, source snapshots, and normalized theorem records.",
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(
        f"audited {output['totals']['case_count']} cases and "
        f"{output['totals']['archive_file_count']} archived files"
    )


if __name__ == "__main__":
    main()
