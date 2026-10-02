#!/usr/bin/env python3
"""Audit the LaTeX-owned fixed-contraction archive and local certificates."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
MANIFEST = HERE / "fixed_contraction_archive_manifest.json"
DEFAULT_OUTPUT = HERE / "certificate_archive_audit.json"


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            value.update(block)
    return value.hexdigest()


def audit_fixed_contraction(latex_root: Path) -> dict:
    manifest = json.loads(MANIFEST.read_text())
    archive = (
        latex_root
        / "papers/dismantlable-ample/research/two_row_fixed_contraction"
    )
    total = 0
    for record in manifest["files"]:
        path = archive / record["path"]
        assert path.stat().st_size == record["bytes"], path
        assert digest(path) == record["sha256"], path
        total += record["bytes"]
    verification = json.loads((archive / "verification.json").read_text())
    assert verification["status"] == "PASS"
    assert len(verification["pairs"]) == 6
    assert all(item["models"] == 2 for item in verification["pairs"])
    return {
        "archive_file_count": len(manifest["files"]),
        "archive_bytes": total,
        "all_archive_hashes_match": True,
        "presentations": 6,
        "models_per_presentation": 2,
        "status": "all complete exclusions independently DRAT-checked",
    }


def audit_local_certificates() -> dict:
    block = HERE / "papers/dismantlable-ample/research"
    source = block / "block_row_completion_source.json"
    replay = block / "block_row_completion.json"
    assert replay.read_bytes() == source.read_bytes()
    block_report = json.loads(replay.read_text())
    assert block_report["order"] == 684
    assert block_report["unique_corner"] == 2167
    assert block_report["minimum_degree"] == 3
    assert len(block_report["column_classes"]) == 13
    assert len(block_report["negative_overlap"]) == 267

    overlap = HERE / "papers/ample-corners/evidence/two_row_overlap_envelope"
    overlap_verification = json.loads((overlap / "verification.json").read_text())
    assert overlap_verification["status"] == "PASS"
    assert overlap_verification["base_order"] == 48
    assert overlap_verification["envelope_order"] == 96
    assert overlap_verification["formula_reconstructed"]
    assert overlap_verification["proof_replayed"]
    overlap_report = json.loads((overlap / "report.json").read_text())
    assert digest(overlap / "search.cnf") == overlap_report["cnf_sha256"]
    assert digest(overlap / "search.drat.gz") == "bb436bd511ea54af1c09b8e8f14bbb2f7c30c62f8d0fc76ad8952ccd8903e408"
    return {
        "block_row_completion": {
            "replay_matches_source_byte_for_byte": True,
            "order": 684,
            "unique_corner": 2167,
            "minimum_degree": 3,
            "signed_column_classes": 13,
            "negative_overlap_order": 267,
        },
        "two_row_overlap_envelope": {
            "base_order": 48,
            "envelope_order": 96,
            "variables": overlap_report["variables"],
            "clauses": overlap_report["clauses"],
            "formula_reconstructed": True,
            "proof_replayed": True,
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--latex-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    report = {
        "schema": "partial-cubes-block-row-and-two-row-rigidity-audit-v1",
        "fixed_contraction": audit_fixed_contraction(args.latex_root),
        "local_certificates": audit_local_certificates(),
        "ownership": "The database owns the complete block-row replay, complete overlap-envelope certificate, compact fixed-contraction reports, the archive manifest, and this audit. The 23 MB fixed-contraction CNF, DRAT, and log archive remains LaTeX-owned.",
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(
        "PASS: audited "
        f"{report['fixed_contraction']['archive_file_count']} fixed-contraction "
        "files and both local certificates"
    )


if __name__ == "__main__":
    main()
