#!/usr/bin/env python3
"""Validate human-owned migration queues without changing their decisions."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

from common import DEFAULT_CONFIG, load_config, load_json


def unique(records: list[dict], field: str, name: str, errors: list[str]) -> None:
    values = [record[field] for record in records]
    duplicates = sorted(value for value, count in Counter(values).items() if count > 1)
    if duplicates:
        errors.append(f"Duplicate {name} {field} values: {duplicates}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--migration-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    config = load_config(args.config)
    root = args.migration_dir
    errors: list[str] = []

    crosswalk = load_json(root / "class_crosswalk.json")
    occurrences = crosswalk["occurrences"]
    expected_occurrences = config["expected_counts"]["survey_nodes"] + config["expected_counts"]["minor_families"]
    if len(occurrences) != expected_occurrences:
        errors.append(f"Class crosswalk has {len(occurrences)} occurrences; expected {expected_occurrences}")
    unique(occurrences, "occurrence_id", "crosswalk", errors)
    allowed = set(crosswalk["allowed_dispositions"])
    invalid = sorted({record["disposition"] for record in occurrences} - allowed)
    if invalid:
        errors.append(f"Invalid crosswalk dispositions: {invalid}")
    for record in occurrences:
        if record["disposition"] != "unresolved" and not record.get("decision_basis"):
            errors.append(f"Resolved crosswalk occurrence lacks decision_basis: {record['occurrence_id']}")

    queue = load_json(root / "claim_queue.json")
    claims = queue["claims"]
    unique(claims, "candidate_id", "claim", errors)
    if len(claims) != queue["counts"]["candidates"]:
        errors.append("Claim queue count does not match claims array")
    actual_types = dict(Counter(claim["candidate_type"] for claim in claims))
    if actual_types != queue["counts"]["by_type"]:
        errors.append(f"Claim type counts are stale: {actual_types}")

    conflicts = load_json(root / "conflicts.json")
    unique(conflicts["items"], "id", "conflict", errors)
    invalid_conflicts = sorted(
        {item["disposition"] for item in conflicts["items"]}
        - set(conflicts["allowed_dispositions"])
    )
    if invalid_conflicts:
        errors.append(f"Invalid conflict dispositions: {invalid_conflicts}")

    consumers = load_json(root / "latex_consumers.json")
    papers = consumers["documents"]
    unique(papers, "paper", "LaTeX consumer", errors)
    configured = set(config["paper_workspaces"])
    inventoried = {paper["paper"] for paper in papers}
    if configured != inventoried:
        errors.append(
            f"LaTeX consumer inventory mismatch; missing={sorted(configured - inventoried)}, "
            f"extra={sorted(inventoried - configured)}"
        )

    manifest = load_json(root / "source_manifest.json")
    if manifest["missing_references"]:
        errors.append(f"Source manifest has {len(manifest['missing_references'])} missing references")

    if errors:
        print("Migration review queue validation failed:")
        for error in errors:
            print(f"- {error}")
        return 1
    print(
        f"Validated {len(occurrences)} class occurrences, {len(claims)} claim candidates, "
        f"{len(conflicts['items'])} conflicts, and {len(papers)} LaTeX documents."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
