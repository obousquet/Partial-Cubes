#!/usr/bin/env python3
"""Validate human-owned migration queues without changing their decisions."""

from __future__ import annotations

import argparse
import re
from collections import Counter
from pathlib import Path

from common import DEFAULT_CONFIG, load_config, load_json, sha256_file


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

    canonical = load_json(root / "canonical_classes.json")
    classes = canonical["classes"]
    if len(classes) != canonical["count"]:
        errors.append("Canonical class count does not match classes array")
    unique(classes, "allocated_id", "canonical class", errors)
    unique(classes, "short_name", "canonical class", errors)
    allocated = sorted(record["allocated_id"] for record in classes)
    if allocated != list(range(1, len(classes) + 1)):
        errors.append("Canonical class IDs are not contiguous from 1")
    invalid_short_names = sorted(
        record["short_name"]
        for record in classes
        if not re.fullmatch(r"[a-z][a-z0-9_]*", record["short_name"])
    )
    if invalid_short_names:
        errors.append(f"Invalid canonical short names: {invalid_short_names}")
    canonical_names = {record["short_name"] for record in classes}
    missing_identities = sorted(
        {record["canonical_short_name"] for record in occurrences} - canonical_names
    )
    if missing_identities:
        errors.append(f"Crosswalk references unallocated identities: {missing_identities}")
    if crosswalk["counts"].get("unresolved") != 0:
        errors.append(f"Crosswalk still has {crosswalk['counts'].get('unresolved')} unresolved occurrences")

    identity_decisions = load_json(root / "identity_decisions.json")
    if len(identity_decisions["overlaps"]) != 9:
        errors.append("Expected nine reviewed overlap decisions")

    data_class_dir = root.parent / "data" / "classes"
    database_classes = {
        record["short_name"]: record
        for path in data_class_dir.glob("[0-9]*.json")
        for record in [load_json(path)]
    }
    allocated_classes = {record["short_name"]: record for record in classes}
    unexpected_database_classes = sorted(set(database_classes) - set(allocated_classes))
    if unexpected_database_classes:
        errors.append(f"Database classes lack reserved identities: {unexpected_database_classes}")
    for short_name, database_record in database_classes.items():
        allocation = allocated_classes[short_name]
        if allocation["allocated_id"] != database_record["id"]:
            errors.append(f"Allocated ID mismatch for database class {short_name}")
        if allocation.get("promotion_state") != "database":
            errors.append(f"Database class is not marked promoted in allocation: {short_name}")
        if allocation.get("database_record") != f"#classes/{short_name}":
            errors.append(f"Database pointer is missing from allocation: {short_name}")

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
    excluded_directories = set(config.get("excluded_source_directories", []))
    inadmissible_paths = sorted(
        f"{entry['repository']}:{entry['path']}"
        for entry in manifest["files"]
        if excluded_directories.intersection(Path(entry["path"]).parts)
    )
    if inadmissible_paths:
        errors.append(f"Source manifest includes excluded research/archive paths: {inadmissible_paths}")
    archive = load_json(root / "source_archive.json")
    archived_keys = {(entry["repository"], entry["path"]) for entry in archive["files"]}
    required_archive_keys = {
        (entry["repository"], entry["path"])
        for entry in manifest["files"]
        if entry["git_state"] != "tracked"
    }
    if archived_keys != required_archive_keys:
        errors.append(
            f"Dirty source archive mismatch; missing={sorted(required_archive_keys - archived_keys)}, "
            f"extra={sorted(archived_keys - required_archive_keys)}"
        )
    for entry in archive["files"]:
        archive_path = root / entry["archive_path"]
        if not archive_path.is_file() or sha256_file(archive_path) != entry["sha256"]:
            errors.append(f"Dirty source archive is missing or corrupt: {entry['archive_path']}")

    if errors:
        print("Migration review queue validation failed:")
        for error in errors:
            print(f"- {error}")
        return 1
    print(
        f"Validated {len(occurrences)} class occurrences into {len(classes)} identities, "
        f"{len(claims)} claim candidates, "
        f"{len(conflicts['items'])} conflicts, and {len(papers)} LaTeX documents."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
