#!/usr/bin/env python3
"""Validate candidate reconciliation totals and write a migration coverage report."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

from common import DEFAULT_CONFIG, emit_json, load_config, load_json


TABLES = (
    "classes",
    "characterizations",
    "characterization_equivalences",
    "relations",
    "operations",
    "operation_results",
    "invariants",
    "invariant_values",
    "obstruction_families",
    "obstruction_bases",
    "examples",
    "results",
)


def audit(config: dict, generated_dir: Path, data_dir: Path) -> dict:
    migration_dir = generated_dir.parent
    survey = load_json(generated_dir / "survey_atlas.json")
    minor = load_json(generated_dir / "minor_inventory.json")
    labels = load_json(generated_dir / "latex_labels.json")
    expected = config["expected_counts"]
    actual = {
        "survey_nodes": survey["counts"]["class_candidates"],
        "survey_containments": survey["counts"]["containment_candidates"],
        "minor_families": minor["counts"]["family_candidates"],
        "minor_relations": minor["counts"]["relation_candidates"],
        "minor_sources": minor["counts"]["source_packets"],
        "minor_questions": minor["counts"]["excluded_questions"],
    }
    errors = [
        f"{key}: expected {value}, found {actual[key]}"
        for key, value in expected.items()
        if actual.get(key) != value
    ]
    if errors:
        raise SystemExit("Coverage gate failed:\n" + "\n".join(errors))

    relation_statuses = Counter(record["status"] for record in minor["relations"])
    relation_types = Counter(record.get("relation", "containment") for record in minor["relations"])
    database_counts = {
        table: len(list((data_dir / table).glob("[0-9]*.json")))
        for table in TABLES
    }
    crosswalk_path = migration_dir / "class_crosswalk.json"
    canonical_path = migration_dir / "canonical_classes.json"
    crosswalk = load_json(crosswalk_path) if crosswalk_path.is_file() else None
    canonical = load_json(canonical_path) if canonical_path.is_file() else None
    disposition_counts = (
        Counter(record["disposition"] for record in crosswalk["occurrences"])
        if crosswalk
        else Counter({"unresolved": actual["survey_nodes"] + actual["minor_families"]})
    )
    unresolved = disposition_counts.get("unresolved", 0)
    return {
        "version": 1,
        "kind": "migration_coverage",
        "phase": 1 if canonical and unresolved == 0 else 0,
        "expected_seed_counts": expected,
        "actual_seed_counts": actual,
        "candidate_totals": {
            "class_occurrences": actual["survey_nodes"] + actual["minor_families"],
            "relation_occurrences": actual["survey_containments"] + actual["minor_relations"],
            "closure_cells": actual["minor_families"] * 3,
            "latex_labels": labels["counts"]["labels"],
        },
        "minor_relation_statuses": dict(sorted(relation_statuses.items())),
        "minor_relation_types": dict(sorted(relation_types.items())),
        "source_packet_fingerprints": {
            "matching_inventory_fingerprint": actual["minor_sources"] - minor["counts"]["stale_source_packets"],
            "mismatching_inventory_fingerprint": minor["counts"]["stale_source_packets"],
            "stale_ids": minor["stale_source_packet_ids"],
        },
        "dispositions": {
            "candidate_class_occurrences": actual["survey_nodes"] + actual["minor_families"],
            "class_occurrences_by_disposition": dict(sorted(disposition_counts.items())),
            "reviewed_canonical_identities": canonical["count"] if canonical else 0,
            "promoted_classes": database_counts["classes"],
            "excluded_broad_questions": actual["minor_questions"],
            "unresolved_class_occurrences": unresolved,
        },
        "database_record_counts": database_counts,
        "gate": "passed",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--generated-dir", type=Path, default=Path(__file__).with_name("generated"))
    parser.add_argument("--data-dir", type=Path, default=Path(__file__).resolve().parents[1] / "data")
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("coverage.json"))
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    emit_json(
        audit(load_config(args.config), args.generated_dir.resolve(), args.data_dir.resolve()),
        args.output.resolve(),
        args.check,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
