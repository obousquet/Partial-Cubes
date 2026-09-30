#!/usr/bin/env python3
"""Extract minor-atlas candidates and verify its embedded source fingerprints."""

from __future__ import annotations

import argparse
from pathlib import Path

from common import DEFAULT_CONFIG, emit_json, load_config, load_json, resolve_source_path, sha256_file


def fingerprint(path: Path, expected: str | None) -> dict:
    if not path.is_file():
        return {
            "exists": False,
            "expected_sha256": expected,
            "actual_sha256": None,
            "matches": False,
        }
    actual = sha256_file(path)
    return {
        "exists": True,
        "expected_sha256": expected,
        "actual_sha256": actual,
        "matches": expected == actual if expected else None,
    }


def extract(latex_root: Path, config: dict) -> dict:
    inventory_rel = config["structured_sources"]["minor_inventory"]
    inventory_path = latex_root / inventory_rel
    raw = load_json(inventory_path)
    expected = config["expected_counts"]
    for key, expected_key in (
        ("families", "minor_families"),
        ("relations", "minor_relations"),
        ("sources", "minor_sources"),
        ("questions", "minor_questions"),
    ):
        if len(raw[key]) != expected[expected_key]:
            raise SystemExit(f"Expected {expected[expected_key]} {key}, found {len(raw[key])}")

    families = []
    for index, record in enumerate(raw["families"]):
        candidate = dict(record)
        candidate.update(
            {
                "source_locator": f"PartialCubes:{inventory_rel}#/families/{index}",
                "migration_state": "candidate",
            }
        )
        families.append(candidate)

    relations = []
    for index, record in enumerate(raw["relations"], start=1):
        candidate = dict(record)
        candidate.update(
            {
                "candidate_id": f"minor-relation-{index:03d}",
                "source_locator": f"PartialCubes:{inventory_rel}#/relations/{index - 1}",
                "migration_state": "candidate",
            }
        )
        relations.append(candidate)

    sources = []
    for index, record in enumerate(raw["sources"]):
        candidate = dict(record)
        source_path = resolve_source_path(latex_root, record["path"])
        candidate["fingerprint_check"] = fingerprint(source_path, record.get("sha256"))
        label_checks = {}
        for label, label_source in record.get("label_sources", {}).items():
            label_path = resolve_source_path(latex_root, label_source["path"])
            label_checks[label] = {
                "path": label_source["path"],
                **fingerprint(label_path, label_source.get("sha256")),
            }
        candidate["label_fingerprint_checks"] = label_checks
        candidate["source_locator"] = f"PartialCubes:{inventory_rel}#/sources/{index}"
        sources.append(candidate)

    questions = []
    for index, record in enumerate(raw["questions"]):
        candidate = dict(record)
        candidate.update(
            {
                "source_locator": f"PartialCubes:{inventory_rel}#/questions/{index}",
                "migration_disposition": "excluded_broad_question",
                "disposition_reason": "Broad research questions remain prose in the LaTeX workspace.",
            }
        )
        questions.append(candidate)

    stale_sources = [
        source["id"]
        for source in sources
        if not source["fingerprint_check"]["matches"]
        or any(not check["matches"] for check in source["label_fingerprint_checks"].values())
    ]

    return {
        "version": 1,
        "kind": "minor_inventory_candidates",
        "source": {
            "path": inventory_rel,
            "sha256": sha256_file(inventory_path),
            "title": raw.get("title"),
            "date": raw.get("date"),
            "scope": raw.get("scope"),
            "conventions": raw.get("conventions"),
        },
        "counts": {
            "family_candidates": len(families),
            "relation_candidates": len(relations),
            "source_packets": len(sources),
            "excluded_questions": len(questions),
            "stale_source_packets": len(stale_sources),
        },
        "stale_source_packet_ids": sorted(stale_sources),
        "families": families,
        "relations": relations,
        "sources": sources,
        "questions": questions,
        "boundary_study": raw.get("boundary_study"),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--latex-root", type=Path)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name("generated") / "minor_inventory.json",
    )
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    config = load_config(args.config)
    latex_root = (args.latex_root or Path(config["latex_repository"]).expanduser()).resolve()
    emit_json(extract(latex_root, config), args.output.resolve(), args.check)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
