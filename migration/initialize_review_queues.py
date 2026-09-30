#!/usr/bin/env python3
"""Initialize human-owned crosswalk and claim queues from generated candidates."""

from __future__ import annotations

import argparse
from pathlib import Path

from common import canonical_json, load_json


EXACT_OVERLAPS = {"ample", "daisy", "learning", "maximum", "median", "treelike"}
ALIAS_OVERLAPS = {
    "damp": ("dismantlable_ample", "peelable"),
    "peelable": ("dismantlable_ample", "damp"),
    "rper": ("cover_recursive_peripheral", "cover"),
    "cover": ("cover_recursive_peripheral", "rper"),
    "realample": ("realizable_ample", "realizable"),
    "realizable": ("realizable_ample", "realample"),
}


def overlap_hint(source_id: str) -> dict | None:
    if source_id in EXACT_OVERLAPS:
        return {"group": source_id, "kind": "exact_source_id", "counterpart_source_id": source_id}
    if source_id in ALIAS_OVERLAPS:
        group, counterpart = ALIAS_OVERLAPS[source_id]
        return {"group": group, "kind": "possible_alias", "counterpart_source_id": counterpart}
    return None


def build_crosswalk(survey: dict, minor: dict) -> dict:
    occurrences = []
    for record in survey["classes"]:
        occurrences.append(
            {
                "occurrence_id": f"survey:{record['source_id']}",
                "source": "survey_atlas",
                "source_id": record["source_id"],
                "display_tex": record["label_tex"],
                "style": record["style"],
                "source_locator": record["source_locator"],
                "overlap_hint": overlap_hint(record["source_id"]),
                "disposition": "unresolved",
                "canonical_short_name": None,
                "decision_basis": None,
            }
        )
    for record in minor["families"]:
        occurrences.append(
            {
                "occurrence_id": f"minor:{record['id']}",
                "source": "minor_inventory",
                "source_id": record["id"],
                "display_name": record["name"],
                "symbol": record["symbol"],
                "source_locator": record["source_locator"],
                "overlap_hint": overlap_hint(record["id"]),
                "disposition": "unresolved",
                "canonical_short_name": None,
                "decision_basis": None,
            }
        )
    return {
        "version": 1,
        "kind": "human_reviewed_class_crosswalk",
        "warning": "Do not regenerate over human decisions. Extraction hints are not equivalence decisions.",
        "allowed_dispositions": [
            "canonical",
            "same_as",
            "parameterized_member",
            "related_not_equal",
            "presentation_only",
            "unresolved",
            "excluded",
        ],
        "counts": {"occurrences": len(occurrences), "unresolved": len(occurrences)},
        "occurrences": occurrences,
    }


def build_claim_queue(survey: dict, minor: dict) -> dict:
    claims = []
    for record in survey["containments"]:
        claims.append(
            {
                "candidate_id": record["candidate_id"],
                "candidate_type": "class_relation",
                "raw_relation_type": "inclusion",
                "subject_source_id": record["subclass_source_id"],
                "object_source_id": record["superclass_source_id"],
                "source_locator": record["source_locator"],
                "review_state": "unreviewed",
                "disposition": None,
            }
        )
    for record in minor["relations"]:
        claims.append(
            {
                "candidate_id": record["candidate_id"],
                "candidate_type": "class_relation",
                "raw_relation_type": record.get("relation", "inclusion"),
                "subject_source_id": record["source"],
                "object_source_id": record["target"],
                "raw_status": record["status"],
                "raw_strict": record.get("strict"),
                "reference_source_id": record["reference"],
                "source_locator": record["source_locator"],
                "review_state": "unreviewed",
                "disposition": None,
            }
        )
    for family in minor["families"]:
        for operation in ("p", "c", "pc"):
            claims.append(
                {
                    "candidate_id": f"minor-closure-{family['id']}-{operation}",
                    "candidate_type": "operation_result",
                    "class_source_id": family["id"],
                    "operation_source_id": operation,
                    "raw_value": family[operation],
                    "source_ids": family["sources"],
                    "source_locator": family["source_locator"],
                    "review_state": "unreviewed",
                    "disposition": None,
                }
            )
        claims.append(
            {
                "candidate_id": f"minor-obstruction-summary-{family['id']}",
                "candidate_type": "obstruction_summary",
                "class_source_id": family["id"],
                "raw_statement": family["basis"],
                "source_ids": family["sources"],
                "source_locator": family["source_locator"],
                "review_state": "requires_atomic_decomposition",
                "disposition": None,
            }
        )
        for closure_kind in ("c_closure", "pc_closure"):
            claims.append(
                {
                    "candidate_id": f"minor-{closure_kind.replace('_', '-')}-{family['id']}",
                    "candidate_type": "closure_description",
                    "class_source_id": family["id"],
                    "closure_kind": closure_kind,
                    "raw_statement": family[closure_kind],
                    "source_ids": family["sources"],
                    "source_locator": family["source_locator"],
                    "review_state": "requires_atomic_decomposition",
                    "disposition": None,
                }
            )
    type_counts: dict[str, int] = {}
    for claim in claims:
        type_counts[claim["candidate_type"]] = type_counts.get(claim["candidate_type"], 0) + 1
    return {
        "version": 1,
        "kind": "human_reviewed_claim_queue",
        "warning": "Do not regenerate over human decisions. Raw status and diagram assertions are not proof review.",
        "counts": {"candidates": len(claims), "by_type": type_counts},
        "claims": claims,
    }


def write_new(path: Path, data: dict) -> None:
    if path.exists():
        raise SystemExit(f"Refusing to overwrite human-owned queue: {path}")
    path.write_text(canonical_json(data), encoding="utf-8")
    print(f"Initialized: {path}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--generated-dir", type=Path, default=Path(__file__).with_name("generated"))
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    survey = load_json(args.generated_dir / "survey_atlas.json")
    minor = load_json(args.generated_dir / "minor_inventory.json")
    write_new(args.output_dir / "class_crosswalk.json", build_crosswalk(survey, minor))
    write_new(args.output_dir / "claim_queue.json", build_claim_queue(survey, minor))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
