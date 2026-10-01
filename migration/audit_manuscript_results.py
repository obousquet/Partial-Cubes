#!/usr/bin/env python3
"""Audit theorem-like LaTeX labels against database proof-source locators."""

from __future__ import annotations

import argparse
import json
import re
from collections import Counter, defaultdict
from pathlib import Path

from common import emit_json, load_json


RESULT_PREFIXES = {"thm", "prop", "cor", "lem", "ex", "conj", "prob"}
NON_RESULT_PREFIXES = {"app", "def", "eq", "fig", "item", "q", "rem", "sec", "tab"}
RESULT_ENVIRONMENTS = {
    "theorem",
    "maintheorem",
    "proposition",
    "corollary",
    "lemma",
    "example",
    "conjecture",
    "problem",
}


def iter_record_files(data_dir: Path):
    for path in sorted(data_dir.glob("*/*.json")):
        if path.name == "schema.json" or path.parent.name == "assets":
            continue
        yield path


def proof_sources(data_dir: Path) -> list[dict]:
    sources = []
    for path in iter_record_files(data_dir):
        try:
            record = json.loads(path.read_text())
        except json.JSONDecodeError:
            continue
        value = record.get("proof_source")
        if isinstance(value, str) and value:
            sources.append(
                {
                    "path": path.relative_to(data_dir.parent).as_posix(),
                    "proof_source": value,
                }
            )
    return sources


def workspace(path: str) -> str:
    parts = Path(path).parts
    if len(parts) >= 2 and parts[0] == "papers":
        return parts[1]
    return parts[0] if parts else "unknown"


def audit(label_index: dict, data_dir: Path) -> dict:
    grouped: dict[str, list[dict]] = defaultdict(list)
    for file_entry in label_index["files"]:
        for label in file_entry["labels"]:
            if (
                label["prefix"] in RESULT_PREFIXES
                or (
                    label["prefix"] not in NON_RESULT_PREFIXES
                    and (label.get("environment") or "").lower()
                    in RESULT_ENVIRONMENTS
                )
            ):
                grouped[label["label"]].append(
                    {
                        "path": file_entry["path"],
                        "line": label["line"],
                        "prefix": label["prefix"],
                        "environment": label.get("environment"),
                    }
                )

    sources = proof_sources(data_dir)
    claims = []
    workspace_counts: dict[str, Counter] = defaultdict(Counter)
    for label, occurrences in sorted(grouped.items()):
        pattern = re.compile(r"#" + re.escape(label) + r"(?![A-Za-z0-9:_.-])")
        covered_by = [
            source["path"]
            for source in sources
            if pattern.search(source["proof_source"])
        ]
        workspaces = sorted({workspace(item["path"]) for item in occurrences})
        status = "covered" if covered_by else "uncovered"
        for name in workspaces:
            workspace_counts[name]["total"] += 1
            workspace_counts[name][status] += 1
        claims.append(
            {
                "label": label,
                "prefix": occurrences[0]["prefix"],
                "environment": occurrences[0]["environment"],
                "status": status,
                "workspaces": workspaces,
                "occurrences": [
                    {"path": item["path"], "line": item["line"]}
                    for item in occurrences
                ],
                "covered_by": covered_by,
            }
        )

    status_counts = Counter(claim["status"] for claim in claims)
    return {
        "version": 1,
        "kind": "manuscript_result_audit",
        "scope": "Unique theorem-, proposition-, corollary-, lemma-, example-, conjecture-, and problem-like labels in the preserved current paper sources. Coverage requires an exact #label locator in a database proof_source field; open questions are excluded.",
        "counts": {
            "unique_result_labels": len(claims),
            "label_occurrences": sum(len(claim["occurrences"]) for claim in claims),
            "covered": status_counts["covered"],
            "uncovered": status_counts["uncovered"],
            "proof_source_records": len(sources),
        },
        "by_workspace": {
            name: dict(sorted(counts.items()))
            for name, counts in sorted(workspace_counts.items())
        },
        "claims": claims,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--label-index",
        type=Path,
        default=Path(__file__).with_name("generated") / "latex_labels.json",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "data",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name("manuscript_result_audit.json"),
    )
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    emit_json(
        audit(load_json(args.label_index.resolve()), args.data_dir.resolve()),
        args.output.resolve(),
        args.check,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
