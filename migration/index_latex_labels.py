#!/usr/bin/env python3
"""Index stable labels in paper TeX sources without interpreting claims."""

from __future__ import annotations

import argparse
import re
from collections import Counter
from pathlib import Path

from common import DEFAULT_CONFIG, SourceSnapshot, emit_json, load_config


TOKEN_RE = re.compile(
    r"\\begin\{(?P<begin>[^{}]+)\}|"
    r"\\end\{(?P<end>[^{}]+)\}|"
    r"\\(?P<section>part|chapter|section|subsection|subsubsection)\*?\{(?P<title>[^{}]*)\}|"
    r"\\label\s*\{(?P<label>[^{}]+)\}"
)


def index_file(path: Path, latex_root: Path, snapshot: SourceSnapshot) -> dict:
    text = snapshot.read_text(path)
    environments: list[str] = []
    current_section: dict[str, str] | None = None
    labels = []
    for match in TOKEN_RE.finditer(text):
        if match.group("begin"):
            environments.append(match.group("begin"))
        elif match.group("end"):
            environment = match.group("end")
            if environment in environments:
                reverse_index = environments[::-1].index(environment)
                del environments[len(environments) - reverse_index - 1 :]
        elif match.group("section"):
            current_section = {
                "level": match.group("section"),
                "title_tex": match.group("title"),
            }
        elif match.group("label"):
            label = match.group("label")
            labels.append(
                {
                    "label": label,
                    "prefix": label.split(":", 1)[0] if ":" in label else "unprefixed",
                    "line": text.count("\n", 0, match.start()) + 1,
                    "environment": environments[-1] if environments else None,
                    "section": current_section,
                }
            )
    return {
        "path": path.relative_to(latex_root).as_posix(),
        "sha256": snapshot.recorded_sha256(path),
        "labels": labels,
    }


def extract(latex_root: Path, config: dict, snapshot: SourceSnapshot) -> dict:
    if snapshot.manifest:
        files = [
            snapshot.repositories[entry["repository"]] / entry["path"]
            for entry in snapshot.manifest["files"]
            if "paper_tex" in entry["roles"] and Path(entry["path"]).suffix == ".tex"
        ]
    else:
        paper_root = latex_root / "papers"
        configured = set(config["paper_workspaces"])
        excluded = set(config.get("excluded_source_directories", []))
        files = [
            path
            for path in sorted(paper_root.glob("**/*.tex"))
            if path.relative_to(paper_root).parts[0] in configured
            and not excluded.intersection(path.relative_to(paper_root).parts[1:-1])
        ]
    indexed = [index_file(path, latex_root, snapshot) for path in sorted(files)]
    all_labels = [
        {"path": entry["path"], **label}
        for entry in indexed
        for label in entry["labels"]
    ]
    counts = Counter(label["prefix"] for label in all_labels)
    duplicate_map: dict[str, list[dict]] = {}
    for label in all_labels:
        duplicate_map.setdefault(label["label"], []).append(
            {"path": label["path"], "line": label["line"]}
        )
    duplicates = {
        key: occurrences
        for key, occurrences in sorted(duplicate_map.items())
        if len(occurrences) > 1
    }
    return {
        "version": 1,
        "kind": "latex_label_index",
        "scope": "Current TeX in configured paper workspaces; research/ and archive/ trees are excluded from admissible proof-source indexing.",
        "counts": {
            "tex_files": len(indexed),
            "labels": len(all_labels),
            "duplicate_label_names": len(duplicates),
            "labels_by_prefix": dict(sorted(counts.items())),
        },
        "duplicates": duplicates,
        "files": indexed,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--latex-root", type=Path)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name("generated") / "latex_labels.json",
    )
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    config = load_config(args.config)
    latex_root = (args.latex_root or Path(config["latex_repository"]).expanduser()).resolve()
    snapshot = SourceSnapshot(latex_root, config)
    emit_json(extract(latex_root, config, snapshot), args.output.resolve(), args.check)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
