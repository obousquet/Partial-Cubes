#!/usr/bin/env python3
"""Extract survey class and containment candidates without promoting records."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

from common import DEFAULT_CONFIG, SourceSnapshot, emit_json, line_number, load_config


NODE_RE = re.compile(
    r"\\node\[(?P<options>[^\]]+)\]\s*"
    r"\((?P<id>[^)]+)\)\s*at\s*\((?P<position>[^)]+)\)\s*"
    r"\{(?P<label>.*?)\};",
    re.DOTALL,
)


def compact_tex(text: str) -> str:
    return re.sub(r"\s+", " ", text.strip())


def extract(latex_root: Path, config: dict, snapshot: SourceSnapshot) -> dict:
    sources = config["structured_sources"]
    atlas_rel = sources["survey_atlas"]
    containments_rel = sources["survey_containments"]
    atlas_path = latex_root / atlas_rel
    containments_path = latex_root / containments_rel

    source_data = json.loads(snapshot.read_text(containments_path))
    containments = source_data["containments"]
    endpoints = {
        record[key]
        for record in containments
        for key in ("subclass", "superclass")
    }

    atlas_text = snapshot.read_text(atlas_path)
    parsed_nodes: dict[str, dict] = {}
    presentation_nodes: list[str] = []
    for match in NODE_RE.finditer(atlas_text):
        node_id = match.group("id").strip()
        options = [part.strip() for part in match.group("options").split(",")]
        if node_id not in endpoints:
            presentation_nodes.append(node_id)
            continue
        if node_id in parsed_nodes:
            raise SystemExit(f"Duplicate survey node: {node_id}")
        coordinates = [part.strip() for part in match.group("position").split(",")]
        if len(coordinates) != 2:
            raise SystemExit(f"Cannot parse coordinates for survey node {node_id}")
        parsed_nodes[node_id] = {
            "source_id": node_id,
            "style": options[0],
            "tikz_options": options,
            "position": {"x": coordinates[0], "y": coordinates[1]},
            "label_tex": compact_tex(match.group("label")),
            "source_locator": f"PartialCubes:{atlas_rel}:{line_number(atlas_text, match.start())}",
            "migration_state": "candidate",
        }

    missing = sorted(endpoints - parsed_nodes.keys())
    extra = sorted(parsed_nodes.keys() - endpoints)
    expected = config["expected_counts"]
    if missing or extra:
        raise SystemExit(f"Survey node mismatch; missing={missing}, extra={extra}")
    if len(parsed_nodes) != expected["survey_nodes"]:
        raise SystemExit(f"Expected {expected['survey_nodes']} survey nodes, found {len(parsed_nodes)}")
    if len(containments) != expected["survey_containments"]:
        raise SystemExit(
            f"Expected {expected['survey_containments']} survey containments, found {len(containments)}"
        )

    relation_candidates = []
    for index, record in enumerate(containments, start=1):
        relation_candidates.append(
            {
                "candidate_id": f"survey-containment-{index:03d}",
                "subclass_source_id": record["subclass"],
                "superclass_source_id": record["superclass"],
                "source_locator": f"PartialCubes:{containments_rel}#/containments/{index - 1}",
                "migration_state": "candidate",
            }
        )

    return {
        "version": 1,
        "kind": "survey_atlas_candidates",
        "source": {
            "atlas": atlas_rel,
            "atlas_sha256": snapshot.recorded_sha256(atlas_path),
            "containments": containments_rel,
            "containments_sha256": snapshot.recorded_sha256(containments_path),
            "schema_version": source_data.get("schema_version"),
            "description": source_data.get("description"),
        },
        "counts": {
            "class_candidates": len(parsed_nodes),
            "containment_candidates": len(relation_candidates),
            "presentation_only_nodes": len(presentation_nodes),
        },
        "presentation_only_nodes": sorted(presentation_nodes),
        "classes": [parsed_nodes[key] for key in sorted(parsed_nodes)],
        "containments": relation_candidates,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--latex-root", type=Path)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name("generated") / "survey_atlas.json",
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
