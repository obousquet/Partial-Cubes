#!/usr/bin/env python3
"""Check hierarchy reduction, visual semantics, and generated graph edges."""

from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from pathlib import Path


def has_path(source: str, target: str, adjacency: dict[str, set[str]]) -> bool:
    pending = [source]
    seen = {source}
    while pending:
        current = pending.pop()
        for successor in adjacency.get(current, set()):
            if successor == target:
                return True
            if successor not in seen:
                seen.add(successor)
                pending.append(successor)
    return False


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--site-dir", type=Path, default=Path("docs"))
    parser.add_argument("--math-database-dir", type=Path, default=Path("../math_database"))
    parser.add_argument(
        "--minor-atlas-inventory",
        type=Path,
        default=Path("migration/generated/minor_inventory.json"),
    )
    parser.add_argument(
        "--minor-atlas-graph-map",
        type=Path,
        default=Path("migration/minor_atlas_graph_classes.json"),
    )
    args = parser.parse_args()

    data_dir = args.data_dir.resolve()
    site_dir = args.site_dir.resolve()
    sys.path.insert(0, str(args.math_database_dir.resolve()))
    sys.path.insert(0, str(data_dir))
    import load_utils  # type: ignore[import-not-found]
    import make_graph  # type: ignore[import-not-found]

    cache = load_utils.get_table_entries_cache(str(data_dir))
    relations = cache.get_table_entries("relations")
    classes = cache.get_table_entries("classes")
    class_refs = {f'#classes/{entry["short_name"]}' for entry in classes}
    aliases = {
        reference: f'#classes/{entry["short_name"]}'
        for entry in classes
        for reference in (f'#classes/{entry["short_name"]}', f'#classes/{entry["id"]}')
    }
    established = [
        relation
        for relation in relations
        if relation.get("status") == "established"
        and relation.get("relation_type") in {"inclusion", "strict_inclusion"}
    ]

    hierarchy = make_graph.generate_class_hierarchy(cache)
    rendered_refs = {
        edge["ref"]: edge
        for edge in hierarchy["edges"]
        if edge.get("ref", "").startswith("#relations/")
    }
    cover_pairs: dict[str, tuple[str, str]] = {}
    adjacency: dict[str, set[str]] = defaultdict(set)
    for relation in established:
        ref = f'#relations/{relation["short_name"]}'
        if ref not in rendered_refs:
            continue
        subject = aliases.get(relation["subject_id"])
        object_ = aliases.get(relation["object_id"])
        if subject and object_:
            cover_pairs[ref] = (subject, object_)
            adjacency[subject].add(object_)

    errors: list[str] = []
    for relation in established:
        subject = aliases.get(relation["subject_id"])
        object_ = aliases.get(relation["object_id"])
        if not subject or not object_:
            errors.append(f"unresolved established inclusion: {relation['short_name']}")
        elif not has_path(subject, object_, adjacency):
            errors.append(f"established inclusion is absent from the Hasse reachability: {relation['short_name']}")

    for ref, (subject, object_) in cover_pairs.items():
        reduced = {node: set(successors) for node, successors in adjacency.items()}
        reduced[subject].discard(object_)
        if has_path(subject, object_, reduced):
            errors.append(f"redundant inclusion was rendered: {ref}")

    by_short_name = {relation["short_name"]: relation for relation in relations}
    for ref, edge in rendered_refs.items():
        relation = by_short_name.get(ref.removeprefix("#relations/"))
        if not relation:
            errors.append(f"rendered relation does not exist: {ref}")
            continue
        if relation["relation_type"] == "strict_inclusion":
            if edge.get("color") != make_graph.STRICT_INCLUSION_COLOR or edge.get("arrowhead") != "normal":
                errors.append(f"strict inclusion lacks its visual encoding: {ref}")
        elif relation["relation_type"] == "inclusion":
            if edge.get("color") != make_graph.INCLUSION_COLOR or edge.get("arrowhead") != "vee":
                errors.append(f"inclusion lacks its visual encoding: {ref}")

    nodes = {node["id"]: node for node in hierarchy["nodes"]}
    if set(nodes) != class_refs:
        errors.append("hierarchy nodes do not match the classes table")
    for entry in classes:
        node = nodes.get(f'#classes/{entry["short_name"]}', {})
        expected = make_graph._closure_style(entry)
        if any(node.get(key) != value for key, value in expected.items()):
            errors.append(f"closure profile is not encoded on node: {entry['short_name']}")

    # The minor atlas is the main Ample-class landscape.  Every atlas family
    # that is a subclass of Ample must have a real database node and an
    # established inclusion path to Ample.  The two excluded atlas rows are
    # explicitly outside that down-set.
    inventory = json.loads(args.minor_atlas_inventory.read_text(encoding="utf-8"))
    graph_map = json.loads(args.minor_atlas_graph_map.read_text(encoding="utf-8"))
    excluded_atlas_rows = {"minimal-nonample", "treelike"}
    expected_atlas_ids = {
        entry["id"] for entry in inventory["families"]
        if entry["id"] not in excluded_atlas_rows
    }
    mapped_atlas_ids = {entry["source_id"] for entry in graph_map["classes"]}
    if mapped_atlas_ids != expected_atlas_ids:
        missing = sorted(expected_atlas_ids - mapped_atlas_ids)
        extra = sorted(mapped_atlas_ids - expected_atlas_ids)
        errors.append(f"minor-atlas graph map mismatch; missing={missing}, extra={extra}")
    if len(mapped_atlas_ids) != len(graph_map["classes"]):
        errors.append("minor-atlas graph map contains duplicate source ids")
    ample_ref = "#classes/ample"
    for mapping in graph_map["classes"]:
        reference = aliases.get(mapping["class_id"])
        if not reference or reference not in nodes:
            errors.append(f"minor-atlas class is absent from hierarchy: {mapping['source_id']}")
        elif reference != ample_ref and not has_path(reference, ample_ref, adjacency):
            errors.append(
                f"minor-atlas Ample subclass lacks Hasse reachability to Ample: {mapping['source_id']}"
            )

    closure = make_graph.generate_minor_closure_map(cache)
    closure_refs = {edge.get("ref") for edge in closure["edges"]}
    for result in cache.get_table_entries("operation_results"):
        if result.get("status") == "established" and result.get("result_type") == "closure_equals":
            ref = f'#operation_results/{result["short_name"]}'
            if ref not in closure_refs:
                errors.append(f"exact closure arrow is absent: {ref}")

    for name, graph in (("class-hierarchy", hierarchy), ("minor-closure", closure)):
        path = site_dir / "graphs" / f"{name}.html"
        if not path.is_file():
            errors.append(f"missing generated graph: {path}")
            continue
        text = path.read_text(encoding="utf-8")
        if "strict digraph" in text:
            errors.append(f"{path}: strict DOT graph would collapse parallel facts")
        edge_count = text.count('id="graph-edge-')
        if edge_count != len(graph["edges"]):
            errors.append(
                f"{path}: generated {edge_count} DOT edges for {len(graph['edges'])} graph facts"
            )

    if errors:
        print("Graph check failed:", *errors, sep="\n- ", file=sys.stderr)
        return 1
    print(
        f"Graph check passed: {len(cover_pairs)} Hasse covers represent "
        f"{len(established)} stored inclusions; {len(mapped_atlas_ids)} minor-atlas Ample-family "
        f"nodes covered; {len(closure['edges'])} closure-map arrows rendered."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
