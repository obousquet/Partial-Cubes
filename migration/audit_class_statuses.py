#!/usr/bin/env python3
"""Audit closure epistemic status and every class's position relative to Partial Cubes."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict, deque
from pathlib import Path


def load_records(data_dir: Path, table: str) -> list[dict]:
    return [
        json.loads(path.read_text(encoding="utf-8"))
        for path in sorted((data_dir / table).glob("[0-9]*.json"))
    ]


def emit(payload: dict, path: Path, check: bool) -> None:
    text = json.dumps(payload, indent=2, ensure_ascii=False) + "\n"
    if check:
        if not path.is_file() or path.read_text(encoding="utf-8") != text:
            raise SystemExit(f"Stale: {path}")
        print(f"Current: {path}")
    else:
        path.write_text(text, encoding="utf-8")
        print(f"Wrote: {path}")


def shortest_path(source: str, target: str, adjacency: dict[str, list[tuple[str, str]]]):
    pending = deque([source])
    previous: dict[str, tuple[str, str] | None] = {source: None}
    while pending:
        current = pending.popleft()
        if current == target:
            break
        for successor, relation_ref in adjacency.get(current, []):
            if successor not in previous:
                previous[successor] = (current, relation_ref)
                pending.append(successor)
    if target not in previous:
        return None
    nodes, relations = [target], []
    while nodes[-1] != source:
        prior, relation_ref = previous[nodes[-1]]
        nodes.append(prior)
        relations.append(relation_ref)
    return {"class_path": list(reversed(nodes)), "relation_path": list(reversed(relations))}


def audit(data_dir: Path) -> tuple[dict, dict]:
    classes = load_records(data_dir, "classes")
    characterizations = load_records(data_dir, "characterizations")
    relations = load_records(data_dir, "relations")
    operation_results = load_records(data_dir, "operation_results")
    operations = load_records(data_dir, "operations")

    class_aliases = {}
    for entry in classes:
        ref = f'#classes/{entry["short_name"]}'
        class_aliases[ref] = ref
        class_aliases[f'#classes/{entry["id"]}'] = ref
    operation_names = {f'#operations/{entry["short_name"]}': entry["short_name"] for entry in operations}
    operation_names.update({f'#operations/{entry["id"]}': entry["short_name"] for entry in operations})
    characterization_sources = {
        f'#characterizations/{entry["short_name"]}': entry.get("proof_source")
        for entry in characterizations
    }
    characterization_sources.update({
        f'#characterizations/{entry["id"]}': entry.get("proof_source")
        for entry in characterizations
    })

    result_refs: dict[tuple[str, str], list[str]] = defaultdict(list)
    for result in operation_results:
        class_ref = class_aliases.get(result.get("class_id"))
        operation_name = operation_names.get(result.get("operation_id"))
        if class_ref and operation_name in {"projection", "conditioning", "pc_minor"}:
            result_refs[(class_ref, operation_name)].append(
                f'#operation_results/{result["short_name"]}'
            )

    closure_rows = []
    prefixes = (("p", "projection"), ("c", "conditioning"), ("pc", "pc_minor"))
    for entry in sorted(classes, key=lambda row: row["id"]):
        class_ref = f'#classes/{entry["short_name"]}'
        row = {
            "class_id": class_ref,
            "primary_source": characterization_sources.get(entry["primary_characterization_id"]),
        }
        for prefix, operation_name in prefixes:
            row[f"{prefix}_closure_status"] = entry[f"{prefix}_closure_status"]
            refs = sorted(result_refs.get((class_ref, operation_name), []))
            if refs:
                row[f"{prefix}_operation_result_ids"] = refs
        closure_rows.append(row)
    closure = {
        "version": 1,
        "kind": "class_closure_status_audit",
        "generated_at": "2026-10-01",
        "semantics": {
            "closed": "Established closure theorem, backed by an operation result and cached Boolean true.",
            "not_closed": "Established counterexample, backed by an operation result and cached Boolean false.",
            "open": "An active source explicitly leaves the closure question open.",
            "not_assessed": "No source-backed closure determination has yet been promoted; this does not assert that the problem is mathematically open.",
            "not_applicable": "The operation has no canonical meaning at this class scope."
        },
        "classes": closure_rows,
    }

    adjacency: dict[str, list[tuple[str, str]]] = defaultdict(list)
    negative: dict[str, list[str]] = defaultdict(list)
    for relation in relations:
        if relation.get("status") != "established":
            continue
        subject = class_aliases.get(relation.get("subject_id"))
        object_ = class_aliases.get(relation.get("object_id"))
        if not subject or not object_:
            continue
        relation_ref = f'#relations/{relation["short_name"]}'
        if relation.get("relation_type") in {"inclusion", "strict_inclusion", "equality"}:
            adjacency[subject].append((object_, relation_ref))
            if relation.get("relation_type") == "equality":
                adjacency[object_].append((subject, relation_ref))
        if object_ == "#classes/partial_cubes" and relation.get("relation_type") in {"noncontainment", "incomparable"}:
            negative[subject].append(relation_ref)

    partial = "#classes/partial_cubes"
    containment_rows = []
    unresolved = []
    for entry in sorted(classes, key=lambda row: row["id"]):
        ref = f'#classes/{entry["short_name"]}'
        if ref == partial:
            row = {"class_id": ref, "disposition": "reference_class"}
        elif (path := shortest_path(ref, partial, adjacency)) is not None:
            row = {"class_id": ref, "disposition": "contained", **path}
        elif (path := shortest_path(partial, ref, adjacency)) is not None:
            row = {"class_id": ref, "disposition": "superclass", **path}
        elif negative.get(ref):
            row = {"class_id": ref, "disposition": "not_contained", "negative_relation_ids": sorted(negative[ref])}
        else:
            row = {"class_id": ref, "disposition": "unresolved"}
            unresolved.append(ref)
        containment_rows.append(row)
    containment = {
        "version": 1,
        "kind": "partial_cube_containment_audit",
        "generated_at": "2026-10-01",
        "counts": {
            key: sum(row["disposition"] == key for row in containment_rows)
            for key in ("reference_class", "contained", "superclass", "not_contained", "unresolved")
        },
        "unresolved_class_ids": unresolved,
        "classes": containment_rows,
    }
    return closure, containment


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--migration-dir", type=Path, default=Path("migration"))
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    closure, containment = audit(args.data_dir.resolve())
    migration_dir = args.migration_dir.resolve()
    emit(closure, migration_dir / "closure_status_audit.json", args.check)
    emit(containment, migration_dir / "partial_cube_containment_audit.json", args.check)
    if containment["unresolved_class_ids"]:
        raise SystemExit("Unresolved containment positions: " + ", ".join(containment["unresolved_class_ids"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
