#!/usr/bin/env python3
"""Inventory BibTeX candidates and report key or work-identity collisions."""

from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

from common import DEFAULT_CONFIG, SourceSnapshot, emit_json, load_config, sha256_text


ENTRY_START_RE = re.compile(r"@(?P<type>[A-Za-z]+)\s*(?P<open>[{(])", re.MULTILINE)


def bib_entries(text: str) -> list[dict]:
    entries = []
    position = 0
    while True:
        match = ENTRY_START_RE.search(text, position)
        if not match:
            break
        opening = match.group("open")
        closing = "}" if opening == "{" else ")"
        depth = 1
        cursor = match.end()
        escaped = False
        while cursor < len(text) and depth:
            character = text[cursor]
            if escaped:
                escaped = False
            elif character == "\\":
                escaped = True
            elif character == opening:
                depth += 1
            elif character == closing:
                depth -= 1
            cursor += 1
        if depth:
            raise ValueError(f"Unterminated BibTeX entry at offset {match.start()}")
        raw = text[match.start() : cursor]
        body = text[match.end() : cursor - 1]
        entry_type = match.group("type").lower()
        comma = body.find(",")
        key = body[:comma].strip() if comma >= 0 else ""
        if entry_type not in {"comment", "preamble", "string"} and key:
            entries.append(
                {
                    "entry_type": entry_type,
                    "key": key,
                    "body": body[comma + 1 :],
                    "raw_sha256": sha256_text(raw.strip()),
                }
            )
        position = cursor
    return entries


def field(body: str, name: str) -> str | None:
    match = re.search(
        rf"(?is)(?:^|,)\s*{re.escape(name)}\s*=\s*(?:\{{((?:[^{{}}]|\{{[^{{}}]*\}})*)\}}|\"([^\"]*)\")",
        body,
    )
    if not match:
        return None
    return re.sub(r"\s+", " ", next(value for value in match.groups() if value is not None).strip())


def normalize_doi(value: str | None) -> str | None:
    if not value:
        return None
    normalized = value.lower().strip()
    normalized = re.sub(r"^https?://(?:dx\.)?doi\.org/", "", normalized)
    return normalized


def normalize_title(value: str | None) -> str | None:
    if not value:
        return None
    value = re.sub(r"\\[A-Za-z]+\*?(?:\[[^\]]*\])?", "", value)
    value = value.replace("{", "").replace("}", "")
    normalized = re.sub(r"[^a-z0-9]+", " ", value.lower()).strip()
    return normalized or None


def normalize_eprint(value: str | None) -> str | None:
    if not value:
        return None
    return re.sub(r"v\d+$", "", value.lower().strip())


def audit(latex_root: Path, config: dict, snapshot: SourceSnapshot) -> dict:
    occurrences = []
    source_counts = {}
    for relative in config["bibliographies"]:
        path = latex_root / relative
        entries = bib_entries(snapshot.read_text(path))
        source_counts[relative] = len(entries)
        for entry in entries:
            doi = normalize_doi(field(entry["body"], "doi"))
            eprint = normalize_eprint(field(entry["body"], "eprint"))
            title = normalize_title(field(entry["body"], "title"))
            identity = f"doi:{doi}" if doi else f"arxiv:{eprint.lower()}" if eprint else f"title:{title}" if title else None
            identities = [
                token
                for token in (
                    f"doi:{doi}" if doi else None,
                    f"arxiv:{eprint}" if eprint else None,
                    f"title:{title}" if title else None,
                )
                if token
            ]
            occurrences.append(
                {
                    "source": relative,
                    "key": entry["key"],
                    "entry_type": entry["entry_type"],
                    "identity": identity,
                    "identities": identities,
                    "doi": doi,
                    "eprint": eprint,
                    "normalized_title": title,
                    "raw_sha256": entry["raw_sha256"],
                    "migration_state": "candidate",
                }
            )

    by_key: dict[str, list[dict]] = defaultdict(list)
    by_identity: dict[str, list[dict]] = defaultdict(list)
    for occurrence in occurrences:
        by_key[occurrence["key"]].append(occurrence)
        if occurrence["identity"]:
            by_identity[occurrence["identity"]].append(occurrence)

    repeated_keys = {}
    conflicting_keys = {}
    for key, items in sorted(by_key.items()):
        if len(items) < 2:
            continue
        summary = [
            {"source": item["source"], "identity": item["identity"], "raw_sha256": item["raw_sha256"]}
            for item in items
        ]
        repeated_keys[key] = summary
        common_identities = set(items[0]["identities"])
        for item in items[1:]:
            common_identities.intersection_update(item["identities"])
        if not common_identities and len({item["raw_sha256"] for item in items}) > 1:
            conflicting_keys[key] = summary

    duplicate_identities = {
        identity: [{"source": item["source"], "key": item["key"]} for item in items]
        for identity, items in sorted(by_identity.items())
        if len({(item["source"], item["key"]) for item in items}) > 1
    }
    return {
        "version": 1,
        "kind": "bibliography_migration_queue",
        "counts": {
            "source_files": len(source_counts),
            "entry_occurrences": len(occurrences),
            "unique_keys": len(by_key),
            "repeated_keys": len(repeated_keys),
            "conflicting_keys": len(conflicting_keys),
            "duplicate_identities": len(duplicate_identities),
        },
        "entries_by_source": source_counts,
        "conflicting_keys": conflicting_keys,
        "repeated_keys": repeated_keys,
        "duplicate_identities": duplicate_identities,
        "occurrences": occurrences,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--latex-root", type=Path)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name("bibliography_queue.json"),
    )
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    config = load_config(args.config)
    latex_root = (args.latex_root or Path(config["latex_repository"]).expanduser()).resolve()
    snapshot = SourceSnapshot(latex_root, config)
    emit_json(audit(latex_root, config, snapshot), args.output.resolve(), args.check)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
