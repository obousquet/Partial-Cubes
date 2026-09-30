#!/usr/bin/env python3
"""Byte-archive manifest files that cannot be recovered from a clean Git tree."""

from __future__ import annotations

import argparse
from pathlib import Path

from common import DEFAULT_CONFIG, MIGRATION_DIR, canonical_json, load_config, load_json, sha256_bytes


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--manifest", type=Path, default=MIGRATION_DIR / "source_manifest.json")
    parser.add_argument("--latex-root", type=Path)
    parser.add_argument("--output", type=Path, default=MIGRATION_DIR / "source_archive.json")
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit(f"Refusing to overwrite source archive index: {args.output}")

    config = load_config(args.config)
    manifest = load_json(args.manifest)
    latex_root = (args.latex_root or Path(config["latex_repository"]).expanduser()).resolve()
    roots = {"PartialCubes": latex_root}
    roots.update(
        {
            name: (latex_root / relative).resolve()
            for name, relative in config["external_repositories"].items()
        }
    )
    archive_root = MIGRATION_DIR / "snapshots" / manifest["snapshot_id"]
    archived = []
    for entry in manifest["files"]:
        if entry["git_state"] == "tracked":
            continue
        source = roots[entry["repository"]] / entry["path"]
        content = source.read_bytes()
        actual = sha256_bytes(content)
        if actual != entry["sha256"]:
            raise SystemExit(
                f"Source changed after manifest creation: {entry['repository']}:{entry['path']}"
            )
        destination = archive_root / entry["repository"] / entry["path"]
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(content)
        archived.append(
            {
                "repository": entry["repository"],
                "path": entry["path"],
                "sha256": entry["sha256"],
                "size_bytes": entry["size_bytes"],
                "git_state": entry["git_state"],
                "archive_path": destination.relative_to(MIGRATION_DIR).as_posix(),
            }
        )
    result = {
        "version": 1,
        "kind": "dirty_source_byte_archive",
        "snapshot_id": manifest["snapshot_id"],
        "source_manifest": args.manifest.relative_to(MIGRATION_DIR).as_posix(),
        "count": len(archived),
        "files": archived,
    }
    args.output.write_text(canonical_json(result), encoding="utf-8")
    print(f"Archived {len(archived)} dirty source files under {archive_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
