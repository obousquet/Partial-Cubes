#!/usr/bin/env python3
"""Check that pinned migration inputs still match their recorded snapshot."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

from common import DEFAULT_CONFIG, MIGRATION_DIR, git_head, load_config, load_json, sha256_bytes, sha256_file


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--manifest", type=Path, default=Path(__file__).with_name("source_manifest.json"))
    parser.add_argument("--latex-root", type=Path)
    parser.add_argument("--archive-index", type=Path, default=MIGRATION_DIR / "source_archive.json")
    args = parser.parse_args()

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

    archive = load_json(args.archive_index) if args.archive_index.is_file() else {"files": []}
    archived = {
        (entry["repository"], entry["path"]): MIGRATION_DIR / entry["archive_path"]
        for entry in archive["files"]
    }
    heads = {entry["name"]: entry["head"] for entry in manifest["repositories"]}
    head_drift = []
    for repository in manifest["repositories"]:
        name = repository["name"]
        actual = git_head(roots[name])
        if actual != repository["head"]:
            head_drift.append(f"{name}: HEAD {repository['head']} -> {actual}")
    live_drift = []
    unrecoverable = []
    for entry in manifest["files"]:
        key = (entry["repository"], entry["path"])
        path = roots[entry["repository"]] / entry["path"]
        if path.is_file() and sha256_file(path) == entry["sha256"]:
            continue
        live_drift.append(f"{entry['repository']}:{entry['path']}")
        archive_path = archived.get(key)
        if archive_path and archive_path.is_file() and sha256_file(archive_path) == entry["sha256"]:
            continue
        if entry["git_state"] == "tracked":
            result = subprocess.run(
                ["git", "-C", str(roots[entry["repository"]]), "show", f"{heads[entry['repository']]}:{entry['path']}"],
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            if result.returncode == 0 and sha256_bytes(result.stdout) == entry["sha256"]:
                continue
        unrecoverable.append(f"{entry['repository']}:{entry['path']}")

    if unrecoverable:
        print("Frozen source bytes are not recoverable:")
        for item in unrecoverable:
            print(f"- {item}")
        return 1
    print(
        f"Source snapshot is recoverable: {len(manifest['files'])} files in "
        f"{len(manifest['repositories'])} repositories; {len(archive['files'])} dirty files byte-archived."
    )
    if head_drift:
        print(f"Repository HEAD drift recorded: {len(head_drift)}")
    if live_drift:
        print(f"Live-file drift safely isolated from snapshot: {len(live_drift)}")
    if manifest["missing_references"]:
        print(f"Recorded missing source references: {len(manifest['missing_references'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
