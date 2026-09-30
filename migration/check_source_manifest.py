#!/usr/bin/env python3
"""Check that pinned migration inputs still match their recorded snapshot."""

from __future__ import annotations

import argparse
from pathlib import Path

from common import DEFAULT_CONFIG, git_head, load_config, load_json, sha256_file


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--manifest", type=Path, default=Path(__file__).with_name("source_manifest.json"))
    parser.add_argument("--latex-root", type=Path)
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

    drift = []
    for repository in manifest["repositories"]:
        name = repository["name"]
        actual = git_head(roots[name])
        if actual != repository["head"]:
            drift.append(f"{name}: HEAD {repository['head']} -> {actual}")
    for entry in manifest["files"]:
        path = roots[entry["repository"]] / entry["path"]
        if not path.is_file():
            drift.append(f"{entry['repository']}:{entry['path']}: missing")
            continue
        actual = sha256_file(path)
        if actual != entry["sha256"]:
            drift.append(
                f"{entry['repository']}:{entry['path']}: sha256 {entry['sha256']} -> {actual}"
            )

    if drift:
        print("Source snapshot drift detected:")
        for item in drift:
            print(f"- {item}")
        return 1
    print(
        f"Source snapshot is current: {len(manifest['files'])} files in "
        f"{len(manifest['repositories'])} repositories."
    )
    if manifest["missing_references"]:
        print(f"Recorded missing source references: {len(manifest['missing_references'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
