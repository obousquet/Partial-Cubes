#!/usr/bin/env python3
"""Build the exact, repository-aware source snapshot for a migration batch."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

from common import (
    DEFAULT_CONFIG,
    emit_json,
    git_head,
    git_status_lines,
    latex_labels,
    load_config,
    load_json,
    logical_repo_path,
    path_git_state,
    resolve_source_path,
    sha256_file,
)


def repositories_for(latex_root: Path, config: dict) -> dict[str, Path]:
    repositories = {"PartialCubes": latex_root.resolve()}
    for name, relative in config["external_repositories"].items():
        repositories[name] = (latex_root / relative).resolve()
    return repositories


def collect_paths(latex_root: Path, config: dict) -> tuple[dict[Path, set[str]], list[dict]]:
    paths: dict[Path, set[str]] = defaultdict(set)
    configured_papers = set(config["paper_workspaces"])
    for path in sorted((latex_root / "papers").glob("**/*.tex")):
        if path.relative_to(latex_root / "papers").parts[0] in configured_papers:
            paths[path.resolve()].add("paper_tex")
    for relative in config["bibliographies"]:
        paths[(latex_root / relative).resolve()].add("bibliography")
    for role, relative in config["structured_sources"].items():
        paths[(latex_root / relative).resolve()].add(role)

    inventory_path = latex_root / config["structured_sources"]["minor_inventory"]
    inventory = load_json(inventory_path)
    missing = []
    for packet in inventory["sources"]:
        for field in ("path", "pdf"):
            value = packet.get(field)
            if not value:
                continue
            path = resolve_source_path(latex_root, value)
            if path.is_file():
                paths[path].add(f"minor_source_packet:{packet['id']}:{field}")
            else:
                missing.append(
                    {
                        "source_packet": packet["id"],
                        "field": field,
                        "configured_path": value,
                        "resolved_path": str(path),
                    }
                )
        for label, source in packet.get("label_sources", {}).items():
            path = resolve_source_path(latex_root, source["path"])
            if path.is_file():
                paths[path].add(f"minor_label_source:{packet['id']}:{label}")
            else:
                missing.append(
                    {
                        "source_packet": packet["id"],
                        "field": f"label_source:{label}",
                        "configured_path": source["path"],
                        "resolved_path": str(path),
                    }
                )
    return paths, missing


def build(latex_root: Path, config: dict) -> dict:
    repositories = repositories_for(latex_root, config)
    for name, root in repositories.items():
        if not (root / ".git").exists():
            raise SystemExit(f"Configured repository does not exist: {name} at {root}")

    paths, missing = collect_paths(latex_root, config)
    files = []
    repo_file_counts: dict[str, int] = defaultdict(int)
    for path in sorted(paths, key=str):
        repository, relative = logical_repo_path(path, repositories)
        root = repositories[repository]
        entry = {
            "repository": repository,
            "path": relative,
            "roles": sorted(paths[path]),
            "sha256": sha256_file(path),
            "size_bytes": path.stat().st_size,
            "git_state": path_git_state(root, relative),
        }
        if path.suffix == ".tex":
            entry["labels"] = latex_labels(path)
        files.append(entry)
        repo_file_counts[repository] += 1

    repository_entries = []
    for name, root in repositories.items():
        status = git_status_lines(root)
        repository_entries.append(
            {
                "name": name,
                "relative_to_latex_root": "." if name == "PartialCubes" else config["external_repositories"][name],
                "head": git_head(root),
                "dirty": bool(status),
                "dirty_path_count": len(status),
                "dirty_paths": status,
                "manifest_file_count": repo_file_counts.get(name, 0),
            }
        )

    return {
        "version": 1,
        "kind": "migration_source_manifest",
        "snapshot_date": config["snapshot_date"],
        "latex_root_logical": "~/latex/PartialCubes",
        "selection_rule": "Configured paper TeX and bibliographies, structured seed files, and every file or PDF named by a minor-atlas source packet.",
        "repositories": repository_entries,
        "counts": {
            "repositories": len(repository_entries),
            "files": len(files),
            "missing_references": len(missing),
        },
        "missing_references": missing,
        "files": files,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--latex-root", type=Path)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name("source_manifest.json"),
    )
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    config = load_config(args.config)
    latex_root = (args.latex_root or Path(config["latex_repository"]).expanduser()).resolve()
    emit_json(build(latex_root, config), args.output.resolve(), args.check)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
