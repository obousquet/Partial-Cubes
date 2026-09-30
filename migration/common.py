#!/usr/bin/env python3
"""Shared, standard-library helpers for migration tools."""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any


MIGRATION_DIR = Path(__file__).resolve().parent
DEFAULT_CONFIG = MIGRATION_DIR / "config.json"


def load_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def load_config(path: Path = DEFAULT_CONFIG) -> dict[str, Any]:
    return load_json(path)


def canonical_json(data: Any) -> str:
    return json.dumps(data, indent=2, ensure_ascii=False) + "\n"


def emit_json(data: Any, output: Path, check: bool = False) -> None:
    content = canonical_json(data)
    if check:
        if not output.is_file() or output.read_text(encoding="utf-8") != content:
            raise SystemExit(f"Generated migration artifact is stale: {output}")
        print(f"Current: {output}")
        return
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(content, encoding="utf-8")
    print(f"Wrote: {output}")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def run_git(root: Path, *args: str, check: bool = True) -> str:
    result = subprocess.run(
        ["git", "-C", str(root), *args],
        check=check,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    return result.stdout.strip()


def git_head(root: Path) -> str:
    return run_git(root, "rev-parse", "HEAD")


def git_status_lines(root: Path) -> list[str]:
    output = run_git(root, "status", "--porcelain=v1", "--untracked-files=all")
    return output.splitlines() if output else []


def path_git_state(root: Path, relative_path: str) -> str:
    tracked = subprocess.run(
        ["git", "-C", str(root), "ls-files", "--error-unmatch", "--", relative_path],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    ).returncode == 0
    if not tracked:
        return "untracked"
    unstaged = subprocess.run(
        ["git", "-C", str(root), "diff", "--quiet", "--", relative_path]
    ).returncode != 0
    staged = subprocess.run(
        ["git", "-C", str(root), "diff", "--cached", "--quiet", "--", relative_path]
    ).returncode != 0
    if staged and unstaged:
        return "staged_and_modified"
    if staged:
        return "staged"
    if unstaged:
        return "modified"
    return "tracked"


def latex_labels(path: Path) -> list[str]:
    """Return labels in source order without trying to interpret TeX."""
    import re

    text = path.read_text(encoding="utf-8", errors="replace")
    return re.findall(r"\\label\s*\{([^{}]+)\}", text)


def line_number(text: str, offset: int) -> int:
    return text.count("\n", 0, offset) + 1


def logical_repo_path(path: Path, repositories: dict[str, Path]) -> tuple[str, str]:
    resolved = path.resolve()
    matches: list[tuple[int, str, Path]] = []
    for name, root in repositories.items():
        try:
            resolved.relative_to(root.resolve())
        except ValueError:
            continue
        matches.append((len(root.resolve().parts), name, root.resolve()))
    if not matches:
        raise ValueError(f"Path is outside configured repositories: {path}")
    _, name, root = max(matches)
    return name, resolved.relative_to(root).as_posix()


def resolve_source_path(latex_root: Path, source_path: str) -> Path:
    return (latex_root / source_path).resolve()


class SourceSnapshot:
    """Read a pinned file from the live tree, byte archive, or recorded Git tree."""

    def __init__(
        self,
        latex_root: Path,
        config: dict[str, Any],
        manifest_path: Path | None = None,
        archive_index_path: Path | None = None,
    ) -> None:
        self.latex_root = latex_root.resolve()
        self.config = config
        self.repositories = {"PartialCubes": self.latex_root}
        self.repositories.update(
            {
                name: (self.latex_root / relative).resolve()
                for name, relative in config["external_repositories"].items()
            }
        )
        self.manifest_path = manifest_path or MIGRATION_DIR / "source_manifest.json"
        self.manifest = load_json(self.manifest_path) if self.manifest_path.is_file() else None
        self.entries = {}
        self.heads = {}
        if self.manifest:
            self.entries = {
                (entry["repository"], entry["path"]): entry
                for entry in self.manifest["files"]
            }
            self.heads = {
                entry["name"]: entry["head"]
                for entry in self.manifest["repositories"]
            }
        archive_path = archive_index_path or MIGRATION_DIR / "source_archive.json"
        self.archive_index = load_json(archive_path) if archive_path.is_file() else None
        self.archives = {}
        if self.archive_index:
            self.archives = {
                (entry["repository"], entry["path"]): MIGRATION_DIR / entry["archive_path"]
                for entry in self.archive_index["files"]
            }

    def key(self, path: Path) -> tuple[str, str]:
        return logical_repo_path(path, self.repositories)

    def read_bytes(self, path: Path) -> bytes:
        path = path.resolve()
        if not self.manifest:
            return path.read_bytes()
        key = self.key(path)
        entry = self.entries.get(key)
        if not entry:
            raise FileNotFoundError(f"Path is not in the frozen source manifest: {key[0]}:{key[1]}")
        if path.is_file():
            live = path.read_bytes()
            if sha256_bytes(live) == entry["sha256"]:
                return live
        archive = self.archives.get(key)
        if archive and archive.is_file():
            content = archive.read_bytes()
            if sha256_bytes(content) == entry["sha256"]:
                return content
        if entry["git_state"] == "tracked":
            result = subprocess.run(
                ["git", "-C", str(self.repositories[key[0]]), "show", f"{self.heads[key[0]]}:{key[1]}"],
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            if result.returncode == 0 and sha256_bytes(result.stdout) == entry["sha256"]:
                return result.stdout
        raise RuntimeError(f"Frozen source is not recoverable: {key[0]}:{key[1]}")

    def read_text(self, path: Path) -> str:
        return self.read_bytes(path).decode("utf-8", errors="replace")

    def recorded_sha256(self, path: Path) -> str:
        if not self.manifest:
            return sha256_file(path)
        key = self.key(path.resolve())
        if key not in self.entries:
            raise FileNotFoundError(f"Path is not in the frozen source manifest: {key[0]}:{key[1]}")
        return self.entries[key]["sha256"]
