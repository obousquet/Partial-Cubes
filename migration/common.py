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
