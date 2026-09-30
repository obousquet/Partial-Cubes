#!/usr/bin/env python3
"""Normalize generated HTML so committed Pages output has stable whitespace."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--site-dir", type=Path, default=Path("docs"))
    args = parser.parse_args()
    changed = 0
    for path in sorted(args.site_dir.rglob("*.html")):
        content = path.read_text(encoding="utf-8")
        lines = [line.rstrip() for line in content.splitlines()]
        while lines and not lines[-1]:
            lines.pop()
        normalized = "\n".join(lines) + "\n"
        if normalized != content:
            path.write_text(normalized, encoding="utf-8")
            changed += 1
    print(f"Normalized {changed} generated HTML files.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
