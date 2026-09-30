#!/usr/bin/env python3
"""Check that a complete, portable static site was generated."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


TABLES = {
    "classes",
    "characterizations",
    "characterization_equivalences",
    "relations",
    "operations",
    "operation_results",
    "invariants",
    "invariant_values",
    "obstruction_families",
    "obstruction_bases",
    "examples",
    "results",
}
GRAPHS = {"class-hierarchy", "minor-closure"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--site-dir", type=Path, default=Path("docs"))
    args = parser.parse_args()
    site_dir = args.site_dir.resolve()
    errors = []

    expected = {
        site_dir / "index.html",
        site_dir / "styles.css",
        site_dir / "bibliography.html",
        site_dir / ".nojekyll",
        *(site_dir / table / "index.html" for table in TABLES),
        *(site_dir / "graphs" / f"{graph}.html" for graph in GRAPHS),
    }
    for path in sorted(expected):
        if not path.is_file():
            errors.append(f"missing generated file: {path}")
        elif path.name != ".nojekyll" and path.stat().st_size == 0:
            errors.append(f"empty generated file: {path}")

    if (site_dir / "index.html").is_file():
        index = (site_dir / "index.html").read_text(encoding="utf-8")
        if "Partial Cubes" not in index:
            errors.append("main index does not contain the site title")
        for table in TABLES:
            if f"{table}/" not in index:
                errors.append(f"main index has no link to {table}")

    for path in site_dir.rglob("*.html") if site_dir.exists() else []:
        text = path.read_text(encoding="utf-8")
        for forbidden in ("file://", "/home/ec2-user/", "~/code/", "~/latex/"):
            if forbidden in text:
                errors.append(f"{path}: generated HTML exposes local path {forbidden!r}")

    if errors:
        print("Generated site check failed:", *errors, sep="\n- ", file=sys.stderr)
        return 1
    print(f"Generated site contains {len(TABLES)} tables and {len(GRAPHS)} graph views.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
