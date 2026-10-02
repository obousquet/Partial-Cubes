#!/usr/bin/env python3
"""Exhaustive small-cube regression for the signed-facet shell theorem.

This verifier is deliberately independent of the campaign search code.  It
checks every proper row in Q_d, d <= max_dimension, against the structural
claims used in Lemma ``damp-row-facet-shell-core``:

* the canonical signed-facet shell decomposition;
* equivalence of ampleness and corner peelability with the residual core;
* the corner formula and safe deletion formula when the residual is empty.

Run nontrivial instances on the designated worker.
"""

from __future__ import annotations

import argparse
import json
from functools import lru_cache


def projection_count(row: frozenset[int], dimension: int, mask: int) -> int:
    coordinates = [e for e in range(dimension) if mask & (1 << e)]
    traces = {
        sum(((word >> e) & 1) << j for j, e in enumerate(coordinates))
        for word in row
    }
    return len(traces)


@lru_cache(maxsize=None)
def is_ample(dimension: int, row_tuple: tuple[int, ...]) -> bool:
    row = frozenset(row_tuple)
    shattered = sum(
        projection_count(row, dimension, mask) == 1 << bin(mask).count("1")
        for mask in range(1 << dimension)
    )
    return shattered == len(row)


@lru_cache(maxsize=None)
def is_corner_peelable(dimension: int, row_tuple: tuple[int, ...]) -> bool:
    if len(row_tuple) <= 1:
        return True
    for word in row_tuple:
        remainder = tuple(x for x in row_tuple if x != word)
        if is_ample(dimension, remainder) and is_corner_peelable(
            dimension, remainder
        ):
            return True
    return False


def canonical_core(
    row: frozenset[int], dimension: int
) -> tuple[list[tuple[int, int]], frozenset[int], frozenset[int]]:
    full_facets: list[tuple[int, int]] = []
    for e in range(dimension):
        for bit in (0, 1):
            facet = {
                word
                for word in range(1 << dimension)
                if ((word >> e) & 1) == bit
            }
            if facet <= row:
                full_facets.append((e, bit))

    opposite_face = frozenset(
        word
        for word in range(1 << dimension)
        if all(((word >> e) & 1) == 1 - bit for e, bit in full_facets)
    )
    return full_facets, opposite_face, row & opposite_face


def suppress_fixed_coordinates(
    row: frozenset[int], dimension: int, full_facets: list[tuple[int, int]]
) -> tuple[int, tuple[int, ...]]:
    fixed = {e for e, _ in full_facets}
    free = [e for e in range(dimension) if e not in fixed]
    projected = tuple(
        sorted(
            sum(((word >> e) & 1) << j for j, e in enumerate(free))
            for word in row
        )
    )
    return len(free), projected


def verify(max_dimension: int) -> dict[str, object]:
    checked_by_dimension: dict[str, int] = {}
    for dimension in range(1, max_dimension + 1):
        checked = 0
        universe = frozenset(range(1 << dimension))
        for indicator in range(1 << (1 << dimension)):
            row = frozenset(
                word
                for word in universe
                if indicator & (1 << word)
            )
            if row == universe:
                continue

            full_facets, opposite_face, core = canonical_core(row, dimension)
            shell = universe - opposite_face
            assert row == shell | core
            assert not (shell & core)

            residual_dimension, residual_tuple = suppress_fixed_coordinates(
                core, dimension, full_facets
            )
            row_tuple = tuple(sorted(row))
            row_ample = is_ample(dimension, row_tuple)
            core_ample = is_ample(residual_dimension, residual_tuple)
            assert row_ample == core_ample

            row_peelable = (
                is_corner_peelable(dimension, row_tuple) if row_ample else False
            )
            core_peelable = (
                is_corner_peelable(residual_dimension, residual_tuple)
                if core_ample
                else False
            )
            assert row_peelable == core_peelable

            if not core and row:
                for word in row:
                    incident_full_facets = sum(
                        ((word >> e) & 1) == bit for e, bit in full_facets
                    )
                    remainder = tuple(sorted(row - {word}))
                    is_corner = is_ample(dimension, remainder)
                    assert is_corner == (incident_full_facets == 1)
                    if is_corner:
                        assert is_corner_peelable(dimension, remainder)

            checked += 1
        checked_by_dimension[str(dimension)] = checked

    return {
        "schema": "damp-terminal-facet-partition-regression-v1",
        "status": "ok",
        "max_dimension": max_dimension,
        "proper_rows_by_dimension": checked_by_dimension,
        "proper_rows_checked": sum(checked_by_dimension.values()),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-dimension", type=int, default=4)
    args = parser.parse_args()
    print(json.dumps(verify(args.max_dimension), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
