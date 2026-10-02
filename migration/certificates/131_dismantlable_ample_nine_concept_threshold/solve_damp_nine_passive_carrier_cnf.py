#!/usr/bin/env python3
"""Certify the passive-coordinate-safe nine-facet carrier reduction.

The formula does not enumerate lower or upper cube families.  It encodes
only consequences of the two-shore trace-defect theorem:

* every doubly blocked residual vertex of degree two marks an isolated
  relative support;
* a doubly blocked degree-three vertex marks either an isolated support or
  a minimal--maximal support cover;
* every compatible lower--upper carrier intersection is a full cube in the
  residual graph; and
* double blocking means that every residual vertex has lower and upper
  carrier depth at least two.

The active graph classification first reduces thirteen coordinate-coloured
nine-vertex graphs to the three profiles below.  This script proves that
even those profiles have no abstract carrier system.  Passive parts of the
nine supports are represented on nine formal coordinates.  This loses no
generality: the inclusion preorder of any nine sets has a faithful
down-set representation on at most nine coordinates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

from z3 import (
    And,
    BitVec,
    BitVecVal,
    Bool,
    Extract,
    Goal,
    If,
    Implies,
    Int,
    Not,
    Or,
    Solver,
    Sum,
    Then,
    ULE,
    unknown,
)


CAMPAIGN = Path(__file__).resolve().parent
PROFILES = {
    "5.1": (0, 2, 4, 6, 7, 8, 10, 11, 15),
    "8.1": (0, 2, 4, 6, 8, 10, 11, 14, 15),
    "9.1": (0, 2, 3, 4, 6, 7, 8, 10, 11),
}
ACTIVE_DIMENSION = 4
SUPPORT_COUNT = 9


def residual_degrees(vertices: tuple[int, ...]) -> dict[int, tuple[int, ...]]:
    vertex_set = set(vertices)
    return {
        vertex: tuple(
            coordinate
            for coordinate in range(ACTIVE_DIMENSION)
            if vertex ^ (1 << coordinate) in vertex_set
        )
        for vertex in vertices
    }


def contained_active_cubes(
    vertices: tuple[int, ...],
) -> tuple[tuple[int, int], ...]:
    vertex_set = set(vertices)
    cubes = []
    for free in range(1, 1 << ACTIVE_DIMENSION):
        for base in range(1 << ACTIVE_DIMENSION):
            if base & free:
                continue
            submask = free
            while (base ^ submask) in vertex_set:
                if submask == 0:
                    cubes.append((free, base))
                    break
                submask = (submask - 1) & free
    return tuple(cubes)


def build_formula(profile: str) -> tuple[Solver, dict[str, object]]:
    vertices = PROFILES[profile]
    vertex_index = {vertex: index for index, vertex in enumerate(vertices)}
    directions = residual_degrees(vertices)
    degree_two = tuple(v for v in vertices if len(directions[v]) == 2)
    degree_three = tuple(v for v in vertices if len(directions[v]) == 3)
    cubes = contained_active_cubes(vertices)
    width = ACTIVE_DIMENSION
    count = SUPPORT_COUNT

    solver = Solver()
    active = [BitVec(f"active_{index}", width) for index in range(count)]
    passive = [
        [Bool(f"passive_{index}_{coordinate}") for coordinate in range(count)]
        for index in range(count)
    ]
    lower_anchor = [Int(f"lower_anchor_{index}") for index in range(count)]
    upper_anchor = [Int(f"upper_anchor_{index}") for index in range(count)]
    lower_word = [BitVec(f"lower_word_{index}", width) for index in range(count)]
    upper_word = [BitVec(f"upper_word_{index}", width) for index in range(count)]

    for anchor in lower_anchor + upper_anchor:
        solver.add(anchor >= 0, anchor < count)

    def active_subset(first: int, second: int):
        return (active[first] & ~active[second]) == 0

    def passive_subset(first: int, second: int):
        return And(
            *(
                Implies(passive[first][coordinate], passive[second][coordinate])
                for coordinate in range(count)
            )
        )

    def support_subset(first: int, second: int):
        return And(active_subset(first, second), passive_subset(first, second))

    def passive_equal(first: int, second: int):
        return And(
            *(
                passive[first][coordinate] == passive[second][coordinate]
                for coordinate in range(count)
            )
        )

    for first in range(count):
        for second in range(first + 1, count):
            solver.add(
                Or(
                    active[first] != active[second],
                    *(
                        passive[first][coordinate]
                        != passive[second][coordinate]
                        for coordinate in range(count)
                    ),
                )
            )

    minimal = [
        And(
            *(Not(support_subset(other, index)) for other in range(count) if other != index)
        )
        for index in range(count)
    ]
    maximal = [
        And(
            *(Not(support_subset(index, other)) for other in range(count) if other != index)
        )
        for index in range(count)
    ]

    for index in range(count):
        solver.add(
            Or(
                *(
                    And(lower_anchor[index] == anchor, lower_word[index] == vertex)
                    for anchor, vertex in enumerate(vertices)
                )
            )
        )
        solver.add(
            Or(
                *(
                    And(upper_anchor[index] == anchor, upper_word[index] == vertex)
                    for anchor, vertex in enumerate(vertices)
                )
            )
        )

    def lower_carrier(index: int, vertex: int):
        return ((BitVecVal(vertex, width) ^ lower_word[index]) & active[index]) == 0

    def upper_carrier(index: int, vertex: int):
        return ((BitVecVal(vertex, width) ^ upper_word[index]) & ~active[index]) == 0

    # The degree-two marker theorem supplies distinct isolated supports.
    for index, vertex in enumerate(degree_two):
        solver.add(minimal[index], maximal[index])
        solver.add(
            lower_anchor[index] == vertex_index[vertex],
            upper_anchor[index] == vertex_index[vertex],
        )
        solver.add(
            Sum(
                *(
                    If(Extract(coordinate, coordinate, active[index]) == 1, 1, 0)
                    for coordinate in directions[vertex]
                )
            )
            == 1
        )

    # A blocked singleton move has carrier depth at least two on each shore.
    for vertex in vertices:
        solver.add(
            Sum(
                *(
                    If(And(minimal[index], lower_carrier(index, vertex)), 1, 0)
                    for index in range(count)
                )
            )
            >= 2
        )
        solver.add(
            Sum(
                *(
                    If(And(maximal[index], upper_carrier(index, vertex)), 1, 0)
                    for index in range(count)
                )
            )
            >= 2
        )

    def active_cube(free, base):
        return Or(
            *(And(free == cube_free, (base & ~free) == cube_base) for cube_free, cube_base in cubes)
        )

    # Compatible minimal--maximal traces produce a full residual cube.  Since
    # residual concepts agree on passive coordinates, the maximal support's
    # passive part must be contained in the minimal support's passive part.
    for first in range(count):
        for second in range(count):
            if first == second:
                continue
            compatible = (
                (lower_word[first] ^ upper_word[second])
                & (active[first] & ~active[second])
            ) == 0
            free = active[second] & ~active[first]
            base = (
                (lower_word[first] & active[first])
                | (upper_word[second] & ~active[second])
            )
            solver.add(
                Implies(
                    And(minimal[first], maximal[second], compatible),
                    And(active_cube(free, base), passive_subset(second, first)),
                )
            )

    # The degree-three local dichotomy: isolated marker, or a support cover
    # with one active coordinate and no passive-coordinate difference.
    for vertex in degree_three:
        alternatives = []
        for first in range(count):
            alternatives.append(
                And(
                    minimal[first],
                    maximal[first],
                    lower_anchor[first] == vertex_index[vertex],
                    upper_anchor[first] == vertex_index[vertex],
                )
            )
            for second in range(count):
                if first == second:
                    continue
                difference = active[second] & ~active[first]
                one_coordinate = And(
                    difference != 0,
                    (difference & (difference - 1)) == 0,
                )
                alternatives.append(
                    And(
                        minimal[first],
                        maximal[second],
                        support_subset(first, second),
                        passive_equal(first, second),
                        one_coordinate,
                        lower_carrier(first, vertex),
                        upper_carrier(second, vertex),
                    )
                )
        solver.add(Or(*alternatives))

    # The support indices not fixed by degree-two markers are interchangeable.
    for index in range(len(degree_two), count - 1):
        solver.add(ULE(active[index], active[index + 1]))

    metadata = {
        "profile": profile,
        "active_dimension": ACTIVE_DIMENSION,
        "vertices": list(vertices),
        "degree_two_vertices": list(degree_two),
        "degree_three_vertices": list(degree_three),
        "contained_nontrivial_active_cubes": [list(cube) for cube in cubes],
        "support_count": count,
        "passive_preorder_coordinate_count": count,
        "assertion_count": len(solver.assertions()),
    }
    return solver, metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--profile", choices=sorted(PROFILES), required=True)
    parser.add_argument("--timeout-ms", type=int, default=300_000)
    parser.add_argument("--dimacs", type=Path)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()

    solver, metadata = build_formula(arguments.profile)
    solver.set(timeout=arguments.timeout_ms)
    smt2 = solver.to_smt2()
    metadata["smt2_sha256"] = hashlib.sha256(smt2.encode()).hexdigest()

    if arguments.dimacs is not None:
        goal = Goal()
        goal.add(*solver.assertions())
        compiled = Then("simplify", "card2bv", "bit-blast", "tseitin-cnf")(goal)
        if len(compiled) != 1:
            raise AssertionError(f"expected one CNF goal, got {len(compiled)}")
        dimacs = compiled[0].dimacs()
        arguments.dimacs.write_text(dimacs)
        metadata["dimacs_path"] = arguments.dimacs.name
        metadata["dimacs_sha256"] = hashlib.sha256(dimacs.encode()).hexdigest()
        metadata["dimacs_clause_count"] = len(compiled[0])

    started = time.monotonic()
    status = solver.check()
    metadata["solve_seconds"] = time.monotonic() - started
    metadata["status"] = str(status)
    metadata["reason_unknown"] = solver.reason_unknown() if status == unknown else None
    canonical = json.dumps(metadata, sort_keys=True, separators=(",", ":"))
    metadata["mathematical_digest"] = hashlib.sha256(canonical.encode()).hexdigest()

    rendered = json.dumps(metadata, indent=2, sort_keys=True) + "\n"
    if arguments.output is not None:
        arguments.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
