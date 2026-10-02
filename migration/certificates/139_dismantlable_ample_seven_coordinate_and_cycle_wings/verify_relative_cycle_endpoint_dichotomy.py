#!/usr/bin/env python3
"""Self-contained replay of the relative C6 endpoint counterexample."""

from __future__ import annotations

import json
from pathlib import Path


OUTPUT = Path(__file__).with_name("relative_cycle_endpoint_dichotomy.json")


def shadow(family: set[int], dimension: int) -> set[int]:
    shattered = set()
    for support in range(1 << dimension):
        traces = {x & support for x in family}
        submasks = set()
        submask = support
        while True:
            submasks.add(submask)
            if submask == 0:
                break
            submask = (submask - 1) & support
        if traces == submasks:
            shattered.add(support)
    return shattered


def ample(family: set[int], dimension: int) -> bool:
    return len(shadow(family, dimension)) == len(family)


def corner(family: set[int], vertex: int, dimension: int) -> bool:
    return vertex in family and ample(family - {vertex}, dimension)


def main() -> None:
    upper = set(range(7)) | {8, 10, 11}
    lower = {0, 8, 10, 11}
    residual, designated, dimension = upper - lower, 1, 4
    assert ample(upper, dimension)
    assert ample(lower, dimension)
    cycle = [1, 3, 2, 6, 4, 5]
    edges = {
        tuple(sorted((x, x ^ (1 << coordinate))))
        for x in residual
        for coordinate in range(dimension)
        if x ^ (1 << coordinate) in residual
    }
    assert edges == {
        tuple(sorted((cycle[i], cycle[(i + 1) % 6])))
        for i in range(6)
    }

    first = lower | {designated}
    assert {1, 11} <= first and not {3, 9} & first
    assert not ample(first, dimension)
    last = upper - {designated}
    assert {3, 5} <= last and not {1, 7} & last
    assert not ample(last, dimension)
    assert not corner(upper, designated, dimension)

    deletion = [5, 6, 4, 1, 3, 2]
    remaining = set(upper)
    for vertex in deletion:
        assert vertex in remaining - lower and corner(remaining, vertex, dimension)
        remaining.remove(vertex)
        assert ample(remaining, dimension)
    assert remaining == lower

    supports = shadow(upper, dimension) - shadow(lower, dimension)
    isolated = [
        support
        for support in supports
        if all(
            support == other or support & other not in (support, other)
            for other in supports
        )
    ]
    markers = {}
    complement = set(range(1 << dimension)) - upper
    for support in isolated:
        other = ((1 << dimension) - 1) ^ support
        candidates = [
            x
            for x in residual
            if x & support not in {h & support for h in lower}
            and x & other not in {y & other for y in complement}
        ]
        assert len(candidates) == 1
        markers[support] = candidates[0]
    assert supports == {3, 4, 5, 6, 9, 10}
    assert markers == {3: 1, 9: 3, 10: 2}
    available = [
        x
        for x in sorted(residual)
        if ample(lower | {x}, dimension) or corner(upper, x, dimension)
    ]
    assert len(available) >= len(supports) - len(isolated)

    output = {
        "status": "PASS",
        "A": sorted(upper),
        "H": sorted(lower),
        "designated_vertex": designated,
        "exclusive_cycle": cycle,
        "relative_deletion_order": deletion,
        "first_failure": {"pair": [1, 11], "missing_intermediates": [3, 9]},
        "last_failure": {"pair": [3, 5], "missing_intermediates": [1, 7]},
        "relative_supports": sorted(supports),
        "isolated_markers": markers,
        "available_vertices": available,
        "nonisolated_support_count": len(supports) - len(isolated),
        "scope": "Refutes the degree-two vertexwise first-or-last dichotomy only. The pair is relatively peelable and supplies no Damp obstruction or counterexample to arbitrary-cycle relative peeling.",
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n")
    print("PASS: C6 exclusive graph; vertex 1 fails both endpoints; relative peeling replayed")


if __name__ == "__main__":
    main()
