#!/usr/bin/env python3
"""Exact checks for nested geometric substitution and the corner deletion."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json

from verify_signed_pattern_recognition import shattered


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "data/examples/assets"
RECONCILIATION = ASSETS / "maximum_parent_reconciliation_source.json"
OUTPUT = ASSETS / "geometric_substitution_database_verification.json"
EXPECTED_RECONCILIATION = "6769f2952b81b90c808deef4a3c859b8617fb1cd3752eda45a68c6c03d020386"


def vc_dimension(concepts, n: int):
    return max((support.bit_count() for support in shattered(concepts, n)), default=-1)


def is_ample(concepts, n: int):
    return len(concepts) == len(shattered(concepts, n))


def nested_outer(base, upper, x_count: int):
    return frozenset(base | {x | (1 << x_count) for x in upper})


def nested_substitution(base, upper, test, x_count: int, y_count: int):
    answer = set()
    for x in base:
        fibers = range(1 << y_count) if x in upper else test
        answer.update(x | (y << x_count) for y in fibers)
    return frozenset(answer)


def replay_nested_ampleness(counts: Counter):
    x_count = y_count = 2
    x_cube = range(1 << x_count)
    y_cube = range(1 << y_count)
    tests = [
        frozenset(y for y in y_cube if mask >> y & 1)
        for mask in range(1, (1 << (1 << y_count)) - 1)
    ]
    for base_mask in range(1, 1 << (1 << x_count)):
        base = frozenset(x for x in x_cube if base_mask >> x & 1)
        for upper_mask in range(1, 1 << (1 << x_count)):
            upper = frozenset(x for x in x_cube if upper_mask >> x & 1)
            if not upper < base:
                continue
            outer = nested_outer(base, upper, x_count)
            if not is_ample(outer, x_count + 1):
                continue
            for test in tests:
                if not is_ample(test, y_count):
                    continue
                output = nested_substitution(base, upper, test, x_count, y_count)
                assert is_ample(output, x_count + y_count)
                # A test word outside T recovers the upper shore; one in T
                # recovers the base shore, and x in B\S recovers T.
                inside = next(iter(test))
                outside = next(iter(set(y_cube) - set(test)))
                section_inside = {w & 3 for w in output if w >> x_count == inside}
                section_outside = {w & 3 for w in output if w >> x_count == outside}
                assert section_inside == set(base)
                assert section_outside == set(upper)
                x = next(iter(set(base) - set(upper)))
                recovered_test = {w >> x_count for w in output if w & 3 == x}
                assert recovered_test == set(test)
                counts["nested_ample_substitutions"] += 1


E26 = frozenset({
    0, 1, 2, 3, 16, 17, 18, 19, 20, 22, 24, 25, 28,
    32, 33, 34, 35, 36, 37, 40, 42, 44, 48, 52, 56, 60,
})
D27 = E26 | {49}


def guard_vectors():
    d = Fraction(1, 100)
    guards = []
    for triple in combinations(range(4), 3):
        vector = [Fraction(0) for _ in range(6)]
        for coordinate in triple:
            vector[coordinate] = 1
        guards.append(tuple(vector))
    guards.extend([
        (0, 1, 0, 0, d, 2 * d),
        (1, 1, 0, 0, d, 2 * d),
        (0, 0, 1, 0, -d, -d),
        (0, 0, 0, 1, -d, -3 * d),
        (0, 0, 1, 1, -d, -2 * d),
        (1, 0, 1, 0, d, 0),
        (1, 0, 1, 0, 0, -d),
        (0, 1, 0, 1, d, 0),
        (0, 1, 0, 1, 0, -d),
        (0, 1, 1, 0, -d, 0),
        (0, 1, 1, 0, 0, d),
        (1, 0, 0, 1, -d, 0),
        (1, 0, 0, 1, 0, d),
    ])
    return tuple(tuple(Fraction(c) for c in guard) for guard in guards)


RHO = {
    (0, 0): Fraction(1), (0, 1): Fraction(1), (0, 2): Fraction(1), (0, 3): Fraction(1),
    (1, 0): Fraction(1), (1, 1): Fraction(1), (1, 2): Fraction(1), (1, 3): Fraction(1),
    (2, 0): Fraction(1), (2, 1): Fraction(1), (2, 2): Fraction(3),
    (3, 0): Fraction(1), (3, 1): Fraction(1), (3, 2): Fraction(3),
    (4, 1): Fraction(2), (4, 2): Fraction(1, 2), (4, 3): Fraction(1),
    (8, 1): Fraction(4), (8, 2): Fraction(1), (8, 3): Fraction(1),
    (12, 1): Fraction(4), (12, 2): Fraction(1, 2), (12, 3): Fraction(1),
    (5, 2): Fraction(1, 2), (10, 2): Fraction(5, 2),
    (6, 1): Fraction(3, 2), (9, 1): Fraction(4),
}


def replay_corner_deletion(counts: Counter):
    assert len(E26) == 26 and len(D27) == 27 and D27 - E26 == {49}
    assert is_ample(E26, 6) and is_ample(D27, 6)
    assert vc_dimension(E26, 6) == vc_dimension(D27, 6) == 3

    w = 49
    cube_supports = []
    for support in range(1 << 6):
        vertices = {w ^ sub for sub in range(1 << 6) if sub & ~support == 0}
        if vertices <= D27:
            cube_supports.append(support)
    maximal = {s for s in cube_supports if not any(s != t and s & ~t == 0 for t in cube_supports)}
    assert maximal == {49}  # free coordinates {a1,i,j}

    guards = guard_vectors()
    eps = Fraction(1, 10000)
    for word in D27:
        ordinary = word & 15
        controls = (word >> 4) & 3
        rho = RHO[(ordinary, controls)]
        point = [eps if ordinary >> j & 1 else Fraction(-1) for j in range(4)]
        point.append(rho if controls & 1 else -rho)
        point.append(Fraction(1) if controls & 2 else Fraction(-1))
        assert all(Fraction(-2) < point[j] < 1 for j in range(4))
        assert -5 < point[4] < 5 and -2 < point[5] < 2
        assert all(sum(a * z for a, z in zip(guard, point)) < 0 for guard in guards)
        counts["corner_extension_strict_witnesses"] += 1

    for word in set(range(64)) - set(D27):
        signs = tuple(1 if word >> j & 1 else -1 for j in range(6))
        excluding = [guard for guard in guards if all(a * sign >= 0 for a, sign in zip(guard, signs))]
        assert excluding
        counts["corner_extension_forbidden_orthants"] += 1

    counts["corner_deletion_unique_maximal_cubes"] = len(maximal)


def main():
    assert sha256(RECONCILIATION.read_bytes()).hexdigest() == EXPECTED_RECONCILIATION
    reconciliation = json.loads(RECONCILIATION.read_text())
    required = {
        "thm:nested-geometric-substitution",
        "prop:realizable-complement",
        "prop:geometric-corner-deletion",
        "rem:geometric-nonrealizability-not-nonmembership",
    }
    entries = {entry["label"]: entry for entry in reconciliation["statements"]}
    assert required <= entries.keys()
    assert all(entries[label]["reconciliation"] == "scoped_audit_located" for label in required)

    counts = Counter()
    replay_nested_ampleness(counts)
    replay_corner_deletion(counts)
    output = {
        "date": "2026-10-01",
        "status": "passed",
        "arithmetic": "exhaustive small nested substitutions and exact rational inequalities",
        "scope": "Checks every qualified nested ample substitution on two outer and two test coordinates; verifies the E26/D27 sizes, ampleness, VC dimensions, unique corner cube, all 27 rational realization witnesses, and an explicit excluding guard for every forbidden orthant. The general separator-cone, polar-duality, and ratio-cycle proofs are database-owned analytic arguments.",
        "source_sha256": {RECONCILIATION.name: EXPECTED_RECONCILIATION},
        "source_statement_labels": sorted(required),
        "counts": dict(counts),
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"status": "passed", "counts": dict(counts)}))


if __name__ == "__main__":
    main()
