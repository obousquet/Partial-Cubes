#!/usr/bin/env python3
"""Exhaustive small checks for maximum sequences and Boolean substitution."""

from __future__ import annotations

from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from math import comb
from pathlib import Path
import json

from verify_signed_pattern_recognition import (
    literal_pattern,
    minimal_generators,
    patterns,
    shattered,
)


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "data/examples/assets"
RECONCILIATION = ASSETS / "maximum_parent_reconciliation_source.json"
OUTPUT = ASSETS / "boolean_substitution_database_verification.json"
EXPECTED_RECONCILIATION = "6769f2952b81b90c808deef4a3c859b8617fb1cd3752eda45a68c6c03d020386"


def vc_dimension(concepts, n: int):
    return max((support.bit_count() for support in shattered(concepts, n)), default=-1)


def is_maximum(concepts, n: int, rank: int):
    return (
        len(concepts) == sum(comb(n, i) for i in range(rank + 1))
        and vc_dimension(concepts, n) == rank
    )


def is_ample(concepts, n: int):
    return len(concepts) == len(shattered(concepts, n))


def delete_bit(word: int, coordinate: int):
    lower = word & ((1 << coordinate) - 1)
    upper = word >> (coordinate + 1)
    return lower | (upper << coordinate)


def section(concepts, n: int, coordinate: int, value: int):
    return frozenset(
        delete_bit(word, coordinate)
        for word in concepts
        if ((word >> coordinate) & 1) == value
    )


@lru_cache(maxsize=None)
def recursively_peripheral(concepts_tuple, n: int):
    concepts = frozenset(concepts_tuple)
    if not concepts:
        return True
    if n == 0:
        return concepts == {0}
    for coordinate in range(n):
        zero = section(concepts, n, coordinate, 0)
        one = section(concepts, n, coordinate, 1)
        if not (zero <= one or one <= zero):
            continue
        if recursively_peripheral(tuple(sorted(zero)), n - 1) and recursively_peripheral(tuple(sorted(one)), n - 1):
            return True
    return False


def substitute(outer, test, nx: int, ny: int):
    answer = set()
    for x in range(1 << nx):
        for y in range(1 << ny):
            old_value = 0 if y in test else 1
            if (x | (old_value << nx)) in outer:
                answer.add(x | (y << nx))
    return frozenset(answer)


def shifted_generators(concepts, n: int, offset: int = 0):
    return {
        frozenset((coordinate + offset, value) for coordinate, value in literal_pattern(pattern))
        for pattern in patterns(concepts, n)
    }


def ideal_product(left, right):
    return {a | b for a in left for b in right}


def predicted_shattering(c0, c1, test, nx: int, ny: int):
    full_test_supports = range(1 << ny)
    test_complement = frozenset(set(range(1 << ny)) - set(test))
    answer = set()
    for u in shattered(c0 & c1, nx):
        answer.update(u | (v << nx) for v in full_test_supports)
    for u in shattered(c0, nx):
        answer.update(u | (v << nx) for v in shattered(test, ny))
    for u in shattered(c1, nx):
        answer.update(u | (v << nx) for v in shattered(test_complement, ny))
    return answer


def replay_maximum_sequences(counts: Counter):
    n = 3
    universe = range(1 << n)
    maximum_by_rank = {}
    for rank in range(n + 1):
        size = sum(comb(n, i) for i in range(rank + 1))
        classes = [
            frozenset(candidate)
            for candidate in combinations(universe, size)
            if is_maximum(candidate, n, rank)
        ]
        maximum_by_rank[rank] = classes
        counts["maximum_q3_classes"] += len(classes)

    for rank, classes in maximum_by_rank.items():
        for cls in classes:
            if rank > 0:
                assert any(candidate < cls for candidate in maximum_by_rank[rank - 1])
                counts["maximum_predecessors_found"] += 1
            if rank < n:
                assert any(cls < candidate for candidate in maximum_by_rank[rank + 1])
                counts["maximum_successors_found"] += 1

    # Hamming-weight chains are explicit maximum sequences; their level sums
    # are Hamming balls on the disjoint union.
    for total_rank in range(6):
        level_sum = {
            x | (y << 2)
            for x in range(4)
            for y in range(8)
            if x.bit_count() + y.bit_count() <= total_rank
        }
        assert is_maximum(level_sum, 5, total_rank)
        counts["level_sum_ranks"] += 1

    # Verify every one-step delay of the weight chain on Q2.
    base = {
        rank: frozenset(word for word in range(4) if word.bit_count() <= rank)
        for rank in range(-1, 4)
    }
    base[-1] = frozenset()
    base[3] = frozenset(range(4))
    for threshold in range(3):
        delayed = []
        for rank in range(4):
            if rank <= threshold:
                cls = set(base[rank]) | {word | 4 for word in base[rank - 1]}
            else:
                cls = set(base[rank - 1]) | {word | 4 for word in base[rank]}
            cls = frozenset(cls)
            assert is_maximum(cls, 3, rank)
            delayed.append(cls)
            counts["delayed_sequence_levels"] += 1
        assert all(delayed[i] <= delayed[i + 1] for i in range(3))


def replay_closure_basics(counts: Counter):
    ample_q2 = []
    for family_mask in range(1 << 4):
        concepts = frozenset(word for word in range(4) if family_mask >> word & 1)
        if is_ample(concepts, 2):
            ample_q2.append(concepts)
            complement = frozenset(set(range(4)) - set(concepts))
            assert is_ample(complement, 2)
            counts["ample_complements_q2"] += 1
    for left in ample_q2:
        for right in ample_q2:
            product_class = frozenset(x | (y << 2) for x in left for y in right)
            assert is_ample(product_class, 4)
            counts["ample_products_q2"] += 1


def replay_substitution(counts: Counter):
    nx = ny = 2
    tests = [
        frozenset(word for word in range(4) if family_mask >> word & 1)
        for family_mask in range(1, 15)
    ]
    for outer_mask in range(1, 1 << 8):
        outer = frozenset(word for word in range(8) if outer_mask >> word & 1)
        c0 = section(outer, 3, 2, 0)
        c1 = section(outer, 3, 2, 1)
        if not c0 or not c1 or c0 == c1:
            continue
        p = c0 | c1
        for test in tests:
            test_complement = frozenset(set(range(4)) - set(test))
            result = substitute(outer, test, nx, ny)

            ample_inputs = is_ample(outer, 3) and is_ample(test, 2)
            assert is_ample(result, 4) == ample_inputs
            counts["substitution_ample_equivalences"] += 1

            rper_inputs = recursively_peripheral(tuple(sorted(outer)), 3) and recursively_peripheral(tuple(sorted(test)), 2)
            assert recursively_peripheral(tuple(sorted(result)), 4) == rper_inputs
            counts["substitution_recursive_equivalences"] += 1

            if ample_inputs:
                assert shattered(result, 4) == predicted_shattering(c0, c1, test, nx, ny)
                counts["substitution_shattering_formulas"] += 1

            q_result = shifted_generators(result, 4)
            q_p = shifted_generators(p, nx)
            q_c0 = shifted_generators(c0, nx)
            q_c1 = shifted_generators(c1, nx)
            q_test = shifted_generators(test, ny, nx)
            q_test_complement = shifted_generators(test_complement, ny, nx)
            predicted_ideal = minimal_generators(
                q_p | ideal_product(q_c0, q_test_complement) | ideal_product(q_c1, q_test)
            )
            assert q_result == predicted_ideal
            counts["substitution_signed_ideal_formulas"] += 1


def main():
    assert sha256(RECONCILIATION.read_bytes()).hexdigest() == EXPECTED_RECONCILIATION
    reconciliation = json.loads(RECONCILIATION.read_text())
    required = {
        "def:peripheral",
        "lem:peripheral-closures",
        "lem:sequences",
        "thm:substitution",
        "cor:substitution-factor",
    }
    entries = {entry["label"]: entry for entry in reconciliation["statements"]}
    assert required <= entries.keys()
    assert all(entries[label]["reconciliation"] == "scoped_audit_located" for label in required)

    counts = Counter()
    replay_maximum_sequences(counts)
    replay_closure_basics(counts)
    replay_substitution(counts)
    output = {
        "date": "2026-10-01",
        "status": "passed",
        "arithmetic": "exhaustive finite binary families, exact shattering sets, and squarefree monomial generators",
        "scope": "Checks predecessors and successors for every maximum class on Q3, all ranks of a Q2-by-Q3 level sum, every one-step delay of the Q2 weight chain, all Q2 ample complements and ordered products, and every qualified substitution with a three-coordinate outer class and two-coordinate nonconstant test. The substitution replay checks Ample and recursive-peripheral equivalence, exact shattering formulas, and signed-ideal identities. General maximum-parent constructions use database-owned proofs.",
        "source_sha256": {RECONCILIATION.name: EXPECTED_RECONCILIATION},
        "source_statement_labels": sorted(required),
        "counts": dict(counts),
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"status": "passed", "counts": dict(counts)}))


if __name__ == "__main__":
    main()
