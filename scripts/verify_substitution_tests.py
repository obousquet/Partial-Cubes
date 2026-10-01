#!/usr/bin/env python3
"""Exact finite checks for conjunction and universal substitution tests."""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import product
from math import comb
from pathlib import Path
import json

from verify_signed_pattern_recognition import patterns, shattered


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "data/examples/assets"
RECONCILIATION = ASSETS / "maximum_parent_reconciliation_source.json"
OUTPUT = ASSETS / "substitution_tests_database_verification.json"
EXPECTED_RECONCILIATION = "6769f2952b81b90c808deef4a3c859b8617fb1cd3752eda45a68c6c03d020386"


def vc_dimension(concepts, n: int):
    return max((support.bit_count() for support in shattered(concepts, n)), default=-1)


def is_maximum(concepts, n: int, rank: int):
    return len(concepts) == sum(comb(n, j) for j in range(rank + 1)) and vc_dimension(concepts, n) == rank


def is_ample(concepts, n: int):
    return len(concepts) == len(shattered(concepts, n))


def section(concepts, n: int, coordinate: int, value: int):
    low = (1 << coordinate) - 1
    return frozenset(
        (word & low) | ((word >> 1) & ~low)
        for word in concepts
        if ((word >> coordinate) & 1) == value
    )


def projection(concepts, coordinate: int):
    low = (1 << coordinate) - 1
    return frozenset((word & low) | ((word >> 1) & ~low) for word in concepts)


def full_reduction(concepts, n: int, coordinate: int):
    return section(concepts, n, coordinate, 0) & section(concepts, n, coordinate, 1)


def complement(concepts, n: int):
    return frozenset(set(range(1 << n)) - set(concepts))


def conjunction_substitution(concepts, x_count: int, k: int):
    answer = set()
    old_bit = 1 << x_count
    all_one = (1 << k) - 1
    for old in concepts:
        x = old & (old_bit - 1)
        old_value = bool(old & old_bit)
        for z in range(1 << k):
            if (z == all_one) == old_value:
                answer.add(x | (z << x_count))
    return frozenset(answer)


def replay_conjunction(counts: Counter):
    x_count = 2
    n = x_count + 1
    for family_mask in range(1, 1 << (1 << n)):
        concepts = frozenset(w for w in range(1 << n) if family_mask >> w & 1)
        a = section(concepts, n, x_count, 0)
        b = section(concepts, n, x_count, 1)
        for k in (2, 3):
            output = conjunction_substitution(concepts, x_count, k)
            actual = shattered(output, x_count + k)
            predicted = set()
            for sx in range(1 << x_count):
                for tz in range(1 << k):
                    if tz == 0:
                        ok = sx in shattered(a | b, x_count)
                    elif tz != (1 << k) - 1:
                        ok = sx in shattered(a, x_count)
                    else:
                        old_support = sx | (1 << x_count)
                        ok = old_support in shattered(concepts, n)
                    if ok:
                        predicted.add(sx | (tz << x_count))
            assert actual == predicted
            counts["conjunction_shattering_instances"] += 1

            rank = vc_dimension(concepts, n)
            if rank >= 1 and is_maximum(concepts, n, rank):
                assert vc_dimension(output, x_count + k) == rank + k - 1
                counts["maximum_conjunction_rank_instances"] += 1

    # Exact top-rank boundary behavior for every missing word and k=2,3.
    for k in (2, 3):
        full = frozenset(range(1 << n))
        assert conjunction_substitution(full, x_count, k) == frozenset(range(1 << (x_count + k)))
        for missing in range(1 << n):
            concepts = full - {missing}
            output = conjunction_substitution(concepts, x_count, k)
            if missing >> x_count:
                assert is_maximum(output, x_count + k, n + k - 2)
            else:
                assert not is_maximum(output, x_count + k, n + k - 2)
                missing_patterns = patterns(output, x_count + k)
                visible = missing & ((1 << x_count) - 1)
                for new_coordinate in range(k):
                    support = ((1 << x_count) - 1) | (1 << (x_count + new_coordinate))
                    bits = visible
                    assert (support, bits) in missing_patterns
            counts["conjunction_top_rank_boundaries"] += 1


def monotone_truth_tables(k: int):
    for mask in range(1, (1 << (1 << k)) - 1):
        values = tuple((mask >> word) & 1 for word in range(1 << k))
        if all(not values[x] or values[y] for x in range(1 << k) for y in range(1 << k) if x & ~y == 0):
            yield values


def inclusion_minimal(words):
    return tuple(w for w in words if not any(v != w and v & ~w == 0 for v in words))


def three_section(test, k: int):
    test = set(test)
    answer = set()
    for y in range(1 << k):
        if y in test:
            answer.add(y)  # controls 00
        answer.add(y | (1 << k))  # controls 01 in (u,v) order
        if y not in test:
            answer.add(y | (3 << k))  # controls 11
    return frozenset(answer)


def guards_hold(z, u, v, minimal_ones, minimal_zeros):
    return (
        u - v < 0
        and all(sum(z[j] for j in range(len(z)) if support >> j & 1) - v < 0 for support in minimal_ones)
        and all(u - sum(z[j] for j in range(len(z)) if support >> j & 1) < 0 for support in minimal_zeros)
    )


def replay_monotone_cones(counts: Counter):
    k = 3
    all_bits = (1 << k) - 1
    for values in monotone_truth_tables(k):
        zero_set = frozenset(w for w, value in enumerate(values) if value == 0)
        one_set = set(range(1 << k)) - set(zero_set)
        minimal_ones = inclusion_minimal(one_set)
        # Zero-coordinate sets of words in the zero fiber, ordered by inclusion.
        zero_supports = {all_bits ^ w for w in zero_set}
        minimal_zeros = inclusion_minimal(zero_supports)
        h = three_section(zero_set, k)
        assert len(h) == 1 << (k + 1)
        assert is_ample(h, k + 2)

        for y in zero_set:
            z = tuple(1 if y >> j & 1 else -(k + 1) for j in range(k))
            b_sums = [sum(z[j] for j in range(k) if support >> j & 1) for support in minimal_zeros]
            u, v = min([-2] + [value - 1 for value in b_sums]), -1
            assert guards_hold(z, u, v, minimal_ones, minimal_zeros)
            assert u < 0 and v < 0
        for y in one_set:
            z = tuple(k + 1 if y >> j & 1 else -1 for j in range(k))
            a_sums = [sum(z[j] for j in range(k) if support >> j & 1) for support in minimal_ones]
            u, v = 1, max([2] + [value + 1 for value in a_sums])
            assert guards_hold(z, u, v, minimal_ones, minimal_zeros)
            assert u > 0 and v > 0
        for y in range(1 << k):
            z = tuple(1 if y >> j & 1 else -1 for j in range(k))
            assert guards_hold(z, -(k + 1), k + 1, minimal_ones, minimal_zeros)
        counts["nonconstant_monotone_q3_tests"] += 1
        counts["monotone_three_section_words"] += len(h)


def replay_nonmonotone_path(counts: Counter):
    test = frozenset({0b000, 0b001, 0b011, 0b111})  # displayed order (1,2,3), bit 0 first
    test_complement = complement(test, 3)
    values = (-3, -2, -1, 1, 2, 3)
    first = set()
    second = set()
    for z1, z2, z3 in product(values, repeat=3):
        word = (z1 > 0) | ((z2 > 0) << 1) | ((z3 > 0) << 2)
        if z1 > z2 > z3:
            first.add(word)
        if z3 > -z2 > z1:
            second.add(word)
    assert first == set(test)
    assert second == set(test_complement)
    assert is_ample(test, 3)
    for reflection in range(8):
        reflected = {word ^ reflection for word in test}
        assert not all(y in reflected for x in reflected for y in range(8) if y & ~x == 0)
    counts["nonmonotone_path_words"] = len(test)
    counts["nonmonotone_reflections_rejected"] = 8


def replay_universal_identities(counts: Counter):
    for k in range(1, 4):
        cube = set(range(1 << k))
        for family_mask in range(1 << (1 << k)):
            test = frozenset(w for w in cube if family_mask >> w & 1)
            h = three_section(test, k)
            transformed = {
                (word & ((1 << k) - 1))
                | ((((1 - ((word >> (k + 1)) & 1)) << 0) | ((1 - ((word >> k) & 1)) << 1)) << k)
                for word in h
            }
            assert transformed == three_section(complement(test, k), k)
            for coordinate in range(k):
                for value in (0, 1):
                    assert section(h, k + 2, coordinate, value) == three_section(section(test, k, coordinate, value), k - 1)
                assert full_reduction(test, k, coordinate) == complement(
                    projection(complement(test, k), coordinate), k - 1
                )
            counts["universal_test_identity_families"] += 1


def main():
    assert sha256(RECONCILIATION.read_bytes()).hexdigest() == EXPECTED_RECONCILIATION
    reconciliation = json.loads(RECONCILIATION.read_text())
    required = {
        "thm:and-substitution",
        "cor:read-once-substitution",
        "prop:and-shattering",
        "prop:substitution-test",
        "thm:monotone-substitution",
        "prop:complementary-realizations",
        "def:universal-tests",
        "prop:universal-test-closure",
        "prop:geometric-test-closure",
    }
    entries = {entry["label"]: entry for entry in reconciliation["statements"]}
    assert required <= entries.keys()
    assert all(entries[label]["reconciliation"] == "scoped_audit_located" for label in required)

    counts = Counter()
    replay_conjunction(counts)
    replay_monotone_cones(counts)
    replay_nonmonotone_path(counts)
    replay_universal_identities(counts)
    output = {
        "date": "2026-10-01",
        "status": "passed",
        "arithmetic": "exhaustive finite Boolean families, exact shattering sets, integer cone witnesses, and set identities",
        "scope": "Checks the conjunction shattering law for every nonempty class on three old coordinates at arities two and three, all maximum-rank consequences and top-rank boundaries there, every nonconstant monotone test on Q3 with explicit cone witnesses, the nonmonotone path realization, and the three-section section/complement/reduction identities for every test through Q3. General parent constructions use the database-owned proofs.",
        "source_sha256": {RECONCILIATION.name: EXPECTED_RECONCILIATION},
        "source_statement_labels": sorted(required),
        "counts": dict(counts),
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"status": "passed", "counts": dict(counts)}))


if __name__ == "__main__":
    main()
