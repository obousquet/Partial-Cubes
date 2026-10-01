#!/usr/bin/env python3
"""Exact finite checks for coordinate products and missing-pattern parents."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from math import comb
from pathlib import Path
import json


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "data/examples/assets"
RECONCILIATION = ASSETS / "maximum_parent_reconciliation_source.json"
ROUND954 = ASSETS / "maximum_parent_round954_source.json"
OUTPUT = ASSETS / "maximum_parent_products_database_verification.json"

EXPECTED_SHA256 = {
    RECONCILIATION.name: "6769f2952b81b90c808deef4a3c859b8617fb1cd3752eda45a68c6c03d020386",
    ROUND954.name: "40aea16f70bf750bd521eee285d189209bddccef35c8995aa3ce6f2d3a7eb756",
}


def submasks(mask: int):
    current = mask
    while True:
        yield current
        if current == 0:
            break
        current = (current - 1) & mask


def shattered(concepts, n: int):
    concepts = set(concepts)
    return {
        support
        for support in range(1 << n)
        if {concept & support for concept in concepts} == set(submasks(support))
    }


def strongly_shattered(concepts, n: int):
    concepts = set(concepts)
    answer = set()
    full = (1 << n) - 1
    for support in range(1 << n):
        outside = full ^ support
        if any(
            all((base | trace) in concepts for trace in submasks(support))
            for base in submasks(outside)
        ):
            answer.add(support)
    return answer


def is_ample(concepts, n: int):
    return shattered(concepts, n) == strongly_shattered(concepts, n)


def missing(concepts, n: int):
    concepts = tuple(concepts)
    result = set()
    for support in range(1, 1 << n):
        traces = {concept & support for concept in concepts}
        for bits in submasks(support):
            if bits not in traces and all(
                any(
                    concept & (support ^ (1 << coordinate))
                    == bits & (support ^ (1 << coordinate))
                    for concept in concepts
                )
                for coordinate in range(n)
                if support >> coordinate & 1
            ):
                result.add((support, bits))
    return result


def dot(first, second):
    return sum(x * y for x, y in zip(first, second))


def strict_linear(rows_, n: int):
    size = len(rows_)
    initial = [
        (tuple(map(F, row)), tuple(F(i == j) for j in range(size)))
        for i, row in enumerate(rows_)
    ]

    def solve(system, dimension):
        compact = {}
        for row, weights in system:
            if not any(row):
                return None, weights
            scale = next(abs(value) for value in row if value)
            compact.setdefault(
                tuple(value / scale for value in row),
                tuple(value / scale for value in weights),
            )
        system = list(compact.items())
        if not system:
            return (F(0),) * dimension, None
        positive = [(row, weights) for row, weights in system if row[-1] > 0]
        negative = [(row, weights) for row, weights in system if row[-1] < 0]
        reduced = [(row[:-1], weights) for row, weights in system if row[-1] == 0]
        for first, first_weights in positive:
            for second, second_weights in negative:
                reduced.append(
                    (
                        tuple(
                            x / first[-1] - y / second[-1]
                            for x, y in zip(first[:-1], second[:-1])
                        ),
                        tuple(
                            x / first[-1] - y / second[-1]
                            for x, y in zip(first_weights, second_weights)
                        ),
                    )
                )
        point, dual = solve(reduced, dimension - 1)
        if point is None:
            return None, dual
        lower = [-dot(row[:-1], point) / row[-1] for row, _ in positive]
        upper = [-dot(row[:-1], point) / row[-1] for row, _ in negative]
        low, high = max(lower, default=None), min(upper, default=None)
        if low is not None and high is not None:
            assert low < high
            last = (low + high) / 2
        elif low is not None:
            last = low + 1
        elif high is not None:
            last = high - 1
        else:
            last = F(0)
        return point + (last,), None

    point, dual = solve(initial, n)
    if point is not None:
        assert all(dot(row, point) > 0 for row in rows_)
    return point, dual


def unit_separator_witnesses(concepts, n: int):
    patterns = missing(concepts, n)
    witnesses = {}
    for concept in concepts:
        rows = [
            tuple(F(1 if concept >> i & 1 else -1) if i == j else F(0) for i in range(n))
            for j in range(n)
        ]
        for support, bits in patterns:
            rows.append(
                tuple(
                    F(-(1 if bits >> i & 1 else -1)) if support >> i & 1 else F(0)
                    for i in range(n)
                )
            )
        witness, _ = strict_linear(rows, n)
        assert witness is not None
        witnesses[concept] = witness
    return patterns, witnesses


def one_coordinate_product(first, second):
    """Local coordinates are (private, shared); output is (private, private, shared)."""
    return {
        (left & 1) | ((right & 1) << 1) | (((left >> 1) & 1) << 2)
        for left in first
        for right in second
        if (left >> 1) & 1 == (right >> 1) & 1
    }


def replay_one_coordinate(counts):
    classes = [
        {word for word in range(4) if mask >> word & 1}
        for mask in range(16)
    ]
    ample = [concepts for concepts in classes if is_ample(concepts, 2)]
    for first, second in product(ample, repeat=2):
        assert is_ample(one_coordinate_product(first, second), 3)
        counts["one_coordinate_ample_pairs"] += 1
    counts["one_coordinate_ample_inputs"] = len(ample)


def fibers(concepts, shared_word: int):
    return {word >> 2 for word in concepts if word & 3 == shared_word}


def replay_tree_products(counts):
    tree = (0, 1, 3)
    universe = tuple(shared | (private << 2) for shared in tree for private in (0, 1))
    classes = []
    for mask in range(1 << len(universe)):
        concepts = {word for i, word in enumerate(universe) if mask >> i & 1}
        if is_ample(concepts, 3):
            classes.append(concepts)
    edges = ((0, 1, 0), (1, 3, 1))
    for first, second in product(classes, repeat=2):
        joined = {
            shared | (left_private << 2) | (right_private << 3)
            for shared in tree
            for left_private in fibers(first, shared)
            for right_private in fibers(second, shared)
        }
        expected = set()
        for shared in tree:
            first_sh = shattered(fibers(first, shared), 1)
            second_sh = shattered(fibers(second, shared), 1)
            for left_support in first_sh:
                for right_support in second_sh:
                    expected.add((left_support << 2) | (right_support << 3))
        for left, right, direction in edges:
            first_edge = fibers(first, left) & fibers(first, right)
            second_edge = fibers(second, left) & fibers(second, right)
            for left_support in shattered(first_edge, 1):
                for right_support in shattered(second_edge, 1):
                    expected.add((1 << direction) | (left_support << 2) | (right_support << 3))
        assert shattered(joined, 4) == expected
        assert is_ample(joined, 4)
        counts["tree_product_pairs"] += 1
    counts["tree_product_inputs"] = len(classes)


def bitword(text: str):
    return sum((character == "1") << i for i, character in enumerate(text))


def replay_opposite_triple(counts):
    lowest = {bitword("000")}
    middle = {bitword(word) for word in ("000", "100", "101", "111")}
    highest = set(range(8)) - {bitword("001")}
    sequence = (lowest, middle, highest, set(range(8)))
    level = {
        word: next(index for index, concepts in enumerate(sequence) if word in concepts)
        for word in range(8)
    }
    visible = set()
    for first in range(8):
        for second in range(8):
            if first in lowest and second in highest:
                visible.add(first | (second << 3) | (bitword("00") << 6))
            if first in middle and second in middle:
                visible.add(first | (second << 3) | (bitword("10") << 6))
            if first in highest and second in lowest:
                visible.add(first | (second << 3) | (bitword("11") << 6))
    assert len(visible) == 30 and is_ample(visible, 8)
    visible_sh = shattered(visible, 8)
    assert max(mask.bit_count() for mask in visible_sh) == 2

    for coordinate in range(8):
        zero = {
            (word & ((1 << coordinate) - 1)) | ((word >> (coordinate + 1)) << coordinate)
            for word in visible
            if not word >> coordinate & 1
        }
        one = {
            (word & ((1 << coordinate) - 1)) | ((word >> (coordinate + 1)) << coordinate)
            for word in visible
            if word >> coordinate & 1
        }
        assert not zero <= one and not one <= zero
        counts["nonnested_coordinates"] += 1

    full = set(range(8))
    punctured = full - {bitword("010")}
    left_path = {bitword(word) for word in ("100", "110", "101", "001")}
    right_path = {bitword(word) for word in ("000", "100", "001", "011")}
    parent = set()
    rank = 2
    for first in range(8):
        for second in range(8):
            p, q = level[first], level[second]
            if p + q <= rank - 3:
                control = full
            elif p + q == rank - 2:
                control = punctured
            elif p + q == rank - 1:
                control = left_path if p in (1, 2) else right_path
            elif p + q == rank:
                selected = {0: bitword("000"), 1: bitword("100"), 2: bitword("110")}.get(p, bitword("001"))
                control = {selected}
            else:
                control = set()
            parent.update(first | (second << 3) | (word << 6) for word in control)
    assert len(parent) == sum(comb(9, i) for i in range(3)) == 46
    assert max(mask.bit_count() for mask in shattered(parent, 9)) == 2
    marked = {
        (word & ((1 << 8) - 1))
        for word in parent
        if not word >> 8 & 1
    }
    assert marked == visible
    counts["opposite_triple_parent_words"] = len(parent)
    counts["opposite_triple_visible_words"] = len(visible)


def replay_large_patterns(counts):
    for mask in range(1, 1 << 8):
        concepts = {word for word in range(8) if mask >> word & 1}
        if not is_ample(concepts, 3):
            continue
        patterns = missing(concepts, 3)
        if sum(support.bit_count() >= 3 for support, _ in patterns) > 2:
            continue
        _, witnesses = unit_separator_witnesses(concepts, 3)
        counts["small_ample_classes_with_unit_separators"] += 1
        counts["small_separator_witnesses"] += len(witnesses)

    down = set(range(8)) - {bitword("111")}
    anchor = bitword("110")
    agreeing = {
        first | (second << 3)
        for first in range(8)
        for second in range(8)
        if (first in down and second == anchor) or (first == anchor and second in down)
    }
    assert len(agreeing) == 13 and is_ample(agreeing, 6)
    patterns, witnesses = unit_separator_witnesses(agreeing, 6)
    large = [(support, bits) for support, bits in patterns if support.bit_count() >= 3]
    small = [(support, bits) for support, bits in patterns if support.bit_count() == 2]
    assert len(patterns) == 11 and len(large) == 2 and len(small) == 9
    assert large[0][0] & large[1][0] == 0
    comparisons = 0
    for concept, witness in witnesses.items():
        for support, bits in patterns:
            value = sum(
                (1 if bits >> i & 1 else -1) * witness[i]
                for i in range(6)
                if support >> i & 1
            )
            assert value < 0
            comparisons += 1
    for reflection in range(1 << 6):
        reflected = {word ^ reflection for word in agreeing}
        assert any((first & second) not in reflected for first in reflected for second in reflected)
    counts["agreeing_block_concepts"] = len(agreeing)
    counts["agreeing_block_missing_patterns"] = len(patterns)
    counts["agreeing_block_separator_comparisons"] = comparisons
    counts["agreeing_block_reflections_rejected"] = 1 << 6


def replay_scope_witness():
    first = {bitword(word) for word in ("00", "10", "11")}
    second = {bitword(word) for word in ("00", "01", "11")}
    intersection = first & second
    assert is_ample(first, 2) and is_ample(second, 2)
    assert intersection == {bitword("00"), bitword("11")}
    assert not is_ample(intersection, 2)
    return sorted(intersection)


def main():
    hashes = {
        path.name: sha256(path.read_bytes()).hexdigest()
        for path in (RECONCILIATION, ROUND954)
    }
    assert hashes == EXPECTED_SHA256
    reconciliation = json.loads(RECONCILIATION.read_text())
    round954 = json.loads(ROUND954.read_text())
    assert reconciliation["source_sha256"] == "c1e8df40a3004200cc83bd4c1c3f5831588568bba301cb767d351ca39ec89f53"
    assert reconciliation["counts"] == {"scoped_audit_located": 126, "unreconciled": 0}
    required = {
        "thm:coordinate-product",
        "rem:coordinate-product-scope",
        "thm:tree-product-ample",
        "thm:path-triple-parent",
        "ex:path-nonperipheral",
        "prop:path-cumulative-parents",
        "cor:monotone-path-product",
        "thm:two-large-patterns",
        "thm:agreeing-large-patterns",
        "ex:agreeing-nonintersection",
    }
    entries = {entry["label"]: entry for entry in reconciliation["statements"]}
    assert required <= entries.keys()
    assert all(entries[label]["reconciliation"] == "scoped_audit_located" for label in required)
    assert round954["status"] == "PASS"
    assert {
        "thm:two-large-patterns",
        "thm:agreeing-large-patterns",
    } <= round954["latex"]["labels"].keys()

    counts = Counter()
    replay_one_coordinate(counts)
    replay_tree_products(counts)
    replay_opposite_triple(counts)
    replay_large_patterns(counts)
    intersection = replay_scope_witness()
    output = {
        "date": "2026-10-01",
        "status": "passed",
        "arithmetic": "exact finite sets and fractions.Fraction Fourier--Motzkin elimination",
        "scope": "Bounded independent replay of one-coordinate and tree product ampleness, the nonample two-coordinate witness, the 30-word opposite-triple class and its explicit 46-word rank-two maximum parent, and unit signed separators for small ample and agreeing-block classes. Universal maximum-parent constructions use database-owned proofs.",
        "source_sha256": hashes,
        "source_statement_labels": sorted(required),
        "counts": dict(counts),
        "matched_path_intersection": intersection,
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"status": "passed", "counts": dict(counts)}))


if __name__ == "__main__":
    main()
