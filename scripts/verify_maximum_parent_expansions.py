#!/usr/bin/env python3
"""Exact examples and bounded checks for maximum-parent expansions."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from math import comb
from pathlib import Path
import json

from verify_maximum_parent_products import (
    RECONCILIATION,
    bitword,
    is_ample,
    missing,
    shattered,
)


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "data/examples/assets"
OUTPUT = ASSETS / "maximum_parent_expansions_database_verification.json"
EXPECTED_RECONCILIATION = "6769f2952b81b90c808deef4a3c859b8617fb1cd3752eda45a68c6c03d020386"


def is_maximum(concepts, n: int, rank: int):
    return (
        len(concepts) == sum(comb(n, i) for i in range(rank + 1))
        and max((support.bit_count() for support in shattered(concepts, n)), default=-1) == rank
    )


def is_learning_space(concepts, n: int):
    concepts = set(concepts)
    if not concepts or 0 not in concepts:
        return False
    if any((first | second) not in concepts for first in concepts for second in concepts):
        return False
    return all(
        concept == 0
        or any(concept ^ (1 << i) in concepts for i in range(n) if concept >> i & 1)
        for concept in concepts
    )


def replay_square_parent(counts):
    sequence = (
        {0},
        {word for word in range(4) if word.bit_count() <= 1},
        set(range(4)),
    )
    level = {
        word: next(index for index, concepts in enumerate(sequence) if word in concepts)
        for word in range(4)
    }
    full = set(range(16))
    tree_x = {bitword(word) for word in ("0000", "0100", "0110", "1110", "0111")}
    tree_y = {bitword(word) for word in ("0100", "1100", "0110", "0010", "0111")}
    kx, ky = full - tree_x, full - tree_y
    k = full - {bitword("0111")}
    stars = {
        (2, 0): {bitword(word) for word in ("1101", "1100", "0101", "1001", "1111")},
        (1, 1): {bitword(word) for word in ("1001", "1000", "1011", "0001", "1101")},
        (0, 2): {bitword(word) for word in ("0001", "0000", "1001", "0101", "0011")},
        "other": {bitword(word) for word in ("1011", "1010", "0011", "1111", "1001")},
    }
    singleton = bitword("1001")
    rank = 3
    parent = set()
    for first, second in product(range(4), repeat=2):
        first_level, second_level = level[first], level[second]
        total = first_level + second_level
        if total <= rank - 4:
            control = full
        elif total == rank - 3:
            control = k
        elif total == rank - 2:
            control = kx if first_level >= 1 else ky
        elif total == rank - 1:
            control = stars.get((first_level, second_level), stars["other"])
        elif total == rank:
            control = {singleton}
        else:
            control = set()
        parent.update(first | (second << 2) | (word << 4) for word in control)
    assert is_maximum(parent, 8, 3) and len(parent) == 93
    visible = {word & 63 for word in parent if not word & (3 << 6)}
    expected = set()
    for first, second in product(range(4), repeat=2):
        if first in sequence[0] and second in sequence[2]:
            expected.add(first | (second << 2) | (bitword("00") << 4))
        if first in sequence[1] and second in sequence[1]:
            expected.add(first | (second << 2) | (bitword("10") << 4))
        if first in sequence[2] and second in sequence[0]:
            expected.add(first | (second << 2) | (bitword("11") << 4))
        if first in sequence[0] and second in sequence[0]:
            expected.add(first | (second << 2) | (bitword("01") << 4))
    assert visible == expected and len(visible) == 18
    counts["square_parent_words"] = len(parent)
    counts["square_section_words"] = len(visible)


def replay_monotone_optimum(counts):
    section = {0, 1, 2, 3, 4, 5, 6, 8}
    upper = {0, 1, 2, 4, 8, 9, 10, 12}
    parent = section | {16 + word for word in upper}
    assert is_ample(section, 4)
    patterns = missing(section, 4)
    assert max(support.bit_count() for support in shattered(section, 4)) == 2
    assert min(support.bit_count() for support, _ in patterns) == 2
    assert is_maximum(parent, 5, 2)
    assert {word for word in parent if not word & 16} == section
    counts["monotone_section_words"] = len(section)
    counts["monotone_parent_words"] = len(parent)
    counts["monotone_missing_patterns"] = len(patterns)


def replay_learning_spaces(counts):
    n = 3
    for mask in range(1 << (1 << n)):
        concepts = {word for word in range(1 << n) if mask >> word & 1}
        if not is_learning_space(concepts, n):
            continue
        assert is_ample(concepts, n)
        patterns = missing(concepts, n)
        assert all(bits.bit_count() == 1 for _, bits in patterns)
        epsilon = F(1, n + 1)
        magnitude = 3 * (n + 1)
        for concept in concepts:
            order = []
            current = concept
            while current:
                coordinate = next(i for i in range(n) if current >> i & 1 and current ^ (1 << i) in concepts)
                order.append(coordinate)
                current ^= 1 << coordinate
            order.reverse()
            witness = [F(-1)] * n
            for index, coordinate in enumerate(order):
                witness[coordinate] = F(magnitude ** (len(order) - index - 1))
            for support, bits in patterns:
                root = bits.bit_length() - 1
                prerequisites = [i for i in range(n) if support >> i & 1 and i != root]
                value = witness[root] - epsilon * sum(witness[i] for i in prerequisites)
                assert value < 0
                counts["learning_inequality_checks"] += 1
            counts["learning_witnesses"] += 1
        counts["learning_spaces_q3"] += 1


def replay_nested_learning(counts):
    def allowed(word):
        c_rule = not word & 4 or bool(word & 3)
        d_rule = not word & 8 or bool(word & 6)
        return c_rule and d_rule

    base = {word for word in range(16) if allowed(word)}
    site = {0, 1, 3, 7, 15}
    assert site <= base and is_learning_space(base, 4) and is_learning_space(site, 4)
    reflected_expansion = site | {word | 16 for word in base}
    assert is_learning_space(reflected_expansion, 5)
    assert 1 in site and 7 in site and 5 in base and 5 not in site
    complement = set(range(16)) - base
    assert len(complement) == 4 and all(not word & 2 for word in complement)
    site_complement = set(range(16)) - site
    face_trace = {word for word in site_complement if not word & 2}
    assert len(face_trace) == 6
    counts["nested_learning_base_words"] = len(base)
    counts["nested_learning_site_words"] = len(site)
    counts["nested_learning_expansion_words"] = len(reflected_expansion)


def main():
    assert sha256(RECONCILIATION.read_bytes()).hexdigest() == EXPECTED_RECONCILIATION
    reconciliation = json.loads(RECONCILIATION.read_text())
    required = {
        "thm:square-triple-parent",
        "thm:block-gluing",
        "cor:face-expansion",
        "cor:dual-convex-expansion",
        "lem:strict-realization",
        "cor:convex-realization",
        "thm:monotone-parameters",
        "ex:monotone-optimal-parent",
        "thm:learning-maxmin",
        "cor:nested-learning",
    }
    entries = {entry["label"]: entry for entry in reconciliation["statements"]}
    assert required <= entries.keys()
    assert all(entries[label]["reconciliation"] == "scoped_audit_located" for label in required)
    counts = Counter()
    replay_square_parent(counts)
    replay_monotone_optimum(counts)
    replay_learning_spaces(counts)
    replay_nested_learning(counts)
    output = {
        "date": "2026-10-01",
        "status": "passed",
        "arithmetic": "exact finite sets and fractions.Fraction",
        "scope": "Independent construction of the rank-three square parent and optimal downset parent, exhaustive prerequisite-inequality checks for all learning spaces on three coordinates, and the nested nonconvex learning-space expansion. Universal gluing and realization theorems use database-owned proofs.",
        "source_sha256": {RECONCILIATION.name: EXPECTED_RECONCILIATION},
        "source_statement_labels": sorted(required),
        "counts": dict(counts),
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"status": "passed", "counts": dict(counts)}))


if __name__ == "__main__":
    main()
