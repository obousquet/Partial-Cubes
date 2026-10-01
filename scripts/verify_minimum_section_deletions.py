#!/usr/bin/env python3
"""Exact examples for minimum-section deletion and downset filtering."""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from math import comb
from pathlib import Path
import json

from verify_maximum_parent_products import RECONCILIATION, bitword, is_ample, shattered


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "data/examples/assets"
OUTPUT = ASSETS / "minimum_section_deletions_database_verification.json"
EXPECTED_RECONCILIATION = "6769f2952b81b90c808deef4a3c859b8617fb1cd3752eda45a68c6c03d020386"


def is_maximum(concepts, n: int, rank: int):
    return len(concepts) == sum(comb(n, i) for i in range(rank + 1)) and max(
        support.bit_count() for support in shattered(concepts, n)
    ) == rank


def section(concepts, block, pattern):
    block_mask = sum(1 << i for i in block)
    pattern_word = sum(int(value) << i for i, value in zip(block, pattern))
    return {word & ~block_mask for word in concepts if word & block_mask == pattern_word}


def reduction(concepts, n: int, block):
    block_mask = sum(1 << i for i in block)
    outside = ((1 << n) - 1) ^ block_mask
    return {
        word
        for word in range(1 << n)
        if not word & block_mask
        and all((word | assignment) in concepts for assignment in range(1 << n) if not assignment & outside)
    }


def replay_hamming_ball(counts):
    parent = {word for word in range(32) if word.bit_count() <= 3}
    block = (0, 1)
    minimum = section(parent, block, (1, 1))
    assert len(minimum) == 4 and minimum == reduction(parent, 5, block)
    for pattern in ((0, 0), (1, 0), (0, 1)):
        assert minimum <= section(parent, block, pattern)
    deleted_words = {word for word in parent if word & 3 == 3}
    residual = parent - deleted_words
    assert len(parent) == 26 and len(residual) == 22 and is_ample(residual, 5)
    assert 3 + 1 - len(block) == 2
    counts["hamming_parent_words"] = len(parent)
    counts["hamming_deleted_words"] = len(deleted_words)
    counts["hamming_residual_words"] = len(residual)


def replay_forest_example(counts):
    path = {bitword("1" * prefix + "0" * (7 - prefix)) for prefix in range(8)}
    maximum = set(range(128)) - path
    assert is_maximum(maximum, 7, 5) and len(maximum) == 120
    specifications = [
        ((2, 3, 4, 5, 6), "00000", "0100000"),
        ((0, 3, 4, 5, 6), "10000", "1010000"),
        ((0, 1, 2, 3, 4), "11111", "1111101"),
    ]
    deleted = set()
    patterns = []
    for block, pattern_text, word_text in specifications:
        pattern = tuple(map(int, pattern_text))
        current = section(maximum, block, pattern)
        assert current == reduction(maximum, 7, block) and len(current) == 1
        word = bitword(word_text)
        assert word in maximum
        assert word & ~sum(1 << i for i in block) in current
        deleted.add(word)
        patterns.append((block, pattern))
    residual = maximum - deleted
    assert len(residual) == 117 and is_ample(residual, 7)
    expected = {
        support
        for support in range(1 << 7)
        if support.bit_count() <= 5
        and all(sum(1 << i for i in block) & support != sum(1 << i for i in block) for block, _ in patterns)
    }
    assert shattered(residual, 7) == expected
    disagreements = set()
    for i, (first_block, first_pattern) in enumerate(patterns):
        first = dict(zip(first_block, first_pattern))
        for j, (second_block, second_pattern) in enumerate(patterns[:i]):
            second = dict(zip(second_block, second_pattern))
            if any(first[k] != second[k] for k in first.keys() & second.keys()):
                disagreements.add((j, i))
    assert disagreements == {(0, 2), (1, 2)}
    assert sum(max(5 + 1 - len(block), 0) for block, _ in patterns) == 3
    counts["forest_maximum_words"] = len(maximum)
    counts["forest_residual_words"] = len(residual)
    counts["forest_shattered_supports"] = len(expected)
    counts["forest_disagreement_edges"] = len(disagreements)


def is_downset(concepts):
    return all(all(subset in concepts for subset in range(word + 1) if subset & ~word == 0) for word in concepts)


def replay_downset_filters(counts):
    downsets = []
    for mask in range(1, 1 << 8):
        concepts = {word for word in range(8) if mask >> word & 1}
        if is_downset(concepts):
            downsets.append(concepts)
    for concepts in downsets:
        for condition in downsets:
            filtered = concepts & condition
            expected = {
                support
                for support in shattered(concepts, 3)
                if support in condition
            }
            assert filtered and shattered(filtered, 3) == expected and is_ample(filtered, 3)
            counts["downset_filter_pairs"] += 1
    counts["downsets_q3"] = len(downsets)


def main():
    assert sha256(RECONCILIATION.read_bytes()).hexdigest() == EXPECTED_RECONCILIATION
    reconciliation = json.loads(RECONCILIATION.read_text())
    required = {
        "thm:avoidance",
        "cor:cube-addition",
        "thm:minimum-section-deletion",
        "thm:one-missing-block",
        "thm:compatible-minimum-deletion",
        "thm:forest-minimum-deletion",
        "thm:monotone-filter",
    }
    entries = {entry["label"]: entry for entry in reconciliation["statements"]}
    assert required <= entries.keys()
    assert all(entries[label]["reconciliation"] == "scoped_audit_located" for label in required)
    counts = Counter()
    replay_hamming_ball(counts)
    replay_forest_example(counts)
    replay_downset_filters(counts)
    output = {
        "date": "2026-10-01",
        "status": "passed",
        "arithmetic": "exact finite sets",
        "scope": "Independent checks of both promoted deletion examples, the exact forest shattering family and disagreement graph, and every pair of downset filters on three coordinates. Universal parent constructions use database-owned proofs.",
        "source_sha256": {RECONCILIATION.name: EXPECTED_RECONCILIATION},
        "source_statement_labels": sorted(required),
        "counts": dict(counts),
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"status": "passed", "counts": dict(counts)}))


if __name__ == "__main__":
    main()
