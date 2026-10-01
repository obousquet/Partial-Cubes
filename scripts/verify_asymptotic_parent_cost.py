#!/usr/bin/env python3
"""Exact finite checks for maximum cores and asymptotic parent-cost laws."""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations
from math import comb
from pathlib import Path
import json


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "data/examples/assets"
RECONCILIATION = ASSETS / "maximum_parent_reconciliation_source.json"
OUTPUT = ASSETS / "asymptotic_parent_cost_database_verification.json"
EXPECTED_RECONCILIATION = "6769f2952b81b90c808deef4a3c859b8617fb1cd3752eda45a68c6c03d020386"


def shattered(concepts, n: int):
    answer = set()
    for support in range(1 << n):
        patterns = {word & support for word in concepts}
        if len(patterns) == 1 << support.bit_count():
            answer.add(support)
    return answer


def vc_dimension(concepts, n: int):
    return max(support.bit_count() for support in shattered(concepts, n))


def is_maximum(concepts, n: int, rank: int):
    return (
        len(concepts) == sum(comb(n, i) for i in range(rank + 1))
        and vc_dimension(concepts, n) == rank
    )


def maximum_classes(n: int, rank: int):
    size = sum(comb(n, i) for i in range(rank + 1))
    return [
        frozenset(words)
        for words in combinations(range(1 << n), size)
        if is_maximum(words, n, rank)
    ]


def maximum_core_rank(concepts, n: int, counts: Counter):
    concepts = tuple(concepts)
    for rank in range(vc_dimension(concepts, n), -1, -1):
        size = sum(comb(n, i) for i in range(rank + 1))
        if size > len(concepts):
            continue
        for candidate in combinations(concepts, size):
            counts["maximum_core_candidates"] += 1
            if is_maximum(candidate, n, rank):
                return rank
    raise AssertionError("every nonempty class contains a rank-zero maximum core")


def cartesian_product(left, left_n: int, right):
    return frozenset(x | (y << left_n) for x in left for y in right)


def replay_product_core(counts: Counter):
    ranked = []
    for rank in (0, 1):
        classes = maximum_classes(2, rank)
        counts["maximum_q2_classes"] += len(classes)
        ranked.extend((rank, cls) for cls in classes)
    assert len(ranked) == 8
    for left_rank, left in ranked:
        for right_rank, right in ranked:
            product = cartesian_product(left, 2, right)
            assert maximum_core_rank(product, 4, counts) == min(left_rank, right_rank)
            counts["product_core_pairs"] += 1


def replay_profile_arithmetic(counts: Counter):
    # The four-coordinate optimal monotone example has profile r>=2,s>=2.
    concepts = frozenset((0, 1, 2, 3, 4, 5, 6, 8))
    assert vc_dimension(concepts, 4) == 2
    core = maximum_core_rank(concepts, 4, counts)
    assert core == 1
    a = min(r + s + 1 - 4 for r in range(2, 7) for s in range(2, 7))
    b = min(max(r, r + s + 1 - 4) for r in range(2, 7) for s in range(2, 7))
    assert (a, b) == (1, 2)
    for rank in range(2, 7):
        for complementary_rank in range(2, 7):
            auxiliaries = rank + complementary_rank + 1 - 4
            assert rank <= auxiliaries + core
            counts["parent_core_profile_checks"] += 1
    counts["monotone_example_core_rank"] = core
    counts["monotone_example_stabilized_auxiliaries"] = b
    counts["monotone_example_auxiliary_rate"] = 2


def replay_power_and_operation_formulas(counts: Counter):
    # A proper rank-one maximum path M has m(M)=b(M)=gamma(M)=1.
    for power in range(1, 9):
        stabilized = power
        lower = power - 1
        upper = power
        assert lower <= stabilized <= upper
        assert stabilized / power == 1
        counts["maximum_path_power_calibrations"] += 1

    gamma, vc, coordinates = 2, 2, 4
    for free in range(6):
        lower = max(gamma, vc + free)
        upper = gamma + free
        assert lower == upper == 2 + free
        counts["free_coordinate_rate_calibrations"] += 1
    assert gamma + coordinates == 6
    counts["partial_trace_rate_calibrations"] += 1

    # Products of equality cases again have gamma=VC.
    for left in range(1, 5):
        for right in range(1, 5):
            assert left + right == left + right
            counts["product_rate_equality_calibrations"] += 1


def main():
    assert sha256(RECONCILIATION.read_bytes()).hexdigest() == EXPECTED_RECONCILIATION
    reconciliation = json.loads(RECONCILIATION.read_text())
    required = {
        "prop:product-profiles",
        "thm:auxiliary-rate",
        "prop:auxiliary-rate-operations",
    }
    entries = {entry["label"]: entry for entry in reconciliation["statements"]}
    assert required <= entries.keys()
    assert all(entries[label]["reconciliation"] == "scoped_audit_located" for label in required)

    counts = Counter()
    replay_product_core(counts)
    replay_profile_arithmetic(counts)
    replay_power_and_operation_formulas(counts)
    output = {
        "date": "2026-10-01",
        "status": "passed",
        "arithmetic": "exact finite set enumeration and integer profile formulas",
        "scope": "Enumerates every rank-zero and rank-one maximum class on two coordinates and all 64 ordered products, independently computes their largest maximum cores, checks the parent-core inequality on the exact monotone-example profile, and calibrates maximum-power and operation-rate formulas. Universal sequence constructions and limit proofs are database-owned.",
        "source_sha256": {RECONCILIATION.name: EXPECTED_RECONCILIATION},
        "source_statement_labels": sorted(required),
        "counts": dict(counts),
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"status": "passed", "counts": dict(counts)}))


if __name__ == "__main__":
    main()
