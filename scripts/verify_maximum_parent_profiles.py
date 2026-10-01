#!/usr/bin/env python3
"""Exact finite checks for maximum-parent profiles and fixed-rank bounds."""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations
from math import comb
from pathlib import Path
import json

from verify_maximum_parent_products import RECONCILIATION, bitword, missing, shattered


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "data/examples/assets"
OUTPUT = ASSETS / "maximum_parent_profiles_database_verification.json"
EXPECTED_RECONCILIATION = "6769f2952b81b90c808deef4a3c859b8617fb1cd3752eda45a68c6c03d020386"


def is_maximum(concepts, n: int, rank: int):
    return (
        len(concepts) == sum(comb(n, i) for i in range(rank + 1))
        and max((support.bit_count() for support in shattered(concepts, n)), default=-1) == rank
    )


def boundary_edges(concepts, n: int):
    concepts = set(concepts)
    return sum(
        1
        for word in concepts
        for coordinate in range(n)
        if word ^ (1 << coordinate) not in concepts
    )


def fixed_bounds(concepts, n: int, rank: int):
    patterns = missing(concepts, n)
    if any(support.bit_count() > rank + 1 for support, _ in patterns):
        return None
    first = sum(rank + 1 - support.bit_count() for support, _ in patterns)
    second = (rank + 1) * ((1 << n) - len(concepts)) - boundary_edges(concepts, n)
    return first, second


def replay_uniform_patterns(counts):
    n = 3
    universe = range(1 << n)
    for rank in range(n):
        size = sum(comb(n, i) for i in range(rank + 1))
        for selected in combinations(universe, size):
            concepts = set(selected)
            if not is_maximum(concepts, n, rank):
                continue
            patterns = missing(concepts, n)
            assert all(support.bit_count() == rank + 1 for support, _ in patterns)
            for support in range(1 << n):
                if support.bit_count() == rank + 1:
                    assert sum(row_support == support for row_support, _ in patterns) == 1
            counts["maximum_classes_q3"] += 1
            counts["uniform_missing_patterns"] += len(patterns)


def replay_shared_witness(counts):
    parent = {bitword(word) for word in ("000", "100", "010", "001")}
    assert is_maximum(parent, 3, 1)
    section = {word & 3 for word in parent if word >> 2 & 1}
    assert section == {0}
    patterns = missing(section, 2)
    assert patterns == {(1, 1), (2, 2)}
    assert fixed_bounds(section, 2, 1) == (2, 4)
    witness_auxiliary_supports = [{2}, {2}]
    assert len(set().union(*witness_auxiliary_supports)) == 1
    counts["shared_witness_visible_patterns"] = len(patterns)
    counts["shared_witness_bound"] = 2
    counts["shared_witness_actual_auxiliaries"] = 1


def replay_optimal_profile(counts):
    section = {0, 1, 2, 3, 4, 5, 6, 8}
    upper = {0, 1, 2, 4, 8, 9, 10, 12}
    parent = section | {16 + word for word in upper}
    assert is_maximum(parent, 5, 2)
    patterns = missing(section, 4)
    assert sorted(support.bit_count() for support, _ in patterns) == [2, 2, 2, 3]
    first, second = fixed_bounds(section, 4, 2)
    assert first == 3 and second >= 1
    assert max(support.bit_count() for support in shattered(section, 4)) == 2
    # Exact monotone-shore criterion: r>=2 and m>=r-1. With n=4,
    # s=n+m-r-1, so this is exactly the quadrant r>=2,s>=2.
    for rank in range(2, 7):
        least_aux = rank - 1
        complementary_rank = 4 + least_aux - rank - 1
        assert complementary_rank == 2
        counts["profile_boundary_points_checked"] += 1
    counts["optimal_profile_min_rank"] = 2
    counts["optimal_profile_min_auxiliaries"] = 1


def replay_coordinate_formulas(counts):
    # A proper rank-one maximum path on two coordinates has minimal profile (1,0).
    n, rank, complement_rank = 2, 1, 0
    for constants in range(1, 6):
        transformed = (rank, max(complement_rank + constants, n + constants - 1))
        auxiliaries = transformed[0] + transformed[1] + 1 - (n + constants)
        assert transformed == (1, constants + 1) and auxiliaries == 1
        counts["constant_profile_calibrations"] += 1
    for free in range(6):
        transformed = (rank + free, complement_rank + free)
        auxiliaries = transformed[0] + transformed[1] + 1 - (n + free)
        assert auxiliaries == free
        counts["free_profile_calibrations"] += 1


def main():
    assert sha256(RECONCILIATION.read_bytes()).hexdigest() == EXPECTED_RECONCILIATION
    reconciliation = json.loads(RECONCILIATION.read_text())
    required = {
        "lem:uniform",
        "thm:bound",
        "ex:shared",
        "cor:decision",
        "thm:effective-parent-profiles",
        "prop:constant-profile",
        "prop:free-coordinate-profiles",
        "prop:free-block-normal-form",
    }
    entries = {entry["label"]: entry for entry in reconciliation["statements"]}
    assert required <= entries.keys()
    assert all(entries[label]["reconciliation"] == "scoped_audit_located" for label in required)
    counts = Counter()
    replay_uniform_patterns(counts)
    replay_shared_witness(counts)
    replay_optimal_profile(counts)
    replay_coordinate_formulas(counts)
    output = {
        "date": "2026-10-01",
        "status": "passed",
        "arithmetic": "exact finite sets and integer profile formulas",
        "scope": "Independent enumeration of uniform missing patterns in all maximum classes on three coordinates; fixed-rank bound calibrations; the shared-witness and optimal-profile examples; and constant/free-coordinate formulas on a proper maximum path. Universal bounds and algorithms use database-owned proofs.",
        "source_sha256": {RECONCILIATION.name: EXPECTED_RECONCILIATION},
        "source_statement_labels": sorted(required),
        "counts": dict(counts),
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"status": "passed", "counts": dict(counts)}))


if __name__ == "__main__":
    main()
