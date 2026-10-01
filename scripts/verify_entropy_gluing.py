#!/usr/bin/env python3
"""Replay exact entropy-balance and common-section gluing certificates."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path
import json


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "data/examples/assets"
BALANCE = ASSETS / "entropy_balance_source.json"
GLUING = ASSETS / "entropy_gluing_source.json"
MARGIN = ASSETS / "margin_incompatible_extensions_verification.json"
OUTPUT = ASSETS / "entropy_gluing_database_verification.json"

EXPECTED_SHA256 = {
    BALANCE.name: "3b42a3dc917d51f343ba596fda7ef48cca44f5b1f461c597a50f8d9af82a6f38",
    GLUING.name: "b3d6663b5119ca62f2324a75311fda3d39130fda9bad8fa17c8b01e653a4b1d0",
}


def missing(concepts, n: int):
    concepts = tuple(concepts)
    result = []
    for support in range(1, 1 << n):
        traces = {concept & support for concept in concepts}
        bits = support
        while True:
            if bits not in traces and all(
                any(
                    concept & (support ^ (1 << coordinate))
                    == bits & (support ^ (1 << coordinate))
                    for concept in concepts
                )
                for coordinate in range(n)
                if support >> coordinate & 1
            ):
                result.append((support, bits))
            if bits == 0:
                break
            bits = (bits - 1) & support
    return set(result)


def dot(first, second):
    return sum(x * y for x, y in zip(first, second))


def strict_linear(rows_, n: int):
    """Solve A x > 0 by exact Fourier--Motzkin elimination."""
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
        low = max(lower, default=None)
        high = min(upper, default=None)
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


def ratio(coefficients, magnitudes, concept: int):
    result = F(0)
    for (support, bits), row in coefficients.items():
        favorable = sum(
            value * magnitudes[i]
            for i, value in row.items()
            if (concept ^ bits) >> i & 1
        )
        unfavorable = sum(
            value * magnitudes[i]
            for i, value in row.items()
            if not (concept ^ bits) >> i & 1
        )
        assert favorable > 0
        result = max(result, unfavorable / favorable)
    return result


def solve_threshold(concepts, n: int, coefficients, threshold):
    witnesses = {}
    for concept in concepts:
        rows = [tuple(F(i == j) for j in range(n)) for i in range(n)]
        for (support, bits), entry in coefficients.items():
            rows.append(
                tuple(
                    entry.get(i, F(0))
                    * (threshold if (concept ^ bits) >> i & 1 else -1)
                    for i in range(n)
                )
            )
        vector, _ = strict_linear(rows, n)
        if vector is None:
            return None
        assert ratio(coefficients, vector, concept) < threshold
        scale = max(vector, default=F(1))
        witnesses[concept] = tuple(value / scale for value in vector)
    return witnesses


def lift(word: int, block: int):
    return (word & 3) | (((word >> 2) & 1) << block)


def replay_gluing_thresholds(source, counts):
    records = []
    cases = [
        ((0, 1, 2), (1, 2), (0, 1, 2)),
        ((0, 1, 2), (0, 1), (0, 2)),
        ((0, 1, 2), (), (0, 1, 2)),
        ((0, 3), (0, 3), (0, 3)),
        ((0, 3), (0,), (3,)),
        ((0, 3), (), (0, 3)),
        ((0, 1, 2, 3), (0, 3), (0, 1, 2, 3)),
        ((0, 1, 2, 3), (0, 1, 2), (0, 1, 3)),
    ]
    for section, top1, top2 in cases:
        first = tuple(sorted(set(section) | {c + 4 for c in top1}))
        second = tuple(sorted(set(section) | {c + 4 for c in top2}))
        union = tuple(sorted({lift(c, 2) for c in first} | {lift(c, 3) for c in second}))
        rows1 = set(missing(first, 3))
        rows2 = {tuple(lift(value, 3) for value in row) for row in missing(second, 3)}
        cross = {(12, 12)} if top1 and top2 else set()
        assert rows1 | rows2 | cross == missing(union, 4)
        counts["row_partitions"] += 1
        for variant in (0, 1):
            coefficients = {
                row: {
                    i: F(1 + (row[0] + 2 * row[1] + i + variant) % 5)
                    for i in range(4)
                    if row[0] >> i & 1
                }
                for row in missing(union, 4)
            }
            coefficients = {
                row: {i: value / sum(entry.values()) for i, value in entry.items()}
                for row, entry in coefficients.items()
            }
            if cross:
                coefficients[(12, 12)] = {
                    2: F(1 if variant == 0 else 100, 101),
                    3: F(100 if variant == 0 else 1, 101),
                }
            coefficients1 = {row: coefficients[row] for row in missing(first, 3)}
            coefficients2 = {
                row: {
                    (2 if i == 3 else i): value
                    for i, value in coefficients[tuple(lift(value, 3) for value in row)].items()
                }
                for row in missing(second, 3)
            }
            for threshold in map(F, ("1/4", "3/4", "1", "5/4", "2", "4")):
                witnesses1 = solve_threshold(first, 3, coefficients1, threshold)
                witnesses2 = solve_threshold(second, 3, coefficients2, threshold)
                union_witnesses = solve_threshold(union, 4, coefficients, threshold)
                feasible = witnesses1 is not None and witnesses2 is not None
                assert (union_witnesses is not None) == feasible
                counts["threshold_cases"] += 1
                if feasible:
                    private_masses = [
                        sum(value for i, value in row.items() if i >= 2)
                        for pattern, row in coefficients.items()
                        if pattern[0] & 12
                    ]
                    beta = min(
                        [F(1)]
                        + private_masses
                        + [value for pattern in cross for value in coefficients[pattern].values()]
                    )
                    epsilon = min(F(1, 2), threshold * beta / 2)
                    for concept in union:
                        if concept & 4:
                            vector = tuple(epsilon * value for value in witnesses1[concept]) + (F(1),)
                        elif concept & 8:
                            local = (concept & 3) | 4
                            witness = witnesses2[local]
                            vector = (
                                epsilon * witness[0],
                                epsilon * witness[1],
                                F(1),
                                epsilon * witness[2],
                            )
                        else:
                            vector = (
                                epsilon * witnesses1[concept][0],
                                epsilon * witnesses1[concept][1],
                                F(1),
                                F(1),
                            )
                        assert ratio(coefficients, vector, concept) < threshold
                        counts["constructed_witnesses"] += 1
                records.append(
                    {
                        "D": list(section),
                        "top1": list(top1),
                        "top2": list(top2),
                        "variant": variant,
                        "t": str(threshold),
                        "feasible": feasible,
                    }
                )
    assert records == source["threshold_checks"]


def decode_certificate(packet):
    concepts = tuple(packet["words"])
    rows = {
        (row["support"], row["bits"]): {
            int(i): F(value) for i, value in row["coefficients"].items()
        }
        for row in packet["rows"]
    }
    witnesses = {
        int(concept): tuple(map(F, vector))
        for concept, vector in packet["witnesses"].items()
    }
    assert set(rows) == missing(concepts, packet["n"])
    return concepts, rows, witnesses


def dominance(packet):
    concepts, rows, witnesses = decode_certificate(packet)
    comparisons = 0
    for concept, pattern in product(concepts, rows):
        favorable = [
            coefficient * abs(witnesses[concept][i])
            for i, coefficient in rows[pattern].items()
            if (concept ^ pattern[1]) >> i & 1
        ]
        unfavorable = [
            coefficient * abs(witnesses[concept][i])
            for i, coefficient in rows[pattern].items()
            if not (concept ^ pattern[1]) >> i & 1
        ]
        assert max(favorable) > max(unfavorable, default=F(0))
        comparisons += 1
    return comparisons


def cycle_rank(patterns, coordinates):
    vertices = [("p", pattern) for pattern in patterns] + [("i", i) for i in coordinates]
    parent = {vertex: vertex for vertex in vertices}

    def root(vertex):
        while parent[vertex] != vertex:
            vertex = parent[vertex]
        return vertex

    edges = 0
    for pattern in patterns:
        for coordinate in coordinates:
            if pattern[0] >> coordinate & 1:
                edges += 1
                first, second = root(("p", pattern)), root(("i", coordinate))
                parent[first] = second
    components = len({root(vertex) for vertex in vertices})
    return edges - len(vertices) + components


def drop(word: int, coordinate: int):
    return (word & ((1 << coordinate) - 1)) | ((word >> (coordinate + 1)) << coordinate)


def deterministic_array(concepts, n: int):
    return {
        pattern: {
            i: F(1 + (pattern[0] + pattern[1] + i) % 4)
            for i in range(n)
            if pattern[0] >> i & 1
        }
        for pattern in missing(concepts, n)
    }


def minor_rows(concepts, n: int, coefficients, coordinate: int, fixed, last: bool):
    kept = [
        concept
        for concept in concepts
        if fixed is None or concept >> coordinate & 1 == fixed
    ]
    if not kept:
        return None
    image = tuple(sorted({drop(concept, coordinate) for concept in kept}))
    result = {}
    for target in missing(image, n - 1):
        choices = [
            pattern
            for pattern in sorted(coefficients)
            if (drop(pattern[0], coordinate), drop(pattern[1], coordinate)) == target
            and (
                not pattern[0] >> coordinate & 1
                or (fixed is not None and pattern[1] >> coordinate & 1 == fixed)
            )
        ]
        assert choices
        pattern = choices[-1 if last else 0]
        row = {
            i if i < coordinate else i - 1: value
            for i, value in coefficients[pattern].items()
            if i != coordinate
        }
        total = sum(row.values())
        result[target] = {i: value / total for i, value in row.items()}
    return image, result


def replay_balance_calibrations(source, counts):
    for record in source["minor_maps"]:
        concepts = tuple(record["C"])
        coefficients = deterministic_array(concepts, 2)
        reduced = minor_rows(
            concepts,
            2,
            coefficients,
            record["coordinate"],
            record["fixed"],
            record["last"],
        )
        assert reduced is not None
        image, rows = reduced
        magnitudes = {
            concept: tuple(F(1 + (3 * concept + 2 * i) % 5) for i in range(2))
            for concept in concepts
        }
        values = []
        for target in image:
            lift = next(
                concept
                for concept in concepts
                if drop(concept, record["coordinate"]) == target
                and (
                    record["fixed"] is None
                    or concept >> record["coordinate"] & 1 == record["fixed"]
                )
            )
            vector = tuple(
                value
                for i, value in enumerate(magnitudes[lift])
                if i != record["coordinate"]
            )
            values.append(str(ratio(rows, vector, target)))
            counts["minor_concepts"] += 1
        assert values == record["ratios"]
        counts["minor_maps"] += 1

    feature_balance = Counter()
    for row in source["opposite_pair_balance"]:
        mass = F(row["mass"])
        for i in row["bad"]:
            feature_balance[("a", tuple(row["p"]), i)] += mass
            feature_balance[("z", row["c"], i)] += mass
        for i in row["good"]:
            feature_balance[("a", tuple(row["p"]), i)] -= mass
            feature_balance[("z", row["c"], i)] -= mass
    assert all(value == 0 for value in feature_balance.values())
    counts["opposite_pair_balance_rows"] = len(source["opposite_pair_balance"])

    for record in source["branching_balances"]:
        p = F(record["p"])
        masses = tuple(map(F, record["lambda_values"]))
        assert masses == (F(1, 2), p / 2, (1 - p) / 2)
        differences = ((-p, -(1 - p), F(1)), (F(1), F(0), F(-1)), (F(0), F(1), F(-1)))
        assert sum(masses) == 1
        assert all(sum(masses[r] * differences[r][i] for r in range(3)) == 0 for i in range(3))
        assert record["entropy_value"] == "-H(p)/2"
        counts["branching_balances"] += 1


def main():
    actual_hashes = {
        path.name: sha256(path.read_bytes()).hexdigest() for path in (BALANCE, GLUING)
    }
    assert actual_hashes == EXPECTED_SHA256
    balance_source = json.loads(BALANCE.read_text())
    gluing_source = json.loads(GLUING.read_text())
    margin = json.loads(MARGIN.read_text())
    assert balance_source["status"] == gluing_source["status"] == margin["status"] == "passed"

    counts = Counter()
    replay_balance_calibrations(balance_source, counts)
    replay_gluing_thresholds(gluing_source, counts)

    counts["dominance_comparisons"] = sum(
        dominance(margin["certificates"][name]) for name in ("A", "B")
    )
    assert counts["dominance_comparisons"] == gluing_source["counts"]["dominance_comparisons"]

    first, _, _ = decode_certificate(margin["certificates"]["A"])
    second, _, _ = decode_certificate(margin["certificates"]["B"])
    union = tuple(sorted(set(first) | {(concept & ~8) | ((concept & 8) << 3) for concept in second}))
    assert union == tuple(margin["union"]["words"]) and len(union) == 33
    all_rows = {
        pattern: {i: F(1) for i in range(7) if pattern[0] >> i & 1}
        for pattern in missing(union, 7)
    }
    assert all(not (support ^ bits) & 79 for support, bits in all_rows)
    upper = []
    for epsilon in map(F, ("1/10", "1/100", "1/1000")):
        absent = 2 + 5 * epsilon
        worst = max(
            ratio(
                all_rows,
                tuple(
                    F(1) if i in (4, 5) else epsilon if concept >> i & 1 else absent
                    for i in range(7)
                ),
                concept,
            )
            for concept in union
        )
        assert worst <= 1 + 5 * epsilon
        upper.append({"eps": str(epsilon), "max_ratio": str(worst), "bound": str(1 + 5 * epsilon)})
        counts["union_upper_witnesses"] += len(union)
    assert upper == gluing_source["union_upper_arrays"]

    boundary_a = {
        (tuple(row["pattern"]), row["coordinate"]): F(row["value"])
        for row in gluing_source["boundary_A"]
    }
    boundary_b = {
        (tuple(row["pattern"]), row["coordinate"]): F(row["value"])
        for row in gluing_source["boundary_B"]
    }
    assert set(boundary_a) == set(boundary_b)
    assert all(boundary_a[key] + boundary_b[key] == 0 for key in boundary_a)
    for pattern, _ in boundary_a:
        assert sum(value for (other, _), value in boundary_a.items() if other == pattern) == 0
    for coordinate in (0, 1, 2, 4, 5):
        assert sum(value for (_, i), value in boundary_a.items() if i == coordinate) == 0
    assert sum(map(abs, boundary_a.values())) <= 2

    shared = {pattern for pattern in missing(first, 6) if not pattern[0] & 8}
    rank = cycle_rank(shared, (0, 1, 2, 4, 5))
    assert rank == gluing_source["cycle_rank"] == 7
    assert dict(counts)["threshold_cases"] == gluing_source["counts"]["threshold_cases"]
    assert dict(counts)["constructed_witnesses"] == gluing_source["counts"]["constructed_witnesses"]
    assert dict(counts)["row_partitions"] == gluing_source["counts"]["row_partitions"]

    output = {
        "date": "2026-10-01",
        "status": "passed",
        "arithmetic": "fractions.Fraction; exact Fourier--Motzkin elimination; no floating-point optimizer",
        "scope": "Independent replay of retained bounded entropy calibrations, all common-section threshold reconstructions, strict-dominance certificates, opposite section circulation, and union upper arrays. Universal duality, composition, and forest laws use database-owned proofs.",
        "source_sha256": actual_hashes,
        "counts": dict(counts),
        "cycle_rank": rank,
        "boundary_entries_per_piece": len(boundary_a),
        "opposite_boundaries": True,
        "global_mass_per_piece": gluing_source["global_mass_per_piece"],
        "global_entropy": gluing_source["global_entropy"],
        "union_upper_arrays": upper,
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"status": "passed", "counts": dict(counts), "cycle_rank": rank}))


if __name__ == "__main__":
    main()
