#!/usr/bin/env python3
"""Replay the exact bounded evidence for the two-control graph construction.

The retained source packet contains the published finite audit summary.  This
database-owned verifier independently reconstructs every tested class, missing
pattern, strict rational feasibility decision, elementary minor, and bouquet
comparison.  The general all-length theorems remain justified by the written
proof recorded in the database.
"""

from __future__ import annotations

from collections import Counter, deque
from fractions import Fraction as F
from functools import cache
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "data/examples/assets"
SOURCE = ASSETS / "two_control_graphs_source.json"
OUTPUT = ASSETS / "two_control_graphs_verification.json"
SOURCE_SHA256 = "46ba60957e32630d9e810d0c935054da8ebcaba7e1399c9f2d88a6655236c389"


def sign(word: int, coordinate: int) -> int:
    return 1 if word >> coordinate & 1 else -1


@cache
def missing(words: tuple[int, ...], n: int) -> tuple[tuple[int, int], ...]:
    patterns = []
    for support in range(1, 1 << n):
        traces = {word & support for word in words}
        bits = support
        while True:
            if bits not in traces and all(
                any(
                    word & (support ^ (1 << coordinate))
                    == bits & (support ^ (1 << coordinate))
                    for word in words
                )
                for coordinate in range(n)
                if support >> coordinate & 1
            ):
                patterns.append((support, bits))
            if bits == 0:
                break
            bits = (bits - 1) & support
    return tuple(sorted(patterns))


@cache
def ample(words: tuple[int, ...], n: int) -> bool:
    shattered = sum(
        len({word & support for word in words}) == 1 << support.bit_count()
        for support in range(1 << n)
    )
    return shattered == len(words)


def vc_dimension(words: tuple[int, ...], n: int) -> int:
    return max(
        support.bit_count()
        for support in range(1 << n)
        if len({word & support for word in words}) == 1 << support.bit_count()
    )


def isometric(words: tuple[int, ...], n: int) -> bool:
    family = set(words)
    for start in words:
        distance = {start: 0}
        queue = deque([start])
        while queue:
            word = queue.popleft()
            for coordinate in range(n):
                neighbor = word ^ (1 << coordinate)
                if neighbor in family and neighbor not in distance:
                    distance[neighbor] = distance[word] + 1
                    queue.append(neighbor)
        if any(distance.get(word) != (word ^ start).bit_count() for word in words):
            return False
    return True


def family(m: int, arcs: list[tuple[int, int]], left=None):
    left = set(range(0, m, 2)) if left is None else set(left)
    assert all((source in left) != (target in left) for source, target in arcs)
    edges = {frozenset(edge) for edge in arcs}
    assert len(edges) == len(arcs) and all(len(edge) == 2 for edge in edges)
    words = set(range(4))
    order = [0, 1, 2, 3]
    expected = []
    for vertex in range(m):
        excluded = 0 if vertex in left else 3
        private = 1 << (vertex + 2)
        expected.append((private | 3, private | excluded))
        layer = [
            private | control
            for control in sorted(
                (control for control in range(4) if control != excluded),
                key=lambda control: ((control ^ (3 ^ excluded)).bit_count(), control),
            )
        ]
        words.update(layer)
        order.extend(layer)
    for source, target in arcs:
        support = (1 << (source + 2)) | (1 << (target + 2))
        control = 2 if source in left else 1
        words.add(support | control)
        order.append(support | control)
        for coordinate in range(2):
            expected.append(
                (
                    support | (1 << coordinate),
                    support | ((control ^ (1 << coordinate)) & (1 << coordinate)),
                )
            )
    for first, second in combinations(range(m), 2):
        if frozenset((first, second)) not in edges:
            support = (1 << (first + 2)) | (1 << (second + 2))
            expected.append((support, support))
    return tuple(sorted(words)), tuple(sorted(expected)), order


def ranking(m: int, arcs: list[tuple[int, int]]):
    incoming = [0] * m
    adjacency = [[] for _ in range(m)]
    for source, target in arcs:
        incoming[target] += 1
        adjacency[source].append(target)
    queue = deque(vertex for vertex in range(m) if incoming[vertex] == 0)
    order = []
    while queue:
        vertex = queue.popleft()
        order.append(vertex)
        for target in adjacency[vertex]:
            incoming[target] -= 1
            if incoming[target] == 0:
                queue.append(target)
    return {vertex: index for index, vertex in enumerate(order)} if len(order) == m else None


def rows(m: int, patterns, ratios, unequal=False):
    result = {}
    for index, pattern in enumerate(patterns):
        raw = {
            coordinate: F(1 + (coordinate + index) % 4 if unequal else 1)
            for coordinate in range(m + 2)
            if pattern[0] >> coordinate & 1
        }
        if pattern[0] & 3 == 3 and (pattern[0] >> 2).bit_count() == 1:
            vertex = (pattern[0] >> 2).bit_length() - 1
            raw[0] = F(ratios[vertex])
            raw[1] = F(1)
        total = sum(raw.values())
        result[pattern] = {coordinate: value / total for coordinate, value in raw.items()}
    return result


def valid(coefficients, concept: int, vector) -> bool:
    return all(sign(concept, i) * value > 0 for i, value in enumerate(vector)) and all(
        sum(sign(pattern[1], i) * value * vector[i] for i, value in row.items()) < 0
        for pattern, row in coefficients.items()
    )


def dot(first, second):
    return sum(x * y for x, y in zip(first, second))


def strict_linear(rows_, n: int):
    """Solve A x > 0 exactly, returning a rational point or Farkas multiplier."""
    size = len(rows_)
    initial = [
        (
            tuple(map(F, row)),
            tuple(F(int(first == second)) for second in range(size)),
        )
        for first, row in enumerate(rows_)
    ]

    def solve(system, dimension):
        compact = {}
        for row, weights in system:
            if not any(row):
                return None, weights
            scale = next(abs(value) for value in row if value)
            normalized = tuple(value / scale for value in row)
            compact.setdefault(normalized, tuple(value / scale for value in weights))
        system = list(compact.items())
        if not system:
            return (F(0),) * dimension, None
        assert dimension > 0
        positive = [(row, weight) for row, weight in system if row[-1] > 0]
        negative = [(row, weight) for row, weight in system if row[-1] < 0]
        reduced = [(row[:-1], weight) for row, weight in system if row[-1] == 0]
        for first, first_weight in positive:
            for second, second_weight in negative:
                reduced.append(
                    (
                        tuple(
                            x / first[-1] - y / second[-1]
                            for x, y in zip(first[:-1], second[:-1])
                        ),
                        tuple(
                            x / first[-1] - y / second[-1]
                            for x, y in zip(first_weight, second_weight)
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
    else:
        assert all(value >= 0 for value in dual) and any(dual)
        assert all(
            sum(value * row[i] for value, row in zip(dual, rows_)) == 0
            for i in range(n)
        )
    return point, dual


def orthant_system(n: int, coefficients, concept: int):
    system = [tuple(F(int(i == j)) for i in range(n)) for j in range(n)]
    system.extend(
        tuple(
            -sign(pattern[1], i) * sign(concept, i) * row.get(i, 0)
            for i in range(n)
        )
        for pattern, row in sorted(coefficients.items())
    )
    return tuple(system)


def witnesses(m: int, arcs, words, coefficients, left=None):
    left = set(range(0, m, 2)) if left is None else set(left)
    vertex_patterns = {
        vertex: (
            (1 << (vertex + 2)) | 3,
            (1 << (vertex + 2)) | (0 if vertex in left else 3),
        )
        for vertex in range(m)
    }
    ratios = {
        vertex: coefficients[vertex_patterns[vertex]][0]
        / coefficients[vertex_patterns[vertex]][1]
        for vertex in range(m)
    }
    assert all(ratios[source] < ratios[target] for source, target in arcs)
    result = {}
    for concept in words:
        active = [vertex for vertex in range(m) if concept >> (vertex + 2) & 1]
        vector = [F(sign(concept, i)) for i in range(m + 2)]
        if len(active) == 1:
            pattern = vertex_patterns[active[0]]
            agree = [i for i in range(2) if sign(pattern[1], i) == sign(concept, i)]
            disagree = [i for i in range(2) if i not in agree]
            large = (1 + sum(coefficients[pattern][i] for i in agree)) / sum(
                coefficients[pattern][i] for i in disagree
            )
            for i in disagree:
                vector[i] *= large
        elif len(active) == 2:
            source, target = next(
                (source, target)
                for source, target in arcs
                if {source, target} == set(active)
            )
            midpoint = (ratios[source] + ratios[target]) / 2
            vector[:2] = [F(-1), midpoint] if source in left else [F(1), -midpoint]
        relevant = [
            pattern
            for pattern in coefficients
            if all(i < 2 or i - 2 in active for i in coefficients[pattern])
        ]
        epsilon = F(1)
        for pattern in relevant:
            base = sum(
                sign(pattern[1], i) * coefficients[pattern][i] * vector[i]
                for i in coefficients[pattern]
                if i < 2
            )
            denominator = sum(
                coefficients[pattern][i] for i in coefficients[pattern] if i >= 2
            )
            assert base < 0 and denominator > 0
            epsilon = min(epsilon, -base / (2 * denominator))
        for vertex in active:
            vector[vertex + 2] = epsilon
        large = F(1)
        for pattern, row in coefficients.items():
            negative = sum(
                row[i] for i in row if i >= 2 and i - 2 not in active
            )
            if negative:
                rest = sum(
                    sign(pattern[1], i) * row[i] * vector[i]
                    for i in row
                    if i < 2 or i - 2 in active
                )
                large = max(large, (max(F(0), rest) + 1) / negative)
        for vertex in range(m):
            if vertex not in active:
                vector[vertex + 2] = -large
        assert valid(coefficients, concept, vector)
        result[concept] = tuple(vector)
    return result


def check_lp(words, n: int, coefficients, counts: Counter) -> bool:
    feasible = True
    for concept in words:
        system = orthant_system(n, coefficients, concept)
        point, dual = strict_linear(system, n)
        counts["orthant_decisions"] += 1
        if point is not None:
            signed = tuple(sign(concept, i) * point[i] for i in range(n))
            assert valid(coefficients, concept, signed)
            counts["primal_checks"] += 1
        else:
            assert all(value >= 0 for value in dual) and any(dual)
            assert all(
                sum(value * row[i] for value, row in zip(dual, system)) == 0
                for i in range(n)
            )
            counts["dual_checks"] += 1
            feasible = False
    return feasible


def learning(words, n: int) -> bool:
    family_ = set(words)
    return (
        0 in family_
        and all(first | second in family_ for first in words for second in words)
        and all(
            word == 0
            or any(
                word ^ (1 << coordinate) in family_
                for coordinate in range(n)
                if word >> coordinate & 1
            )
            for word in words
        )
    )


def ratios_to_rows(words, n: int):
    result = {}
    for pattern in missing(words, n):
        assert pattern[1].bit_count() == 1
        root = pattern[1].bit_length() - 1
        others = [
            coordinate
            for coordinate in range(n)
            if pattern[0] >> coordinate & 1 and coordinate != root
        ]
        ratios = {coordinate: F(1, n + 1) for coordinate in others}
        denominator = 1 + sum(ratios.values())
        result[pattern] = {
            root: F(1, denominator),
            **{coordinate: value / denominator for coordinate, value in ratios.items()},
        }
    return result


def feasible_order(words, concept: int, n: int):
    family_ = set(words)
    current = concept
    reverse = []
    while current:
        coordinate = next(
            coordinate
            for coordinate in range(n)
            if current >> coordinate & 1 and current ^ (1 << coordinate) in family_
        )
        reverse.append(coordinate)
        current ^= 1 << coordinate
    return tuple(reversed(reverse))


def construct_learning_witnesses(words, n: int, coefficients):
    normalized = []
    root_vector = tuple(F(1) for _ in range(n))
    for pattern, values in sorted(coefficients.items()):
        root = pattern[1].bit_length() - 1
        ratios = {
            coordinate: value / values[root]
            for coordinate, value in values.items()
            if coordinate != root
        }
        assert root_vector[root] > sum(
            ratios[coordinate] * root_vector[coordinate] for coordinate in ratios
        )
        normalized.append((root, ratios))
    result = {}
    for concept in words:
        present = {coordinate for coordinate in range(n) if concept >> coordinate & 1}
        positive = {}
        for coordinate in feasible_order(words, concept, n):
            bounds = [F(1)]
            for root, ratios in normalized:
                if root == coordinate:
                    previous = sum(
                        ratios[j] * positive[j] for j in ratios if j in positive
                    )
                    assert previous > 0
                    bounds.append(previous / 2)
            positive[coordinate] = min(bounds)
        bounds = [F(1)]
        for root, ratios in normalized:
            if root in present:
                gap = (
                    sum(ratios[j] * positive[j] for j in ratios if j in present)
                    - positive[root]
                )
                tail = sum(
                    ratios[j] * root_vector[j] for j in ratios if j not in present
                )
                assert gap > 0
                bounds.append(gap / (2 * (1 + tail)))
        scale = min(bounds)
        vector = tuple(
            positive[i] if i in present else -scale * root_vector[i]
            for i in range(n)
        )
        assert valid(coefficients, concept, vector)
        result[concept] = vector
    return result


def deleted(word: int, coordinate: int) -> int:
    return (word & ((1 << coordinate) - 1)) | ((word >> (coordinate + 1)) << coordinate)


def minor(words, coordinate: int, condition=None):
    return tuple(
        sorted(
            {
                deleted(word, coordinate)
                for word in words
                if condition is None or word >> coordinate & 1 == condition
            }
        )
    )


def certify_minor(words, n: int, counts: Counter):
    root = next(
        (
            root
            for root in words
            if learning(tuple(sorted(word ^ root for word in words)), n)
        ),
        None,
    )
    if root is not None:
        reflected = tuple(sorted(word ^ root for word in words))
        coefficients = ratios_to_rows(reflected, n)
        result = construct_learning_witnesses(reflected, n, coefficients)
        assert len(result) == len(reflected)
        counts["minor_rational_witnesses"] += len(result)
        return "learning-space reflection"
    coefficients = {
        pattern: {
            coordinate: F(1, pattern[0].bit_count())
            for coordinate in range(n)
            if pattern[0] >> coordinate & 1
        }
        for pattern in missing(words, n)
    }
    assert check_lp(words, n, coefficients, counts)
    counts["minor_rational_witnesses"] += len(words)
    return "uniform-row exact LP"


def replay():
    assert sha256(SOURCE.read_bytes()).hexdigest() == SOURCE_SHA256
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    assert source["status"] == "passed"

    counts = Counter()
    graph_records = []
    inputs = [(0, []), (1, []), (3, [(0, 1), (1, 2)]), (3, [(0, 1), (2, 1)])]
    for orientation in range(16):
        cycle_edges = [(i, (i + 1) % 4) for i in range(4)]
        inputs.append(
            (
                4,
                [
                    edge if orientation >> i & 1 else edge[::-1]
                    for i, edge in enumerate(cycle_edges)
                ],
            )
        )
    for m, arcs in inputs:
        words, expected, order = family(m, arcs)
        n = m + 2
        assert missing(words, n) == expected
        assert ample(words, n) and isometric(words, n)
        assert vc_dimension(words, n) == 2
        assert all(ample(tuple(sorted(order[:size])), n) for size in range(1, len(order) + 1))
        reflected = {word ^ 1 for word in words}
        assert all(first & second in reflected for first in reflected for second in reflected)
        rank = ranking(m, arcs)
        choices = [
            {vertex: F(1) for vertex in range(m)},
            {vertex: F(2**vertex) for vertex in range(m)},
        ]
        if rank is not None:
            choices.append({vertex: F(2 ** rank[vertex]) for vertex in range(m)})
        for index, ratios in enumerate(choices):
            coefficients = rows(m, expected, ratios, unequal=bool(index % 2))
            expected_feasible = all(ratios[source_] < ratios[target] for source_, target in arcs)
            assert check_lp(words, n, coefficients, counts) == expected_feasible
            if expected_feasible:
                counts["constructed_full_witnesses"] += len(
                    witnesses(m, arcs, words, coefficients)
                )
            counts["row_families"] += 1
        graph_records.append(
            {
                "vertices": m,
                "arcs": [list(arc) for arc in arcs],
                "concepts": len(words),
                "acyclic": rank is not None,
            }
        )

    minor_records = []
    for m in (4, 6):
        arcs = [(vertex, (vertex + 1) % m) for vertex in range(m)]
        words, expected, order = family(m, arcs)
        n = m + 2
        assert len(words) == 4 * m + 4 and missing(words, n) == expected
        assert ample(words, n) and isometric(words, n) and vc_dimension(words, n) == 2
        reflected = {word ^ 1 for word in words}
        assert all(first & second in reflected for first in reflected for second in reflected)
        assert all(ample(tuple(sorted(order[:size])), n) for size in range(1, len(order) + 1))
        for coordinate in range(n):
            for condition in (None, 0, 1):
                reduced = minor(words, coordinate, condition)
                assert reduced and ample(reduced, n - 1) and isometric(reduced, n - 1)
                if coordinate >= 2 and condition in (None, 0):
                    removed = coordinate - 2
                    kept = [vertex for vertex in range(m) if vertex != removed]
                    index = {vertex: position for position, vertex in enumerate(kept)}
                    reduced_arcs = [
                        (index[source_], index[target])
                        for source_, target in arcs
                        if removed not in (source_, target)
                    ]
                    left = {index[vertex] for vertex in kept if vertex % 2 == 0}
                    built, patterns, _ = family(m - 1, reduced_arcs, left)
                    assert built == reduced
                    rank = ranking(m - 1, reduced_arcs)
                    assert rank is not None
                    coefficients = rows(
                        m - 1,
                        patterns,
                        {vertex: F(2 ** rank[vertex]) for vertex in rank},
                        unequal=True,
                    )
                    counts["minor_rational_witnesses"] += len(
                        witnesses(m - 1, reduced_arcs, reduced, coefficients, left)
                    )
                    assert check_lp(reduced, n - 1, coefficients, counts)
                    method = "acyclic graph: ranked ratios and exact LP"
                else:
                    method = certify_minor(reduced, n - 1, counts)
                minor_records.append(
                    {
                        "cycle": m,
                        "coordinate": coordinate,
                        "condition": condition,
                        "concepts": len(reduced),
                        "certificate": method,
                    }
                )
                counts["elementary_minors"] += 1

    words, _, _ = family(4, [(vertex, (vertex + 1) % 4) for vertex in range(4)])
    mapping = {0: 0, 2: 1, 3: 2, 1: 3}
    transformed = {
        sum(((word >> (vertex + 2)) & 1) << mapping[vertex] for vertex in range(4))
        | (((word & 3) ^ 1) << 4)
        for word in words
    }
    bouquet = set()
    for support in range(16):
        controls = (
            [0, 1, 2, 3]
            if support == 0
            else [0, 2, 3]
            if support & 12 == 0
            else [0, 1, 3]
            if support & 3 == 0
            else [0]
            if support in (5, 10)
            else [3]
            if support in (6, 9)
            else []
        )
        bouquet.update(support | (control << 4) for control in controls)
    deleted_words = {
        support | (control << 4)
        for support, controls in [(3, [0, 2, 3]), (12, [0, 1, 3])]
        for control in controls
    }
    assert len(bouquet) == 26 and transformed < bouquet
    assert bouquet - transformed == deleted_words
    counts["source_bouquet_deleted_concepts"] = 6

    assert dict(counts) == source["counts"]
    assert graph_records == source["graphs"]
    assert minor_records == source["minors"]
    return {
        "status": "passed",
        "source_sha256": SOURCE_SHA256,
        "counts": dict(counts),
        "graphs_replayed": len(graph_records),
        "cyclic_square_orientations": sum(
            not record["acyclic"] for record in graph_records if record["vertices"] == 4
        ),
        "elementary_minors_replayed": len(minor_records),
        "cycle_lengths_replayed": [4, 6],
        "bouquet_comparison": {
            "bouquet_concepts": len(bouquet),
            "d4_concepts": len(transformed),
            "deleted_concepts": len(deleted_words),
        },
        "scope": source["limits"],
    }


def main():
    report = replay()
    OUTPUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
