#!/usr/bin/env python3
"""Replay exact source-sign and scaling-indicator certificates."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import json


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "data/examples/assets"
SIGN_SOURCE = ASSETS / "source_sign_decomposition_source.json"
SCALING_SOURCE = ASSETS / "coefficient_scaling_source.json"
OUTPUT = ASSETS / "source_scaling_invariants_verification.json"


def sign(word: int, coordinate: int) -> int:
    return 1 if word >> coordinate & 1 else -1


def minimal_missing(words: tuple[int, ...], n: int) -> set[tuple[int, int]]:
    family = set(words)
    out: set[tuple[int, int]] = set()
    for support in range(1, 1 << n):
        bits = support
        while True:
            if not any((word & support) == bits for word in family) and all(
                any(
                    (word & (support ^ (1 << coordinate)))
                    == (bits & ~(1 << coordinate))
                    for word in family
                )
                for coordinate in range(n)
                if support >> coordinate & 1
            ):
                out.add((support, bits))
            if bits == 0:
                break
            bits = (bits - 1) & support
    return out


def decode_rows(encoded):
    return {
        (record["support"], record["bits"]): {
            int(coordinate): Fraction(value)
            for coordinate, value in record["coefficients"]
        }
        for record in encoded
    }


def valid(rows, concept: int, vector: tuple[Fraction, ...]) -> bool:
    if not all(sign(concept, i) * vector[i] > 0 for i in range(len(vector))):
        return False
    return all(
        sum(sign(bits, i) * coefficient * vector[i] for i, coefficient in row.items()) < 0
        for (_, bits), row in rows.items()
    )


def orthant_system(n: int, rows, concept: int):
    system = [
        tuple(Fraction(int(i == j)) for i in range(n))
        for j in range(n)
    ]
    system.extend(
        tuple(
            -sign(bits, i) * sign(concept, i) * row.get(i, 0)
            for i in range(n)
        )
        for (_, bits), row in sorted(rows.items())
    )
    return tuple(system)


def check_dual(system, encoded) -> None:
    weights = tuple(Fraction(value) for value in encoded)
    assert len(weights) == len(system)
    assert any(weights) and all(value >= 0 for value in weights)
    assert all(
        sum(weight * row[i] for weight, row in zip(weights, system)) == 0
        for i in range(len(system[0]))
    )


def graph_sources(words: tuple[int, ...], n: int) -> tuple[int, ...]:
    family = set(words)
    return tuple(
        word
        for word in words
        if not any(
            word ^ (1 << coordinate) in family
            for coordinate in range(n)
            if word >> coordinate & 1
        )
    )


def positive_core(coordinates, rows):
    coordinates = set(coordinates)
    rows = list(rows)
    while True:
        roots = {bits.bit_length() - 1 for _, bits in rows}
        free = coordinates - roots
        if not free:
            return tuple(sorted(coordinates)), tuple(rows)
        rows = [
            pattern
            for pattern in rows
            if not any(
                (pattern[0] ^ pattern[1]) >> coordinate & 1
                for coordinate in free
            )
        ]
        coordinates -= free


def acyclic(nodes, arcs) -> bool:
    successors = {node: set() for node in nodes}
    indegree = {node: 0 for node in nodes}
    for first, second in arcs:
        if second not in successors[first]:
            successors[first].add(second)
            indegree[second] += 1
    pending = [node for node in nodes if indegree[node] == 0]
    seen = 0
    while pending:
        node = pending.pop()
        seen += 1
        for successor in successors[node]:
            indegree[successor] -= 1
            if indegree[successor] == 0:
                pending.append(successor)
    return seen == len(nodes)


def side_bits(words: tuple[int, ...], n: int) -> tuple[int, int]:
    patterns = minimal_missing(words, n)
    assert all(bits.bit_count() <= 1 for _, bits in patterns)
    epsilon_zero = epsilon_infinity = 0
    for source in graph_sources(words, n):
        positive_coordinates = tuple(i for i in range(n) if source >> i & 1)
        positive_rows = tuple(
            (support & source, bits)
            for support, bits in patterns
            if bits & source
        )
        core, _ = positive_core(positive_coordinates, positive_rows)
        epsilon_zero = max(epsilon_zero, int(bool(core)))

        negative_coordinates = tuple(i for i in range(n) if not source >> i & 1)
        negative_rows = tuple(
            (support, bits)
            for support, bits in patterns
            if bits and not bits & source and not (support ^ bits) & source
        )
        arcs = tuple(
            (bits.bit_length() - 1, coordinate)
            for support, bits in negative_rows
            for coordinate in negative_coordinates
            if (support ^ bits) >> coordinate & 1
        )
        epsilon_infinity = max(
            epsilon_infinity,
            int(not acyclic(negative_coordinates, arcs)),
        )
    return epsilon_zero, epsilon_infinity


def reflected_eta(words: tuple[int, ...], n: int):
    values = []
    for reflection in range(1 << n):
        reflected = tuple(sorted(word ^ reflection for word in words))
        family = set(reflected)
        if all(first | second in family for first in family for second in family):
            values.append((reflection, side_bits(reflected, n)))
    if not values:
        return values, ("infinity", "infinity")
    return values, tuple(min(bits[i] for _, bits in values) for i in range(2))


def replay_sign_report(report) -> dict[str, object]:
    expected_counts = {
        "retained_certificates_replayed": 2057,
        "positive_core_extensions": 1921,
        "split_systems": 2057,
        "joined_witnesses": 1804,
        "retained_row_families": 561,
        "example_orthant_decisions": 2754,
        "example_arrays": 162,
        "example_constructed_witnesses": 2448,
    }
    assert report["status"] == "passed" and report["counts"] == expected_counts
    geometry = report["geometry"]
    words = tuple(geometry["concepts"])
    assert words == (7, 15, 30, 31, 35, 39, 45, 47, 50, 51, 55, 57, 58, 59, 61, 62, 63)
    assert all(first | second in set(words) for first in words for second in words)
    patterns = minimal_missing(words, 6)
    assert patterns == {
        (record["support"], record["bits"])
        for record in geometry["minimal_patterns"]
    }
    assert tuple(sorted(graph_sources(words, 6))) == tuple(geometry["minimal_states"])
    assert geometry["vc_dimension"] == 2

    endpoint_statuses = {}
    for record in report["log_convexity_counterexample"]:
        rows = decode_rows(record["rows"])
        assert set(rows) == patterns
        assert all(
            all(value > 0 for value in row.values()) and sum(row.values()) == 1
            for row in rows.values()
        )
        if "witnesses" in record:
            witnesses = {
                int(word): tuple(Fraction(value) for value in vector)
                for word, vector in record["witnesses"].items()
            }
            assert set(witnesses) == set(words)
            assert all(valid(rows, word, vector) for word, vector in witnesses.items())
            endpoint_statuses[record["name"]] = True
        else:
            check_dual(
                orthant_system(6, rows, record["excluded_state"]),
                record["dual"],
            )
            endpoint_statuses[record["name"]] = False
    assert endpoint_statuses == {
        "left": True,
        "right": True,
        "geometric_midpoint": False,
    }
    return {
        "counts": expected_counts,
        "minimal_patterns": len(patterns),
        "endpoint_statuses": endpoint_statuses,
    }


def replay_scaling_report(report) -> dict[str, object]:
    expected_counts = {
        "marker_learning_parents": 107,
        "small_tail_witnesses": 390,
        "large_tail_witnesses": 379,
        "retained_classes": 141,
        "elementary_minor_type_checks": 1159,
        "minor_row_equivariance_checks": 3038,
        "deleted_positive_root_rows": 184,
        "fixed_tail_witnesses": 33,
        "ray_class_decisions": 12,
        "product_type_checks": 4,
    }
    assert report["status"] == "passed" and report["counts"] == expected_counts
    calibrations = {}
    for record in report["reflection_examples"]:
        words = tuple(record["concepts"])
        orientations, eta = reflected_eta(words, record["n"])
        encoded = [
            (item["reflection"], tuple(item["bits"]))
            for item in record["orientations"]
        ]
        assert orientations == encoded
        assert eta == tuple(record["eta"])
        calibrations[record["name"]] = list(eta)
    assert calibrations == {
        "gain17": [1, 0],
        "path_learning5": [0, 1],
        "daisy_star": [0, 0],
    }

    certificate_count = 0
    for record in report["fixed_tail_certificates"]:
        rows = decode_rows(record["rows"])
        witnesses = {
            int(word): tuple(Fraction(value) for value in vector)
            for word, vector in record["witnesses"].items()
        }
        assert set(witnesses) == set(record["concepts"])
        assert all(valid(rows, word, vector) for word, vector in witnesses.items())
        certificate_count += len(witnesses)
    assert certificate_count == expected_counts["fixed_tail_witnesses"]

    for item in report["mixed_ray"]:
        scale = Fraction(item["scale"])
        assert item["gain_class"] == (2 * scale * scale > 1)
        assert item["learning_class"] == (scale * scale < 4)
        assert item["product"] == (item["gain_class"] and item["learning_class"])
    assert report["mixed_ray_exact_interval"] == "1/sqrt(2) < scale < 2"
    return {
        "counts": expected_counts,
        "reflection_calibrations": calibrations,
        "fixed_tail_witnesses_replayed": certificate_count,
        "mixed_ray_decisions_replayed": len(report["mixed_ray"]),
    }


def run() -> None:
    sign_report = json.loads(SIGN_SOURCE.read_text(encoding="utf-8"))
    scaling_report = json.loads(SCALING_SOURCE.read_text(encoding="utf-8"))
    result = {
        "status": "passed",
        "source_reports": {
            SIGN_SOURCE.name: sha256(SIGN_SOURCE.read_bytes()).hexdigest(),
            SCALING_SOURCE.name: sha256(SCALING_SOURCE.read_bytes()).hexdigest(),
        },
        "source_sign_replay": replay_sign_report(sign_report),
        "scaling_replay": replay_scaling_report(scaling_report),
        "scope": "Exact replay of the 17-state endpoint witnesses and midpoint dual, fixed-tail witnesses, reflected indicator calibrations, and mixed-ray decisions. Corpus-wide counts are integrity-checked; universal topology and operation laws use the written proofs.",
        "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    run()
