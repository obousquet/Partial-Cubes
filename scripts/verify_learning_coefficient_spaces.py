#!/usr/bin/env python3
"""Replay exact learning-space coefficient certificates stored by the database."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import json


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "data/examples/assets"
SOURCE_REPORT = ASSETS / "learning_space_coefficients_source.json"
OUTPUT = ASSETS / "learning_coefficient_spaces_verification.json"


def sign(word: int, coordinate: int) -> int:
    return 1 if word >> coordinate & 1 else -1


def dot(first, second):
    return sum(x * y for x, y in zip(first, second))


def minimal_missing(words: tuple[int, ...], n: int) -> set[tuple[int, int]]:
    present = set(words)
    patterns: set[tuple[int, int]] = set()
    for support in range(1, 1 << n):
        bits = support
        while True:
            realized = any((word & support) == bits for word in present)
            all_deletions_realized = all(
                any(
                    (word & (support ^ (1 << coordinate)))
                    == (bits & ~(1 << coordinate))
                    for word in present
                )
                for coordinate in range(n)
                if support >> coordinate & 1
            )
            if not realized and all_deletions_realized:
                patterns.add((support, bits))
            if bits == 0:
                break
            bits = (bits - 1) & support
    return patterns


def is_learning_space(words: tuple[int, ...], n: int) -> bool:
    family = set(words)
    return (
        0 in family
        and all(first | second in family for first in family for second in family)
        and all(
            word == 0
            or any(
                word ^ (1 << coordinate) in family
                for coordinate in range(n)
                if word >> coordinate & 1
            )
            for word in family
        )
    )


def decode_rows(rows):
    return {
        (row["support"], row["bits"]): {
            int(coordinate): Fraction(value)
            for coordinate, value in row["coefficients"].items()
        }
        for row in rows
    }


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
        for (support, bits), row in sorted(rows.items())
    )
    return tuple(system)


def check_dual(system, encoded_weights) -> None:
    weights = tuple(Fraction(value) for value in encoded_weights)
    assert len(weights) == len(system)
    assert any(weights) and all(value >= 0 for value in weights)
    dimension = len(system[0]) if system else 0
    assert all(
        sum(weight * row[i] for weight, row in zip(weights, system)) == 0
        for i in range(dimension)
    )


def check_full_certificate(words, n, rows, record) -> None:
    margin = Fraction(record["margin"])
    assert margin > 0
    witnesses = {
        int(word): tuple(Fraction(value) for value in vector)
        for word, vector in record["witnesses"].items()
    }
    assert set(witnesses) == set(words)
    for word, vector in witnesses.items():
        assert len(vector) == n
        assert all(sign(word, i) * vector[i] >= margin for i in range(n))
        for (_, bits), row in rows.items():
            separator = sum(sign(bits, i) * coefficient * vector[i] for i, coefficient in row.items())
            assert separator <= -margin
    root = tuple(Fraction(value) for value in record["root_vector"])
    assert len(root) == n
    assert all(dot(row, root) > 0 for row in orthant_system(n, rows, 0))


def run() -> None:
    report = json.loads(SOURCE_REPORT.read_text(encoding="utf-8"))
    expected_coverage = {
        "nonempty_classes_screened": 274,
        "learning_spaces_through_Q3": 44,
        "prescribed_row_families": 181,
        "independent_orthant_decisions": 755,
        "full_positive_certificates": 178,
        "root_infeasibility_certificates": 3,
    }
    assert report["status"] == "passed"
    assert report["coverage"] == expected_coverage
    assert len(report["cases"]) == expected_coverage["prescribed_row_families"]

    feasible = infeasible = orthant_decisions = 0
    cyclic_statuses = {}
    for record in report["cases"]:
        words = tuple(record["words"])
        n = record["n"]
        rows = decode_rows(record["rows"])
        assert is_learning_space(words, n)
        assert set(rows) == minimal_missing(words, n)
        assert all(
            set(row) == {i for i in range(n) if support >> i & 1}
            and all(value > 0 for value in row.values())
            and sum(row.values()) == 1
            for (support, _), row in rows.items()
        )
        statuses = record["independent_orthants_feasible"]
        assert len(statuses) == len(words)
        assert record["root_feasible"] == all(statuses)
        orthant_decisions += len(statuses)
        if record["root_feasible"]:
            check_full_certificate(words, n, rows, record)
            feasible += 1
        else:
            check_dual(
                orthant_system(n, rows, 0),
                record["root_infeasibility_multipliers"],
            )
            infeasible += 1
        if record["name"].startswith("cyclic_"):
            cyclic_statuses[record["name"]] = record["root_feasible"]

    assert feasible == 178 and infeasible == 3 and orthant_decisions == 755
    assert cyclic_statuses == {
        "cyclic_1/2_1/2": True,
        "cyclic_100_1/1000": True,
        "cyclic_1_1": False,
        "cyclic_2_2": False,
        "cyclic_10_1/10": False,
    }

    control = report["inaccessible_control"]
    control_rows = decode_rows(control["rows"])
    control_root = tuple(Fraction(value) for value in control["root_vector"])
    assert all(dot(row, control_root) > 0 for row in orthant_system(2, control_rows, 0))
    check_dual(
        orthant_system(2, control_rows, control["failed_concept"]),
        control["infeasibility_multipliers"],
    )

    result = {
        "status": "passed",
        "source_report_sha256": sha256(SOURCE_REPORT.read_bytes()).hexdigest(),
        "coverage": expected_coverage,
        "serialized_cases_replayed": len(report["cases"]),
        "cyclic_domain_checks": cyclic_statuses,
        "inaccessible_control_replayed": True,
        "scope": "Exact replay of every serialized primal and dual certificate in the preserved source report. The general criterion and Euclidean topology use the written proofs.",
        "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    run()
