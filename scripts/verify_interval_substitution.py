#!/usr/bin/env python3
"""Replay retained coordinate-interval and Boolean-substitution certificates."""

from __future__ import annotations

from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import json


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "data/examples/assets"
COORDINATE = ASSETS / "coordinate_interval_source.json"
SUBSTITUTION = ASSETS / "dominance_substitution_source.json"
STANDALONE = ASSETS / "dominance_substitution_replay_source.json"
INTEGRATION = ASSETS / "interval_substitution_integration_source.json"
OUTPUT = ASSETS / "interval_substitution_verification.json"

EXPECTED_SHA256 = {
    COORDINATE.name: "013e97f53f5d17724854c00ec1b223974af75e85ca7afb6160fc3431b5435823",
    SUBSTITUTION.name: "bf9044a8701bb363dadfd61130a1c907febf1fcd20fd7f37620c12b9fc7f1da0",
    STANDALONE.name: "db801a02dc5a2414f02a22ea7fcc976086a875bc2b33badb5e7c4a8bd02a88db",
    INTEGRATION.name: "6e7cb7636f17788bab736bb2f744439b4e1f7e23de2d7cd8694ce83accf71c5a",
}

SPECTRAL_MAXIMUM = (1, 3, 7, 14, 15, 16, 17, 19, 23, 24, 25, 27, 28, 29, 30, 31)


def missing(concepts, n: int):
    concepts = tuple(concepts)
    answer = []
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
                answer.append((support, bits))
            if bits == 0:
                break
            bits = (bits - 1) & support
    return set(answer)


def interval_ratio(row, pattern, concept: int, magnitudes):
    favorable = sum(
        value * magnitudes[coordinate]
        for coordinate, value in row.items()
        if (pattern[1] ^ concept) >> coordinate & 1
    )
    unfavorable = sum(
        value * magnitudes[coordinate]
        for coordinate, value in row.items()
        if not (pattern[1] ^ concept) >> coordinate & 1
    )
    assert favorable > 0
    return unfavorable / favorable


def replay_coordinate(data):
    assert data["status"] == "passed"
    assert tuple(data["concepts"]) == SPECTRAL_MAXIMUM
    rows = {
        (record["support"], record["bits"]): {
            int(coordinate): F(value)
            for coordinate, value in record["coefficients"].items()
        }
        for record in data["rows"]
    }
    assert set(rows) == missing(SPECTRAL_MAXIMUM, 5)
    assert all(value > 0 for row in rows.values() for value in row.values())
    centers = {
        int(concept): tuple(map(F, values))
        for concept, values in data["centers"].items()
    }
    assert set(centers) == set(SPECTRAL_MAXIMUM)
    endpoint_factor = F(data["endpoint_factor"])
    threshold = 1 / endpoint_factor
    assert endpoint_factor > 1
    maximum = F(0)
    endpoint_count = row_count = 0
    for concept, center in centers.items():
        assert all(value > 0 for value in center)
        assert all(
            interval_ratio(row, pattern, concept, center) < threshold
            for pattern, row in rows.items()
        )
        for coordinate in range(5):
            for factor in (1 / endpoint_factor, endpoint_factor):
                endpoint = tuple(
                    value * factor if index == coordinate else value
                    for index, value in enumerate(center)
                )
                for pattern, row in rows.items():
                    ratio = interval_ratio(row, pattern, concept, endpoint)
                    assert ratio <= endpoint_factor * interval_ratio(
                        row, pattern, concept, center
                    )
                    assert ratio < 1
                    maximum = max(maximum, ratio)
                    row_count += 1
                endpoint_count += 1
    assert all(concept ^ 4 not in SPECTRAL_MAXIMUM for concept in (14, 16, 1))
    assert {(7, 2), (14, 4), (28, 8)} <= set(rows)
    assert len(SPECTRAL_MAXIMUM) == 16
    for support in range(32):
        traces = {concept & support for concept in SPECTRAL_MAXIMUM}
        if support.bit_count() == 2:
            assert len(traces) == 4
        if support.bit_count() == 3:
            assert len(traces) < 8
    checks = {
        "centers": len(centers),
        "endpoints": endpoint_count,
        "endpoint_row_checks": row_count,
        "missing_patterns": len(rows),
        "maximum_endpoint_ratio": str(maximum),
        "interval_ratio": str(endpoint_factor * endpoint_factor),
    }
    assert checks == data["checks"]
    return checks


def check_dominance(concepts, n: int, rows, witnesses, minimal=True):
    assert set(witnesses) == set(concepts)
    if minimal:
        assert set(rows) == missing(concepts, n)
    comparisons = 0
    for concept in concepts:
        assert len(witnesses[concept]) == n
        for pattern, entries in rows.items():
            assert set(entries) == {
                coordinate
                for coordinate in range(n)
                if pattern[0] >> coordinate & 1
            }
            favorable = [
                weight + witnesses[concept][coordinate]
                for coordinate, weight in entries.items()
                if (concept ^ pattern[1]) >> coordinate & 1
            ]
            unfavorable = [
                weight + witnesses[concept][coordinate]
                for coordinate, weight in entries.items()
                if not (concept ^ pattern[1]) >> coordinate & 1
            ]
            assert favorable
            if unfavorable:
                assert max(favorable) > max(unfavorable)
            comparisons += 1
    return comparisons


def decode(record):
    rows = {
        tuple(row["pattern"]): {
            int(coordinate): F(value)
            for coordinate, value in row["entries"].items()
        }
        for row in record["rows"]
    }
    witnesses = {
        int(concept): tuple(map(F, values))
        for concept, values in record["witnesses"].items()
    }
    check_dominance(record["concepts"], record["n"], rows, witnesses)
    return tuple(record["concepts"]), record["n"], rows, witnesses


def drop(word: int, coordinate: int) -> int:
    return (word & ((1 << coordinate) - 1)) | (
        (word >> (coordinate + 1)) << coordinate
    )


def replay_substitution(data):
    assert data["status"] == "passed" and len(data["cases"]) == 5
    totals = {
        "cases": 0,
        "output_concepts": 0,
        "output_rows": 0,
        "output_comparisons": 0,
        "extended_comparisons": 0,
    }
    for case in data["cases"]:
        outer, n, _, _ = decode(case["outer"])
        test, k, _, _ = decode(case["test"])
        complement, complement_n, _, _ = decode(case["test_complement"])
        assert complement_n == k
        assert set(complement) == set(range(1 << k)) - set(test)
        output, output_n, output_rows, _ = decode(case["output"])
        coordinate = case["coordinate"]
        expected = {
            drop(concept, coordinate) | (word << (n - 1))
            for concept in outer
            for word in range(1 << k)
            if concept >> coordinate & 1 == (0 if word in test else 1)
        }
        assert output_n == n - 1 + k and set(output) == expected

        extended = case["extended"]
        extended_rows = {
            tuple(row["pattern"]): {
                int(index): F(value) for index, value in row["entries"].items()
            }
            for row in extended["rows"]
        }
        extended_witnesses = {
            int(concept): tuple(map(F, values))
            for concept, values in extended["witnesses"].items()
        }
        check_dominance(
            extended["concepts"],
            extended["n"],
            extended_rows,
            extended_witnesses,
            minimal=False,
        )
        assert set(extended["concepts"]) == {
            concept
            for concept in range(1 << extended["n"])
            if all(concept & pattern[0] != pattern[1] for pattern in extended_rows)
        }
        assert {
            concept & ((1 << output_n) - 1) for concept in extended["concepts"]
        } == set(output)
        totals["cases"] += 1
        totals["output_concepts"] += len(output)
        totals["output_rows"] += len(output_rows)
        totals["output_comparisons"] += len(output) * len(output_rows)
        totals["extended_comparisons"] += len(extended["concepts"]) * len(
            extended_rows
        )
    assert totals == data["checks"]
    return totals


def replay():
    for path in (COORDINATE, SUBSTITUTION, STANDALONE, INTEGRATION):
        assert sha256(path.read_bytes()).hexdigest() == EXPECTED_SHA256[path.name]
    coordinate = json.loads(COORDINATE.read_text(encoding="utf-8"))
    substitution = json.loads(SUBSTITUTION.read_text(encoding="utf-8"))
    standalone = json.loads(STANDALONE.read_text(encoding="utf-8"))
    integration = json.loads(INTEGRATION.read_text(encoding="utf-8"))
    coordinate_checks = replay_coordinate(coordinate)
    substitution_checks = replay_substitution(substitution)
    assert standalone["status"] == "replayed"
    assert standalone["checks"] == substitution_checks
    assert standalone["sha256"]["dominance_substitution_verification.json"] == EXPECTED_SHA256[SUBSTITUTION.name]
    assert integration["status"] == "passed"
    assert integration["sha256"]["coordinate_substitution.tex"] == "6b22e8aee0c585e885a9a2e5f5f7fb1069f8ed1016ed3e2e32f2fc13a50050b3"
    assert integration["new_diagnostics_regenerated"] == 2
    assert integration["prior_diagnostics_retained_by_hash"] == 20
    return {
        "status": "passed",
        "source_sha256": EXPECTED_SHA256,
        "coordinate_interval_checks": coordinate_checks,
        "substitution_checks": substitution_checks,
        "standalone_replay_matched": True,
        "integration": {
            "pages": integration["pages"],
            "labels": integration["labels"],
            "citations": integration["citations"],
            "prior_labels_preserved": integration["prior_labels_preserved"],
        },
        "scope": {
            "coordinate_interval": coordinate["scope"],
            "substitution": substitution["scope"],
            "integration": integration["scope"],
        },
    }


def main():
    report = replay()
    OUTPUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
