#!/usr/bin/env python3
"""Replay the exact finite evidence for the source-spectral invariant."""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import json


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "data/examples/assets"
SPECTRA = ASSETS / "source_spectra_source.json"
SPECTRA_REPLAY = ASSETS / "source_spectra_replay_source.json"
ZERO = ASSETS / "spectral_zero_source.json"
ZERO_REPLAY = ASSETS / "spectral_zero_replay_source.json"
OUTPUT = ASSETS / "spectral_invariant_verification.json"


def minimal_missing(words: tuple[int, ...], n: int) -> set[tuple[int, int]]:
    family = set(words)
    answer: set[tuple[int, int]] = set()
    for support in range(1, 1 << n):
        bits = support
        while True:
            if not any(word & support == bits for word in family) and all(
                any(
                    word & (support ^ (1 << coordinate))
                    == bits & ~(1 << coordinate)
                    for word in family
                )
                for coordinate in range(n)
                if support >> coordinate & 1
            ):
                answer.add((support, bits))
            if bits == 0:
                break
            bits = (bits - 1) & support
    return answer


def sources(words: tuple[int, ...], n: int) -> tuple[int, ...]:
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


def determinant(matrix) -> Fraction:
    work = [list(map(Fraction, row)) for row in matrix]
    result = Fraction(1)
    for index in range(len(work)):
        pivot = next(
            (row for row in range(index, len(work)) if work[row][index]),
            None,
        )
        if pivot is None:
            return Fraction(0)
        if pivot != index:
            work[index], work[pivot] = work[pivot], work[index]
            result = -result
        value = work[index][index]
        result *= value
        for row in range(index + 1, len(work)):
            ratio = work[row][index] / value
            for column in range(index + 1, len(work)):
                work[row][column] -= ratio * work[index][column]
    return result


def radius_compare_one(matrix) -> int:
    """Return -1, 0, or 1 according as a nonnegative matrix has radius <, =, or > 1."""
    n = len(matrix)
    for size in range(1, n + 1):
        for indices in combinations(range(n), size):
            minor = tuple(
                tuple(Fraction(i == j) - matrix[i][j] for j in indices)
                for i in indices
            )
            if determinant(minor) < 0:
                return 1
    full = tuple(
        tuple(Fraction(i == j) - matrix[i][j] for j in range(n))
        for i in range(n)
    )
    return 0 if determinant(full) == 0 else -1


def row_choices(rows, positive: bool):
    n = len(rows)
    if positive:
        for mask in range(1, 1 << n):
            indices = tuple(i for i in range(n) if mask >> i & 1)
            choices = [
                [
                    row
                    for row in rows[i]
                    if all(not row[j] for j in range(n) if not mask >> j & 1)
                ]
                for i in indices
            ]
            if all(choices):
                for selected in product(*choices):
                    yield tuple(
                        tuple(row[j] for j in indices) for row in selected
                    )
    else:
        choices = [options or [(Fraction(0),) * n] for options in rows]
        yield from product(*choices)


def source_rows(n: int, coefficients, concept: int, positive: bool):
    indices = [
        coordinate
        for coordinate in range(n)
        if bool(concept >> coordinate & 1) == positive
    ]
    rows = [[] for _ in indices]
    for (support, bits), entries in coefficients.items():
        if not bits:
            continue
        root = bits.bit_length() - 1
        applies = bool(bits & concept) if positive else not bool(support & concept)
        if applies:
            rows[indices.index(root)].append(
                tuple(
                    entries.get(coordinate, Fraction(0)) / entries[root]
                    if coordinate != root
                    else Fraction(0)
                    for coordinate in indices
                )
            )
    return rows


def square_scaled(matrix, squared_radius: Fraction):
    n = len(matrix)
    return tuple(
        tuple(
            sum(matrix[i][k] * matrix[k][j] for k in range(n))
            / squared_radius
            for j in range(n)
        )
        for i in range(n)
    )


def parsed_rows(records):
    return {
        tuple(record["pattern"]): {
            int(coordinate): Fraction(value)
            for coordinate, value in record["entries"].items()
        }
        for record in records
    }


def coefficient_array(words: tuple[int, ...], n: int, ratio: Fraction):
    result = {}
    for support, bits in minimal_missing(words, n):
        root = bits.bit_length() - 1 if bits else None
        result[(support, bits)] = {
            coordinate: (
                Fraction(1)
                if root is None or coordinate == root
                else ratio
            )
            for coordinate in range(n)
            if support >> coordinate & 1
        }
    return result


def dominance_gap(coefficients, concept: int, magnitudes) -> Fraction:
    gaps = []
    for (_, bits), row in coefficients.items():
        favorable = [
            value * magnitudes[coordinate]
            for coordinate, value in row.items()
            if bool(concept >> coordinate & 1) != bool(bits >> coordinate & 1)
        ]
        unfavorable = [
            value * magnitudes[coordinate]
            for coordinate, value in row.items()
            if bool(concept >> coordinate & 1) == bool(bits >> coordinate & 1)
        ]
        assert favorable
        gap = max(favorable) - max(unfavorable, default=Fraction(0))
        assert gap > 0
        gaps.append(gap)
    return min(gaps, default=Fraction(1))


def replay_spectra(report) -> dict[str, object]:
    expected_counts = {
        "abstract_side_tests": 24,
        "policy_matrices": 329,
        "source_side_tests": 80,
        "full_orthant_decisions": 196,
        "minor_array_maps": 192,
        "separate_side_implications": 384,
        "nonvacuous_minor_side_implications": 276,
        "minor_maps_from_infeasible_arrays": 108,
        "block_product_side_tests": 12,
    }
    assert report["status"] == "passed" and report["counts"] == expected_counts
    assert len(report["abstract"]) == 12
    assert len(report["minor_maps"]) == 192
    assert len(report["block_products"]) == 12

    maximum = next(item for item in report["classes"] if item["name"] == "maximum16")
    words = tuple(maximum["concepts"])
    n = maximum["n"]
    coefficients = parsed_rows(maximum["base_rows"])
    assert set(coefficients) == minimal_missing(words, n)
    assert all(first | second in set(words) for first in words for second in words)
    union_closed_reflections = [
        reflection
        for reflection in range(1 << n)
        if all(
            (first ^ reflection) | (second ^ reflection)
            in {word ^ reflection for word in words}
            for first in words
            for second in words
        )
    ]
    assert union_closed_reflections == [0]

    extrema = []
    for positive, target in ((True, Fraction(9, 8)), (False, Fraction(9, 16))):
        comparisons = []
        for concept in sources(words, n):
            for matrix in row_choices(source_rows(n, coefficients, concept, positive), positive):
                comparisons.append(radius_compare_one(square_scaled(matrix, target)))
        assert comparisons and 0 in comparisons
        assert all(value >= 0 for value in comparisons) if positive else all(
            value <= 0 for value in comparisons
        )
        extrema.append(
            {
                "positive": positive,
                "squared_extremum": str(target),
                "matrices": len(comparisons),
                "comparison_counts": {str(k): v for k, v in Counter(comparisons).items()},
            }
        )
    assert extrema == report["maximum16_extrema"]
    assert report["maximum16_calibration_ratio_squared"] == "1/2"
    return {
        "counts": expected_counts,
        "maximum16_sources": list(sources(words, n)),
        "maximum16_union_closed_reflections": union_closed_reflections,
        "maximum16_extrema": extrema,
    }


def replay_zero(report) -> dict[str, object]:
    assert report["status"] == "passed"
    assert report["counts"] == {
        "full_dominance_witnesses": 47,
        "power_calibrations": 4,
        "balanced_obstructions": 2,
    }
    gain = (7, 15, 30, 31, 35, 39, 45, 47, 50, 51, 55, 57, 58, 59, 61, 62, 63)
    learning = tuple(
        sorted(
            {
                ((1 << initial) - 1) | (31 ^ ((1 << final) - 1))
                for initial in range(6)
                for final in range(6)
            }
        )
    )
    cases = {
        "gain17": (gain, 6, Fraction(2)),
        "learning16": (learning, 5, Fraction(1, 2)),
        "zero_coordinates": ((0,), 0, Fraction(1)),
        "zero_state": ((0,), 2, Fraction(1)),
        "one_state": ((3,), 2, Fraction(1)),
        "full_square": ((0, 1, 2, 3), 2, Fraction(1)),
    }
    checked = 0
    for name, (words, n, ratio) in cases.items():
        coefficients = coefficient_array(tuple(words), n, ratio)
        for witness in report["reconstruction"][name]["full_witnesses"]:
            magnitudes = tuple(Fraction(value) for value in witness["magnitudes"])
            gap = dominance_gap(coefficients, witness["concept"], magnitudes)
            assert gap == Fraction(witness["minimum_gap"])
            checked += 1

    repaired_words = tuple(report["direct_extension_failure"]["concepts"])
    repaired_coefficients = coefficient_array(repaired_words, 3, Fraction(1, 2))
    for witness in report["repaired_extension"]["full_witnesses"]:
        magnitudes = tuple(Fraction(value) for value in witness["magnitudes"])
        gap = dominance_gap(repaired_coefficients, witness["concept"], magnitudes)
        assert gap == Fraction(witness["minimum_gap"])
        checked += 1
    assert checked == 47

    for calibration in report["power_calibrations"]:
        exponent = calibration["t"]
        assert Fraction(calibration["alpha_G_squared"]) == 2 * Fraction(4) ** exponent
        assert calibration["beta_G"] == 0
        assert calibration["alpha_L"] == "infinity"
        assert Fraction(calibration["beta_L_squared"]) == 2 * Fraction(4) ** (-exponent)
        assert Fraction(calibration["product_ratio"]) == Fraction(4) ** (-exponent)

    for obstruction in report["maximum_balanced_obstructions"]:
        total = Counter()
        for comparison in obstruction:
            total.update(comparison)
        assert all(value == 0 for value in total.values())
    return {
        "counts": report["counts"],
        "dominance_witnesses_replayed": checked,
        "product_ratios": [item["product_ratio"] for item in report["power_calibrations"]],
    }


def run() -> None:
    spectra = json.loads(SPECTRA.read_text(encoding="utf-8"))
    spectra_replay = json.loads(SPECTRA_REPLAY.read_text(encoding="utf-8"))
    zero = json.loads(ZERO.read_text(encoding="utf-8"))
    zero_replay = json.loads(ZERO_REPLAY.read_text(encoding="utf-8"))
    assert spectra_replay["status"] == "passed"
    assert spectra_replay["report_bytes_match"] is True
    assert spectra_replay["report_sha256"] == sha256(SPECTRA.read_bytes()).hexdigest()
    assert zero_replay["status"] == "passed"
    assert zero_replay["report_byte_match"] is True
    assert zero_replay["report_sha256"] == sha256(ZERO.read_bytes()).hexdigest()
    result = {
        "status": "passed",
        "source_reports": {
            path.name: sha256(path.read_bytes()).hexdigest()
            for path in (SPECTRA, SPECTRA_REPLAY, ZERO, ZERO_REPLAY)
        },
        "source_spectra_replay": replay_spectra(spectra),
        "spectral_zero_replay": replay_zero(zero),
        "scope": "Exact replay of the maximum-class spectral extrema, its unique union-closed orientation, all 47 strict-dominance witnesses, four nonattainment calibrations, and both balanced obstructions. Complete source-report counts and cold-replay hashes are integrity-checked; universal laws use the written proofs.",
        "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    run()
