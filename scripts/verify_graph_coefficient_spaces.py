#!/usr/bin/env python3
"""Replay exact finite evidence for graph-model coefficient spaces."""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "data/examples/assets"
SPHERES = ASSETS / "coefficient_cycle_spheres_source.json"
PRIVATE = ASSETS / "private_coordinate_products_source.json"
ARRANGEMENT = ASSETS / "cycle_arrangement_source.json"
UNION = ASSETS / "union_closed_graph_spaces_source.json"
OUTPUT = ASSETS / "graph_coefficient_spaces_verification.json"


def traces(words: tuple[int, ...], n: int) -> list[set[int]]:
    return [{word & support for word in words} for support in range(1 << n)]


def minimal_missing(words: tuple[int, ...], n: int) -> set[tuple[int, int]]:
    realized = traces(words, n)
    answer = set()
    for support in range(1, 1 << n):
        bits = support
        while True:
            if bits not in realized[support] and all(
                bits & ~(1 << coordinate) in realized[support ^ (1 << coordinate)]
                for coordinate in range(n)
                if support >> coordinate & 1
            ):
                answer.add((support, bits))
            if bits == 0:
                break
            bits = (bits - 1) & support
    return answer


def vc_dimension(words: tuple[int, ...], n: int) -> int:
    realized = traces(words, n)
    return max(
        support.bit_count()
        for support in range(1 << n)
        if len(realized[support]) == 1 << support.bit_count()
    )


def decode_rows(records):
    rows = {}
    for record in records:
        encoded = record["coefficients"]
        entries = encoded.items() if isinstance(encoded, dict) else encoded
        rows[(record["support"], record["bits"])] = {
            int(coordinate): Fraction(value) for coordinate, value in entries
        }
    return rows


def valid(rows, concept: int, vector: tuple[Fraction, ...]) -> bool:
    if not all(
        (1 if concept >> coordinate & 1 else -1) * value > 0
        for coordinate, value in enumerate(vector)
    ):
        return False
    return all(
        sum(
            (1 if bits >> coordinate & 1 else -1)
            * coefficient
            * vector[coordinate]
            for coordinate, coefficient in row.items()
        )
        < 0
        for (_, bits), row in rows.items()
    )


def gf2_rank(rows: list[int]) -> int:
    pivots: dict[int, int] = {}
    for value in rows:
        while value:
            pivot = value.bit_length() - 1
            if pivot not in pivots:
                pivots[pivot] = value
                break
            value ^= pivots[pivot]
    return len(pivots)


def nerve_data(facets: list[list[int]]):
    faces: dict[int, set[tuple[int, ...]]] = {}
    for facet in facets:
        for size in range(1, len(facet) + 1):
            faces.setdefault(size - 1, set()).update(combinations(facet, size))
    dimension = max(faces)
    ranks = [0]
    for dim in range(1, dimension + 1):
        lower = {face: index for index, face in enumerate(sorted(faces[dim - 1]))}
        columns = []
        for face in sorted(faces[dim]):
            column = 0
            for omitted in range(len(face)):
                boundary = face[:omitted] + face[omitted + 1 :]
                column ^= 1 << lower[boundary]
            columns.append(column)
        ranks.append(gf2_rank(columns))
    ranks.append(0)
    f_vector = [len(faces[dim]) for dim in range(dimension + 1)]
    betti = [
        f_vector[dim] - ranks[dim] - ranks[dim + 1]
        for dim in range(dimension + 1)
    ]
    return f_vector, ranks, betti


def replay_spheres(report) -> dict[str, object]:
    assert report["status"] == "passed" and len(report["cases"]) == 2
    witness_count = 0
    cases = []
    for record in report["cases"]:
        m = record["cycle_length"]
        n = record["n"]
        words = tuple(record["words"])
        patterns = minimal_missing(words, n)
        assert len(words) == record["cardinality"] == 8 * m + 8
        assert vc_dimension(words, n) == record["VC"] == 3
        assert patterns == {tuple(pattern) for pattern in record["missing_patterns"]}
        assert max(support.bit_count() for support, _ in patterns) == 4
        dimension = sum(support.bit_count() - 1 for support, _ in patterns)
        assert dimension == 8 * m + m * (m - 1) // 2
        certificate = record["certificate"]
        rows = decode_rows(certificate["rows"])
        assert set(rows) == patterns
        witnesses = {
            int(word): tuple(Fraction(value) for value in vector)
            for word, vector in certificate["witnesses"].items()
        }
        assert set(witnesses) == set(words)
        assert all(valid(rows, word, vector) for word, vector in witnesses.items())
        witness_count += len(witnesses)
        cases.append(
            {
                "cycle_length": m,
                "coordinates": n,
                "concepts": len(words),
                "patterns": len(patterns),
                "coefficient_dimension": dimension,
            }
        )
    nerve = report["m4_nerve"]
    f_vector, ranks, betti = nerve_data(nerve["facets"])
    assert f_vector == nerve["f_vector"]
    assert ranks == nerve["boundary_ranks"]
    assert betti == nerve["betti_F2"]
    return {
        "cases": cases,
        "separator_witnesses_replayed": witness_count,
        "m4_nerve": {
            "f_vector": f_vector,
            "boundary_ranks": ranks,
            "betti_F2": betti,
        },
    }


def replay_private(report) -> dict[str, object]:
    expected = {
        "full_reconstructed_witnesses": 786,
        "private_array_changes": 12,
        "radial_inverse_checks": 20,
    }
    assert report["status"] == "passed" and report["counts"] == expected
    assert [item["name"] for item in report["classes"]] == [
        "bipartite_cycle_4",
        "bipartite_cycle_6",
        "component_one",
        "triangle_section",
    ]
    for item in report["radial_sections"]:
        source_a, source_b, target_a, target_b = map(Fraction, item["base"])
        p, q = map(Fraction, item["section"])
        assert p > 0 and q > 0
        assert source_a * p + source_b * q > 1
        assert target_a * p + target_b * q < 1
    return {"counts": expected, "radial_sections_replayed": len(report["radial_sections"])}


def replay_arrangement(report) -> dict[str, object]:
    assert report["status"] == "passed" and len(report["cases"]) == 21
    parameters = feasible = centers = 0
    for case in report["cases"]:
        n = case["n"]
        edges = [tuple(edge) for edge in case["edges"]]
        cycles = case["cycles"]
        for item in case["parameters"]:
            d = item["d"]
            expected = not any(len({d[vertex] for vertex in cycle}) == 1 for cycle in cycles)
            assert item["feasible"] is expected
            parameters += 1
            if expected:
                feasible += 1
                eta = Fraction(item["eta"])
                center = tuple(Fraction(value) for value in item["center"])
                assert len(center) == n and eta > 0
                assert all(
                    center[target] - center[source]
                    > -abs(Fraction(d[target]) - Fraction(d[source]))
                    for source, target in edges
                )
                centers += 1
    summary = report["summary"]
    assert parameters == summary["parameter_points"] == 2305
    assert feasible == summary["feasible"] == 2275
    assert summary["rational_inverse_checks"] == 9100
    actual = report["actual_class"]
    words = tuple(actual["words"])
    assert vc_dimension(words, actual["n"]) == actual["VC"] == 3
    assert set(actual["shattered"]) == set(actual["strongly_shattered"])
    return {
        "digraph_cases": len(report["cases"]),
        "parameters_replayed": parameters,
        "feasible_centers_replayed": centers,
        "reported_radial_inverse_checks": summary["rational_inverse_checks"],
    }


def replay_union(report) -> dict[str, object]:
    expected = {
        "complete_pattern_lists": 2,
        "arc_conditioning_checks": 5,
        "full_orthant_decisions": 536,
        "full_array_decisions": 7,
        "feasible_arrays": 4,
        "infeasible_arrays": 3,
        "rational_section_checks": 4,
    }
    assert report["status"] == "passed" and report["counts"] == expected
    structures = []
    for record in report["structures"]:
        words = tuple(record["concepts"])
        n = record["coordinates"]
        patterns = minimal_missing(words, n)
        assert len(words) == record["cardinality"]
        assert all(first | second in set(words) for first in words for second in words)
        assert len(patterns) == record["patterns"]
        assert vc_dimension(words, n) == record["vc_dimension"] == 4
        assert sum(support.bit_count() - 1 for support, _ in patterns) == record["coefficient_dimension"]
        assert record["fixed_orientation_side_indicators"] == [1, 1]
        structures.append(
            {
                "vertices": record["vertices"],
                "arcs": len(record["arcs"]),
                "coordinates": n,
                "concepts": len(words),
                "patterns": len(patterns),
                "coefficient_dimension": record["coefficient_dimension"],
            }
        )
    witness_count = 0
    for certificate in report["certificates"]:
        if not certificate["actual"]:
            continue
        rows = decode_rows(certificate["rows"])
        witnesses = {
            int(word): tuple(Fraction(value) for value in vector)
            for word, vector in certificate["witnesses"].items()
        }
        assert set(witnesses) == set(certificate["concepts"])
        assert all(valid(rows, word, vector) for word, vector in witnesses.items())
        witness_count += len(witnesses)
    assert witness_count == 332
    return {
        "counts": expected,
        "structures": structures,
        "feasible_array_witnesses_replayed": witness_count,
    }


def run() -> None:
    reports = {
        path.name: json.loads(path.read_text(encoding="utf-8"))
        for path in (SPHERES, PRIVATE, ARRANGEMENT, UNION)
    }
    result = {
        "status": "passed",
        "source_reports": {
            path.name: sha256(path.read_bytes()).hexdigest()
            for path in (SPHERES, PRIVATE, ARRANGEMENT, UNION)
        },
        "sphere_replay": replay_spheres(reports[SPHERES.name]),
        "private_coordinate_replay": replay_private(reports[PRIVATE.name]),
        "arrangement_replay": replay_arrangement(reports[ARRANGEMENT.name]),
        "union_closed_replay": replay_union(reports[UNION.name]),
        "scope": "Exact replay of the m=4,6 sphere certificates and nerve homology, 2,305 arrangement parameters and centers, four feedback sections, both union-closed structures, and 332 feasible-array witnesses. Complete report counts are checked; general product and topology theorems use the written proofs.",
        "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    OUTPUT.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    run()
