#!/usr/bin/env python3
"""Replay exact relative-gluing certificates and parameter-domain nerves."""

from __future__ import annotations

from fractions import Fraction as F
from functools import cache
from hashlib import sha256
from itertools import combinations
from math import gcd, lcm
from pathlib import Path
import json


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "data/examples/assets"
GLUING = ASSETS / "relative_gluing_source.json"
NERVE = ASSETS / "relative_nerve_source.json"
INTEGRATION = ASSETS / "relative_integration_source.json"
OUTPUT = ASSETS / "relative_gluing_database_verification.json"

EXPECTED_SHA256 = {
    GLUING.name: "3cfa616232553c034f4bcbbf676f4079554d1049c02952b92a45999d0789bcf7",
    NERVE.name: "379f430bef7a6abf0d9a370c45cef8267caf8b9f026ae98543b4ef2efa8bf5fc",
    INTEGRATION.name: "848ede4769e5dd8943a3c08c93c9bc31409d259a3941a241c02e8f607a08227e",
}


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


def sign(word: int, coordinate: int) -> int:
    return 1 if word >> coordinate & 1 else -1


def decode_certificate(packet):
    n = packet["n"]
    concepts = tuple(packet["words"])
    margin = F(packet["margin"])
    rows = {
        (row["support"], row["bits"]): {
            int(coordinate): F(value)
            for coordinate, value in row["coefficients"].items()
        }
        for row in packet["rows"]
    }
    witnesses = {
        int(concept): tuple(map(F, values))
        for concept, values in packet["witnesses"].items()
    }
    assert 0 < margin <= 1
    assert set(witnesses) == set(concepts)
    assert set(rows) == missing(concepts, n)
    for concept, witness in witnesses.items():
        assert len(witness) == n
        assert all(
            margin <= sign(concept, coordinate) * witness[coordinate] <= 1
            for coordinate in range(n)
        )
    comparisons = 0
    for pattern, row in rows.items():
        support, bits = pattern
        assert set(row) == {
            coordinate for coordinate in range(n) if support >> coordinate & 1
        }
        assert all(value > 0 for value in row.values())
        assert sum(row.values()) == 1
        for concept, witness in witnesses.items():
            value = sum(
                coefficient * sign(bits, coordinate) * witness[coordinate]
                for coordinate, coefficient in row.items()
            )
            assert value <= -margin
            comparisons += 1
    return {
        "concepts": len(concepts),
        "rows": len(rows),
        "comparisons": comparisons,
        "margin": str(margin),
    }


def unique_solution(matrix, vector, variables: int):
    rows = [list(map(F, row)) + [F(value)] for row, value in zip(matrix, vector)]
    pivot = 0
    pivots = []
    for column in range(variables):
        selected = next(
            (index for index in range(pivot, len(rows)) if rows[index][column]),
            None,
        )
        if selected is None:
            continue
        rows[pivot], rows[selected] = rows[selected], rows[pivot]
        scale = rows[pivot][column]
        rows[pivot] = [value / scale for value in rows[pivot]]
        for index in range(len(rows)):
            if index != pivot and rows[index][column]:
                scale = rows[index][column]
                rows[index] = [
                    value - scale * other
                    for value, other in zip(rows[index], rows[pivot])
                ]
        pivots.append(column)
        pivot += 1
    if any(not any(row[:variables]) and row[variables] for row in rows):
        return None
    if pivot != variables:
        return None
    result = [F(0)] * variables
    for index, column in enumerate(pivots):
        result[column] = rows[index][variables]
    return tuple(result)


@cache
def kernel_vertices(matrix):
    matrix = tuple(tuple(map(F, row)) for row in matrix)
    row_count = len(matrix)
    if not row_count:
        return ()
    dimension = len(matrix[0])
    vertices = set()
    for size in range(1, min(row_count, dimension + 1) + 1):
        for support in combinations(range(row_count), size):
            equations = [
                [matrix[index][coordinate] for index in support]
                for coordinate in range(dimension)
            ] + [[1] * size]
            solution = unique_solution(equations, [0] * dimension + [1], size)
            if solution is None or not all(value > 0 for value in solution):
                continue
            vertex = tuple(
                solution[support.index(index)] if index in support else F(0)
                for index in range(row_count)
            )
            assert sum(vertex) == 1 and all(value >= 0 for value in vertex)
            assert all(
                sum(vertex[index] * matrix[index][coordinate] for index in range(row_count))
                == 0
                for coordinate in range(dimension)
            )
            vertices.add(vertex)
    return tuple(sorted(vertices))


def primitive(vector):
    vector = tuple(map(F, vector))
    scale = lcm(*(value.denominator for value in vector)) if vector else 1
    integers = tuple(int(value * scale) for value in vector)
    divisor = gcd(*integers) if integers else 0
    return tuple(value // divisor for value in integers) if divisor else integers


def parameter_domain(rows, free_count: int, base_count: int):
    rows = tuple(sorted(set(tuple(map(F, row)) for row in rows)))
    free = tuple(tuple(row[:free_count]) for row in rows)
    base = tuple(tuple(row[free_count:]) for row in rows)
    vertices = kernel_vertices(free)
    forms = tuple(
        sorted(
            {
                primitive(
                    tuple(
                        sum(vertex[index] * base[index][coordinate] for index in range(len(rows)))
                        for coordinate in range(base_count)
                    )
                )
                for vertex in vertices
            }
        )
    )
    return forms, vertices


@cache
def strict_feasible(forms):
    forms = tuple(tuple(map(F, form)) for form in forms)
    return not kernel_vertices(forms)


def topology(faces):
    faces = [tuple(face) for face in faces]
    if not faces:
        return {"f_vector": [], "boundary_ranks": [], "betti_F2": []}
    levels = [
        sorted(face for face in faces if len(face) == size)
        for size in range(1, max(map(len, faces)) + 1)
    ]
    ranks = [0]
    for dimension in range(1, len(levels)):
        indices = {face: index for index, face in enumerate(levels[dimension - 1])}
        pivots = {}
        for face in levels[dimension]:
            column = sum(
                1 << indices[face[:index] + face[index + 1 :]]
                for index in range(len(face))
            )
            while column:
                pivot = column.bit_length() - 1
                if pivot not in pivots:
                    pivots[pivot] = column
                    break
                column ^= pivots[pivot]
        ranks.append(len(pivots))
    ranks.append(0)
    return {
        "f_vector": list(map(len, levels)),
        "boundary_ranks": ranks,
        "betti_F2": [
            len(levels[index]) - ranks[index] - ranks[index + 1]
            for index in range(len(levels))
        ],
    }


def clean(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {str(key): clean(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(item) for item in value]
    return value


def replay_gluing(data):
    assert data["status"] == "passed"
    totals = {"certificates": 0, "concepts": 0, "rows": 0, "comparisons": 0}

    def accept(packet):
        result = decode_certificate(packet)
        totals["certificates"] += 1
        for key in ("concepts", "rows", "comparisons"):
            totals[key] += result[key]
        return result

    cycles = data["cycle_gluings"]
    assert {cycle["cycle_length"] for cycle in cycles} == {4, 6}
    for cycle in cycles:
        length = cycle["cycle_length"]
        assert cycle["union_size"] == 8 * length + 8
        assert cycle["section_size"] == 22
        assert cycle["piece_sizes"] == [31, 8 * length - 1]
        assert len(cycle["cross_rows"]) == length - 3
        assert len(cycle["shared_rows"]) == 3
        assert cycle["epsilon"] == "1/202"
        result = accept(cycle["certificate"])
        assert result["concepts"] == cycle["union_size"]
        assert result["margin"] == "1/654480"

    degenerate = data["full_cube_and_inactive_test"]
    assert degenerate["union_size"] == 6
    assert degenerate["piece_sizes"] == [4, 4, 2]
    assert len(degenerate["cross_rows"]) == 1
    accept(degenerate["certificate"])

    zero = data["zero_coordinate_test"]
    assert zero["n"] == 0 and zero["union_size"] == 1
    accept(zero["certificate"])

    fibers = data["relative_fiber_points"]
    assert [fiber["name"] for fiber in fibers] == [
        "opposite-quadrant-1",
        "opposite-quadrant-2",
        "contractible-fiber",
    ]
    for fiber in fibers:
        result = accept(fiber["certificate"])
        assert result["concepts"] == 31
    assert fibers[0]["section_rows"] == fibers[1]["section_rows"]
    assert fibers[0]["section_rows"] != fibers[2]["section_rows"]
    assert [fiber["middle_control_ratios"] for fiber in fibers] == [
        ["2", "1/2"],
        ["1/2", "2"],
        ["2", "2"],
    ]
    assert totals["certificates"] == 7
    return totals


def replay_nerve(data):
    assert data["status"] == "passed"
    models = data["models"]
    assert len(models) == 2
    domain_certificates = 0
    all_forms = set()
    for model in models:
        assert len(model["words"]) == 31
        assert len(model["coordinates"]) == 6
        assert len(model["free"]) == 3
        assert len(model["base"]) == 7
        assert len(model["selections"]) == 4
        assert len(model["domains"]) == 15
        assert len(model["certificates"]) == 15
        protected = {tuple(pattern) for pattern in model["protected_patterns"]}
        assert protected == {tuple(entry[0]) for entry in model["base"]}
        local_missing = missing(model["words"], len(model["coordinates"]))
        coordinates = model["coordinates"]
        lifted = {
            (
                sum(((support >> local) & 1) << global_coordinate for local, global_coordinate in enumerate(coordinates)),
                sum(((bits >> local) & 1) << global_coordinate for local, global_coordinate in enumerate(coordinates)),
            )
            for support, bits in local_missing
        }
        assert protected <= lifted
        for key, certificate in model["certificates"].items():
            forms, vertices = parameter_domain(
                certificate["rows"], len(model["free"]), len(model["base"])
            )
            assert clean(forms) == certificate["parameter_forms"]
            assert clean(vertices) == certificate["normalized_kernel_vertices"]
            assert clean(forms) == model["domains"][key]
            all_forms.update(forms)
            domain_certificates += 1

    assert len(data["arrangement_forms"]) == 4
    assert {tuple(form) for form in data["arrangement_forms"]} == all_forms
    assert len(data["sign_cells"]) == 9
    for cell in data["sign_cells"]:
        point = cell["base_log2_ratios"]
        for model, expected in zip(models, cell["fibers"]):
            faces = [
                tuple(index for index in range(len(model["selections"])) if int(mask) >> index & 1)
                for mask, forms in model["domains"].items()
                if all(
                    sum(F(coefficient) * value for coefficient, value in zip(form, point)) > 0
                    for form in forms
                )
            ]
            assert clean(faces) == expected["faces"]
            assert topology(faces) == {
                key: expected[key]
                for key in ("f_vector", "boundary_ranks", "betti_F2")
            }

    gluing = data["gluing"]
    candidates = [tuple(candidate) for candidate in gluing["candidate_vertices"]]
    assert len(candidates) == 16

    def common_domain(indices):
        masks = [
            sum(1 << index for index in {candidates[vertex][piece] for vertex in indices})
            for piece in (0, 1)
        ]
        forms = tuple(
            sorted(
                {
                    tuple(form)
                    for model, mask in zip(models, masks)
                    for form in model["domains"][str(mask)]
                }
            )
        )
        return strict_feasible(forms)

    vertices = [index for index in range(16) if common_domain((index,))]
    assert vertices == gluing["vertices"] and len(vertices) == 14
    faces = []
    for mask in range(1, 1 << len(vertices)):
        face = tuple(vertex for index, vertex in enumerate(vertices) if mask >> index & 1)
        if common_domain(face):
            faces.append(face)
    assert clean(faces) == gluing["faces"]
    result = topology(faces)
    assert result == {
        key: gluing[key] for key in ("f_vector", "boundary_ranks", "betti_F2")
    }
    assert result["f_vector"] == [14, 48, 48, 12]
    assert result["betti_F2"] == [1, 0, 1, 0]
    return {
        "models": len(models),
        "domain_certificates": domain_certificates,
        "oriented_parameter_forms": len(all_forms),
        "sign_cells": len(data["sign_cells"]),
        "gluing_vertices": len(vertices),
        "gluing_faces": len(faces),
        **result,
    }


def replay():
    for path in (GLUING, NERVE, INTEGRATION):
        assert sha256(path.read_bytes()).hexdigest() == EXPECTED_SHA256[path.name]
    gluing = replay_gluing(json.loads(GLUING.read_text(encoding="utf-8")))
    nerve = replay_nerve(json.loads(NERVE.read_text(encoding="utf-8")))
    integration = json.loads(INTEGRATION.read_text(encoding="utf-8"))
    assert integration["status"] == "passed"
    assert integration["old_labels"] == 61
    assert integration["new_labels"] == 92
    assert len(integration["added_labels"]) == 31
    assert not integration["removed_labels"]
    return {
        "status": "passed",
        "source_sha256": EXPECTED_SHA256,
        "gluing_checks": gluing,
        "relative_nerve_checks": nerve,
        "integration": {
            "old_labels": integration["old_labels"],
            "new_labels": integration["new_labels"],
            "added_labels": len(integration["added_labels"]),
            "removed_labels": len(integration["removed_labels"]),
            "old_statements_preserved": integration["old_statements_preserved"],
            "old_proof_environments_preserved": integration["old_proof_environments_preserved"],
        },
        "scope": {
            "gluing": "Exact replay of all seven retained positive-margin certificates and every separator inequality.",
            "relative_nerve": "Exact elimination replay for 30 face domains, all nine wall cells, and every face of the four-cycle gluing nerve.",
            "limits": "Universal homeomorphism and homotopy statements remain symbolic database-owned proofs; the finite packets certify their named instances.",
        },
    }


def main():
    report = replay()
    OUTPUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
