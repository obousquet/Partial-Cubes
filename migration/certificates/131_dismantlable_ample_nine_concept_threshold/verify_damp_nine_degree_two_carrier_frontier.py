#!/usr/bin/env python3
"""Reduce the nine-facet frontier using degree-two carrier geometry.

This verifier keeps passive coordinates symbolically through the inclusion
preorder of their support parts.  It proves that only profiles 5.1, 8.1,
and 9.1 survive the degree-two isolated-support and cross-carrier laws.
"""

from __future__ import annotations

import hashlib
import json

import enumerate_damp_nine_colored_embeddings as colored_frontier


EXPECTED_SURVIVORS = ("5.1", "8.1", "9.1")


def profile_witness(row: dict[str, object]) -> tuple[dict[int, int] | None, int]:
    vertices = tuple(sorted(int(word) for word in row["embedding"]))
    vertex_set = set(vertices)
    dimension = int(row["dimension"])
    full = (1 << dimension) - 1
    directions = {
        vertex: tuple(
            coordinate
            for coordinate in range(dimension)
            if vertex ^ (1 << coordinate) in vertex_set
        )
        for vertex in vertices
    }
    markers = tuple(vertex for vertex in vertices if len(directions[vertex]) == 2)
    domains = {
        vertex: tuple(
            support
            for support in range(1 << dimension)
            if sum(
                (support >> coordinate) & 1
                for coordinate in directions[vertex]
            )
            == 1
        )
        for vertex in markers
    }

    def directed_cross(
        lower_marker: int,
        lower_support: int,
        upper_marker: int,
        upper_support: int,
    ) -> tuple[bool, bool]:
        overlap = lower_support & (full ^ upper_support)
        if (lower_marker ^ upper_marker) & overlap:
            return True, False
        free = upper_support & (full ^ lower_support)
        if not free:
            return False, False
        base = (
            (lower_marker & lower_support)
            | (upper_marker & (full ^ upper_support))
        )
        submask = free
        while (base ^ submask) in vertex_set:
            if submask == 0:
                # The full cross cube has no passive direction, so the
                # upper passive support must be contained in the lower one.
                return True, True
            submask = (submask - 1) & free
        return False, False

    def consistent(assignment: dict[int, int]) -> bool:
        assigned = tuple(assignment)
        passive_subset = {(vertex, vertex) for vertex in assigned}
        for first_index, first in enumerate(assigned):
            for second in assigned[first_index + 1 :]:
                first_support = assignment[first]
                second_support = assignment[second]
                if (first ^ second).bit_count() == 1:
                    coordinate = (first ^ second).bit_length() - 1
                    if (
                        (first_support >> coordinate) & 1
                    ) == ((second_support >> coordinate) & 1):
                        return False

                valid, inclusion = directed_cross(
                    first, first_support, second, second_support
                )
                if not valid:
                    return False
                if inclusion:
                    passive_subset.add((second, first))

                valid, inclusion = directed_cross(
                    second, second_support, first, first_support
                )
                if not valid:
                    return False
                if inclusion:
                    passive_subset.add((first, second))

        changed = True
        while changed:
            changed = False
            for first in assigned:
                for middle in assigned:
                    if (first, middle) not in passive_subset:
                        continue
                    for second in assigned:
                        if (
                            (middle, second) in passive_subset
                            and (first, second) not in passive_subset
                        ):
                            passive_subset.add((first, second))
                            changed = True

        # Distinct isolated supports are incomparable.  An active inclusion
        # together with a forced passive inclusion would violate this.
        for first_index, first in enumerate(assigned):
            for second in assigned[first_index + 1 :]:
                first_support = assignment[first]
                second_support = assignment[second]
                if (
                    not (first_support & (full ^ second_support))
                    and (first, second) in passive_subset
                ):
                    return False
                if (
                    not (second_support & (full ^ first_support))
                    and (second, first) in passive_subset
                ):
                    return False
        return True

    nodes = 0

    def search(assignment: dict[int, int]) -> dict[int, int] | None:
        nonlocal nodes
        nodes += 1
        if len(assignment) == len(markers):
            return dict(assignment)

        best_marker = None
        best_domain = None
        for marker in markers:
            if marker in assignment:
                continue
            feasible = []
            for support in domains[marker]:
                assignment[marker] = support
                if consistent(assignment):
                    feasible.append(support)
                del assignment[marker]
            if not feasible:
                return None
            if best_domain is None or len(feasible) < len(best_domain):
                best_marker = marker
                best_domain = feasible

        assert best_marker is not None and best_domain is not None
        for support in best_domain:
            assignment[best_marker] = support
            witness = search(assignment)
            if witness is not None:
                return witness
            del assignment[best_marker]
        return None

    return search({}), nodes


def verify() -> dict[str, object]:
    frontier = colored_frontier.enumerate_colored_types()
    rows = []
    for graph in frontier["graphs"]:
        for profile in graph["colored_types"]:
            witness, nodes = profile_witness(profile)
            rows.append(
                {
                    "profile": profile["colored_type"],
                    "dimension": profile["dimension"],
                    "degree_two_carrier_status": (
                        "SURVIVES" if witness is not None else "EXCLUDED"
                    ),
                    "search_nodes": nodes,
                    "witness_active_supports": (
                        {str(key): value for key, value in sorted(witness.items())}
                        if witness is not None
                        else None
                    ),
                }
            )

    survivors = tuple(
        row["profile"]
        for row in rows
        if row["degree_two_carrier_status"] == "SURVIVES"
    )
    if survivors != EXPECTED_SURVIVORS:
        raise AssertionError((survivors, EXPECTED_SURVIVORS))

    record: dict[str, object] = {
        "status": "VERIFIED",
        "graph_type_count": frontier["graph_type_count"],
        "coordinate_coloured_profile_count": frontier[
            "colored_embedding_type_count"
        ],
        "excluded_profile_count": len(rows) - len(survivors),
        "surviving_profiles": list(survivors),
        "profiles": rows,
    }
    canonical = json.dumps(record, sort_keys=True, separators=(",", ":"))
    record["mathematical_digest"] = hashlib.sha256(canonical.encode()).hexdigest()
    return record


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
