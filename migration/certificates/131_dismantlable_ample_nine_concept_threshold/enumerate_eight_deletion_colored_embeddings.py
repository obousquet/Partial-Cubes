#!/usr/bin/env python3
"""Enumerate color-respecting embeddings of the seven D8 graph types.

The uncolored induced-cube graph does not determine the coordinate labels on
its edges.  A coordinate of an embedding is a normalized vertex cut whose
crossing graph edges form a matching.  A color-respecting cube embedding is
therefore exactly a set of such cuts which

* partitions the graph edges, and
* separates every graph nonedge at least twice.

This script enumerates every such cut system and quotients by graph
automorphisms, cut complementation, and coordinate permutation.  Constant
coordinates are omitted.  Thus its output classifies active cube embeddings,
not merely abstract graph types.
"""

from __future__ import annotations

from itertools import permutations

from enumerate_eight_deletion_graphs import (
    N,
    PAIRS,
    degrees,
    edges_from_mask,
)
from diagnose_eight_facet_active_states import graph_types


FULL_VERTEX_MASK = (1 << N) - 1


def graph_automorphisms(mask: int) -> tuple[tuple[int, ...], ...]:
    edge_set = set(edges_from_mask(mask))
    degree_sequence = degrees(mask)
    result = []
    for permutation in permutations(range(N)):
        if any(
            degree_sequence[vertex]
            != degree_sequence[permutation[vertex]]
            for vertex in range(N)
        ):
            continue
        image = {
            tuple(sorted((permutation[first], permutation[second])))
            for first, second in edge_set
        }
        if image == edge_set:
            result.append(permutation)
    return tuple(result)


def normalized_image(side: int, permutation: tuple[int, ...]) -> int:
    image = sum(
        1 << permutation[vertex]
        for vertex in range(N)
        if (side >> vertex) & 1
    )
    if image & 1:
        image ^= FULL_VERTEX_MASK
    return image


def cut_system_canonical(
    sides: tuple[int, ...],
    automorphisms: tuple[tuple[int, ...], ...],
) -> tuple[int, ...]:
    return min(
        tuple(
            sorted(
                normalized_image(side, permutation)
                for side in sides
            )
        )
        for permutation in automorphisms
    )


def admissible_cuts(mask: int):
    edges = edges_from_mask(mask)
    cuts = []
    for side_without_zero in range(1, 1 << (N - 1)):
        side = side_without_zero << 1
        cut_edges = 0
        used_vertices = 0
        valid = True
        for edge_number, (first, second) in enumerate(edges):
            if not (((side >> first) ^ (side >> second)) & 1):
                continue
            endpoints = (1 << first) | (1 << second)
            if used_vertices & endpoints:
                valid = False
                break
            used_vertices |= endpoints
            cut_edges |= 1 << edge_number
        if valid and cut_edges:
            cuts.append((cut_edges, side))
    return tuple(cuts)


def separates_nonedges(
    sides: tuple[int, ...], mask: int
) -> bool:
    edge_set = set(edges_from_mask(mask))
    return all(
        sum(
            ((side >> first) ^ (side >> second)) & 1
            for side in sides
        )
        >= 2
        for first, second in PAIRS
        if (first, second) not in edge_set
    )


def embedding_from_sides(sides: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        sum(
            1 << coordinate
            for coordinate, side in enumerate(sides)
            if (side >> vertex) & 1
        )
        for vertex in range(N)
    )


def colored_embedding_types(mask: int):
    edges = edges_from_mask(mask)
    all_edges = (1 << len(edges)) - 1
    cuts = admissible_cuts(mask)
    cuts_by_edge = [[] for _ in edges]
    for cut_number, (cut_edges, _side) in enumerate(cuts):
        for edge_number in range(len(edges)):
            if (cut_edges >> edge_number) & 1:
                cuts_by_edge[edge_number].append(cut_number)

    raw_systems = []

    def exact_cover(covered: int, chosen: tuple[int, ...]) -> None:
        if covered == all_edges:
            sides = tuple(sorted(cuts[index][1] for index in chosen))
            if separates_nonedges(sides, mask):
                raw_systems.append(sides)
            return
        missing = all_edges ^ covered
        edge_number = (missing & -missing).bit_length() - 1
        for cut_number in cuts_by_edge[edge_number]:
            cut_edges = cuts[cut_number][0]
            if cut_edges & covered:
                continue
            exact_cover(covered | cut_edges, chosen + (cut_number,))

    exact_cover(0, ())
    automorphisms = graph_automorphisms(mask)
    representatives = {}
    for sides in raw_systems:
        canonical = cut_system_canonical(sides, automorphisms)
        representatives.setdefault(canonical, 0)
        representatives[canonical] += 1
    return automorphisms, raw_systems, representatives


def bitset_text(value: int, width: int) -> str:
    return format(value, f"0{width}b")


def main() -> int:
    total_types = 0
    for type_number, (_key, mask, _minimum_embedding) in enumerate(
        graph_types(), start=1
    ):
        automorphisms, raw_systems, representatives = (
            colored_embedding_types(mask)
        )
        total_types += len(representatives)
        dimension_counts = {}
        for sides in representatives:
            dimension_counts[len(sides)] = (
                dimension_counts.get(len(sides), 0) + 1
            )
        print(
            f"type {type_number}: graph_automorphisms={len(automorphisms)} "
            f"raw_cut_systems={len(raw_systems)} "
            f"colored_types={len(representatives)} "
            f"dimensions={dimension_counts}"
        )
        for number, (sides, orbit_size) in enumerate(
            sorted(
                representatives.items(),
                key=lambda item: (len(item[0]), item[0]),
            ),
            start=1,
        ):
            embedding = embedding_from_sides(sides)
            cut_sizes = tuple(
                sum(
                    ((side >> first) ^ (side >> second)) & 1
                    for first, second in edges_from_mask(mask)
                )
                for side in sides
            )
            print(
                f"  {type_number}.{number}: dimension={len(sides)} "
                f"cut_sizes={cut_sizes} orbit={orbit_size} "
                f"embedding="
                + "{"
                + ",".join(
                    bitset_text(vertex, len(sides))
                    for vertex in embedding
                )
                + "}"
            )
    print({"colored_embedding_types": total_types})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
