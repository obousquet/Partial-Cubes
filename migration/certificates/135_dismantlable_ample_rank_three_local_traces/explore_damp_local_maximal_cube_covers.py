#!/usr/bin/env python3
"""Enumerate irredundant local trace covers of maximal cubes of rank <= 3.

If ``C`` is a maximal cube of a cornerless cube family, the nonempty
proper traces ``C intersect C'`` over the other maximal cubes cover every
vertex of ``C``.  When the VC dimension is at most three, ``C`` is one of
``Q_0,Q_1,Q_2,Q_3``.  This script enumerates the irredundant covers of
``Q_r`` by proper subcubes, modulo the full automorphism group of ``Q_r``.

The output is a conjecture diagnostic and a finite local-pattern catalogue;
it is not used as a proof certificate for the dimension-free trace-cover
lemma.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path


OUTPUT = Path(__file__).with_name("damp_local_maximal_cube_covers.json")


def proper_subcubes(rank: int) -> tuple[int, ...]:
    """Return the vertex masks of all nonempty proper subcubes of Q_rank."""

    answer = set()
    for fixed_mask in range(1, 1 << rank):
        value_mask = fixed_mask
        while True:
            vertices = 0
            for vertex in range(1 << rank):
                if vertex & fixed_mask == value_mask:
                    vertices |= 1 << vertex
            answer.add(vertices)
            if value_mask == 0:
                break
            value_mask = (value_mask - 1) & fixed_mask
    return tuple(sorted(answer, key=lambda mask: (-mask.bit_count(), mask)))


def cube_automorphisms(rank: int) -> tuple[tuple[int, ...], ...]:
    """Return all vertex permutations induced by coordinate permutations/flips."""

    maps = []
    for permutation in itertools.permutations(range(rank)):
        for flip in range(1 << rank):
            image = []
            for vertex in range(1 << rank):
                target = flip
                for target_coordinate, source_coordinate in enumerate(permutation):
                    if vertex & (1 << source_coordinate):
                        target ^= 1 << target_coordinate
                image.append(target)
            maps.append(tuple(image))
    return tuple(maps)


def transform_mask(mask: int, vertex_map: tuple[int, ...]) -> int:
    answer = 0
    for vertex, image in enumerate(vertex_map):
        if mask & (1 << vertex):
            answer |= 1 << image
    return answer


def canonical_cover(
    cover: frozenset[int], automorphisms: tuple[tuple[int, ...], ...]
) -> tuple[int, ...]:
    return min(
        tuple(sorted(transform_mask(mask, automorphism) for mask in cover))
        for automorphism in automorphisms
    )


def irredundant_covers(rank: int) -> set[frozenset[int]]:
    """Enumerate covers in which every selected subcube has a private vertex."""

    universe = (1 << (1 << rank)) - 1
    candidates = proper_subcubes(rank)
    through_vertex = tuple(
        tuple(index for index, mask in enumerate(candidates) if mask & (1 << vertex))
        for vertex in range(1 << rank)
    )
    visited: set[frozenset[int]] = set()
    answers: set[frozenset[int]] = set()

    def search(chosen_indices: frozenset[int], covered: int) -> None:
        if chosen_indices in visited:
            return
        visited.add(chosen_indices)
        if covered == universe:
            cover = frozenset(candidates[index] for index in chosen_indices)
            for index in chosen_indices:
                other_union = 0
                for other in chosen_indices:
                    if other != index:
                        other_union |= candidates[other]
                if candidates[index] & ~other_union == 0:
                    return
            answers.add(cover)
            return

        uncovered = universe & ~covered
        vertex_bit = uncovered & -uncovered
        vertex = vertex_bit.bit_length() - 1
        for index in through_vertex[vertex]:
            if index in chosen_indices:
                continue
            candidate = candidates[index]
            # An irredundant cover cannot contain comparable subcubes.
            if any(
                candidate & candidates[other] in (candidate, candidates[other])
                for other in chosen_indices
            ):
                continue
            search(chosen_indices | {index}, covered | candidate)

    search(frozenset(), 0)
    return answers


def trace_word(mask: int, rank: int) -> str:
    symbols = []
    vertices = [vertex for vertex in range(1 << rank) if mask & (1 << vertex)]
    for coordinate in range(rank):
        values = {(vertex >> coordinate) & 1 for vertex in vertices}
        symbols.append("*" if len(values) == 2 else str(next(iter(values))))
    return "".join(symbols)


def is_chordal(adjacency: tuple[frozenset[int], ...]) -> bool:
    count = len(adjacency)
    unnumbered = set(range(count))
    weights = [0] * count
    order = []
    while unnumbered:
        vertex = min(unnumbered, key=lambda item: (-weights[item], item))
        unnumbered.remove(vertex)
        order.append(vertex)
        for neighbor in adjacency[vertex] & unnumbered:
            weights[neighbor] += 1
    position = {vertex: index for index, vertex in enumerate(order)}
    for vertex in order:
        earlier = [
            neighbor
            for neighbor in adjacency[vertex]
            if position[neighbor] < position[vertex]
        ]
        if len(earlier) < 2:
            continue
        parent = max(earlier, key=position.__getitem__)
        if any(
            neighbor != parent and neighbor not in adjacency[parent]
            for neighbor in earlier
        ):
            return False
    return True


def cover_profile(cover: tuple[int, ...], rank: int) -> dict[str, object]:
    adjacency = [set() for _ in cover]
    for left, right in itertools.combinations(range(len(cover)), 2):
        if cover[left] & cover[right]:
            adjacency[left].add(right)
            adjacency[right].add(left)
    adjacency_tuple = tuple(frozenset(row) for row in adjacency)
    edge_count = sum(map(len, adjacency_tuple)) // 2
    dimensions = sorted(
        (mask.bit_count().bit_length() - 1 for mask in cover), reverse=True
    )
    common_intersection = cover[0]
    for mask in cover[1:]:
        common_intersection &= mask
    return {
        "cover_size": len(cover),
        "trace_dimensions": dimensions,
        "trace_words": sorted(trace_word(mask, rank) for mask in cover),
        "intersection_edge_count": edge_count,
        "disjoint_pair_count": len(cover) * (len(cover) - 1) // 2 - edge_count,
        "intersection_degree_sequence": sorted(
            (len(row) for row in adjacency_tuple), reverse=True
        ),
        "intersection_graph_chordal": is_chordal(adjacency_tuple),
        "common_intersection_size": common_intersection.bit_count(),
    }


def rank_catalogue(rank: int) -> dict[str, object]:
    raw = irredundant_covers(rank)
    automorphisms = cube_automorphisms(rank)
    orbit_counts: Counter[tuple[int, ...]] = Counter(
        canonical_cover(cover, automorphisms) for cover in raw
    )
    representatives = sorted(
        orbit_counts,
        key=lambda cover: (
            len(cover),
            tuple(-mask.bit_count() for mask in cover),
            cover,
        ),
    )
    profiles = [cover_profile(cover, rank) for cover in representatives]
    by_size = Counter(profile["cover_size"] for profile in profiles)
    by_dimension_signature = Counter(
        tuple(profile["trace_dimensions"]) for profile in profiles
    )
    by_disjoint_pairs = Counter(profile["disjoint_pair_count"] for profile in profiles)
    return {
        "rank": rank,
        "proper_subcube_count": len(proper_subcubes(rank)),
        "raw_irredundant_cover_count": len(raw),
        "automorphism_orbit_count": len(representatives),
        "minimum_cover_size": min(by_size),
        "maximum_cover_size": max(by_size),
        "minimum_disjoint_pair_count": min(by_disjoint_pairs),
        "all_common_intersections_empty": all(
            profile["common_intersection_size"] == 0 for profile in profiles
        ),
        "orbit_count_by_cover_size": {
            str(key): value for key, value in sorted(by_size.items())
        },
        "orbit_count_by_trace_dimensions": {
            ",".join(map(str, signature)): count
            for signature, count in sorted(by_dimension_signature.items())
        },
        "orbit_count_by_disjoint_pair_count": {
            str(key): value for key, value in sorted(by_disjoint_pairs.items())
        },
        "orbits": [
            {
                "orbit_size": orbit_counts[representative],
                **profile,
            }
            for representative, profile in zip(representatives, profiles)
        ],
    }


def main() -> None:
    catalogues = [rank_catalogue(rank) for rank in range(1, 4)]
    payload = {
        "schema": "damp-local-maximal-cube-covers-v1",
        "scope": (
            "Irredundant covers of Q_r by nonempty proper subcubes, modulo "
            "coordinate permutations and flips, for 1 <= r <= 3."
        ),
        "catalogues": catalogues,
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    payload["sha256_without_digest"] = hashlib.sha256(canonical.encode()).hexdigest()
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    for catalogue in catalogues:
        print(
            "rank={rank} raw={raw_irredundant_cover_count} "
            "orbits={automorphism_orbit_count} sizes={minimum_cover_size}.."
            "{maximum_cover_size} min_disjoint_pairs={minimum_disjoint_pair_count}".format(
                **catalogue
            )
        )
    print(f"wrote {OUTPUT}")


if __name__ == "__main__":
    main()
