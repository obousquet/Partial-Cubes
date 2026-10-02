#!/usr/bin/env python3
"""Audit conditioned active states for the seven eight-facet graph types.

For every graph whose minimum cube embedding has at most four coordinates,
the script enumerates all conditioned lower classes, retains those with
nonnegative active fine numerator, reduces them by affine cube symmetries
preserving the relative set, and tests interval partitions and initial
facets.  The five-coordinate barbell is handled separately by the symbolic
boundary solver.
"""

from __future__ import annotations

from collections import Counter
from itertools import permutations

from enumerate_eight_deletion_graphs import (
    N,
    canonical_matrix,
    degrees,
    edges_from_mask,
    induced_cube_embedding,
    matrix_graph_mask,
    relevant_matrices,
)


def graph_types():
    abstract = {}
    for left in range(1, N // 2 + 1):
        right = N - left
        for matrix in relevant_matrices(left, right):
            canonical = canonical_matrix(matrix, left, right)
            key = (left, right, canonical)
            abstract.setdefault(
                key, matrix_graph_mask(canonical, left, right)
            )
    result = []
    for key, mask in abstract.items():
        embedding = induced_cube_embedding(mask)
        if embedding is not None:
            result.append((key, mask, embedding))
    return sorted(
        result,
        key=lambda item: (
            len(edges_from_mask(item[1])),
            sorted(degrees(item[1]), reverse=True),
            item[0],
        ),
    )


def projection_blocks(
    relative: tuple[int, ...], support: int
) -> tuple[tuple[int, ...], ...]:
    blocks = {}
    for concept in relative:
        blocks.setdefault(concept & support, []).append(concept)
    return tuple(tuple(block) for block in blocks.values())


def eligible_blocks(
    relative: tuple[int, ...],
    lower: frozenset[int],
    dimension: int,
) -> dict[int, tuple[tuple[int, ...], ...]]:
    result = {}
    for support in range(1 << dimension):
        traces = {concept & support for concept in lower}
        result[support] = tuple(
            block
            for block in projection_blocks(relative, support)
            if (block[0] & support) not in traces
        )
    return result


def submasks(mask: int):
    submask = mask
    while True:
        yield submask
        if submask == 0:
            return
        submask = (submask - 1) & mask


def active_fine_h(
    blocks: dict[int, tuple[tuple[int, ...], ...]], dimension: int
) -> tuple[int, ...]:
    return tuple(
        sum(
            (-1) ** (support.bit_count() - submask.bit_count())
            * len(blocks[submask])
            for submask in submasks(support)
        )
        for support in range(1 << dimension)
    )


def permute_bits(concept: int, permutation: tuple[int, ...]) -> int:
    return sum(
        ((concept >> old) & 1) << new
        for new, old in enumerate(permutation)
    )


def stabilizer(
    relative: tuple[int, ...], dimension: int
) -> tuple[tuple[int, tuple[int, ...]], ...]:
    relative_set = set(relative)
    result = []
    for permutation in permutations(range(dimension)):
        permuted = tuple(permute_bits(concept, permutation) for concept in relative)
        for translation in range(1 << dimension):
            if {concept ^ translation for concept in permuted} == relative_set:
                result.append((translation, permutation))
    return tuple(result)


def canonical_lower(
    lower: frozenset[int],
    automorphisms: tuple[tuple[int, tuple[int, ...]], ...],
) -> tuple[int, ...]:
    return min(
        tuple(
            sorted(
                permute_bits(concept, permutation) ^ translation
                for concept in lower
            )
        )
        for translation, permutation in automorphisms
    )


def interval_partitions(
    relative: tuple[int, ...],
    blocks: dict[int, tuple[tuple[int, ...], ...]],
    active_h: tuple[int, ...],
) -> tuple[tuple[int, ...], ...]:
    index = {concept: position for position, concept in enumerate(relative)}
    lower_supports = tuple(
        support
        for support, multiplicity in enumerate(active_h)
        for _ in range(multiplicity)
    )
    assert len(lower_supports) == len(relative)
    result = []
    for assignment in set(permutations(lower_supports)):
        valid = True
        for support, support_blocks in blocks.items():
            selected = {
                position
                for position, root_support in enumerate(assignment)
                if root_support & ~support == 0
            }
            if any(
                len(selected & {index[concept] for concept in block}) != 1
                for block in support_blocks
            ):
                valid = False
                break
        if valid:
            result.append(assignment)
    return tuple(sorted(result))


def initial_count(
    relative: tuple[int, ...],
    lower: frozenset[int],
    assignment: tuple[int, ...],
    dimension: int,
) -> int:
    return sum(
        assignment[index]
        & ~sum(
            1 << coordinate
            for coordinate in range(dimension)
            if concept ^ (1 << coordinate) in lower
        )
        == 0
        for index, concept in enumerate(relative)
    )


def active_teaching_roots(
    relative: tuple[int, ...],
    lower: frozenset[int],
    dimension: int,
) -> tuple[int, ...]:
    result = []
    for concept in relative:
        teaching_coordinates = {
            coordinate
            for coordinate in range(dimension)
            if concept ^ (1 << coordinate) in lower
        }
        if all(
            any(((concept ^ lower_concept) >> coordinate) & 1
                for coordinate in teaching_coordinates)
            for lower_concept in lower
        ):
            result.append(concept)
    return tuple(result)


def main() -> int:
    failures = 0
    for number, (_key, mask, embedding) in enumerate(graph_types(), start=1):
        dimension = max(embedding).bit_length()
        relative = tuple(embedding)
        if dimension > 4:
            print(
                f"type {number}: dimension={dimension}; "
                "deferred to symbolic five-coordinate audit"
            )
            continue
        holes = tuple(
            concept for concept in range(1 << dimension)
            if concept not in relative
        )
        nonnegative = []
        data = {}
        for choice in range(1 << len(holes)):
            lower = frozenset(
                concept
                for index, concept in enumerate(holes)
                if (choice >> index) & 1
            )
            blocks = eligible_blocks(relative, lower, dimension)
            fine_h = active_fine_h(blocks, dimension)
            if min(fine_h) >= 0:
                nonnegative.append(lower)
                data[lower] = (blocks, fine_h)
        automorphisms = stabilizer(relative, dimension)
        orbit_sizes = Counter(
            canonical_lower(lower, automorphisms) for lower in nonnegative
        )
        empty_partition_orbits = 0
        no_initial_orbits = 0
        rootless_orbits = 0
        profiles = Counter()
        for representative in sorted(
            orbit_sizes, key=lambda item: (len(item), item)
        ):
            lower = frozenset(representative)
            blocks, fine_h = data[lower]
            partitions = interval_partitions(relative, blocks, fine_h)
            teaching = active_teaching_roots(relative, lower, dimension)
            minimum_initial = (
                min(
                    initial_count(
                        relative, lower, assignment, dimension
                    )
                    for assignment in partitions
                )
                if partitions else 0
            )
            empty_partition_orbits += not partitions
            no_initial_orbits += bool(partitions) and minimum_initial == 0
            rootless_orbits += not teaching
            profiles[
                (
                    len(lower),
                    len(partitions),
                    minimum_initial,
                    len(teaching),
                )
            ] += 1
        failures += empty_partition_orbits + no_initial_orbits
        print(
            f"type {number}: dimension={dimension} holes={len(holes)} "
            f"automorphisms={len(automorphisms)} "
            f"nonnegative_states={len(nonnegative)} "
            f"orbits={len(orbit_sizes)} rootless_orbits={rootless_orbits} "
            f"empty_partition_orbits={empty_partition_orbits} "
            f"no_initial_orbits={no_initial_orbits}"
        )
        print(f"  profiles={dict(sorted(profiles.items()))}")
    print({"failures": failures})
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
