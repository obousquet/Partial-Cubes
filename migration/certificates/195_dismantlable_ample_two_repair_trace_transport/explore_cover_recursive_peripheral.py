"""Explore the cover-recursively peripheral partial cubes.

Mathematical model
------------------
A finite partial cube is represented by an irredundant well-graded family
``F`` contained in ``{0,1}^d``.  For a coordinate ``i``, let ``F_i^0`` and
``F_i^1`` be its two projected sections.  The coordinate is peripheral
exactly when these sections are nested.  Contracting the coordinate gives
their union; the smaller section is the site of the inverse peripheral
expansion.

The class ``R_per`` is tested by the exact recursion

    K_1 is in R_per;
    F is in R_per iff some peripheral coordinate has both its contraction
    and its expansion site in R_per.

For comparison, tree-likeness only recursively tests the contraction and
does not require the expansion site to be tree-like.  The script exhaustively
enumerates cube-labelled partial cubes through a requested isometric
dimension, quotients them by coordinate permutations and flips, checks the
immediate pc-minors, and tests Cartesian products in the enumerated range.

This is a finite diagnostic, not a proof of the general closure theorems.
Run with

    uv run python explore_cover_recursive_peripheral.py --max-dimension 4
"""

from __future__ import annotations

import argparse
import itertools
from collections import deque
from functools import lru_cache
from typing import Iterable


def vertices(family: int, dimension: int) -> tuple[int, ...]:
    """Return the represented cube vertices in increasing order."""

    return tuple(v for v in range(1 << dimension) if family >> v & 1)


def family_from_vertices(members: Iterable[int]) -> int:
    return sum(1 << vertex for vertex in members)


def project_vertex(vertex: int, coordinate: int) -> int:
    low_mask = (1 << coordinate) - 1
    return (vertex & low_mask) | ((vertex >> (coordinate + 1)) << coordinate)


def section(family: int, dimension: int, coordinate: int, value: int) -> int:
    return family_from_vertices(
        [
            project_vertex(vertex, coordinate)
            for vertex in vertices(family, dimension)
            if (vertex >> coordinate) & 1 == value
        ]
    )


def irredundant_projection(family: int, dimension: int) -> tuple[int, int]:
    """Delete all coordinates that are constant on the family."""

    members = vertices(family, dimension)
    if not members:
        raise ValueError("the empty family does not represent a partial cube")
    varying = [
        coordinate
        for coordinate in range(dimension)
        if len({(vertex >> coordinate) & 1 for vertex in members}) == 2
    ]
    projected = family_from_vertices(
        [
            sum(
                ((vertex >> old_coordinate) & 1) << new_coordinate
                for new_coordinate, old_coordinate in enumerate(varying)
            )
            for vertex in members
        ]
    )
    return projected, len(varying)


def transform_vertex(
    vertex: int, permutation: tuple[int, ...], flip: int
) -> int:
    transformed = 0
    for new_coordinate, old_coordinate in enumerate(permutation):
        bit = ((vertex >> old_coordinate) & 1) ^ ((flip >> old_coordinate) & 1)
        transformed |= bit << new_coordinate
    return transformed


@lru_cache(maxsize=None)
def canonical_family(family: int, dimension: int) -> tuple[int, int]:
    family, dimension = irredundant_projection(family, dimension)
    if dimension == 0:
        return 1, 0
    members = vertices(family, dimension)
    best: int | None = None
    for permutation in itertools.permutations(range(dimension)):
        for flip in range(1 << dimension):
            transformed = family_from_vertices(
                [transform_vertex(vertex, permutation, flip) for vertex in members]
            )
            if best is None or transformed < best:
                best = transformed
    assert best is not None
    return best, dimension


@lru_cache(maxsize=None)
def is_irredundant(family: int, dimension: int) -> bool:
    members = vertices(family, dimension)
    return all(
        any(not ((vertex >> coordinate) & 1) for vertex in members)
        and any((vertex >> coordinate) & 1 for vertex in members)
        for coordinate in range(dimension)
    )


@lru_cache(maxsize=None)
def is_partial_cube_family(family: int, dimension: int) -> bool:
    """Check that the one-inclusion graph realizes Hamming distance."""

    members = vertices(family, dimension)
    if not members:
        return False
    member_set = set(members)
    for source in members:
        distances = {source: 0}
        pending = deque([source])
        while pending:
            vertex = pending.popleft()
            for coordinate in range(dimension):
                neighbor = vertex ^ (1 << coordinate)
                if neighbor in member_set and neighbor not in distances:
                    distances[neighbor] = distances[vertex] + 1
                    pending.append(neighbor)
        if len(distances) != len(members):
            return False
        if any(
            distances[target] != (source ^ target).bit_count()
            for target in members
        ):
            return False
    return True


def nested(left: int, right: int) -> bool:
    return left & ~right == 0


@lru_cache(maxsize=None)
def peripheral_decompositions(
    family: int, dimension: int
) -> tuple[tuple[int, int, int, int], ...]:
    """Return ``(coordinate, base, site, lower_dimension)`` decompositions."""

    family, dimension = irredundant_projection(family, dimension)
    decompositions: list[tuple[int, int, int, int]] = []
    for coordinate in range(dimension):
        zero = section(family, dimension, coordinate, 0)
        one = section(family, dimension, coordinate, 1)
        if nested(zero, one):
            decompositions.append((coordinate, one, zero, dimension - 1))
        elif nested(one, zero):
            decompositions.append((coordinate, zero, one, dimension - 1))
    return tuple(decompositions)


@lru_cache(maxsize=None)
def is_cover_recursive_peripheral(family: int, dimension: int) -> bool:
    family, dimension = irredundant_projection(family, dimension)
    if family.bit_count() == 1:
        return True
    for _, base, site, lower_dimension in peripheral_decompositions(
        family, dimension
    ):
        assert is_partial_cube_family(base, lower_dimension)
        assert is_partial_cube_family(site, lower_dimension)
        if is_cover_recursive_peripheral(
            base, lower_dimension
        ) and is_cover_recursive_peripheral(site, lower_dimension):
            return True
    return False


@lru_cache(maxsize=None)
def is_tree_like(family: int, dimension: int) -> bool:
    family, dimension = irredundant_projection(family, dimension)
    if family.bit_count() == 1:
        return True
    return any(
        is_tree_like(base, lower_dimension)
        for _, base, _, lower_dimension in peripheral_decompositions(
            family, dimension
        )
    )


def immediate_pc_minors(family: int, dimension: int) -> set[tuple[int, int]]:
    """Return elementary halfspace restrictions and Theta-contractions."""

    minors: set[tuple[int, int]] = set()
    for coordinate in range(dimension):
        zero = section(family, dimension, coordinate, 0)
        one = section(family, dimension, coordinate, 1)
        for minor in (zero, one, zero | one):
            minors.add(canonical_family(minor, dimension - 1))
    minors.discard(canonical_family(family, dimension))
    return minors


@lru_cache(maxsize=None)
def contains_pc_minor(
    family: int, dimension: int, target_family: int, target_dimension: int
) -> bool:
    family, dimension = canonical_family(family, dimension)
    target_family, target_dimension = canonical_family(
        target_family, target_dimension
    )
    if (family, dimension) == (target_family, target_dimension):
        return True
    if dimension <= target_dimension:
        return False
    return any(
        contains_pc_minor(
            minor, minor_dimension, target_family, target_dimension
        )
        for minor, minor_dimension in immediate_pc_minors(family, dimension)
    )


def cartesian_product(
    left: int, left_dimension: int, right: int, right_dimension: int
) -> tuple[int, int]:
    product = family_from_vertices(
        [
            left_vertex | (right_vertex << left_dimension)
            for left_vertex in vertices(left, left_dimension)
            for right_vertex in vertices(right, right_dimension)
        ]
    )
    return product, left_dimension + right_dimension


@lru_cache(maxsize=None)
def is_gated_subfamily(subfamily: int, family: int, dimension: int) -> bool:
    """Check gatedness inside a represented partial cube."""

    submembers = vertices(subfamily, dimension)
    if not submembers or subfamily & ~family:
        return False
    for vertex in vertices(family, dimension):
        gates = [
            candidate
            for candidate in submembers
            if all(
                (vertex ^ target).bit_count()
                == (vertex ^ candidate).bit_count()
                + (candidate ^ target).bit_count()
                for target in submembers
            )
        ]
        if len(gates) != 1:
            return False
    return True


def antipodally_punctured_cube(dimension: int) -> int:
    """Return ``Q_d - {0^d,1^d}``, for ``d >= 3``."""

    top = (1 << dimension) - 1
    return family_from_vertices(range(1, top))


def singly_punctured_cube(dimension: int) -> int:
    """Return ``Q_d - {1^d}``."""

    top = (1 << dimension) - 1
    return family_from_vertices(range(top))


def peripheral_cube_envelope(site: int, dimension: int) -> tuple[int, int]:
    """Peripherally expand ``Q_d`` along the represented site."""

    cube = family_from_vertices(range(1 << dimension))
    envelope = cube | family_from_vertices(
        [vertex | (1 << dimension) for vertex in vertices(site, dimension)]
    )
    return envelope, dimension + 1


def enumerate_partial_cubes(max_dimension: int) -> dict[int, list[int]]:
    """Enumerate cube-symmetry classes with irredundant dimension at most four."""

    catalogue: dict[int, list[int]] = {0: [1]}
    for dimension in range(1, max_dimension + 1):
        canonical: set[int] = set()
        # Every cube-symmetry orbit has a representative containing vertex zero.
        for family in range(1, 1 << (1 << dimension), 2):
            if not is_irredundant(family, dimension):
                continue
            if not is_partial_cube_family(family, dimension):
                continue
            canonical.add(canonical_family(family, dimension)[0])
        catalogue[dimension] = sorted(canonical)
    return catalogue


def bitstrings(family: int, dimension: int) -> str:
    return "{" + ",".join(
        format(vertex, f"0{dimension}b") for vertex in vertices(family, dimension)
    ) + "}"


def run(max_dimension: int) -> None:
    if max_dimension > 4:
        raise ValueError("exhaustive subset enumeration is bounded at dimension four")

    print("=== Cover-recursive peripheral finite diagnostic ===")
    catalogue = enumerate_partial_cubes(max_dimension)

    print("\nClass counts up to cube symmetry")
    print("dim   partial cubes   tree-like   R_per   periphery-free")
    for dimension, families in catalogue.items():
        tree_like = sum(is_tree_like(family, dimension) for family in families)
        recursive = sum(
            is_cover_recursive_peripheral(family, dimension)
            for family in families
        )
        periphery_free = sum(
            not peripheral_decompositions(family, dimension)
            for family in families
            if family.bit_count() > 1
        )
        print(
            f"{dimension:>3} {len(families):>15} {tree_like:>11} "
            f"{recursive:>7} {periphery_free:>17}"
        )

    print("\nImmediate pc-minor closure checks")
    checked_recursive = 0
    checked_tree_contractions = 0
    recursive_failures: list[tuple[int, int, int, int]] = []
    tree_contraction_failures: list[tuple[int, int, int, int]] = []
    for dimension, families in catalogue.items():
        for family in families:
            if is_cover_recursive_peripheral(family, dimension):
                for minor, minor_dimension in immediate_pc_minors(
                    family, dimension
                ):
                    checked_recursive += 1
                    if not is_cover_recursive_peripheral(minor, minor_dimension):
                        recursive_failures.append(
                            (family, dimension, minor, minor_dimension)
                        )
            if is_tree_like(family, dimension):
                for coordinate in range(dimension):
                    contraction = section(
                        family, dimension, coordinate, 0
                    ) | section(family, dimension, coordinate, 1)
                    contraction, contraction_dimension = canonical_family(
                        contraction, dimension - 1
                    )
                    checked_tree_contractions += 1
                    if not is_tree_like(contraction, contraction_dimension):
                        tree_contraction_failures.append(
                            (
                                family,
                                dimension,
                                contraction,
                                contraction_dimension,
                            )
                        )
    print(
        f"R_per immediate minors checked: {checked_recursive}; "
        f"failures: {len(recursive_failures)}"
    )
    print(
        f"tree-like contractions checked: {checked_tree_contractions}; "
        f"failures: {len(tree_contraction_failures)}"
    )

    c6 = family_from_vertices([1, 2, 3, 4, 5, 6])
    c6, c6_dimension = canonical_family(c6, 3)
    print("\nC6 pc-minor exclusion check")
    c6_failures: list[tuple[int, int]] = []
    for dimension, families in catalogue.items():
        for family in families:
            if is_cover_recursive_peripheral(family, dimension) and contains_pc_minor(
                family, dimension, c6, c6_dimension
            ):
                c6_failures.append((family, dimension))
    print(f"R_per members containing C6 as a pc-minor: {len(c6_failures)}")

    print("\nTree-like and C6-free comparison")
    print("dim   tree-like and C6-free   outside R_per")
    for dimension, families in catalogue.items():
        comparison = [
            family
            for family in families
            if is_tree_like(family, dimension)
            and not contains_pc_minor(family, dimension, c6, c6_dimension)
        ]
        outside = [
            family
            for family in comparison
            if not is_cover_recursive_peripheral(family, dimension)
        ]
        print(f"{dimension:>3} {len(comparison):>25} {len(outside):>15}")

    print("\npc-minor-minimal nonmembers in the enumerated range")
    minimal_nonmembers: list[tuple[int, int]] = []
    for dimension, families in catalogue.items():
        for family in families:
            if is_cover_recursive_peripheral(family, dimension):
                continue
            if all(
                is_cover_recursive_peripheral(minor, minor_dimension)
                for minor, minor_dimension in immediate_pc_minors(
                    family, dimension
                )
            ):
                minimal_nonmembers.append((family, dimension))
    for index, (family, dimension) in enumerate(minimal_nonmembers, start=1):
        print(
            f"{index:>2}. dim={dimension}, order={family.bit_count()}, "
            f"periphery_free={not peripheral_decompositions(family, dimension)}, "
            f"vertices={bitstrings(family, dimension)}"
        )
    if not minimal_nonmembers:
        print("none")

    print("\nCartesian product checks with total dimension in range")
    product_checks = 0
    product_failures: list[tuple[int, int, int, int]] = []
    for left_dimension, left_families in catalogue.items():
        for right_dimension, right_families in catalogue.items():
            if left_dimension + right_dimension > max_dimension:
                continue
            for left in left_families:
                for right in right_families:
                    product_checks += 1
                    product, product_dimension = cartesian_product(
                        left, left_dimension, right, right_dimension
                    )
                    expected = is_cover_recursive_peripheral(
                        left, left_dimension
                    ) and is_cover_recursive_peripheral(right, right_dimension)
                    actual = is_cover_recursive_peripheral(
                        product, product_dimension
                    )
                    if expected != actual:
                        product_failures.append(
                            (left, left_dimension, right, right_dimension)
                        )
    print(f"factorwise product tests: {product_checks}; failures: {len(product_failures)}")

    print("\nGated-amalgam counterexample search in the enumerated range")
    gated_amalgam_counterexamples: list[
        tuple[int, int, int, int]
    ] = []
    gated_subgraphs_tested = 0
    for dimension, families in catalogue.items():
        for family in families:
            if is_cover_recursive_peripheral(family, dimension):
                continue
            proper_gated_recursive: list[int] = []
            subfamily = family
            while subfamily:
                if subfamily != family:
                    gated_subgraphs_tested += 1
                    if is_gated_subfamily(subfamily, family, dimension):
                        normalized, normalized_dimension = canonical_family(
                            subfamily, dimension
                        )
                        if is_cover_recursive_peripheral(
                            normalized, normalized_dimension
                        ):
                            proper_gated_recursive.append(subfamily)
                subfamily = (subfamily - 1) & family
            found = False
            for left_index, left in enumerate(proper_gated_recursive):
                for right in proper_gated_recursive[left_index:]:
                    if left & right and left | right == family:
                        gated_amalgam_counterexamples.append(
                            (family, dimension, left, right)
                        )
                        found = True
                        break
                if found:
                    break
    print(
        f"proper subfamilies tested for gatedness: {gated_subgraphs_tested}; "
        f"non-R_per gated amalgams of two R_per pieces: "
        f"{len(gated_amalgam_counterexamples)}"
    )
    for family, dimension, left, right in gated_amalgam_counterexamples:
        print(
            f"counterexample dim={dimension}, order={family.bit_count()}: "
            f"G={bitstrings(family, dimension)}, "
            f"G0={bitstrings(left, dimension)}, "
            f"G1={bitstrings(right, dimension)}"
        )

    print("\nStructured antipodal-puncture family")
    print(
        " n   order   partial cube   periphery-free   R_per   "
        "all immediate minors R_per   C6-minor"
    )
    for dimension in range(3, 7):
        family = antipodally_punctured_cube(dimension)
        immediate = immediate_pc_minors(family, dimension)
        has_c6 = contains_pc_minor(family, dimension, c6, c6_dimension)
        print(
            f"{dimension:>2} {family.bit_count():>7} "
            f"{str(is_partial_cube_family(family, dimension)):>14} "
            f"{str(not peripheral_decompositions(family, dimension)):>16} "
            f"{str(is_cover_recursive_peripheral(family, dimension)):>7} "
            f"{str(all(is_cover_recursive_peripheral(minor, minor_dimension) for minor, minor_dimension in immediate)):>28} "
            f"{str(has_c6):>10}"
        )
        assert is_cover_recursive_peripheral(
            singly_punctured_cube(dimension - 1), dimension - 1
        )

    figure_one = antipodally_punctured_cube(4)
    envelope, envelope_dimension = peripheral_cube_envelope(figure_one, 4)
    print("\nAlmost-median tree-like envelope diagnostic")
    print(
        f"site=Q4 minus antipodes; envelope order={envelope.bit_count()}, "
        f"dimension={envelope_dimension}, "
        f"tree_like={is_tree_like(envelope, envelope_dimension)}, "
        f"R_per={is_cover_recursive_peripheral(envelope, envelope_dimension)}, "
        f"C6_minor={contains_pc_minor(envelope, envelope_dimension, c6, c6_dimension)}"
    )

    assert not recursive_failures
    assert not tree_contraction_failures
    assert not c6_failures
    assert not product_failures


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-dimension", type=int, default=4)
    args = parser.parse_args()
    run(args.max_dimension)


if __name__ == "__main__":
    main()
