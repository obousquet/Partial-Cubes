"""Search for an explicit peripheral-recursion obstruction beyond ``Q_n^{--}``.

The class ``R_per`` consists of the one-vertex partial cube and all partial
cubes admitting a peripheral coordinate whose contraction base and duplicated
site both belong to ``R_per``.  Its forbidden pc-minors are therefore the
pc-minor-minimal nonmembers, equivalently the minimal periphery-free partial
cubes.

This script follows pc-minors of Hall's 299-concept maximum class ``C_H``.
The concepts below are copied from ``anc/CH.txt`` in the source archive of

    Chalopin--Chepoi--Moran--Warmuth,
    Unlabeled sample compression schemes and corner peelings for ample and
    maximum classes, JCSS 127 (2022), 1--28,
    https://arxiv.org/src/1812.02099

The ancillary file has SHA256
``0367726b9934753a89711b42750d59f8a327860263383740ab85e6b6b27cc682``.

Every reported witness is independently checked to be a partial cube, ample,
periphery-free, outside ``R_per``, and to have every immediate restriction and
contraction in ``R_per``.  Thus finding a witness is a finite computer proof
for that explicitly listed family; claims about global smallestness require a
separate exhaustive census.

Run conservatively with, for example,

    nice -n 10 uv run python explore_smallest_additional_peripheral_obstruction.py
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import itertools
import json
import random
import resource
import time
from collections import Counter, defaultdict, deque
from functools import lru_cache
from pathlib import Path
from typing import Iterable

try:
    import pynauty
except ImportError:  # Optional exact accelerator; the Python fallback remains.
    pynauty = None

from explore_cover_recursive_peripheral import (
    enumerate_partial_cubes,
    is_gated_subfamily,
)


HALL_CONCEPTS = tuple(
    int(token)
    for token in """
0 1 2 3 16 17 19 48 49 51 55 59 63 64 65 67 80 81 83 112 113 115
119 123 127 131 147 179 187 191 192 193 195 208 209 211 240 241
243 244 245 247 251 252 253 255 256 272 304 305 307 311 319 368
369 371 375 383 447 496 497 499 500 501 503 508 509 511 563 567
571 575 703 767 768 772 776 780 784 788 796 800 804 808 812 816
817 819 820 821 823 824 825 827 828 829 831 880 881 883 884 885
887 892 893 895 959 1008 1009 1011 1012 1013 1015 1020 1021 1023
1024 1025 1026 1027 1040 1072 1088 1089 1091 1104 1136 1216 1217
1219 1232 1233 1235 1264 1265 1267 1268 1269 1271 1276 1277 1279
1280 1296 1328 1392 1520 1524 1532 1533 1535 1791 1792 1796 1800
1804 1808 1812 1820 1840 1844 1852 1904 1908 1916 2032 2036 2044
2045 2047 2048 2050 2051 2115 2179 2243 2304 2816 2824 2828 2844
2860 2876 3072 3073 3074 3075 3088 3120 3136 3137 3139 3152 3184
3200 3201 3202 3203 3264 3265 3266 3267 3268 3269 3270 3271 3276
3277 3278 3279 3280 3281 3283 3284 3285 3287 3292 3293 3295 3311
3312 3313 3315 3316 3317 3319 3324 3325 3327 3328 3344 3376 3392
3408 3440 3520 3524 3532 3533 3535 3536 3540 3548 3549 3551 3567
3568 3572 3580 3581 3583 3788 3789 3790 3791 3823 3839 3840 3844
3848 3852 3856 3860 3868 3884 3888 3892 3900 3904 3908 3916 3920
3924 3932 3948 3952 3956 3964 4032 4036 4044 4045 4046 4047 4048
4052 4060 4061 4062 4063 4076 4078 4079 4080 4084 4092 4093 4094
4095
""".split()
)


def family_from_vertices(members: Iterable[int]) -> int:
    family = 0
    for vertex in members:
        family |= 1 << vertex
    return family


def vertices(family: int, dimension: int) -> tuple[int, ...]:
    return tuple(vertex for vertex in range(1 << dimension) if family >> vertex & 1)


def project_vertex(vertex: int, coordinate: int) -> int:
    lower = (1 << coordinate) - 1
    return (vertex & lower) | ((vertex >> (coordinate + 1)) << coordinate)


@lru_cache(maxsize=None)
def section(family: int, dimension: int, coordinate: int, value: int) -> int:
    return family_from_vertices(
        project_vertex(vertex, coordinate)
        for vertex in vertices(family, dimension)
        if (vertex >> coordinate) & 1 == value
    )


@lru_cache(maxsize=None)
def irredundant_projection(family: int, dimension: int) -> tuple[int, int]:
    members = vertices(family, dimension)
    if not members:
        raise ValueError("the empty family is not a partial cube")
    varying = tuple(
        coordinate
        for coordinate in range(dimension)
        if any((vertex >> coordinate) & 1 for vertex in members)
        and any(not ((vertex >> coordinate) & 1) for vertex in members)
    )
    projected = family_from_vertices(
        sum(
            ((vertex >> old_coordinate) & 1) << new_coordinate
            for new_coordinate, old_coordinate in enumerate(varying)
        )
        for vertex in members
    )
    return projected, len(varying)


def irredundant_rooted_projection(
    family: int, dimension: int, root: int
) -> tuple[int, int, int]:
    """Remove constant coordinates while transporting a distinguished root."""

    members = vertices(family, dimension)
    if root not in members:
        raise ValueError("the root is not a member of the family")
    varying = tuple(
        coordinate
        for coordinate in range(dimension)
        if any((vertex >> coordinate) & 1 for vertex in members)
        and any(not ((vertex >> coordinate) & 1) for vertex in members)
    )

    def project_varying(vertex: int) -> int:
        return sum(
            ((vertex >> old_coordinate) & 1) << new_coordinate
            for new_coordinate, old_coordinate in enumerate(varying)
        )

    return (
        family_from_vertices(project_varying(vertex) for vertex in members),
        len(varying),
        project_varying(root),
    )


def nested(left: int, right: int) -> bool:
    return left & ~right == 0


@lru_cache(maxsize=None)
def peripheral_decompositions(
    family: int, dimension: int
) -> tuple[tuple[int, int, int, int], ...]:
    family, dimension = irredundant_projection(family, dimension)
    answer: list[tuple[int, int, int, int]] = []
    for coordinate in range(dimension):
        zero = section(family, dimension, coordinate, 0)
        one = section(family, dimension, coordinate, 1)
        if nested(zero, one):
            answer.append((coordinate, one, zero, dimension - 1))
        elif nested(one, zero):
            answer.append((coordinate, zero, one, dimension - 1))
    return tuple(answer)


@lru_cache(maxsize=None)
def is_r_per(family: int, dimension: int) -> bool:
    family, dimension = irredundant_projection(family, dimension)
    if family.bit_count() == 1:
        return True
    return any(
        is_r_per(base, lower_dimension) and is_r_per(site, lower_dimension)
        for _, base, site, lower_dimension in peripheral_decompositions(
            family, dimension
        )
    )


def r_per_certificate(family: int, dimension: int) -> object:
    """Return one recursive peripheral-decomposition certificate."""

    family, dimension = irredundant_projection(family, dimension)
    if family.bit_count() == 1:
        return "K1"
    for coordinate, base, site, lower_dimension in peripheral_decompositions(
        family, dimension
    ):
        if is_r_per(base, lower_dimension) and is_r_per(
            site, lower_dimension
        ):
            return {
                "coordinate": coordinate,
                "base": r_per_certificate(base, lower_dimension),
                "site": r_per_certificate(site, lower_dimension),
            }
    raise ValueError("the supplied family is not in R_per")


@lru_cache(maxsize=None)
def is_rooted_r_per(
    family: int, dimension: int, root: int
) -> bool:
    """Whether ``R_per`` has a certificate avoiding a fixed root.

    At every recursion node the chosen peripheral halfspace must avoid the
    current image of the root.  Such a certificate can be executed inside a
    one-vertex amalgam without the other factor spoiling the periphery.
    """

    family, dimension, root = irredundant_rooted_projection(
        family, dimension, root
    )
    if family.bit_count() == 1:
        return True
    for coordinate in range(dimension):
        zero = section(family, dimension, coordinate, 0)
        one = section(family, dimension, coordinate, 1)
        root_value = (root >> coordinate) & 1
        root_projection = project_vertex(root, coordinate)
        for peripheral_value, site, base in (
            (0, zero, one),
            (1, one, zero),
        ):
            if root_value == peripheral_value or not nested(site, base):
                continue
            if is_r_per(site, dimension - 1) and is_rooted_r_per(
                base, dimension - 1, root_projection
            ):
                return True
    return False


def rooted_r_per_certificate(
    family: int, dimension: int, root: int
) -> object:
    """Emit one root-avoiding recursive peripheral certificate."""

    family, dimension, root = irredundant_rooted_projection(
        family, dimension, root
    )
    if family.bit_count() == 1:
        return "root"
    for coordinate in range(dimension):
        zero = section(family, dimension, coordinate, 0)
        one = section(family, dimension, coordinate, 1)
        root_value = (root >> coordinate) & 1
        root_projection = project_vertex(root, coordinate)
        for peripheral_value, site, base in (
            (0, zero, one),
            (1, one, zero),
        ):
            if root_value == peripheral_value or not nested(site, base):
                continue
            if is_r_per(site, dimension - 1) and is_rooted_r_per(
                base, dimension - 1, root_projection
            ):
                return {
                    "coordinate": coordinate,
                    "peripheral_value": peripheral_value,
                    "base": rooted_r_per_certificate(
                        base, dimension - 1, root_projection
                    ),
                    "site": r_per_certificate(site, dimension - 1),
                }
    raise ValueError("the supplied rooted family has no avoiding certificate")


def root_avoiding_peripheral_coordinates(
    family: int, dimension: int, root: int
) -> tuple[int, ...]:
    """Coordinates whose peripheral semicube omits the distinguished root."""

    family, dimension, root = irredundant_rooted_projection(
        family, dimension, root
    )
    answer: list[int] = []
    for coordinate in range(dimension):
        root_value = (root >> coordinate) & 1
        root_side = section(
            family, dimension, coordinate, root_value
        )
        opposite_side = section(
            family, dimension, coordinate, 1 - root_value
        )
        if nested(opposite_side, root_side):
            answer.append(coordinate)
    return tuple(answer)


def raw_root_avoiding_peripheral_coordinates(
    family: int,
    dimension: int,
    root: int,
) -> tuple[int, ...]:
    """Available coordinates without deleting constants or renumbering.

    This common-coordinate version is needed for a contained pair
    ``site <= base``: a coordinate can be constant in one member while
    varying in the other, so projecting the two members separately would
    destroy their coordinate alignment.
    """

    if root not in vertices(family, dimension):
        raise ValueError("the root is not a member of the family")
    return tuple(
        coordinate
        for coordinate in range(dimension)
        if nested(
            section(
                family,
                dimension,
                coordinate,
                1 - ((root >> coordinate) & 1),
            ),
            section(
                family,
                dimension,
                coordinate,
                (root >> coordinate) & 1,
            ),
        )
    )


@lru_cache(maxsize=None)
def contract_coordinate(
    family: int,
    dimension: int,
    coordinate: int,
) -> int:
    """Contract one cube coordinate while retaining the common labelling."""

    return (
        section(family, dimension, coordinate, 0)
        | section(family, dimension, coordinate, 1)
    )


@lru_cache(maxsize=None)
def has_root_elimination_order(
    family: int,
    dimension: int,
    root: int,
) -> bool:
    """Test rooted recursion using available-coordinate contractions only."""

    family, dimension, root = irredundant_rooted_projection(
        family,
        dimension,
        root,
    )
    if not is_r_per(family, dimension):
        return False
    if family.bit_count() == 1:
        return True
    return any(
        has_root_elimination_order(
            contract_coordinate(family, dimension, coordinate),
            dimension - 1,
            project_vertex(root, coordinate),
        )
        for coordinate in raw_root_avoiding_peripheral_coordinates(
            family,
            dimension,
            root,
        )
    )


@lru_cache(maxsize=None)
def is_synchronizable_reverse_pair(
    base: int,
    site: int,
    dimension: int,
    root: int,
) -> bool:
    """Whether a reverse peripheral pair admits joint rooted elimination."""

    if site & ~base:
        raise ValueError("the reverse-expansion site must lie in its base")
    if root not in vertices(site, dimension):
        raise ValueError("the root must lie in the reverse-expansion site")
    if not (
        is_r_per(base, dimension)
        and is_r_per(site, dimension)
    ):
        return False
    if base == site:
        return is_rooted_r_per(base, dimension, root)

    common_available = (
        set(
            raw_root_avoiding_peripheral_coordinates(
                base,
                dimension,
                root,
            )
        )
        & set(
            raw_root_avoiding_peripheral_coordinates(
                site,
                dimension,
                root,
            )
        )
    )
    return any(
        is_synchronizable_reverse_pair(
            contract_coordinate(base, dimension, coordinate),
            contract_coordinate(site, dimension, coordinate),
            dimension - 1,
            project_vertex(root, coordinate),
        )
        for coordinate in common_available
    )


@lru_cache(maxsize=None)
def common_availability_closure(
    base: int,
    site: int,
    dimension: int,
    root: int,
) -> tuple[int, int, int, int, int]:
    """Greedily contract the unique maximal common-feasible coordinate set.

    The return value is
    ``(contracted_mask, final_base, final_site, final_dimension, final_root)``.
    The mask uses the original coordinate labels.  Availability is permanent
    under other contractions, so the maximal contracted set and the two final
    projections are independent of the greedy tie-breaking order.
    """

    if site & ~base:
        raise ValueError("the reverse-expansion site must lie in its base")
    if root not in vertices(site, dimension):
        raise ValueError("the root must lie in the reverse-expansion site")

    current_base = base
    current_site = site
    current_dimension = dimension
    current_root = root
    labels = list(range(dimension))
    contracted_mask = 0

    while True:
        common_available = sorted(
            set(
                raw_root_avoiding_peripheral_coordinates(
                    current_base,
                    current_dimension,
                    current_root,
                )
            )
            & set(
                raw_root_avoiding_peripheral_coordinates(
                    current_site,
                    current_dimension,
                    current_root,
                )
            ),
            key=lambda coordinate: labels[coordinate],
        )
        if not common_available:
            return (
                contracted_mask,
                current_base,
                current_site,
                current_dimension,
                current_root,
            )

        coordinate = common_available[0]
        contracted_mask |= 1 << labels[coordinate]
        current_base = contract_coordinate(
            current_base,
            current_dimension,
            coordinate,
        )
        current_site = contract_coordinate(
            current_site,
            current_dimension,
            coordinate,
        )
        current_root = project_vertex(current_root, coordinate)
        current_dimension -= 1
        del labels[coordinate]


@lru_cache(maxsize=None)
def common_availability_closure_avoiding(
    base: int,
    site: int,
    dimension: int,
    root: int,
    forbidden_mask: int,
) -> tuple[int, int, int, int, int]:
    """Greedily close while forbidding selected original coordinates."""

    if site & ~base:
        raise ValueError("the reverse-expansion site must lie in its base")
    if root not in vertices(site, dimension):
        raise ValueError("the root must lie in the reverse-expansion site")
    if forbidden_mask & ~((1 << dimension) - 1):
        raise ValueError("the forbidden mask uses an unknown coordinate")

    current_base = base
    current_site = site
    current_dimension = dimension
    current_root = root
    labels = list(range(dimension))
    contracted_mask = 0

    while True:
        common_available = sorted(
            (
                coordinate
                for coordinate in (
                    set(
                        raw_root_avoiding_peripheral_coordinates(
                            current_base,
                            current_dimension,
                            current_root,
                        )
                    )
                    & set(
                        raw_root_avoiding_peripheral_coordinates(
                            current_site,
                            current_dimension,
                            current_root,
                        )
                    )
                )
                if not (forbidden_mask & (1 << labels[coordinate]))
            ),
            key=lambda coordinate: labels[coordinate],
        )
        if not common_available:
            return (
                contracted_mask,
                current_base,
                current_site,
                current_dimension,
                current_root,
            )

        coordinate = common_available[0]
        contracted_mask |= 1 << labels[coordinate]
        current_base = contract_coordinate(
            current_base,
            current_dimension,
            coordinate,
        )
        current_site = contract_coordinate(
            current_site,
            current_dimension,
            coordinate,
        )
        current_root = project_vertex(current_root, coordinate)
        current_dimension -= 1
        del labels[coordinate]


@lru_cache(maxsize=None)
def common_availability_antimatroid(
    base: int,
    site: int,
    dimension: int,
    root: int,
) -> tuple[int, ...]:
    """Return all coordinate sets admitting a common-availability order."""

    if site & ~base:
        raise ValueError("the reverse-expansion site must lie in its base")
    if root not in vertices(site, dimension):
        raise ValueError("the root must lie in the reverse-expansion site")

    states: dict[
        int,
        tuple[int, int, int, int, tuple[int, ...]],
    ] = {
        0: (
            base,
            site,
            dimension,
            root,
            tuple(range(dimension)),
        )
    }
    queue = deque([0])
    while queue:
        contracted_mask = queue.popleft()
        (
            current_base,
            current_site,
            current_dimension,
            current_root,
            labels,
        ) = states[contracted_mask]
        common_available = (
            set(
                raw_root_avoiding_peripheral_coordinates(
                    current_base,
                    current_dimension,
                    current_root,
                )
            )
            & set(
                raw_root_avoiding_peripheral_coordinates(
                    current_site,
                    current_dimension,
                    current_root,
                )
            )
        )
        for coordinate in sorted(
            common_available,
            key=lambda current: labels[current],
        ):
            original_coordinate = labels[coordinate]
            next_mask = contracted_mask | (1 << original_coordinate)
            if next_mask in states:
                continue
            next_labels = (
                labels[:coordinate] + labels[coordinate + 1 :]
            )
            states[next_mask] = (
                contract_coordinate(
                    current_base,
                    current_dimension,
                    coordinate,
                ),
                contract_coordinate(
                    current_site,
                    current_dimension,
                    coordinate,
                ),
                current_dimension - 1,
                project_vertex(current_root, coordinate),
                next_labels,
            )
            queue.append(next_mask)
    return tuple(sorted(states))


@lru_cache(maxsize=None)
def is_synchronizable_by_common_closure(
    base: int,
    site: int,
    dimension: int,
    root: int,
) -> bool:
    """Test synchronization by the deterministic common-availability closure."""

    if not (
        is_r_per(base, dimension)
        and is_r_per(site, dimension)
    ):
        return False
    contracted_mask = common_availability_closure(
        base,
        site,
        dimension,
        root,
    )[0]
    return contracted_mask == (1 << dimension) - 1


def is_q_activated(
    base: int,
    site: int,
    dimension: int,
    root: int,
    q: int,
) -> bool:
    """Whether the marked coordinate belongs to the common closure."""

    if q < 0 or q >= dimension:
        raise ValueError("the activation coordinate is outside the cube")
    return bool(
        common_availability_closure(
            base,
            site,
            dimension,
            root,
        )[0]
        & (1 << q)
    )


def is_q_activated_by_q_free_fibers(
    base: int,
    site: int,
    dimension: int,
    root: int,
    q: int,
) -> bool:
    """Test marked activation after the maximal ``q``-free closure."""

    if q < 0 or q >= dimension:
        raise ValueError("the activation coordinate is outside the cube")
    (
        contracted_mask,
        final_base,
        final_site,
        final_dimension,
        final_root,
    ) = common_availability_closure_avoiding(
        base,
        site,
        dimension,
        root,
        1 << q,
    )
    projected_q = q - (
        contracted_mask & ((1 << q) - 1)
    ).bit_count()
    if projected_q < 0 or projected_q >= final_dimension:
        raise AssertionError("the q-free closure contracted q")
    return (
        projected_q
        in raw_root_avoiding_peripheral_coordinates(
            final_base,
            final_dimension,
            final_root,
        )
        and projected_q
        in raw_root_avoiding_peripheral_coordinates(
            final_site,
            final_dimension,
            final_root,
        )
    )


@lru_cache(maxsize=None)
def q_defect_support_system(
    family: int,
    dimension: int,
    root: int,
    q: int,
) -> tuple[tuple[int, tuple[int, ...]], ...]:
    """Return every labelled ``q`` defect and its minimal repair supports.

    The two ``q``-sections live in a cube of dimension ``dimension - 1``.
    For an opposite-side projection ``u`` absent from the root-side section,
    a contraction set repairs ``u`` exactly when it contains the support of
    ``u xor v`` for some root-side projection ``v``.  Each returned row is
    ``(u, clause)``, where ``clause`` is the antichain of inclusion-minimal
    repair supports for the labelled defect ``u``.
    """

    if q < 0 or q >= dimension:
        raise ValueError("the repair coordinate is outside the cube")
    if root not in vertices(family, dimension):
        raise ValueError("the root is not a member of the family")
    root_value = (root >> q) & 1
    root_side = section(
        family,
        dimension,
        q,
        root_value,
    )
    opposite_side = section(
        family,
        dimension,
        q,
        1 - root_value,
    )
    root_vertices = vertices(root_side, dimension - 1)
    root_vertex_set = set(root_vertices)
    defect_system: list[tuple[int, tuple[int, ...]]] = []
    for opposite_vertex in vertices(
        opposite_side,
        dimension - 1,
    ):
        if opposite_vertex in root_vertex_set:
            continue
        candidate_supports = {
            opposite_vertex ^ root_vertex
            for root_vertex in root_vertices
        }
        minimal_supports = tuple(
            sorted(
                support
                for support in candidate_supports
                if not any(
                    other != support
                    and not (other & ~support)
                    for other in candidate_supports
                )
            )
        )
        if not minimal_supports or 0 in minimal_supports:
            raise AssertionError("an unmatched q-fiber has no proper clause")
        defect_system.append((opposite_vertex, minimal_supports))
    return tuple(defect_system)


@lru_cache(maxsize=None)
def minimal_q_repair_supports(
    family: int,
    dimension: int,
    root: int,
    q: int,
) -> tuple[tuple[int, ...], ...]:
    """Return the unlabelled multiset of minimal ``q``-repair clauses."""

    return tuple(
        sorted(
            clause
            for _, clause in q_defect_support_system(
                family,
                dimension,
                root,
                q,
            )
        )
    )


def inclusion_minimal_masks(masks: Iterable[int]) -> tuple[int, ...]:
    """Return the inclusion-minimal members of a finite mask family."""

    distinct = set(masks)
    return tuple(
        sorted(
            mask
            for mask in distinct
            if not any(
                other != mask
                and not (other & ~mask)
                for other in distinct
            )
        )
    )


def contract_q_defect_support_system(
    defect_system: tuple[tuple[int, tuple[int, ...]], ...],
    fiber_dimension: int,
    coordinate: int,
) -> tuple[tuple[int, tuple[int, ...]], ...]:
    """Compute the signed defect system after contracting one coordinate."""

    if coordinate < 0 or coordinate >= fiber_dimension:
        raise ValueError("the contracted fiber coordinate is outside the cube")
    singleton = 1 << coordinate
    transformed: dict[int, tuple[int, ...]] = {}
    for defect, clause in defect_system:
        if singleton in clause:
            continue
        projected_defect = project_vertex(defect, coordinate)
        projected_clause = inclusion_minimal_masks(
            project_vertex(support, coordinate)
            for support in clause
        )
        if not projected_clause or 0 in projected_clause:
            raise AssertionError(
                "a surviving contracted defect has an empty repair support"
            )
        previous = transformed.setdefault(
            projected_defect,
            projected_clause,
        )
        if previous != projected_clause:
            raise AssertionError(
                "two lifts of one contracted defect have different clauses"
            )
    return tuple(sorted(transformed.items()))


def restrict_q_defect_support_system_to_root(
    defect_system: tuple[tuple[int, tuple[int, ...]], ...],
    fiber_dimension: int,
    coordinate: int,
    root_value: int,
) -> tuple[tuple[int, tuple[int, ...]], ...]:
    """Compute the signed defect system in a root-side restriction."""

    if coordinate < 0 or coordinate >= fiber_dimension:
        raise ValueError("the restricted fiber coordinate is outside the cube")
    if root_value not in (0, 1):
        raise ValueError("a root coordinate value must be zero or one")
    restricted: list[tuple[int, tuple[int, ...]]] = []
    singleton = 1 << coordinate
    for defect, clause in defect_system:
        if ((defect >> coordinate) & 1) != root_value:
            continue
        projected_clause = tuple(
            project_vertex(support, coordinate)
            for support in clause
            if not (support & singleton)
        )
        if not projected_clause:
            raise AssertionError(
                "a root-side defect has no root-side repair support"
            )
        restricted.append(
            (
                project_vertex(defect, coordinate),
                projected_clause,
            )
        )
    return tuple(restricted)


@lru_cache(maxsize=None)
def rooted_defect_atlas(
    family: int,
    dimension: int,
    root: int,
) -> tuple[tuple[tuple[int, tuple[int, ...]], ...], ...]:
    """Return the signed defect-support system in every coordinate."""

    return tuple(
        q_defect_support_system(
            family,
            dimension,
            root,
            coordinate,
        )
        for coordinate in range(dimension)
    )


def contract_rooted_defect_atlas(
    atlas: tuple[tuple[tuple[int, tuple[int, ...]], ...], ...],
    dimension: int,
    coordinate: int,
) -> tuple[tuple[tuple[int, tuple[int, ...]], ...], ...]:
    """Apply a coordinate contraction to every surviving atlas entry."""

    if len(atlas) != dimension:
        raise ValueError("the defect atlas has the wrong dimension")
    if coordinate < 0 or coordinate >= dimension:
        raise ValueError("the contracted coordinate is outside the atlas")
    return tuple(
        contract_q_defect_support_system(
            atlas[marked_coordinate],
            dimension - 1,
            coordinate - (coordinate > marked_coordinate),
        )
        for marked_coordinate in range(dimension)
        if marked_coordinate != coordinate
    )


def restrict_rooted_defect_atlas_to_root(
    atlas: tuple[tuple[tuple[int, tuple[int, ...]], ...], ...],
    dimension: int,
    root: int,
    coordinate: int,
) -> tuple[tuple[tuple[int, tuple[int, ...]], ...], ...]:
    """Apply a root-side restriction to every surviving atlas entry."""

    if len(atlas) != dimension:
        raise ValueError("the defect atlas has the wrong dimension")
    if coordinate < 0 or coordinate >= dimension:
        raise ValueError("the restricted coordinate is outside the atlas")
    root_value = (root >> coordinate) & 1
    return tuple(
        restrict_q_defect_support_system_to_root(
            atlas[marked_coordinate],
            dimension - 1,
            coordinate - (coordinate > marked_coordinate),
            root_value,
        )
        for marked_coordinate in range(dimension)
        if marked_coordinate != coordinate
    )


def common_availability_closure_from_defect_atlases(
    base_atlas: tuple[
        tuple[tuple[int, tuple[int, ...]], ...],
        ...,
    ],
    site_atlas: tuple[
        tuple[tuple[int, tuple[int, ...]], ...],
        ...,
    ],
    dimension: int,
    forbidden_mask: int,
) -> int:
    """Compute the common closure using only two rooted defect atlases."""

    if len(base_atlas) != dimension or len(site_atlas) != dimension:
        raise ValueError("the defect atlases have the wrong dimension")
    if forbidden_mask & ~((1 << dimension) - 1):
        raise ValueError("the forbidden mask uses an unknown coordinate")
    current_base = base_atlas
    current_site = site_atlas
    current_dimension = dimension
    labels = list(range(dimension))
    contracted_mask = 0
    while True:
        common_available = [
            coordinate
            for coordinate in range(current_dimension)
            if not current_base[coordinate]
            and not current_site[coordinate]
            and not (forbidden_mask & (1 << labels[coordinate]))
        ]
        if not common_available:
            return contracted_mask
        coordinate = min(
            common_available,
            key=lambda current: labels[current],
        )
        contracted_mask |= 1 << labels[coordinate]
        current_base = contract_rooted_defect_atlas(
            current_base,
            current_dimension,
            coordinate,
        )
        current_site = contract_rooted_defect_atlas(
            current_site,
            current_dimension,
            coordinate,
        )
        current_dimension -= 1
        del labels[coordinate]


def insert_vertex_coordinate(
    vertex: int,
    coordinate: int,
    value: int,
) -> int:
    """Insert one coordinate bit into a cube vertex."""

    if coordinate < 0 or value not in (0, 1):
        raise ValueError("invalid inserted coordinate")
    lower_mask = (1 << coordinate) - 1
    return (
        (vertex & lower_mask)
        | (value << coordinate)
        | ((vertex & ~lower_mask) << 1)
    )


def defect_atlas_horn_constraints(
    atlas: tuple[tuple[tuple[int, tuple[int, ...]], ...], ...],
    dimension: int,
    root: int,
) -> tuple[int, int, tuple[tuple[int, int], ...]]:
    """Encode exact atlas realization as forced/forbidden Horn data."""

    if len(atlas) != dimension:
        raise ValueError("the defect atlas has the wrong dimension")
    if root < 0 or root >= 1 << dimension:
        raise ValueError("the root is outside the cube")
    forced_mask = 1 << root
    forbidden_mask = 0
    implications: set[tuple[int, int]] = set()
    for coordinate, defect_system in enumerate(atlas):
        root_value = (root >> coordinate) & 1
        opposite_value = 1 - root_value
        defects = dict(defect_system)
        fiber_dimension = dimension - 1
        for defect, clause in defect_system:
            forced_mask |= 1 << insert_vertex_coordinate(
                defect,
                coordinate,
                opposite_value,
            )
            for support in clause:
                forced_mask |= 1 << insert_vertex_coordinate(
                    defect ^ support,
                    coordinate,
                    root_value,
                )

        for root_projection in range(1 << fiber_dimension):
            allowed = all(
                any(
                    not (
                        support
                        & ~(defect ^ root_projection)
                    )
                    for support in clause
                )
                for defect, clause in defect_system
            )
            root_vertex = insert_vertex_coordinate(
                root_projection,
                coordinate,
                root_value,
            )
            if not allowed:
                forbidden_mask |= 1 << root_vertex
            if root_projection not in defects:
                opposite_vertex = insert_vertex_coordinate(
                    root_projection,
                    coordinate,
                    opposite_value,
                )
                implications.add((opposite_vertex, root_vertex))
    return forced_mask, forbidden_mask, tuple(sorted(implications))


def horn_implication_closure(
    seed_mask: int,
    implications: tuple[tuple[int, int], ...],
) -> int:
    """Close a vertex family under unary Horn implications."""

    closure = seed_mask
    while True:
        next_closure = closure
        for source, target in implications:
            if closure & (1 << source):
                next_closure |= 1 << target
        if next_closure == closure:
            return closure
        closure = next_closure


def family_satisfies_atlas_horn_constraints(
    family: int,
    forced_mask: int,
    forbidden_mask: int,
    implications: tuple[tuple[int, int], ...],
) -> bool:
    """Test one family against an atlas's exact Horn realization system."""

    return (
        not (forced_mask & ~family)
        and not (forbidden_mask & family)
        and all(
            not (family & (1 << source))
            or bool(family & (1 << target))
            for source, target in implications
        )
    )


def maximal_safe_horn_assignment(
    variable_count: int,
    forbidden_mask: int,
    implications: tuple[tuple[int, int], ...],
) -> int:
    """Return all variables whose implication closure avoids a forbidden one."""

    safe_mask = 0
    for variable in range(variable_count):
        closure = horn_implication_closure(
            1 << variable,
            implications,
        )
        if not (closure & forbidden_mask):
            safe_mask |= 1 << variable
    return safe_mask


def contained_pair_atlas_horn_constraints(
    base_atlas: tuple[
        tuple[tuple[int, tuple[int, ...]], ...],
        ...,
    ],
    site_atlas: tuple[
        tuple[tuple[int, tuple[int, ...]], ...],
        ...,
    ],
    dimension: int,
    root: int,
) -> tuple[int, int, tuple[tuple[int, int], ...]]:
    """Encode simultaneous atlas realization with site contained in base."""

    (
        base_forced,
        base_forbidden,
        base_implications,
    ) = defect_atlas_horn_constraints(
        base_atlas,
        dimension,
        root,
    )
    (
        site_forced,
        site_forbidden,
        site_implications,
    ) = defect_atlas_horn_constraints(
        site_atlas,
        dimension,
        root,
    )
    vertex_count = 1 << dimension
    forced_mask = base_forced | (site_forced << vertex_count)
    forbidden_mask = (
        base_forbidden
        | (site_forbidden << vertex_count)
    )
    implications = set(base_implications)
    implications.update(
        (source + vertex_count, target + vertex_count)
        for source, target in site_implications
    )
    implications.update(
        (vertex + vertex_count, vertex)
        for vertex in range(vertex_count)
    )
    return forced_mask, forbidden_mask, tuple(sorted(implications))


def q_repair_supports_cover(
    repair_supports: tuple[tuple[int, ...], ...],
    contraction_mask: int,
) -> bool:
    """Whether a contraction mask satisfies every repair-support clause."""

    return all(
        any(not (support & ~contraction_mask) for support in clause)
        for clause in repair_supports
    )


def contract_original_coordinate_mask(
    family: int,
    dimension: int,
    root: int,
    contraction_mask: int,
) -> tuple[int, int, int]:
    """Contract a set named in the original coordinate system."""

    if contraction_mask & ~((1 << dimension) - 1):
        raise ValueError("the contraction mask uses an unknown coordinate")
    current_family = family
    current_root = root
    current_dimension = dimension
    for coordinate in reversed(range(dimension)):
        if not (contraction_mask & (1 << coordinate)):
            continue
        current_family = contract_coordinate(
            current_family,
            current_dimension,
            coordinate,
        )
        current_root = project_vertex(current_root, coordinate)
        current_dimension -= 1
    return current_family, current_dimension, current_root


def marked_fiber_state_signature(
    base: int,
    site: int,
    dimension: int,
    root: int,
    q: int,
    mode: str,
) -> tuple[object, ...]:
    """Encode q-free feasible sets and q-fiber repair incidence."""

    if mode not in {"closure", "logical", "multiset", "incidence"}:
        raise ValueError(f"unknown marked fiber signature mode: {mode}")
    q_free_feasible = tuple(
        project_vertex(mask, q)
        for mask in common_availability_antimatroid(
            base,
            site,
            dimension,
            root,
        )
        if not (mask & (1 << q))
    )

    def family_payload(family: int) -> object:
        clauses = minimal_q_repair_supports(
            family,
            dimension,
            root,
            q,
        )
        if mode == "closure":
            return ()
        if mode == "logical":
            return tuple(sorted(set(clauses)))
        if mode == "multiset":
            return clauses
        root_value = (root >> q) & 1
        return (
            section(
                family,
                dimension,
                q,
                root_value,
            ).bit_count(),
            section(
                family,
                dimension,
                q,
                1 - root_value,
            ).bit_count(),
            clauses,
        )

    return (
        q,
        q_free_feasible,
        family_payload(base),
        family_payload(site),
    )


def marked_fiber_witness_signature(
    base: int,
    site: int,
    dimension: int,
    q: int,
    mode: str,
) -> tuple[object, ...]:
    """Fiber data for a marked pair and every forced minor retaining q."""

    initial = marked_fiber_state_signature(
        base,
        site,
        dimension,
        0,
        q,
        mode,
    )
    forced_rows: list[tuple[object, object]] = []
    for coordinate in range(dimension):
        if coordinate == q:
            continue
        projected_q = q - (q > coordinate)
        forced_rows.append(
            (
                marked_fiber_state_signature(
                    contract_coordinate(
                        base,
                        dimension,
                        coordinate,
                    ),
                    contract_coordinate(
                        site,
                        dimension,
                        coordinate,
                    ),
                    dimension - 1,
                    0,
                    projected_q,
                    mode,
                ),
                marked_fiber_state_signature(
                    section(
                        base,
                        dimension,
                        coordinate,
                        0,
                    ),
                    section(
                        site,
                        dimension,
                        coordinate,
                        0,
                    ),
                    dimension - 1,
                    0,
                    projected_q,
                    mode,
                ),
            )
        )
    return initial, tuple(forced_rows)


@lru_cache(maxsize=None)
def canonical_marked_fiber_witness_signature(
    base: int,
    site: int,
    dimension: int,
    q: int,
    mode: str,
) -> tuple[object, ...]:
    """Canonical marked fiber signature under coordinate permutations."""

    return min(
        marked_fiber_witness_signature(
            permute_family_by_coordinate_order(
                base,
                dimension,
                order,
            ),
            permute_family_by_coordinate_order(
                site,
                dimension,
                order,
            ),
            dimension,
            order.index(q),
            mode,
        )
        for order in itertools.permutations(range(dimension))
    )


def q_section_orders(
    family: int,
    dimension: int,
    root: int,
    q: int,
) -> tuple[int, int]:
    """Return the root-side and opposite-side ``q``-section orders."""

    root_value = (root >> q) & 1
    return (
        section(
            family,
            dimension,
            q,
            root_value,
        ).bit_count(),
        section(
            family,
            dimension,
            q,
            1 - root_value,
        ).bit_count(),
    )


def marked_q_free_closure_payload(
    base: int,
    site: int,
    dimension: int,
    root: int,
    q: int,
    include_section_orders: bool,
) -> tuple[object, ...]:
    """Encode the maximal q-free closure, optionally with section orders."""

    contracted_mask = common_availability_closure_avoiding(
        base,
        site,
        dimension,
        root,
        1 << q,
    )[0]
    payload: tuple[object, ...] = (
        project_vertex(contracted_mask, q),
    )
    if include_section_orders:
        payload += (
            q_section_orders(
                base,
                dimension,
                root,
                q,
            ),
            q_section_orders(
                site,
                dimension,
                root,
                q,
            ),
        )
    return payload


def marked_signed_repair_signature(
    base: int,
    site: int,
    dimension: int,
    q: int,
    mode: str,
) -> tuple[object, ...]:
    """Encode a signed initial repair system and exact minor certificates."""

    if mode not in {
        "signed",
        "signed_orders",
        "signed_section_profile",
        "certificate",
        "certificate_initial_orders",
        "certificate_orders",
    }:
        raise ValueError(f"unknown signed repair signature mode: {mode}")
    include_certificates = mode in {
        "certificate",
        "certificate_initial_orders",
        "certificate_orders",
    }
    include_initial_orders = mode in {
        "signed_orders",
        "signed_section_profile",
        "certificate_initial_orders",
        "certificate_orders",
    }
    include_forced_orders = mode in {
        "signed_section_profile",
        "certificate_orders",
    }
    initial: tuple[object, ...] = (
        q_defect_support_system(
            base,
            dimension,
            0,
            q,
        ),
        q_defect_support_system(
            site,
            dimension,
            0,
            q,
        ),
    )
    if include_initial_orders:
        initial += (
            q_section_orders(base, dimension, 0, q),
            q_section_orders(site, dimension, 0, q),
        )
    if not include_certificates and not include_forced_orders:
        return q, initial

    if include_certificates:
        initial += (
            marked_q_free_closure_payload(
                base,
                site,
                dimension,
                0,
                q,
                False,
            ),
        )
    forced_rows: list[tuple[object, object]] = []
    for coordinate in range(dimension):
        if coordinate == q:
            continue
        projected_q = q - (q > coordinate)
        contracted_base = contract_coordinate(
            base,
            dimension,
            coordinate,
        )
        contracted_site = contract_coordinate(
            site,
            dimension,
            coordinate,
        )
        restricted_base = section(
            base,
            dimension,
            coordinate,
            0,
        )
        restricted_site = section(
            site,
            dimension,
            coordinate,
            0,
        )
        if not include_certificates:
            forced_rows.append(
                (
                    (
                        q_section_orders(
                            contracted_base,
                            dimension - 1,
                            0,
                            projected_q,
                        ),
                        q_section_orders(
                            contracted_site,
                            dimension - 1,
                            0,
                            projected_q,
                        ),
                    ),
                    (
                        q_section_orders(
                            restricted_base,
                            dimension - 1,
                            0,
                            projected_q,
                        ),
                        q_section_orders(
                            restricted_site,
                            dimension - 1,
                            0,
                            projected_q,
                        ),
                    ),
                )
            )
            continue
        forced_rows.append(
            (
                marked_q_free_closure_payload(
                    contracted_base,
                    contracted_site,
                    dimension - 1,
                    0,
                    projected_q,
                    include_forced_orders,
                ),
                marked_q_free_closure_payload(
                    restricted_base,
                    restricted_site,
                    dimension - 1,
                    0,
                    projected_q,
                    include_forced_orders,
                ),
            )
        )
    return q, initial, tuple(forced_rows)


@lru_cache(maxsize=None)
def canonical_marked_signed_repair_signature(
    base: int,
    site: int,
    dimension: int,
    q: int,
    mode: str,
) -> tuple[object, ...]:
    """Canonicalize a signed repair signature under marked permutations."""

    return min(
        marked_signed_repair_signature(
            permute_family_by_coordinate_order(
                base,
                dimension,
                order,
            ),
            permute_family_by_coordinate_order(
                site,
                dimension,
                order,
            ),
            dimension,
            order.index(q),
            mode,
        )
        for order in itertools.permutations(range(dimension))
    )


def marked_defect_atlas_signature(
    base: int,
    site: int,
    dimension: int,
    q: int,
    mode: str,
) -> tuple[object, ...]:
    """Encode the full rooted defect atlases of a marked pair."""

    if mode not in {"atlas", "atlas_orders", "atlas_section_orders"}:
        raise ValueError(f"unknown defect-atlas signature mode: {mode}")
    signature: tuple[object, ...] = (
        q,
        rooted_defect_atlas(base, dimension, 0),
        rooted_defect_atlas(site, dimension, 0),
    )
    if mode == "atlas_orders":
        signature += (
            base.bit_count(),
            site.bit_count(),
        )
    if mode == "atlas_section_orders":
        signature += (
            tuple(
                q_section_orders(
                    base,
                    dimension,
                    0,
                    coordinate,
                )
                for coordinate in range(dimension)
            ),
            tuple(
                q_section_orders(
                    site,
                    dimension,
                    0,
                    coordinate,
                )
                for coordinate in range(dimension)
            ),
        )
    return signature


@lru_cache(maxsize=None)
def canonical_marked_defect_atlas_signature(
    base: int,
    site: int,
    dimension: int,
    q: int,
    mode: str,
) -> tuple[object, ...]:
    """Canonicalize a marked defect-atlas signature."""

    return min(
        marked_defect_atlas_signature(
            permute_family_by_coordinate_order(
                base,
                dimension,
                order,
            ),
            permute_family_by_coordinate_order(
                site,
                dimension,
                order,
            ),
            dimension,
            order.index(q),
            mode,
        )
        for order in itertools.permutations(range(dimension))
    )


def is_rooted_q_activation_obstruction(
    base: int,
    site: int,
    dimension: int,
    root: int,
    q: int,
) -> bool:
    """Test rooted-minor minimal failure to activate ``q``.

    The simultaneous rooted minors of a contained pair are obtained by
    contracting a coordinate or restricting both members to the side
    containing the projected root.  Activation is permanent under either
    operation when the marked coordinate is retained.  It is therefore
    enough to test the two immediate minors in every coordinate other than
    ``q``.
    """

    if (
        q < 0
        or q >= dimension
        or site & ~base
        or root not in vertices(site, dimension)
    ):
        return False

    def activated(
        current_base: int,
        current_site: int,
        current_dimension: int,
        current_root: int,
        current_q: int,
    ) -> bool:
        closure_test = is_q_activated(
            current_base,
            current_site,
            current_dimension,
            current_root,
            current_q,
        )
        fiber_test = is_q_activated_by_q_free_fibers(
            current_base,
            current_site,
            current_dimension,
            current_root,
            current_q,
        )
        if closure_test != fiber_test:
            raise AssertionError(
                "the closure and q-free fiber activation tests disagree"
            )
        return closure_test

    if activated(base, site, dimension, root, q):
        return False

    for coordinate in range(dimension):
        if coordinate == q:
            continue
        root_value = (root >> coordinate) & 1
        projected_root = project_vertex(root, coordinate)
        projected_q = q - (q > coordinate)
        if not activated(
            contract_coordinate(base, dimension, coordinate),
            contract_coordinate(site, dimension, coordinate),
            dimension - 1,
            projected_root,
            projected_q,
        ):
            return False
        if not activated(
            section(base, dimension, coordinate, root_value),
            section(site, dimension, coordinate, root_value),
            dimension - 1,
            projected_root,
            projected_q,
        ):
            return False
    return True


def is_critical_reverse_pair(
    base: int,
    site: int,
    dimension: int,
    root: int,
) -> bool:
    """Test the exact pair conditions for a root-critical reverse lift."""

    if site & ~base or root not in vertices(site, dimension):
        return False
    if base == site:
        return False
    if not (
        is_rooted_r_per(base, dimension, root)
        and is_rooted_r_per(site, dimension, root)
    ):
        return False
    if set(
        raw_root_avoiding_peripheral_coordinates(
            base,
            dimension,
            root,
        )
    ) & set(
        raw_root_avoiding_peripheral_coordinates(
            site,
            dimension,
            root,
        )
    ):
        return False

    for coordinate in range(dimension):
        root_value = (root >> coordinate) & 1
        projected_root = project_vertex(root, coordinate)
        if not is_synchronizable_reverse_pair(
            section(base, dimension, coordinate, root_value),
            section(site, dimension, coordinate, root_value),
            dimension - 1,
            projected_root,
        ):
            return False
        if not is_synchronizable_reverse_pair(
            contract_coordinate(base, dimension, coordinate),
            contract_coordinate(site, dimension, coordinate),
            dimension - 1,
            projected_root,
        ):
            return False
    return True


def is_critical_reverse_pair_by_closure(
    base: int,
    site: int,
    dimension: int,
    root: int,
) -> bool:
    """Test root-criticality using deterministic common closures only."""

    if site & ~base or root not in vertices(site, dimension):
        return False
    if base == site:
        return False
    if not (
        is_rooted_r_per(base, dimension, root)
        and is_rooted_r_per(site, dimension, root)
    ):
        return False
    if set(
        raw_root_avoiding_peripheral_coordinates(
            base,
            dimension,
            root,
        )
    ) & set(
        raw_root_avoiding_peripheral_coordinates(
            site,
            dimension,
            root,
        )
    ):
        return False

    for coordinate in range(dimension):
        root_value = (root >> coordinate) & 1
        projected_root = project_vertex(root, coordinate)
        if not is_synchronizable_by_common_closure(
            section(base, dimension, coordinate, root_value),
            section(site, dimension, coordinate, root_value),
            dimension - 1,
            projected_root,
        ):
            return False
        if not is_synchronizable_by_common_closure(
            contract_coordinate(base, dimension, coordinate),
            contract_coordinate(site, dimension, coordinate),
            dimension - 1,
            projected_root,
        ):
            return False
    return True


def endpoint_rooted_three_vertex_site_coordinates(
    site: int,
    dimension: int,
    root: int,
) -> tuple[int, int] | None:
    """Recognize an endpoint-rooted three-vertex path.

    The returned pair ``(p,q)`` means that the path, starting at ``root``,
    first crosses coordinate ``p`` and then coordinate ``q``.  A path rooted
    at its middle vertex returns ``None``.
    """

    members = vertices(site, dimension)
    if len(members) != 3 or root not in members:
        return None
    root_neighbors = [
        vertex
        for vertex in members
        if vertex != root and (vertex ^ root).bit_count() == 1
    ]
    if len(root_neighbors) != 1:
        return None
    near = root_neighbors[0]
    far = next(
        vertex
        for vertex in members
        if vertex not in (root, near)
    )
    far_difference = far ^ root
    first_difference = near ^ root
    if (
        far_difference.bit_count() != 2
        or not (far_difference & first_difference)
        or (far ^ near).bit_count() != 1
    ):
        return None
    p = first_difference.bit_length() - 1
    q = (far_difference ^ first_difference).bit_length() - 1
    return p, q


def is_three_vertex_collar_critical_by_lower_activation(
    base: int,
    site: int,
    dimension: int,
    root: int,
) -> bool:
    """Test the exact lower-dimensional collar criterion.

    For an endpoint-rooted path site, let ``p`` be its first coordinate and
    ``q`` its second.  Splitting the base along ``p`` gives a rooted-side
    family ``C`` and an opposite family ``T <= C``.  The upper pair is
    critical exactly when ``C`` is rooted-good, ``T`` is ``R_per``, the
    lower pair has disjoint initial availability, and it is a rooted
    ``q``-activation obstruction.  Rooted-goodness of every root-side
    restriction of ``C`` follows automatically from rooted-minor permanence.
    """

    path_coordinates = endpoint_rooted_three_vertex_site_coordinates(
        site,
        dimension,
        root,
    )
    if (
        path_coordinates is None
        or site & ~base
        or base == site
    ):
        return False
    p, q = path_coordinates
    root_value = (root >> p) & 1
    lower_root = project_vertex(root, p)
    lower_dimension = dimension - 1
    lower_q = q - (q > p)
    rooted_side = section(
        base,
        dimension,
        p,
        root_value,
    )
    opposite_side = section(
        base,
        dimension,
        p,
        1 - root_value,
    )
    if not (
        nested(opposite_side, rooted_side)
        and is_rooted_r_per(
            rooted_side,
            lower_dimension,
            lower_root,
        )
        and is_r_per(opposite_side, lower_dimension)
    ):
        return False
    if set(
        raw_root_avoiding_peripheral_coordinates(
            rooted_side,
            lower_dimension,
            lower_root,
        )
    ) & set(
        raw_root_avoiding_peripheral_coordinates(
            opposite_side,
            lower_dimension,
            lower_root,
        )
    ):
        return False

    return is_rooted_q_activation_obstruction(
        rooted_side,
        opposite_side,
        lower_dimension,
        lower_root,
        lower_q,
    )


def successful_synchronizing_moves(
    base: int,
    site: int,
    dimension: int,
    root: int,
) -> tuple[int, bool]:
    """Return successful first old moves and whether equality terminates."""

    terminal = (
        base == site
        and is_rooted_r_per(base, dimension, root)
    )
    common_available = (
        set(
            raw_root_avoiding_peripheral_coordinates(
                base,
                dimension,
                root,
            )
        )
        & set(
            raw_root_avoiding_peripheral_coordinates(
                site,
                dimension,
                root,
            )
        )
    )
    successful_mask = 0
    for coordinate in common_available:
        if is_synchronizable_reverse_pair(
            contract_coordinate(base, dimension, coordinate),
            contract_coordinate(site, dimension, coordinate),
            dimension - 1,
            project_vertex(root, coordinate),
        ):
            successful_mask |= 1 << coordinate
    return successful_mask, terminal


def permute_family_by_coordinate_order(
    family: int,
    dimension: int,
    order: tuple[int, ...],
) -> int:
    """Relabel old coordinates according to ``new -> old`` order."""

    return family_from_vertices(
        sum(
            ((vertex >> old_coordinate) & 1) << new_coordinate
            for new_coordinate, old_coordinate in enumerate(order)
        )
        for vertex in vertices(family, dimension)
    )


def permute_coordinate_mask(
    mask: int,
    order: tuple[int, ...],
) -> int:
    """Relabel a coordinate subset according to ``new -> old`` order."""

    return sum(
        ((mask >> old_coordinate) & 1) << new_coordinate
        for new_coordinate, old_coordinate in enumerate(order)
    )


def lift_coordinate_mask(
    mask: int,
    dimension: int,
    omitted_coordinate: int,
) -> int:
    """Lift a mask from ``E minus omitted_coordinate`` back to ``E``."""

    old_labels = tuple(
        coordinate
        for coordinate in range(dimension)
        if coordinate != omitted_coordinate
    )
    return sum(
        ((mask >> new_coordinate) & 1) << old_coordinate
        for new_coordinate, old_coordinate in enumerate(old_labels)
    )


@lru_cache(maxsize=None)
def canonical_reverse_pair(
    base: int,
    site: int,
    dimension: int,
) -> tuple[int, int]:
    """Canonicalize a rooted-at-zero contained pair under permutations."""

    return min(
        (
            permute_family_by_coordinate_order(
                base,
                dimension,
                order,
            ),
            permute_family_by_coordinate_order(
                site,
                dimension,
                order,
            ),
        )
        for order in itertools.permutations(range(dimension))
    )


@lru_cache(maxsize=None)
def canonical_marked_reverse_pair(
    base: int,
    site: int,
    dimension: int,
    q: int,
) -> tuple[int, int, int]:
    """Canonicalize a rooted-at-zero contained pair with marked coordinate."""

    if q < 0 or q >= dimension:
        raise ValueError("the marked coordinate is outside the cube")
    return min(
        (
            permute_family_by_coordinate_order(
                base,
                dimension,
                order,
            ),
            permute_family_by_coordinate_order(
                site,
                dimension,
                order,
            ),
            order.index(q),
        )
        for order in itertools.permutations(range(dimension))
    )


@lru_cache(maxsize=None)
def canonical_marked_pair_extension(
    base: int,
    site: int,
    extended_site: int,
    dimension: int,
    q: int,
) -> tuple[int, int, int, int]:
    """Canonicalize a marked pair together with one larger site."""

    if site & ~extended_site or extended_site & ~base:
        raise ValueError(
            "the marked extension must satisfy site <= extended site <= base"
        )
    if q < 0 or q >= dimension:
        raise ValueError("the marked coordinate is outside the cube")
    return min(
        (
            permute_family_by_coordinate_order(
                base,
                dimension,
                order,
            ),
            permute_family_by_coordinate_order(
                site,
                dimension,
                order,
            ),
            permute_family_by_coordinate_order(
                extended_site,
                dimension,
                order,
            ),
            order.index(q),
        )
        for order in itertools.permutations(range(dimension))
    )


@lru_cache(maxsize=None)
def canonical_activation_antimatroid_signature(
    base: int,
    site: int,
    dimension: int,
) -> tuple[
    tuple[int, ...],
    tuple[tuple[tuple[int, ...], tuple[int, ...]], ...],
]:
    """Canonical feasible-set antimatroids for a pair and its forced minors."""

    initial_feasible = common_availability_antimatroid(
        base,
        site,
        dimension,
        0,
    )
    rows: list[tuple[tuple[int, ...], tuple[int, ...]]] = []
    for forced_coordinate in range(dimension):
        contraction_feasible = tuple(
            sorted(
                (1 << forced_coordinate)
                | lift_coordinate_mask(
                    mask,
                    dimension,
                    forced_coordinate,
                )
                for mask in common_availability_antimatroid(
                    contract_coordinate(
                        base,
                        dimension,
                        forced_coordinate,
                    ),
                    contract_coordinate(
                        site,
                        dimension,
                        forced_coordinate,
                    ),
                    dimension - 1,
                    0,
                )
            )
        )
        restriction_feasible = tuple(
            sorted(
                lift_coordinate_mask(
                    mask,
                    dimension,
                    forced_coordinate,
                )
                for mask in common_availability_antimatroid(
                    section(base, dimension, forced_coordinate, 0),
                    section(site, dimension, forced_coordinate, 0),
                    dimension - 1,
                    0,
                )
            )
        )
        rows.append((contraction_feasible, restriction_feasible))

    return min(
        (
            tuple(
                sorted(
                    permute_coordinate_mask(mask, order)
                    for mask in initial_feasible
                )
            ),
            tuple(
                (
                    tuple(
                        sorted(
                            permute_coordinate_mask(mask, order)
                            for mask in rows[old_coordinate][0]
                        )
                    ),
                    tuple(
                        sorted(
                            permute_coordinate_mask(mask, order)
                            for mask in rows[old_coordinate][1]
                        )
                    ),
                )
                for old_coordinate in order
            ),
        )
        for order in itertools.permutations(range(dimension))
    )


@lru_cache(maxsize=None)
def canonical_repair_signature(
    base: int,
    site: int,
    dimension: int,
) -> tuple[int, int, tuple[tuple[int, bool, int, bool], ...]]:
    """Canonical first-successful-move data for a critical reverse pair."""

    rows: list[tuple[int, bool, int, bool]] = []
    for forced_coordinate in range(dimension):
        old_labels = tuple(
            coordinate
            for coordinate in range(dimension)
            if coordinate != forced_coordinate
        )

        def lift_mask(mask: int) -> int:
            return sum(
                ((mask >> new_coordinate) & 1) << old_coordinate
                for new_coordinate, old_coordinate in enumerate(
                    old_labels
                )
            )

        contraction_moves, contraction_terminal = (
            successful_synchronizing_moves(
                contract_coordinate(
                    base,
                    dimension,
                    forced_coordinate,
                ),
                contract_coordinate(
                    site,
                    dimension,
                    forced_coordinate,
                ),
                dimension - 1,
                0,
            )
        )
        root_base = section(
            base,
            dimension,
            forced_coordinate,
            0,
        )
        root_site = section(
            site,
            dimension,
            forced_coordinate,
            0,
        )
        restriction_moves, restriction_terminal = (
            successful_synchronizing_moves(
                root_base,
                root_site,
                dimension - 1,
                0,
            )
        )
        rows.append(
            (
                lift_mask(contraction_moves),
                contraction_terminal,
                lift_mask(restriction_moves),
                restriction_terminal,
            )
        )

    base_available = sum(
        1 << coordinate
        for coordinate in raw_root_avoiding_peripheral_coordinates(
            base,
            dimension,
            0,
        )
    )
    site_available = sum(
        1 << coordinate
        for coordinate in raw_root_avoiding_peripheral_coordinates(
            site,
            dimension,
            0,
        )
    )
    return min(
        (
            permute_coordinate_mask(base_available, order),
            permute_coordinate_mask(site_available, order),
            tuple(
                (
                    permute_coordinate_mask(rows[old_coordinate][0], order),
                    rows[old_coordinate][1],
                    permute_coordinate_mask(rows[old_coordinate][2], order),
                    rows[old_coordinate][3],
                )
                for old_coordinate in order
            ),
        )
        for order in itertools.permutations(range(dimension))
    )


def is_root_locked(family: int, dimension: int, root: int) -> bool:
    """Whether a nontrivial class member has no periphery avoiding its root."""

    family, dimension, root = irredundant_rooted_projection(
        family, dimension, root
    )
    return (
        family.bit_count() > 1
        and is_r_per(family, dimension)
        and not root_avoiding_peripheral_coordinates(
            family, dimension, root
        )
    )


def is_root_critical(family: int, dimension: int, root: int) -> bool:
    """Whether the pointed factor yields minimal obstructions under wedging.

    Besides being root-locked, every root-side restriction and every
    contraction must have a root-avoiding recursive certificate.
    """

    family, dimension, root = irredundant_rooted_projection(
        family, dimension, root
    )
    if not is_root_locked(family, dimension, root):
        return False
    for coordinate in range(dimension):
        root_value = (root >> coordinate) & 1
        root_side = section(
            family, dimension, coordinate, root_value
        )
        contraction = (
            section(family, dimension, coordinate, 0)
            | section(family, dimension, coordinate, 1)
        )
        projected_root = project_vertex(root, coordinate)
        if not is_rooted_r_per(
            root_side, dimension - 1, projected_root
        ):
            return False
        if not is_rooted_r_per(
            contraction, dimension - 1, projected_root
        ):
            return False
    return True


def canonical_pointed_family(
    family: int, dimension: int, root: int
) -> int:
    """Canonical mask of a pointed family, with its root transported to zero."""

    family, dimension, root = irredundant_rooted_projection(
        family, dimension, root
    )
    translated = tuple(
        vertex ^ root for vertex in vertices(family, dimension)
    )
    return min(
        family_from_vertices(vertex_map[vertex] for vertex in translated)
        for vertex_map in coordinate_permutation_maps(dimension)
    )


def one_vertex_wedge(
    left: int,
    left_dimension: int,
    left_root: int,
    right: int,
    right_dimension: int,
    right_root: int,
) -> tuple[int, int]:
    """Embed two pointed families on disjoint blocks and identify their roots."""

    left, left_dimension, left_root = irredundant_rooted_projection(
        left, left_dimension, left_root
    )
    right, right_dimension, right_root = irredundant_rooted_projection(
        right, right_dimension, right_root
    )
    left_members = {
        vertex ^ left_root for vertex in vertices(left, left_dimension)
    }
    right_members = {
        (vertex ^ right_root) << left_dimension
        for vertex in vertices(right, right_dimension)
    }
    return (
        family_from_vertices(left_members | right_members),
        left_dimension + right_dimension,
    )


def root_locked_factor() -> int:
    """The rooted 11-vertex factor used in the gated-amalgam obstruction."""

    return family_from_vertices(
        (0, 4, 5, 7, 8, 10, 11, 12, 13, 14, 15)
    )


def rooted_wedge_obstruction() -> int:
    """Two root-locked factors wedged at zero in disjoint coordinates."""

    factor_members = set(vertices(root_locked_factor(), 4))
    return family_from_vertices(
        factor_members | {vertex << 4 for vertex in factor_members}
    )


def immediate_minor_isomorphism_groups(
    family: int, dimension: int
) -> list[dict[str, object]]:
    groups: dict[tuple[int, int], list[str]] = {}
    for minor, minor_dimension, operation in immediate_pc_minors_with_operations(
        family, dimension
    ):
        canonical, _ = cube_symmetry_profile(minor, minor_dimension)
        groups.setdefault((canonical, minor_dimension), []).append(operation)
    answer: list[dict[str, object]] = []
    for (canonical, minor_dimension), operations in sorted(
        groups.items(), key=lambda item: (item[0][0].bit_count(), item[0][0])
    ):
        answer.append(
            {
                "order": canonical.bit_count(),
                "operations": tuple(operations),
                "vertices": bitstrings(canonical, minor_dimension),
                "R_certificate": r_per_certificate(
                    canonical, minor_dimension
                ),
            }
        )
    return answer


def immediate_pc_minors_with_operations(
    family: int, dimension: int
) -> list[tuple[int, int, str]]:
    """Return distinct normalized elementary pc-minors and one operation."""

    seen: set[tuple[int, int]] = set()
    answer: list[tuple[int, int, str]] = []
    for coordinate in range(dimension):
        zero = section(family, dimension, coordinate, 0)
        one = section(family, dimension, coordinate, 1)
        for name, minor in (("r0", zero), ("r1", one), ("c", zero | one)):
            normalized = irredundant_projection(minor, dimension - 1)
            if normalized not in seen:
                seen.add(normalized)
                answer.append((*normalized, f"{name}@{coordinate}"))
    return answer


@lru_cache(maxsize=100_000)
def is_partial_cube_family(family: int, dimension: int) -> bool:
    """Check equality of graph and Hamming distances for every vertex pair."""

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


def shattering_profile(family: int, dimension: int) -> tuple[int, int, bool]:
    """Return numbers shattered/strongly shattered and their exact equality."""

    members = vertices(family, dimension)
    shattered: set[int] = set()
    strongly_shattered: set[int] = set()
    universe_mask = (1 << dimension) - 1
    for coordinates in range(1 << dimension):
        cube_order = 1 << coordinates.bit_count()
        if len({vertex & coordinates for vertex in members}) == cube_order:
            shattered.add(coordinates)
        complement = universe_mask ^ coordinates
        fiber_sizes = Counter(vertex & complement for vertex in members)
        if fiber_sizes and max(fiber_sizes.values()) == cube_order:
            strongly_shattered.add(coordinates)
    return (
        len(shattered),
        len(strongly_shattered),
        shattered == strongly_shattered,
    )


def shattered_rank_histogram(
    family: int,
    dimension: int,
) -> dict[int, int]:
    """Count shattered coordinate sets by cardinality."""
    members = vertices(family, dimension)
    histogram: Counter[int] = Counter()
    for coordinates in range(1 << dimension):
        cube_order = 1 << coordinates.bit_count()
        if len({vertex & coordinates for vertex in members}) == cube_order:
            histogram[coordinates.bit_count()] += 1
    return dict(sorted(histogram.items()))


def is_ample_family(family: int, dimension: int) -> bool:
    """Test ampleness via the exact Sandwich-equality criterion.

    A family is shattering-extremal precisely when the number of shattered
    coordinate sets equals its cardinality.  Counting only shattered sets is
    substantially cheaper than materializing both shattering complexes.
    """

    members = vertices(family, dimension)
    shattered_count = 0
    for coordinates in range(1 << dimension):
        cube_order = 1 << coordinates.bit_count()
        if len({vertex & coordinates for vertex in members}) == cube_order:
            shattered_count += 1
            if shattered_count > len(members):
                return False
    return shattered_count == len(members)


@lru_cache(maxsize=50_000)
def shattered_coordinate_sets_mask(family: int, dimension: int) -> int:
    """Encode the shattered coordinate sets as a bitset."""

    members = vertices(family, dimension)
    shattered = 0
    for coordinates in range(1 << dimension):
        if len({vertex & coordinates for vertex in members}) == (
            1 << coordinates.bit_count()
        ):
            shattered |= 1 << coordinates
    return shattered


def ample_expansion_of_ample_base(
    zero: int,
    one: int,
    base_dimension: int,
) -> bool:
    """Test ampleness of an expansion whose contraction is already ample.

    If ``C`` is the lift of sections ``zero`` and ``one`` and
    ``B = zero | one`` is ample, then the shattered sets of ``C`` split into
    those not using the new coordinate (the shattered sets of ``B``) and
    those using it (the sets shattered by both sections).  Since
    ``|C| = |B| + |zero & one|``, Sandwich equality gives the criterion below.
    """

    common_shattered = (
        shattered_coordinate_sets_mask(zero, base_dimension)
        & shattered_coordinate_sets_mask(one, base_dimension)
    )
    return common_shattered.bit_count() == (zero & one).bit_count()


def ample_predicate_regression(
    max_base_dimension: int = 4,
    order_limit: int = 10,
) -> None:
    """Compare Sandwich-count and full-complex ampleness predicates."""

    census = enumerate_partial_cubes(max_base_dimension)
    tested: set[tuple[int, int]] = set()
    for dimension in range(max_base_dimension + 1):
        tested.update((family, dimension) for family in census[dimension])
    for base_dimension in range(1, max_base_dimension + 1):
        tested.update(
            (family, base_dimension + 1)
            for family in iter_expansion_candidates(
                census[base_dimension],
                base_dimension,
                order_limit,
            )
        )
    failures = [
        (family, dimension)
        for family, dimension in tested
        if is_ample_family(family, dimension)
        != shattering_profile(family, dimension)[2]
    ]
    expansion_tested = 0
    expansion_failures: list[tuple[int, int]] = []
    for family, dimension in tested:
        if dimension == 0:
            continue
        zero = section(family, dimension, dimension - 1, 0)
        one = section(family, dimension, dimension - 1, 1)
        if not zero or not one:
            continue
        base_dimension = dimension - 1
        if not is_ample_family(zero | one, base_dimension):
            continue
        expansion_tested += 1
        if ample_expansion_of_ample_base(
            zero, one, base_dimension
        ) != is_ample_family(family, dimension):
            expansion_failures.append((family, dimension))
    print(
        "\nAmple predicate regression\n"
        f"  max_base_dimension={max_base_dimension}; "
        f"order_limit={order_limit}; tested={len(tested)}; "
        f"failures={len(failures)}; "
        f"ample-base-expansions={expansion_tested}; "
        f"expansion_failures={len(expansion_failures)}",
        flush=True,
    )
    if failures:
        family, dimension = failures[0]
        raise AssertionError(
            "Sandwich-count ampleness mismatch: "
            f"dimension={dimension}, vertices={bitstrings(family, dimension)}"
        )
    if expansion_failures:
        family, dimension = expansion_failures[0]
        raise AssertionError(
            "ample-base expansion criterion mismatch: "
            f"dimension={dimension}, vertices={bitstrings(family, dimension)}"
        )


def shattered_size_histogram(family: int, dimension: int) -> tuple[int, ...]:
    members = vertices(family, dimension)
    histogram = [0] * (dimension + 1)
    for coordinates in range(1 << dimension):
        if len({vertex & coordinates for vertex in members}) == (
            1 << coordinates.bit_count()
        ):
            histogram[coordinates.bit_count()] += 1
    return tuple(histogram)


def crossing_pairs(
    family: int, dimension: int
) -> tuple[tuple[int, int], ...]:
    """Coordinate pairs shattered by the family.

    For an ample partial cube these are exactly the pairs of crossing
    ``Theta``-classes, since shattered and strongly shattered sets agree.
    """

    members = vertices(family, dimension)
    answer: list[tuple[int, int]] = []
    for left in range(dimension):
        for right in range(left + 1, dimension):
            coordinates = (1 << left) | (1 << right)
            if len({vertex & coordinates for vertex in members}) == 4:
                answer.append((left, right))
    return tuple(answer)


def noncrossing_pairs(
    family: int, dimension: int
) -> tuple[tuple[int, int], ...]:
    crossed = set(crossing_pairs(family, dimension))
    return tuple(
        (left, right)
        for left in range(dimension)
        for right in range(left + 1, dimension)
        if (left, right) not in crossed
    )


def periphery_free_semicube_witnesses(
    family: int, dimension: int, left: int, right: int
) -> tuple[str, ...]:
    """Semicube restrictions at a pair that remain periphery-free."""

    answer: list[str] = []
    for coordinate in (left, right):
        for value in (0, 1):
            restricted, lower_dimension = irredundant_projection(
                section(family, dimension, coordinate, value),
                dimension - 1,
            )
            if (
                restricted.bit_count() > 1
                and not peripheral_decompositions(
                    restricted, lower_dimension
                )
            ):
                answer.append(f"r{value}@{coordinate}")
    return tuple(answer)


def maximal_cubes(family: int, dimension: int) -> tuple[frozenset[int], ...]:
    """Enumerate inclusion-maximal cube faces contained in the family."""

    members = set(vertices(family, dimension))
    cubes: set[frozenset[int]] = set()
    full = (1 << dimension) - 1
    for support in range(1 << dimension):
        complement = full ^ support
        anchor = complement
        subanchor = anchor
        while True:
            face = frozenset(
                subanchor | varying
                for varying in range(1 << dimension)
                if varying & ~support == 0
            )
            if face <= members:
                cubes.add(face)
            if subanchor == 0:
                break
            subanchor = (subanchor - 1) & anchor
    return tuple(
        cube
        for cube in cubes
        if not any(cube < other for other in cubes)
    )


@lru_cache(maxsize=None)
def ccmw_corners(family: int, dimension: int) -> tuple[int, ...]:
    cubes = maximal_cubes(family, dimension)
    return tuple(
        vertex
        for vertex in vertices(family, dimension)
        if sum(vertex in cube for cube in cubes) == 1
    )


@lru_cache(maxsize=None)
def has_corner_peeling(family: int, dimension: int) -> bool:
    family, dimension = irredundant_projection(family, dimension)
    if family.bit_count() == 1:
        return True
    for corner in ccmw_corners(family, dimension):
        residual = family & ~(1 << corner)
        residual, residual_dimension = irredundant_projection(
            residual, dimension
        )
        if shattering_profile(residual, residual_dimension)[2] and has_corner_peeling(
            residual, residual_dimension
        ):
            return True
    return False


@lru_cache(maxsize=None)
def has_corner_peeling_ending_at_root(
    family: int,
    dimension: int,
    root: int,
) -> bool:
    """Whether a corner deletion sequence can leave ``root`` until last."""
    family, dimension, root = irredundant_rooted_projection(
        family,
        dimension,
        root,
    )
    if family.bit_count() == 1:
        return True
    for corner in ccmw_corners(family, dimension):
        if corner == root:
            continue
        residual = family & ~(1 << corner)
        (
            residual,
            residual_dimension,
            residual_root,
        ) = irredundant_rooted_projection(
            residual,
            dimension,
            root,
        )
        if (
            shattering_profile(
                residual,
                residual_dimension,
            )[2]
            and has_corner_peeling_ending_at_root(
                residual,
                residual_dimension,
                residual_root,
            )
        ):
            return True
    return False


@lru_cache(maxsize=None)
def corner_peeling_order_ending_at_root(
    family: int,
    dimension: int,
    root: int,
) -> tuple[int, ...] | None:
    """Return an original-coordinate corner deletion order ending at root."""
    if root not in vertices(family, dimension):
        raise ValueError("the root is not a member of the family")
    if family.bit_count() == 1:
        return (root,)
    for corner in ccmw_corners(family, dimension):
        if corner == root:
            continue
        residual = family & ~(1 << corner)
        if not shattering_profile(residual, dimension)[2]:
            continue
        tail = corner_peeling_order_ending_at_root(
            residual,
            dimension,
            root,
        )
        if tail is not None:
            return (corner,) + tail
    return None


def cube_face_pattern(
    cube: frozenset[int],
    dimension: int,
) -> str:
    """Encode a cube face by ``0``, ``1``, and ``*`` coordinates."""
    pattern: list[str] = []
    for coordinate in reversed(range(dimension)):
        values = {
            (vertex >> coordinate) & 1
            for vertex in cube
        }
        pattern.append(
            "*" if len(values) == 2 else str(next(iter(values)))
        )
    return "".join(pattern)


def corner_peeling_face_certificate(
    family: int,
    dimension: int,
    deletion_order: tuple[int, ...],
) -> tuple[tuple[str, str], ...]:
    """Record the unique maximal cube witnessing each corner deletion."""
    residual = family
    certificate: list[tuple[str, str]] = []
    for corner in deletion_order[:-1]:
        containing = [
            cube
            for cube in maximal_cubes(residual, dimension)
            if corner in cube
        ]
        if len(containing) != 1:
            raise AssertionError(
                "the proposed corner belongs to "
                f"{len(containing)} maximal cubes"
            )
        certificate.append(
            (
                format(corner, f"0{dimension}b"),
                cube_face_pattern(containing[0], dimension),
            )
        )
        residual &= ~(1 << corner)
    if residual != 1 << deletion_order[-1]:
        raise AssertionError("the corner sequence does not end at one vertex")
    return tuple(certificate)


def corner_peeling_certificate(
    family: int, dimension: int
) -> tuple[tuple[int, str], ...] | None:
    """Return one corner deletion sequence in the current coordinates."""

    family, dimension = irredundant_projection(family, dimension)
    if family.bit_count() == 1:
        return ((dimension, bitstrings(family, dimension)),)
    for corner in ccmw_corners(family, dimension):
        residual = family & ~(1 << corner)
        residual, residual_dimension = irredundant_projection(
            residual, dimension
        )
        if not shattering_profile(residual, residual_dimension)[2]:
            continue
        tail = corner_peeling_certificate(residual, residual_dimension)
        if tail is not None:
            return ((dimension, format(corner, f"0{dimension}b")),) + tail
    return None


def transform_family(
    family: int,
    dimension: int,
    permutation: tuple[int, ...],
    flip: int,
) -> int:
    transformed: list[int] = []
    for vertex in vertices(family, dimension):
        image = 0
        for new_coordinate, old_coordinate in enumerate(permutation):
            bit = ((vertex >> old_coordinate) & 1) ^ (
                (flip >> old_coordinate) & 1
            )
            image |= bit << new_coordinate
        transformed.append(image)
    return family_from_vertices(transformed)


@lru_cache(maxsize=None)
def coordinate_permutation_maps(dimension: int) -> tuple[tuple[int, ...], ...]:
    """All vertex maps induced by coordinate permutations (without flips)."""

    maps: list[tuple[int, ...]] = []
    for permutation in itertools.permutations(range(dimension)):
        images: list[int] = []
        for vertex in range(1 << dimension):
            image = 0
            for new_coordinate, old_coordinate in enumerate(permutation):
                image |= ((vertex >> old_coordinate) & 1) << new_coordinate
            images.append(image)
        maps.append(tuple(images))
    return tuple(maps)


@lru_cache(maxsize=None)
def canonical_cube_family(family: int, dimension: int) -> int:
    """Canonicalize under coordinate permutations and flips.

    A representative containing zero can be obtained by translating any one
    member to zero, so only flips indexed by members need be considered.
    """

    members = vertices(family, dimension)
    best: int | None = None
    for origin in members:
        translated = tuple(vertex ^ origin for vertex in members)
        for vertex_map in coordinate_permutation_maps(dimension):
            transformed = family_from_vertices(
                vertex_map[vertex] for vertex in translated
            )
            if best is None or transformed < best:
                best = transformed
    assert best is not None
    return best


def tree_canonical_cube_family(family: int, dimension: int) -> int:
    """Canonical cube representative of a tree partial cube.

    In a tree every edge is its own Theta-class, so graph isomorphisms and
    cube isomorphisms agree.  Rooted-tree codes canonically order the branches;
    assigning one coordinate to each edge in that order then gives a cube
    embedding depending only on the unrooted tree.  At a bicentral symmetry
    the two center-rooted embeddings are compared explicitly.
    """

    members = vertices(family, dimension)
    member_set = set(members)
    adjacency = {
        vertex: tuple(
            vertex ^ (1 << coordinate)
            for coordinate in range(dimension)
            if vertex ^ (1 << coordinate) in member_set
        )
        for vertex in members
    }
    edge_count = sum(map(len, adjacency.values())) // 2
    if edge_count != len(members) - 1:
        raise ValueError("the supplied family is not a tree")

    remaining = set(members)
    while len(remaining) > 2:
        leaves = {
            vertex
            for vertex in remaining
            if sum(neighbor in remaining for neighbor in adjacency[vertex])
            <= 1
        }
        if not leaves:
            raise ValueError("the supplied family is not a tree")
        remaining -= leaves
    centers = tuple(remaining)

    @lru_cache(maxsize=None)
    def rooted_embedding(
        vertex: int, parent: int
    ) -> tuple[str, tuple[int, ...], int]:
        children = [
            rooted_embedding(neighbor, vertex)
            for neighbor in adjacency[vertex]
            if neighbor != parent
        ]
        children.sort(key=lambda item: (item[0], item[2], item[1]))
        code = "(" + "".join(child[0] for child in children) + ")"
        embedded = {0}
        offset = 0
        for _, child_members, child_dimension in children:
            connecting_bit = 1 << offset
            embedded.update(
                connecting_bit | (member << (offset + 1))
                for member in child_members
            )
            offset += child_dimension + 1
        return code, tuple(sorted(embedded)), offset

    rooted_candidates: list[tuple[str, int]] = []
    for center in centers:
        code, embedded_members, embedded_dimension = rooted_embedding(
            center, -1
        )
        if embedded_dimension != dimension:
            raise ValueError("tree dimension does not equal its edge count")
        rooted_candidates.append(
            (code, family_from_vertices(embedded_members))
        )
    least_code = min(code for code, _ in rooted_candidates)
    return min(
        embedded
        for code, embedded in rooted_candidates
        if code == least_code
    )


def tree_canonicalizer_regression(max_dimension: int = 4) -> None:
    """Exhaustively test tree canonicalization under small cube symmetries."""

    census = enumerate_partial_cubes(max_dimension)
    tree_types = 0
    transformed_families = 0
    for dimension in range(max_dimension + 1):
        for family in census[dimension]:
            members = vertices(family, dimension)
            member_set = set(members)
            edge_count = sum(
                (vertex ^ (1 << coordinate)) in member_set
                for vertex in members
                for coordinate in range(dimension)
            ) // 2
            if edge_count != len(members) - 1:
                continue
            tree_types += 1
            expected = tree_canonical_cube_family(family, dimension)
            if canonical_cube_family(expected, dimension) != (
                canonical_cube_family(family, dimension)
            ):
                raise AssertionError(
                    "tree canonicalizer changed the cube-isomorphism type"
                )
            for permutation in itertools.permutations(range(dimension)):
                for flip in range(1 << dimension):
                    transformed = transform_family(
                        family, dimension, permutation, flip
                    )
                    transformed_families += 1
                    if (
                        tree_canonical_cube_family(
                            transformed, dimension
                        )
                        != expected
                    ):
                        raise AssertionError(
                            "tree canonicalizer is not cube-invariant"
                        )
    print(
        "\nTree canonicalizer regression\n"
        f"  max_dimension={max_dimension}; tree_types={tree_types}; "
        f"transformed_families={transformed_families}; failures=0",
        flush=True,
    )


def refined_canonicalizer_regression(max_dimension: int = 4) -> None:
    """Exhaustively test refined canonicalization under cube symmetries."""

    census = enumerate_partial_cubes(max_dimension)
    family_types = 0
    transformed_families = 0
    for dimension in range(max_dimension + 1):
        for family in census[dimension]:
            family_types += 1
            expected = refined_canonical_cube_family(family, dimension)
            if canonical_cube_family(expected, dimension) != (
                canonical_cube_family(family, dimension)
            ):
                raise AssertionError(
                    "refined canonicalizer changed the cube-isomorphism type"
                )
            for permutation in itertools.permutations(range(dimension)):
                for flip in range(1 << dimension):
                    transformed = transform_family(
                        family, dimension, permutation, flip
                    )
                    transformed_families += 1
                    if (
                        refined_canonical_cube_family(
                            transformed, dimension
                        )
                        != expected
                    ):
                        raise AssertionError(
                            "refined canonicalizer is not cube-invariant: "
                            f"dimension={dimension}, family={family:x}, "
                            f"permutation={permutation}, flip={flip:x}, "
                            f"expected={expected:x}, "
                            "actual="
                            f"{refined_canonical_cube_family(transformed, dimension):x}"
                        )
    print(
        "\nRefined canonicalizer regression\n"
        f"  max_dimension={max_dimension}; family_types={family_types}; "
        f"transformed_families={transformed_families}; failures=0",
        flush=True,
    )


def nauty_canonical_cube_family(family: int, dimension: int) -> int:
    """Canonicalize a cube family through a colored incidence gadget.

    Each coordinate is represented by a center adjacent to two indistinguish-
    able literal vertices.  Every concept is adjacent to exactly one literal
    from each pair.  A color-preserving graph isomorphism may permute the
    coordinate centers and swap the two literals of each pair, and these are
    exactly the coordinate permutations and flips of the ambient cube.
    """

    if pynauty is None:
        raise RuntimeError("pynauty is not available")
    members = vertices(family, dimension)
    concept_count = len(members)
    center_start = concept_count
    zero_literal_start = center_start + dimension
    one_literal_start = zero_literal_start + dimension
    graph_order = concept_count + 3 * dimension
    adjacency = {vertex: set() for vertex in range(graph_order)}

    def connect(left: int, right: int) -> None:
        adjacency[left].add(right)
        adjacency[right].add(left)

    for coordinate in range(dimension):
        center = center_start + coordinate
        zero_literal = zero_literal_start + coordinate
        one_literal = one_literal_start + coordinate
        connect(center, zero_literal)
        connect(center, one_literal)
        for concept, member in enumerate(members):
            connect(
                concept,
                one_literal if member & (1 << coordinate) else zero_literal,
            )

    graph = pynauty.Graph(
        number_of_vertices=graph_order,
        directed=False,
        adjacency_dict={
            vertex: sorted(neighbors)
            for vertex, neighbors in adjacency.items()
        },
        vertex_coloring=[
            set(range(concept_count)),
            set(range(center_start, zero_literal_start)),
            set(range(zero_literal_start, graph_order)),
        ],
    )
    canonical_order = pynauty.canon_label(graph)
    canonical_labels = [0] * graph_order
    for new_label, old_vertex in enumerate(canonical_order):
        canonical_labels[old_vertex] = new_label
    permutation = tuple(
        sorted(
            range(dimension),
            key=lambda coordinate: canonical_labels[
                center_start + coordinate
            ],
        )
    )
    flip = 0
    for coordinate in range(dimension):
        if canonical_labels[zero_literal_start + coordinate] > (
            canonical_labels[one_literal_start + coordinate]
        ):
            flip |= 1 << coordinate
    return transform_family(
        family, dimension, permutation, flip
    )


def refined_canonical_cube_family(family: int, dimension: int) -> int:
    """Canonicalize exactly, after invariant root/coordinate refinement.

    Cube isomorphisms preserve Hamming distances and one-inclusion degrees.
    We may therefore restrict possible images of zero to vertices having the
    lexicographically least rooted distance/degree profile.  After translating
    such a vertex to zero, a root-fixing cube isomorphism is a coordinate
    permutation.  The rooted vertex--coordinate incidence refinement below is
    also invariant, so coordinates with distinct stable colors have a forced
    common order; only permutations inside tied stable-color blocks need be
    enumerated.

    This gives an exact canonical form, not a hash: every retained transform is
    a cube automorphism, and every cube isomorphism maps a selected root and
    each signature block to their selected counterparts.
    """

    members = vertices(family, dimension)
    member_set = set(members)
    edge_count = sum(
        (vertex ^ (1 << coordinate)) in member_set
        for vertex in members
        for coordinate in range(dimension)
    ) // 2
    if edge_count == len(members) - 1:
        return tree_canonical_cube_family(family, dimension)
    if pynauty is not None:
        return nauty_canonical_cube_family(family, dimension)

    degrees = {
        vertex: sum(
            (vertex ^ (1 << coordinate)) in member_set
            for coordinate in range(dimension)
        )
        for vertex in members
    }

    def rooted_profile(origin: int) -> tuple[tuple[tuple[int, int], int], ...]:
        return tuple(
            sorted(
                Counter(
                    ((vertex ^ origin).bit_count(), degrees[vertex])
                    for vertex in members
                ).items()
            )
        )

    profiles = {origin: rooted_profile(origin) for origin in members}
    least_profile = min(profiles.values())
    candidate_origins = tuple(
        origin
        for origin in members
        if profiles[origin] == least_profile
    )

    best: int | None = None
    for origin in candidate_origins:
        translated_members = tuple(vertex ^ origin for vertex in members)
        translated_set = set(translated_members)
        translated_degrees = {
            vertex ^ origin: degree for vertex, degree in degrees.items()
        }
        vertex_colors: dict[int, object] = {
            vertex: (vertex.bit_count(), translated_degrees[vertex])
            for vertex in translated_members
        }
        coordinate_colors: dict[int, object] = {}
        for coordinate in range(dimension):
            bit = 1 << coordinate
            coordinate_colors[coordinate] = (
                sum(bool(vertex & bit) for vertex in translated_members),
                sum(
                    (vertex ^ bit) in translated_set
                    for vertex in translated_members
                    if not (vertex & bit)
                ),
            )

        while True:
            old_vertex_partition = {
                frozenset(
                    vertex
                    for vertex, color in vertex_colors.items()
                    if color == old_color
                )
                for old_color in set(vertex_colors.values())
            }
            old_coordinate_partition = {
                frozenset(
                    coordinate
                    for coordinate, color in coordinate_colors.items()
                    if color == old_color
                )
                for old_color in set(coordinate_colors.values())
            }

            vertex_signatures = {
                vertex: (
                    vertex_colors[vertex],
                    tuple(
                        sorted(
                            coordinate_colors[coordinate]
                            for coordinate in range(dimension)
                            if vertex & (1 << coordinate)
                        )
                    ),
                )
                for vertex in translated_members
            }
            vertex_palette = {
                signature: index
                for index, signature in enumerate(
                    sorted(set(vertex_signatures.values()))
                )
            }
            vertex_colors = {
                vertex: vertex_palette[signature]
                for vertex, signature in vertex_signatures.items()
            }

            coordinate_signatures = {
                coordinate: (
                    coordinate_colors[coordinate],
                    tuple(
                        sorted(
                            vertex_colors[vertex]
                            for vertex in translated_members
                            if vertex & (1 << coordinate)
                        )
                    ),
                )
                for coordinate in range(dimension)
            }
            coordinate_palette = {
                signature: index
                for index, signature in enumerate(
                    sorted(set(coordinate_signatures.values()))
                )
            }
            coordinate_colors = {
                coordinate: coordinate_palette[signature]
                for coordinate, signature in coordinate_signatures.items()
            }
            new_vertex_partition = {
                frozenset(
                    vertex
                    for vertex, color in vertex_colors.items()
                    if color == new_color
                )
                for new_color in set(vertex_colors.values())
            }
            new_coordinate_partition = {
                frozenset(
                    coordinate
                    for coordinate, color in coordinate_colors.items()
                    if color == new_color
                )
                for new_color in set(coordinate_colors.values())
            }
            if (
                old_vertex_partition == new_vertex_partition
                and old_coordinate_partition == new_coordinate_partition
            ):
                break

        ordered_groups = tuple(
            tuple(
                coordinate
                for coordinate in range(dimension)
                if coordinate_colors[coordinate] == color
            )
            for color in sorted(set(coordinate_colors.values()))
        )
        for within_group_orders in itertools.product(
            *(
                tuple(itertools.permutations(group))
                for group in ordered_groups
            )
        ):
            permutation = tuple(
                coordinate
                for group_order in within_group_orders
                for coordinate in group_order
            )
            transformed = transform_family(
                family_from_vertices(translated_members),
                dimension,
                permutation,
                0,
            )
            if best is None or transformed < best:
                best = transformed
    assert best is not None
    return best


def cube_symmetry_profile(family: int, dimension: int) -> tuple[int, int]:
    """Return the canonical bit mask and cube-automorphism group order."""

    canonical: int | None = None
    automorphisms = 0
    for permutation in itertools.permutations(range(dimension)):
        for flip in range(1 << dimension):
            transformed = transform_family(
                family, dimension, permutation, flip
            )
            if canonical is None or transformed < canonical:
                canonical = transformed
            if transformed == family:
                automorphisms += 1
    assert canonical is not None
    return canonical, automorphisms


def structural_profile(family: int, dimension: int) -> dict[str, object]:
    members = vertices(family, dimension)
    member_set = set(members)
    degrees = Counter(
        sum(
            (vertex ^ (1 << coordinate)) in member_set
            for coordinate in range(dimension)
        )
        for vertex in members
    )
    cubes = maximal_cubes(family, dimension)
    cube_dimensions = Counter((len(cube).bit_length() - 1) for cube in cubes)
    canonical, automorphisms = cube_symmetry_profile(family, dimension)
    shattered_histogram = shattered_size_histogram(family, dimension)
    vc_dimension = max(
        index for index, count in enumerate(shattered_histogram) if count
    )
    corners = ccmw_corners(family, dimension)
    return {
        "VC_dimension": vc_dimension,
        "shattered_size_histogram": shattered_histogram,
        "degree_histogram": dict(sorted(degrees.items())),
        "maximal_cube_dimension_histogram": dict(
            sorted(cube_dimensions.items())
        ),
        "CCMW_corners": len(corners),
        "corner_peeling": has_corner_peeling(family, dimension),
        "cube_automorphisms": automorphisms,
        "canonical_cube_mask": canonical,
    }


def is_antipodally_punctured_cube(family: int, dimension: int) -> bool:
    """Use uniqueness of irredundant hypercube coordinates for the test."""

    if family.bit_count() != (1 << dimension) - 2:
        return False
    holes = [
        vertex
        for vertex in range(1 << dimension)
        if not (family >> vertex & 1)
    ]
    return len(holes) == 2 and (holes[0] ^ holes[1]).bit_count() == dimension


def pentagonal_maximum_class() -> int:
    """The 16-concept maximum VC-dimension-two class carried by a 5-cycle."""

    cycle_edges = ((0, 1), (1, 2), (2, 4), (4, 3), (3, 0))
    full = (1 << 5) - 1
    concepts = [0]
    concepts.extend(1 << coordinate for coordinate in range(5))
    for left, right in cycle_edges:
        edge = (1 << left) | (1 << right)
        concepts.append(edge)
        concepts.append(full ^ edge)
    return family_from_vertices(concepts)


def hall_22_obstruction() -> int:
    """One 22-vertex minimal obstruction extracted from Hall's class."""

    return family_from_vertices(
        int(token, 2)
        for token in """
000111 001100 001101 001111 010000 010100 010101 010111
011100 011101 011111 101000 101001 101010 101011 101100
101101 101111 110000 110100 111000 111100
""".split()
    )


def bitstrings(family: int, dimension: int) -> str:
    return "{" + ",".join(
        format(vertex, f"0{dimension}b")
        for vertex in vertices(family, dimension)
    ) + "}"


def non_r_children(
    family: int, dimension: int
) -> list[tuple[int, int, str]]:
    return [
        child
        for child in immediate_pc_minors_with_operations(family, dimension)
        if not is_r_per(child[0], child[1])
    ]


def greedy_minimal_obstruction(
    start_family: int,
    start_dimension: int,
    strategy: str,
    rng: random.Random,
) -> tuple[int, int, list[tuple[int, int, str]]]:
    """Descend through non-R immediate minors until every child is in R."""

    family, dimension = irredundant_projection(start_family, start_dimension)
    path: list[tuple[int, int, str]] = [(family, dimension, "start")]
    while True:
        children = non_r_children(family, dimension)
        if not children:
            return family, dimension, path
        if strategy == "small-order":
            family, dimension, operation = min(
                children, key=lambda child: (child[0].bit_count(), child[1], child[0])
            )
        elif strategy == "large-order":
            family, dimension, operation = max(
                children, key=lambda child: (child[0].bit_count(), -child[1], -child[0])
            )
        elif strategy == "few-peripheral":
            family, dimension, operation = min(
                children,
                key=lambda child: (
                    len(peripheral_decompositions(child[0], child[1])),
                    child[0].bit_count(),
                    child[0],
                ),
            )
        elif strategy == "random":
            family, dimension, operation = rng.choice(children)
        else:
            raise ValueError(f"unknown strategy: {strategy}")
        path.append((family, dimension, operation))


def verify_obstruction(family: int, dimension: int) -> dict[str, object]:
    family, dimension = irredundant_projection(family, dimension)
    immediate = immediate_pc_minors_with_operations(family, dimension)
    shattered, strongly_shattered, ample = shattering_profile(family, dimension)
    report: dict[str, object] = {
        "dimension": dimension,
        "order": family.bit_count(),
        "partial_cube": is_partial_cube_family(family, dimension),
        "periphery_free": not peripheral_decompositions(family, dimension),
        "R_per": is_r_per(family, dimension),
        "immediate_minors": len(immediate),
        "all_immediate_minors_R_per": all(
            is_r_per(minor, minor_dimension)
            for minor, minor_dimension, _ in immediate
        ),
        "shattered": shattered,
        "strongly_shattered": strongly_shattered,
        "ample": ample,
        "crossing_pairs": len(crossing_pairs(family, dimension)),
        "complete_crossing_graph": not noncrossing_pairs(
            family, dimension
        ),
        "is_Q_n_minus_minus": is_antipodally_punctured_cube(
            family, dimension
        ),
    }
    return report


def gated_prime_profile(family: int, dimension: int) -> dict[str, object]:
    """Exhaust coordinate faces and test proper gated two-covers.

    Every convex subgraph of a partial cube is obtained by fixing a set of
    hypercube coordinates, so the loop has exactly ``3**dimension`` cases.
    A gated-amalgam decomposition must in particular cover the vertex set.
    """

    family, dimension = irredundant_projection(family, dimension)
    gated: set[int] = set()
    face_count = 0
    for fixed in range(1 << dimension):
        value = fixed
        while True:
            face_count += 1
            face = family_from_vertices(
                vertex
                for vertex in vertices(family, dimension)
                if vertex & fixed == value
            )
            if (
                face
                and face != family
                and is_gated_subfamily(face, family, dimension)
            ):
                gated.add(face)
            if value == 0:
                break
            value = (value - 1) & fixed

    covers = 0
    gated_list = sorted(gated)
    for left_index, left in enumerate(gated_list):
        for right in gated_list[left_index + 1 :]:
            if left & right and left | right == family:
                covers += 1
    return {
        "coordinate_faces": face_count,
        "proper_gated_faces": len(gated),
        "proper_gated_order_histogram": dict(
            sorted(Counter(face.bit_count() for face in gated).items())
        ),
        "proper_intersecting_two_covers": covers,
        "gated_prime": covers == 0,
    }


def components_after_deletion(
    base: int, dimension: int, deleted: int
) -> tuple[int, ...]:
    """Components of the one-inclusion graph induced by ``base - deleted``."""

    unseen = set(vertices(base & ~deleted, dimension))
    components: list[int] = []
    while unseen:
        source = min(unseen)
        unseen.remove(source)
        pending = [source]
        members: list[int] = []
        while pending:
            vertex = pending.pop()
            members.append(vertex)
            for coordinate in range(dimension):
                neighbor = vertex ^ (1 << coordinate)
                if neighbor in unseen:
                    unseen.remove(neighbor)
                    pending.append(neighbor)
        components.append(family_from_vertices(members))
    return tuple(components)


def lift_expansion(zero: int, one: int, base_dimension: int) -> int:
    return zero | family_from_vertices(
        vertex | (1 << base_dimension)
        for vertex in vertices(one, base_dimension)
    )


def search_dimension_five(order_limit: int) -> None:
    """Exhaust expansions of all dimension-four partial cubes up to an order.

    For a coordinate of a dimension-five partial cube, its contraction is a
    dimension-four partial cube ``B``.  Its two projected sections ``A_0`` and
    ``A_1`` are isometric subgraphs covering ``B``.  If ``I=A_0 intersect A_1``,
    no edge of ``B-I`` joins the two sides.  Hence its connected components can
    be assigned independently to the two sections.  Conversely, the exact
    partial-cube check below discards every assignment that is not an expansion.
    Enumerating canonical ``B`` and all such assignments is therefore complete
    up to cube symmetry.
    """

    if order_limit < 6 or order_limit > 32:
        raise ValueError("the dimension-five order limit must lie in [6,32]")
    bases = enumerate_partial_cubes(4)[4]
    intersections_tested = 0
    assignments_tested = 0
    partial_cube_expansions: set[int] = set()
    periphery_free: set[int] = set()
    minimal_obstructions: set[int] = set()
    additional_obstructions: set[int] = set()

    for base in bases:
        base_members = vertices(base, 4)
        maximum_intersection = min(
            len(base_members), order_limit - len(base_members)
        )
        for intersection_order in range(1, maximum_intersection + 1):
            for chosen in itertools.combinations(
                base_members, intersection_order
            ):
                intersection = family_from_vertices(chosen)
                components = components_after_deletion(base, 4, intersection)
                intersections_tested += 1
                # A periphery-free new coordinate needs vertices exclusive to
                # both sections, hence at least two components of B-I.
                if len(components) < 2:
                    continue
                # Swapping the new coordinate exchanges the sections.  Put the
                # first component on side zero to remove exactly that symmetry.
                for assignment in range(1 << (len(components) - 1)):
                    zero = intersection | components[0]
                    one = intersection
                    for index, component in enumerate(components[1:]):
                        if assignment >> index & 1:
                            zero |= component
                        else:
                            one |= component
                    if one == intersection:
                        continue
                    assignments_tested += 1
                    if not is_partial_cube_family(zero, 4):
                        continue
                    if not is_partial_cube_family(one, 4):
                        continue
                    expansion = lift_expansion(zero, one, 4)
                    if expansion.bit_count() > order_limit:
                        continue
                    if not is_partial_cube_family(expansion, 5):
                        continue
                    partial_cube_expansions.add(expansion)
                    if peripheral_decompositions(expansion, 5):
                        continue
                    periphery_free.add(expansion)
                    if not all(
                        is_r_per(minor, minor_dimension)
                        for minor, minor_dimension, _
                        in immediate_pc_minors_with_operations(expansion, 5)
                    ):
                        continue
                    minimal_obstructions.add(expansion)
                    if not is_antipodally_punctured_cube(expansion, 5):
                        additional_obstructions.add(expansion)

    print("\nExact dimension-five expansion census")
    print(
        f"  order limit={order_limit}; dimension-four bases={len(bases)}; "
        f"intersections={intersections_tested}; "
        f"side assignments={assignments_tested}"
    )
    print(
        f"  labelled expansions retained: partial_cubes={len(partial_cube_expansions)}, "
        f"periphery_free={len(periphery_free)}, "
        f"minimal_obstructions={len(minimal_obstructions)}, "
        f"additional_minimal_obstructions={len(additional_obstructions)}"
    )
    if additional_obstructions:
        best = min(
            additional_obstructions,
            key=lambda family: (family.bit_count(), family),
        )
        print(
            f"  smallest additional witness: order={best.bit_count()}, "
            f"vertices={bitstrings(best, 5)}"
        )
        pentagon = pentagonal_maximum_class()
        print(
            "  pentagonal maximum-class model: "
            f"exact_match={best == pentagon}; "
            + ", ".join(
                f"{key}={value}"
                for key, value in verify_obstruction(pentagon, 5).items()
            )
        )
        print(
            "  pentagonal structural profile: "
            + ", ".join(
                f"{key}={value}"
                for key, value in structural_profile(pentagon, 5).items()
            )
        )

    ample_periphery_free = {
        family
        for family in periphery_free
        if shattering_profile(family, 5)[2]
    }
    incomplete = {
        family
        for family in ample_periphery_free
        if noncrossing_pairs(family, 5)
    }
    pair_without_semicube_witnesses: list[
        tuple[int, tuple[int, int]]
    ] = []
    for family in incomplete:
        for pair in noncrossing_pairs(family, 5):
            if not periphery_free_semicube_witnesses(
                family, 5, *pair
            ):
                pair_without_semicube_witnesses.append((family, pair))
    print(
        "  pair-shattering diagnostic: "
        f"ample_periphery_free={len(ample_periphery_free)}, "
        f"incomplete_crossing={len(incomplete)}, "
        "noncrossing_pairs_without_periphery_free_semicube="
        f"{len(pair_without_semicube_witnesses)}"
    )
    if pair_without_semicube_witnesses:
        family, pair = min(
            pair_without_semicube_witnesses,
            key=lambda item: (item[0].bit_count(), item[1], item[0]),
        )
        print(
            "  smallest failure of the semicube-witness proxy: "
            f"order={family.bit_count()}, pair={pair}, "
            f"vertices={bitstrings(family, 5)}"
        )


def iter_expansion_candidates(
    bases: Iterable[int], base_dimension: int, order_limit: int
) -> Iterable[int]:
    """Yield all partial-cube expansions up to ``order_limit``.

    The input needs only one representative of each base isomorphism class:
    every cube isomorphism of a base transports its two expansion sections.
    Different coordinates or covers may yield the same labelled mask; callers
    choose whether to retain those duplicates.
    """

    for base in bases:
        base_members = vertices(base, base_dimension)
        maximum_intersection = min(
            len(base_members), order_limit - len(base_members)
        )
        for intersection_order in range(1, maximum_intersection + 1):
            for chosen in itertools.combinations(
                base_members, intersection_order
            ):
                intersection = family_from_vertices(chosen)
                components = components_after_deletion(
                    base, base_dimension, intersection
                )
                if not components:
                    # Both sections equal the base.
                    expansion = lift_expansion(
                        intersection, intersection, base_dimension
                    )
                    if expansion.bit_count() <= order_limit:
                        yield expansion
                    continue
                # Fix the first component on side zero to quotient the flip of
                # the new coordinate.  Other components may all join it: this
                # includes peripheral expansions whose side one is I.
                for assignment in range(1 << (len(components) - 1)):
                    zero = intersection | components[0]
                    one = intersection
                    for index, component in enumerate(components[1:]):
                        if assignment >> index & 1:
                            zero |= component
                        else:
                            one |= component
                    if not is_partial_cube_family(zero, base_dimension):
                        continue
                    if not is_partial_cube_family(one, base_dimension):
                        continue
                    expansion = lift_expansion(zero, one, base_dimension)
                    if expansion.bit_count() > order_limit:
                        continue
                    # The two isometric sections cover ``base`` and their
                    # exclusive parts are unions of distinct components of
                    # ``base - intersection``.  The partial-cube expansion
                    # theorem therefore makes the lift isometric.  The exact
                    # dimension-six and dimension-seven routines separately
                    # verified this implication by direct distances on every
                    # retained expansion.
                    yield expansion


def expansion_candidates(
    bases: Iterable[int], base_dimension: int, order_limit: int
) -> set[int]:
    """Return the distinct labelled outputs of ``iter_expansion_candidates``."""

    return set(
        iter_expansion_candidates(bases, base_dimension, order_limit)
    )


@lru_cache(maxsize=None)
def coordinate_low_vertex_masks(dimension: int) -> tuple[int, ...]:
    """Return cube-vertex bitmasks selecting value zero in each coordinate."""

    masks: list[int] = []
    vertex_count = 1 << dimension
    for coordinate in range(dimension):
        step = 1 << coordinate
        mask = 0
        for block_start in range(0, vertex_count, 2 * step):
            mask |= ((1 << step) - 1) << block_start
        masks.append(mask)
    return tuple(masks)


def expansion_has_peripheral_coordinate(
    zero: int,
    one: int,
    base_dimension: int,
) -> bool:
    """Test periphery directly on the two projected expansion sections.

    For an old coordinate ``c``, the two lifted ``c``-sections are nested
    exactly when the corresponding sections of both ``zero`` and ``one`` are
    nested in the same direction.  In the cube-vertex bitset, flipping
    coordinate ``c`` shifts a low-section bit by ``2**c`` positions.  Thus
    the two inclusions can be tested directly, without constructing projected
    sections.  The new coordinate is peripheral precisely when the two input
    sections themselves are nested.
    """

    if nested(zero, one) or nested(one, zero):
        return True
    universe = (1 << (1 << base_dimension)) - 1
    for coordinate, low_mask in enumerate(
        coordinate_low_vertex_masks(base_dimension)
    ):
        shift = 1 << coordinate
        high_mask = universe ^ low_mask
        zero_low = zero & low_mask
        zero_high = zero & high_mask
        one_low = one & low_mask
        one_high = one & high_mask
        if (
            nested(zero_low << shift, zero_high)
            and nested(one_low << shift, one_high)
        ) or (
            nested(zero_high >> shift, zero_low)
            and nested(one_high >> shift, one_low)
        ):
            return True
    return False


def expansion_periphery_predicate_regression(
    max_base_dimension: int = 4,
    order_limit: int = 10,
) -> None:
    """Exhaustively compare direct and lifted periphery predicates."""

    if max_base_dimension < 1 or max_base_dimension > 4:
        raise ValueError("the regression supports base dimensions 1 through 4")
    census = enumerate_partial_cubes(max_base_dimension)
    tested = 0
    failures: list[tuple[int, int]] = []
    for base_dimension in range(1, max_base_dimension + 1):
        for expansion in iter_expansion_candidates(
            census[base_dimension],
            base_dimension,
            order_limit,
        ):
            zero = section(
                expansion,
                base_dimension + 1,
                base_dimension,
                0,
            )
            one = section(
                expansion,
                base_dimension + 1,
                base_dimension,
                1,
            )
            direct = expansion_has_peripheral_coordinate(
                zero, one, base_dimension
            )
            lifted = bool(
                peripheral_decompositions(
                    expansion, base_dimension + 1
                )
            )
            tested += 1
            if direct != lifted:
                failures.append((expansion, base_dimension + 1))
    print(
        "\nExpansion-periphery predicate regression\n"
        f"  max_base_dimension={max_base_dimension}; "
        f"order_limit={order_limit}; tested={tested}; "
        f"failures={len(failures)}",
        flush=True,
    )
    if failures:
        family, dimension = failures[0]
        raise AssertionError(
            "direct periphery predicate mismatch: "
            f"dimension={dimension}, vertices={bitstrings(family, dimension)}"
        )


FAMILY_CENSUS_SCHEMA = "partial-cube-family-census-v1"
CANONICALIZER_NAME = (
    f"tree-first-colored-incidence-pynauty-{pynauty.__version__}-v2"
    if pynauty is not None
    else "refined-incidence-tree-python-v2"
)


def family_census_digest(families: Iterable[int]) -> str:
    """Digest a canonical sorted family list independently of JSON layout."""

    payload = "".join(
        f"{family:x}\n" for family in sorted(set(families))
    ).encode("ascii")
    return hashlib.sha256(payload).hexdigest()


def write_family_census_artifact(
    path: str,
    *,
    kind: str,
    dimension: int,
    order_limit: int,
    families: Iterable[int],
    metadata: dict[str, object] | None = None,
) -> None:
    """Write a deterministic, integrity-checked exact-census artifact."""

    ordered = sorted(set(families))
    artifact: dict[str, object] = {
        "schema": FAMILY_CENSUS_SCHEMA,
        "kind": kind,
        "dimension": dimension,
        "order_limit": order_limit,
        "canonicalizer": CANONICALIZER_NAME,
        "family_count": len(ordered),
        "family_digest": family_census_digest(ordered),
        "families_hex": [f"{family:x}" for family in ordered],
    }
    if metadata:
        artifact.update(metadata)
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp")
    temporary.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(output)
    print(
        f"  wrote artifact={output}; families={len(ordered)}; "
        f"digest={artifact['family_digest']}",
        flush=True,
    )


def read_family_census_artifact(path: str) -> tuple[dict[str, object], set[int]]:
    """Read an artifact and verify its schema, count, and family digest."""

    artifact = json.loads(Path(path).read_text(encoding="utf-8"))
    if artifact.get("schema") != FAMILY_CENSUS_SCHEMA:
        raise ValueError(f"{path}: unsupported family-census schema")
    encoded = artifact.get("families_hex")
    if not isinstance(encoded, list) or not all(
        isinstance(item, str) for item in encoded
    ):
        raise ValueError(f"{path}: malformed family list")
    families = {int(item, 16) for item in encoded}
    if len(families) != artifact.get("family_count"):
        raise ValueError(f"{path}: duplicate family or incorrect count")
    if family_census_digest(families) != artifact.get("family_digest"):
        raise ValueError(f"{path}: family digest mismatch")
    dimension = artifact.get("dimension")
    order_limit = artifact.get("order_limit")
    if not isinstance(dimension, int) or not isinstance(order_limit, int):
        raise ValueError(f"{path}: malformed census parameters")
    if any(
        family.bit_count() > order_limit
        or irredundant_projection(family, dimension)[1] != dimension
        for family in families
    ):
        raise ValueError(f"{path}: family outside the declared census bounds")
    return artifact, families


def ample_type_census_chunk(
    source_path: str,
    target_dimension: int,
    order_limit: int,
    base_start: int,
    base_stop: int | None,
    output_path: str,
) -> None:
    """Generate one deterministic chunk of the next ample-type layer.

    Base indices refer to the sorted source families of order strictly below
    ``order_limit``.  Separate chunks can therefore be generated in distinct
    short-lived processes and later merged with an exact interval-cover check.
    """

    source, source_families = read_family_census_artifact(source_path)
    if source["dimension"] != target_dimension - 1:
        raise ValueError("the source artifact has the wrong base dimension")
    if target_dimension < 5 or target_dimension > 14:
        raise ValueError("the target dimension must lie in [5,14]")
    if order_limit < target_dimension + 1 or order_limit > 16:
        raise ValueError("the order bound is incompatible with the dimension")
    bases = sorted(
        family
        for family in source_families
        if family.bit_count() < order_limit
    )
    total_bases = len(bases)
    effective_stop = (
        total_bases if base_stop is None else min(base_stop, total_bases)
    )
    if base_start < 0 or base_start >= effective_stop:
        raise ValueError("the requested base range is empty or reversed")

    started = time.perf_counter()
    ample_labelled: set[int] = set()
    base_dimension = target_dimension - 1
    layer_bit_count = 1 << base_dimension
    zero_layer_mask = (1 << layer_bit_count) - 1
    for local_index, base in enumerate(
        bases[base_start:effective_stop], start=1
    ):
        for family in iter_expansion_candidates(
            (base,), base_dimension, order_limit
        ):
            zero = family & zero_layer_mask
            one = family >> layer_bit_count
            if ample_expansion_of_ample_base(
                zero, one, base_dimension
            ):
                ample_labelled.add(family)
        if local_index % 100 == 0:
            is_partial_cube_family.cache_clear()
            shattered_coordinate_sets_mask.cache_clear()
            section.cache_clear()
            irredundant_projection.cache_clear()
            peripheral_decompositions.cache_clear()
            canonical_cube_family.cache_clear()
            gc.collect()
        if local_index % 10 == 0:
            print(
                f"  generation={base_start + local_index}/"
                f"{total_bases}; labelled={len(ample_labelled)}; "
                "peak_rss_mb="
                f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024:.1f}",
                flush=True,
            )

    ample_types: set[int] = set()
    for family_index, family in enumerate(
        sorted(ample_labelled), start=1
    ):
        ample_types.add(
            refined_canonical_cube_family(family, target_dimension)
        )
        if family_index % 500 == 0:
            gc.collect()
        if family_index % 2000 == 0:
            print(
                f"  canonicalization={family_index}/"
                f"{len(ample_labelled)}; types={len(ample_types)}",
                flush=True,
            )

    histogram = Counter(family.bit_count() for family in ample_types)
    print(
        "\nExact ample-type census chunk\n"
        f"  target_dimension={target_dimension}; order_limit={order_limit}; "
        f"source_types={len(source_families)}; eligible_bases={total_bases}; "
        f"base_range=[{base_start},{effective_stop})\n"
        f"  labelled={len(ample_labelled)}; types={len(ample_types)}; "
        f"order_histogram={dict(sorted(histogram.items()))}; "
        "peak_rss_mb="
        f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024:.1f}; "
        f"runtime_seconds={time.perf_counter() - started:.3f}",
        flush=True,
    )
    write_family_census_artifact(
        output_path,
        kind="ample-type-chunk",
        dimension=target_dimension,
        order_limit=order_limit,
        families=ample_types,
        metadata={
            "source_family_digest": source["family_digest"],
            "source_family_count": len(source_families),
            "eligible_base_count": total_bases,
            "base_start": base_start,
            "base_stop": effective_stop,
        },
    )


def merge_ample_type_census_chunks(
    chunk_paths: Iterable[str],
    output_path: str,
) -> None:
    """Merge chunks only when their base ranges form one exact full cover."""

    chunks: list[tuple[dict[str, object], set[int], str]] = []
    for path in chunk_paths:
        artifact, families = read_family_census_artifact(path)
        if artifact.get("kind") != "ample-type-chunk":
            raise ValueError(f"{path}: not an ample-type chunk")
        chunks.append((artifact, families, path))
    if not chunks:
        raise ValueError("at least one chunk is required")

    reference = chunks[0][0]
    consistency_fields = (
        "dimension",
        "order_limit",
        "canonicalizer",
        "source_family_digest",
        "source_family_count",
        "eligible_base_count",
    )
    for artifact, _, path in chunks[1:]:
        if any(
            artifact.get(field) != reference.get(field)
            for field in consistency_fields
        ):
            raise ValueError(f"{path}: inconsistent chunk metadata")

    ordered_chunks = sorted(
        chunks, key=lambda item: int(item[0]["base_start"])
    )
    expected_start = 0
    merged: set[int] = set()
    ranges: list[list[int]] = []
    for artifact, families, path in ordered_chunks:
        start = int(artifact["base_start"])
        stop = int(artifact["base_stop"])
        if start != expected_start or stop <= start:
            raise ValueError(
                f"{path}: chunk ranges do not form an exact interval cover"
            )
        expected_start = stop
        ranges.append([start, stop])
        merged.update(families)
    if expected_start != int(reference["eligible_base_count"]):
        raise ValueError("the chunks do not cover every eligible base")

    histogram = Counter(family.bit_count() for family in merged)
    print(
        "\nMerged exact ample-type census\n"
        f"  dimension={reference['dimension']}; "
        f"order_limit={reference['order_limit']}; "
        f"chunks={len(ordered_chunks)}; types={len(merged)}; "
        f"order_histogram={dict(sorted(histogram.items()))}",
        flush=True,
    )
    write_family_census_artifact(
        output_path,
        kind="ample-type-census",
        dimension=int(reference["dimension"]),
        order_limit=int(reference["order_limit"]),
        families=merged,
        metadata={
            "complete": True,
            "source_family_digest": reference["source_family_digest"],
            "source_family_count": reference["source_family_count"],
            "eligible_base_count": reference["eligible_base_count"],
            "chunk_ranges": ranges,
        },
    )


def periphery_free_search_chunk(
    source_path: str,
    target_dimension: int,
    order_limit: int,
    base_start: int,
    base_stop: int | None,
    output_path: str,
) -> None:
    """Search one exact base interval for ample periphery-free expansions."""

    source, source_families = read_family_census_artifact(source_path)
    if source.get("kind") != "ample-type-census" or not source.get(
        "complete"
    ):
        raise ValueError("the source must be a complete ample-type census")
    if source["dimension"] != target_dimension - 1:
        raise ValueError("the source artifact has the wrong base dimension")
    if int(source["order_limit"]) < order_limit - 1:
        raise ValueError("the source census does not reach the required order")
    if target_dimension < 6 or target_dimension > 14:
        raise ValueError("the target dimension must lie in [6,14]")
    if order_limit < target_dimension + 1 or order_limit > 18:
        raise ValueError("the order bound is incompatible with the dimension")

    bases = sorted(
        family
        for family in source_families
        if family.bit_count() < order_limit
    )
    total_bases = len(bases)
    effective_stop = (
        total_bases if base_stop is None else min(base_stop, total_bases)
    )
    if base_start < 0 or base_start >= effective_stop:
        raise ValueError("the requested base range is empty or reversed")

    started = time.perf_counter()
    intersections = 0
    assignments = 0
    periphery_free_assignments = 0
    section_isometric = 0
    periphery_free_partial_cubes = 0
    ample_periphery_free = 0
    periphery_free: set[int] = set()
    minimal: set[int] = set()
    base_dimension = target_dimension - 1

    for local_index, base in enumerate(
        bases[base_start:effective_stop], start=1
    ):
        members = vertices(base, base_dimension)
        maximum_intersection = min(
            len(members), order_limit - len(members)
        )
        for intersection_order in range(1, maximum_intersection + 1):
            for chosen in itertools.combinations(
                members, intersection_order
            ):
                intersection = family_from_vertices(chosen)
                components = components_after_deletion(
                    base, base_dimension, intersection
                )
                intersections += 1
                if len(components) < 2:
                    continue
                for assignment in range(1 << (len(components) - 1)):
                    zero = intersection | components[0]
                    one = intersection
                    for index, component in enumerate(components[1:]):
                        if assignment >> index & 1:
                            zero |= component
                        else:
                            one |= component
                    if one == intersection:
                        continue
                    assignments += 1
                    # Periphery is a section-nesting condition and can be
                    # rejected before the more expensive isometry tests.  A
                    # rejected raw cover cannot be a target irrespective of
                    # whether its sides define a valid partial-cube
                    # expansion.
                    if expansion_has_peripheral_coordinate(
                        zero, one, base_dimension
                    ):
                        continue
                    periphery_free_assignments += 1
                    if not is_partial_cube_family(zero, base_dimension):
                        continue
                    if not is_partial_cube_family(one, base_dimension):
                        continue
                    section_isometric += 1
                    # The two isometric sections cover the base and their
                    # exclusive parts are unions of distinct components of
                    # the deleted intersection.  The partial-cube expansion
                    # theorem therefore makes the lift a partial cube; the
                    # exhaustive dimension-six and dimension-seven routines
                    # retain an independent direct-distance cross-check.
                    periphery_free_partial_cubes += 1
                    expansion = lift_expansion(
                        zero, one, base_dimension
                    )
                    if not is_ample_family(
                        expansion, target_dimension
                    ):
                        continue
                    ample_periphery_free += 1
                    canonical = refined_canonical_cube_family(
                        expansion, target_dimension
                    )
                    periphery_free.add(canonical)
                    if all(
                        is_r_per(minor, minor_dimension)
                        for minor, minor_dimension, _
                        in immediate_pc_minors_with_operations(
                            expansion, target_dimension
                        )
                    ):
                        minimal.add(canonical)

        if local_index % 100 == 0:
            is_partial_cube_family.cache_clear()
            section.cache_clear()
            irredundant_projection.cache_clear()
            peripheral_decompositions.cache_clear()
            is_r_per.cache_clear()
            canonical_cube_family.cache_clear()
            gc.collect()
            print(
                f"  search={base_start + local_index}/"
                f"{total_bases}; assignments={assignments}; "
                "periphery_free_assignments="
                f"{periphery_free_assignments}; "
                "periphery_free_partial_cubes="
                f"{periphery_free_partial_cubes}; "
                f"ample_periphery_free={ample_periphery_free}; "
                f"periphery_free_types={len(periphery_free)}; "
                "peak_rss_mb="
                f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024:.1f}",
                flush=True,
            )

    counts = {
        "intersections": intersections,
        "assignments": assignments,
        "periphery_free_assignments": periphery_free_assignments,
        "periphery_free_isometric_section_pairs": section_isometric,
        "periphery_free_partial_cube_expansions":
            periphery_free_partial_cubes,
        "ample_periphery_free_expansions": ample_periphery_free,
    }
    print(
        "\nExact periphery-free search chunk\n"
        f"  target_dimension={target_dimension}; order_limit={order_limit}; "
        f"source_types={len(source_families)}; eligible_bases={total_bases}; "
        f"base_range=[{base_start},{effective_stop})\n"
        f"  counts={counts}; "
        f"periphery_free_types={len(periphery_free)}; "
        f"minimal_obstruction_types={len(minimal)}; "
        "peak_rss_mb="
        f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024:.1f}; "
        f"runtime_seconds={time.perf_counter() - started:.3f}",
        flush=True,
    )
    for family in sorted(
        periphery_free, key=lambda item: (item.bit_count(), item)
    ):
        print(
            "  periphery-free type: "
            f"order={family.bit_count()}, "
            f"minimal={family in minimal}, "
            f"vertices={bitstrings(family, target_dimension)}",
            flush=True,
        )
    write_family_census_artifact(
        output_path,
        kind="periphery-free-search-chunk",
        dimension=target_dimension,
        order_limit=order_limit,
        families=periphery_free,
        metadata={
            "source_family_digest": source["family_digest"],
            "source_family_count": len(source_families),
            "eligible_base_count": total_bases,
            "base_start": base_start,
            "base_stop": effective_stop,
            "counts": counts,
            "minimal_families_hex": [
                f"{family:x}" for family in sorted(minimal)
            ],
        },
    )


def merge_periphery_free_search_chunks(
    chunk_paths: Iterable[str],
    output_path: str,
) -> None:
    """Merge a full interval cover of exact periphery-free search chunks."""

    chunks: list[tuple[dict[str, object], set[int], str]] = []
    for path in chunk_paths:
        artifact, families = read_family_census_artifact(path)
        if artifact.get("kind") != "periphery-free-search-chunk":
            raise ValueError(f"{path}: not a periphery-free search chunk")
        chunks.append((artifact, families, path))
    if not chunks:
        raise ValueError("at least one search chunk is required")

    reference = chunks[0][0]
    consistency_fields = (
        "dimension",
        "order_limit",
        "canonicalizer",
        "source_family_digest",
        "source_family_count",
        "eligible_base_count",
    )
    for artifact, _, path in chunks[1:]:
        if any(
            artifact.get(field) != reference.get(field)
            for field in consistency_fields
        ):
            raise ValueError(f"{path}: inconsistent search metadata")

    ordered_chunks = sorted(
        chunks, key=lambda item: int(item[0]["base_start"])
    )
    expected_start = 0
    merged: set[int] = set()
    minimal: set[int] = set()
    totals: Counter[str] = Counter()
    ranges: list[list[int]] = []
    for artifact, families, path in ordered_chunks:
        start = int(artifact["base_start"])
        stop = int(artifact["base_stop"])
        if start != expected_start or stop <= start:
            raise ValueError(
                f"{path}: search ranges do not form an exact interval cover"
            )
        counts = artifact.get("counts")
        if not isinstance(counts, dict) or not all(
            isinstance(value, int) for value in counts.values()
        ):
            raise ValueError(f"{path}: malformed search counts")
        expected_start = stop
        ranges.append([start, stop])
        totals.update(counts)
        merged.update(families)
        encoded_minimal = artifact.get("minimal_families_hex")
        if not isinstance(encoded_minimal, list):
            raise ValueError(f"{path}: malformed minimal-family list")
        minimal.update(int(item, 16) for item in encoded_minimal)
    if expected_start != int(reference["eligible_base_count"]):
        raise ValueError("the chunks do not cover every eligible base")
    if totals["periphery_free_isometric_section_pairs"] != totals[
        "periphery_free_partial_cube_expansions"
    ]:
        raise AssertionError("the expansion-theorem accounting failed")
    if not minimal.issubset(merged):
        raise AssertionError("a minimal type is absent from the merged output")

    print(
        "\nMerged exact periphery-free search\n"
        f"  dimension={reference['dimension']}; "
        f"order_limit={reference['order_limit']}; "
        f"chunks={len(ordered_chunks)}; counts={dict(totals)}; "
        f"periphery_free_types={len(merged)}; "
        f"minimal_obstruction_types={len(minimal)}",
        flush=True,
    )
    write_family_census_artifact(
        output_path,
        kind="periphery-free-search",
        dimension=int(reference["dimension"]),
        order_limit=int(reference["order_limit"]),
        families=merged,
        metadata={
            "complete": True,
            "source_family_digest": reference["source_family_digest"],
            "source_family_count": reference["source_family_count"],
            "eligible_base_count": reference["eligible_base_count"],
            "chunk_ranges": ranges,
            "counts": dict(totals),
            "minimal_families_hex": [
                f"{family:x}" for family in sorted(minimal)
            ],
        },
    )


def search_all_dimensions_below_order(order_limit: int) -> None:
    """Exhaust ample partial cubes in dimensions at least five up to an order.

    Knauer--Marc's forbidden-pc-minor theorem for ampleness implies that every
    minimal obstruction outside the known ``Q_n^{--}`` family is ample.  Since
    contractions preserve ampleness, it is sound to retain only ample bases at
    every expansion level.
    """

    if order_limit < 6 or order_limit > 18:
        raise ValueError("the global order limit must lie in [6,18]")
    started = time.perf_counter()
    bases = [
        family
        for family in enumerate_partial_cubes(4)[4]
        if family.bit_count() < order_limit
        and shattering_profile(family, 4)[2]
    ]
    print("\nBounded all-dimensional ample expansion census")
    print(
        f"  target order <= {order_limit}; dimension 4 ample bases={len(bases)}",
        flush=True,
    )
    additional: list[tuple[int, int]] = []
    # An irredundant d-dimensional partial cube has at least d+1 vertices.
    for dimension in range(5, order_limit):
        generated = expansion_candidates(bases, dimension - 1, order_limit)
        ample_raw = {
            family
            for family in generated
            if shattering_profile(family, dimension)[2]
        }
        ample = {
            canonical_cube_family(family, dimension)
            for family in ample_raw
        }
        candidates = {
            family
            for family in ample
            if not peripheral_decompositions(family, dimension)
            and all(
                is_r_per(minor, minor_dimension)
                for minor, minor_dimension, _
                in immediate_pc_minors_with_operations(family, dimension)
            )
            and not is_antipodally_punctured_cube(family, dimension)
        }
        additional.extend((family, dimension) for family in candidates)
        print(
            f"  dimension={dimension}: generated={len(generated)}, "
            f"ample_raw={len(ample_raw)}, ample_types={len(ample)}, "
            f"additional_minimal={len(candidates)}",
            flush=True,
        )
        # A base for a further nontrivial expansion must have at most
        # order_limit-1 vertices.
        bases = [
            family for family in ample if family.bit_count() < order_limit
        ]
        if not bases:
            break
    if additional:
        family, dimension = min(
            additional,
            key=lambda item: (item[0].bit_count(), item[1], item[0]),
        )
        print(
            f"  smallest found: dimension={dimension}, "
            f"order={family.bit_count()}, vertices={bitstrings(family, dimension)}"
        )
    else:
        print("  no additional minimal obstruction in the exhausted range")
    print(f"  runtime_seconds={time.perf_counter() - started:.3f}")


def random_ample_deletion_walks(
    dimension: int,
    trials: int,
    seed: int,
    report_below: int,
) -> None:
    """Stress-test sparse ample classes reached by corner deletions.

    The walk begins at the full cube and repeatedly removes a uniformly
    chosen vertex whose deletion remains ample and uses every coordinate.
    This samples only ample classes that admit the displayed reverse
    completion to the cube; it is therefore an adversarial diagnostic, not
    an exhaustive census of all ample classes.
    """

    if dimension < 3 or dimension > 8:
        raise ValueError("random ample walks support dimensions 3 through 8")
    if trials < 1:
        raise ValueError("the number of random walks must be positive")
    if report_below < dimension + 1:
        raise ValueError("the reporting threshold is below irredundant order")

    rng = random.Random(seed)
    full = (1 << (1 << dimension)) - 1
    seen_sparse: set[int] = set()
    periphery_free: set[int] = set()
    incomplete: set[int] = set()
    proxy_failures: list[tuple[int, tuple[int, int]]] = []
    below_pentagon: set[int] = set()
    terminal_sizes: Counter[int] = Counter()

    for _ in range(trials):
        family = full
        while family.bit_count() > dimension + 1:
            removable: list[int] = []
            members = list(vertices(family, dimension))
            rng.shuffle(members)
            for vertex in members:
                candidate = family & ~(1 << vertex)
                projected, projected_dimension = irredundant_projection(
                    candidate, dimension
                )
                if projected_dimension != dimension:
                    continue
                if shattering_profile(projected, dimension)[2]:
                    removable.append(projected)
            if not removable:
                break
            family = rng.choice(removable)
            if family.bit_count() > report_below:
                continue
            seen_sparse.add(family)
            if peripheral_decompositions(family, dimension):
                continue
            periphery_free.add(family)
            if family.bit_count() < 16:
                below_pentagon.add(family)
            pairs = noncrossing_pairs(family, dimension)
            if not pairs:
                continue
            incomplete.add(family)
            for pair in pairs:
                if not periphery_free_semicube_witnesses(
                    family, dimension, *pair
                ):
                    proxy_failures.append((family, pair))
        terminal_sizes[family.bit_count()] += 1

    print("\nRandom ample deletion-walk diagnostic")
    print(
        f"  dimension={dimension}; trials={trials}; seed={seed}; "
        f"report_below={report_below}"
    )
    print(
        f"  sparse_states={len(seen_sparse)}; "
        f"periphery_free={len(periphery_free)}; "
        f"incomplete_crossing={len(incomplete)}; "
        f"order_below_16={len(below_pentagon)}; "
        "noncrossing_pairs_without_periphery_free_semicube="
        f"{len(proxy_failures)}"
    )
    print(f"  terminal_order_histogram={dict(sorted(terminal_sizes.items()))}")
    if below_pentagon:
        family = min(
            below_pentagon, key=lambda item: (item.bit_count(), item)
        )
        print(
            "  smallest periphery-free state below order 16: "
            f"order={family.bit_count()}, "
            f"vertices={bitstrings(family, dimension)}"
        )
    if proxy_failures:
        family, pair = min(
            proxy_failures,
            key=lambda item: (item[0].bit_count(), item[1], item[0]),
        )
        print(
            "  smallest failure of the semicube-witness proxy: "
            f"order={family.bit_count()}, pair={pair}, "
            f"vertices={bitstrings(family, dimension)}"
        )


def random_sparse_dimension_six_expansions(
    trials: int,
    seed: int,
    order_limit: int = 15,
) -> None:
    """Sample sparse dimension-six expansions from exact ample bases.

    All dimension-five bases in the sampling pool are generated exactly as
    expansions of the dimension-four census.  The subsequent choice of the
    overlap and the component assignment is random, so this is not a
    completeness certificate for dimension six.
    """

    if trials < 1:
        raise ValueError("the number of sparse expansion trials must be positive")
    if order_limit < 7 or order_limit > 18:
        raise ValueError("the sparse dimension-six order limit must lie in [7,18]")

    started = time.perf_counter()
    rng = random.Random(seed)
    dimension_four_bases = [
        family
        for family in enumerate_partial_cubes(4)[4]
        if family.bit_count() < order_limit
        and shattering_profile(family, 4)[2]
    ]
    raw_dimension_five = expansion_candidates(
        dimension_four_bases, 4, order_limit
    )
    bases = [
        family
        for family in raw_dimension_five
        if family.bit_count() < order_limit
        and shattering_profile(family, 5)[2]
    ]

    sampled: set[int] = set()
    partial_cubes: set[int] = set()
    ample: set[int] = set()
    periphery_free: set[int] = set()
    minimal: set[int] = set()
    incomplete: set[int] = set()
    proxy_failures: list[tuple[int, tuple[int, int]]] = []

    for _ in range(trials):
        base = rng.choice(bases)
        members = vertices(base, 5)
        maximum_intersection = min(
            len(members), order_limit - len(members)
        )
        intersection_order = rng.randint(1, maximum_intersection)
        intersection = family_from_vertices(
            rng.sample(members, intersection_order)
        )
        components = components_after_deletion(base, 5, intersection)
        if len(components) < 2:
            continue
        zero = intersection | components[0]
        one = intersection
        for component in components[1:]:
            if rng.getrandbits(1):
                zero |= component
            else:
                one |= component
        if one == intersection:
            continue
        expansion = lift_expansion(zero, one, 5)
        sampled.add(expansion)
        if not is_partial_cube_family(zero, 5):
            continue
        if not is_partial_cube_family(one, 5):
            continue
        if not is_partial_cube_family(expansion, 6):
            continue
        partial_cubes.add(expansion)
        if not shattering_profile(expansion, 6)[2]:
            continue
        ample.add(expansion)
        if peripheral_decompositions(expansion, 6):
            continue
        periphery_free.add(expansion)
        if all(
            is_r_per(minor, minor_dimension)
            for minor, minor_dimension, _
            in immediate_pc_minors_with_operations(expansion, 6)
        ):
            minimal.add(expansion)
        pairs = noncrossing_pairs(expansion, 6)
        if not pairs:
            continue
        incomplete.add(expansion)
        for pair in pairs:
            if not periphery_free_semicube_witnesses(
                expansion, 6, *pair
            ):
                proxy_failures.append((expansion, pair))

    print("\nRandom sparse dimension-six expansion diagnostic")
    print(
        f"  trials={trials}; seed={seed}; order_limit={order_limit}; "
        f"d4_ample_bases={len(dimension_four_bases)}; "
        f"d5_ample_labelled_bases={len(bases)}"
    )
    print(
        f"  sampled={len(sampled)}; partial_cubes={len(partial_cubes)}; "
        f"ample={len(ample)}; periphery_free={len(periphery_free)}; "
        f"minimal_obstructions={len(minimal)}; "
        f"incomplete_crossing={len(incomplete)}; "
        "noncrossing_pairs_without_periphery_free_semicube="
        f"{len(proxy_failures)}"
    )
    print(f"  runtime_seconds={time.perf_counter() - started:.3f}")
    if periphery_free:
        family = min(
            periphery_free, key=lambda item: (item.bit_count(), item)
        )
        print(
            "  smallest periphery-free sample: "
            f"order={family.bit_count()}, "
            f"complete_crossing={not noncrossing_pairs(family, 6)}, "
            f"vertices={bitstrings(family, 6)}"
        )
    if minimal:
        family = min(minimal, key=lambda item: (item.bit_count(), item))
        print(
            "  smallest sampled minimal obstruction: "
            f"order={family.bit_count()}, "
            f"complete_crossing={not noncrossing_pairs(family, 6)}, "
            f"vertices={bitstrings(family, 6)}"
        )
    if proxy_failures:
        family, pair = min(
            proxy_failures,
            key=lambda item: (item[0].bit_count(), item[1], item[0]),
        )
        print(
            "  smallest failure of the semicube-witness proxy: "
            f"order={family.bit_count()}, pair={pair}, "
            f"vertices={bitstrings(family, 6)}"
        )


def search_dimension_six_below_order(
    order_limit: int = 15,
    base_start: int = 0,
    base_stop: int | None = None,
    canonicalize_ample: bool = False,
) -> set[int]:
    """Exhaust dimension-six ample expansions through an order bound.

    Every additional minimal obstruction is ample, and contractions preserve
    ampleness.  Hence a dimension-six obstruction of order at most the bound
    is an expansion of one of the ample dimension-five bases generated below.
    The overlap and every bipartition of the components outside it are then
    enumerated exactly, up to swapping the new coordinate.
    """

    if order_limit < 7 or order_limit > 18:
        raise ValueError("the dimension-six order limit must lie in [7,18]")
    if base_start < 0:
        raise ValueError("the first dimension-five base index must be nonnegative")

    started = time.perf_counter()
    dimension_four_bases = [
        family
        for family in enumerate_partial_cubes(4)[4]
        if family.bit_count() < order_limit
        and shattering_profile(family, 4)[2]
    ]
    raw_dimension_five = expansion_candidates(
        dimension_four_bases, 4, order_limit
    )
    print(
        "\nExact dimension-six ample expansion census",
        flush=True,
    )
    raw_dimension_five_count = len(raw_dimension_five)
    print(
        f"  generated dimension-five pool: raw={raw_dimension_five_count}",
        flush=True,
    )
    # Canonicalize as a bounded-memory stream.  The older implementation used
    # the unbounded-cache exhaustive canonicalizer in a set comprehension.
    # Its values are never reused in this one-pass preparation, so retaining
    # every input/result pair only increases the front-end memory footprint.
    # ``refined_canonical_cube_family`` is the exact, regression-tested
    # canonicalizer used by the resumable census pipeline and has no global
    # result cache.
    canonical_dimension_five: set[int] = set()
    eligible_dimension_five = 0
    for raw_index, family in enumerate(
        sorted(raw_dimension_five), start=1
    ):
        if (
            family.bit_count() < order_limit
            and shattering_profile(family, 5)[2]
        ):
            eligible_dimension_five += 1
            canonical_dimension_five.add(
                refined_canonical_cube_family(family, 5)
            )
        if raw_index % 500 == 0:
            canonical_cube_family.cache_clear()
            is_partial_cube_family.cache_clear()
            section.cache_clear()
            irredundant_projection.cache_clear()
            peripheral_decompositions.cache_clear()
            gc.collect()
            print(
                f"  dimension-five preparation={raw_index}/"
                f"{raw_dimension_five_count}; "
                f"eligible={eligible_dimension_five}; "
                f"types={len(canonical_dimension_five)}; "
                "peak_rss_mb="
                f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024:.1f}",
                flush=True,
            )
    dimension_five_bases = sorted(canonical_dimension_five)
    del canonical_dimension_five
    total_bases = len(dimension_five_bases)
    effective_stop = (
        total_bases
        if base_stop is None
        else min(base_stop, total_bases)
    )
    if base_start > effective_stop:
        raise ValueError("the dimension-six base range is empty or reversed")
    selected_bases = dimension_five_bases[base_start:effective_stop]
    complete_coverage = base_start == 0 and effective_stop == total_bases
    print(
        f"  canonical ample bases={total_bases}; "
        f"selected_range=[{base_start},{effective_stop}); "
        f"complete_coverage={complete_coverage}",
        flush=True,
    )
    # The generated pool and its canonicalization caches are no longer needed.
    # Releasing them before the six-dimensional loop is essential on the
    # shared host, whose swap may already be occupied by unrelated services.
    del raw_dimension_five
    canonical_cube_family.cache_clear()
    is_partial_cube_family.cache_clear()
    section.cache_clear()
    irredundant_projection.cache_clear()
    peripheral_decompositions.cache_clear()
    is_r_per.cache_clear()
    gc.collect()

    intersections = 0
    assignments = 0
    section_isometric = 0
    partial_cubes = 0
    ample = 0
    periphery_free: set[int] = set()
    minimal: set[int] = set()
    ample_labelled: set[int] = set()

    for local_index, base in enumerate(selected_bases, start=1):
        base_index = base_start + local_index
        members = vertices(base, 5)
        maximum_intersection = min(
            len(members), order_limit - len(members)
        )
        for intersection_order in range(1, maximum_intersection + 1):
            for chosen in itertools.combinations(
                members, intersection_order
            ):
                intersection = family_from_vertices(chosen)
                components = components_after_deletion(
                    base, 5, intersection
                )
                intersections += 1
                if len(components) < 2:
                    continue
                for assignment in range(1 << (len(components) - 1)):
                    zero = intersection | components[0]
                    one = intersection
                    for index, component in enumerate(components[1:]):
                        if assignment >> index & 1:
                            zero |= component
                        else:
                            one |= component
                    if one == intersection:
                        continue
                    assignments += 1
                    if not is_partial_cube_family(zero, 5):
                        continue
                    if not is_partial_cube_family(one, 5):
                        continue
                    section_isometric += 1
                    expansion = lift_expansion(zero, one, 5)
                    if not is_partial_cube_family(expansion, 6):
                        continue
                    partial_cubes += 1
                    if not shattering_profile(expansion, 6)[2]:
                        continue
                    ample += 1
                    if canonicalize_ample:
                        ample_labelled.add(expansion)
                    if peripheral_decompositions(expansion, 6):
                        continue
                    canonical = canonical_cube_family(expansion, 6)
                    periphery_free.add(canonical)
                    if all(
                        is_r_per(minor, minor_dimension)
                        for minor, minor_dimension, _
                        in immediate_pc_minors_with_operations(
                            expansion, 6
                        )
                    ):
                        minimal.add(canonical)
        # Candidate expansions from different bases share little recursive
        # state. Retaining their unbounded section/periphery caches caused the
        # earlier exact run to exhaust host memory. Clear them in small blocks;
        # the integer counters and canonical witness masks remain exact.
        if local_index % 10 == 0:
            is_partial_cube_family.cache_clear()
            section.cache_clear()
            irredundant_projection.cache_clear()
            peripheral_decompositions.cache_clear()
            is_r_per.cache_clear()
            canonical_cube_family.cache_clear()
            gc.collect()
        if local_index % 25 == 0:
            peak_rss_mb = (
                resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
            )
            print(
                f"  progress bases={base_index}/{total_bases}; "
                f"assignments={assignments}; ample={ample}; "
                f"periphery_free_types={len(periphery_free)}; "
                f"peak_rss_mb={peak_rss_mb:.1f}",
                flush=True,
            )

    print(
        f"  order_limit={order_limit}; "
        f"d4_ample_bases={len(dimension_four_bases)}; "
        f"d5_ample_types={total_bases}; "
        f"selected_bases={len(selected_bases)}; "
        f"complete_coverage={complete_coverage}"
    )
    print(
        f"  intersections={intersections}; assignments={assignments}; "
        f"isometric_section_pairs={section_isometric}; "
        f"partial_cube_expansions={partial_cubes}; "
        f"ample_expansions={ample}; "
        f"periphery_free_types={len(periphery_free)}; "
        f"minimal_obstruction_types={len(minimal)}"
    )
    if section_isometric != partial_cubes:
        raise AssertionError(
            "the expansion theorem cross-check failed: an isometric "
            "two-section cover did not lift to a partial cube"
        )
    ample_types: set[int] = set()
    if canonicalize_ample:
        print(
            "  canonicalizing ample dimension-six families: "
            f"labelled_masks={len(ample_labelled)}",
            flush=True,
        )
        for family_index, family in enumerate(
            sorted(ample_labelled), start=1
        ):
            ample_types.add(refined_canonical_cube_family(family, 6))
            if family_index % 500 == 0:
                gc.collect()
            if family_index % 1000 == 0:
                print(
                    f"  canonicalization progress={family_index}/"
                    f"{len(ample_labelled)}; types={len(ample_types)}",
                    flush=True,
                )
        gc.collect()
        print(
            f"  canonical ample dimension-six types={len(ample_types)}",
            flush=True,
        )
    print(
        "  peak_rss_mb="
        f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024:.1f}"
    )
    print(f"  runtime_seconds={time.perf_counter() - started:.3f}")
    if periphery_free:
        family = min(
            periphery_free, key=lambda item: (item.bit_count(), item)
        )
        print(
            "  smallest periphery-free type: "
            f"order={family.bit_count()}, "
            f"complete_crossing={not noncrossing_pairs(family, 6)}, "
            f"vertices={bitstrings(family, 6)}"
        )
    if minimal:
        family = min(minimal, key=lambda item: (item.bit_count(), item))
        print(
            "  smallest minimal obstruction type: "
            f"order={family.bit_count()}, "
            f"complete_crossing={not noncrossing_pairs(family, 6)}, "
            f"vertices={bitstrings(family, 6)}"
        )
    if complete_coverage and not periphery_free:
        print(
            "  exact conclusion: no periphery-free ample "
            "dimension-six type occurs in the exhausted order range"
        )
    return ample_types


def ample_dimension_six_base_census(order_limit: int = 14) -> set[int]:
    """Return every ample dimension-six type through an order bound.

    Unlike ``search_dimension_six_below_order``, this retains expansions whose
    newly introduced coordinate is peripheral.  Those graphs cannot be the
    current periphery-free target, but they are necessary contraction bases
    for the next-dimensional census.
    """

    if order_limit < 7 or order_limit > 17:
        raise ValueError("the dimension-six base limit must lie in [7,17]")

    started = time.perf_counter()
    dimension_four_bases = [
        family
        for family in enumerate_partial_cubes(4)[4]
        if family.bit_count() < order_limit
        and shattering_profile(family, 4)[2]
    ]
    raw_dimension_five = expansion_candidates(
        dimension_four_bases, 4, order_limit
    )
    dimension_five_bases = {
        refined_canonical_cube_family(family, 5)
        for family in raw_dimension_five
        if family.bit_count() < order_limit
        and shattering_profile(family, 5)[2]
    }
    del raw_dimension_five
    gc.collect()
    print("\nExact ample dimension-six contraction-base census", flush=True)
    print(
        f"  order_limit={order_limit}; "
        f"d4_ample_bases={len(dimension_four_bases)}; "
        f"d5_ample_types={len(dimension_five_bases)}",
        flush=True,
    )

    ample_labelled = {
        family
        for family in iter_expansion_candidates(
            dimension_five_bases, 5, order_limit
        )
        if shattering_profile(family, 6)[2]
    }
    gc.collect()
    print(
        f"  ample labelled masks={len(ample_labelled)}",
        flush=True,
    )
    ample_types: set[int] = set()
    for family_index, family in enumerate(
        sorted(ample_labelled), start=1
    ):
        ample_types.add(refined_canonical_cube_family(family, 6))
        if family_index % 500 == 0:
            gc.collect()
        if family_index % 2000 == 0:
            print(
                f"  canonicalization progress={family_index}/"
                f"{len(ample_labelled)}; types={len(ample_types)}",
                flush=True,
            )
    order_histogram = Counter(
        family.bit_count() for family in ample_types
    )
    print(
        f"  canonical ample dimension-six types={len(ample_types)}; "
        f"order_histogram={dict(sorted(order_histogram.items()))}; "
        "peak_rss_mb="
        f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024:.1f}; "
        f"runtime_seconds={time.perf_counter() - started:.3f}",
        flush=True,
    )
    return ample_types


def ample_type_census(
    target_dimension: int,
    order_limit: int,
) -> set[int]:
    """Enumerate all ample cube-isomorphism types in one bounded dimension.

    Starting from the exact dimension-four census, each layer is generated by
    all partial-cube expansions of the preceding contraction bases.  Types of
    order exactly ``order_limit`` are retained at the final layer but cannot be
    bases for a further nontrivial expansion.
    """

    if target_dimension < 4 or target_dimension > 8:
        raise ValueError("the generic ample census supports dimensions 4 through 8")
    # Dimension five through order twenty is the exact contraction source
    # needed by the constrained dimension-six/order-twenty-one campaign.
    # Higher-dimensional generic censuses retain their earlier guardrail.
    maximum_order = 20 if target_dimension == 5 else 17
    if order_limit < target_dimension + 1 or order_limit > maximum_order:
        raise ValueError("the ample census order bound is incompatible with the dimension")

    started = time.perf_counter()
    initial = {
        family
        for family in enumerate_partial_cubes(4)[4]
        if family.bit_count() <= order_limit
        and shattering_profile(family, 4)[2]
    }
    print("\nExact bounded ample-type census", flush=True)
    print(
        f"  dimension=4: types={len(initial)}; "
        f"order_limit={order_limit}",
        flush=True,
    )
    if target_dimension == 4:
        return initial

    bases = {
        family for family in initial if family.bit_count() < order_limit
    }
    current: set[int] = set()
    for dimension in range(5, target_dimension + 1):
        ample_labelled: set[int] = set()
        ordered_bases = sorted(bases)
        for base_index, base in enumerate(ordered_bases, start=1):
            for family in iter_expansion_candidates(
                (base,), dimension - 1, order_limit
            ):
                if shattering_profile(family, dimension)[2]:
                    ample_labelled.add(family)
            if base_index % 100 == 0:
                print(
                    f"  dimension={dimension} generation="
                    f"{base_index}/{len(ordered_bases)}; "
                    f"ample_labelled={len(ample_labelled)}",
                    flush=True,
                )
        gc.collect()
        current = set()
        for family_index, family in enumerate(
            sorted(ample_labelled), start=1
        ):
            current.add(
                refined_canonical_cube_family(family, dimension)
            )
            if family_index % 1000 == 0:
                gc.collect()
            if family_index % 5000 == 0:
                print(
                    f"  dimension={dimension} canonicalization="
                    f"{family_index}/{len(ample_labelled)}; "
                    f"types={len(current)}",
                    flush=True,
                )
        histogram = Counter(
            family.bit_count() for family in current
        )
        print(
            f"  dimension={dimension}: labelled={len(ample_labelled)}, "
            f"types={len(current)}, "
            f"order_histogram={dict(sorted(histogram.items()))}",
            flush=True,
        )
        del ample_labelled
        gc.collect()
        bases = {
            family
            for family in current
            if family.bit_count() < order_limit
        }
    print(
        "  peak_rss_mb="
        f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024:.1f}; "
        f"runtime_seconds={time.perf_counter() - started:.3f}",
        flush=True,
    )
    return current


def search_dimension_seven_below_order(
    order_limit: int = 15,
    base_start: int = 0,
    base_stop: int | None = None,
) -> None:
    """Exhaust ample dimension-seven candidates through an order bound.

    A dimension-seven candidate of order at most ``N`` contracts to an ample
    dimension-six graph of order at most ``N-1``.  The preceding exact base
    census supplies every such type, including bases obtained through a
    peripheral last coordinate.
    """

    if order_limit < 8 or order_limit > 16:
        raise ValueError("the dimension-seven order limit must lie in [8,16]")
    if base_start < 0:
        raise ValueError("the first dimension-six base index must be nonnegative")

    started = time.perf_counter()
    dimension_six_bases = sorted(
        ample_dimension_six_base_census(order_limit - 1)
    )
    total_bases = len(dimension_six_bases)
    effective_stop = (
        total_bases
        if base_stop is None
        else min(base_stop, total_bases)
    )
    if base_start > effective_stop:
        raise ValueError("the dimension-seven base range is empty or reversed")
    selected_bases = dimension_six_bases[base_start:effective_stop]
    complete_coverage = base_start == 0 and effective_stop == total_bases

    print("\nExact dimension-seven ample expansion census", flush=True)
    print(
        f"  order_limit={order_limit}; d6_ample_types={total_bases}; "
        f"selected_range=[{base_start},{effective_stop}); "
        f"complete_coverage={complete_coverage}",
        flush=True,
    )

    intersections = 0
    assignments = 0
    section_isometric = 0
    partial_cubes = 0
    ample = 0
    periphery_free: set[int] = set()
    minimal: set[int] = set()

    for local_index, base in enumerate(selected_bases, start=1):
        base_index = base_start + local_index
        members = vertices(base, 6)
        maximum_intersection = min(
            len(members), order_limit - len(members)
        )
        for intersection_order in range(1, maximum_intersection + 1):
            for chosen in itertools.combinations(
                members, intersection_order
            ):
                intersection = family_from_vertices(chosen)
                components = components_after_deletion(
                    base, 6, intersection
                )
                intersections += 1
                if len(components) < 2:
                    continue
                for assignment in range(1 << (len(components) - 1)):
                    zero = intersection | components[0]
                    one = intersection
                    for index, component in enumerate(components[1:]):
                        if assignment >> index & 1:
                            zero |= component
                        else:
                            one |= component
                    if one == intersection:
                        continue
                    assignments += 1
                    if not is_partial_cube_family(zero, 6):
                        continue
                    if not is_partial_cube_family(one, 6):
                        continue
                    section_isometric += 1
                    expansion = lift_expansion(zero, one, 6)
                    if not is_partial_cube_family(expansion, 7):
                        continue
                    partial_cubes += 1
                    if not shattering_profile(expansion, 7)[2]:
                        continue
                    ample += 1
                    if peripheral_decompositions(expansion, 7):
                        continue
                    canonical = refined_canonical_cube_family(
                        expansion, 7
                    )
                    periphery_free.add(canonical)
                    if all(
                        is_r_per(minor, minor_dimension)
                        for minor, minor_dimension, _
                        in immediate_pc_minors_with_operations(
                            expansion, 7
                        )
                    ):
                        minimal.add(canonical)

        if local_index % 10 == 0:
            is_partial_cube_family.cache_clear()
            section.cache_clear()
            irredundant_projection.cache_clear()
            peripheral_decompositions.cache_clear()
            is_r_per.cache_clear()
            canonical_cube_family.cache_clear()
            gc.collect()
        if local_index % 100 == 0:
            print(
                f"  progress bases={base_index}/{total_bases}; "
                f"assignments={assignments}; ample={ample}; "
                f"periphery_free_types={len(periphery_free)}; "
                "peak_rss_mb="
                f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024:.1f}",
                flush=True,
            )

    print(
        f"  intersections={intersections}; assignments={assignments}; "
        f"isometric_section_pairs={section_isometric}; "
        f"partial_cube_expansions={partial_cubes}; "
        f"ample_expansions={ample}; "
        f"periphery_free_types={len(periphery_free)}; "
        f"minimal_obstruction_types={len(minimal)}"
    )
    if section_isometric != partial_cubes:
        raise AssertionError(
            "the dimension-seven expansion-theorem cross-check failed"
        )
    print(
        "  peak_rss_mb="
        f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024:.1f}; "
        f"runtime_seconds={time.perf_counter() - started:.3f}"
    )
    for family in sorted(
        periphery_free, key=lambda item: (item.bit_count(), item)
    ):
        print(
            "  periphery-free type: "
            f"order={family.bit_count()}, "
            f"minimal={family in minimal}, "
            f"vertices={bitstrings(family, 7)}"
        )
    if complete_coverage and not periphery_free:
        print(
            "  exact conclusion: no periphery-free ample "
            "dimension-seven type occurs in the exhausted order range"
        )


def rooted_factor_census(
    max_dimension: int,
    max_wedge_dimension: int,
) -> None:
    """Classify small pointed factors and test the rooted-wedge criterion.

    The wedge criterion tested here is

        A wedge B is in R_per
        iff A or B has a root-avoiding recursive certificate.

    The pointed input types are exact up to rooted cube isomorphism because
    ``enumerate_partial_cubes`` is exact up to ordinary cube isomorphism and
    every vertex of each representative is considered as a root.
    """

    if max_dimension < 0 or max_dimension > 4:
        raise ValueError("the rooted-factor census supports dimensions 0 through 4")
    if max_wedge_dimension < 0 or max_wedge_dimension > 8:
        raise ValueError("the wedge-test dimension must lie in [0,8]")

    started = time.perf_counter()
    census = enumerate_partial_cubes(max_dimension)
    r_types: list[tuple[int, int]] = []
    critical_types: list[tuple[int, int]] = []

    print("\nExact rooted-factor census")
    for dimension in range(max_dimension + 1):
        pointed = {
            canonical_pointed_family(family, dimension, root)
            for family in census[dimension]
            for root in vertices(family, dimension)
        }
        members = {
            family
            for family in pointed
            if is_r_per(family, dimension)
        }
        rooted = {
            family
            for family in members
            if is_rooted_r_per(family, dimension, 0)
        }
        root_last_peelable = {
            family
            for family in members
            if has_corner_peeling_ending_at_root(
                family,
                dimension,
                0,
            )
        }
        explicit_root_last_orders = {
            family: corner_peeling_order_ending_at_root(
                family,
                dimension,
                0,
            )
            for family in members
        }
        explicit_root_last_peelable = {
            family
            for family, order in explicit_root_last_orders.items()
            if order is not None
        }
        if root_last_peelable != explicit_root_last_peelable:
            raise AssertionError(
                "projected and original-coordinate root-last "
                "corner-peeling tests disagree"
            )
        locked = {
            family
            for family in members
            if is_root_locked(family, dimension, 0)
        }
        critical = {
            family
            for family in locked
            if is_root_critical(family, dimension, 0)
        }
        r_types.extend((family, dimension) for family in sorted(members))
        critical_types.extend(
            (family, dimension) for family in sorted(critical)
        )
        print(
            f"  dimension={dimension}: pointed_types={len(pointed)}, "
            f"R_per={len(members)}, root_avoiding={len(rooted)}, "
            "root_last_corner_peelable="
            f"{len(root_last_peelable)}, "
            f"root_locked={len(locked)}, root_critical={len(critical)}"
        )
        rooted_not_peelable = sorted(rooted - root_last_peelable)
        peelable_not_rooted = sorted(root_last_peelable - rooted)
        if rooted_not_peelable or peelable_not_rooted:
            print(
                "    root-avoiding/corner-peeling differences: "
                f"rooted_not_peelable={len(rooted_not_peelable)}, "
                f"peelable_not_rooted={len(peelable_not_rooted)}"
            )
        if rooted_not_peelable:
            print(
                "      first rooted_not_peelable: "
                f"order={rooted_not_peelable[0].bit_count()}, "
                "vertices="
                f"{bitstrings(rooted_not_peelable[0], dimension)}"
            )
        if peelable_not_rooted:
            witness = peelable_not_rooted[0]
            peeling = explicit_root_last_orders[witness]
            face_certificate = corner_peeling_face_certificate(
                witness,
                dimension,
                peeling or (),
            )
            print(
                "      first peelable_not_rooted: "
                f"order={witness.bit_count()}, "
                "vertices="
                f"{bitstrings(witness, dimension)}, "
                "root-last deletion_order="
                f"{tuple(format(vertex, f'0{dimension}b') for vertex in peeling or ())}, "
                f"unique_maximal_faces={face_certificate}"
            )
        for family in sorted(
            critical, key=lambda item: (item.bit_count(), item)
        ):
            print(
                "    critical: "
                f"order={family.bit_count()}, "
                f"vertices={bitstrings(family, dimension)}"
            )

    tested_pairs = 0
    rooted_flags = {
        item: is_rooted_r_per(item[0], item[1], 0)
        for item in r_types
    }
    critical_flags = {
        item: is_root_critical(item[0], item[1], 0)
        for item in r_types
    }
    membership_failures: list[
        tuple[int, int, int, int, bool, bool]
    ] = []
    minimality_failures: list[
        tuple[int, int, int, int, bool, bool]
    ] = []
    for left_index, (left, left_dimension) in enumerate(r_types):
        for right, right_dimension in r_types[left_index:]:
            if left_dimension + right_dimension > max_wedge_dimension:
                continue
            wedge, wedge_dimension = one_vertex_wedge(
                left,
                left_dimension,
                0,
                right,
                right_dimension,
                0,
            )
            tested_pairs += 1
            predicted_membership = (
                rooted_flags[(left, left_dimension)]
                or rooted_flags[(right, right_dimension)]
            )
            actual_membership = is_r_per(wedge, wedge_dimension)
            if predicted_membership != actual_membership:
                membership_failures.append(
                    (
                        left,
                        left_dimension,
                        right,
                        right_dimension,
                        predicted_membership,
                        actual_membership,
                    )
                )

            predicted_minimal = (
                critical_flags[(left, left_dimension)]
                and critical_flags[(right, right_dimension)]
            )
            actual_minimal = (
                not actual_membership
                and all(
                    is_r_per(minor, minor_dimension)
                    for minor, minor_dimension, _
                    in immediate_pc_minors_with_operations(
                        wedge, wedge_dimension
                    )
                )
            )
            if predicted_minimal != actual_minimal:
                minimality_failures.append(
                    (
                        left,
                        left_dimension,
                        right,
                        right_dimension,
                        predicted_minimal,
                        actual_minimal,
                    )
                )
        # Recursive certificates of distinct eight-dimensional wedges share
        # little state.  Clear only those unbounded caches periodically so the
        # exact regression remains safe on the shared small-memory host.
        if left_index % 10 == 9:
            is_partial_cube_family.cache_clear()
            section.cache_clear()
            irredundant_projection.cache_clear()
            peripheral_decompositions.cache_clear()
            is_r_per.cache_clear()
            is_rooted_r_per.cache_clear()
            canonical_cube_family.cache_clear()
        if left_index and left_index % 50 == 0:
            print(
                f"  wedge progress: left_types={left_index}/"
                f"{len(r_types)}, tested_pairs={tested_pairs}",
                flush=True,
            )

    print(
        "  rooted-wedge criterion: "
        f"tested_pairs={tested_pairs}, "
        f"max_total_dimension={max_wedge_dimension}, "
        f"membership_failures={len(membership_failures)}, "
        f"minimality_failures={len(minimality_failures)}"
    )
    if membership_failures:
        print(f"  first membership failure={membership_failures[0]}")
    if minimality_failures:
        print(f"  first minimality failure={minimality_failures[0]}")
    if membership_failures or minimality_failures:
        raise AssertionError("the rooted-wedge characterization regression failed")
    if max_dimension == 4:
        expected = (
            canonical_pointed_family(root_locked_factor(), 4, 0),
            4,
        )
        if critical_types != [expected]:
            raise AssertionError(
                "L4 must be the unique root-critical pointed type "
                "through dimension four"
            )
    print(
        f"  critical_pointed_types={len(critical_types)}; "
        f"runtime_seconds={time.perf_counter() - started:.3f}"
    )


def rooted_synchronization_regression(max_dimension: int) -> None:
    """Check elimination and reverse-pair criteria on all small pointed types."""

    if max_dimension < 0 or max_dimension > 4:
        raise ValueError(
            "the rooted synchronization regression supports dimensions 0 through 4"
        )

    started = time.perf_counter()
    census = enumerate_partial_cubes(max_dimension)
    elimination_checks = 0
    reverse_presentations = 0
    critical_pair_presentations = 0
    elimination_failures: list[tuple[int, int, bool, bool]] = []
    synchronization_failures: list[
        tuple[int, int, int, bool, bool]
    ] = []
    closure_failures: list[
        tuple[int, int, int, bool, bool]
    ] = []
    antimatroid_failures: list[
        tuple[int, int, int, int, int, bool, bool]
    ] = []
    criticality_failures: list[
        tuple[int, int, int, bool, bool]
    ] = []
    closure_criticality_failures: list[
        tuple[int, int, int, bool, bool]
    ] = []

    print("\nExact rooted synchronization regression")
    for dimension in range(max_dimension + 1):
        pointed = {
            canonical_pointed_family(family, dimension, root)
            for family in census[dimension]
            for root in vertices(family, dimension)
        }
        members = sorted(
            family
            for family in pointed
            if is_r_per(family, dimension)
        )
        dimension_presentations = 0
        dimension_critical_pairs = 0
        for family in members:
            rooted_actual = is_rooted_r_per(family, dimension, 0)
            rooted_predicted = has_root_elimination_order(
                family,
                dimension,
                0,
            )
            elimination_checks += 1
            if rooted_actual != rooted_predicted:
                elimination_failures.append(
                    (
                        family,
                        dimension,
                        rooted_predicted,
                        rooted_actual,
                    )
                )

            critical_actual = is_root_critical(
                family,
                dimension,
                0,
            )
            for coordinate in range(dimension):
                root_side = section(
                    family,
                    dimension,
                    coordinate,
                    0,
                )
                opposite_side = section(
                    family,
                    dimension,
                    coordinate,
                    1,
                )
                if not nested(root_side, opposite_side):
                    continue
                dimension_presentations += 1
                reverse_presentations += 1
                synchronized = is_synchronizable_reverse_pair(
                    opposite_side,
                    root_side,
                    dimension - 1,
                    0,
                )
                closure_synchronized = (
                    is_synchronizable_by_common_closure(
                        opposite_side,
                        root_side,
                        dimension - 1,
                        0,
                    )
                )
                feasible_sets = common_availability_antimatroid(
                    opposite_side,
                    root_side,
                    dimension - 1,
                    0,
                )
                feasible_union = 0
                for feasible_set in feasible_sets:
                    feasible_union |= feasible_set
                closure_mask = common_availability_closure(
                    opposite_side,
                    root_side,
                    dimension - 1,
                    0,
                )[0]
                antimatroid_synchronized = (
                    (1 << (dimension - 1)) - 1
                ) in feasible_sets
                if (
                    feasible_union != closure_mask
                    or antimatroid_synchronized != synchronized
                ):
                    antimatroid_failures.append(
                        (
                            family,
                            dimension,
                            coordinate,
                            feasible_union,
                            closure_mask,
                            antimatroid_synchronized,
                            synchronized,
                        )
                    )
                if synchronized != rooted_actual:
                    synchronization_failures.append(
                        (
                            family,
                            dimension,
                            coordinate,
                            synchronized,
                            rooted_actual,
                        )
                    )
                if closure_synchronized != synchronized:
                    closure_failures.append(
                        (
                            family,
                            dimension,
                            coordinate,
                            closure_synchronized,
                            synchronized,
                        )
                    )
                critical_predicted = is_critical_reverse_pair(
                    opposite_side,
                    root_side,
                    dimension - 1,
                    0,
                )
                closure_critical_predicted = (
                    is_critical_reverse_pair_by_closure(
                        opposite_side,
                        root_side,
                        dimension - 1,
                        0,
                    )
                )
                if critical_predicted:
                    dimension_critical_pairs += 1
                    critical_pair_presentations += 1
                if critical_predicted != critical_actual:
                    criticality_failures.append(
                        (
                            family,
                            dimension,
                            coordinate,
                            critical_predicted,
                            critical_actual,
                        )
                    )
                if closure_critical_predicted != critical_predicted:
                    closure_criticality_failures.append(
                        (
                            family,
                            dimension,
                            coordinate,
                            closure_critical_predicted,
                            critical_predicted,
                        )
                    )
        print(
            f"  dimension={dimension}: R_per_pointed={len(members)}, "
            f"reverse_presentations={dimension_presentations}, "
            f"critical_pair_presentations={dimension_critical_pairs}"
        )

    print(
        f"  elimination_checks={elimination_checks}; "
        f"reverse_presentations={reverse_presentations}; "
        f"critical_pair_presentations={critical_pair_presentations}; "
        f"elimination_failures={len(elimination_failures)}; "
        "synchronization_failures="
        f"{len(synchronization_failures)}; "
        f"closure_failures={len(closure_failures)}; "
        f"antimatroid_failures={len(antimatroid_failures)}; "
        f"criticality_failures={len(criticality_failures)}; "
        "closure_criticality_failures="
        f"{len(closure_criticality_failures)}; "
        f"runtime_seconds={time.perf_counter() - started:.3f}"
    )
    failures = (
        elimination_failures
        or synchronization_failures
        or closure_failures
        or antimatroid_failures
        or criticality_failures
        or closure_criticality_failures
    )
    if failures:
        print(f"  first failure={failures[0]}")
        raise AssertionError(
            "the rooted synchronization characterization regression failed"
        )


def repair_signature_payload(
    signature: tuple[
        int,
        int,
        tuple[tuple[int, bool, int, bool], ...],
    ],
) -> dict[str, object]:
    """Convert a canonical repair signature to JSON data."""

    base_available, site_available, rows = signature
    dimension = len(rows)

    def coordinates(mask: int) -> list[int]:
        return [
            coordinate
            for coordinate in range(dimension)
            if mask >> coordinate & 1
        ]

    return {
        "base_available": coordinates(base_available),
        "site_available": coordinates(site_available),
        "forced_coordinates": [
            {
                "coordinate": coordinate,
                "contraction_first_moves": coordinates(row[0]),
                "contraction_terminal_equality": row[1],
                "root_restriction_first_moves": coordinates(row[2]),
                "root_restriction_terminal_equality": row[3],
            }
            for coordinate, row in enumerate(rows)
        ],
    }


def activation_antimatroid_signature_payload(
    signature: tuple[
        tuple[int, ...],
        tuple[tuple[tuple[int, ...], tuple[int, ...]], ...],
    ],
) -> dict[str, object]:
    """Convert canonical feasible contraction sets to JSON."""

    initial_feasible, rows = signature
    dimension = len(rows)

    def coordinates(mask: int) -> list[int]:
        return [
            coordinate
            for coordinate in range(dimension)
            if mask >> coordinate & 1
        ]

    def feasible_sets(masks: tuple[int, ...]) -> list[list[int]]:
        return [coordinates(mask) for mask in masks]

    return {
        "initial_feasible_sets": feasible_sets(initial_feasible),
        "forced_coordinates": [
            {
                "coordinate": coordinate,
                "contraction_feasible_sets": feasible_sets(row[0]),
                "root_restriction_feasible_sets": feasible_sets(row[1]),
            }
            for coordinate, row in enumerate(rows)
        ],
    }


@lru_cache(maxsize=None)
def rooted_r_per_families_in_cube(dimension: int) -> tuple[int, ...]:
    """Enumerate every rooted-good family in ``Q_dimension`` with root zero.

    The exact census supplies one cube-isomorphism representative in every
    intrinsic dimension at most ``dimension``.  Choosing the root, translating
    it to zero, and injecting the active coordinates into the ambient cube
    exhausts all labelled rooted embeddings; a set removes automorphic
    duplicates.
    """

    if dimension < 0 or dimension > 4:
        raise ValueError(
            "the complete rooted-family generator supports dimensions 0 through 4"
        )

    census = enumerate_partial_cubes(dimension)
    embedded: set[int] = set()
    for intrinsic_dimension in range(dimension + 1):
        for family in census[intrinsic_dimension]:
            if not is_r_per(family, intrinsic_dimension):
                continue
            for root in vertices(family, intrinsic_dimension):
                translated = tuple(
                    vertex ^ root
                    for vertex in vertices(family, intrinsic_dimension)
                )
                for coordinate_image in itertools.permutations(
                    range(dimension),
                    intrinsic_dimension,
                ):
                    image = family_from_vertices(
                        sum(
                            (
                                (vertex >> old_coordinate) & 1
                            )
                            << coordinate_image[old_coordinate]
                            for old_coordinate in range(intrinsic_dimension)
                        )
                        for vertex in translated
                    )
                    if is_rooted_r_per(image, dimension, 0):
                        embedded.add(image)

    return tuple(sorted(embedded, key=lambda family: (family.bit_count(), family)))


@lru_cache(maxsize=None)
def r_per_families_in_cube(dimension: int) -> tuple[int, ...]:
    """Enumerate every labelled ``R_per`` family embedded in a small cube."""

    if dimension < 0 or dimension > 4:
        raise ValueError(
            "the complete labelled-family generator supports dimensions 0 through 4"
        )

    census = enumerate_partial_cubes(dimension)
    embedded: set[int] = set()
    for intrinsic_dimension in range(dimension + 1):
        for family in census[intrinsic_dimension]:
            if not is_r_per(family, intrinsic_dimension):
                continue
            members = vertices(family, intrinsic_dimension)
            for coordinate_image in itertools.permutations(
                range(dimension),
                intrinsic_dimension,
            ):
                unflipped = tuple(
                    sum(
                        ((vertex >> old_coordinate) & 1)
                        << coordinate_image[old_coordinate]
                        for old_coordinate in range(intrinsic_dimension)
                    )
                    for vertex in members
                )
                for translation in range(1 << dimension):
                    embedded.add(
                        family_from_vertices(
                            vertex ^ translation
                            for vertex in unflipped
                        )
                    )

    return tuple(sorted(embedded, key=lambda family: (family.bit_count(), family)))


@lru_cache(maxsize=None)
def rooted_r_per_families_dimension_five(
    max_order: int = 32,
) -> tuple[tuple[int, ...], int]:
    """Generate every rooted-good family in ``Q_5`` with root zero.

    Every nontrivial rooted-good family has an active root-avoiding peripheral
    coordinate.  Contracting it gives a rooted-good base in ``Q_4`` and its
    opposite section is an ordinary ``R_per`` subfamily of that base.
    Conversely every such contained pair gives a rooted-good peripheral lift.
    Since the four old coordinates are already fully labelled, inserting the
    new coordinate in each of the five positions exhausts all labelled lifts.
    """

    if max_order < 1 or max_order > 32:
        raise ValueError("the rooted Q5 order bound must lie in [1,32]")

    base_dimension = 4
    rooted_bases = rooted_r_per_families_in_cube(base_dimension)
    sites = r_per_families_in_cube(base_dimension)
    generated: set[int] = {1}
    contained_presentations = 0
    for base_index, base in enumerate(rooted_bases, start=1):
        for site in sites:
            if (
                site & ~base
                or base.bit_count() + site.bit_count() > max_order
            ):
                continue
            contained_presentations += 1
            lift = lift_expansion(base, site, base_dimension)
            for new_coordinate in range(base_dimension + 1):
                order = (
                    tuple(range(new_coordinate))
                    + (base_dimension,)
                    + tuple(range(new_coordinate, base_dimension))
                )
                generated.add(
                    permute_family_by_coordinate_order(
                        lift,
                        base_dimension + 1,
                        order,
                    )
                )
        if base_index % 500 == 0:
            print(
                "  rooted Q5 generator progress: "
                f"bases={base_index}/{len(rooted_bases)}; "
                f"contained_presentations={contained_presentations}; "
                f"distinct_families={len(generated)}",
                flush=True,
            )

    return (
        tuple(
            sorted(
                generated,
                key=lambda family: (family.bit_count(), family),
            )
        ),
        contained_presentations,
    )


def rooted_family_dimension_five_inventory(output_path: str) -> None:
    """Write a compact inventory for the exact rooted-good ``Q_5`` generator."""

    started = time.perf_counter()
    rooted_bases = rooted_r_per_families_in_cube(4)
    sites = r_per_families_in_cube(4)
    print(
        "\nExact rooted-good Q5 generator",
        f"  rooted Q4 bases={len(rooted_bases)}; "
        f"labelled R_per Q4 sites={len(sites)}",
        sep="\n",
        flush=True,
    )
    families, contained_presentations = (
        rooted_r_per_families_dimension_five(32)
    )
    irredundant = [
        family
        for family in families
        if irredundant_projection(family, 5)[1] == 5
    ]
    order_histogram = Counter(
        family.bit_count() for family in families
    )
    irredundant_order_histogram = Counter(
        family.bit_count() for family in irredundant
    )
    family_digest = hashlib.sha256(
        (
            "\n".join(format(family, "x") for family in families)
            + "\n"
        ).encode("ascii")
    ).hexdigest()
    payload = {
        "schema_version": 1,
        "kind": "rooted-r-per-q5-inventory",
        "complete": True,
        "dimension": 5,
        "root": 0,
        "rooted_q4_bases": len(rooted_bases),
        "labelled_r_per_q4_sites": len(sites),
        "contained_peripheral_presentations": contained_presentations,
        "rooted_r_per_labelled_families": len(families),
        "irredundant_rooted_r_per_labelled_families": len(irredundant),
        "order_histogram": dict(sorted(order_histogram.items())),
        "irredundant_order_histogram": dict(
            sorted(irredundant_order_histogram.items())
        ),
        "family_digest": family_digest,
    }
    Path(output_path).write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        f"  contained_presentations={contained_presentations}; "
        f"rooted_families={len(families)}; "
        f"irredundant_rooted_families={len(irredundant)}"
    )
    print(f"  order_histogram={dict(sorted(order_histogram.items()))}")
    print(
        "  irredundant_order_histogram="
        f"{dict(sorted(irredundant_order_histogram.items()))}"
    )
    print(f"  family_digest={family_digest}")
    print(
        "  peak_rss_mb="
        f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024:.1f}; "
        f"runtime_seconds={time.perf_counter() - started:.3f}"
    )
    print(f"  wrote={output_path}")


def rooted_critical_pair_profile_dimension_five(
    factor_order_limit: int,
    output_path: str,
) -> None:
    """Search critical reverse pairs on five old coordinates exactly by order.

    A critical pair has a site of order at least three.  A one-vertex site
    has every ambient coordinate available.  A two-vertex rooted partial
    cube is an edge: its varying coordinate has equal projected sections,
    while every constant coordinate has empty opposite section, so again
    every ambient coordinate is available.  Either case intersects the
    nonempty availability set of a nontrivial rooted-good base.  Thus a lift
    of order at most ``factor_order_limit`` has base order at most
    ``factor_order_limit - 3``.
    """

    if factor_order_limit < 8 or factor_order_limit > 17:
        raise ValueError(
            "the exact five-coordinate pair search supports lift orders 8 through 17"
        )

    started = time.perf_counter()
    dimension = 5
    families, generator_presentations = (
        rooted_r_per_families_dimension_five(
            factor_order_limit - 3
        )
    )
    universe_mask = (1 << (1 << dimension)) - 1
    low_masks = coordinate_low_vertex_masks(dimension)

    labelled_base_candidates = [
        family
        for family in families
        if 6 <= family.bit_count() <= factor_order_limit - 3
        and all(
            family & low_mask
            and family & (universe_mask ^ low_mask)
            for low_mask in low_masks
        )
    ]
    canonical_bases: set[int] = set()
    for index, family in enumerate(labelled_base_candidates, start=1):
        canonical_bases.add(
            canonical_pointed_family(family, dimension, 0)
        )
        if index % 25_000 == 0:
            print(
                "  five-coordinate base quotient progress: "
                f"labelled={index}/{len(labelled_base_candidates)}; "
                f"canonical={len(canonical_bases)}",
                flush=True,
            )
    ordered_bases = sorted(
        canonical_bases,
        key=lambda family: (family.bit_count(), family),
    )

    maximum_site_order = factor_order_limit - 6
    rooted_sites = {
        family
        for family in families
        if 3 <= family.bit_count() <= maximum_site_order
    }
    availability_cache: dict[int, int] = {}

    def availability_mask(family: int) -> int:
        answer = availability_cache.get(family)
        if answer is None:
            answer = sum(
                1 << coordinate
                for coordinate in (
                    raw_root_avoiding_peripheral_coordinates(
                        family,
                        dimension,
                        0,
                    )
                )
            )
            availability_cache[family] = answer
        return answer

    subset_positions_tested = 0
    rooted_containments = 0
    availability_disjoint_containments = 0
    availability_disjoint_pairs: set[tuple[int, int]] = set()
    critical_presentations = 0
    critical_pairs: set[tuple[int, int]] = set()
    for base_index, base in enumerate(ordered_bases, start=1):
        base_order = base.bit_count()
        base_available = availability_mask(base)
        site = (base - 1) & base
        while site:
            if site & 1:
                site_order = site.bit_count()
                if (
                    3 <= site_order
                    and base_order + site_order <= factor_order_limit
                ):
                    subset_positions_tested += 1
                    if site in rooted_sites:
                        rooted_containments += 1
                        if not (
                            base_available
                            & availability_mask(site)
                        ):
                            availability_disjoint_containments += 1
                            pair = canonical_reverse_pair(
                                base,
                                site,
                                dimension,
                            )
                            availability_disjoint_pairs.add(pair)
                            if all(
                                is_synchronizable_by_common_closure(
                                    section(
                                        base,
                                        dimension,
                                        coordinate,
                                        0,
                                    ),
                                    section(
                                        site,
                                        dimension,
                                        coordinate,
                                        0,
                                    ),
                                    dimension - 1,
                                    0,
                                )
                                and is_synchronizable_by_common_closure(
                                    contract_coordinate(
                                        base,
                                        dimension,
                                        coordinate,
                                    ),
                                    contract_coordinate(
                                        site,
                                        dimension,
                                        coordinate,
                                    ),
                                    dimension - 1,
                                    0,
                                )
                                for coordinate in range(dimension)
                            ):
                                critical_presentations += 1
                                if not (
                                    is_critical_reverse_pair(
                                        *pair,
                                        dimension,
                                        0,
                                    )
                                    and is_critical_reverse_pair_by_closure(
                                        *pair,
                                        dimension,
                                        0,
                                    )
                                ):
                                    raise AssertionError(
                                        "the bounded direct search disagrees "
                                        "with a critical-pair verifier"
                                    )
                                factor = lift_expansion(
                                    pair[1],
                                    pair[0],
                                    dimension,
                                )
                                if not is_root_critical(
                                    factor,
                                    dimension + 1,
                                    0,
                                ):
                                    raise AssertionError(
                                        "a five-coordinate critical pair did "
                                        "not lift to a root-critical factor"
                                    )
                                critical_pairs.add(pair)
            site = (site - 1) & base
        if base_index % 250 == 0:
            print(
                "  five-coordinate pair search progress: "
                f"bases={base_index}/{len(ordered_bases)}; "
                f"subset_positions={subset_positions_tested}; "
                f"rooted_containments={rooted_containments}; "
                "availability_disjoint="
                f"{availability_disjoint_containments}; "
                f"critical_types={len(critical_pairs)}",
                flush=True,
            )

    pairs = sorted(critical_pairs)
    disjoint_pairs = sorted(availability_disjoint_pairs)
    pair_index = {pair: index for index, pair in enumerate(pairs)}
    repair_by_pair = {
        pair: canonical_repair_signature(*pair, dimension)
        for pair in pairs
    }
    activation_by_pair = {
        pair: canonical_activation_antimatroid_signature(
            *pair,
            dimension,
        )
        for pair in pairs
    }

    def collision_classes(
        signatures: dict[tuple[int, int], object],
    ) -> tuple[int, list[list[tuple[int, int]]]]:
        classes: dict[object, list[tuple[int, int]]] = {}
        for pair in pairs:
            classes.setdefault(
                (
                    pair[0].bit_count(),
                    pair[1].bit_count(),
                    signatures[pair],
                ),
                [],
            ).append(pair)
        collisions = [
            members
            for members in classes.values()
            if len(members) > 1
        ]
        return len(classes), sorted(
            collisions,
            key=lambda members: (
                -len(members),
                tuple(pair_index[pair] for pair in members),
            ),
        )

    ordered_repair_types, ordered_repair_collisions = (
        collision_classes(repair_by_pair)
    )
    ordered_activation_types, ordered_activation_collisions = (
        collision_classes(activation_by_pair)
    )
    pair_digest = hashlib.sha256(
        (
            "\n".join(f"{base:x}:{site:x}" for base, site in pairs)
            + "\n"
        ).encode("ascii")
    ).hexdigest()

    def coordinates(mask: int, labels: tuple[int, ...]) -> list[int]:
        return [
            label
            for coordinate, label in enumerate(labels)
            if mask >> coordinate & 1
        ]

    def forced_minor_closure_profile(
        pair: tuple[int, int],
    ) -> list[dict[str, object]]:
        base, site = pair
        rows: list[dict[str, object]] = []
        for coordinate in range(dimension):
            labels = tuple(
                old_coordinate
                for old_coordinate in range(dimension)
                if old_coordinate != coordinate
            )
            contraction_mask = common_availability_closure(
                contract_coordinate(base, dimension, coordinate),
                contract_coordinate(site, dimension, coordinate),
                dimension - 1,
                0,
            )[0]
            restriction_mask = common_availability_closure(
                section(base, dimension, coordinate, 0),
                section(site, dimension, coordinate, 0),
                dimension - 1,
                0,
            )[0]
            full_mask = (1 << (dimension - 1)) - 1
            rows.append(
                {
                    "coordinate": coordinate,
                    "contraction_closure": coordinates(
                        contraction_mask,
                        labels,
                    ),
                    "contraction_full_rank": (
                        contraction_mask == full_mask
                    ),
                    "root_restriction_closure": coordinates(
                        restriction_mask,
                        labels,
                    ),
                    "root_restriction_full_rank": (
                        restriction_mask == full_mask
                    ),
                }
            )
        return rows

    availability_disjoint_pair_digest = hashlib.sha256(
        (
            "\n".join(
                f"{base:x}:{site:x}"
                for base, site in disjoint_pairs
            )
            + "\n"
        ).encode("ascii")
    ).hexdigest()

    def critical_pair_payload(
        pair: tuple[int, int],
    ) -> dict[str, object]:
        factor = lift_expansion(
            pair[1],
            pair[0],
            dimension,
        )
        factor_dimension = dimension + 1
        wedge, wedge_dimension = one_vertex_wedge(
            root_locked_factor(),
            4,
            0,
            factor,
            factor_dimension,
            0,
        )
        wedge_report = verify_obstruction(wedge, wedge_dimension)
        if not (
            wedge_report["periphery_free"]
            and not wedge_report["R_per"]
            and wedge_report["all_immediate_minors_R_per"]
        ):
            raise AssertionError(
                "a five-coordinate critical pair produced a nonminimal wedge"
            )
        return {
            "index": pair_index[pair],
            "base_mask": format(pair[0], "x"),
            "site_mask": format(pair[1], "x"),
            "base_order": pair[0].bit_count(),
            "site_order": pair[1].bit_count(),
            "factor_order": factor.bit_count(),
            "factor_profile": {
                "dimension": factor_dimension,
                "mask": format(factor, "x"),
                "vertices": [
                    format(vertex, f"0{factor_dimension}b")
                    for vertex in vertices(
                        factor,
                        factor_dimension,
                    )
                ],
                "R_per": is_r_per(factor, factor_dimension),
                "root_locked": is_root_locked(
                    factor,
                    factor_dimension,
                    0,
                ),
                "root_critical": is_root_critical(
                    factor,
                    factor_dimension,
                    0,
                ),
                "root_last_corner_peelable": (
                    has_corner_peeling_ending_at_root(
                        factor,
                        factor_dimension,
                        0,
                    )
                ),
                "peripheral_coordinates": [
                    coordinate
                    for coordinate, _, _, _
                    in peripheral_decompositions(
                        factor,
                        factor_dimension,
                    )
                ],
                "shattered_rank_histogram": (
                    shattered_rank_histogram(
                        factor,
                        factor_dimension,
                    )
                ),
            },
            "wedge_with_L4_profile": wedge_report,
            "canonical_repair_signature": repair_signature_payload(
                repair_by_pair[pair]
            ),
            "canonical_activation_antimatroid": (
                activation_antimatroid_signature_payload(
                    activation_by_pair[pair]
                )
            ),
        }

    payload = {
        "schema_version": 1,
        "kind": "rooted-critical-pair-profile-d5",
        "complete_within_factor_order": True,
        "old_dimension": dimension,
        "factor_order_limit": factor_order_limit,
        "rooted_family_generator_order_limit": factor_order_limit - 3,
        "rooted_family_generator_presentations": generator_presentations,
        "rooted_family_count": len(families),
        "labelled_irredundant_base_candidates": len(
            labelled_base_candidates
        ),
        "canonical_base_candidates": len(ordered_bases),
        "subset_positions_tested": subset_positions_tested,
        "rooted_containments": rooted_containments,
        "availability_disjoint_containments": (
            availability_disjoint_containments
        ),
        "availability_disjoint_pair_type_count": len(disjoint_pairs),
        "availability_disjoint_pair_digest": (
            availability_disjoint_pair_digest
        ),
        "availability_disjoint_pair_types": [
            {
                "index": index,
                "base_mask": format(pair[0], "x"),
                "site_mask": format(pair[1], "x"),
                "base_order": pair[0].bit_count(),
                "site_order": pair[1].bit_count(),
                "factor_order": (
                    pair[0].bit_count() + pair[1].bit_count()
                ),
                "base_available": list(
                    raw_root_avoiding_peripheral_coordinates(
                        pair[0],
                        dimension,
                        0,
                    )
                ),
                "site_available": list(
                    raw_root_avoiding_peripheral_coordinates(
                        pair[1],
                        dimension,
                        0,
                    )
                ),
                "forced_minor_closures": (
                    forced_minor_closure_profile(pair)
                ),
            }
            for index, pair in enumerate(disjoint_pairs)
        ],
        "critical_presentations": critical_presentations,
        "critical_pair_types": len(pairs),
        "pair_digest": pair_digest,
        "order_enriched_repair_signature_types": ordered_repair_types,
        "order_enriched_repair_collision_groups": len(
            ordered_repair_collisions
        ),
        "order_enriched_activation_antimatroid_signature_types": (
            ordered_activation_types
        ),
        "order_enriched_activation_antimatroid_collision_groups": len(
            ordered_activation_collisions
        ),
        "pair_types": [
            critical_pair_payload(pair)
            for pair in pairs
        ],
        "order_enriched_activation_collision_classes": [
            {
                "base_order": members[0][0].bit_count(),
                "site_order": members[0][1].bit_count(),
                "pair_type_indices": [
                    pair_index[pair] for pair in members
                ],
            }
            for members in ordered_activation_collisions
        ],
    }
    Path(output_path).write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print("\nExact rooted critical-pair profile on five old coordinates")
    print(
        f"  factor_order_limit={factor_order_limit}; "
        f"rooted_families={len(families)}; "
        f"canonical_bases={len(ordered_bases)}; "
        f"subset_positions={subset_positions_tested}; "
        f"rooted_containments={rooted_containments}; "
        "availability_disjoint="
        f"{availability_disjoint_containments}; "
        "availability_disjoint_pair_types="
        f"{len(disjoint_pairs)}; "
        f"critical_pair_types={len(pairs)}"
    )
    print(
        "  availability_disjoint_pair_digest="
        f"{availability_disjoint_pair_digest}"
    )
    print(
        "  order_enriched_repair_signatures="
        f"{ordered_repair_types}; "
        "order_enriched_repair_collisions="
        f"{len(ordered_repair_collisions)}"
    )
    print(
        "  order_enriched_activation_antimatroid_signatures="
        f"{ordered_activation_types}; "
        "order_enriched_activation_antimatroid_collisions="
        f"{len(ordered_activation_collisions)}"
    )
    print(f"  pair_digest={pair_digest}")
    print(
        "  peak_rss_mb="
        f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024:.1f}; "
        f"runtime_seconds={time.perf_counter() - started:.3f}"
    )
    print(f"  wrote={output_path}")


def complete_lower_q_activation_obstruction_profile_dimension_four(
    output_path: str,
    site_three_source_path: str | None,
) -> None:
    """Enumerate the marked lower obstructions through four coordinates.

    This is an independent route to the three-vertex-site census.  It scans
    contained pairs ``T <= C`` directly, marks a coordinate ``q`` whose root
    edge lies in ``T``, tests initial availability disjointness, and then
    applies rooted-minor minimality of ``q``-activation.  Adding a fresh
    collar coordinate turns every accepted marked pair into a critical
    three-vertex-site pair.
    """

    started = time.perf_counter()
    maximum_dimension = 4
    profiles: list[dict[str, object]] = []
    target_types: list[tuple[int, int, int]] = []
    target_upper_pairs: list[tuple[int, int]] = []

    def marked_digest(
        pairs: list[tuple[int, int, int]],
    ) -> str:
        return hashlib.sha256(
            (
                "\n".join(
                    f"{base:x}:{site:x}:{q}"
                    for base, site, q in pairs
                )
                + "\n"
            ).encode("ascii")
        ).hexdigest()

    def pair_digest(
        pairs: list[tuple[int, int]],
    ) -> str:
        return hashlib.sha256(
            (
                "\n".join(
                    f"{base:x}:{site:x}"
                    for base, site in pairs
                )
                + "\n"
            ).encode("ascii")
        ).hexdigest()

    for dimension in range(1, maximum_dimension + 1):
        rooted_families = rooted_r_per_families_in_cube(dimension)
        r_per_families = r_per_families_in_cube(dimension)
        rooted_restriction_checks = 0
        for family in rooted_families:
            for coordinate in range(dimension):
                rooted_restriction_checks += 1
                if not is_rooted_r_per(
                    section(family, dimension, coordinate, 0),
                    dimension - 1,
                    0,
                ):
                    raise AssertionError(
                        "rooted-goodness failed root-side restriction "
                        "permanence"
                    )
        irredundant_rooted_sides = tuple(
            family
            for family in rooted_families
            if irredundant_projection(family, dimension)[1] == dimension
        )
        structural_presentations = 0
        initially_locked_presentations = 0
        obstruction_presentations = 0
        marked_types: set[tuple[int, int, int]] = set()

        for rooted_side in irredundant_rooted_sides:
            rooted_available = set(
                raw_root_avoiding_peripheral_coordinates(
                    rooted_side,
                    dimension,
                    0,
                )
            )
            for opposite_side in r_per_families:
                if (
                    opposite_side & ~rooted_side
                    or not (opposite_side & 1)
                ):
                    continue
                opposite_available = set(
                    raw_root_avoiding_peripheral_coordinates(
                        opposite_side,
                        dimension,
                        0,
                    )
                )
                for q in range(dimension):
                    q_neighbor_bit = 1 << (1 << q)
                    if not (opposite_side & q_neighbor_bit):
                        continue
                    structural_presentations += 1
                    if rooted_available & opposite_available:
                        continue
                    initially_locked_presentations += 1
                    if not is_rooted_q_activation_obstruction(
                        rooted_side,
                        opposite_side,
                        dimension,
                        0,
                        q,
                    ):
                        continue
                    obstruction_presentations += 1
                    marked_types.add(
                        canonical_marked_reverse_pair(
                            rooted_side,
                            opposite_side,
                            dimension,
                            q,
                        )
                    )

        ordered_types = sorted(marked_types)
        upper_pairs = sorted(
            {
                canonical_reverse_pair(
                    lift_expansion(
                        rooted_side,
                        opposite_side,
                        dimension,
                    ),
                    family_from_vertices(
                        (
                            0,
                            1 << dimension,
                            (1 << dimension) | (1 << q),
                        )
                    ),
                    dimension + 1,
                )
                for rooted_side, opposite_side, q in ordered_types
            }
        )
        if len(upper_pairs) != len(ordered_types):
            raise AssertionError(
                "distinct marked lower obstructions have isomorphic collars"
            )

        profile = {
            "dimension": dimension,
            "rooted_good_family_count": len(rooted_families),
            "rooted_side_restriction_checks": (
                rooted_restriction_checks
            ),
            "rooted_side_restriction_failures": 0,
            "r_per_family_count": len(r_per_families),
            "irredundant_rooted_side_count": len(
                irredundant_rooted_sides
            ),
            "structural_presentations": structural_presentations,
            "initially_locked_presentations": (
                initially_locked_presentations
            ),
            "obstruction_presentations": obstruction_presentations,
            "canonical_marked_obstruction_types": len(ordered_types),
            "marked_pair_digest": marked_digest(ordered_types),
            "collar_pair_digest": pair_digest(upper_pairs),
        }
        profiles.append(profile)
        print(
            "  lower q-activation progress: "
            f"dimension={dimension}; "
            f"structural={structural_presentations}; "
            f"locked={initially_locked_presentations}; "
            f"obstruction_presentations={obstruction_presentations}; "
            f"types={len(ordered_types)}",
            flush=True,
        )
        if dimension == maximum_dimension:
            target_types = ordered_types
            target_upper_pairs = upper_pairs

    cross_check: dict[str, object] | None = None
    if site_three_source_path is not None:
        source = json.loads(
            Path(site_three_source_path).read_text(encoding="utf-8")
        )
        source_pairs = sorted(
            (
                int(str(row["base_mask"]), 16),
                int(str(row["site_mask"]), 16),
            )
            for row in source["pair_types"]
        )
        if source_pairs != target_upper_pairs:
            raise AssertionError(
                "the direct lower census disagrees with the Q5 collar census"
            )
        expected_presentations = (
            profiles[-1]["obstruction_presentations"]
            * (maximum_dimension + 1)
        )
        if source.get("critical_presentations") != expected_presentations:
            raise AssertionError(
                "the labelled lower and upper presentation counts disagree"
            )
        if source.get("pair_digest") != profiles[-1]["collar_pair_digest"]:
            raise AssertionError(
                "the lower and upper collar digests disagree"
            )
        cross_check = {
            "source_certificate": site_three_source_path,
            "source_pair_digest": source["pair_digest"],
            "source_critical_presentations": source[
                "critical_presentations"
            ],
            "collar_insertions_per_lower_presentation": (
                maximum_dimension + 1
            ),
            "exact_pair_set_match": True,
        }

    def type_payload(
        marked_pair: tuple[int, int, int],
    ) -> dict[str, object]:
        rooted_side, opposite_side, q = marked_pair
        upper_pair = canonical_reverse_pair(
            lift_expansion(
                rooted_side,
                opposite_side,
                maximum_dimension,
            ),
            family_from_vertices(
                (
                    0,
                    1 << maximum_dimension,
                    (1 << maximum_dimension) | (1 << q),
                )
            ),
            maximum_dimension + 1,
        )
        return {
            "rooted_side_mask": format(rooted_side, "x"),
            "opposite_side_mask": format(opposite_side, "x"),
            "activation_coordinate": q,
            "rooted_side_order": rooted_side.bit_count(),
            "opposite_side_order": opposite_side.bit_count(),
            "rooted_side_available": list(
                raw_root_avoiding_peripheral_coordinates(
                    rooted_side,
                    maximum_dimension,
                    0,
                )
            ),
            "opposite_side_available": list(
                raw_root_avoiding_peripheral_coordinates(
                    opposite_side,
                    maximum_dimension,
                    0,
                )
            ),
            "canonical_collar_base_mask": format(upper_pair[0], "x"),
            "canonical_collar_site_mask": format(upper_pair[1], "x"),
        }

    payload = {
        "schema_version": 1,
        "kind": (
            "complete-rooted-q-activation-obstruction-profile-through-d4"
        ),
        "complete_through_dimension": maximum_dimension,
        "profiles": profiles,
        "dimension_four_types": [
            type_payload(marked_pair)
            for marked_pair in target_types
        ],
        "site_three_cross_check": cross_check,
    }
    Path(output_path).write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print("\nComplete lower q-activation obstruction profile")
    print(
        f"  complete_through_dimension={maximum_dimension}; "
        f"dimension_four_types={len(target_types)}; "
        f"marked_pair_digest={profiles[-1]['marked_pair_digest']}; "
        f"collar_pair_digest={profiles[-1]['collar_pair_digest']}"
    )
    print(
        "  peak_rss_mb="
        f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024:.1f}; "
        f"runtime_seconds={time.perf_counter() - started:.3f}"
    )
    print(f"  wrote={output_path}")


def lower_q_fiber_witness_profile_dimension_four(
    input_path: str,
    output_path: str,
) -> None:
    """Profile exact fiber witnesses on the complete marked d=4 census.

    The source certificate supplies all marked rooted q-activation
    obstruction types through coordinate permutations.  Four nested
    signatures record progressively more of the exact q-free activation
    mechanism:

    * ``closure`` keeps only the common q-free feasible-set system;
    * ``logical`` adds the distinct minimal repair clauses of each family;
    * ``multiset`` retains one repair clause for every defective q-fiber;
    * ``incidence`` also records the orders of the two q-sections.

    Every signature includes the initial marked pair and both forced minors
    in each coordinate other than q, then is canonicalized under all marked
    coordinate permutations.
    """

    started = time.perf_counter()
    source = json.loads(Path(input_path).read_text(encoding="utf-8"))
    expected_kind = (
        "complete-rooted-q-activation-obstruction-profile-through-d4"
    )
    if (
        source.get("kind") != expected_kind
        or source.get("complete_through_dimension") != 4
    ):
        raise ValueError(
            "the fiber-witness profile requires the complete marked d=4 "
            "lower-obstruction certificate"
        )
    rows = source.get("dimension_four_types")
    if not isinstance(rows, list) or len(rows) != 174:
        raise ValueError(
            "the complete marked d=4 certificate must contain 174 types"
        )

    dimension = 4
    marked_pairs: list[tuple[int, int, int]] = []
    for row in rows:
        if not isinstance(row, dict):
            raise ValueError("malformed marked-obstruction row")
        pair = (
            int(str(row["rooted_side_mask"]), 16),
            int(str(row["opposite_side_mask"]), 16),
            int(row["activation_coordinate"]),
        )
        if canonical_marked_reverse_pair(*pair[:2], dimension, pair[2]) != pair:
            raise AssertionError(
                "a source marked-obstruction row is not canonical"
            )
        if not is_rooted_q_activation_obstruction(
            pair[0],
            pair[1],
            dimension,
            0,
            pair[2],
        ):
            raise AssertionError(
                "a source row fails the rooted q-activation obstruction test"
            )
        marked_pairs.append(pair)
    if marked_pairs != sorted(set(marked_pairs)):
        raise AssertionError(
            "the source marked-obstruction rows are not sorted and unique"
        )

    fiber_criterion_checks = 0
    for rooted_side, opposite_side, q in marked_pairs:
        states = [
            (rooted_side, opposite_side, dimension, 0, q)
        ]
        for coordinate in range(dimension):
            if coordinate == q:
                continue
            projected_q = q - (q > coordinate)
            states.extend(
                (
                    (
                        contract_coordinate(
                            rooted_side,
                            dimension,
                            coordinate,
                        ),
                        contract_coordinate(
                            opposite_side,
                            dimension,
                            coordinate,
                        ),
                        dimension - 1,
                        0,
                        projected_q,
                    ),
                    (
                        section(
                            rooted_side,
                            dimension,
                            coordinate,
                            0,
                        ),
                        section(
                            opposite_side,
                            dimension,
                            coordinate,
                            0,
                        ),
                        dimension - 1,
                        0,
                        projected_q,
                    ),
                )
            )
        for state_base, state_site, state_dimension, state_root, state_q in (
            states
        ):
            for family in (state_base, state_site):
                repair_supports = minimal_q_repair_supports(
                    family,
                    state_dimension,
                    state_root,
                    state_q,
                )
                for contraction_mask in range(1 << state_dimension):
                    if contraction_mask & (1 << state_q):
                        continue
                    (
                        contracted_family,
                        contracted_dimension,
                        contracted_root,
                    ) = contract_original_coordinate_mask(
                        family,
                        state_dimension,
                        state_root,
                        contraction_mask,
                    )
                    contracted_q = state_q - (
                        contraction_mask & ((1 << state_q) - 1)
                    ).bit_count()
                    direct = contracted_q in (
                        raw_root_avoiding_peripheral_coordinates(
                            contracted_family,
                            contracted_dimension,
                            contracted_root,
                        )
                    )
                    fiber = q_repair_supports_cover(
                        repair_supports,
                        project_vertex(contraction_mask, state_q),
                    )
                    fiber_criterion_checks += 1
                    if direct != fiber:
                        raise AssertionError(
                            "the direct and repair-support fiber criteria "
                            "disagree"
                        )

    def marked_pair_payload(index: int) -> dict[str, object]:
        rooted_side, opposite_side, q = marked_pairs[index]
        return {
            "type_index": index,
            "rooted_side_mask": format(rooted_side, "x"),
            "opposite_side_mask": format(opposite_side, "x"),
            "activation_coordinate": q,
            "rooted_side_order": rooted_side.bit_count(),
            "opposite_side_order": opposite_side.bit_count(),
        }

    mode_payloads: list[dict[str, object]] = []
    modes = ("closure", "logical", "multiset", "incidence")
    for mode in modes:
        classes: dict[
            tuple[object, ...],
            list[int],
        ] = defaultdict(list)
        ordered_signatures: list[tuple[object, ...]] = []
        for index, (rooted_side, opposite_side, q) in enumerate(
            marked_pairs
        ):
            signature = canonical_marked_fiber_witness_signature(
                rooted_side,
                opposite_side,
                dimension,
                q,
                mode,
            )
            ordered_signatures.append(signature)
            classes[signature].append(index)
        collisions = sorted(
            (
                members
                for members in classes.values()
                if len(members) > 1
            ),
            key=lambda members: tuple(members),
        )
        signature_digest = hashlib.sha256(
            (
                json.dumps(
                    ordered_signatures,
                    separators=(",", ":"),
                )
                + "\n"
            ).encode("ascii")
        ).hexdigest()
        mode_payloads.append(
            {
                "mode": mode,
                "signature_types": len(classes),
                "collision_groups": len(collisions),
                "colliding_types": sum(len(group) for group in collisions),
                "largest_collision_group": max(
                    (len(group) for group in collisions),
                    default=1,
                ),
                "signature_digest": signature_digest,
                "collision_classes": [
                    {
                        "type_indices": members,
                        "members": [
                            marked_pair_payload(index)
                            for index in members
                        ],
                    }
                    for members in collisions
                ],
            }
        )
        print(
            "  lower q-fiber witness mode: "
            f"mode={mode}; "
            f"types={len(classes)}; "
            f"collision_groups={len(collisions)}; "
            "colliding_types="
            f"{sum(len(group) for group in collisions)}",
            flush=True,
        )

    source_profile = source["profiles"][-1]
    payload = {
        "schema_version": 1,
        "kind": "rooted-q-fiber-witness-profile-d4",
        "dimension": dimension,
        "source_certificate": input_path,
        "source_marked_pair_digest": source_profile[
            "marked_pair_digest"
        ],
        "source_collar_pair_digest": source_profile[
            "collar_pair_digest"
        ],
        "marked_obstruction_types": len(marked_pairs),
        "fiber_criterion_checks": fiber_criterion_checks,
        "fiber_criterion_failures": 0,
        "signature_modes": mode_payloads,
    }
    Path(output_path).write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print("\nExact lower q-fiber witness profile")
    print(
        f"  marked_obstruction_types={len(marked_pairs)}; "
        f"fiber_criterion_checks={fiber_criterion_checks}; "
        "source_marked_pair_digest="
        f"{source_profile['marked_pair_digest']}"
    )
    print(
        "  peak_rss_mb="
        f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024:.1f}; "
        f"runtime_seconds={time.perf_counter() - started:.3f}"
    )
    print(f"  wrote={output_path}")


def lower_q_signed_repair_functor_profile_dimension_four(
    input_path: str,
    output_path: str,
) -> None:
    """Test the signed repair functor on all marked d=4 obstructions."""

    started = time.perf_counter()
    source = json.loads(Path(input_path).read_text(encoding="utf-8"))
    expected_kind = (
        "complete-rooted-q-activation-obstruction-profile-through-d4"
    )
    if (
        source.get("kind") != expected_kind
        or source.get("complete_through_dimension") != 4
    ):
        raise ValueError(
            "the signed-repair profile requires the complete marked d=4 "
            "lower-obstruction certificate"
        )
    rows = source.get("dimension_four_types")
    if not isinstance(rows, list) or len(rows) != 174:
        raise ValueError(
            "the complete marked d=4 certificate must contain 174 types"
        )

    dimension = 4
    marked_pairs = [
        (
            int(str(row["rooted_side_mask"]), 16),
            int(str(row["opposite_side_mask"]), 16),
            int(row["activation_coordinate"]),
        )
        for row in rows
    ]
    if marked_pairs != sorted(set(marked_pairs)):
        raise AssertionError(
            "the source marked-obstruction rows are not sorted and unique"
        )

    functor_checks = 0
    activation_certificate_checks = 0
    atlas_functor_checks = 0
    atlas_closure_checks = 0
    atlas_horn_family_checks = 0
    atlas_horn_pair_checks = 0
    horn_interval_profiles = {
        "least_marked_valid": 0,
        "greatest_marked_valid": 0,
        "least_single_extension_presentations": 0,
        "least_single_extension_valid": 0,
        "least_single_extension_invalid": 0,
        "original_single_extension_presentations": 0,
        "original_single_extension_valid": 0,
        "original_single_extension_invalid": 0,
        "original_single_deletion_presentations": 0,
        "original_single_deletion_valid": 0,
        "original_single_deletion_invalid": 0,
    }
    horn_interval_failure_counts: Counter[str] = Counter()
    horn_realization_census = {
        "distinct_realizations": 0,
        "valid_realizations": 0,
        "invalid_realizations": 0,
        "intervals_with_upward_violation": 0,
        "principal_valid_intervals": 0,
        "maximum_optional_variables": 0,
        "original_is_unique_minimal_valid": 0,
        "original_is_greatest": 0,
    }
    minimal_valid_realization_histogram: Counter[int] = Counter()
    first_least_extension_failure: dict[str, object] | None = None
    first_original_extension_failure: dict[str, object] | None = None
    first_original_deletion_failure: dict[str, object] | None = None
    first_upward_validity_failure: dict[str, object] | None = None
    first_nonprincipal_valid_interval: dict[str, object] | None = None

    def collar_realization_profile(
        base: int,
        site: int,
        q: int,
    ) -> dict[str, bool]:
        base_partial_cube = is_partial_cube_family(base, dimension)
        site_partial_cube = is_partial_cube_family(site, dimension)
        base_irredundant = (
            irredundant_projection(base, dimension)[1]
            == dimension
        )
        base_rooted_good = (
            base_partial_cube
            and is_rooted_r_per(base, dimension, 0)
        )
        site_r_per = (
            site_partial_cube
            and is_r_per(site, dimension)
        )
        marked_edge = bool(site & (1 << (1 << q)))
        return {
            "base_partial_cube": base_partial_cube,
            "site_partial_cube": site_partial_cube,
            "base_irredundant": base_irredundant,
            "base_rooted_good": base_rooted_good,
            "site_r_per": site_r_per,
            "marked_edge": marked_edge,
            "valid": (
                base_partial_cube
                and site_partial_cube
                and base_irredundant
                and base_rooted_good
                and site_r_per
                and marked_edge
            ),
        }

    for rooted_side, opposite_side, q in marked_pairs:
        initial_systems = {
            family: q_defect_support_system(
                family,
                dimension,
                0,
                q,
            )
            for family in (rooted_side, opposite_side)
        }
        initial_atlases = {
            family: rooted_defect_atlas(
                family,
                dimension,
                0,
            )
            for family in (rooted_side, opposite_side)
        }
        for family in (rooted_side, opposite_side):
            atlas = initial_atlases[family]
            (
                forced_mask,
                forbidden_mask,
                implications,
            ) = defect_atlas_horn_constraints(
                atlas,
                dimension,
                0,
            )
            minimal_family = horn_implication_closure(
                forced_mask,
                implications,
            )
            maximal_family = maximal_safe_horn_assignment(
                1 << dimension,
                forbidden_mask,
                implications,
            )
            checks = (
                not (minimal_family & forbidden_mask),
                family_satisfies_atlas_horn_constraints(
                    family,
                    forced_mask,
                    forbidden_mask,
                    implications,
                ),
                family_satisfies_atlas_horn_constraints(
                    minimal_family,
                    forced_mask,
                    forbidden_mask,
                    implications,
                )
                and rooted_defect_atlas(
                    minimal_family,
                    dimension,
                    0,
                )
                == atlas,
                family_satisfies_atlas_horn_constraints(
                    maximal_family,
                    forced_mask,
                    forbidden_mask,
                    implications,
                )
                and rooted_defect_atlas(
                    maximal_family,
                    dimension,
                    0,
                )
                == atlas,
            )
            atlas_horn_family_checks += len(checks)
            if not all(checks):
                raise AssertionError(
                    "the Horn atlas-realization theorem failed on a "
                    "known family"
                )

        (
            pair_forced,
            pair_forbidden,
            pair_implications,
        ) = contained_pair_atlas_horn_constraints(
            initial_atlases[rooted_side],
            initial_atlases[opposite_side],
            dimension,
            0,
        )
        pair_minimal = horn_implication_closure(
            pair_forced,
            pair_implications,
        )
        vertex_count = 1 << dimension
        vertex_variable_mask = (1 << vertex_count) - 1
        minimal_base = pair_minimal & vertex_variable_mask
        minimal_site = (
            pair_minimal >> vertex_count
        ) & vertex_variable_mask
        original_pair_mask = (
            rooted_side
            | (opposite_side << vertex_count)
        )
        pair_checks = (
            not (pair_minimal & pair_forbidden),
            family_satisfies_atlas_horn_constraints(
                original_pair_mask,
                pair_forced,
                pair_forbidden,
                pair_implications,
            ),
            not (minimal_site & ~minimal_base),
            rooted_defect_atlas(
                minimal_base,
                dimension,
                0,
            )
            == initial_atlases[rooted_side],
            rooted_defect_atlas(
                minimal_site,
                dimension,
                0,
            )
            == initial_atlases[opposite_side],
        )
        atlas_horn_pair_checks += len(pair_checks)
        if not all(pair_checks):
            raise AssertionError(
                "the Horn contained-pair realization theorem failed"
            )

        marked_site_neighbor_variable = (
            vertex_count + (1 << q)
        )
        marked_pair_minimal = horn_implication_closure(
            pair_forced | (1 << marked_site_neighbor_variable),
            pair_implications,
        )
        if marked_pair_minimal & pair_forbidden:
            raise AssertionError(
                "the known marked pair has no Horn realization retaining "
                "its marked edge"
            )
        marked_pair_greatest = maximal_safe_horn_assignment(
            2 * vertex_count,
            pair_forbidden,
            pair_implications,
        )

        def decode_pair(
            assignment: int,
        ) -> tuple[int, int]:
            return (
                assignment & vertex_variable_mask,
                (assignment >> vertex_count) & vertex_variable_mask,
            )

        least_base, least_site = decode_pair(marked_pair_minimal)
        greatest_base, greatest_site = decode_pair(
            marked_pair_greatest
        )
        least_profile = collar_realization_profile(
            least_base,
            least_site,
            q,
        )
        greatest_profile = collar_realization_profile(
            greatest_base,
            greatest_site,
            q,
        )
        horn_interval_profiles["least_marked_valid"] += int(
            least_profile["valid"]
        )
        horn_interval_profiles["greatest_marked_valid"] += int(
            greatest_profile["valid"]
        )

        optional_variables = (
            marked_pair_greatest & ~marked_pair_minimal
        )
        for variable in range(2 * vertex_count):
            if not (optional_variables & (1 << variable)):
                continue
            extension = horn_implication_closure(
                marked_pair_minimal | (1 << variable),
                pair_implications,
            )
            if extension & pair_forbidden:
                raise AssertionError(
                    "a greatest-safe Horn variable produced a forbidden "
                    "extension"
                )
            extension_base, extension_site = decode_pair(extension)
            extension_profile = collar_realization_profile(
                extension_base,
                extension_site,
                q,
            )
            horn_interval_profiles[
                "least_single_extension_presentations"
            ] += 1
            if extension_profile["valid"]:
                horn_interval_profiles[
                    "least_single_extension_valid"
                ] += 1
                continue
            horn_interval_profiles[
                "least_single_extension_invalid"
            ] += 1
            for key, value in extension_profile.items():
                if key != "valid" and not value:
                    horn_interval_failure_counts[key] += 1
            if first_least_extension_failure is None:
                role = (
                    "base"
                    if variable < vertex_count
                    else "site"
                )
                vertex = (
                    variable
                    if variable < vertex_count
                    else variable - vertex_count
                )
                first_least_extension_failure = {
                    "source_rooted_side_mask": format(
                        rooted_side,
                        "x",
                    ),
                    "source_opposite_side_mask": format(
                        opposite_side,
                        "x",
                    ),
                    "activation_coordinate": q,
                    "added_role": role,
                    "added_vertex": format(
                        vertex,
                        f"0{dimension}b",
                    ),
                    "least_base_mask": format(least_base, "x"),
                    "least_site_mask": format(least_site, "x"),
                    "extension_base_mask": format(
                        extension_base,
                        "x",
                    ),
                    "extension_site_mask": format(
                        extension_site,
                        "x",
                    ),
                    "failed_properties": [
                        key
                        for key, value in extension_profile.items()
                        if key != "valid" and not value
                    ],
                }

        original_optional_variables = (
            marked_pair_greatest & ~original_pair_mask
        )
        for variable in range(2 * vertex_count):
            if not (
                original_optional_variables & (1 << variable)
            ):
                continue
            extension = horn_implication_closure(
                original_pair_mask | (1 << variable),
                pair_implications,
            )
            if extension & pair_forbidden:
                raise AssertionError(
                    "a greatest-safe Horn variable produced a forbidden "
                    "extension of the original pair"
                )
            extension_base, extension_site = decode_pair(extension)
            extension_profile = collar_realization_profile(
                extension_base,
                extension_site,
                q,
            )
            horn_interval_profiles[
                "original_single_extension_presentations"
            ] += 1
            if extension_profile["valid"]:
                horn_interval_profiles[
                    "original_single_extension_valid"
                ] += 1
                continue
            horn_interval_profiles[
                "original_single_extension_invalid"
            ] += 1
            if first_original_extension_failure is None:
                role = (
                    "base"
                    if variable < vertex_count
                    else "site"
                )
                vertex = (
                    variable
                    if variable < vertex_count
                    else variable - vertex_count
                )
                first_original_extension_failure = {
                    "source_rooted_side_mask": format(
                        rooted_side,
                        "x",
                    ),
                    "source_opposite_side_mask": format(
                        opposite_side,
                        "x",
                    ),
                    "activation_coordinate": q,
                    "added_role": role,
                    "added_vertex": format(
                        vertex,
                        f"0{dimension}b",
                    ),
                    "extension_base_mask": format(
                        extension_base,
                        "x",
                    ),
                    "extension_site_mask": format(
                        extension_site,
                        "x",
                    ),
                    "failed_properties": [
                        key
                        for key, value in extension_profile.items()
                        if key != "valid" and not value
                    ],
                }

        removable_variables = (
            original_pair_mask
            & ~pair_forced
            & ~(1 << marked_site_neighbor_variable)
        )
        for variable in range(2 * vertex_count):
            if not (removable_variables & (1 << variable)):
                continue
            deletion = original_pair_mask & ~(1 << variable)
            if not family_satisfies_atlas_horn_constraints(
                deletion,
                pair_forced,
                pair_forbidden,
                pair_implications,
            ):
                continue
            deletion_base, deletion_site = decode_pair(deletion)
            deletion_profile = collar_realization_profile(
                deletion_base,
                deletion_site,
                q,
            )
            horn_interval_profiles[
                "original_single_deletion_presentations"
            ] += 1
            if deletion_profile["valid"]:
                horn_interval_profiles[
                    "original_single_deletion_valid"
                ] += 1
                continue
            horn_interval_profiles[
                "original_single_deletion_invalid"
            ] += 1
            if first_original_deletion_failure is None:
                role = (
                    "base"
                    if variable < vertex_count
                    else "site"
                )
                vertex = (
                    variable
                    if variable < vertex_count
                    else variable - vertex_count
                )
                first_original_deletion_failure = {
                    "source_rooted_side_mask": format(
                        rooted_side,
                        "x",
                    ),
                    "source_opposite_side_mask": format(
                        opposite_side,
                        "x",
                    ),
                    "activation_coordinate": q,
                    "deleted_role": role,
                    "deleted_vertex": format(
                        vertex,
                        f"0{dimension}b",
                    ),
                    "deletion_base_mask": format(
                        deletion_base,
                        "x",
                    ),
                    "deletion_site_mask": format(
                        deletion_site,
                        "x",
                    ),
                    "failed_properties": [
                        key
                        for key, value in deletion_profile.items()
                        if key != "valid" and not value
                    ],
                }

        optional_indices = [
            variable
            for variable in range(2 * vertex_count)
            if optional_variables & (1 << variable)
        ]
        horn_realization_census[
            "maximum_optional_variables"
        ] = max(
            horn_realization_census["maximum_optional_variables"],
            len(optional_indices),
        )
        realizations: set[int] = set()
        for subset_mask in range(1 << len(optional_indices)):
            seed = marked_pair_minimal
            for optional_index, variable in enumerate(
                optional_indices
            ):
                if subset_mask & (1 << optional_index):
                    seed |= 1 << variable
            realization = horn_implication_closure(
                seed,
                pair_implications,
            )
            if realization & pair_forbidden:
                raise AssertionError(
                    "an optional subset reached a forbidden Horn variable"
                )
            realizations.add(realization)

        valid_assignments: list[int] = []
        invalid_assignments: list[int] = []
        for realization in sorted(realizations):
            realization_base, realization_site = decode_pair(
                realization
            )
            profile = collar_realization_profile(
                realization_base,
                realization_site,
                q,
            )
            if profile["valid"]:
                valid_assignments.append(realization)
            else:
                invalid_assignments.append(realization)
        horn_realization_census["distinct_realizations"] += len(
            realizations
        )
        horn_realization_census["valid_realizations"] += len(
            valid_assignments
        )
        horn_realization_census["invalid_realizations"] += len(
            invalid_assignments
        )
        minimal_valid_assignments = [
            realization
            for realization in valid_assignments
            if not any(
                other != realization
                and not (other & ~realization)
                for other in valid_assignments
            )
        ]
        minimal_valid_realization_histogram[
            len(minimal_valid_assignments)
        ] += 1
        if (
            len(minimal_valid_assignments) == 1
            and original_pair_mask == minimal_valid_assignments[0]
        ):
            horn_realization_census[
                "original_is_unique_minimal_valid"
            ] += 1
        if original_pair_mask == marked_pair_greatest:
            horn_realization_census[
                "original_is_greatest"
            ] += 1

        upward_failure: tuple[int, int] | None = None
        for valid_assignment in valid_assignments:
            for invalid_assignment in invalid_assignments:
                if not (valid_assignment & ~invalid_assignment):
                    upward_failure = (
                        valid_assignment,
                        invalid_assignment,
                    )
                    break
            if upward_failure is not None:
                break
        if upward_failure is not None:
            horn_realization_census[
                "intervals_with_upward_violation"
            ] += 1
            if first_upward_validity_failure is None:
                valid_base, valid_site = decode_pair(
                    upward_failure[0]
                )
                invalid_base, invalid_site = decode_pair(
                    upward_failure[1]
                )
                first_upward_validity_failure = {
                    "source_rooted_side_mask": format(
                        rooted_side,
                        "x",
                    ),
                    "source_opposite_side_mask": format(
                        opposite_side,
                        "x",
                    ),
                    "activation_coordinate": q,
                    "valid_base_mask": format(valid_base, "x"),
                    "valid_site_mask": format(valid_site, "x"),
                    "larger_invalid_base_mask": format(
                        invalid_base,
                        "x",
                    ),
                    "larger_invalid_site_mask": format(
                        invalid_site,
                        "x",
                    ),
                }

        principal = False
        if valid_assignments:
            valid_core = valid_assignments[0]
            for realization in valid_assignments[1:]:
                valid_core &= realization
            principal = (
                valid_core in realizations
                and valid_core in valid_assignments
                and all(
                    bool(valid_core & ~realization)
                    or realization in valid_assignments
                    for realization in realizations
                )
            )
        if principal:
            horn_realization_census[
                "principal_valid_intervals"
            ] += 1
        elif first_nonprincipal_valid_interval is None:
            first_nonprincipal_valid_interval = {
                "source_rooted_side_mask": format(
                    rooted_side,
                    "x",
                ),
                "source_opposite_side_mask": format(
                    opposite_side,
                    "x",
                ),
                "activation_coordinate": q,
                "realizations": len(realizations),
                "valid_realizations": len(valid_assignments),
                "minimal_valid_realizations": len(
                    minimal_valid_assignments
                ),
            }
        state_rows: list[
            tuple[int, int, int, int, bool]
        ] = [
            (
                rooted_side,
                opposite_side,
                dimension,
                q,
                False,
            )
        ]
        for coordinate in range(dimension):
            if coordinate == q:
                continue
            projected_q = q - (q > coordinate)
            fiber_coordinate = coordinate - (coordinate > q)
            for family in (rooted_side, opposite_side):
                initial_system = initial_systems[family]
                initial_atlas = initial_atlases[family]
                contracted_family = contract_coordinate(
                    family,
                    dimension,
                    coordinate,
                )
                contracted_direct = q_defect_support_system(
                    contracted_family,
                    dimension - 1,
                    0,
                    projected_q,
                )
                contracted_functor = contract_q_defect_support_system(
                    initial_system,
                    dimension - 1,
                    fiber_coordinate,
                )
                functor_checks += 1
                if contracted_direct != contracted_functor:
                    raise AssertionError(
                        "the contraction repair functor disagrees with "
                        "direct recomputation"
                    )

                contracted_atlas_direct = rooted_defect_atlas(
                    contracted_family,
                    dimension - 1,
                    0,
                )
                contracted_atlas_functor = contract_rooted_defect_atlas(
                    initial_atlas,
                    dimension,
                    coordinate,
                )
                atlas_functor_checks += 1
                if contracted_atlas_direct != contracted_atlas_functor:
                    raise AssertionError(
                        "the contraction atlas functor disagrees with "
                        "direct recomputation"
                    )

                restricted_family = section(
                    family,
                    dimension,
                    coordinate,
                    0,
                )
                restricted_direct = q_defect_support_system(
                    restricted_family,
                    dimension - 1,
                    0,
                    projected_q,
                )
                restricted_functor = (
                    restrict_q_defect_support_system_to_root(
                        initial_system,
                        dimension - 1,
                        fiber_coordinate,
                        0,
                    )
                )
                functor_checks += 1
                if restricted_direct != restricted_functor:
                    raise AssertionError(
                        "the restriction repair functor disagrees with "
                        "direct recomputation"
                    )

                restricted_atlas_direct = rooted_defect_atlas(
                    restricted_family,
                    dimension - 1,
                    0,
                )
                restricted_atlas_functor = (
                    restrict_rooted_defect_atlas_to_root(
                        initial_atlas,
                        dimension,
                        0,
                        coordinate,
                    )
                )
                atlas_functor_checks += 1
                if restricted_atlas_direct != restricted_atlas_functor:
                    raise AssertionError(
                        "the restriction atlas functor disagrees with "
                        "direct recomputation"
                    )

            state_rows.extend(
                (
                    (
                        contract_coordinate(
                            rooted_side,
                            dimension,
                            coordinate,
                        ),
                        contract_coordinate(
                            opposite_side,
                            dimension,
                            coordinate,
                        ),
                        dimension - 1,
                        projected_q,
                        True,
                    ),
                    (
                        section(
                            rooted_side,
                            dimension,
                            coordinate,
                            0,
                        ),
                        section(
                            opposite_side,
                            dimension,
                            coordinate,
                            0,
                        ),
                        dimension - 1,
                        projected_q,
                        True,
                    ),
                )
            )

        for (
            state_base,
            state_site,
            state_dimension,
            state_q,
            expected_activation,
        ) in state_rows:
            contracted_mask = common_availability_closure_avoiding(
                state_base,
                state_site,
                state_dimension,
                0,
                1 << state_q,
            )[0]
            atlas_contracted_mask = (
                common_availability_closure_from_defect_atlases(
                    rooted_defect_atlas(
                        state_base,
                        state_dimension,
                        0,
                    ),
                    rooted_defect_atlas(
                        state_site,
                        state_dimension,
                        0,
                    ),
                    state_dimension,
                    1 << state_q,
                )
            )
            atlas_closure_checks += 1
            if contracted_mask != atlas_contracted_mask:
                raise AssertionError(
                    "the defect-atlas closure disagrees with the family "
                    "closure"
                )
            fiber_mask = project_vertex(contracted_mask, state_q)
            systems_cover = all(
                q_repair_supports_cover(
                    tuple(
                        clause
                        for _, clause in q_defect_support_system(
                            family,
                            state_dimension,
                            0,
                            state_q,
                        )
                    ),
                    fiber_mask,
                )
                for family in (state_base, state_site)
            )
            direct_activation = is_q_activated(
                state_base,
                state_site,
                state_dimension,
                0,
                state_q,
            )
            activation_certificate_checks += 1
            if (
                systems_cover != direct_activation
                or direct_activation != expected_activation
            ):
                raise AssertionError(
                    "the signed repair minimality certificate disagrees "
                    "with direct marked activation"
                )

    def marked_pair_payload(index: int) -> dict[str, object]:
        rooted_side, opposite_side, q = marked_pairs[index]
        return {
            "type_index": index,
            "rooted_side_mask": format(rooted_side, "x"),
            "opposite_side_mask": format(opposite_side, "x"),
            "activation_coordinate": q,
            "rooted_side_order": rooted_side.bit_count(),
            "opposite_side_order": opposite_side.bit_count(),
        }

    mode_payloads: list[dict[str, object]] = []
    modes = (
        "signed",
        "signed_orders",
        "signed_section_profile",
        "certificate",
        "certificate_initial_orders",
        "certificate_orders",
    )
    for mode in modes:
        classes: dict[
            tuple[object, ...],
            list[int],
        ] = defaultdict(list)
        ordered_signatures: list[tuple[object, ...]] = []
        for index, (rooted_side, opposite_side, q) in enumerate(
            marked_pairs
        ):
            signature = canonical_marked_signed_repair_signature(
                rooted_side,
                opposite_side,
                dimension,
                q,
                mode,
            )
            ordered_signatures.append(signature)
            classes[signature].append(index)
        collisions = sorted(
            (
                members
                for members in classes.values()
                if len(members) > 1
            ),
            key=lambda members: tuple(members),
        )
        mode_payloads.append(
            {
                "mode": mode,
                "signature_types": len(classes),
                "collision_groups": len(collisions),
                "colliding_types": sum(len(group) for group in collisions),
                "largest_collision_group": max(
                    (len(group) for group in collisions),
                    default=1,
                ),
                "signature_digest": hashlib.sha256(
                    (
                        json.dumps(
                            ordered_signatures,
                            separators=(",", ":"),
                        )
                        + "\n"
                    ).encode("ascii")
                ).hexdigest(),
                "collision_classes": [
                    {
                        "type_indices": members,
                        "members": [
                            marked_pair_payload(index)
                            for index in members
                        ],
                    }
                    for members in collisions
                ],
            }
        )
        print(
            "  signed repair mode: "
            f"mode={mode}; "
            f"types={len(classes)}; "
            f"collision_groups={len(collisions)}; "
            "colliding_types="
            f"{sum(len(group) for group in collisions)}",
            flush=True,
        )

    atlas_mode_payloads: list[dict[str, object]] = []
    for mode in ("atlas", "atlas_orders", "atlas_section_orders"):
        classes: dict[
            tuple[object, ...],
            list[int],
        ] = defaultdict(list)
        ordered_signatures: list[tuple[object, ...]] = []
        for index, (rooted_side, opposite_side, q) in enumerate(
            marked_pairs
        ):
            signature = canonical_marked_defect_atlas_signature(
                rooted_side,
                opposite_side,
                dimension,
                q,
                mode,
            )
            ordered_signatures.append(signature)
            classes[signature].append(index)
        collisions = sorted(
            (
                members
                for members in classes.values()
                if len(members) > 1
            ),
            key=lambda members: tuple(members),
        )
        atlas_mode_payloads.append(
            {
                "mode": mode,
                "signature_types": len(classes),
                "collision_groups": len(collisions),
                "colliding_types": sum(len(group) for group in collisions),
                "largest_collision_group": max(
                    (len(group) for group in collisions),
                    default=1,
                ),
                "signature_digest": hashlib.sha256(
                    (
                        json.dumps(
                            ordered_signatures,
                            separators=(",", ":"),
                        )
                        + "\n"
                    ).encode("ascii")
                ).hexdigest(),
                "collision_classes": [
                    {
                        "type_indices": members,
                        "members": [
                            marked_pair_payload(index)
                            for index in members
                        ],
                    }
                    for members in collisions
                ],
            }
        )
        print(
            "  defect atlas mode: "
            f"mode={mode}; "
            f"types={len(classes)}; "
            f"collision_groups={len(collisions)}; "
            "colliding_types="
            f"{sum(len(group) for group in collisions)}",
            flush=True,
        )

    source_profile = source["profiles"][-1]
    payload = {
        "schema_version": 1,
        "kind": "rooted-q-signed-repair-functor-profile-d4",
        "dimension": dimension,
        "source_certificate": input_path,
        "source_marked_pair_digest": source_profile[
            "marked_pair_digest"
        ],
        "source_collar_pair_digest": source_profile[
            "collar_pair_digest"
        ],
        "marked_obstruction_types": len(marked_pairs),
        "functor_checks": functor_checks,
        "functor_failures": 0,
        "atlas_functor_checks": atlas_functor_checks,
        "atlas_functor_failures": 0,
        "atlas_closure_checks": atlas_closure_checks,
        "atlas_closure_failures": 0,
        "atlas_horn_family_checks": atlas_horn_family_checks,
        "atlas_horn_family_failures": 0,
        "atlas_horn_pair_checks": atlas_horn_pair_checks,
        "atlas_horn_pair_failures": 0,
        "horn_interval_profiles": horn_interval_profiles,
        "horn_interval_failure_counts": dict(
            sorted(horn_interval_failure_counts.items())
        ),
        "first_least_extension_failure": first_least_extension_failure,
        "first_original_extension_failure": (
            first_original_extension_failure
        ),
        "first_original_deletion_failure": (
            first_original_deletion_failure
        ),
        "horn_realization_census": horn_realization_census,
        "minimal_valid_realization_histogram": {
            str(key): value
            for key, value in sorted(
                minimal_valid_realization_histogram.items()
            )
        },
        "first_upward_validity_failure": (
            first_upward_validity_failure
        ),
        "first_nonprincipal_valid_interval": (
            first_nonprincipal_valid_interval
        ),
        "activation_certificate_checks": activation_certificate_checks,
        "activation_certificate_failures": 0,
        "signature_modes": mode_payloads,
        "atlas_signature_modes": atlas_mode_payloads,
    }
    Path(output_path).write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print("\nExact lower signed-repair functor profile")
    print(
        f"  marked_obstruction_types={len(marked_pairs)}; "
        f"functor_checks={functor_checks}; "
        f"atlas_functor_checks={atlas_functor_checks}; "
        f"atlas_closure_checks={atlas_closure_checks}; "
        f"atlas_horn_family_checks={atlas_horn_family_checks}; "
        f"atlas_horn_pair_checks={atlas_horn_pair_checks}; "
        "activation_certificate_checks="
        f"{activation_certificate_checks}"
    )
    print(
        "  Horn interval: "
        + "; ".join(
            f"{key}={value}"
            for key, value in horn_interval_profiles.items()
        )
    )
    print(
        "  Horn interval failures: "
        + ", ".join(
            f"{key}={value}"
            for key, value in sorted(
                horn_interval_failure_counts.items()
            )
        )
    )
    print(
        "  Horn realization census: "
        + "; ".join(
            f"{key}={value}"
            for key, value in horn_realization_census.items()
        )
    )
    print(
        "  minimal valid realization histogram: "
        + ", ".join(
            f"{key}:{value}"
            for key, value in sorted(
                minimal_valid_realization_histogram.items()
            )
        )
    )
    print(
        "  peak_rss_mb="
        f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024:.1f}; "
        f"runtime_seconds={time.perf_counter() - started:.3f}"
    )
    print(f"  wrote={output_path}")


def principal_filter_falsifier_profile(
    output_path: str,
) -> None:
    """Refute atlas-fiber upward closure and audit its surviving proxies.

    The complete small-family part groups every family containing root zero
    by its rooted defect atlas.  The dimension-five part gives a canonical
    marked contained-pair interval in which adding one site vertex restores
    the pentagonal maximum-class obstruction and destroys ``R_per`` without
    changing either atlas.
    """

    started = time.perf_counter()

    def property_profile(
        dimension: int,
        groups: dict[
            tuple[
                tuple[tuple[int, tuple[int, ...]], ...],
                ...,
            ],
            list[int],
        ],
        label: str,
    ) -> dict[str, object]:
        def has_property(family: int) -> bool:
            partial_cube = is_partial_cube_family(
                family,
                dimension,
            )
            if label == "partial_cube":
                return partial_cube
            if label == "partial_cube_irredundant":
                return (
                    partial_cube
                    and irredundant_projection(
                        family,
                        dimension,
                    )[1]
                    == dimension
                )
            if label == "r_per":
                return (
                    partial_cube
                    and is_r_per(family, dimension)
                )
            if label == "rooted_r_per":
                return (
                    partial_cube
                    and is_rooted_r_per(
                        family,
                        dimension,
                        0,
                    )
                )
            raise ValueError(f"unknown atlas-fiber property: {label}")

        fibers_with_valid = 0
        mixed_fibers = 0
        upward_violation_fibers = 0
        nonprincipal_fibers = 0
        first_mixed: dict[str, str] | None = None
        first_upward_violation: dict[str, str] | None = None
        first_nonprincipal: dict[str, object] | None = None
        for families in groups.values():
            valid = [
                family
                for family in families
                if has_property(family)
            ]
            invalid = [
                family
                for family in families
                if not has_property(family)
            ]
            if not valid:
                continue
            fibers_with_valid += 1
            if invalid:
                mixed_fibers += 1
                candidate = {
                    "invalid_mask": format(
                        min(
                            invalid,
                            key=lambda family: (
                                family.bit_count(),
                                family,
                            ),
                        ),
                        "x",
                    ),
                    "valid_mask": format(
                        min(
                            valid,
                            key=lambda family: (
                                family.bit_count(),
                                family,
                            ),
                        ),
                        "x",
                    ),
                }
                if first_mixed is None:
                    first_mixed = candidate

            upward_failure = next(
                (
                    (small, large)
                    for small in sorted(
                        valid,
                        key=lambda family: (
                            family.bit_count(),
                            family,
                        ),
                    )
                    for large in sorted(
                        invalid,
                        key=lambda family: (
                            family.bit_count(),
                            family,
                        ),
                    )
                    if not (small & ~large)
                ),
                None,
            )
            if upward_failure is not None:
                upward_violation_fibers += 1
                if first_upward_violation is None:
                    first_upward_violation = {
                        "valid_mask": format(
                            upward_failure[0],
                            "x",
                        ),
                        "larger_invalid_mask": format(
                            upward_failure[1],
                            "x",
                        ),
                    }

            minimal_valid = [
                family
                for family in valid
                if not any(
                    other != family
                    and not (other & ~family)
                    for other in valid
                )
            ]
            principal = (
                len(minimal_valid) == 1
                and upward_failure is None
            )
            if not principal:
                nonprincipal_fibers += 1
                if first_nonprincipal is None:
                    first_nonprincipal = {
                        "minimal_valid_masks": [
                            format(family, "x")
                            for family in minimal_valid
                        ],
                        "upward_violation": (
                            None
                            if upward_failure is None
                            else {
                                "valid_mask": format(
                                    upward_failure[0],
                                    "x",
                                ),
                                "larger_invalid_mask": format(
                                    upward_failure[1],
                                    "x",
                                ),
                            }
                        ),
                    }
        return {
            "fibers_with_valid_realization": fibers_with_valid,
            "mixed_fibers": mixed_fibers,
            "upward_violation_fibers": upward_violation_fibers,
            "nonprincipal_valid_fibers": nonprincipal_fibers,
            "first_mixed_fiber": first_mixed,
            "first_upward_violation": first_upward_violation,
            "first_nonprincipal_fiber": first_nonprincipal,
        }

    dimension_profiles: list[dict[str, object]] = []
    for dimension in range(5):
        cube_order = 1 << dimension
        groups: dict[
            tuple[
                tuple[tuple[int, tuple[int, ...]], ...],
                ...,
            ],
            list[int],
        ] = defaultdict(list)
        for optional in range(1 << (cube_order - 1)):
            family = 1 | (optional << 1)
            groups[
                rooted_defect_atlas(
                    family,
                    dimension,
                    0,
                )
            ].append(family)
        properties = {
            label: property_profile(
                dimension,
                groups,
                label,
            )
            for label in (
                "partial_cube",
                "partial_cube_irredundant",
                "r_per",
                "rooted_r_per",
            )
        }
        if properties["partial_cube"][
            "upward_violation_fibers"
        ]:
            raise AssertionError(
                "a fixed-atlas extension destroyed partial-cube membership"
            )
        dimension_profiles.append(
            {
                "dimension": dimension,
                "rooted_families": 1 << (cube_order - 1),
                "atlas_fibers": len(groups),
                "largest_atlas_fiber": max(
                    len(families)
                    for families in groups.values()
                ),
                "properties": properties,
            }
        )
        print(
            "  complete rooted-family atlas audit: "
            f"dimension={dimension}; "
            f"families={1 << (cube_order - 1)}; "
            f"fibers={len(groups)}; "
            "partial_cube_upward_failures="
            f"{properties['partial_cube']['upward_violation_fibers']}; "
            "r_per_upward_failures="
            f"{properties['r_per']['upward_violation_fibers']}",
            flush=True,
        )

    dimension = 5
    root = 0
    q = 3
    base = int("18fb3af", 16)
    site = int("8fb3ab", 16)
    invalid_site = int("18fb3ab", 16)
    if canonical_marked_pair_extension(
        base,
        site,
        invalid_site,
        dimension,
        q,
    ) != (base, site, invalid_site, q):
        raise AssertionError(
            "the displayed principal-filter falsifier is not canonical"
        )
    if not (
        canonical_cube_family(
            invalid_site,
            dimension,
        )
        == canonical_cube_family(
            pentagonal_maximum_class(),
            dimension,
        )
    ):
        raise AssertionError(
            "the larger site is not the pentagonal obstruction"
        )

    base_atlas = rooted_defect_atlas(
        base,
        dimension,
        root,
    )
    site_atlas = rooted_defect_atlas(
        site,
        dimension,
        root,
    )
    invalid_site_atlas = rooted_defect_atlas(
        invalid_site,
        dimension,
        root,
    )
    (
        forced_mask,
        forbidden_mask,
        implications,
    ) = contained_pair_atlas_horn_constraints(
        base_atlas,
        site_atlas,
        dimension,
        root,
    )
    vertex_count = 1 << dimension
    vertex_variable_mask = (1 << vertex_count) - 1
    marked_site_neighbor = vertex_count + (1 << q)
    least = horn_implication_closure(
        forced_mask | (1 << marked_site_neighbor),
        implications,
    )
    greatest = maximal_safe_horn_assignment(
        2 * vertex_count,
        forbidden_mask,
        implications,
    )
    optional_variables = greatest & ~least
    optional_indices = [
        variable
        for variable in range(2 * vertex_count)
        if optional_variables & (1 << variable)
    ]

    def decode_pair(assignment: int) -> tuple[int, int]:
        return (
            assignment & vertex_variable_mask,
            (assignment >> vertex_count) & vertex_variable_mask,
        )

    def collar_profile(
        current_base: int,
        current_site: int,
    ) -> dict[str, bool]:
        base_partial_cube = is_partial_cube_family(
            current_base,
            dimension,
        )
        site_partial_cube = is_partial_cube_family(
            current_site,
            dimension,
        )
        base_irredundant = (
            irredundant_projection(
                current_base,
                dimension,
            )[1]
            == dimension
        )
        base_rooted_r_per = (
            base_partial_cube
            and is_rooted_r_per(
                current_base,
                dimension,
                root,
            )
        )
        site_r_per = (
            site_partial_cube
            and is_r_per(current_site, dimension)
        )
        marked_edge = bool(
            current_site
            & (1 << (root ^ (1 << q)))
        )
        return {
            "base_partial_cube": base_partial_cube,
            "site_partial_cube": site_partial_cube,
            "base_irredundant": base_irredundant,
            "base_rooted_r_per": base_rooted_r_per,
            "site_r_per": site_r_per,
            "marked_edge": marked_edge,
            "valid": (
                base_partial_cube
                and site_partial_cube
                and base_irredundant
                and base_rooted_r_per
                and site_r_per
                and marked_edge
            ),
        }

    realizations: set[int] = set()
    for subset_mask in range(1 << len(optional_indices)):
        seed = least
        for index, variable in enumerate(optional_indices):
            if subset_mask & (1 << index):
                seed |= 1 << variable
        realization = horn_implication_closure(
            seed,
            implications,
        )
        if realization & forbidden_mask:
            raise AssertionError(
                "a P5 interval realization reaches a forbidden variable"
            )
        realizations.add(realization)

    interval_rows: list[dict[str, object]] = []
    valid_assignments: list[int] = []
    invalid_assignments: list[int] = []
    for realization in sorted(
        realizations,
        key=lambda assignment: (
            decode_pair(assignment)[0].bit_count()
            + decode_pair(assignment)[1].bit_count(),
            decode_pair(assignment),
        ),
    ):
        current_base, current_site = decode_pair(realization)
        profile = collar_profile(
            current_base,
            current_site,
        )
        (valid_assignments if profile["valid"] else invalid_assignments).append(
            realization
        )
        interval_rows.append(
            {
                "base_mask": format(current_base, "x"),
                "site_mask": format(current_site, "x"),
                "base_order": current_base.bit_count(),
                "site_order": current_site.bit_count(),
                "profile": profile,
            }
        )
    minimal_valid_assignments = [
        realization
        for realization in valid_assignments
        if not any(
            other != realization
            and not (other & ~realization)
            for other in valid_assignments
        )
    ]
    upward_failures = [
        (small, large)
        for small in valid_assignments
        for large in invalid_assignments
        if not (small & ~large)
    ]

    valid_assignment = base | (site << vertex_count)
    invalid_assignment = base | (invalid_site << vertex_count)
    exact_checks = {
        "site_contained_in_invalid_site": not (
            site & ~invalid_site
        ),
        "invalid_site_contained_in_base": not (
            invalid_site & ~base
        ),
        "base_partial_cube": is_partial_cube_family(
            base,
            dimension,
        ),
        "base_irredundant": (
            irredundant_projection(
                base,
                dimension,
            )[1]
            == dimension
        ),
        "base_rooted_r_per": is_rooted_r_per(
            base,
            dimension,
            root,
        ),
        "site_partial_cube": is_partial_cube_family(
            site,
            dimension,
        ),
        "site_r_per": is_r_per(site, dimension),
        "invalid_site_partial_cube": is_partial_cube_family(
            invalid_site,
            dimension,
        ),
        "invalid_site_periphery_free": not peripheral_decompositions(
            invalid_site,
            dimension,
        ),
        "invalid_site_not_r_per": not is_r_per(
            invalid_site,
            dimension,
        ),
        "same_site_atlas": site_atlas == invalid_site_atlas,
        "valid_pair_in_horn_interval": (
            family_satisfies_atlas_horn_constraints(
                valid_assignment,
                forced_mask,
                forbidden_mask,
                implications,
            )
        ),
        "invalid_pair_in_horn_interval": (
            family_satisfies_atlas_horn_constraints(
                invalid_assignment,
                forced_mask,
                forbidden_mask,
                implications,
            )
        ),
        "valid_pair_below_invalid_pair": not (
            valid_assignment & ~invalid_assignment
        ),
        "marked_activation_obstruction": (
            is_rooted_q_activation_obstruction(
                base,
                site,
                dimension,
                root,
                q,
            )
        ),
    }
    if not all(exact_checks.values()):
        raise AssertionError(
            "the canonical P5 principal-filter certificate failed"
        )
    if not upward_failures:
        raise AssertionError(
            "the canonical P5 interval has no upward validity failure"
        )
    interval_profile = {
        "optional_variables": len(optional_indices),
        "distinct_realizations": len(realizations),
        "valid_realizations": len(valid_assignments),
        "invalid_realizations": len(invalid_assignments),
        "minimal_valid_realizations": len(
            minimal_valid_assignments
        ),
        "upward_failure_relations": len(upward_failures),
    }
    expected_interval_profile = {
        "optional_variables": 4,
        "distinct_realizations": 6,
        "valid_realizations": 4,
        "invalid_realizations": 2,
        "minimal_valid_realizations": 1,
        "upward_failure_relations": 5,
    }
    if interval_profile != expected_interval_profile:
        raise AssertionError(
            "the canonical P5 Horn interval profile changed: "
            f"{interval_profile}"
        )
    minimal_base, minimal_site = decode_pair(
        minimal_valid_assignments[0]
    )
    if (minimal_base, minimal_site) != (
        int("8fb3af", 16),
        site,
    ):
        raise AssertionError(
            "the canonical P5 interval has an unexpected admissible minimum"
        )
    for realization in invalid_assignments:
        current_base, current_site = decode_pair(realization)
        profile = collar_profile(current_base, current_site)
        if (
            profile["base_partial_cube"]
            and profile["site_partial_cube"]
            and profile["base_irredundant"]
            and profile["base_rooted_r_per"]
            and profile["marked_edge"]
            and not profile["site_r_per"]
        ):
            continue
        raise AssertionError(
            "a P5 interval realization fails outside site R_per membership"
        )

    added_site_vertex = (invalid_site ^ site).bit_length() - 1
    atlas_extension_witnesses: list[dict[str, object]] = []
    for coordinate, defect_system in enumerate(site_atlas):
        projection = project_vertex(
            added_site_vertex,
            coordinate,
        )
        if added_site_vertex & (1 << coordinate):
            root_mate = added_site_vertex ^ (1 << coordinate)
            if not (site & (1 << root_mate)):
                raise AssertionError(
                    "an opposite-side P5 extension has no root mate"
                )
            atlas_extension_witnesses.append(
                {
                    "coordinate": coordinate,
                    "side": "opposite",
                    "projection": format(
                        projection,
                        f"0{dimension - 1}b",
                    ),
                    "root_mate": format(
                        root_mate,
                        f"0{dimension}b",
                    ),
                }
            )
            continue
        repair_witnesses = []
        for defect, clause in defect_system:
            support = next(
                (
                    candidate
                    for candidate in clause
                    if not (
                        candidate
                        & ~(defect ^ projection)
                    )
                ),
                None,
            )
            if support is None:
                raise AssertionError(
                    "a root-side P5 extension is not atlas-allowed"
                )
            repair_witnesses.append(
                {
                    "defect": format(
                        defect,
                        f"0{dimension - 1}b",
                    ),
                    "support": format(
                        support,
                        f"0{dimension - 1}b",
                    ),
                }
            )
        atlas_extension_witnesses.append(
            {
                "coordinate": coordinate,
                "side": "root",
                "projection": format(
                    projection,
                    f"0{dimension - 1}b",
                ),
                "repair_witnesses": repair_witnesses,
            }
        )

    marked_minor_rows: list[dict[str, object]] = []

    def append_marked_state(
        operation: str,
        current_base: int,
        current_site: int,
        current_dimension: int,
        current_q: int,
    ) -> None:
        closure = common_availability_closure_avoiding(
            current_base,
            current_site,
            current_dimension,
            0,
            1 << current_q,
        )[0]
        marked_minor_rows.append(
            {
                "operation": operation,
                "dimension": current_dimension,
                "activation_coordinate": current_q,
                "base_order": current_base.bit_count(),
                "site_order": current_site.bit_count(),
                "maximal_q_free_closure": format(
                    closure,
                    f"0{current_dimension}b",
                ),
                "activated": is_q_activated(
                    current_base,
                    current_site,
                    current_dimension,
                    0,
                    current_q,
                ),
            }
        )

    append_marked_state(
        "initial",
        base,
        site,
        dimension,
        q,
    )
    for coordinate in range(dimension):
        if coordinate == q:
            continue
        projected_q = q - (q > coordinate)
        append_marked_state(
            f"contract_{coordinate}",
            contract_coordinate(
                base,
                dimension,
                coordinate,
            ),
            contract_coordinate(
                site,
                dimension,
                coordinate,
            ),
            dimension - 1,
            projected_q,
        )
        append_marked_state(
            f"restrict_root_{coordinate}",
            section(
                base,
                dimension,
                coordinate,
                0,
            ),
            section(
                site,
                dimension,
                coordinate,
                0,
            ),
            dimension - 1,
            projected_q,
        )
    if (
        marked_minor_rows[0]["activated"]
        or not all(
            row["activated"]
            for row in marked_minor_rows[1:]
        )
    ):
        raise AssertionError(
            "the P5 marked-minor activation certificate failed"
        )

    payload = {
        "schema_version": 1,
        "kind": "defect-atlas-principal-filter-falsifier",
        "complete_rooted_family_audit_through_dimension": 4,
        "dimension_profiles": dimension_profiles,
        "partial_cube_upward_theorem_regression_failures": 0,
        "canonical_falsifier": {
            "dimension": dimension,
            "root": root,
            "activation_coordinate": q,
            "base_mask": format(base, "x"),
            "valid_site_mask": format(site, "x"),
            "invalid_site_mask": format(invalid_site, "x"),
            "base_order": base.bit_count(),
            "valid_site_order": site.bit_count(),
            "invalid_site_order": invalid_site.bit_count(),
            "added_site_vertex": format(
                added_site_vertex,
                f"0{dimension}b",
            ),
            "invalid_site_family": "pentagonal P5 obstruction",
            "exact_checks": exact_checks,
            "base_rooted_r_per_certificate": (
                rooted_r_per_certificate(
                    base,
                    dimension,
                    root,
                )
            ),
            "valid_site_r_per_certificate": r_per_certificate(
                site,
                dimension,
            ),
            "atlas_extension_witnesses": (
                atlas_extension_witnesses
            ),
            "marked_minor_rows": marked_minor_rows,
            "horn_interval": {
                **interval_profile,
                "rows": interval_rows,
            },
        },
    }
    Path(output_path).write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print("\nDefect-atlas principal-filter falsifier")
    print(
        "  canonical marked triple: "
        f"({base:x},{site:x},{invalid_site:x};q={q})"
    )
    print(
        "  Horn interval: "
        f"optional_variables={len(optional_indices)}; "
        f"realizations={len(realizations)}; "
        f"valid={len(valid_assignments)}; "
        f"invalid={len(invalid_assignments)}; "
        "minimal_valid="
        f"{len(minimal_valid_assignments)}; "
        f"upward_failures={len(upward_failures)}"
    )
    print(
        "  peak_rss_mb="
        f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024:.1f}; "
        f"runtime_seconds={time.perf_counter() - started:.3f}"
    )
    print(f"  wrote={output_path}")


def complete_site_three_critical_pair_profile_dimension_five(
    output_path: str,
) -> None:
    """Exhaust the three-vertex-site branch on five old coordinates.

    The exact rooted-good ``Q_5`` generator supplies every labelled base.
    A three-vertex critical site must be an endpoint-rooted path whose first
    coordinate is the base's unique available coordinate.  The scan compares
    the direct forced-minor predicate with the lower-dimensional collar
    criterion on every such contained presentation.
    """

    started = time.perf_counter()
    dimension = 5
    families, generator_presentations = (
        rooted_r_per_families_dimension_five(32)
    )
    universe = (1 << (1 << dimension)) - 1
    low_masks = coordinate_low_vertex_masks(dimension)
    high_masks = tuple(universe ^ mask for mask in low_masks)

    def root_zero_availability_mask(family: int) -> int:
        answer = 0
        for coordinate, (low_mask, high_mask) in enumerate(
            zip(low_masks, high_masks)
        ):
            shift = 1 << coordinate
            if nested(
                (family & high_mask) >> shift,
                family & low_mask,
            ):
                answer |= 1 << coordinate
        return answer

    irredundant_bases = 0
    unique_availability_bases = 0
    contained_path_presentations = 0
    critical_presentations = 0
    criterion_disagreements: list[dict[str, object]] = []
    critical_pairs: set[tuple[int, int]] = set()

    for family_index, base in enumerate(families, start=1):
        if not all(
            base & low_mask and base & high_mask
            for low_mask, high_mask in zip(low_masks, high_masks)
        ):
            continue
        irredundant_bases += 1
        available_mask = root_zero_availability_mask(base)
        if available_mask.bit_count() != 1:
            continue
        unique_availability_bases += 1
        p = available_mask.bit_length() - 1
        for q in range(dimension):
            if q == p:
                continue
            site = family_from_vertices(
                (0, 1 << p, (1 << p) | (1 << q))
            )
            if site & ~base:
                continue
            contained_path_presentations += 1
            direct = is_critical_reverse_pair_by_closure(
                base,
                site,
                dimension,
                0,
            )
            collar = (
                is_three_vertex_collar_critical_by_lower_activation(
                    base,
                    site,
                    dimension,
                    0,
                )
            )
            if direct != collar:
                criterion_disagreements.append(
                    {
                        "base_mask": format(base, "x"),
                        "site_mask": format(site, "x"),
                        "direct": direct,
                        "collar": collar,
                    }
                )
            if direct:
                critical_presentations += 1
                pair = canonical_reverse_pair(
                    base,
                    site,
                    dimension,
                )
                critical_pairs.add(pair)
        if family_index % 250_000 == 0:
            print(
                "  site-three Q5 progress: "
                f"families={family_index}/{len(families)}; "
                f"irredundant_bases={irredundant_bases}; "
                "unique_availability_bases="
                f"{unique_availability_bases}; "
                f"contained_paths={contained_path_presentations}; "
                f"critical_types={len(critical_pairs)}",
                flush=True,
            )

    if criterion_disagreements:
        raise AssertionError(
            "the direct and collar criteria disagree on "
            f"{len(criterion_disagreements)} Q5 presentations"
        )

    ordered_pairs = sorted(critical_pairs)
    pair_digest = hashlib.sha256(
        (
            "\n".join(
                f"{base:x}:{site:x}"
                for base, site in ordered_pairs
            )
            + "\n"
        ).encode("ascii")
    ).hexdigest()

    def pair_payload(
        pair: tuple[int, int],
    ) -> dict[str, object]:
        base, site = pair
        path_coordinates = (
            endpoint_rooted_three_vertex_site_coordinates(
                site,
                dimension,
                0,
            )
        )
        if path_coordinates is None:
            raise AssertionError("a critical site is not endpoint-rooted")
        p, q = path_coordinates
        lower_labels = tuple(
            coordinate
            for coordinate in range(dimension)
            if coordinate != p
        )
        lower_q = lower_labels.index(q)
        rooted_side = section(base, dimension, p, 0)
        opposite_side = section(base, dimension, p, 1)
        lower_dimension = dimension - 1

        def labelled_availability(family: int) -> list[int]:
            return [
                lower_labels[coordinate]
                for coordinate in (
                    raw_root_avoiding_peripheral_coordinates(
                        family,
                        lower_dimension,
                        0,
                    )
                )
            ]

        activation_rows: list[dict[str, object]] = []
        for coordinate, label in enumerate(lower_labels):
            row: dict[str, object] = {
                "coordinate": label,
                "rooted_side_restriction_rooted_good": (
                    is_rooted_r_per(
                        section(
                            rooted_side,
                            lower_dimension,
                            coordinate,
                            0,
                        ),
                        lower_dimension - 1,
                        0,
                    )
                ),
            }
            if coordinate != lower_q:
                projected_q = lower_q - (lower_q > coordinate)
                contraction_closure = common_availability_closure(
                    contract_coordinate(
                        rooted_side,
                        lower_dimension,
                        coordinate,
                    ),
                    contract_coordinate(
                        opposite_side,
                        lower_dimension,
                        coordinate,
                    ),
                    lower_dimension - 1,
                    0,
                )[0]
                restriction_closure = common_availability_closure(
                    section(
                        rooted_side,
                        lower_dimension,
                        coordinate,
                        0,
                    ),
                    section(
                        opposite_side,
                        lower_dimension,
                        coordinate,
                        0,
                    ),
                    lower_dimension - 1,
                    0,
                )[0]
                row.update(
                    {
                        "q_activated_after_contraction": bool(
                            contraction_closure
                            & (1 << projected_q)
                        ),
                        "q_activated_after_root_restriction": bool(
                            restriction_closure
                            & (1 << projected_q)
                        ),
                    }
                )
            activation_rows.append(row)

        lower_initial_closure = common_availability_closure(
            rooted_side,
            opposite_side,
            lower_dimension,
            0,
        )[0]
        return {
            "base_mask": format(base, "x"),
            "site_mask": format(site, "x"),
            "base_order": base.bit_count(),
            "site_order": site.bit_count(),
            "factor_order": base.bit_count() + site.bit_count(),
            "collar_coordinate": p,
            "activation_coordinate": q,
            "lower_rooted_side_mask": format(rooted_side, "x"),
            "lower_opposite_side_mask": format(opposite_side, "x"),
            "lower_rooted_side_order": rooted_side.bit_count(),
            "lower_opposite_side_order": opposite_side.bit_count(),
            "lower_rooted_side_available": (
                labelled_availability(rooted_side)
            ),
            "lower_opposite_side_available": (
                labelled_availability(opposite_side)
            ),
            "lower_initial_closure": [
                lower_labels[coordinate]
                for coordinate in range(lower_dimension)
                if lower_initial_closure >> coordinate & 1
            ],
            "lower_forced_minor_activation": activation_rows,
        }

    payload = {
        "schema_version": 1,
        "kind": "complete-rooted-site-three-critical-pair-profile-d5",
        "complete": True,
        "old_dimension": dimension,
        "rooted_family_count": len(families),
        "rooted_family_generator_presentations": generator_presentations,
        "irredundant_base_count": irredundant_bases,
        "unique_availability_base_count": unique_availability_bases,
        "contained_endpoint_path_presentations": (
            contained_path_presentations
        ),
        "criterion_disagreements": len(criterion_disagreements),
        "critical_presentations": critical_presentations,
        "critical_pair_types": len(ordered_pairs),
        "critical_base_order_histogram": dict(
            sorted(
                Counter(
                    base.bit_count()
                    for base, _site in ordered_pairs
                ).items()
            )
        ),
        "pair_digest": pair_digest,
        "pair_types": [
            pair_payload(pair)
            for pair in ordered_pairs
        ],
    }
    Path(output_path).write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print("\nComplete three-vertex-site profile on five old coordinates")
    print(
        f"  rooted_families={len(families)}; "
        f"irredundant_bases={irredundant_bases}; "
        f"unique_availability_bases={unique_availability_bases}; "
        f"contained_paths={contained_path_presentations}"
    )
    print(
        f"  critical_presentations={critical_presentations}; "
        f"critical_pair_types={len(ordered_pairs)}; "
        f"criterion_disagreements={len(criterion_disagreements)}"
    )
    print(f"  pair_digest={pair_digest}")
    print(
        "  peak_rss_mb="
        f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024:.1f}; "
        f"runtime_seconds={time.perf_counter() - started:.3f}"
    )
    print(f"  wrote={output_path}")


def site_three_activation_antimatroid_profile(
    input_path: str,
    output_path: str,
) -> None:
    """Profile activation invariants on the complete Q5 site-three branch."""

    started = time.perf_counter()
    source = json.loads(Path(input_path).read_text(encoding="utf-8"))
    if not isinstance(source, dict) or source.get("kind") != (
        "complete-rooted-site-three-critical-pair-profile-d5"
    ):
        raise ValueError(
            "the input is not a complete Q5 site-three certificate"
        )
    if not source.get("complete") or source.get("old_dimension") != 5:
        raise ValueError("the site-three source must be complete in Q5")
    raw_rows = source.get("pair_types")
    if not isinstance(raw_rows, list):
        raise ValueError("the site-three source has no pair_types array")

    dimension = 5
    pairs: list[tuple[int, int]] = []
    rows_by_pair: dict[tuple[int, int], dict[str, object]] = {}
    activation_by_pair: dict[tuple[int, int], object] = {}
    repair_by_pair: dict[tuple[int, int], object] = {}
    for row in raw_rows:
        if not isinstance(row, dict):
            raise ValueError("every site-three pair row must be an object")
        pair = (
            int(str(row["base_mask"]), 16),
            int(str(row["site_mask"]), 16),
        )
        if canonical_reverse_pair(*pair, dimension) != pair:
            raise AssertionError("a site-three source pair is not canonical")
        if not (
            is_critical_reverse_pair_by_closure(
                *pair,
                dimension,
                0,
            )
            and is_three_vertex_collar_critical_by_lower_activation(
                *pair,
                dimension,
                0,
            )
        ):
            raise AssertionError(
                "a site-three source pair failed a criticality predicate"
            )
        pairs.append(pair)
        rows_by_pair[pair] = row
        activation_by_pair[pair] = (
            canonical_activation_antimatroid_signature(
                *pair,
                dimension,
            )
        )
        repair_by_pair[pair] = canonical_repair_signature(
            *pair,
            dimension,
        )

    pair_digest = hashlib.sha256(
        (
            "\n".join(
                f"{base:x}:{site:x}"
                for base, site in pairs
            )
            + "\n"
        ).encode("ascii")
    ).hexdigest()
    if pair_digest != source.get("pair_digest"):
        raise AssertionError("the site-three source pair digest changed")
    pair_index = {
        pair: index
        for index, pair in enumerate(pairs)
    }

    def signature_groups(
        signatures: dict[tuple[int, int], object],
        extra_key,
    ) -> dict[object, list[tuple[int, int]]]:
        groups: dict[object, list[tuple[int, int]]] = {}
        for pair in pairs:
            groups.setdefault(
                (
                    pair[0].bit_count(),
                    pair[1].bit_count(),
                    extra_key(pair),
                    signatures[pair],
                ),
                [],
            ).append(pair)
        return groups

    activation_groups = signature_groups(
        activation_by_pair,
        lambda _pair: (),
    )
    repair_groups = signature_groups(
        repair_by_pair,
        lambda _pair: (),
    )
    collar_activation_groups = signature_groups(
        activation_by_pair,
        lambda pair: (
            int(rows_by_pair[pair]["lower_rooted_side_order"]),
            int(rows_by_pair[pair]["lower_opposite_side_order"]),
        ),
    )
    activation_collisions = [
        group
        for group in activation_groups.values()
        if len(group) > 1
    ]
    repair_collisions = [
        group
        for group in repair_groups.values()
        if len(group) > 1
    ]
    collar_activation_collisions = [
        group
        for group in collar_activation_groups.values()
        if len(group) > 1
    ]

    signature_payloads = [
        activation_antimatroid_signature_payload(
            activation_by_pair[pair]
        )
        for pair in pairs
    ]
    signature_digest = hashlib.sha256(
        (
            "\n".join(
                json.dumps(
                    payload,
                    sort_keys=True,
                    separators=(",", ":"),
                )
                for payload in signature_payloads
            )
            + "\n"
        ).encode("utf-8")
    ).hexdigest()

    def collision_payload(
        group: list[tuple[int, int]],
    ) -> dict[str, object]:
        answer: dict[str, object] = {
            "base_order": group[0][0].bit_count(),
            "site_order": group[0][1].bit_count(),
            "factor_order": (
                group[0][0].bit_count()
                + group[0][1].bit_count()
            ),
            "pair_type_indices": [
                pair_index[pair]
                for pair in group
            ],
            "pairs": [
                {
                    "base_mask": format(pair[0], "x"),
                    "site_mask": format(pair[1], "x"),
                    "lower_rooted_side_order": int(
                        rows_by_pair[pair][
                            "lower_rooted_side_order"
                        ]
                    ),
                    "lower_opposite_side_order": int(
                        rows_by_pair[pair][
                            "lower_opposite_side_order"
                        ]
                    ),
                }
                for pair in group
            ],
        }
        if len(group) == 2:
            left_vertices = set(vertices(group[0][0], dimension))
            right_vertices = set(vertices(group[1][0], dimension))
            answer["base_symmetric_difference"] = {
                "first_only": [
                    format(vertex, f"0{dimension}b")
                    for vertex in sorted(
                        left_vertices - right_vertices
                    )
                ],
                "second_only": [
                    format(vertex, f"0{dimension}b")
                    for vertex in sorted(
                        right_vertices - left_vertices
                    )
                ],
            }
        return answer

    payload = {
        "schema_version": 1,
        "kind": (
            "complete-rooted-site-three-activation-profile-d5"
        ),
        "source_certificate": input_path,
        "source_pair_digest": pair_digest,
        "old_dimension": dimension,
        "critical_pair_types": len(pairs),
        "order_enriched_repair_signature_types": len(repair_groups),
        "order_enriched_repair_collision_groups": len(
            repair_collisions
        ),
        "order_enriched_activation_antimatroid_signature_types": len(
            activation_groups
        ),
        "order_enriched_activation_antimatroid_collision_groups": len(
            activation_collisions
        ),
        "collar_order_enriched_activation_antimatroid_signature_types": (
            len(collar_activation_groups)
        ),
        "collar_order_enriched_activation_antimatroid_collision_groups": (
            len(collar_activation_collisions)
        ),
        "activation_antimatroid_signature_digest": signature_digest,
        "order_enriched_activation_collision_classes": [
            {
                **collision_payload(group),
                "canonical_activation_antimatroid": (
                    activation_antimatroid_signature_payload(
                        activation_by_pair[group[0]]
                    )
                ),
            }
            for group in sorted(
                activation_collisions,
                key=lambda members: (
                    members[0][0].bit_count()
                    + members[0][1].bit_count(),
                    tuple(
                        pair_index[pair]
                        for pair in members
                    ),
                ),
            )
        ],
        "order_enriched_repair_collision_classes": [
            {
                **collision_payload(group),
                "canonical_repair_signature": (
                    repair_signature_payload(
                        repair_by_pair[group[0]]
                    )
                ),
            }
            for group in sorted(
                repair_collisions,
                key=lambda members: (
                    members[0][0].bit_count()
                    + members[0][1].bit_count(),
                    tuple(
                        pair_index[pair]
                        for pair in members
                    ),
                ),
            )
        ],
    }
    Path(output_path).write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print("\nComplete Q5 site-three activation profile")
    print(
        f"  critical_pair_types={len(pairs)}; "
        "order_enriched_activation_signatures="
        f"{len(activation_groups)}; "
        "activation_collision_groups="
        f"{len(activation_collisions)}; "
        "collar_order_enriched_activation_signatures="
        f"{len(collar_activation_groups)}; "
        "collar_activation_collision_groups="
        f"{len(collar_activation_collisions)}"
    )
    print(
        "  order_enriched_repair_signatures="
        f"{len(repair_groups)}; "
        f"repair_collision_groups={len(repair_collisions)}"
    )
    print(f"  source_pair_digest={pair_digest}")
    print(f"  activation_antimatroid_signature_digest={signature_digest}")
    print(
        "  peak_rss_mb="
        f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024:.1f}; "
        f"runtime_seconds={time.perf_counter() - started:.3f}"
    )
    print(f"  wrote={output_path}")


def complete_rooted_critical_pair_profile_dimension_four(
    output_path: str,
) -> None:
    """Exhaust all critical reverse pairs on four old coordinates.

    This bypasses five-dimensional expansion generation.  The base is
    required to use all four old coordinates; the site may use a proper
    subset.  The rooted-family generator is complete up to labelled cube
    embeddings, and canonical reverse pairs quotient the remaining
    coordinate permutations.
    """

    started = time.perf_counter()
    dimension = 4
    rooted_families = rooted_r_per_families_in_cube(dimension)
    bases = tuple(
        family
        for family in rooted_families
        if irredundant_projection(family, dimension)[1] == dimension
    )

    containment_presentations = 0
    critical_presentations = 0
    critical_pairs: set[tuple[int, int]] = set()
    for base in bases:
        for site in rooted_families:
            if site == base or site & ~base:
                continue
            containment_presentations += 1
            if not is_critical_reverse_pair(
                base,
                site,
                dimension,
                0,
            ):
                continue
            critical_presentations += 1
            pair = canonical_reverse_pair(base, site, dimension)
            factor = lift_expansion(site, base, dimension)
            if not is_root_critical(factor, dimension + 1, 0):
                raise AssertionError(
                    "a direct critical pair did not lift to a root-critical factor"
                )
            critical_pairs.add(pair)

    pairs = sorted(critical_pairs)
    pair_index = {pair: index for index, pair in enumerate(pairs)}
    repair_by_pair = {
        pair: canonical_repair_signature(*pair, dimension)
        for pair in pairs
    }
    activation_by_pair = {
        pair: canonical_activation_antimatroid_signature(*pair, dimension)
        for pair in pairs
    }

    def signature_classes(
        signatures: dict[tuple[int, int], object],
        *,
        with_orders: bool,
    ) -> list[list[tuple[int, int]]]:
        classes: dict[object, list[tuple[int, int]]] = {}
        for pair in pairs:
            key: object = signatures[pair]
            if with_orders:
                key = (
                    pair[0].bit_count(),
                    pair[1].bit_count(),
                    key,
                )
            classes.setdefault(key, []).append(pair)
        return sorted(
            classes.values(),
            key=lambda members: (
                -len(members),
                tuple(pair_index[pair] for pair in members),
            ),
        )

    repair_classes = signature_classes(
        repair_by_pair,
        with_orders=False,
    )
    ordered_repair_classes = signature_classes(
        repair_by_pair,
        with_orders=True,
    )
    activation_classes = signature_classes(
        activation_by_pair,
        with_orders=False,
    )
    ordered_activation_classes = signature_classes(
        activation_by_pair,
        with_orders=True,
    )
    repair_collisions = [
        members for members in repair_classes if len(members) > 1
    ]
    ordered_repair_collisions = [
        members for members in ordered_repair_classes if len(members) > 1
    ]
    activation_collisions = [
        members for members in activation_classes if len(members) > 1
    ]
    ordered_activation_collisions = [
        members
        for members in ordered_activation_classes
        if len(members) > 1
    ]

    pair_digest = hashlib.sha256(
        (
            "\n".join(f"{base:x}:{site:x}" for base, site in pairs)
            + "\n"
        ).encode("ascii")
    ).hexdigest()
    bounded_pairs = [
        pair
        for pair in pairs
        if pair[0].bit_count() + pair[1].bit_count() <= 22
    ]
    bounded_digest = hashlib.sha256(
        (
            "\n".join(
                f"{base:x}:{site:x}" for base, site in bounded_pairs
            )
            + "\n"
        ).encode("ascii")
    ).hexdigest()
    expected_order_22_digest = (
        "4e9f3db6e7ea70ea964093c96ae56ef287732a13a183d73f1f57389cbabf39b7"
    )
    if bounded_digest != expected_order_22_digest:
        raise AssertionError(
            "the direct pair census disagrees with the independent "
            "order-22 expansion census: "
            f"expected={expected_order_22_digest}, observed={bounded_digest}"
        )

    factor_order_histogram = Counter(
        base.bit_count() + site.bit_count()
        for base, site in pairs
    )
    signature_digest = hashlib.sha256(
        (
            "\n".join(
                json.dumps(
                    activation_antimatroid_signature_payload(
                        activation_by_pair[pair]
                    ),
                    sort_keys=True,
                    separators=(",", ":"),
                )
                for pair in pairs
            )
            + "\n"
        ).encode("utf-8")
    ).hexdigest()
    observed_complete = (
        len(rooted_families),
        len(bases),
        containment_presentations,
        critical_presentations,
        len(pairs),
        len(repair_classes),
        len(ordered_repair_classes),
        len(ordered_repair_collisions),
        len(activation_classes),
        len(ordered_activation_classes),
        len(ordered_activation_collisions),
        pair_digest,
        signature_digest,
    )
    expected_complete = (
        2753,
        2532,
        342630,
        1152,
        56,
        42,
        53,
        3,
        44,
        56,
        0,
        "e34e706c161f70bc451ac700663369cc51d230a00832941d84465db32eb4d747",
        "1e3df177f52b5663248d9028bc726e43bf07376204bebb29a08488df4becddc6",
    )
    if observed_complete != expected_complete:
        raise AssertionError(
            "the complete four-coordinate critical-pair regression changed: "
            f"expected={expected_complete}, observed={observed_complete}"
        )
    payload = {
        "schema_version": 1,
        "kind": "complete-rooted-critical-pair-profile",
        "old_dimension": dimension,
        "complete": True,
        "rooted_r_per_labelled_families": len(rooted_families),
        "irredundant_base_labelled_families": len(bases),
        "proper_containment_presentations": containment_presentations,
        "critical_labelled_presentations": critical_presentations,
        "critical_pair_types": len(pairs),
        "factor_order_histogram": dict(sorted(factor_order_histogram.items())),
        "pair_digest": pair_digest,
        "order_22_pair_types": len(bounded_pairs),
        "order_22_pair_digest": bounded_digest,
        "repair_signature_types": len(repair_classes),
        "repair_signature_collision_groups": len(repair_collisions),
        "order_enriched_repair_signature_types": len(
            ordered_repair_classes
        ),
        "order_enriched_repair_collision_groups": len(
            ordered_repair_collisions
        ),
        "activation_antimatroid_signature_types": len(
            activation_classes
        ),
        "activation_antimatroid_collision_groups": len(
            activation_collisions
        ),
        "order_enriched_activation_antimatroid_signature_types": len(
            ordered_activation_classes
        ),
        "order_enriched_activation_antimatroid_collision_groups": len(
            ordered_activation_collisions
        ),
        "activation_antimatroid_signature_digest": signature_digest,
        "pair_types": [
            {
                "index": pair_index[pair],
                "base_mask": format(pair[0], "x"),
                "site_mask": format(pair[1], "x"),
                "base_order": pair[0].bit_count(),
                "site_order": pair[1].bit_count(),
                "factor_order": (
                    pair[0].bit_count() + pair[1].bit_count()
                ),
                "canonical_repair_signature": repair_signature_payload(
                    repair_by_pair[pair]
                ),
                "canonical_activation_antimatroid": (
                    activation_antimatroid_signature_payload(
                        activation_by_pair[pair]
                    )
                ),
            }
            for pair in pairs
        ],
        "order_enriched_activation_collision_classes": [
            {
                "base_order": members[0][0].bit_count(),
                "site_order": members[0][1].bit_count(),
                "pair_type_indices": [
                    pair_index[pair] for pair in members
                ],
                "canonical_activation_antimatroid": (
                    activation_antimatroid_signature_payload(
                        activation_by_pair[members[0]]
                    )
                ),
            }
            for members in ordered_activation_collisions
        ],
    }
    Path(output_path).write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    print("\nComplete rooted critical-pair profile on four old coordinates")
    print(
        f"  rooted_families={len(rooted_families)}; "
        f"irredundant_bases={len(bases)}; "
        f"proper_containments={containment_presentations}; "
        f"critical_presentations={critical_presentations}; "
        f"critical_pair_types={len(pairs)}"
    )
    print(
        f"  repair_signatures={len(repair_classes)}; "
        "order_enriched_repair_signatures="
        f"{len(ordered_repair_classes)}; "
        "order_enriched_repair_collisions="
        f"{len(ordered_repair_collisions)}"
    )
    print(
        "  activation_antimatroid_signatures="
        f"{len(activation_classes)}; "
        "order_enriched_activation_antimatroid_signatures="
        f"{len(ordered_activation_classes)}; "
        "order_enriched_activation_antimatroid_collisions="
        f"{len(ordered_activation_collisions)}"
    )
    print(f"  factor_order_histogram={dict(sorted(factor_order_histogram.items()))}")
    print(f"  pair_digest={pair_digest}")
    print(f"  activation_antimatroid_signature_digest={signature_digest}")
    print(
        "  peak_rss_mb="
        f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024:.1f}; "
        f"runtime_seconds={time.perf_counter() - started:.3f}"
    )
    print(f"  wrote={output_path}")


def rooted_repair_profile_dimension_five(
    order_limit: int,
    output_path: str | None,
) -> None:
    """Profile first-step repair data for exact dimension-five critical pairs."""

    if order_limit < 6 or order_limit > 22:
        raise ValueError(
            "the dimension-five repair profile supports orders 6 through 22"
        )

    started = time.perf_counter()
    dimension_four_bases = [
        family
        for family in enumerate_partial_cubes(4)[4]
        if family.bit_count() < order_limit
        and shattering_profile(family, 4)[2]
    ]
    raw = expansion_candidates(
        dimension_four_bases,
        4,
        order_limit,
    )
    r_types = {
        canonical_cube_family(family, 5)
        for family in raw
        if shattering_profile(family, 5)[2]
        and is_r_per(family, 5)
    }
    pointed = {
        canonical_pointed_family(family, 5, root)
        for family in r_types
        for root in vertices(family, 5)
    }
    critical = sorted(
        (
            family
            for family in pointed
            if is_root_critical(family, 5, 0)
        ),
        key=lambda family: (family.bit_count(), family),
    )

    occurrences_by_pair: dict[
        tuple[int, int],
        list[dict[str, object]],
    ] = {}
    signature_by_pair: dict[
        tuple[int, int],
        tuple[int, int, tuple[tuple[int, bool, int, bool], ...]],
    ] = {}
    presentation_count = 0
    for factor in critical:
        for peripheral_coordinate, base, site, lower_dimension in (
            peripheral_decompositions(factor, 5)
        ):
            if 0 not in vertices(site, lower_dimension):
                raise AssertionError(
                    "a root-critical periphery omitted its root"
                )
            if not is_critical_reverse_pair(
                base,
                site,
                lower_dimension,
                0,
            ):
                raise AssertionError(
                    "a root-critical presentation failed the exact "
                    "critical-pair test"
                )
            pair = canonical_reverse_pair(
                base,
                site,
                lower_dimension,
            )
            signature = canonical_repair_signature(
                pair[0],
                pair[1],
                lower_dimension,
            )
            previous = signature_by_pair.setdefault(pair, signature)
            if previous != signature:
                raise AssertionError(
                    "one pair type has two repair signatures"
                )
            occurrences_by_pair.setdefault(pair, []).append(
                {
                    "factor_mask": format(factor, "x"),
                    "factor_order": factor.bit_count(),
                    "peripheral_coordinate": peripheral_coordinate,
                }
            )
            presentation_count += 1

    ordered_pairs = sorted(occurrences_by_pair)
    pair_index = {
        pair: index for index, pair in enumerate(ordered_pairs)
    }
    pairs_by_signature: dict[
        tuple[int, int, tuple[tuple[int, bool, int, bool], ...]],
        list[tuple[int, int]],
    ] = {}
    for pair in ordered_pairs:
        pairs_by_signature.setdefault(
            signature_by_pair[pair],
            [],
        ).append(pair)
    collision_groups = [
        pairs
        for pairs in pairs_by_signature.values()
        if len(pairs) > 1
    ]
    largest_collision = max(
        (len(pairs) for pairs in collision_groups),
        default=0,
    )
    pairs_by_ordered_signature: dict[
        tuple[
            int,
            int,
            tuple[
                int,
                int,
                tuple[tuple[int, bool, int, bool], ...],
            ],
        ],
        list[tuple[int, int]],
    ] = {}
    for pair in ordered_pairs:
        pairs_by_ordered_signature.setdefault(
            (
                pair[0].bit_count(),
                pair[1].bit_count(),
                signature_by_pair[pair],
            ),
            [],
        ).append(pair)
    ordered_collision_groups = [
        pairs
        for pairs in pairs_by_ordered_signature.values()
        if len(pairs) > 1
    ]
    largest_ordered_collision = max(
        (len(pairs) for pairs in ordered_collision_groups),
        default=0,
    )

    pair_digest_lines = [
        f"{base:x}:{site:x}" for base, site in ordered_pairs
    ]
    pair_digest = hashlib.sha256(
        ("\n".join(pair_digest_lines) + "\n").encode("ascii")
    ).hexdigest()
    payload: dict[str, object] = {
        "schema_version": 1,
        "kind": "rooted-critical-repair-profile",
        "dimension": 5,
        "order_limit": order_limit,
        "max_factor_order": order_limit,
        "old_dimension": 4,
        "critical_factor_types": len(critical),
        "peripheral_presentations": presentation_count,
        "critical_pair_types": len(ordered_pairs),
        "repair_signature_types": len(pairs_by_signature),
        "repair_signature_collision_groups": len(collision_groups),
        "largest_collision_group": largest_collision,
        "order_enriched_repair_signature_types": len(
            pairs_by_ordered_signature
        ),
        "order_enriched_collision_groups": len(
            ordered_collision_groups
        ),
        "largest_order_enriched_collision_group": (
            largest_ordered_collision
        ),
        "critical_pair_digest": pair_digest,
        "pair_types": [
            {
                "index": pair_index[pair],
                "base_mask": format(pair[0], "x"),
                "site_mask": format(pair[1], "x"),
                "base_order": pair[0].bit_count(),
                "site_order": pair[1].bit_count(),
                "canonical_repair_signature": repair_signature_payload(
                    signature_by_pair[pair]
                ),
                "occurrences": sorted(
                    occurrences_by_pair[pair],
                    key=lambda item: (
                        int(str(item["factor_mask"]), 16),
                        int(item["peripheral_coordinate"]),
                    ),
                ),
            }
            for pair in ordered_pairs
        ],
        "collision_groups": [
            {
                "pair_type_indices": [
                    pair_index[pair] for pair in pairs
                ],
                "canonical_repair_signature": repair_signature_payload(
                    signature_by_pair[pairs[0]]
                ),
            }
            for pairs in sorted(
                collision_groups,
                key=lambda group: (
                    -len(group),
                    tuple(pair_index[pair] for pair in group),
                ),
            )
        ],
        "order_enriched_collision_groups": [
            {
                "base_order": pairs[0][0].bit_count(),
                "site_order": pairs[0][1].bit_count(),
                "pair_type_indices": [
                    pair_index[pair] for pair in pairs
                ],
                "canonical_repair_signature": repair_signature_payload(
                    signature_by_pair[pairs[0]]
                ),
            }
            for pairs in sorted(
                ordered_collision_groups,
                key=lambda group: (
                    -len(group),
                    tuple(pair_index[pair] for pair in group),
                ),
            )
        ],
    }

    if order_limit == 16:
        expected = (14, 18, 16, 8, 4, 6, 14, 2, 2)
        observed = (
            len(critical),
            presentation_count,
            len(ordered_pairs),
            len(pairs_by_signature),
            len(collision_groups),
            largest_collision,
            len(pairs_by_ordered_signature),
            len(ordered_collision_groups),
            largest_ordered_collision,
        )
        if observed != expected:
            raise AssertionError(
                "dimension-five repair-profile regression changed: "
                f"expected={expected}, observed={observed}"
            )
    if order_limit == 21:
        expected = (28, 41, 37, 25, 5, 8, 34, 3, 2)
        observed = (
            len(critical),
            presentation_count,
            len(ordered_pairs),
            len(pairs_by_signature),
            len(collision_groups),
            largest_collision,
            len(pairs_by_ordered_signature),
            len(ordered_collision_groups),
            largest_ordered_collision,
        )
        if observed != expected:
            raise AssertionError(
                "dimension-five order-21 repair-profile regression changed: "
                f"expected={expected}, observed={observed}"
            )

    if output_path is not None:
        Path(output_path).write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print("\nExact dimension-five rooted repair profile")
    print(
        f"  order_limit={order_limit}; "
        f"critical_factors={len(critical)}; "
        f"peripheral_presentations={presentation_count}; "
        f"critical_pair_types={len(ordered_pairs)}; "
        f"repair_signatures={len(pairs_by_signature)}; "
        f"collision_groups={len(collision_groups)}; "
        f"largest_collision_group={largest_collision}; "
        "order_enriched_signatures="
        f"{len(pairs_by_ordered_signature)}; "
        "order_enriched_collision_groups="
        f"{len(ordered_collision_groups)}; "
        "largest_order_enriched_collision_group="
        f"{largest_ordered_collision}"
    )
    print(f"  critical_pair_digest={pair_digest}")
    print(
        "  peak_rss_mb="
        f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024:.1f}; "
        f"runtime_seconds={time.perf_counter() - started:.3f}"
    )
    if output_path is not None:
        print(f"  wrote={output_path}")


def rooted_activation_antimatroid_profile(
    input_path: str,
    output_path: str | None,
) -> None:
    """Profile full common-availability antimatroids of critical pairs."""

    started = time.perf_counter()
    source = json.loads(Path(input_path).read_text(encoding="utf-8"))
    if not isinstance(source, dict):
        raise ValueError("the rooted repair certificate must be a JSON object")
    if source.get("kind") != "rooted-critical-repair-profile":
        raise ValueError("the input is not a rooted critical-pair certificate")
    dimension = int(source["old_dimension"])
    if dimension != 4:
        raise ValueError(
            "the activation-antimatroid profile currently expects dimension four"
        )
    raw_pair_types = source.get("pair_types")
    if not isinstance(raw_pair_types, list):
        raise ValueError("the rooted repair certificate has no pair_types array")

    pairs: list[tuple[int, int]] = []
    pair_rows: list[dict[str, object]] = []
    signature_by_pair: dict[tuple[int, int], object] = {}
    full_mask = (1 << dimension) - 1
    for expected_index, item in enumerate(raw_pair_types):
        if not isinstance(item, dict) or item.get("index") != expected_index:
            raise ValueError("pair type indices must be consecutive from zero")
        base = int(str(item["base_mask"]), 16)
        site = int(str(item["site_mask"]), 16)
        pair = (base, site)
        if canonical_reverse_pair(base, site, dimension) != pair:
            raise AssertionError("a source pair is not canonically labelled")
        if not (
            is_critical_reverse_pair(base, site, dimension, 0)
            and is_critical_reverse_pair_by_closure(
                base,
                site,
                dimension,
                0,
            )
        ):
            raise AssertionError(
                "a certified pair failed recursive or closure criticality"
            )

        signature = canonical_activation_antimatroid_signature(
            base,
            site,
            dimension,
        )
        initial_feasible, forced_rows = signature
        if initial_feasible != (0,):
            raise AssertionError(
                "a critical pair has a nonempty initially feasible set"
            )
        for coordinate, row in enumerate(forced_rows):
            if full_mask not in row[0]:
                raise AssertionError(
                    "a forced contraction antimatroid is not full rank"
                )
            surviving_mask = full_mask ^ (1 << coordinate)
            if surviving_mask not in row[1]:
                raise AssertionError(
                    "a forced root-restriction antimatroid is not full rank"
                )

        pairs.append(pair)
        signature_by_pair[pair] = signature
        pair_rows.append(
            {
                "index": expected_index,
                "base_mask": format(base, "x"),
                "site_mask": format(site, "x"),
                "base_order": base.bit_count(),
                "site_order": site.bit_count(),
                "canonical_activation_antimatroid": (
                    activation_antimatroid_signature_payload(signature)
                ),
            }
        )

    pair_digest = hashlib.sha256(
        (
            "\n".join(f"{base:x}:{site:x}" for base, site in pairs)
            + "\n"
        ).encode("ascii")
    ).hexdigest()
    if pair_digest != source.get("critical_pair_digest"):
        raise AssertionError("the source critical-pair digest changed")

    pairs_by_signature: dict[object, list[tuple[int, int]]] = {}
    pairs_by_ordered_signature: dict[
        object,
        list[tuple[int, int]],
    ] = {}
    for pair in pairs:
        signature = signature_by_pair[pair]
        pairs_by_signature.setdefault(signature, []).append(pair)
        pairs_by_ordered_signature.setdefault(
            (
                pair[0].bit_count(),
                pair[1].bit_count(),
                signature,
            ),
            [],
        ).append(pair)

    collisions = [
        group
        for group in pairs_by_signature.values()
        if len(group) > 1
    ]
    ordered_collisions = [
        group
        for group in pairs_by_ordered_signature.values()
        if len(group) > 1
    ]
    pair_index = {pair: index for index, pair in enumerate(pairs)}
    signature_payloads = [
        activation_antimatroid_signature_payload(
            signature_by_pair[pair]
        )
        for pair in pairs
    ]
    signature_digest = hashlib.sha256(
        (
            "\n".join(
                json.dumps(
                    payload,
                    sort_keys=True,
                    separators=(",", ":"),
                )
                for payload in signature_payloads
            )
            + "\n"
        ).encode("utf-8")
    ).hexdigest()

    observed = (
        len(pairs),
        len(pairs_by_signature),
        len(collisions),
        max((len(group) for group in collisions), default=0),
        len(pairs_by_ordered_signature),
        len(ordered_collisions),
    )
    expected = (16, 8, 5, 5, 16, 0)
    if pair_digest == (
        "8b027112c96da224be5edb59f3e34529b85fef12a3b81b3ab9cdae6b21d88e4f"
    ) and observed != expected:
        raise AssertionError(
            "dimension-five activation-antimatroid regression changed: "
            f"expected={expected}, observed={observed}"
        )
    expected_order_21 = (37, 25, 6, 6, 37, 0)
    if pair_digest == (
        "f6471f050a0cd552ab259603dccea09c013b335d9aaa411e6891e16c997bbe0e"
    ) and observed != expected_order_21:
        raise AssertionError(
            "dimension-five order-21 activation-antimatroid regression "
            f"changed: expected={expected_order_21}, observed={observed}"
        )

    payload = {
        "schema_version": 1,
        "kind": "rooted-critical-activation-antimatroid-profile",
        "source_certificate": input_path,
        "source_pair_digest": pair_digest,
        "dimension": int(source["dimension"]),
        "old_dimension": dimension,
        "max_factor_order": int(source["max_factor_order"]),
        "critical_pair_types": len(pairs),
        "activation_antimatroid_signature_types": len(
            pairs_by_signature
        ),
        "activation_antimatroid_collision_groups": len(collisions),
        "largest_activation_antimatroid_collision_group": max(
            (len(group) for group in collisions),
            default=0,
        ),
        "order_enriched_activation_antimatroid_signature_types": len(
            pairs_by_ordered_signature
        ),
        "order_enriched_activation_antimatroid_collision_groups": len(
            ordered_collisions
        ),
        "activation_antimatroid_signature_digest": signature_digest,
        "pair_types": pair_rows,
        "collision_classes": [
            {
                "pair_type_indices": [
                    pair_index[pair] for pair in group
                ],
                "canonical_activation_antimatroid": (
                    activation_antimatroid_signature_payload(
                        signature_by_pair[group[0]]
                    )
                ),
            }
            for group in sorted(
                collisions,
                key=lambda members: (
                    -len(members),
                    tuple(pair_index[pair] for pair in members),
                ),
            )
        ],
        "order_enriched_collision_classes": [
            {
                "base_order": group[0][0].bit_count(),
                "site_order": group[0][1].bit_count(),
                "pair_type_indices": [
                    pair_index[pair] for pair in group
                ],
                "canonical_activation_antimatroid": (
                    activation_antimatroid_signature_payload(
                        signature_by_pair[group[0]]
                    )
                ),
            }
            for group in sorted(
                ordered_collisions,
                key=lambda members: (
                    -len(members),
                    tuple(pair_index[pair] for pair in members),
                ),
            )
        ],
    }
    if output_path is not None:
        Path(output_path).write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

    print("\nExact dimension-five activation-antimatroid profile")
    print(
        f"  critical_pair_types={len(pairs)}; "
        "activation_antimatroid_signatures="
        f"{len(pairs_by_signature)}; "
        "collision_groups="
        f"{len(collisions)}; "
        "order_enriched_activation_antimatroid_signatures="
        f"{len(pairs_by_ordered_signature)}; "
        "order_enriched_collision_groups="
        f"{len(ordered_collisions)}"
    )
    print(f"  source_pair_digest={pair_digest}")
    print(f"  activation_antimatroid_signature_digest={signature_digest}")
    print(
        "  peak_rss_mb="
        f"{resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024:.1f}; "
        f"runtime_seconds={time.perf_counter() - started:.3f}"
    )
    if output_path is not None:
        print(f"  wrote={output_path}")


def rooted_factor_dimension_five_census(order_limit: int) -> None:
    """Find all root-critical dimension-five factors through an order bound.

    Every member of ``R_per`` is ample and contractions preserve ampleness, so
    it is complete to expand the exact ample dimension-four bases.  All roots
    of every resulting cube-isomorphism type are then retained before pointed
    canonicalization.
    """

    if order_limit < 6 or order_limit > 21:
        raise ValueError("the dimension-five rooted census supports orders 6 through 21")

    started = time.perf_counter()
    dimension_four_bases = [
        family
        for family in enumerate_partial_cubes(4)[4]
        if family.bit_count() < order_limit
        and shattering_profile(family, 4)[2]
    ]
    raw = expansion_candidates(
        dimension_four_bases, 4, order_limit
    )
    r_types = {
        canonical_cube_family(family, 5)
        for family in raw
        if shattering_profile(family, 5)[2]
        and is_r_per(family, 5)
    }
    pointed = {
        canonical_pointed_family(family, 5, root)
        for family in r_types
        for root in vertices(family, 5)
    }
    critical = {
        family
        for family in pointed
        if is_root_critical(family, 5, 0)
    }

    print("\nExact dimension-five rooted-factor census")
    print(
        f"  order_limit={order_limit}; "
        f"d4_ample_bases={len(dimension_four_bases)}; "
        f"raw_expansions={len(raw)}; R_per_types={len(r_types)}; "
        f"pointed_R_per_types={len(pointed)}; "
        f"root_critical_types={len(critical)}"
    )
    factor = root_locked_factor()
    wedge_orders: Counter[int] = Counter()
    for family in sorted(
        critical, key=lambda item: (item.bit_count(), item)
    ):
        peripheries = peripheral_decompositions(family, 5)
        root_degree = sum(
            bool(family >> (1 << coordinate) & 1)
            for coordinate in range(5)
        )
        root_last_peelable = has_corner_peeling_ending_at_root(
            family,
            5,
            0,
        )
        root_last_order = corner_peeling_order_ending_at_root(
            family,
            5,
            0,
        )
        if root_last_peelable != (root_last_order is not None):
            raise AssertionError(
                "projected and original-coordinate root-last "
                "corner-peeling tests disagree"
            )
        shattered_histogram = shattered_rank_histogram(
            family,
            5,
        )
        locked_roots = tuple(
            root
            for root in vertices(family, 5)
            if is_root_locked(family, 5, root)
        )
        critical_roots = tuple(
            root
            for root in locked_roots
            if is_root_critical(family, 5, root)
        )
        wedge, wedge_dimension = one_vertex_wedge(
            factor, 4, 0, family, 5, 0
        )
        wedge_report = verify_obstruction(wedge, wedge_dimension)
        if not (
            wedge_report["periphery_free"]
            and not wedge_report["R_per"]
            and wedge_report["all_immediate_minors_R_per"]
        ):
            raise AssertionError("a root-critical factor produced a nonminimal wedge")
        wedge_orders[wedge_report["order"]] += 1
        print(
            "  critical: "
            f"order={family.bit_count()}, root_degree={root_degree}, "
            "root_last_corner_peelable="
            f"{root_last_peelable}, "
            "root_last_deletion_order="
            f"{tuple(format(vertex, '05b') for vertex in root_last_order or ())}, "
            f"shattered_rank_histogram={shattered_histogram}, "
            "root_locked_vertices="
            f"{tuple(format(root, '05b') for root in locked_roots)}, "
            "root_critical_vertices="
            f"{tuple(format(root, '05b') for root in critical_roots)}, "
            f"peripheral_coordinates={len(peripheries)}, "
            "peripheral_site_orders="
            f"{tuple(site.bit_count() for _, _, site, _ in peripheries)}, "
            f"vertices={bitstrings(family, 5)}"
        )
    print(
        "  wedges with L4: "
        f"verified={len(critical)}, "
        f"order_histogram={dict(sorted(wedge_orders.items()))}"
    )
    print(f"  runtime_seconds={time.perf_counter() - started:.3f}")


def run(trials: int, seed: int, dimension_five_limit: int | None) -> None:
    if len(HALL_CONCEPTS) != 299 or len(set(HALL_CONCEPTS)) != 299:
        raise AssertionError("the embedded Hall class must have 299 concepts")
    hall = family_from_vertices(HALL_CONCEPTS)
    hall_dimension = 12
    rng = random.Random(seed)

    print("=== Smallest additional peripheral-obstruction search ===")
    hall_shattered, hall_strong, hall_ample = shattering_profile(
        hall, hall_dimension
    )
    print(
        "Hall class: "
        f"dimension={hall_dimension}, order={hall.bit_count()}, "
        f"partial_cube={is_partial_cube_family(hall, hall_dimension)}, "
        f"shattered={hall_shattered}, strongly_shattered={hall_strong}, "
        f"ample={hall_ample}, R_per={is_r_per(hall, hall_dimension)}"
    )
    print(
        "Hall shattered-size histogram: "
        f"{shattered_size_histogram(hall, hall_dimension)}"
    )

    searches = ["small-order", "large-order", "few-peripheral"] + [
        "random"
    ] * trials
    witnesses: dict[tuple[int, int], tuple[str, list[tuple[int, int, str]]]] = {}
    for index, strategy in enumerate(searches, start=1):
        witness, dimension, path = greedy_minimal_obstruction(
            hall, hall_dimension, strategy, rng
        )
        witnesses.setdefault((witness, dimension), (strategy, path))
        print(
            f"descent {index:>3}/{len(searches)} ({strategy}): "
            f"terminal dimension={dimension}, order={witness.bit_count()}, "
            f"steps={len(path) - 1}"
        )

    ranked = sorted(
        witnesses.items(),
        key=lambda item: (item[0][1], item[0][0].bit_count(), item[0][0]),
    )
    print(f"\nDistinct labelled terminal witnesses: {len(ranked)}")
    for number, ((family, dimension), (strategy, path)) in enumerate(
        ranked, start=1
    ):
        report = verify_obstruction(family, dimension)
        print(f"\nWitness {number} (first found by {strategy})")
        print("  " + ", ".join(f"{key}={value}" for key, value in report.items()))
        print("  descent path (operation, dimension, order):")
        print(
            "    "
            + " -> ".join(
                f"{operation}:d{step_dimension}/n{step_family.bit_count()}"
                for step_family, step_dimension, operation in path
            )
        )
        print(f"  vertices={bitstrings(family, dimension)}")
        assert report["partial_cube"]
        assert report["periphery_free"]
        assert not report["R_per"]
        assert report["all_immediate_minors_R_per"]
        assert report["ample"]
        assert not report["is_Q_n_minus_minus"]

    pentagon = pentagonal_maximum_class()
    print("\nPentagonal 16-vertex obstruction")
    print(
        "  gated profile: "
        + ", ".join(
            f"{key}={value}"
            for key, value in gated_prime_profile(pentagon, 5).items()
        )
    )

    benchmark = hall_22_obstruction()
    print("\nHall-derived 22-vertex benchmark")
    print(
        "  "
        + ", ".join(
            f"{key}={value}"
            for key, value in verify_obstruction(benchmark, 6).items()
        )
    )
    print(
        "  structural profile: "
        + ", ".join(
            f"{key}={value}"
            for key, value in structural_profile(benchmark, 6).items()
        )
    )
    print(
        "  gated profile: "
        + ", ".join(
            f"{key}={value}"
            for key, value in gated_prime_profile(benchmark, 6).items()
        )
    )

    factor = root_locked_factor()
    wedge = rooted_wedge_obstruction()
    print("\nRooted-wedge 21-vertex obstruction")
    print(
        "  factor: "
        f"vertices={bitstrings(factor, 4)}, "
        f"R_per={is_r_per(factor, 4)}, "
        f"root_avoiding_at_zero={is_rooted_r_per(factor, 4, 0)}"
    )
    rooted_minor_rows = []
    for coordinate in range(4):
        zero = section(factor, 4, coordinate, 0)
        one = section(factor, 4, coordinate, 1)
        contraction = zero | one
        rooted_minor_rows.append(
            (
                coordinate,
                is_rooted_r_per(zero, 3, 0),
                is_rooted_r_per(contraction, 3, 0),
                is_r_per(one, 3),
            )
        )
    print(f"  rooted factor-minor rows={rooted_minor_rows}")
    print(
        "  obstruction: "
        + ", ".join(
            f"{key}={value}"
            for key, value in verify_obstruction(wedge, 8).items()
        )
    )
    print(
        "  gated profile: "
        + ", ".join(
            f"{key}={value}"
            for key, value in gated_prime_profile(wedge, 8).items()
        )
    )

    cache = is_r_per.cache_info()
    print(
        "\nR_per cache: "
        f"hits={cache.hits}, misses={cache.misses}, size={cache.currsize}"
    )
    if dimension_five_limit is not None:
        search_dimension_five(dimension_five_limit)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--random-trials",
        type=int,
        default=8,
        help="number of additional randomized Hall-minor descents",
    )
    parser.add_argument(
        "--skip-main-search",
        action="store_true",
        help="skip the Hall, pentagonal, and rooted-wedge benchmark report",
    )
    parser.add_argument(
        "--global-order-limit",
        type=int,
        default=None,
        help="exhaust ample expansions in every dimension through this order",
    )
    parser.add_argument("--seed", type=int, default=6601)
    parser.add_argument(
        "--dimension-five-limit",
        type=int,
        default=None,
        help="exhaust dimension-five expansions through this vertex order",
    )
    parser.add_argument(
        "--ample-walk-dimension",
        type=int,
        default=None,
        help="run random ample deletion walks in this cube dimension",
    )
    parser.add_argument(
        "--ample-walk-trials",
        type=int,
        default=0,
        help="number of random ample deletion walks",
    )
    parser.add_argument(
        "--ample-walk-report-below",
        type=int,
        default=24,
        help="record random-walk states at or below this order",
    )
    parser.add_argument(
        "--sparse-d6-trials",
        type=int,
        default=0,
        help="sample this many sparse dimension-six expansions",
    )
    parser.add_argument(
        "--exact-d6-order-limit",
        type=int,
        default=None,
        help="exhaust ample dimension-six expansions through this order",
    )
    parser.add_argument(
        "--exact-d6-base-start",
        type=int,
        default=0,
        help="zero-based first canonical dimension-five base in the exact census",
    )
    parser.add_argument(
        "--exact-d6-base-stop",
        type=int,
        default=None,
        help="exclusive canonical-base stop index for a bounded exact census chunk",
    )
    parser.add_argument(
        "--exact-d6-canonicalize-ample",
        action="store_true",
        help="also cube-canonicalize all ample dimension-six outputs",
    )
    parser.add_argument(
        "--ample-d6-base-order-limit",
        type=int,
        default=None,
        help="enumerate all ample dimension-six contraction bases through this order",
    )
    parser.add_argument(
        "--ample-d6-output",
        type=str,
        default=None,
        help="write the exact dimension-six base census to this JSON artifact",
    )
    parser.add_argument(
        "--exact-d7-order-limit",
        type=int,
        default=None,
        help="exhaust ample dimension-seven expansions through this order",
    )
    parser.add_argument(
        "--exact-d7-base-start",
        type=int,
        default=0,
        help="zero-based first dimension-six base in the exact dimension-seven census",
    )
    parser.add_argument(
        "--exact-d7-base-stop",
        type=int,
        default=None,
        help="exclusive dimension-six-base stop index for an exact dimension-seven chunk",
    )
    parser.add_argument(
        "--ample-type-dimension",
        type=int,
        default=None,
        help="enumerate all ample types in this exact dimension",
    )
    parser.add_argument(
        "--ample-type-order-limit",
        type=int,
        default=15,
        help="vertex-order bound for the generic ample-type census",
    )
    parser.add_argument(
        "--ample-type-source",
        type=str,
        default=None,
        help="read exact preceding-dimension types from this census artifact",
    )
    parser.add_argument(
        "--ample-type-base-start",
        type=int,
        default=0,
        help="first eligible source-base index for a resumable ample-type chunk",
    )
    parser.add_argument(
        "--ample-type-base-stop",
        type=int,
        default=None,
        help="exclusive eligible source-base index for a resumable ample-type chunk",
    )
    parser.add_argument(
        "--ample-type-output",
        type=str,
        default=None,
        help="write a resumable ample-type chunk to this JSON artifact",
    )
    parser.add_argument(
        "--merge-ample-type-chunks",
        nargs="+",
        default=None,
        help="merge these exact-cover ample-type chunk artifacts",
    )
    parser.add_argument(
        "--merge-ample-type-output",
        type=str,
        default=None,
        help="write the merged complete ample-type census to this artifact",
    )
    parser.add_argument(
        "--periphery-free-source",
        type=str,
        default=None,
        help="read exact ample contraction bases from this census artifact",
    )
    parser.add_argument(
        "--periphery-free-dimension",
        type=int,
        default=None,
        help="target dimension for an exact periphery-free search chunk",
    )
    parser.add_argument(
        "--periphery-free-order-limit",
        type=int,
        default=15,
        help="vertex-order bound for an exact periphery-free search chunk",
    )
    parser.add_argument(
        "--periphery-free-base-start",
        type=int,
        default=0,
        help="first source-base index for an exact periphery-free search chunk",
    )
    parser.add_argument(
        "--periphery-free-base-stop",
        type=int,
        default=None,
        help="exclusive source-base index for an exact periphery-free search chunk",
    )
    parser.add_argument(
        "--periphery-free-output",
        type=str,
        default=None,
        help="write an exact periphery-free search chunk to this artifact",
    )
    parser.add_argument(
        "--merge-periphery-free-chunks",
        nargs="+",
        default=None,
        help="merge these exact-cover periphery-free search chunks",
    )
    parser.add_argument(
        "--merge-periphery-free-output",
        type=str,
        default=None,
        help="write the merged complete periphery-free search artifact",
    )
    parser.add_argument(
        "--expansion-periphery-regression",
        action="store_true",
        help="exhaustively compare direct and lifted periphery predicates",
    )
    parser.add_argument(
        "--tree-canonicalizer-regression",
        action="store_true",
        help="exhaustively test the exact tree canonicalizer",
    )
    parser.add_argument(
        "--refined-canonicalizer-regression",
        action="store_true",
        help="exhaustively test the exact refined canonicalizer",
    )
    parser.add_argument(
        "--ample-predicate-regression",
        action="store_true",
        help="compare the Sandwich-count and full ampleness predicates",
    )
    parser.add_argument(
        "--rooted-factor-max-dimension",
        type=int,
        default=None,
        help="classify pointed R_per factors through this dimension",
    )
    parser.add_argument(
        "--rooted-wedge-test-max-dimension",
        type=int,
        default=8,
        help="test pointed-factor wedges through this total dimension",
    )
    parser.add_argument(
        "--rooted-factor-d5-order-limit",
        type=int,
        default=None,
        help="classify dimension-five root-critical factors through this order",
    )
    parser.add_argument(
        "--rooted-synchronization-max-dimension",
        type=int,
        default=None,
        help="verify rooted elimination and reverse-pair synchronization",
    )
    parser.add_argument(
        "--rooted-repair-d5-order-limit",
        type=int,
        default=None,
        help="profile first-step repair signatures of critical reverse pairs",
    )
    parser.add_argument(
        "--rooted-repair-output",
        type=str,
        default=None,
        help="write the exact rooted repair profile to this JSON artifact",
    )
    parser.add_argument(
        "--rooted-activation-profile-input",
        type=str,
        default=None,
        help="read exact critical pairs and profile full common antimatroids",
    )
    parser.add_argument(
        "--rooted-activation-output",
        type=str,
        default=None,
        help="write the rooted activation-antimatroid profile to JSON",
    )
    parser.add_argument(
        "--complete-rooted-pair-d4-output",
        type=str,
        default=None,
        help=(
            "exhaust all critical reverse pairs on four old coordinates "
            "and write their activation profile"
        ),
    )
    parser.add_argument(
        "--rooted-family-d5-inventory-output",
        type=str,
        default=None,
        help=(
            "generate every rooted-good family in Q5 and write a compact "
            "state-space inventory"
        ),
    )
    parser.add_argument(
        "--rooted-pair-d5-factor-order-limit",
        type=int,
        default=None,
        help=(
            "search critical reverse pairs on five old coordinates through "
            "this lift order"
        ),
    )
    parser.add_argument(
        "--rooted-pair-d5-output",
        type=str,
        default=None,
        help="write the exact bounded five-coordinate pair profile",
    )
    parser.add_argument(
        "--complete-site-three-pair-d5-output",
        type=str,
        default=None,
        help=(
            "exhaust the three-vertex-site critical-pair branch on five "
            "old coordinates"
        ),
    )
    parser.add_argument(
        "--complete-lower-q-obstruction-d4-output",
        type=str,
        default=None,
        help=(
            "directly exhaust marked lower q-activation obstructions "
            "through four coordinates"
        ),
    )
    parser.add_argument(
        "--lower-q-obstruction-site-three-source",
        type=str,
        default=None,
        help=(
            "cross-check the direct lower census against this complete "
            "Q5 site-three certificate"
        ),
    )
    parser.add_argument(
        "--lower-q-fiber-input",
        type=str,
        default=None,
        help=(
            "read the complete marked d=4 lower-obstruction certificate "
            "and profile its exact q-fiber witness systems"
        ),
    )
    parser.add_argument(
        "--lower-q-fiber-output",
        type=str,
        default=None,
        help="write the exact marked d=4 q-fiber witness profile",
    )
    parser.add_argument(
        "--lower-q-functor-input",
        type=str,
        default=None,
        help=(
            "read the complete marked d=4 lower-obstruction certificate "
            "and test the signed repair functor"
        ),
    )
    parser.add_argument(
        "--lower-q-functor-output",
        type=str,
        default=None,
        help="write the exact marked d=4 signed-repair functor profile",
    )
    parser.add_argument(
        "--principal-filter-falsifier-output",
        type=str,
        default=None,
        help=(
            "audit complete small defect-atlas fibers and write the "
            "canonical P5 principal-filter falsifier"
        ),
    )
    parser.add_argument(
        "--site-three-activation-input",
        type=str,
        default=None,
        help=(
            "read a complete Q5 site-three pair certificate and profile "
            "its activation invariants"
        ),
    )
    parser.add_argument(
        "--site-three-activation-output",
        type=str,
        default=None,
        help="write the complete Q5 site-three activation profile",
    )
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    if not arguments.skip_main_search:
        run(
            arguments.random_trials,
            arguments.seed,
            arguments.dimension_five_limit,
        )
    if arguments.global_order_limit is not None:
        search_all_dimensions_below_order(arguments.global_order_limit)
    if arguments.ample_walk_dimension is not None:
        random_ample_deletion_walks(
            arguments.ample_walk_dimension,
            arguments.ample_walk_trials,
            arguments.seed,
            arguments.ample_walk_report_below,
        )
    if arguments.sparse_d6_trials:
        random_sparse_dimension_six_expansions(
            arguments.sparse_d6_trials,
            arguments.seed,
        )
    if arguments.exact_d6_order_limit is not None:
        search_dimension_six_below_order(
            arguments.exact_d6_order_limit,
            arguments.exact_d6_base_start,
            arguments.exact_d6_base_stop,
            arguments.exact_d6_canonicalize_ample,
        )
    if arguments.ample_d6_base_order_limit is not None:
        dimension_six_types = ample_dimension_six_base_census(
            arguments.ample_d6_base_order_limit
        )
        if arguments.ample_d6_output is not None:
            write_family_census_artifact(
                arguments.ample_d6_output,
                kind="ample-type-census",
                dimension=6,
                order_limit=arguments.ample_d6_base_order_limit,
                families=dimension_six_types,
                metadata={"complete": True},
            )
    if arguments.exact_d7_order_limit is not None:
        search_dimension_seven_below_order(
            arguments.exact_d7_order_limit,
            arguments.exact_d7_base_start,
            arguments.exact_d7_base_stop,
        )
    if arguments.ample_type_dimension is not None:
        if arguments.ample_type_source is None:
            ample_types = ample_type_census(
                arguments.ample_type_dimension,
                arguments.ample_type_order_limit,
            )
            if arguments.ample_type_output is not None:
                write_family_census_artifact(
                    arguments.ample_type_output,
                    kind="ample-type-census",
                    dimension=arguments.ample_type_dimension,
                    order_limit=arguments.ample_type_order_limit,
                    families=ample_types,
                    metadata={"complete": True},
                )
        else:
            if arguments.ample_type_output is None:
                raise ValueError(
                    "--ample-type-output is required with --ample-type-source"
                )
            ample_type_census_chunk(
                arguments.ample_type_source,
                arguments.ample_type_dimension,
                arguments.ample_type_order_limit,
                arguments.ample_type_base_start,
                arguments.ample_type_base_stop,
                arguments.ample_type_output,
            )
    if arguments.merge_ample_type_chunks is not None:
        if arguments.merge_ample_type_output is None:
            raise ValueError(
                "--merge-ample-type-output is required when merging chunks"
            )
        merge_ample_type_census_chunks(
            arguments.merge_ample_type_chunks,
            arguments.merge_ample_type_output,
        )
    if arguments.periphery_free_dimension is not None:
        if (
            arguments.periphery_free_source is None
            or arguments.periphery_free_output is None
        ):
            raise ValueError(
                "--periphery-free-source and --periphery-free-output "
                "are required for a search chunk"
            )
        periphery_free_search_chunk(
            arguments.periphery_free_source,
            arguments.periphery_free_dimension,
            arguments.periphery_free_order_limit,
            arguments.periphery_free_base_start,
            arguments.periphery_free_base_stop,
            arguments.periphery_free_output,
        )
    if arguments.merge_periphery_free_chunks is not None:
        if arguments.merge_periphery_free_output is None:
            raise ValueError(
                "--merge-periphery-free-output is required when merging chunks"
            )
        merge_periphery_free_search_chunks(
            arguments.merge_periphery_free_chunks,
            arguments.merge_periphery_free_output,
        )
    if arguments.expansion_periphery_regression:
        expansion_periphery_predicate_regression()
    if arguments.tree_canonicalizer_regression:
        tree_canonicalizer_regression()
    if arguments.refined_canonicalizer_regression:
        refined_canonicalizer_regression()
    if arguments.ample_predicate_regression:
        ample_predicate_regression()
    if arguments.rooted_factor_max_dimension is not None:
        rooted_factor_census(
            arguments.rooted_factor_max_dimension,
            arguments.rooted_wedge_test_max_dimension,
        )
    if arguments.rooted_factor_d5_order_limit is not None:
        rooted_factor_dimension_five_census(
            arguments.rooted_factor_d5_order_limit
        )
    if arguments.rooted_synchronization_max_dimension is not None:
        rooted_synchronization_regression(
            arguments.rooted_synchronization_max_dimension
        )
    if arguments.rooted_repair_d5_order_limit is not None:
        rooted_repair_profile_dimension_five(
            arguments.rooted_repair_d5_order_limit,
            arguments.rooted_repair_output,
        )
    if arguments.rooted_activation_profile_input is not None:
        rooted_activation_antimatroid_profile(
            arguments.rooted_activation_profile_input,
            arguments.rooted_activation_output,
        )
    if arguments.complete_rooted_pair_d4_output is not None:
        complete_rooted_critical_pair_profile_dimension_four(
            arguments.complete_rooted_pair_d4_output,
        )
    if arguments.rooted_family_d5_inventory_output is not None:
        rooted_family_dimension_five_inventory(
            arguments.rooted_family_d5_inventory_output,
        )
    if arguments.rooted_pair_d5_factor_order_limit is not None:
        if arguments.rooted_pair_d5_output is None:
            raise ValueError(
                "--rooted-pair-d5-output is required with "
                "--rooted-pair-d5-factor-order-limit"
            )
        rooted_critical_pair_profile_dimension_five(
            arguments.rooted_pair_d5_factor_order_limit,
            arguments.rooted_pair_d5_output,
        )
    if arguments.complete_site_three_pair_d5_output is not None:
        complete_site_three_critical_pair_profile_dimension_five(
            arguments.complete_site_three_pair_d5_output,
        )
    if arguments.complete_lower_q_obstruction_d4_output is not None:
        complete_lower_q_activation_obstruction_profile_dimension_four(
            arguments.complete_lower_q_obstruction_d4_output,
            arguments.lower_q_obstruction_site_three_source,
        )
    if arguments.lower_q_fiber_input is not None:
        if arguments.lower_q_fiber_output is None:
            raise ValueError(
                "--lower-q-fiber-output is required with "
                "--lower-q-fiber-input"
            )
        lower_q_fiber_witness_profile_dimension_four(
            arguments.lower_q_fiber_input,
            arguments.lower_q_fiber_output,
        )
    if arguments.lower_q_functor_input is not None:
        if arguments.lower_q_functor_output is None:
            raise ValueError(
                "--lower-q-functor-output is required with "
                "--lower-q-functor-input"
            )
        lower_q_signed_repair_functor_profile_dimension_four(
            arguments.lower_q_functor_input,
            arguments.lower_q_functor_output,
        )
    if arguments.principal_filter_falsifier_output is not None:
        principal_filter_falsifier_profile(
            arguments.principal_filter_falsifier_output,
        )
    if arguments.site_three_activation_input is not None:
        if arguments.site_three_activation_output is None:
            raise ValueError(
                "--site-three-activation-output is required with "
                "--site-three-activation-input"
            )
        site_three_activation_antimatroid_profile(
            arguments.site_three_activation_input,
            arguments.site_three_activation_output,
        )
