#!/usr/bin/env python3
"""Audit how local one-hole certificates pass through Hall contractions.

For a restriction, local minimal nonfaces are obtained simply by discarding
the root supports that contain the restricted coordinate.  A contraction is
different: two side traces can jointly cover a punctured cube although
neither side contains that punctured cube.  Such a *switched hole* is
invisible if one only transports the root one-hole pieces independently.

This script classifies every canonical one-hole certificate in every
coordinate contraction of Hall's 299-concept obstruction as direct or
switched.  For switched certificates it records the two-colour punctured-cube
profile, verifies that the double-fibre set separates the two exclusive
regions, and checks that every component outside the double-fibre set has a
constant side colour.

The output is a bounded structural diagnostic, not a peelability certificate.
Run nontrivial instances on the worker.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import sys
from collections import Counter, deque
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from explore_smallest_additional_peripheral_obstruction import HALL_CONCEPTS


OUTPUT = (
    ROOT
    / "forbidden_pc_minor_campaign"
    / "damp_hall_contraction_hole_transport.json"
)


def digest(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def suppress(vertex: int, coordinate: int) -> int:
    lower = (1 << coordinate) - 1
    return (vertex & lower) | ((vertex >> (coordinate + 1)) << coordinate)


def insert_bit(vertex: int, coordinate: int, value: int) -> int:
    lower = (1 << coordinate) - 1
    return (
        (vertex & lower)
        | (value << coordinate)
        | ((vertex & ~lower) << 1)
    )


def lift_support(support: int, coordinate: int) -> int:
    return insert_bit(support, coordinate, 0)


def submasks(mask: int) -> tuple[int, ...]:
    answer = []
    submask = mask
    while True:
        answer.append(submask)
        if submask == 0:
            return tuple(answer)
        submask = (submask - 1) & mask


def local_minimal_nonfaces(
    family: frozenset[int], vertex: int, dimension: int
) -> tuple[int, ...]:
    """Return the supports of the canonical local one-hole certificates."""

    answer = []
    # Hall's class and all its elementary images have VC dimension three.
    # A minimal nonface then has order at most four.
    for size in range(2, min(4, dimension) + 1):
        for coordinates in itertools.combinations(range(dimension), size):
            support = sum(1 << item for item in coordinates)
            if vertex ^ support in family:
                continue
            if all(
                vertex ^ proper in family
                for proper in submasks(support)
                if proper != support
            ):
                answer.append(support)
    return tuple(sorted(answer))


def all_local_holes(
    family: frozenset[int], dimension: int
) -> dict[int, tuple[int, ...]]:
    return {
        vertex: local_minimal_nonfaces(family, vertex, dimension)
        for vertex in family
    }


def section(
    family: frozenset[int], dimension: int, coordinate: int, value: int
) -> frozenset[int]:
    return frozenset(
        suppress(vertex, coordinate)
        for vertex in family
        if (vertex >> coordinate) & 1 == value
    )


def contraction(
    family: frozenset[int], dimension: int, coordinate: int
) -> frozenset[int]:
    return frozenset(suppress(vertex, coordinate) for vertex in family)


def components_outside_separator(
    support: int,
    double_fibres: frozenset[int],
) -> tuple[frozenset[int], ...]:
    """Components of the punctured support cube after removing double fibres."""

    tip = support
    remaining = set(submasks(support)) - {tip} - set(double_fibres)
    answer = []
    while remaining:
        root = min(remaining)
        component = {root}
        remaining.remove(root)
        pending = deque([root])
        while pending:
            word = pending.popleft()
            for coordinate in range(support.bit_length()):
                if not (support >> coordinate) & 1:
                    continue
                neighbor = word ^ (1 << coordinate)
                if neighbor in remaining:
                    remaining.remove(neighbor)
                    component.add(neighbor)
                    pending.append(neighbor)
        answer.append(frozenset(component))
    return tuple(answer)


def switched_profile(
    family: frozenset[int],
    projected_vertex: int,
    projected_support: int,
    coordinate: int,
) -> dict[str, object]:
    root_support = lift_support(projected_support, coordinate)
    side_sets: dict[int, tuple[int, ...]] = {}
    for projected_flip in submasks(projected_support):
        root_flip = lift_support(projected_flip, coordinate)
        sides = tuple(
            side
            for side in (0, 1)
            if insert_bit(projected_vertex, coordinate, side) ^ root_flip
            in family
        )
        side_sets[projected_flip] = sides

    if side_sets[projected_support]:
        raise AssertionError("the projected missing tip has a surviving lift")
    if any(
        not sides
        for flip, sides in side_sets.items()
        if flip != projected_support
    ):
        raise AssertionError("the projected punctured cube is incomplete")

    double_fibres = frozenset(
        flip for flip, sides in side_sets.items() if len(sides) == 2
    )
    exclusive = {
        side: frozenset(
            flip for flip, sides in side_sets.items() if sides == (side,)
        )
        for side in (0, 1)
    }
    crossing_edges = []
    for left in exclusive[0]:
        for bit in range(projected_support.bit_length()):
            if not (projected_support >> bit) & 1:
                continue
            right = left ^ (1 << bit)
            if right in exclusive[1]:
                crossing_edges.append((left, right))
    if crossing_edges:
        raise AssertionError(("exclusive regions touch", crossing_edges))

    components = components_outside_separator(projected_support, double_fibres)
    component_colours = []
    for component in components:
        colours = {
            side
            for side in (0, 1)
            if component & exclusive[side]
        }
        if len(colours) != 1:
            raise AssertionError(("non-monochromatic component", component, colours))
        component_colours.append(next(iter(colours)))

    return {
        "support_order": bin(projected_support).count("1"),
        "root_support": root_support,
        "double_fibre_order": len(double_fibres),
        "exclusive_orders": [len(exclusive[0]), len(exclusive[1])],
        "component_count": len(components),
        "component_orders": sorted(map(len, components)),
        "component_colour_counts": {
            str(side): component_colours.count(side) for side in (0, 1)
        },
        "antipode_side_count": len(side_sets[0]),
    }


def hall_audit() -> dict[str, object]:
    family = frozenset(HALL_CONCEPTS)
    dimension = 12
    root_holes = all_local_holes(family, dimension)
    if any(not holes for holes in root_holes.values()):
        raise AssertionError("Hall's root unexpectedly has a corner")

    restriction_rows = []
    contraction_rows = []
    all_switched_profiles = []
    for coordinate in range(dimension):
        for value in (0, 1):
            restricted = section(family, dimension, coordinate, value)
            restricted_holes = all_local_holes(restricted, dimension - 1)
            predicted = {
                suppress(vertex, coordinate): tuple(
                    sorted(suppress(support, coordinate) for support in holes)
                )
                for vertex, holes in root_holes.items()
                if (vertex >> coordinate) & 1 == value
                and any(not ((support >> coordinate) & 1) for support in holes)
            }
            predicted = {
                vertex: tuple(
                    sorted(
                        suppress(support, coordinate)
                        for support in root_holes[
                            insert_bit(vertex, coordinate, value)
                        ]
                        if not ((support >> coordinate) & 1)
                    )
                )
                for vertex in restricted
            }
            if predicted != restricted_holes:
                disagreements = [
                    {
                        "vertex": vertex,
                        "predicted": predicted.get(vertex),
                        "actual": restricted_holes.get(vertex),
                    }
                    for vertex in sorted(set(predicted) | set(restricted_holes))
                    if predicted.get(vertex) != restricted_holes.get(vertex)
                ]
                raise AssertionError(
                    (
                        "restriction transport failed",
                        coordinate,
                        value,
                        disagreements[:5],
                    )
                )
            restriction_rows.append(
                {
                    "coordinate": coordinate,
                    "value": value,
                    "order": len(restricted),
                    "canonical_hole_count": sum(map(len, restricted_holes.values())),
                    "corner_count": sum(not holes for holes in restricted_holes.values()),
                    "transport_exact": True,
                }
            )

        projected = contraction(family, dimension, coordinate)
        projected_holes = all_local_holes(projected, dimension - 1)
        direct_count = 0
        switched = []
        for projected_vertex, supports in projected_holes.items():
            lifts = tuple(
                insert_bit(projected_vertex, coordinate, side)
                for side in (0, 1)
                if insert_bit(projected_vertex, coordinate, side) in family
            )
            for support in supports:
                root_support = lift_support(support, coordinate)
                direct_sides = [
                    (lift >> coordinate) & 1
                    for lift in lifts
                    if root_support in root_holes[lift]
                ]
                if direct_sides:
                    direct_count += 1
                else:
                    profile = switched_profile(
                        family, projected_vertex, support, coordinate
                    )
                    switched.append(profile)
                    all_switched_profiles.append(profile)
        contraction_rows.append(
            {
                "coordinate": coordinate,
                "order": len(projected),
                "canonical_hole_count": sum(map(len, projected_holes.values())),
                "corner_count": sum(not holes for holes in projected_holes.values()),
                "direct_hole_count": direct_count,
                "switched_hole_count": len(switched),
                "switched_profile_counts": [
                    {"profile": list(profile), "multiplicity": multiplicity}
                    for profile, multiplicity in sorted(
                        Counter(
                            (
                                row["support_order"],
                                row["double_fibre_order"],
                                *row["exclusive_orders"],
                                row["component_count"],
                                *row["component_orders"],
                            )
                            for row in switched
                        ).items()
                    )
                ],
            }
        )

    return {
        "dimension": dimension,
        "order": len(family),
        "root_canonical_hole_count": sum(map(len, root_holes.values())),
        "restrictions": restriction_rows,
        "contractions": contraction_rows,
        "totals": {
            "restriction_holes": sum(
                row["canonical_hole_count"] for row in restriction_rows
            ),
            "contraction_holes": sum(
                row["canonical_hole_count"] for row in contraction_rows
            ),
            "direct_contraction_holes": sum(
                row["direct_hole_count"] for row in contraction_rows
            ),
            "switched_contraction_holes": len(all_switched_profiles),
        },
    }


def small_switched_adversary() -> dict[str, object]:
    """The four-vertex path over a punctured square is the first switch."""

    # Coordinates 0,1 form the projected square; coordinate 2 is contracted.
    family = frozenset({0b000, 0b001, 0b100, 0b110})
    projected = contraction(family, 3, 2)
    holes = all_local_holes(projected, 2)
    if holes[0] != (0b11,):
        raise AssertionError(holes)
    root_holes = all_local_holes(family, 3)
    if 0b011 in root_holes[0b000] or 0b011 in root_holes[0b100]:
        raise AssertionError(root_holes)
    profile = switched_profile(family, 0, 0b11, 2)
    return {
        "family": sorted(family),
        "projected_family": sorted(projected),
        "projected_hole": {"antipode": 0, "support": 3},
        "profile": profile,
    }


def main() -> None:
    payload = {
        "schema": "damp-hall-contraction-hole-transport-v1",
        "scope": (
            "Canonical one-hole transport through all restrictions and "
            "contractions of Hall's 299-concept obstruction."
        ),
        "smallest_switched_adversary": small_switched_adversary(),
        "hall": hall_audit(),
    }
    payload["sha256_without_digest"] = digest(payload)
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    totals = payload["hall"]["totals"]
    print(f"restriction_holes={totals['restriction_holes']}")
    print(f"contraction_holes={totals['contraction_holes']}")
    print(f"direct_contraction_holes={totals['direct_contraction_holes']}")
    print(f"switched_contraction_holes={totals['switched_contraction_holes']}")
    print(f"sha256_without_digest={payload['sha256_without_digest']}")
    print(f"wrote {OUTPUT}")


if __name__ == "__main__":
    main()
