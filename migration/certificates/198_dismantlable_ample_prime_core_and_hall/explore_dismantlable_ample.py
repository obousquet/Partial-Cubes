"""Audit corner-peelable (dismantlable) ample classes.

For an ample family, a vertex is a corner precisely when it belongs to a
unique maximal cube.  Equivalently, if ``D(v)`` is the set of coordinates
on edges incident with ``v``, then every vertex obtained from ``v`` by
toggling an arbitrary subset of ``D(v)`` belongs to the family.

The script uses this local test to:

* verify Hall's 299-concept maximum class has no corner;
* inspect all of its elementary coordinate sections and contractions; and
* scan the existing exact ample-family censuses through order 14.

The finite output is evidence about small obstructions.  Pc-minor closure of
the class is proved symbolically in the survey and does not depend on this
audit.

Run conservatively with

    nice -n 10 uv run python explore_dismantlable_ample.py
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from functools import lru_cache
from pathlib import Path

from explore_smallest_additional_peripheral_obstruction import (
    HALL_CONCEPTS,
    ample_type_census,
    enumerate_partial_cubes,
    family_from_vertices,
    irredundant_projection,
    is_ample_family,
    project_vertex,
    section,
    vertices,
)


DEFAULT_CENSUS_PATTERN = (
    "computations/cover_recursive_peripheral/ample_d{dimension}_o14.json"
)


@lru_cache(maxsize=None)
def corners(family: int, dimension: int) -> tuple[int, ...]:
    """Return the CCMW corners of an ample family.

    The incident directions at ``v`` generate the smallest cube containing
    ``v`` and all its neighbors.  The vertex is in a unique maximal cube
    exactly when this whole generated cube is present.
    """

    members = set(vertices(family, dimension))
    answer: list[int] = []
    for vertex in members:
        directions = tuple(
            coordinate
            for coordinate in range(dimension)
            if (vertex ^ (1 << coordinate)) in members
        )
        generated_face = {
            vertex ^ sum(1 << directions[index] for index in range(len(directions))
                         if subset & (1 << index))
            for subset in range(1 << len(directions))
        }
        if generated_face <= members:
            answer.append(vertex)
    return tuple(sorted(answer))


@lru_cache(maxsize=None)
def is_dismantlable(family: int, dimension: int) -> bool:
    """Decide whether an ample family has a complete corner peeling."""

    family, dimension = irredundant_projection(family, dimension)
    if family.bit_count() == 1:
        return True
    return any(
        is_dismantlable(
            *irredundant_projection(family & ~(1 << corner), dimension)
        )
        for corner in corners(family, dimension)
    )


def contraction(family: int, dimension: int, coordinate: int) -> int:
    """Project away one coordinate, i.e. contract its Theta-class."""

    return family_from_vertices(
        project_vertex(vertex, coordinate)
        for vertex in vertices(family, dimension)
    )


def hall_audit() -> dict[str, object]:
    """Audit Hall's class and all elementary coordinate pc-minors."""

    hall = family_from_vertices(HALL_CONCEPTS)
    dimension = 12
    section_results: list[dict[str, object]] = []
    contraction_results: list[dict[str, object]] = []

    for coordinate in range(dimension):
        for value in (0, 1):
            minor, lower_dimension = irredundant_projection(
                section(hall, dimension, coordinate, value),
                dimension - 1,
            )
            section_results.append(
                {
                    "coordinate": coordinate,
                    "value": value,
                    "dimension": lower_dimension,
                    "order": minor.bit_count(),
                    "corners": len(corners(minor, lower_dimension)),
                    "dismantlable": is_dismantlable(minor, lower_dimension),
                }
            )

        minor, lower_dimension = irredundant_projection(
            contraction(hall, dimension, coordinate),
            dimension - 1,
        )
        contraction_results.append(
            {
                "coordinate": coordinate,
                "dimension": lower_dimension,
                "order": minor.bit_count(),
                "corners": len(corners(minor, lower_dimension)),
                "dismantlable": is_dismantlable(minor, lower_dimension),
            }
        )

    return {
        "dimension": dimension,
        "order": hall.bit_count(),
        "corners": len(corners(hall, dimension)),
        "dismantlable": is_dismantlable(hall, dimension),
        "sections": section_results,
        "contractions": contraction_results,
    }


def census_families(dimension: int) -> tuple[int, ...]:
    """Load or generate the exact ample types in one dimension."""

    if dimension <= 3:
        return tuple(
            family
            for family in enumerate_partial_cubes(3)[dimension]
            if is_ample_family(family, dimension)
        )
    if dimension <= 5:
        return tuple(sorted(ample_type_census(dimension, 14)))
    path = Path(DEFAULT_CENSUS_PATTERN.format(dimension=dimension))
    payload = json.loads(path.read_text())
    return tuple(int(value, 16) for value in payload["families_hex"])


def census_audit(first_dimension: int, last_dimension: int) -> dict[str, object]:
    """Scan the repo's exact ample-family censuses."""

    rows: list[dict[str, object]] = []
    total = 0
    failures: list[dict[str, object]] = []
    for dimension in range(first_dimension, last_dimension + 1):
        families = census_families(dimension)
        order_histogram = Counter(family.bit_count() for family in families)
        nondismantlable: list[int] = []
        for family in families:
            if not is_dismantlable(family, dimension):
                nondismantlable.append(family)
        total += len(families)
        row = {
            "dimension": dimension,
            "families": len(families),
            "order_histogram": dict(sorted(order_histogram.items())),
            "nondismantlable": len(nondismantlable),
        }
        rows.append(row)
        for family in nondismantlable:
            failures.append(
                {
                    "dimension": dimension,
                    "order": family.bit_count(),
                    "family_hex": format(family, "x"),
                }
            )
        print(
            f"dimension={dimension}: families={len(families)}; "
            f"nondismantlable={len(nondismantlable)}",
            flush=True,
        )

    return {
        "first_dimension": first_dimension,
        "last_dimension": last_dimension,
        "order_limit": 14,
        "families": total,
        "rows": rows,
        "nondismantlable": failures,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--first-dimension", type=int, default=0)
    parser.add_argument("--last-dimension", type=int, default=11)
    parser.add_argument(
        "--skip-census",
        action="store_true",
        help="audit only Hall's class and its elementary minors",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "computations/dismantlable_ample/"
            "dismantlable_ample_audit.json"
        ),
    )
    args = parser.parse_args()

    hall = hall_audit()
    print(
        "Hall class: "
        f"order={hall['order']}; corners={hall['corners']}; "
        f"dismantlable={hall['dismantlable']}",
        flush=True,
    )
    print(
        "Hall elementary minors: "
        f"sections_dismantlable="
        f"{sum(item['dismantlable'] for item in hall['sections'])}/"
        f"{len(hall['sections'])}; "
        f"contractions_dismantlable="
        f"{sum(item['dismantlable'] for item in hall['contractions'])}/"
        f"{len(hall['contractions'])}",
        flush=True,
    )

    census = None
    if not args.skip_census:
        census = census_audit(
            args.first_dimension,
            args.last_dimension,
        )

    payload = {
        "schema": "dismantlable-ample-audit-v1",
        "hall": hall,
        "census": census,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
