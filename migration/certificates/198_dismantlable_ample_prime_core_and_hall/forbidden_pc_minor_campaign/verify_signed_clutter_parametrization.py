#!/usr/bin/env python3
"""Bounded regressions for the signed-clutter theorems in the paper.

This verifier is deliberately not an enumeration route.  It checks the
proved symbolic formulas against all ample cube families through dimension
four and against Hall's already-certified source neighborhood.  It does not
claim to classify higher-dimensional families or mutation components.
"""

from __future__ import annotations

import hashlib
import json
import sys
from collections import Counter, defaultdict
from itertools import combinations, permutations, product
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from explore_smallest_additional_peripheral_obstruction import (
    HALL_CONCEPTS,
    is_ample_family,
    shattered_coordinate_sets_mask,
)


HALL_SOURCE = (
    ROOT
    / "forbidden_pc_minor_campaign"
    / "computations"
    / "hall_corner_exchange"
    / "source_H.json"
)
OUTPUT = ROOT / "forbidden_pc_minor_campaign" / "signed_clutter_regression.json"


def members(family: int, dimension: int) -> tuple[int, ...]:
    return tuple(x for x in range(1 << dimension) if family & (1 << x))


def subsets(support: int):
    current = support
    while True:
        yield current
        if current == 0:
            return
        current = (current - 1) & support


def minimal_nonfaces(shadow: int, dimension: int) -> tuple[int, ...]:
    return tuple(
        support
        for support in range(1, 1 << dimension)
        if not shadow & (1 << support)
        and all(
            shadow & (1 << proper)
            for proper in subsets(support)
            if proper != support
        )
    )


def support_traces(support: int, dimension: int) -> tuple[int, ...]:
    return tuple(
        trace
        for trace in range(1 << dimension)
        if trace & ~support == 0
    )


def signing(
    family: int,
    shadow: int,
    dimension: int,
) -> dict[int, int]:
    points = members(family, dimension)
    result: dict[int, int] = {}
    for support in minimal_nonfaces(shadow, dimension):
        omitted = set(support_traces(support, dimension)) - {
            x & support for x in points
        }
        if len(omitted) != 1:
            raise AssertionError(
                f"support {support} has {len(omitted)} omitted traces"
            )
        result[support] = omitted.pop()
    return result


def satisfying_family(
    signs: dict[int, int],
    dimension: int,
) -> int:
    family = 0
    for x in range(1 << dimension):
        if all(x & support != trace for support, trace in signs.items()):
            family |= 1 << x
    return family


def first_loss_is_bijective(
    order: tuple[int, ...],
    shadow: int,
    dimension: int,
) -> bool:
    deaths: list[int] = []
    for support in range(1 << dimension):
        if not shadow & (1 << support):
            continue
        last = {
            trace: max(
                index
                for index, x in enumerate(order, start=1)
                if x & support == trace
            )
            for trace in support_traces(support, dimension)
        }
        deaths.append(min(last.values()))
    return sorted(deaths) == list(range(1, len(order) + 1))


def suffixes_are_ample(order: tuple[int, ...], dimension: int) -> bool:
    suffix = 0
    ample = True
    for x in reversed(order):
        suffix |= 1 << x
        ample &= is_ample_family(suffix, dimension)
    return ample


def low_dimension_profile() -> dict[str, object]:
    dimensions: dict[str, object] = {}
    for dimension in range(1, 5):
        ample_by_shadow: dict[int, list[int]] = defaultdict(list)
        for family in range(1, 1 << (1 << dimension)):
            if is_ample_family(family, dimension):
                shadow = shattered_coordinate_sets_mask(family, dimension)
                ample_by_shadow[shadow].append(family)

        signed_presentations = 0
        admissible_signings = 0
        exchange_tests = 0
        first_loss_orders = 0
        one_clause_pivots = 0
        wide_one_clause_pivots = 0

        for shadow, families in ample_by_shadow.items():
            shadows_size = shadow.bit_count()
            signs_by_family: dict[int, dict[int, int]] = {}
            for family in families:
                signs = signing(family, shadow, dimension)
                if satisfying_family(signs, dimension) != family:
                    raise AssertionError("signing does not recover ample family")
                signs_by_family[family] = signs

                points = members(family, dimension)
                if dimension <= 3:
                    orders = permutations(points)
                else:
                    orders = iter((points, tuple(reversed(points))))
                for order in orders:
                    first_loss_orders += 1
                    if first_loss_is_bijective(
                        order, shadow, dimension
                    ) != suffixes_are_ample(order, dimension):
                        raise AssertionError("trace-death criterion failed")

            nonfaces = minimal_nonfaces(shadow, dimension)
            domains = [support_traces(support, dimension) for support in nonfaces]
            for values in product(*domains):
                signed_presentations += 1
                signs = dict(zip(nonfaces, values))
                family = satisfying_family(signs, dimension)
                if family.bit_count() == shadows_size:
                    admissible_signings += 1
                    if not is_ample_family(family, dimension):
                        raise AssertionError("admissible signing is not ample")
                    if shattered_coordinate_sets_mask(
                        family, dimension
                    ) != shadow:
                        raise AssertionError("admissible signing has wrong shadow")

            for family in families:
                signs = signs_by_family[family]
                points = members(family, dimension)
                for added in range(1 << dimension):
                    if family & (1 << added):
                        continue
                    violated = {
                        support
                        for support, trace in signs.items()
                        if added & support == trace
                    }
                    extension = family | (1 << added)
                    if (len(violated) == 1) != is_ample_family(
                        extension, dimension
                    ):
                        raise AssertionError("one-clause extension law failed")

                    for removed in points:
                        private = {
                            support
                            for support in signs
                            if added & support != removed & support
                            and sum(
                                x & support == removed & support
                                for x in points
                            )
                            == 1
                        }
                        exchanged = extension & ~(1 << removed)
                        actual = is_ample_family(exchanged, dimension) and (
                            shattered_coordinate_sets_mask(
                                exchanged, dimension
                            )
                            == shadow
                        )
                        exchange_tests += 1
                        if (violated == private) != actual:
                            raise AssertionError("signed exchange law failed")

            for first, second in combinations(families, 2):
                changed_clauses = sum(
                    signs_by_family[first][support]
                    != signs_by_family[second][support]
                    for support in signs_by_family[first]
                )
                if changed_clauses == 1:
                    one_clause_pivots += 1
                    if (first ^ second).bit_count() > 2:
                        wide_one_clause_pivots += 1

        dimensions[str(dimension)] = {
            "ample_families": sum(map(len, ample_by_shadow.values())),
            "shattered_complexes": len(ample_by_shadow),
            "signed_presentations": signed_presentations,
            "admissible_signings": admissible_signings,
            "exchange_tests": exchange_tests,
            "first_loss_orders": first_loss_orders,
            "one_clause_pivots": one_clause_pivots,
            "wide_one_clause_pivots": wide_one_clause_pivots,
        }
    return dimensions


def hall_profile() -> dict[str, object]:
    payload = json.loads(HALL_SOURCE.read_text(encoding="utf-8"))
    recorded_extensions = {row["added"] for row in payload["extensions"]}
    hall = set(HALL_CONCEPTS)
    dimension = 12
    support_size = 4
    supports = tuple(
        sum(1 << coordinate for coordinate in subset)
        for subset in combinations(range(dimension), support_size)
    )
    signs = {
        support: (
            set(support_traces(support, dimension))
            - {x & support for x in hall}
        ).pop()
        for support in supports
    }
    one_clause_violators = {
        x
        for x in range(1 << dimension)
        if x not in hall
        and sum(x & support == trace for support, trace in signs.items()) == 1
    }
    if one_clause_violators != recorded_extensions:
        raise AssertionError("Hall one-clause violators differ from extensions")
    return {
        "dimension": dimension,
        "order": len(hall),
        "minimal_nonfaces": len(supports),
        "outside_concepts": (1 << dimension) - len(hall),
        "one_clause_violators": len(one_clause_violators),
        "recorded_ample_extensions": len(recorded_extensions),
        "source_profile_sha256": payload["mathematical_profile_sha256"],
    }


def main() -> int:
    mathematics = {
        "schema": "damp-signed-clutter-regression-v1",
        "scope": (
            "all ample families through dimension four plus Hall's "
            "already-certified source neighborhood; regression only"
        ),
        "low_dimensions": low_dimension_profile(),
        "hall": hall_profile(),
    }
    digest = hashlib.sha256(
        json.dumps(mathematics, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    OUTPUT.write_text(
        json.dumps(
            {**mathematics, "mathematical_profile_sha256": digest},
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"wrote {OUTPUT}")
    print(f"mathematical_profile_sha256={digest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
