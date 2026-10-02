#!/usr/bin/env python3
"""Verify the local degree-three two-shore transversal lemma.

For two disjoint nonempty ample classes L,U in Q_3 that both avoid zero,
the lemma asks for disjoint inclusion-minimal transversals of the
inclusion-minimal words of L and U.  Their union leaves at most one of the
three coordinates.  In the carrier proof this yields either an isolated
relative support or a minimal--maximal support cover.
"""

from __future__ import annotations

import hashlib
import json


DIMENSION = 3
AMBIENT = tuple(range(1 << DIMENSION))


def shattered_supports(family: frozenset[int]) -> frozenset[int]:
    return frozenset(
        support
        for support in AMBIENT
        if len({word & support for word in family}) == 1 << support.bit_count()
    )


def is_ample(family: frozenset[int]) -> bool:
    return len(family) == len(shattered_supports(family))


def inclusion_minimal_words(family: frozenset[int]) -> frozenset[int]:
    return frozenset(
        word
        for word in family
        if not any(
            other != word and (other & word) == other for other in family
        )
    )


def minimal_transversals(clutter: frozenset[int]) -> tuple[int, ...]:
    transversals = []
    for candidate in range(1, 1 << DIMENSION):
        if not all(candidate & edge for edge in clutter):
            continue
        if any(
            all((candidate ^ (1 << coordinate)) & edge for edge in clutter)
            for coordinate in range(DIMENSION)
            if candidate & (1 << coordinate)
        ):
            continue
        transversals.append(candidate)
    return tuple(transversals)


def verify() -> dict[str, object]:
    classes = tuple(
        family
        for mask in range(1, 1 << len(AMBIENT))
        if not (mask & 1)
        for family in [
            frozenset(word for word in AMBIENT if mask & (1 << word))
        ]
        if is_ample(family)
    )

    pair_count = 0
    leftover_histogram = {0: 0, 1: 0}
    for lower in classes:
        for upper in classes:
            if lower & upper:
                continue
            pair_count += 1
            lower_blockers = minimal_transversals(
                inclusion_minimal_words(lower)
            )
            upper_blockers = minimal_transversals(
                inclusion_minimal_words(upper)
            )
            witnesses = [
                (first, second)
                for first in lower_blockers
                for second in upper_blockers
                if not (first & second)
            ]
            if not witnesses:
                raise AssertionError(
                    {
                        "lower": sorted(lower),
                        "upper": sorted(upper),
                        "lower_blockers": lower_blockers,
                        "upper_blockers": upper_blockers,
                    }
                )
            leftover = min(
                DIMENSION - (first | second).bit_count()
                for first, second in witnesses
            )
            if leftover not in leftover_histogram:
                raise AssertionError((lower, upper, witnesses, leftover))
            leftover_histogram[leftover] += 1

    record: dict[str, object] = {
        "status": "VERIFIED",
        "dimension": DIMENSION,
        "nonempty_zero_free_ample_classes": len(classes),
        "ordered_disjoint_pairs": pair_count,
        "minimum_leftover_coordinate_histogram": leftover_histogram,
        "claim": (
            "Every ordered disjoint pair has disjoint minimal transversals; "
            "their union leaves zero or one coordinate."
        ),
    }
    canonical = json.dumps(record, sort_keys=True, separators=(",", ":"))
    record["mathematical_digest"] = hashlib.sha256(canonical.encode()).hexdigest()
    return record


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
