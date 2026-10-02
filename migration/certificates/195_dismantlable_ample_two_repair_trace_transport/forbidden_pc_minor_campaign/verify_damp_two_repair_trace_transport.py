#!/usr/bin/env python3
"""Regression for the signed trace transport of a two-repair lift.

Let ``S`` be ample and let ``A=S+{a}``, ``B=S+{b}`` be ample
one-concept extensions.  If their unique new shattered supports ``P_a``
and ``P_b`` are distinct, the opposite-layer lift

    X = (B x {0}) union (A x {1})

has a completely explicit shattered complex, signed-minimal-nonface
presentation, and trace table.  Moreover ``S,a,b`` are recovered from
``(X,e)`` by intersecting its two sections.

This script checks those symbolic formulas on every eligible ordered repair
pair of every ample family through dimension four.  It is a bounded
regression for the proof, not a proof of the unbounded theorem.
"""

from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from pathlib import Path

from verify_signed_clutter_parametrization import (
    is_ample_family,
    members,
    minimal_nonfaces,
    satisfying_family,
    shattered_coordinate_sets_mask,
    signing,
    support_traces,
)


ROOT = Path(__file__).resolve().parent.parent
OUTPUT = (
    ROOT
    / "forbidden_pc_minor_campaign"
    / "damp_two_repair_trace_transport_regression.json"
)


def one_new_support(
    base_shadow: int,
    extension_shadow: int,
) -> int | None:
    """Return the unique new shattered support, if there is one."""

    difference = extension_shadow & ~base_shadow
    if difference.bit_count() != 1:
        return None
    return difference.bit_length() - 1


def repair_support(base: int, tip: int, dimension: int) -> int:
    """Return the old directions in which ``tip`` has a neighbour in base."""

    support = 0
    for coordinate in range(dimension):
        if base & (1 << (tip ^ (1 << coordinate))):
            support |= 1 << coordinate
    return support


def opposite_layer_lift(
    lower: int,
    upper: int,
    dimension: int,
) -> int:
    """Place ``lower`` in layer zero and ``upper`` in layer one."""

    outer = 1 << dimension
    result = 0
    for word in members(lower, dimension):
        result |= 1 << word
    for word in members(upper, dimension):
        result |= 1 << (word | outer)
    return result


def sections(lift: int, dimension: int) -> tuple[int, int]:
    """Recover the projected zero and one sections in the outer coordinate."""

    outer = 1 << dimension
    lower = 0
    upper = 0
    for word in members(lift, dimension + 1):
        if word & outer:
            upper |= 1 << (word ^ outer)
        else:
            lower |= 1 << word
    return lower, upper


def expected_lift_shadow(
    base_shadow: int,
    union_shadow: int,
    dimension: int,
) -> int:
    """Return Sigma^+ union {D+e:D in Sigma}."""

    outer = 1 << dimension
    result = union_shadow
    for support in range(1 << dimension):
        if base_shadow & (1 << support):
            result |= 1 << (support | outer)
    return result


def fibre(
    family: int,
    dimension: int,
    support: int,
    trace: int,
) -> frozenset[int]:
    """Return one trace fibre as a set of cube words."""

    return frozenset(
        word
        for word in members(family, dimension)
        if word & support == trace
    )


def verify_trace_table(
    base_shadow: int,
    union_shadow: int,
    lower: int,
    upper: int,
    lift: int,
    dimension: int,
) -> int:
    """Check every old and vertical trace-fibre transport formula."""

    outer = 1 << dimension
    checks = 0

    for support in range(1 << dimension):
        if not union_shadow & (1 << support):
            continue
        for trace in support_traces(support, dimension):
            expected = fibre(lower, dimension, support, trace) | frozenset(
                word | outer
                for word in fibre(upper, dimension, support, trace)
            )
            actual = fibre(lift, dimension + 1, support, trace)
            if actual != expected:
                raise AssertionError("old-support trace transport failed")
            checks += 1

    for support in range(1 << dimension):
        if not base_shadow & (1 << support):
            continue
        vertical = support | outer
        for trace in support_traces(support, dimension):
            expected_lower = fibre(lower, dimension, support, trace)
            actual_lower = fibre(lift, dimension + 1, vertical, trace)
            if actual_lower != expected_lower:
                raise AssertionError("lower vertical trace transport failed")
            checks += 1

            expected_upper = frozenset(
                word | outer
                for word in fibre(upper, dimension, support, trace)
            )
            actual_upper = fibre(
                lift,
                dimension + 1,
                vertical,
                trace | outer,
            )
            if actual_upper != expected_upper:
                raise AssertionError("upper vertical trace transport failed")
            checks += 1

    return checks


def verify_distinct_pair(
    base: int,
    first_tip: int,
    second_tip: int,
    dimension: int,
    base_shadow: int,
    first_shadow: int,
    second_shadow: int,
    first_support: int,
    second_support: int,
) -> dict[str, int]:
    """Check the complete theorem packet for one ordered repair pair."""

    # The paper uses A=S+{a} in the upper layer and B=S+{b} below.
    upper = base | (1 << first_tip)
    lower = base | (1 << second_tip)
    union = upper | lower
    lift = opposite_layer_lift(lower, upper, dimension)
    outer = 1 << dimension

    union_shadow = shattered_coordinate_sets_mask(union, dimension)
    expected_union_shadow = (
        base_shadow
        | (1 << first_support)
        | (1 << second_support)
    )
    if union_shadow != expected_union_shadow:
        raise AssertionError("the two-tip union has the wrong shadow")

    lift_shadow = shattered_coordinate_sets_mask(lift, dimension + 1)
    predicted_shadow = expected_lift_shadow(
        base_shadow,
        union_shadow,
        dimension,
    )
    if lift_shadow != predicted_shadow:
        raise AssertionError("the lift has the wrong shadow")
    if not is_ample_family(lift, dimension + 1):
        raise AssertionError("distinct repair supports did not give an ample lift")

    actual_nonfaces = set(minimal_nonfaces(lift_shadow, dimension + 1))
    predicted_nonfaces = set(minimal_nonfaces(union_shadow, dimension)) | {
        first_support | outer,
        second_support | outer,
    }
    if actual_nonfaces != predicted_nonfaces:
        raise AssertionError("the lift minimal nonfaces have the wrong form")

    lift_signs = signing(lift, lift_shadow, dimension + 1)
    union_signs = signing(union, union_shadow, dimension)
    predicted_signs = dict(union_signs)
    predicted_signs[first_support | outer] = first_tip & first_support
    predicted_signs[second_support | outer] = (
        (second_tip & second_support) | outer
    )
    if lift_signs != predicted_signs:
        raise AssertionError("the two vertical signed clauses are wrong")
    if satisfying_family(predicted_signs, dimension + 1) != lift:
        raise AssertionError("the transported signing does not recover the lift")

    recovered_lower, recovered_upper = sections(lift, dimension)
    if recovered_lower != lower or recovered_upper != upper:
        raise AssertionError("section recovery failed")
    recovered_base = recovered_lower & recovered_upper
    if recovered_base != base:
        raise AssertionError("the paired core was not recovered by intersection")
    recovered_first = members(recovered_upper & ~recovered_lower, dimension)
    recovered_second = members(recovered_lower & ~recovered_upper, dimension)
    if recovered_first != (first_tip,) or recovered_second != (second_tip,):
        raise AssertionError("the two repair tips were not recovered")
    if repair_support(recovered_base, first_tip, dimension) != first_support:
        raise AssertionError("the first support was not recovered locally")
    if repair_support(recovered_base, second_tip, dimension) != second_support:
        raise AssertionError("the second support was not recovered locally")

    trace_checks = verify_trace_table(
        base_shadow,
        union_shadow,
        lower,
        upper,
        lift,
        dimension,
    )
    return {
        "signed_clauses_checked": len(predicted_signs),
        "trace_fibres_checked": trace_checks,
    }


def dimension_profile(dimension: int) -> dict[str, int]:
    """Check every eligible ordered pair over every ample base."""

    counts: defaultdict[str, int] = defaultdict(int)
    ambient = 1 << dimension
    for base in range(1, 1 << ambient):
        if not is_ample_family(base, dimension):
            continue
        counts["ample_bases"] += 1
        base_shadow = shattered_coordinate_sets_mask(base, dimension)
        repairs: list[tuple[int, int, int]] = []
        for tip in range(ambient):
            if base & (1 << tip):
                continue
            extension = base | (1 << tip)
            if not is_ample_family(extension, dimension):
                continue
            extension_shadow = shattered_coordinate_sets_mask(
                extension,
                dimension,
            )
            support = one_new_support(base_shadow, extension_shadow)
            if support is None:
                raise AssertionError("an ample one-tip extension has no unique support")
            if support != repair_support(base, tip, dimension):
                raise AssertionError("the local repair-support formula failed")
            repairs.append((tip, extension_shadow, support))
        counts["ample_one_tip_extensions"] += len(repairs)

        for first_tip, first_shadow, first_support in repairs:
            for second_tip, second_shadow, second_support in repairs:
                if first_tip == second_tip:
                    continue
                counts["ordered_repair_pairs"] += 1
                upper = base | (1 << first_tip)
                lower = base | (1 << second_tip)
                lift = opposite_layer_lift(lower, upper, dimension)
                actual_ample = is_ample_family(lift, dimension + 1)
                predicted_ample = first_support != second_support
                if actual_ample != predicted_ample:
                    raise AssertionError("support inequality is not the ampleness test")
                if not predicted_ample:
                    counts["equal_support_pairs"] += 1
                    continue
                counts["distinct_support_pairs"] += 1
                checked = verify_distinct_pair(
                    base,
                    first_tip,
                    second_tip,
                    dimension,
                    base_shadow,
                    first_shadow,
                    second_shadow,
                    first_support,
                    second_support,
                )
                counts["signed_clauses_checked"] += checked[
                    "signed_clauses_checked"
                ]
                counts["trace_fibres_checked"] += checked[
                    "trace_fibres_checked"
                ]
                counts["canonical_recoveries"] += 1

    return dict(sorted(counts.items()))


def main() -> int:
    mathematics = {
        "schema": "damp-two-repair-trace-transport-regression-v1",
        "scope": (
            "every eligible ordered pair of ample one-concept extensions "
            "over every ample base through dimension four; regression only"
        ),
        "dimensions": {
            str(dimension): dimension_profile(dimension)
            for dimension in range(1, 5)
        },
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
    print(json.dumps(mathematics["dimensions"], indent=2, sort_keys=True))
    print(f"wrote {OUTPUT}")
    print(f"mathematical_profile_sha256={digest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
