#!/usr/bin/env python3
"""Verify fixed controls for the repair-support strong-projection theorem.

Check every coordinate subset in six already recorded constructions, with
no peeling search. Stored positive orders also give explicit peelings of
all twelve one-coordinate overlaps of the 291-concept seed. The uniform
theorem and the conditional three-tip exclusion are proved in the manuscript.
"""
from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path

from verify_damp_layer_midpoint_cores import corner, shadow, subsets

ROOT = Path(__file__).resolve().parent


def strong(family, mask):
    """Suppress coordinates in mask, retaining their now-zero bit positions."""
    counts = Counter(x & ~mask for x in family)
    return {x for x, count in counts.items() if count == 1 << mask.bit_count()}


def incident(family, tip, dimension):
    return sum(1 << e for e in range(dimension) if tip ^ (1 << e) in family)


def check_order(family, order, dimension):
    remaining = set(family)
    assert len(order) == len(remaining) and set(order) == remaining
    for x in order:
        assert corner(remaining, x, dimension), (x, len(remaining))
        remaining.remove(x)
    assert not remaining


def projected_order(family, order, mask, dimension):
    projected = strong(family, mask)
    first = {}
    for index, x in enumerate(order):
        first.setdefault(x & ~mask, index)
    answer = sorted(projected, key=first.__getitem__)
    check_order(projected, answer, dimension)
    return answer


def check_formula(name, family, tip, dimension, support, order=None):
    assert tip not in family
    assert incident(family, tip, dimension) == support
    assert all(tip ^ u in family for u in subsets(support) if u)
    # Absence of the tip's trace certifies that support is not shattered.
    assert all((x & support) != (tip & support) for x in family)
    extended = family | {tip}
    if order is not None:
        check_order(extended, order, dimension)
    changed = 0
    for mask in range(1 << dimension):
        small = strong(family, mask)
        large = strong(extended, mask)
        if mask & ~support:
            assert small == large
        else:
            word = tip & ~mask
            assert word not in small and large == small | {word}
            assert incident(small, word, dimension) == support & ~mask
            assert corner(large, word, dimension)
            changed += 1
    assert changed == 1 << support.bit_count()
    return dict(name=name, base_order=len(family), dimension=dimension,
                tip=tip, acquired_support=support,
                coordinate_subsets_checked=1 << dimension,
                changed_projections=changed, peeling_checked=order is not None)


def main():
    source_path = ROOT / 'damp_fixed_291_repair_creation.json'
    source = json.loads(source_path.read_text())
    base = set(source['base'])
    assert len(base) == 291 and len(shadow(base, 12)) == 291
    assert not any(corner(base, x, 12) for x in base)
    by_tip = {row['tip']: row for row in source['extensions']}
    repairs = [91, 157, 319]
    for tip in repairs:
        row = by_tip[tip]
        assert row['status'] == 'peelable'
        check_order(base | {tip}, row['peeling'], 12)
    intersection = 4095
    for tip in repairs:
        intersection &= by_tip[tip]['support']
    assert intersection == 0
    overlaps = []
    for e in range(12):
        mask = 1 << e
        tip = next(t for t in repairs if not by_tip[t]['support'] & mask)
        assert strong(base, mask) == strong(base | {tip}, mask)
        order = projected_order(base | {tip}, by_tip[tip]['peeling'], mask, 12)
        overlaps.append(dict(coordinate=e, repair_used=tip,
                             order=len(order), peeling=order))

    cases = [check_formula('original resolving addition', base, 91, 12,
                           by_tip[91]['support'], by_tip[91]['peeling']),
             check_formula('recorded rank-three addition', base, 2696, 12, 2312)]
    one_corner = base | {2696}
    assert [x for x in one_corner if corner(one_corner, x, 12)] == [2696]
    for tip in [2184, 3720]:
        gained = 2312 | (2696 ^ tip)
        cases.append(check_formula('recorded adjacent pair', one_corner,
                                   tip, 12, gained))

    tower_projections = []
    for r in [1, 2]:
        layers = range(1 << r)
        layer_mask = ((1 << r) - 1) << 12
        family = {x | (t << 12) for x in base for t in layers}
        family |= {91 | (t << 12) for t in layers if t != 0}
        family |= {157 | (t << 12) for t in layers if t != (1 << r) - 1}
        order = [157 | (t << 12) for t in reversed(range((1 << r) - 1))]
        order += [x | (t << 12) for x in by_tip[91]['peeling'] for t in layers]
        cases.append(check_formula('terminal tower, layers=' + str(r), family,
                                   91, 12 + r, 102 | layer_mask, order))
        # A fixed nonpeelable projection: it equals the cornerless base.
        assert strong(family, layer_mask) == base
        assert strong(family | {91}, layer_mask) == base | {91}
        tower_projections.append(dict(layers=r, projection_mask=layer_mask,
                                      base_projection_order=291,
                                      inherited_repair_support=102))

    report = dict(scope=__doc__.strip(),
                  source_file_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
                  cases=cases, base_overlap_peelings=overlaps,
                  nonpeelable_tower_projections=tower_projections,
                  all_checks_passed=True)
    (ROOT / 'damp_repair_strong_projections_verification.json').write_text(
        json.dumps(report, indent=2) + '\n')
    print(json.dumps(dict(all_checks_passed=True, cases=len(cases),
                          coordinate_subsets_checked=sum(c['coordinate_subsets_checked']
                                                         for c in cases),
                          base_overlap_peelings=len(overlaps)), indent=2))


if __name__ == '__main__':
    main()
