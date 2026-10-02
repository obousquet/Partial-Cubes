#!/usr/bin/env python3
"""Check two fixed nonpeelable controls for the low-support addition theorem.

Attach a path to the certified 291-concept seed and complete its square;
then complete one edge of the missing square fibre in the product of
that seed with a three-vertex path. These are fixed controls, not a search
for new obstructions. Nonpeelability is certified by complete deletion
state graphs. The uniform conclusions rely on the manuscript proofs.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from verify_damp_singleton_recovery_implications import shadow
from verify_damp_terminal_clause_recovery import corner

ROOT = Path(__file__).resolve().parent
DIMENSION = 14


def support(family, tip):
    return sum(1 << f for f in range(DIMENSION) if tip ^ (1 << f) in family)


def negative_certificate(family):
    initial = set(family)
    states = {}

    def visit(remaining):
        key = frozenset(remaining)
        if key in states:
            return
        assert remaining, 'A peeling was found: the proposed control is false.'
        assert len(states) < 100, 'Fixed-control state budget exceeded.'
        choices = [x for x in sorted(remaining) if corner(remaining, x, DIMENSION)]
        states[key] = choices
        for x in choices:
            visit(remaining - {x})

    visit(initial)
    rows = []
    for remaining, choices in states.items():
        assert all(frozenset(remaining - {x}) in states for x in choices)
        rows.append({'removed': sorted(initial - remaining), 'corners': choices})
    rows.sort(key=lambda row: (len(row['removed']), row['removed']))
    return {'order': len(initial), 'status': 'nonpeelable', 'states': rows}


def main():
    source_path = ROOT / 'damp_fixed_291_repair_creation.json'
    source = json.loads(source_path.read_text())
    base = set(source['base'])
    assert len(base) == len(shadow(base, 4095)) == 291
    assert not any(corner(base, x, 12) for x in base)
    model_path = ROOT / 'damp_hall_carrier_defect_min8.json'
    assert hashlib.sha256(model_path.read_bytes()).hexdigest() == source['source_model_sha256']
    g = min(base)
    h = next(g ^ (1 << f) for f in range(12) if g ^ (1 << f) in base)
    p, q = 1 << 12, 1 << 13
    # A path attached at g. Adding a completes its square; the opposite
    # vertex z becomes a new corner and is interchangeable with a.
    path = base | {g | p, g | p | q}
    a, z = g | q, g | p
    square = path | {a}
    assert len(shadow(path, 16383)) == len(path)
    assert len(shadow(square, 16383)) == len(square)
    assert support(path, a) == p | q
    assert not corner(path, z, DIMENSION) and corner(square, z, DIMENSION)
    swapped = lambda x: z if x == a else a if x == z else x
    assert {swapped(x) for x in square - {z}} == path
    assert {x for x in square if (x ^ a).bit_count() == 1} == {
        x for x in square if (x ^ z).bit_count() == 1}
    minor_images = []
    for f in range(DIMENSION):
        bit = 1 << f
        for side in [None, 0, 1]:
            small = {x & ~bit for x in path
                     if side is None or (x >> f & 1) == side}
            large = {x & ~bit for x in square
                     if side is None or (x >> f & 1) == side}
            gain = large - small
            assert small <= large and len(gain) <= 1
            acquired = None
            if gain:
                tip = next(iter(gain))
                acquired = support(small, tip)
                assert acquired.bit_count() <= 2
                assert corner(large, tip, DIMENSION)
            minor_images.append({'coordinate': f, 'side': side,
                                 'image_gain': sorted(gain),
                                 'support': acquired})
    # Product with a punctured square, followed by two adjacent completions.
    product = {x | t for x in base for t in [0, p, q]}
    c, b = g | p | q, h | p | q
    one = product | {c}
    two = one | {b}
    assert not any(corner(product, x, DIMENSION) for x in product)
    assert len(shadow(product, 16383)) == len(product)
    assert len(shadow(one, 16383)) == len(one)
    assert len(shadow(two, 16383)) == len(two)
    assert support(product, c) == support(product, b) == p | q
    assert support(one, b) == p | q | (g ^ h)
    assert [x for x in sorted(one) if corner(one, x, DIMENSION)] == [c]
    controls = {
        'attached_path': negative_certificate(path),
        'completed_attached_square': negative_certificate(square),
        'one_square_completion': negative_certificate(one),
        'two_adjacent_square_completions': negative_certificate(two),
    }
    result = {
        'scope': __doc__.strip(),
        'source_file_sha256': hashlib.sha256(source_path.read_bytes()).hexdigest(),
        'source_model_sha256': source['source_model_sha256'],
        'base_order': len(base), 'fixed_base_edge': [g, h],
        'new_coordinates': [12, 13],
        'attached_square_tip': a, 'new_opposite_corner': z,
        'elementary_image_checks': minor_images,
        'product_tips': [c, b], 'controls': controls,
        'all_checks_passed': True,
    }
    (ROOT / 'damp_rank_two_additions_verification.json').write_text(
        json.dumps(result, indent=2) + '\n')
    print(json.dumps({
        'all_checks_passed': True,
        'controls': {name: {'order': row['order'], 'deletion_states': len(row['states'])}
                     for name, row in controls.items()},
    }, indent=2))


if __name__ == '__main__':
    main()
