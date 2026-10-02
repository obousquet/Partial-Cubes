#!/usr/bin/env python3
"""Verify separated cooperative additions using two fixed Hall copies.

The controls have bridge lengths zero and three. Both additions together
peel, whereas each individual addition retains a proper cornerless Hall
minor. The completed cubes have disjoint supports and are disjoint. This
refutes geometric isolation as a substitute for pc-minimality; it is not
a counterexample with a forbidden base. No family census is performed.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from verify_damp_singleton_recovery_implications import shadow
from verify_damp_terminal_clause_recovery import corner, subsets

ROOT = Path(__file__).resolve().parent
BLOCK = 12
MASK = (1 << BLOCK) - 1


def support(family, x, dimension):
    return sum(1 << e for e in range(dimension) if x ^ (1 << e) in family)


def check_ampleness(family, left, right, bridge, length):
    """Use exact block projections and the absence of mixed shattered pairs.

    Every support meeting two blocks contains such a pair. The bridge
    projection is a path and hence shatters only empty/singleton supports.
    Thus the entire shattered complex is determined without enumerating
    an ambient cube of dimension 24+length.
    """
    assert {x & MASK for x in family} == left
    assert {(x >> BLOCK) & MASK for x in family} == right
    assert {x >> (2 * BLOCK) for x in family} == set(bridge)
    blocks = [set(range(BLOCK)), set(range(BLOCK, 2 * BLOCK)),
              set(range(2 * BLOCK, 2 * BLOCK + length))]
    checked_pairs = 0
    for i, first in enumerate(blocks):
        for second in blocks[i + 1:]:
            for e in first:
                for f in second:
                    pair = (1 << e) | (1 << f)
                    assert len({x & pair for x in family}) < 4
                    checked_pairs += 1
    left_shadow = shadow(left, MASK)
    right_shadow = shadow(right, MASK)
    bridge_shadow = shadow(set(bridge), (1 << length) - 1)
    assert len(bridge_shadow) == length + 1
    shadow_order = len(left_shadow) + len(right_shadow) - 1 + length
    assert shadow_order == len(family)
    return {'shadow_order': shadow_order, 'mixed_pairs_checked': checked_pairs}


def main():
    proof_path = ROOT / 'damp_terminal_clause_recovery_verification.json'
    proof = json.loads(proof_path.read_text())
    source_path = ROOT / proof['source_file']
    assert hashlib.sha256(source_path.read_bytes()).hexdigest() == proof['source_file_sha256']
    source = json.loads(source_path.read_text())
    bits = int(source['source_family_hex'], 16)
    hall = {x for x in range(1 << BLOCK) if bits >> x & 1}
    record = next(row for row in proof['repair_proofs'] if row['tip'] == 183)
    order = record['peeling']
    assert order[-1] == MASK and MASK in hall
    remaining = hall | {183}
    for x in order:
        assert x in remaining and corner(remaining, x, BLOCK)
        remaining.remove(x)
    assert not remaining
    base = {x ^ MASK for x in hall}
    tip = 183 ^ MASK
    acquired = sum(1 << e for e in record['support'])
    assert len(base) == len(shadow(base, MASK)) == 299
    assert not any(corner(base, x, BLOCK) for x in base)
    assert (tip & ~acquired).bit_count() == 4
    records = []
    for length in [0, 3]:
        dimension = 2 * BLOCK + length
        bridge = [(1 << k) - 1 for k in range(length + 1)]
        far = ((1 << length) - 1) << (2 * BLOCK)
        family = (set(base) | {(x << BLOCK) | far for x in base}
                  | {t << (2 * BLOCK) for t in bridge})
        a, b = tip, (tip << BLOCK) | far
        both = family | {a, b}
        assert len(family) == 597 + length
        assert not any(corner(family, x, dimension) for x in family)
        controls = []
        for label, added in [('base', set()), ('left', {a}),
                             ('right', {b}), ('both', {a, b})]:
            current = family | added
            left = base | ({tip} if a in added else set())
            right = base | ({tip} if b in added else set())
            controls.append({'label': label, 'order': len(current),
                             **check_ampleness(current, left, right, bridge, length)})
        assert support(family, a, dimension) == acquired
        assert support(family, b, dimension) == acquired << BLOCK
        cube_a = {a ^ u for u in subsets(acquired)}
        cube_b = {b ^ u for u in subsets(acquired << BLOCK)}
        assert cube_a <= both and cube_b <= both and not cube_a & cube_b
        cube_distance = min((x ^ y).bit_count() for x in cube_a for y in cube_b)
        assert cube_distance == 8 + length
        # Retain the opposite Hall copy as a proper signed restriction.
        left_augmented = family | {a}
        right_minor = {((x >> BLOCK) & MASK) for x in left_augmented
                       if x & MASK == 0 and x >> (2 * BLOCK) == (1 << length) - 1}
        right_augmented = family | {b}
        left_minor = {x & MASK for x in right_augmented
                      if (x >> BLOCK) & MASK == 0 and x >> (2 * BLOCK) == 0}
        assert right_minor == left_minor == base
        # Peel each completed copy to its prescribed attachment vertex,
        # then peel the remaining path from its left endpoint.
        peeling = ([x ^ MASK for x in order[:-1]]
                   + [((x ^ MASK) << BLOCK) | far for x in order[:-1]]
                   + [t << (2 * BLOCK) for t in bridge])
        remaining = set(both)
        for x in peeling:
            assert x in remaining and corner(remaining, x, dimension)
            remaining.remove(x)
        assert not remaining
        records.append({'bridge_length': length, 'dimension': dimension,
                        'base_order': len(family), 'tips': [a, b],
                        'supports': [acquired, acquired << BLOCK],
                        'cube_distance': cube_distance,
                        'tip_distance': (a ^ b).bit_count(),
                        'ampleness_controls': controls,
                        'base_cornerless': True,
                        'single_additions_have_proper_Hall_minor': True,
                        'both_additions_peel': peeling})
    result = {'scope': __doc__.strip(),
              'source_file_sha256': proof['source_file_sha256'],
              'peeling_source_sha256': hashlib.sha256(proof_path.read_bytes()).hexdigest(),
              'original_tip': 183, 'attachment_vertex': MASK,
              'records': records, 'all_checks_passed': True}
    (ROOT / 'damp_separated_cooperative_repairs_verification.json').write_text(
        json.dumps(result, indent=2) + '\n')
    print(json.dumps({'all_checks_passed': True, 'records': [
        {key: row[key] for key in ['bridge_length', 'base_order', 'cube_distance', 'tip_distance']}
        for row in records]}, indent=2))


if __name__ == '__main__':
    main()
