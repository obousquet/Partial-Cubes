#!/usr/bin/env python3
"""Fixed local checks for adjacent square-completion invariance.

Attach a certified cornerless nonpeelable class to a small punctured cube,
leaving the base cornered. Check the first-deletion reductions, commuting
second deletions, and opposite-edge isomorphism without a peeling search.
The attached class remains a proper pc-minor; these are not new forbidden
witnesses. The manuscript supplies the uniform induction proof.
"""
from __future__ import annotations

from collections import Counter
import hashlib
import json
from pathlib import Path

from verify_damp_layer_midpoint_cores import corner, shadow
from verify_damp_rank_three_corner_neighbors import check_block_shadow

ROOT = Path(__file__).resolve().parent
CUBE = set(range(8))
TIPS = {0, 4}
OPPOSITE = {3, 7}
NEAR = {1, 2, 5, 6}


def corners(family, dimension):
    return {x for x in family if corner(family, x, dimension)}


def support(family, tip, dimension):
    return sum(1 << e for e in range(dimension) if tip ^ (1 << e) in family)


def check_reductions(base, dimension):
    full = base | TIPS
    first = corners(full, dimension)
    rows = []
    for z in sorted(first):
        after = full - {z}
        row = dict(first=z)
        if z in TIPS:
            other = next(iter(TIPS - {z}))
            assert support(base, other, dimension) == 3
            row['case'] = 'one square completion remains'
        elif z not in CUBE:
            assert corner(base, z, dimension)
            assert support(base - {z}, 0, dimension) == 3
            assert support(base - {z}, 4, dimension) == 3
            row['case'] = 'smaller instance with unchanged punctured cubes'
        elif z in NEAR:
            assert corner(base, z, dimension)
            small = base - {z}
            retained_tip = 4 if z < 4 else 0
            last_tip = 4 ^ retained_tip
            middle = small | {retained_tip}
            assert support(small, retained_tip, dimension) == 3
            assert corner(middle, retained_tip, dimension)
            assert support(middle, last_tip, dimension).bit_count() == 2
            assert corner(after, last_tip, dimension)
            row.update(case='two support-two additions after a common corner',
                       addition_order=[retained_tip, last_tip])
        else:
            assert z in OPPOSITE
            second = corners(after, dimension)
            assert not second & TIPS
            second_rows = []
            for w in sorted(second):
                if w not in CUBE or w in NEAR:
                    assert corner(full, w, dimension)
                    assert corner(full - {w}, z, dimension)
                    assert corner(base, w, dimension)
                    second_rows.append(dict(second=w, case='commutes to first'))
                else:
                    assert w in OPPOSITE - {z}
                    assert not support(full, w, dimension) & ~7
                    assert not support(full, z, dimension) & ~7
                    switch = {0: 3, 3: 0, 4: 7, 7: 4}
                    image = lambda x: switch.get(x, x)
                    remaining = after - {w}
                    assert {image(x) for x in remaining} == base
                    # Verify the claimed graph isomorphism on both edge sets.
                    edges = {tuple(sorted((x, x ^ (1 << e))))
                             for x in remaining for e in range(dimension)
                             if x ^ (1 << e) in remaining}
                    mapped = {tuple(sorted((image(x), image(y)))) for x, y in edges}
                    base_edges = {tuple(sorted((x, x ^ (1 << e))))
                                  for x in base for e in range(dimension)
                                  if x ^ (1 << e) in base}
                    assert mapped == base_edges
                    second_rows.append(dict(second=w, case='isomorphic base'))
            row.update(case='opposite edge first', second_deletions=second_rows)
        rows.append(row)
    return rows


def check_minor_images(base, dimension):
    counts = Counter()
    full = base | TIPS
    for e in range(dimension):
        bit = 1 << e
        for side in [None, 0, 1]:
            def project(family):
                return {x & ~bit for x in family
                        if side is None or (x >> e & 1) == side}
            small, large = project(base), project(full)
            gain = large - small
            assert small <= large and len(gain) <= 2
            if not gain:
                kind = 'unchanged'
            elif len(gain) == 1:
                tip = next(iter(gain))
                assert support(small, tip, dimension).bit_count() <= 2
                assert corner(large, tip, dimension)
                kind = 'one support-at-most-two addition'
            else:
                c, u = sorted(gain)
                assert (c ^ u).bit_count() == 1
                common = support(small, c, dimension)
                assert common == support(small, u, dimension)
                assert 1 <= common.bit_count() <= 2
                assert corner(small | {c}, c, dimension)
                assert corner(small | {u}, u, dimension)
                assert corner(large, c, dimension) and corner(large, u, dimension)
                if common.bit_count() == 1:
                    assert support(small | {c}, u, dimension).bit_count() == 2
                kind = 'adjacent pair with common support ' + str(common.bit_count())
            counts[kind] += 1
    assert sum(counts.values()) == 3 * dimension
    return dict(counts)


def main():
    source_path = ROOT / 'damp_fixed_291_repair_creation.json'
    source = json.loads(source_path.read_text())
    block = set(source['base'])
    block_shadow = shadow(block, 12)
    assert 0 in block and len(block) == len(block_shadow) == 291
    assert not corners(block, 12)
    controls = []
    for anchor, leaf in [(1, False), (3, False), (1, True)]:
        local = CUBE - TIPS
        local_dimension = 4 if leaf else 3
        if leaf:
            local |= {10}
        base = local | {anchor | (x << local_dimension) for x in block}
        dimension = local_dimension + 12
        assert corners(base, dimension)
        assert support(base, 0, dimension) == support(base, 4, dimension) == 3
        assert {x & 3 for x in base} == {1, 2, 3}
        for tips in [set(), {0}, {4}, TIPS]:
            family = base | tips
            check_block_shadow(family, local | tips, local_dimension,
                               [anchor], block, block_shadow)
        # Fix all local coordinates to the attachment anchor.
        assert {x >> local_dimension for x in base
                if x & ((1 << local_dimension) - 1) == anchor} == block
        controls.append(dict(attachment=anchor, outside_leaf=leaf,
                             dimension=dimension, base_order=len(base),
                             base_corners=sorted(corners(base, dimension)),
                             reductions=check_reductions(base, dimension),
                             elementary_minor_images=check_minor_images(base, dimension)))
    cases = {r['case'] for c in controls for r in c['reductions']}
    second_cases = {r['case'] for c in controls for first in c['reductions']
                    for r in first.get('second_deletions', [])}
    assert len(cases) == 4
    assert second_cases == {'commutes to first', 'isomorphic base'}
    assert any(r['second'] not in CUBE for c in controls
               for first in c['reductions'] for r in first.get('second_deletions', []))
    report = dict(scope=__doc__.strip(),
                  source_file_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
                  controls=controls, all_checks_passed=True)
    (ROOT / 'damp_square_pair_invariance_verification.json').write_text(
        json.dumps(report, indent=2) + '\n')
    print(json.dumps(dict(all_checks_passed=True, controls=len(controls),
                          first_deletion_cases=len(cases),
                          elementary_minor_images=sum(3 * c['dimension'] for c in controls),
                          second_deletion_cases=sorted(second_cases)), indent=2))


if __name__ == '__main__':
    main()
