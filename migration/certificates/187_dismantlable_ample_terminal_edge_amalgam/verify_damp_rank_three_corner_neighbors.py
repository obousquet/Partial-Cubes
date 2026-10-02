#!/usr/bin/env python3
"""Fixed checks for the two forced deletions at a rank-three neighbor.

The two constructed controls attach copies of the fixed cornerless ample
291 family at selected vertices of a small ample family. They have a proper
nonpeelable minor and are deliberately not forbidden witnesses. One tests
the forced two-deletion prefix; the other tests its companion obstruction.
Also replay the 15 rank-three neighbors of the nine already recorded
one-corner extensions. No peeling search or expanded seed atlas is used.
"""
from __future__ import annotations

import hashlib
import json
from itertools import combinations
from pathlib import Path

from verify_damp_layer_midpoint_cores import corner, shadow, subsets

ROOT = Path(__file__).resolve().parent


def directions(family, word, dimension):
    return sum(1 << e for e in range(dimension) if word ^ (1 << e) in family)


def corners(family, dimension):
    return [v for v in sorted(family) if corner(family, v, dimension)]


def check_local(family, c, u, dimension):
    assert corners(family, dimension) == [c]
    I = directions(family, c, dimension)
    f = c ^ u
    assert f.bit_count() == 1 and not f & I
    P = directions(family, u, dimension)
    assert P.bit_count() == 3 and P & f
    J = P ^ f
    assert J & I == J and J != I
    assert all(u ^ s in family for s in subsets(P) if s)
    upper = family | {u}
    x, y = u ^ J, c ^ J
    cs = corners(upper, dimension)
    expected = [u] + ([x] if not directions(upper, x, dimension) & ~P else [])
    assert cs == sorted(expected)
    after = None
    outside = directions(family, y, dimension) & ~I
    if x in cs:
        after = corners(upper - {x}, dimension)
        assert after == ([y] if outside == f else [])
    return {'corner': c, 'corner_support': I, 'tip': u, 'acquired_support': P,
            'first_old_vertex': x, 'second_old_vertex': y, 'upper_corners': cs,
            'corners_after_first_old_deletion': after,
            'outside_directions_at_second_vertex': outside}


def attached_family(local, local_dimension, anchors, block):
    result = set(local)
    for k, a in enumerate(anchors):
        shift = local_dimension + 12 * k
        result.update(a | (v << shift) for v in block)
    return result, local_dimension + 12 * len(anchors)


def check_block_shadow(family, local, local_dimension, anchors, block, block_shadow):
    """Exact shattering count using projections and unshattered mixed pairs."""
    local_mask = (1 << local_dimension) - 1
    assert {v & local_mask for v in family} == local
    blocks = [list(range(local_dimension))]
    for k in range(len(anchors)):
        shift = local_dimension + 12 * k
        assert {(v >> shift) & 4095 for v in family} == block
        blocks.append(list(range(shift, shift + 12)))
    pair_count = 0
    for first, second in combinations(blocks, 2):
        for e in first:
            for f in second:
                mask = (1 << e) | (1 << f)
                assert len({v & mask for v in family}) < 4
                pair_count += 1
    total = len(shadow(local, local_dimension)) + len(anchors) * (len(block_shadow) - 1)
    assert total == len(family)
    return {'shattered_support_count': total, 'unshattered_mixed_pairs': pair_count}


def check_state_certificate(initial, dimension, rows):
    by_removed = {frozenset(r['removed']): r for r in rows}
    assert len(by_removed) == len(rows) and frozenset() in by_removed
    for removed, row in by_removed.items():
        family = initial - removed
        assert family and removed <= initial
        assert corners(family, dimension) == row['corners']
        for c in row['corners']:
            assert removed | {c} in by_removed


def main():
    source_path = ROOT / 'damp_fixed_291_repair_creation.json'
    source = json.loads(source_path.read_text())
    block = set(source['base'])
    block_shadow = shadow(block, 12)
    assert 0 in block and len(block_shadow) == len(block) == 291
    assert not corners(block, 12)
    models = []
    for with_companion in (False, True):
        c, u, x, y = 0, 8, 11, 3
        local = set(range(8)) | {9, 10, 11}
        local_dimension = 4
        if with_companion:
            local |= {16 + v for v in range(1, 8)}
            local_dimension = 5
        assert len(shadow(local, local_dimension)) == len(local)
        anchors = sorted((set(corners(local, local_dimension)) - {c}) | {4})
        family, dimension = attached_family(local, local_dimension, anchors, block)
        ampleness = check_block_shadow(family, local, local_dimension, anchors,
                                      block, block_shadow)
        assert corners(family, dimension) == [c]
        assert not corners(family - {c}, dimension)
        case = check_local(family, c, u, dimension)
        assert case['upper_corners'] == [u, x]
        upper = family | {u}
        check_block_shadow(upper, local | {u}, local_dimension, anchors, block, block_shadow)
        rows = [{'removed': [], 'corners': [u, x]},
                {'removed': [u], 'corners': [c]},
                {'removed': sorted([u, c]), 'corners': []}]
        if with_companion:
            assert case['corners_after_first_old_deletion'] == []
            assert case['outside_directions_at_second_vertex'] == 8 + 16
            rows.append({'removed': [x], 'corners': []})
            predecessor = family - {c}
            assert directions(predecessor, 16, dimension) == 7
            assert all(16 ^ s in predecessor for s in subsets(7) if s)
            assert len({v & 7 for v in predecessor}) == 7
        else:
            assert case['corners_after_first_old_deletion'] == [y]
            assert case['outside_directions_at_second_vertex'] == 8
            rows.extend([{'removed': [x], 'corners': [y]},
                         {'removed': sorted([x, y]), 'corners': []}])
        check_state_certificate(upper, dimension, rows)
        # Fix the local word and all other attachment blocks to zero:
        # the selected attachment copy is a proper restriction of the family.
        shift = local_dimension
        mask = 4095 << shift
        restricted = {(v >> shift) & 4095 for v in family
                      if v & ~mask == anchors[0]}
        assert restricted == block and len(family) > len(block)
        models.append({'has_companion': with_companion, 'order': len(family),
                       'dimension': dimension, 'local_family': sorted(local),
                       'attachment_anchors': anchors, **ampleness, **case,
                       'negative_certificate': rows,
                       'has_proper_cornerless_291_minor': True})

    # Recompute the already stored rank-three comparisons without searching
    # for a peeling. This checks the formula when the first old deletion is blocked.
    old_path = ROOT / 'computations/round767_recorded_rank_three_neighbors.json'
    old = json.loads(old_path.read_text())
    replayed = []
    for row in old['rows']:
        family = block | {row['c']}
        assert len(shadow(family, 12)) == len(family)
        assert len(shadow(family | {row['u']}, 12)) == len(family) + 1
        result = check_local(family, row['c'], row['u'], 12)
        assert result['upper_corners'] == row['corners']
        replayed.append(result)
    assert len(replayed) == 15
    report = {'schema': 'damp-rank-three-corner-neighbors-verification-v1',
              'source_sha256': hashlib.sha256(source_path.read_bytes()).hexdigest(),
              'models': models, 'recorded_rank_three_neighbors': replayed,
              'all_checks_passed': True,
              'scope': 'Local proof controls only. The constructed families have a '
                       'proper nonpeelable minor and are not forbidden witnesses. '
                       'The theorem is symbolic; these checks are not a census proof.'}
    (ROOT / 'damp_rank_three_corner_neighbors_verification.json').write_text(
        json.dumps(report, indent=2) + '\n')
    print({'all_checks_passed': True, 'models': [
        {k: r[k] for k in ['has_companion', 'order', 'dimension', 'upper_corners',
                           'corners_after_first_old_deletion']} for r in models],
           'recorded_rank_three_neighbors_checked': len(replayed)})


if __name__ == '__main__':
    main()
