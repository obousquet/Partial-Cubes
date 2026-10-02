#!/usr/bin/env python3
"""Fixed controls for the VC-three companion exclusion.

Attach the stored cornerless 291 family at every local corner. These
controls have a proper nonpeelable minor and are not forbidden witnesses.
Enumerate at most 32 deletion states per control, then independently
validate the resulting complete negative certificate. No seed search.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from verify_damp_layer_midpoint_cores import shadow, subsets
from verify_damp_rank_three_corner_neighbors import (
    attached_family, check_block_shadow, check_state_certificate,
    corners, directions,
)

ROOT = Path(__file__).resolve().parent


def negative_certificate(initial, dimension, limit=32):
    pending = [frozenset()]
    seen = set(pending)
    rows = []
    while pending:
        removed = pending.pop()
        family = initial - removed
        assert family, 'A complete peeling contradicts the control contract'
        cs = corners(family, dimension)
        rows.append({'removed': sorted(removed), 'corners': cs})
        for c in cs:
            nxt = removed | {c}
            if nxt not in seen:
                seen.add(nxt)
                assert len(seen) <= limit, 'Fixed diagnostic state cap exceeded'
                pending.append(nxt)
    rows.sort(key=lambda row: (len(row['removed']), row['removed']))
    check_state_certificate(initial, dimension, rows)
    return rows


def main():
    source_path = ROOT / 'damp_fixed_291_repair_creation.json'
    block = set(json.loads(source_path.read_text())['base'])
    block_shadow = shadow(block, 12)
    assert 0 in block and len(block_shadow) == len(block) == 291
    assert max(map(int.bit_count, block_shadow)) == 3
    assert not corners(block, 12)
    models = []
    for escape in (False, True):
        local = {t | s for t in (0, 8, 16) for s in range(1, 8)}
        local_dimension = 5
        if escape:
            local |= {36, 37, 38, 39, 44, 45, 46, 47}
            local_dimension = 6
        local_shadow = shadow(local, local_dimension)
        assert len(local_shadow) == len(local)
        vc = max(map(int.bit_count, local_shadow))
        assert vc == 3 + int(escape)
        anchors = corners(local, local_dimension)
        base, dimension = attached_family(local, local_dimension, anchors, block)
        assert not corners(base, dimension)
        checks = {}
        for tips in ((), (0,), (8,), (16,), (0, 8)):
            name = ','.join(map(str, tips)) or 'base'
            small, large = local | set(tips), base | set(tips)
            small_shadow = shadow(small, local_dimension)
            gained = small_shadow - local_shadow
            assert gained == (set() if not tips else {7, 15} if len(tips) == 2 else {7})
            checks[name] = check_block_shadow(
                large, small, local_dimension, anchors, block, block_shadow)
            checks[name]['gained_supports'] = sorted(gained)
            if len(tips) == 1:
                assert directions(base, tips[0], dimension) == 7
                assert all(tips[0] ^ s in base for s in subsets(7) if s)
                assert corners(large, dimension) == list(tips)
        upper = base | {0, 8}
        rows = negative_certificate(upper, dimension)
        first = corners(upper, dimension)
        if not escape:
            assert anchors == [9, 10, 12, 17, 18, 20]
            assert len(base) == 1761 and dimension == 77
            assert first == [0, 8, 11, 13, 14, 15] and len(rows) == 15
            far = {11, 13, 14, 15}
            for row in rows:
                removed = set(row['removed'])
                if removed & far:
                    assert removed <= far
                    remaining = upper - removed
                    for x in row['corners']:
                        assert x in far
                        assert not directions(upper, x, dimension) & ~15
                    # The only new four-support is lost at the first deletion.
                    local_remaining = (local | {0, 8}) - removed
                    assert max(map(int.bit_count, shadow(local_remaining, 5))) == 3
        else:
            assert first == [0, 8, 11]
            assert corners(upper - {11}, dimension) == [15]
            support = directions(upper - {11}, 15, dimension)
            assert support == 43 and support.bit_count() == 4
            assert 47 in upper and 47 not in set(range(16))
        # Fix the local word and all other attachment blocks: obtain K.
        shift = local_dimension
        mask = 4095 << shift
        assert {(x >> shift) & 4095 for x in base
                if x & ~mask == anchors[0]} == block
        models.append({'vc_dimension': vc, 'order': len(base),
                       'dimension': dimension, 'local_family': sorted(local),
                       'attachment_anchors': anchors, 'ampleness': checks,
                       'pair_completion_corners': first,
                       'negative_certificate': rows,
                       'has_proper_cornerless_291_minor': True,
                       'outside_neighbor_after_first_deletion': 47 if escape else None})
    report = {'schema': 'damp-vc-three-companion-verification-v1',
              'source_sha256': hashlib.sha256(source_path.read_bytes()).hexdigest(),
              'verifier_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'models': models, 'all_checks_passed': True,
              'scope': 'The theorem has a symbolic proof. Both fixed controls are '
                       'nonminimal. The VC-four control exhibits a local outside '
                       'deletion, not a peelable completion or a counterexample '
                       'to nonpeelability without a dimension bound.'}
    (ROOT / 'damp_vc_three_companion_verification.json').write_text(
        json.dumps(report, indent=2) + '\n')
    print({'all_checks_passed': True, 'models': [
        {k: r[k] for k in ('vc_dimension', 'order', 'dimension', 'attachment_anchors',
                          'pair_completion_corners')} | {'states': len(r['negative_certificate'])}
        for r in models]})


if __name__ == '__main__':
    main()
