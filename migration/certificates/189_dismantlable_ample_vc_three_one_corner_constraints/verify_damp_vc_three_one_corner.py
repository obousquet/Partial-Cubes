#!/usr/bin/env python3
"""Fixed controls for the general VC-three one-corner exclusion.

Two incomparable rank-three neighbor supports either stop before the
old cube antipode, or the cube at that antipode blocks the other neighbor's
second deletion. Check both mechanisms on nonminimal ample controls.
The uniform theorem is proved symbolically in the manuscript.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from verify_damp_layer_midpoint_cores import shadow, subsets
from verify_damp_rank_three_corner_neighbors import (
    attached_family, check_block_shadow, check_local, corners, directions,
)
from verify_damp_vc_three_companion import negative_certificate

ROOT = Path(__file__).resolve().parent


def unused_direction_control(block, block_shadow):
    """Both three-deletion prefixes can survive with one unused direction.

    The fixed VC-four control is still negative because the attached block
    is a proper restriction. Only these two prefixes are checked here.
    """
    local = set(range(16)) | {17, 18, 19, 34, 36, 38} | set(range(72, 80))
    local_shadow = shadow(local, 7)
    assert len(local_shadow) == len(local) == 30
    assert max(map(int.bit_count, local_shadow)) == 4
    anchors = [x for x in corners(local, 7) if x]
    source, dimension = attached_family(local, 7, anchors, block)
    assert corners(source, dimension) == [0]
    check_block_shadow(source, local, 7, anchors, block, block_shadow)
    rows = []
    for tip, J in [(16, 3), (32, 6)]:
        full = source | {tip}
        check_block_shadow(full, local | {tip}, 7, anchors, block, block_shadow)
        assert shadow(local | {tip}, 7) - local_shadow == {tip | J}
        case = check_local(source, 0, tip, dimension)
        assert case['corners_after_first_old_deletion'] == [J]
        after_two = full - {tip ^ J, J}
        third = J | 8
        assert corners(after_two, dimension) == [third]
        assert directions(after_two, third, dimension) == 71
        assert 8 == 15 & ~(3 | 6)
        forced_cube = {third ^ U for U in subsets(71)}
        assert forced_cube <= source
        # This cube has the unused direction set to one; both old second
        # vertices have it set to zero and escape its forced outside edges.
        assert not {3, 6} & forced_cube
        mask = 4095 << 7
        assert {(x >> 7) & 4095 for x in full if x & ~mask == anchors[0]} == block
        rows.append({**case, 'third_deletion': third, 'unused_direction_mask': 8,
                     'third_cube_support': 71, 'has_proper_cornerless_291_minor': True})
    return {'vc_dimension': 4, 'order': len(source), 'dimension': dimension,
            'local_family': sorted(local), 'attachment_anchors': anchors,
            'neighbors': rows, 'scope': 'Two valid three-deletion prefixes, '
            'not two resolving additions; each extension has a proper nonpeelable minor.'}


def main():
    source_path = ROOT / 'damp_fixed_291_repair_creation.json'
    block = set(json.loads(source_path.read_text())['base'])
    block_shadow = shadow(block, 12)
    assert len(block) == len(block_shadow) == 291 and 0 in block
    assert max(map(int.bit_count, block_shadow)) == 3
    assert not corners(block, 12)
    models = []
    for third_cube in (False, True):
        local = set(range(8)) | {9, 10, 11, 18, 20, 22}
        local_dimension = 5
        if third_cube:
            local |= {36, 37, 38, 39}
            local_dimension = 6
        local_shadow = shadow(local, local_dimension)
        assert len(local_shadow) == len(local)
        assert max(map(int.bit_count, local_shadow)) == 3
        anchors = [x for x in corners(local, local_dimension) if x != 0]
        source, dimension = attached_family(local, local_dimension, anchors, block)
        assert corners(source, dimension) == [0]
        # The selected block remains a proper coordinate restriction.
        mask = 4095 << local_dimension
        assert {(x >> local_dimension) & 4095 for x in source
                if x & ~mask == anchors[0]} == block
        ampleness = check_block_shadow(source, local, local_dimension, anchors,
                                      block, block_shadow)
        rows = []
        for tip, J in [(8, 3), (16, 6)]:
            full = source | {tip}
            gain = shadow(local | {tip}, local_dimension) - local_shadow
            assert gain == {tip | J}
            check_block_shadow(full, local | {tip}, local_dimension, anchors,
                               block, block_shadow)
            case = check_local(source, 0, tip, dimension)
            assert case['upper_corners'] == [tip, tip ^ J]
            certificate = negative_certificate(full, dimension, limit=32)
            after_two = full - {tip ^ J, J}
            if third_cube and tip == 16:
                assert case['corners_after_first_old_deletion'] == []
                assert directions(source, J, dimension) & ~7 == 48
                third = None
            else:
                assert case['corners_after_first_old_deletion'] == [J]
                third = corners(after_two, dimension)
                assert third == ([7] if third_cube else [])
            if third_cube and tip == 8:
                assert directions(source, 7, dimension) & ~7 == 32
                assert directions(after_two, 7, dimension) == 35
                forced_cube = {7 ^ U for U in subsets(35)}
                assert forced_cube <= source
                # The second neighbor's opposite old vertex acquires an
                # outside neighbor besides its prescribed mate 22.
                assert 6 in forced_cube and 38 in forced_cube
                assert 22 in source and 38 != 22
                case['forced_third_cube'] = sorted(forced_cube)
                case['other_second_vertex'] = 6
                case['blocking_outside_neighbor'] = 38
            rows.append({**case, 'corners_after_two_forced_deletions': third,
                         'negative_certificate': certificate})
        models.append({'has_third_cube': third_cube, 'order': len(source),
                       'dimension': dimension, 'vc_dimension': 3,
                       'local_family': sorted(local), 'attachment_anchors': anchors,
                       'ampleness': ampleness, 'neighbors': rows,
                       'has_proper_cornerless_291_minor': True})
    unused = unused_direction_control(block, block_shadow)
    report = {'schema': 'damp-vc-three-one-corner-verification-v1',
              'source_sha256': hashlib.sha256(source_path.read_bytes()).hexdigest(),
              'verifier_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'models': models, 'unused_direction_control': unused,
              'all_checks_passed': True,
              'scope': 'Fixed nonminimal controls, not new forbidden witnesses '
                       'or a classification by finite enumeration. The first '
                       'tests a blocked third deletion; the second tests the '
                       'conflict forced by an available third deletion.'}
    (ROOT / 'damp_vc_three_one_corner_verification.json').write_text(
        json.dumps(report, indent=2) + '\n')
    print({'all_checks_passed': True, 'models': [
        {k: m[k] for k in ['has_third_cube', 'order', 'dimension']} |
        {'certificate_states': [len(r['negative_certificate']) for r in m['neighbors']]}
        for m in models], 'unused_direction_control':
        {k: unused[k] for k in ['vc_dimension', 'order', 'dimension']}})


if __name__ == '__main__':
    main()
