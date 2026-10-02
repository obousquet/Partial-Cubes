#!/usr/bin/env python3
"""Reconstruct and check the existing 294-concept one-corner positive class.

No search is performed. The stored full peeling proves that one initial
corner does not imply nonpeelability. Since that corner must be deleted
first, a peeling cannot be required to end at any prescribed concept.
"""
from __future__ import annotations

import hashlib
import json
from itertools import combinations
from pathlib import Path

from verify_damp_layer_midpoint_cores import corner, shadow, subsets
from verify_damp_rank_three_corner_neighbors import corners, directions

ROOT = Path(__file__).resolve().parent


def reconstruct():
    owner = json.loads((ROOT / 'damp_hall_owner_forced_shadow65.json').read_text())['mathematics']
    source = json.loads((ROOT / 'damp_hall_forced_core_relative_defect1_f61_o67.json').read_text())['mathematics']
    hall = json.loads((ROOT / 'computations/hall_corner_exchange/source_H.json').read_text())
    bits = int(hall['source_family_hex'], 16)
    projection = {x >> 1 for x in range(4096) if bits >> x & 1}
    assert len(projection) == 232
    overlap = set(source['rows'][0]['relative_corner_completion_members'])
    red = set(owner['witness']['red_supports'])
    lower, upper = set(overlap), set(overlap)
    for coords in combinations(range(11), 3):
        support = sum(1 << e for e in coords)
        contained = []
        for anchor in {x & ~support for x in projection}:
            cube = {anchor | t for t in subsets(support)}
            if cube <= projection:
                contained.append(cube)
        assert len(contained) == 1
        (lower if support in red else upper).update(contained[0])
    assert lower & upper == overlap and lower | upper == projection
    assert (len(lower), len(upper), len(overlap)) == (89, 205, 62)
    for family in [lower, upper, overlap]:
        assert len(shadow(family, 11)) == len(family)
    return lower | {x | 2048 for x in upper}


def main():
    path = ROOT / 'damp_one_corner_peelable_control.json'
    data = json.loads(path.read_text())
    for name, digest in data['source_files'].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest
    source = reconstruct()
    assert source == set(data['family']) and len(source) == data['order'] == 294
    sh = shadow(source, 12)
    assert len(sh) == len(source) and max(map(int.bit_count, sh)) == 3
    only = data['only_corner']
    assert only == 1894 and corners(source, 12) == [only]
    support = directions(source, only, 12)
    assert support == 2305
    order = data['peeling']
    assert len(order) == len(source) and set(order) == source
    assert order[0] == only and order[-1] != only
    remaining = set(source)
    deletion_supports = []
    for x in order:
        assert x in remaining and corner(remaining, x, 12)
        deletion_supports.append(directions(remaining, x, 12))
        remaining.remove(x)
    assert not remaining
    # Inspect only the nine missing neighbors of the already fixed corner.
    # No ample neighbor exists, so this is not a two-neighbor repair witness.
    adjacent_checks = []
    for e in range(12):
        tip = only ^ (1 << e)
        if tip in source:
            continue
        neighbor_support = directions(source, tip, 12)
        missing = [t for t in subsets(neighbor_support) if t and tip ^ t not in source]
        gains_new_support = neighbor_support not in sh
        assert missing or not gains_new_support
        adjacent_checks.append({'tip': tip, 'incident_support': neighbor_support,
                                'missing_cube_offsets': missing,
                                'support_already_shattered': not gains_new_support})
    assert len(adjacent_checks) == 9
    report = {'schema': 'damp-one-corner-peelable-control-verification-v1',
              'certificate_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
              'verifier_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'order': 294, 'dimension': 12, 'vc_dimension': 3,
              'only_corner': only, 'corner_support': support,
              'corner_cube': sorted(only ^ t for t in subsets(support)),
              'peeling_steps_checked': len(order), 'deletion_supports': deletion_supports,
              'adjacent_missing_words_checked': adjacent_checks,
              'prescribed_last_vertex_1894_is_impossible': True,
              'all_checks_passed': True,
              'scope': 'A positive class, not a new forbidden witness. This refutes '
                       'universal prescribed-endpoint peeling and the implication '
                       'from one corner to nonpeelability. It does not decide the '
                       'VC-four two-resolving-neighbor problem.'}
    (ROOT / 'damp_one_corner_peelable_control_verification.json').write_text(
        json.dumps(report, indent=2) + '\n')
    print({k: report[k] for k in ['all_checks_passed', 'order', 'dimension',
                                  'vc_dimension', 'only_corner', 'peeling_steps_checked']})


if __name__ == '__main__':
    main()
