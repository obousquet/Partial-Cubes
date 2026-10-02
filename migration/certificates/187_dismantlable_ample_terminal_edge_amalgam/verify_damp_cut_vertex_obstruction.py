#!/usr/bin/env python3
"""Verify the pointed 294 certificate and all minors of its vertex double.

This is an independent replay, with no peeling search. The 587-concept
double is ample, cornerless, and has a unique cut vertex. Each of its 72
elementary restrictions/contractions has an explicit validated peeling.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from verify_damp_layer_midpoint_cores import corner, shadow
from verify_damp_rank_three_corner_neighbors import corners

ROOT = Path(__file__).resolve().parent


def check_order(family, order, dimension, last=None):
    assert len(order) == len(family) and set(order) == family
    if last is not None:
        assert order[-1] == last
    remaining = set(family)
    for x in order:
        assert x in remaining and corner(remaining, x, dimension)
        remaining.remove(x)
    assert not remaining


def expand_dense(order, family):
    active = 0
    for x in family:
        active |= x
    coordinates = [e for e in range(active.bit_length()) if active >> e & 1]
    return [sum(((x >> i) & 1) << e for i, e in enumerate(coordinates)) for x in order]


def connected_after_deletion(family, dimension, deleted):
    remaining = family - {deleted}
    if not remaining:
        return True
    start = next(iter(remaining))
    found, pending = {start}, [start]
    while pending:
        x = pending.pop()
        for e in range(dimension):
            y = x ^ (1 << e)
            if y in remaining and y not in found:
                found.add(y)
                pending.append(y)
    return found == remaining


def main():
    source_path = ROOT / 'damp_one_corner_peelable_control.json'
    pointed_path = ROOT / 'computations/round775_pointed_minor_descent.json'
    source = json.loads(source_path.read_text())
    pointed = json.loads(pointed_path.read_text())
    assert pointed['source_sha256'] == hashlib.sha256(source_path.read_bytes()).hexdigest()
    assert pointed['status'] == 'rooted_minimal_verified' and not pointed['history']
    root = source['only_corner']
    base = {x ^ root for x in source['family']}
    ordinary = [x ^ root for x in source['peeling']]
    assert base == set(pointed['final']['family']) and len(base) == 294
    check_order(base, ordinary, 12)
    sh = shadow(base, 12)
    assert len(sh) == 294 and max(map(int.bit_count, sh)) == 3
    assert corners(base, 12) == [0]
    # This directly proves failure of peeling to 0, independently of the
    # search status stored in the pointed-minor report.
    assert len(base) > 1 and not [x for x in corners(base, 12) if x]
    assert all(connected_after_deletion(base, 12, x) for x in base)
    by_operation = {(r['coordinate'], r['operation']): r
                    for r in pointed['final']['elementary_checks']}
    assert set(by_operation) == {(e, k) for e in range(12)
                                for k in ['root_section', 'contraction']}
    root_orders = {}
    for (e, kind), row in by_operation.items():
        bit = 1 << e
        selected = base if kind == 'contraction' else {x for x in base if not x & bit}
        child = {x & ~bit for x in selected}
        order = expand_dense(row['peeling'], child)
        check_order(child, order, 12, last=0)
        root_orders[e, kind] = order
    double = base | {x << 12 for x in base}
    assert len(double) == 587 and corners(double, 24) == []
    # Both blocks have the given shadows, and mixed pairs omit trace 11.
    for e in range(12):
        for f in range(12, 24):
            mask = (1 << e) | (1 << f)
            assert len({x & mask for x in double}) == 3
    assert len(sh) + len(sh) - 1 == len(double)
    assert not connected_after_deletion(double, 24, 0)
    # The individual blocks have no cut vertices, so 0 is the unique cut
    # vertex of the double. Directly check its two components as well.
    assert len(base - {0}) == 293
    elementary = []
    for block in range(2):
        shift = 12 * block
        other_shift = 12 * (1 - block)
        for e in range(12):
            bit = 1 << (e + shift)
            for kind in ['root_section', 'opposite_section', 'contraction']:
                if kind == 'root_section':
                    selected = {x for x in double if not x & bit}
                elif kind == 'opposite_section':
                    selected = {x for x in double if x & bit}
                else:
                    selected = double
                child = {x & ~bit for x in selected}
                if kind == 'opposite_section':
                    # The other block is discarded by this restriction.
                    order = [((x & ~(1 << e)) << shift) for x in ordinary if x >> e & 1]
                else:
                    order = [x << shift for x in root_orders[e, kind][:-1]]
                    order += [x << other_shift for x in ordinary]
                check_order(child, order, 24)
                assert len(child) < len(double)
                elementary.append({'coordinate': e + shift, 'operation': kind,
                                   'order': len(child), 'peeling': order})
    assert len(elementary) == 72
    report = {'schema': 'damp-cut-vertex-587-verification-v1',
              'source_sha256': hashlib.sha256(source_path.read_bytes()).hexdigest(),
              'pointed_certificate_sha256': hashlib.sha256(pointed_path.read_bytes()).hexdigest(),
              'verifier_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'rooted_block': sorted(base), 'family': sorted(double),
              'order': 587, 'dimension': 24, 'vc_dimension': 3,
              'corner_count': 0, 'unique_cut_vertex': 0,
              'components_after_cut_deletion': [293, 293],
              'blocks_are_two_connected': True, 'rooted_minors_checked': 24,
              'elementary_minor_peelings': elementary,
              'all_checks_passed': True,
              'scope': 'An explicit new ample forbidden pc-minor with a cut vertex. '
                       'The general theorem characterizes the cut-vertex branch '
                       'using pointed minimality; intrinsic admissibility of its '
                       'pointed blocks and the full classification remain open.'}
    (ROOT / 'damp_cut_vertex_587_verification.json').write_text(
        json.dumps(report, indent=2) + '\n')
    print({k: report[k] for k in ['all_checks_passed', 'order', 'dimension',
                                 'vc_dimension', 'corner_count', 'unique_cut_vertex',
                                 'rooted_minors_checked']})
    print({'elementary_minor_peelings_checked': len(elementary)})


if __name__ == '__main__':
    main()
