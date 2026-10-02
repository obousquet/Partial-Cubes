#!/usr/bin/env python3
"""Fixed controls for one-sided midpoint pruning and nonempty-core wings.

No search is performed. The full cores and the relative pruning systems
keep their original unique-midpoint constraints fixed. The uniform wing
bounds use the manuscript proof and the relative Yang forest/eight-facet
theorems, not extrapolation from these controls.
"""
from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from pathlib import Path

from verify_damp_layer_midpoint_cores import corner, midpoint_core, shadow

ROOT = Path(__file__).resolve().parent


def relative_core(upper, lower, dimension):
    assert lower <= upper
    triples = midpoint_core(upper, dimension)[2]
    by_root = defaultdict(list)
    for w, u, v in sorted(triples):
        by_root[w].append((u, v))
    remaining = upper - lower
    pruning = []
    while remaining:
        present = lower | remaining
        removable = sorted(w for w in remaining
                            if not any(u in present and v in present
                                       for u, v in by_root[w]))
        if not removable:
            break
        for w in removable:
            present = lower | remaining
            assert not any(u in present and v in present for u, v in by_root[w])
            remaining.remove(w)
            pruning.append(w)
    present = lower | remaining
    witnesses = [[w, *next((u, v) for u, v in by_root[w]
                           if u in present and v in present)]
                 for w in sorted(remaining)]
    return remaining, {'relative_core': sorted(remaining),
                       'fixed_boundary': sorted(lower),
                       'pruning_order': pruning, 'witnesses': witnesses}


def cycle_witness(family, dimension):
    parent = {}
    depth = {}

    def visit(v, previous):
        parent[v] = previous
        depth[v] = 0 if previous is None else depth[previous] + 1
        for u in sorted(v ^ (1 << e) for e in range(dimension)
                        if v ^ (1 << e) in family):
            if u == previous:
                continue
            if u not in depth:
                found = visit(u, v)
                if found:
                    return found
            elif depth[u] < depth[v]:
                path = [v]
                while path[-1] != u:
                    path.append(parent[path[-1]])
                return path
        return None

    for v in sorted(family):
        if v not in depth:
            cycle = visit(v, None)
            if cycle:
                assert len(cycle) >= 4 and len(set(cycle)) == len(cycle)
                assert all((cycle[i] ^ cycle[(i + 1) % len(cycle)]).bit_count() == 1
                           for i in range(len(cycle)))
                return cycle
    return None


def check_splits(family, dimension):
    whole_core = midpoint_core(family, dimension)[0]
    rows = []
    for e in range(dimension):
        bit = 1 << e
        sides = [{v & ~bit for v in family if (v >> e & 1) == i}
                 for i in [0, 1]]
        assert all(sides)
        common = sides[0] & sides[1]
        assert common
        assert all(not midpoint_core(side, dimension)[0] for side in sides)
        checks = []
        for i, side in enumerate(sides):
            wing = side - common
            relative, certificate = relative_core(side, common, dimension)
            projection = {v & ~bit for v in whole_core if (v >> e & 1) == i}
            assert projection & wing <= relative
            cycle = cycle_witness(wing, dimension)
            if whole_core:
                assert relative and cycle and len(wing) >= 9
            checks.append({'wing_order': len(wing), 'cycle': cycle, **certificate})
        if any(not c['relative_core'] for c in checks):
            assert not whole_core
        rows.append({'coordinate': e, 'overlap_order': len(common), 'sides': checks})
    return whole_core, rows


def verify_peeling(family, order, dimension):
    remaining = set(family)
    for v in order:
        assert v in remaining and corner(remaining, v, dimension)
        remaining.remove(v)
    assert not remaining


def main():
    hall_path = ROOT / 'computations/hall_corner_exchange/source_H.json'
    source = json.loads(hall_path.read_text())
    bits = int(source['source_family_hex'], 16)
    hall = {v for v in range(4096) if bits >> v & 1}
    assert len(hall) == len(shadow(hall, 12)) == 299
    assert not any(corner(hall, v, 12) for v in hall)
    extended = hall | {373}
    assert len(shadow(extended, 12)) == 300
    assert [v for v in sorted(extended) if corner(extended, v, 12)] == [373, 1011]
    cases = []
    for name, family, expected in [('Hall', hall, 158), ('Hall_plus_373', extended, 147)]:
        core, splits = check_splits(family, 12)
        assert len(core) == expected
        cases.append({'name': name, 'order': len(family), 'core_order': len(core),
                      'splits': splits})

    empty_path = ROOT / 'damp_empty_core_counterexample.json'
    empty = set(json.loads(empty_path.read_text())['counterexample'])
    core, splits = check_splits(empty, 12)
    assert not core
    assert [len(c['relative_core']) for c in splits[0]['sides']] == [54, 17]
    cases.append({'name': 'empty_core_292', 'order': len(empty), 'core_order': 0,
                  'splits': splits})

    # A genuinely one-sided control: the zero wing cannot start peeling,
    # whereas the one section already equals the overlap.
    lower_section = {v >> 1 for v in hall if not v & 1}
    upper_section = {v >> 1 for v in hall if v & 1}
    common = lower_section & upper_section
    peel_path = ROOT / 'computations/round765_peripheral_section_peelings.json'
    orders = json.loads(peel_path.read_text())
    verify_peeling(lower_section, orders['section'], 11)
    verify_peeling(common, orders['overlap'], 11)
    assert len(lower_section) == 157 and len(common) == 67
    section_corners = [v for v in sorted(lower_section) if corner(lower_section, v, 11)]
    assert len(section_corners) == 9 and set(section_corners) <= common
    left_relative = relative_core(lower_section, common, 11)[0]
    assert len(left_relative) == 52
    control = lower_section | {v | 2048 for v in common}
    assert len(control) == len(shadow(control, 12)) == 224
    assert not midpoint_core(control, 12)[0]
    order = [v | 2048 for v in orders['overlap']] + orders['section']
    verify_peeling(control, order, 12)
    result = {
        'scope': __doc__.strip(),
        'source_hashes': {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                          for p in [hall_path, empty_path, peel_path]},
        'cases': cases,
        'one_sided_control': {'order': 224, 'section_orders': [157, 67],
                              'relative_core_orders': [52, 0], 'core_order': 0,
                              'all_zero_section_corners_in_overlap': section_corners,
                              'peeling': order},
        'all_checks_passed': True,
    }
    (ROOT / 'damp_relative_midpoint_wings_verification.json').write_text(
        json.dumps(result, indent=2) + '\n')
    print(json.dumps({'all_checks_passed': True,
                      'nonempty_core_coordinate_checks': 24,
                      'cyclic_wings_checked': 48,
                      'empty_example_coordinate_zero_relative_cores': [54, 17],
                      'one_sided_control_order': 224,
                      'one_sided_control_relative_core_orders': [52, 0]}, indent=2))


if __name__ == '__main__':
    main()
