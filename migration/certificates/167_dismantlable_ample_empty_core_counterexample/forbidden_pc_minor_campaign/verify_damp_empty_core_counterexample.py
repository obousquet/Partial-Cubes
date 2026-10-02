#!/usr/bin/env python3
"""Check the fixed 292-concept counterexample to the unique-midpoint criterion.

This is a certificate checker, with no peeling search. It recomputes the
complete shattered complexes, corner sets, static rooted-pair cores, and
all 72 supplied elementary-minor peelings of the example and predecessor.
The empty-core pruning order and one rooted-pair witness per member of
the activated core are written to the verification report.
"""
from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def subsets(mask):
    current = mask
    while True:
        yield current
        if not current:
            break
        current = (current - 1) & mask


def incident(family, v):
    return sum(1 << e for e in range(12) if v ^ (1 << e) in family)


def is_corner(family, v):
    return all(v ^ mask in family for mask in subsets(incident(family, v)))


def shattered(family):
    return {mask for mask in range(4096)
            if len({v & mask for v in family}) == 1 << mask.bit_count()}


def static_core(family):
    """Enumerate distance-two endpoint pairs; keep original constraints fixed."""
    constraints = []
    words = sorted(family)
    for i, u in enumerate(words):
        for v in words[i + 1:]:
            difference = u ^ v
            if difference.bit_count() != 2:
                continue
            bit = difference & -difference
            midpoints = [w for w in [u ^ bit, v ^ bit] if w in family]
            if len(midpoints) == 1:
                constraints.append((midpoints[0], u, v))
    by_root = defaultdict(list)
    for w, u, v in constraints:
        by_root[w].append((u, v))
    remaining = set(family)
    pruning = []
    while remaining:
        removable = sorted(w for w in remaining
                           if not any(u in remaining and v in remaining
                                      for u, v in by_root[w]))
        if not removable:
            break
        for w in removable:
            assert not any(u in remaining and v in remaining for u, v in by_root[w])
            remaining.remove(w)
            pruning.append(w)
    witnesses = [[w, *next((u, v) for u, v in by_root[w]
                           if u in remaining and v in remaining)]
                 for w in sorted(remaining)]
    return constraints, remaining, pruning, witnesses


def verify_order(family, order):
    remaining = set(family)
    for v in order:
        assert v in remaining and is_corner(remaining, v)
        remaining.remove(v)
    assert not remaining


def main():
    path = ROOT / 'damp_empty_core_counterexample.json'
    data = json.loads(path.read_text())
    source_path = ROOT / data['source_file']
    assert hashlib.sha256(source_path.read_bytes()).hexdigest() == data['source_file_sha256']
    source = json.loads(source_path.read_text())
    model_path = ROOT / 'damp_hall_carrier_defect_min8.json'
    assert hashlib.sha256(model_path.read_bytes()).hexdigest() == data['source_model_sha256']
    base = set(source['base'])
    predecessor = (base - set(data['removed_from_original'])) | set(data['added_to_original'])
    assert predecessor == set(data['predecessor']) and len(predecessor) == 291
    corner = data['corner']
    family = predecessor | {corner}
    assert family == set(data['counterexample']) and len(family) == 292
    assert data['removed_from_original'] == [96, 126, 394, 1922, 1934]
    assert data['added_to_original'] == [553, 763, 2968, 4008, 4026]
    shadows = [shattered(x) for x in [predecessor, family]]
    assert [len(s) for s in shadows] == [291, 292]
    assert not any(is_corner(predecessor, v) for v in predecessor)
    assert [v for v in sorted(family) if is_corner(family, v)] == [3959]
    acquired = incident(family, corner)
    assert acquired == 2192 and shadows[1] - shadows[0] == {acquired}
    constraints, core, pruning, witnesses = static_core(family)
    old_constraints, old_core, old_pruning, old_witnesses = static_core(predecessor)
    assert not core and len(pruning) == 292
    assert len(old_core) == 76 and len(old_witnesses) == 76
    # A second formulation enumerates rooted neighbor pairs directly.
    direct = set()
    for w in family:
        neighbors = sorted(v for v in family if (v ^ w).bit_count() == 1)
        for i, u in enumerate(neighbors):
            for v in neighbors[i + 1:]:
                if u ^ v ^ w not in family:
                    direct.add((w, u, v))
    assert direct == set(constraints)
    minor_counts = {}
    for record in data['elementary_minor_peelings']:
        current = set(record['family'])
        assert current == (family if record['label'] == 'counterexample' else predecessor)
        assert len(record['minor_checks']) == 36
        seen = set()
        for row in record['minor_checks']:
            e, side = row['coordinate'], row['side']
            seen.add((e, side))
            image = {v & ~(1 << e) for v in current
                     if side is None or (v >> e & 1) == side}
            assert len(image) == row['order'] and row['status'] == 'peelable'
            verify_order(image, row['peeling'])
        assert seen == {(e, side) for e in range(12) for side in [None, 0, 1]}
        assert all({v >> e & 1 for v in current} == {0, 1} for e in range(12))
        minor_counts[record['label']] = len(seen)
    assert minor_counts == {'counterexample': 36, 'cornerless_predecessor': 36}
    # Discovery provenance is checked without repeating the search.
    initial = base | set(data['discovery']['maximal_compatible_tips'])
    assert len(shattered(initial)) == len(initial) == 305
    for v in data['discovery']['corner_deletion_prefix']:
        assert v in initial and is_corner(initial, v)
        initial.remove(v)
    assert initial == family
    result = {
        'schema': 'damp-empty-core-counterexample-verification-v1',
        'witness_file_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
        'dimension': 12, 'order': len(family), 'shattered_order': len(shadows[1]),
        'only_corner': corner, 'corner_support': [e for e in range(12) if acquired >> e & 1],
        'unique_midpoint_constraint_count': len(constraints),
        'empty_core_pruning_order': pruning,
        'predecessor_order': len(predecessor), 'predecessor_cornerless': True,
        'predecessor_core_order': len(old_core),
        'predecessor_core': sorted(old_core),
        'predecessor_core_witnesses': old_witnesses,
        'elementary_minor_peelings_checked': minor_counts,
        'counterexample_is_forbidden': True,
        'predecessor_is_forbidden': True,
        'all_checks_passed': True,
    }
    (ROOT / 'damp_empty_core_counterexample_verification.json').write_text(
        json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ['empty_core_pruning_order', 'predecessor_core',
                                   'predecessor_core_witnesses']}, indent=2))


if __name__ == '__main__':
    main()
