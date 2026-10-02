#!/usr/bin/env python3
"""Check the fixed controls for the layerwise unique-midpoint core formula.

No search is performed. The symbolic layer-core and terminal tower theorems
give the uniform conclusions; these controls independently enumerate rooted
pairs, validate stored repair peelings, and check complete static prunings.
"""
from __future__ import annotations

import hashlib
import json
from collections import defaultdict, deque
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def subsets(mask):
    sub = mask
    while True:
        yield sub
        if not sub:
            return
        sub = (sub - 1) & mask


def corner(family, word, dimension):
    directions = sum(1 << e for e in range(dimension)
                     if word ^ (1 << e) in family)
    return all(word ^ sub in family for sub in subsets(directions))


def shadow(family, dimension):
    return {mask for mask in range(1 << dimension)
            if len({v & mask for v in family}) == 1 << mask.bit_count()}


def midpoint_core(family, dimension):
    triples = []
    by_root = defaultdict(int)
    by_head = defaultdict(list)
    for w in sorted(family):
        neighbors = sorted(w ^ (1 << e) for e in range(dimension)
                           if w ^ (1 << e) in family)
        for i, u in enumerate(neighbors):
            for v in neighbors[i + 1:]:
                if u ^ v ^ w not in family:
                    k = len(triples)
                    triples.append((w, u, v))
                    by_root[w] += 1
                    by_head[u].append(k)
                    by_head[v].append(k)
    remaining = set(family)
    active = [True] * len(triples)
    pending = deque(w for w in sorted(family) if not by_root[w])
    order = []
    while pending:
        w = pending.popleft()
        if w not in remaining:
            continue
        assert not by_root[w]
        remaining.remove(w)
        order.append(w)
        for k in by_head[w]:
            if not active[k]:
                continue
            active[k] = False
            root = triples[k][0]
            by_root[root] -= 1
            if root in remaining and not by_root[root]:
                pending.append(root)
    # Recheck the pruning independently against the original triples.
    by_root_pairs = defaultdict(list)
    for w, u, v in triples:
        by_root_pairs[w].append((u, v))
    checked = set(family)
    for w in order:
        assert not any(u in checked and v in checked for u, v in by_root_pairs[w])
        checked.remove(w)
    assert checked == remaining
    witnesses = []
    for w in sorted(remaining):
        pair = next((u, v) for u, v in by_root_pairs[w]
                    if u in remaining and v in remaining)
        witnesses.append([w, *pair])
    return remaining, {
        'constraint_count': len(triples), 'core': sorted(remaining),
        'static_pruning_order': order, 'core_witnesses': witnesses,
    }, set(triples)


def assemble(seed, rows, old_dimension, layer_dimension):
    return ({s | (t << old_dimension) for s in seed
             for t in range(1 << layer_dimension)}
            | {a | (t << old_dimension) for a, row in rows.items() for t in row})


def check_formula(seed, rows, old_dimension, layer_dimension):
    for row in rows.values():
        assert not midpoint_core(row, layer_dimension)[0]
    family = assemble(seed, rows, old_dimension, layer_dimension)
    actual, certificate, triples = midpoint_core(family, old_dimension + layer_dimension)
    expected = set()
    layer_sizes = []
    for t in range(1 << layer_dimension):
        section = seed | {a for a, row in rows.items() if t in row}
        k = midpoint_core(section, old_dimension)[0]
        expected.update(v | (t << old_dimension) for v in k)
        layer_sizes.append(len(k))
    assert actual == expected
    assert len(shadow(family, old_dimension + layer_dimension)) == len(family)
    return family, {
        'order': len(family), 'dimension': old_dimension + layer_dimension,
        'rows': {str(a): sorted(row) for a, row in rows.items()},
        'layer_core_orders': layer_sizes, **certificate,
    }, triples


def main():
    source_path = ROOT / 'damp_fixed_291_repair_creation.json'
    source = json.loads(source_path.read_text())
    seed = set(source['base'])
    model = ROOT / 'damp_hall_carrier_defect_min8.json'
    assert hashlib.sha256(model.read_bytes()).hexdigest() == source['source_model_sha256']
    assert len(seed) == len(shadow(seed, 12)) == 291
    assert not any(corner(seed, v, 12) for v in seed)
    seed_core = midpoint_core(seed, 12)[0]
    assert len(seed_core) == 76
    repairs = {}
    for a, support in [(91, 102), (157, 390)]:
        row = next(r for r in source['extensions'] if r['tip'] == a)
        assert row['support'] == support and row['status'] == 'peelable'
        assert shadow(seed | {a}, 12) - shadow(seed, 12) == {support}
        remaining = seed | {a}
        for v in row['peeling']:
            assert v in remaining and corner(remaining, v, 12)
            remaining.remove(v)
        assert not remaining
        repairs[str(a)] = {'support': support, 'peeling': row['peeling']}

    controls = []
    _, row, _ = check_formula(seed, {}, 12, 1)
    assert len(row['core']) == 152
    controls.append({'name': 'unrepaired_product', **row})
    for r in [1, 2]:
        full = (1 << r) - 1
        rows = {91: set(range(1 << r)) - {0},
                157: set(range(1 << r)) - {full}}
        family, row, triples = check_formula(seed, rows, 12, r)
        assert len(family) == 293 * (1 << r) - 2 and not row['core']
        corners = sorted(v for v in family if corner(family, v, 12 + r))
        assert len(corners) == 2 * r
        if r == 1:
            # Independent endpoint-pair enumeration checks completeness.
            direct = set()
            words = sorted(family)
            for i, u in enumerate(words):
                for v in words[i + 1:]:
                    delta = u ^ v
                    if delta.bit_count() != 2:
                        continue
                    bit = delta & -delta
                    mids = [w for w in [u ^ bit, v ^ bit] if w in family]
                    if len(mids) == 1:
                        direct.add((mids[0], u, v))
            assert direct == triples
        deletions = []
        for v in corners:
            a, t = v & 4095, v >> 12
            reduced = {b: row - {t} if b == a else set(row)
                       for b, row in rows.items()}
            reduced = {b: row for b, row in reduced.items() if row}
            child, child_row, _ = check_formula(seed, reduced, 12, r)
            assert child == family - {v}
            assert len(child_row['core']) == (76 if r == 1 else 0)
            deletions.append({'deleted': v, 'core_order': len(child_row['core'])})
        controls.append({'name': f'antipodal_tower_{r}', 'corners': corners,
                         'corner_deletions': deletions, **row})

    # The empty-row-core hypothesis matters: all horizontal sections here
    # are vertices or edges, while the single occurrence row has a core.
    family = assemble({0}, {1: seed}, 1, 12)
    actual, row, _ = midpoint_core(family, 13)
    assert actual == {1 | (v << 1) for v in seed_core}
    assert len(shadow(family, 13)) == len(family) == 4387
    controls.append({'name': 'nonempty_row_core', 'order': len(family),
                     'all_horizontal_sections_peel': True, **row})
    result = {
        'scope': __doc__.strip(),
        'source_file_sha256': hashlib.sha256(source_path.read_bytes()).hexdigest(),
        'source_model_sha256': source['source_model_sha256'],
        'seed_order': 291, 'seed_core_order': 76, 'checked_repairs': repairs,
        'controls': controls, 'all_checks_passed': True,
        'uniform_forbidden_and_terminal_status': 'proved by manuscript assembly theorems',
    }
    (ROOT / 'damp_layer_midpoint_cores_verification.json').write_text(
        json.dumps(result, indent=2) + '\n')
    print(json.dumps({'all_checks_passed': True,
                      'controls': [{'name': r['name'], 'order': r['order'],
                                    'core_order': len(r['core'])} for r in controls]}, indent=2))


if __name__ == '__main__':
    main()
