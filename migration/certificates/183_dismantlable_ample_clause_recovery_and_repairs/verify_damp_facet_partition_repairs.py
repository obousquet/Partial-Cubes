#!/usr/bin/env python3
"""Verify zero-, one-, and two-repair Hall assemblies and nested flattening.

Completeness of the repair list follows from the nonpeelable strong
projection theorem. This checks every possible single-word fibre
completion and supplies explicit peelings for every listed repair.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from verify_damp_terminal_clause_recovery import corner, subsets

ROOT = Path(__file__).resolve().parent


def build(seed, rows, layer_mask):
    family = {x | t for x in seed for t in subsets(layer_mask)}
    for tip, row in rows.items():
        family.update(tip | t for t in row)
    return family


def checked_repair_peeling(seed, rows, layer_mask, tip, old_order, dimension):
    layers = set(subsets(layer_mask))
    missing = layers - rows[tip]
    assert len(missing) == 1
    z = tip | next(iter(missing))
    remaining = build(seed, rows, layer_mask) | {z}
    order = []
    for other, row in rows.items():
        if other == tip:
            continue
        row = set(row)
        while row:
            t = next(t for t in sorted(row) if corner(row, t, dimension))
            word = other | t
            assert word in remaining and corner(remaining, word, dimension)
            remaining.remove(word)
            order.append(word)
            row.remove(t)
    assert set(old_order) == seed | {tip} and len(old_order) == len(seed) + 1
    for word in old_order:
        for t in sorted(layers):
            lifted = word | t
            assert lifted in remaining and corner(remaining, lifted, dimension)
            remaining.remove(lifted)
            order.append(lifted)
    assert not remaining
    return {'added_word': z, 'peeling': order}


def main():
    source_path = ROOT / 'computations/hall_corner_exchange/source_H.json'
    source = json.loads(source_path.read_text())
    bits = int(source['source_family_hex'], 16)
    seed = {x for x in range(1 << 12) if bits >> x & 1}
    assert len(seed) == 299 and source['source_label'] == 'H'
    assert not any(corner(seed, x, 12) for x in seed)
    for support in subsets((1 << 12) - 1):
        if support.bit_count() == 4:
            assert len({x & support for x in seed}) == 15
    old_certificate_path = ROOT / 'damp_terminal_clause_recovery_verification.json'
    old_certificate = json.loads(old_certificate_path.read_text())
    assert old_certificate['source_mathematical_profile_sha256'] == source['mathematical_profile_sha256']
    orders = {r['tip']: r['peeling'] for r in old_certificate['repair_proofs']}
    acquired = set()
    for tip, old_order in orders.items():
        support = sum(1 << e for e in range(12) if tip ^ (1 << e) in seed)
        assert support.bit_count() == 4 and support not in acquired
        assert all(tip ^ sub in seed for sub in subsets(support) if sub)
        assert tip & support not in {x & support for x in seed}
        acquired.add(support)
        remaining = seed | {tip}
        for word in old_order:
            assert word in remaining and corner(remaining, word, 12)
            remaining.remove(word)
        assert not remaining

    a, b, c = 183, 315, 1836
    p, q, r = 1 << 12, 1 << 13, 1 << 14
    pq, pqr = set(subsets(p | q)), set(subsets(p | q | r))
    models = [
        ('two_owners', {a: pq - {0}, b: pq - {p | q}}, p | q, 2),
        ('three_owners_one_full_support',
         {a: pq - {p | q}, b: {p, p | q}, c: {q, p | q}}, p | q, 1),
        ('three_owners_no_full_support',
         {a: {t for t in pqr if not t & p or not t & q},
          b: {t for t in pqr if t & q or not t & r},
          c: {t for t in pqr if t & p or t & r}}, p | q | r, 0),
    ]
    reports = []
    for name, rows, mask, expected in models:
        family = build(seed, rows, mask)
        layers = set(subsets(mask))
        # Enumerate all old words: an addition must complete one of these
        # full layer fibres to change the nonpeelable strong projection.
        completions = []
        for old in range(1 << 12):
            missing = {t for t in layers if old | t not in family}
            if len(missing) == 1:
                completions.append(old | next(iter(missing)))
        expected_tips = [tip for tip, row in rows.items() if len(layers - row) == 1]
        assert len(completions) == len(expected_tips) == expected
        peelings = [checked_repair_peeling(seed, rows, mask, tip, orders[tip], mask.bit_length())
                    for tip in expected_tips]
        assert sorted(completions) == sorted(item['added_word'] for item in peelings)
        reports.append({'name': name, 'order': len(family), 'dimension': mask.bit_length(),
                        'complete_repair_list': sorted(completions), 'peelings': peelings})

    inner = build(seed, {a: {p}, b: {0}}, p)
    nested = build(inner, {a: {q}, b | p: {0}}, q)
    flattened = build(seed, {a: pq - {0}, b: pq - {p | q}}, p | q)
    assert nested == flattened and len(nested) == 1202
    result = {'scope': __doc__.strip(),
              'source_mathematical_profile_sha256': source['mathematical_profile_sha256'],
              'source_file_sha256': hashlib.sha256(source_path.read_bytes()).hexdigest(),
              'seed_repair_certificate_sha256': hashlib.sha256(old_certificate_path.read_bytes()).hexdigest(),
              'models': reports, 'nested_flattening_order': len(nested),
              'all_checks_passed': True}
    output = ROOT / 'damp_facet_partition_repairs_verification.json'
    output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'models': [{k: v for k, v in item.items() if k != 'peelings'} for item in reports],
                      'nested_flattening_order': len(nested), 'all_checks_passed': True}, indent=2))


if __name__ == '__main__':
    main()
