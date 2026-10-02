#!/usr/bin/env python3
"""Fixed regression for flattening arbitrary proper-row common-seed covers.

The inner cover has a four-word path row, no balanced coordinates, exactly
two resolving additions, and two obstruction-preserving corner predecessors
that also have no balanced coordinates. Its outer two-repair lift flattens
over Hall and descends by eight corners to a terminal antipodal tower.
Obstruction and terminality claims use the proved common-seed criteria;
this script checks their finite inputs and constructions, not a census.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from explore_damp_forced_layer_closure import (
    balanced_coordinates, check_family, signed_clauses,
)
from verify_damp_facet_partition_repairs import build, checked_repair_peeling
from verify_damp_singleton_recovery_implications import shadow, subsets
from verify_damp_terminal_clause_recovery import corner

ROOT = Path(__file__).resolve().parent


def full_facets(rows, mask):
    return {(e, sign): [a for a, row in rows.items()
                       if {t for t in subsets(mask) if bool(t & e) == bool(sign)} <= row]
            for e in [1 << f for f in range(mask.bit_length()) if mask & (1 << f)]
            for sign in [0, 1]}


def peel_row(row, dimension):
    remaining = set(row)
    order = []
    while remaining:
        word = next(x for x in sorted(remaining) if corner(remaining, x, dimension))
        remaining.remove(word)
        order.append(word)
    return order


def main():
    source_path = ROOT / 'computations/hall_corner_exchange/source_H.json'
    source = json.loads(source_path.read_text())
    bits = int(source['source_family_hex'], 16)
    hall = {x for x in range(4096) if bits >> x & 1}
    assert len(hall) == 299 and not any(corner(hall, x, 12) for x in hall)
    certificate_path = ROOT / 'damp_terminal_clause_recovery_verification.json'
    certificate = json.loads(certificate_path.read_text())
    assert certificate['source_mathematical_profile_sha256'] == source['mathematical_profile_sha256']
    orders = {r['tip']: r['peeling'] for r in certificate['repair_proofs']}
    a, b, d = 183, 315, 1836
    gains = {tip: sum(1 << f for f in range(12) if tip ^ (1 << f) in hall)
             for tip in [a, b, d]}
    assert len(set(gains.values())) == 3
    for tip in [a, b, d]:
        assert gains[tip].bit_count() == 4
        remaining = hall | {tip}
        assert len(orders[tip]) == len(remaining)
        for x in orders[tip]:
            assert x in remaining and corner(remaining, x, 12)
            remaining.remove(x)
        assert not remaining

    hall_clauses = {}
    for j in subsets(4095):
        if j.bit_count() == 4:
            missing = set(subsets(j)) - {x & j for x in hall}
            assert len(missing) == 1
            hall_clauses[j] = missing.pop()

    def clauses_for(rows, mask):
        clauses = dict(hall_clauses)
        family = build(hall, rows, mask)
        for tip, row in rows.items():
            assert row and len(shadow(row, mask)) == len(row)
            assert set(row) < set(subsets(mask))
            clauses.pop(gains[tip])
            for j, trace in signed_clauses(row, mask).items():
                clauses[gains[tip] | j] = (tip & gains[tip]) | trace
        # At most three four-set gains create no new old minimal five-set.
        # The replicated-repair theorem gives completeness of this list.
        for j, trace in clauses.items():
            assert set(subsets(j)) - {x & j for x in family} == {trace}
        return clauses

    p, q, r, e = 1 << 12, 1 << 13, 1 << 14, 1 << 15
    inner_mask = p | q | r
    cube = set(subsets(inner_mask))
    path = {0, p, p | q, inner_mask}
    rows = {a: cube - {0}, b: cube - {inner_mask}, d: path}
    seed = build(hall, rows, inner_mask)
    assert len(seed) == 2410
    inner_clauses = clauses_for(rows, inner_mask)
    inner_degrees = [sum(bool(j & (1 << f)) for j in inner_clauses) for f in range(15)]
    assert inner_degrees == [165 + 2 * bool(gains[d] & (1 << f)) for f in range(12)] + [4] * 3
    assert balanced_coordinates(inner_clauses, 4095 | inner_mask) == 0
    assert all(full_facets(rows, inner_mask).values())
    assert not any(d in owners for owners in full_facets(rows, inner_mask).values())
    completions = []
    for old in range(4096):
        missing = {t for t in cube if old | t not in seed}
        if len(missing) == 1:
            completions.append(old | next(iter(missing)))
    assert completions == [a, b | inner_mask]
    repair_peelings = [checked_repair_peeling(hall, rows, inner_mask, tip, orders[tip], 15)
                       for tip in [a, b]]
    assert sorted(rp['added_word'] for rp in repair_peelings) == completions

    predecessor_records = []
    all_corners = {x for x in seed if corner(seed, x, 15)}
    expected = {a | t for t in [p, q, r]} | {
        b | (inner_mask ^ t) for t in [p, q, r]} | {d, d | inner_mask}
    assert all_corners == expected
    for word in sorted(all_corners):
        tip, t = word & 4095, word & inner_mask
        shortened = {label: set(row) for label, row in rows.items()}
        shortened[tip].remove(t)
        # Every shortened row peels, so facet coverage is the exact criterion.
        for row in shortened.values():
            peel_row(row, 15)
        preserves = all(full_facets(shortened, inner_mask).values())
        assert preserves == (tip == d)
        if preserves:
            clauses = clauses_for(shortened, inner_mask)
            degrees = [sum(bool(j & (1 << f)) for j in clauses) for f in range(15)]
            assert degrees == [165 + bool(gains[d] & (1 << f)) for f in range(12)] + [3] * 3
            assert balanced_coordinates(clauses, 4095 | inner_mask) == 0
            predecessor_records.append({'deleted_corner': word, 'order': len(seed) - 1,
                                        'clause_degrees': degrees, 'balanced_coordinates': 0})
    assert len(predecessor_records) == 2
    print('Inner cover: 2410 concepts, two resolving additions, two preserving predecessors, all three balanced sets empty.', flush=True)

    outer_rows = {a: {0}, b | inner_mask: {e}}
    outer = build(seed, outer_rows, e)
    flat_rows = {a: {t | v for t in rows[a] for v in [0, e]} | {0},
                 b: {t | v for t in rows[b] for v in [0, e]} | {inner_mask | e},
                 d: {t | v for t in path for v in [0, e]}}
    total_mask = inner_mask | e
    assert outer == build(hall, flat_rows, total_mask) and len(outer) == 4822
    assert flat_rows[a] == set(subsets(total_mask)) - {e}
    assert flat_rows[b] == set(subsets(total_mask)) - {inner_mask}
    outer_clauses = clauses_for(flat_rows, total_mask)
    balanced, successes, _, bad = check_family(outer, 4095 | total_mask, outer_clauses)
    assert balanced == e and successes == {e} and bad is None
    assert {x & ~e for x in outer if not x & e and x | e in outer} == seed

    path_order = peel_row(flat_rows[d], 16)
    current = set(outer)
    current_rows = {tip: set(row) for tip, row in flat_rows.items()}
    descent = []
    for t in path_order:
        word = d | t
        assert corner(current, word, 16)
        current.remove(word)
        current_rows[d].remove(t)
        assert len(shadow(current_rows[d], total_mask)) == len(current_rows[d])
        assert all(full_facets(current_rows, total_mask).values())
        assert current == build(hall, current_rows, total_mask)
        descent.append({'deleted_corner': word, 'remaining_order': len(current)})
    del current_rows[d]
    assert len(descent) == 8 and len(current) == 4814
    assert all(len(owners) == 1 for owners in full_facets(current_rows, total_mask).values())
    assert current_rows == {a: set(subsets(total_mask)) - {e},
                            b: set(subsets(total_mask)) - {inner_mask}}
    terminal_inner = build(hall, {a: rows[a], b: rows[b]}, inner_mask)
    assert current == build(terminal_inner, outer_rows, e)
    endpoint_clauses = clauses_for(current_rows, total_mask)
    _, endpoint_successes, _, bad = check_family(current, 4095 | total_mask, endpoint_clauses)
    assert endpoint_successes == set(subsets(total_mask)) - {0} and bad is None
    result = {'scope': __doc__.strip(),
              'source_file_sha256': hashlib.sha256(source_path.read_bytes()).hexdigest(),
              'source_mathematical_profile_sha256': source['mathematical_profile_sha256'],
              'repair_certificate_sha256': hashlib.sha256(certificate_path.read_bytes()).hexdigest(),
              'inner_order': len(seed), 'inner_path_row': sorted(path),
              'inner_clause_degrees': inner_degrees,
              'inner_balanced_coordinates': 0, 'inner_corners': sorted(all_corners),
              'complete_resolving_additions': completions, 'repair_peelings': repair_peelings,
              'preserving_corner_predecessors': predecessor_records,
              'outer_order': len(outer), 'outer_successful_layers': sorted(successes),
              'flattened_path_row': sorted(flat_rows[d]), 'corner_descent': descent,
              'terminal_endpoint_order': len(current),
              'endpoint_successful_layers': sorted(endpoint_successes),
              'all_checks_passed': True}
    (ROOT / 'damp_iterated_common_seed_covers_verification.json').write_text(
        json.dumps(result, indent=2) + '\n')
    print('Outer cover: 4822 concepts, unique layer set; eight verified corner deletions reach terminal order 4814.', flush=True)


if __name__ == '__main__':
    main()
