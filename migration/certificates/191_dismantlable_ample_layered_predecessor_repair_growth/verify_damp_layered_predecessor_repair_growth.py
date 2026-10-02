#!/usr/bin/env python3
"""Verify the sharp 600-to-601 repair-growth example and its outer descent.

All inputs are the fixed Hall seed and three already certified repair
tips. No new obstruction search is performed. The theorem supplies
completeness of resolving-fibre tests and pc-minimality; this script
checks the constructions, complete candidate lists, explicit repair
peelings, acquired supports, clause degrees, and the terminal endpoint.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from verify_damp_facet_partition_repairs import build, checked_repair_peeling
from verify_damp_terminal_clause_recovery import corner, subsets

ROOT = Path(__file__).resolve().parent


def completions(family, mask):
    layers = set(subsets(mask))
    result = []
    for old in range(1 << 12):
        missing = {t for t in layers if old | t not in family}
        if len(missing) == 1:
            result.append(old | next(iter(missing)))
    return sorted(result)


def clause_degrees(family, supports, dimension):
    # Each listed support really has exactly one omitted trace. The
    # replicated-repair theorem gives completeness of the support list.
    for support in supports:
        assert len({x & support for x in family}) == (1 << support.bit_count()) - 1
    return [sum(bool(s & (1 << f)) for s in supports) for f in range(dimension)]


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
    p, q = 1 << 12, 1 << 13
    acquired = {tip: sum(1 << e for e in range(12) if tip ^ (1 << e) in hall)
                for tip in [a, b, d]}
    assert len(set(acquired.values())) == 3
    for tip, support in acquired.items():
        assert support.bit_count() == 4
        assert all(tip ^ t in hall for t in subsets(support) if t)
    rows_h = {a: {p}, b: {0}}
    rows_s = {**rows_h, d: {0}}
    h = build(hall, rows_h, p)
    s = build(hall, rows_s, p)
    assert len(h) == 600 and s == h | {d} and len(s) == 601
    assert {x for x in h if corner(h, x, 13)} == {a | p, b}
    assert {x for x in s if corner(s, x, 13)} == {a | p, b, d}
    assert completions(h, p) == [a, b | p]
    assert completions(s, p) == [a, b | p, d | p]
    reports = []
    for name, family, rows in [('H', h, rows_h), ('S', s, rows_s)]:
        repairs = [checked_repair_peeling(hall, rows, p, tip, orders[tip], 13)
                   for tip in rows]
        assert sorted(r['added_word'] for r in repairs) == completions(family, p)
        reports.append({'name': name, 'order': len(family),
                        'complete_resolving_tips': completions(family, p),
                        'repair_peelings': repairs})
    assert (d ^ (d | p)).bit_count() == 1
    assert all((d ^ old).bit_count() > 1 for old in completions(h, p))
    assert {x for x in s if not x & p and x | p in s} == hall

    old_supports = [j for j in subsets(4095)
                    if j.bit_count() == 4 and j not in set(acquired.values())]
    s_clauses = old_supports + [support | p for support in acquired.values()]
    s_degrees = clause_degrees(s, s_clauses, 13)
    assert s_degrees == [165] * 12 + [3]
    # No coordinate can satisfy the degree-two balanced-clause condition.
    assert 2 not in s_degrees

    u, v = d | p, b | p
    outer = build(s, {u: {0}, v: {q}}, q)
    assert len(outer) == 1204
    deleted = d | q
    assert not corner(outer, d, 14)
    assert corner(outer, deleted, 14)
    after = outer - {deleted}
    final_rows = {a: {p, p | q}, b: {0, q, p | q}, d: {0, p}}
    assert after == build(hall, final_rows, p | q) and len(after) == 1203
    facet_owners = []
    for e in [p, q]:
        for sign in [0, 1]:
            facet = {t for t in subsets(p | q) if bool(t & e) == bool(sign)}
            owners = [tip for tip, row in final_rows.items() if facet <= row]
            assert len(owners) == 1
            facet_owners.append({'coordinate': e.bit_length() - 1, 'sign': sign,
                                 'owner': owners[0]})
    outer_clauses = old_supports + [acquired[a] | p, acquired[b] | p | q,
                                   acquired[d] | p | q]
    outer_degrees = clause_degrees(outer, outer_clauses, 14)
    assert outer_degrees == [165] * 12 + [3, 2]
    endpoint_clauses = old_supports + [acquired[a] | p, acquired[b] | p | q,
                                      acquired[d] | q]
    endpoint_degrees = clause_degrees(after, endpoint_clauses, 14)
    assert endpoint_degrees == [165] * 12 + [2, 2]
    result = {'scope': __doc__.strip(),
              'source_file_sha256': hashlib.sha256(source_path.read_bytes()).hexdigest(),
              'source_mathematical_profile_sha256': source['mathematical_profile_sha256'],
              'repair_certificate_sha256': hashlib.sha256(certificate_path.read_bytes()).hexdigest(),
              'models': reports,
              'added_seed_corner': d, 'unique_new_resolving_tip': d | p,
              'seed_clause_degrees': s_degrees,
              'outer_order': len(outer), 'outer_deleted_corner': deleted,
              'outer_clause_degrees': outer_degrees, 'endpoint_order': len(after),
              'endpoint_clause_degrees': endpoint_degrees, 'endpoint_facet_owners': facet_owners,
              'all_checks_passed': True}
    output = ROOT / 'damp_layered_predecessor_repair_growth_verification.json'
    output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'models'}, indent=2))


if __name__ == '__main__':
    main()
