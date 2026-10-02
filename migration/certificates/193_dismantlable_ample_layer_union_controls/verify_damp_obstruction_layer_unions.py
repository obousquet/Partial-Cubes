#!/usr/bin/env python3
"""Fixed regressions for union closure and the full layer-recovery dichotomy.

This is not an obstruction census or a proof of the general theorem.
The Hall seed and its certified repairs give the obstruction inputs by the
assembly theorems. We check all candidate layer sets on six fixed inputs,
the disjoint and overlapping reductions, signed seed-coordinate inheritance,
the sublayer fibre formula, the thin-coordinate and repair-count dichotomy,
and the peelable six-word counterexample to unrestricted union closure.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

from explore_damp_forced_layer_closure import (
    balanced_coordinates, check_family, fibres, signed_clauses, successful,
)
from verify_damp_singleton_recovery_implications import subsets

ROOT = Path(__file__).resolve().parent


def main():
    source_path = ROOT / 'computations/hall_corner_exchange/source_H.json'
    source = json.loads(source_path.read_text())
    bits = int(source['source_family_hex'], 16)
    hall = {x for x in range(4096) if bits >> x & 1}
    assert len(hall) == 299
    hall_clauses = {}
    for j in subsets(4095):
        if j.bit_count() == 4:
            missing = set(subsets(j)) - {x & j for x in hall}
            assert len(missing) == 1
            hall_clauses[j] = missing.pop()
    a, b, d = 183, 315, 1836
    gains = {tip: sum(1 << i for i in range(12) if tip ^ (1 << i) in hall)
             for tip in [a, b, d]}
    assert len(set(gains.values())) == 3
    for tip, gain in gains.items():
        assert gain.bit_count() == 4 and hall_clauses[gain] == tip & gain
        assert all(tip ^ t in hall for t in subsets(gain) if t)

    def model(rows, layers):
        """The replicated-repair formula; all retained rows miss one face."""
        cube = set(subsets(layers))
        rows = {tip: row for tip, row in rows.items() if row}
        family = {x | t for x in hall for t in cube}
        clauses = dict(hall_clauses)
        for tip, row in rows.items():
            assert row < cube
            missing = cube - row
            first = min(missing)
            varying = 0
            for t in missing:
                varying |= first ^ t
            support = layers ^ varying
            trace = first & support
            assert support and missing == {t for t in cube if t & support == trace}
            gain = gains[tip]
            clauses.pop(gain)
            clauses[gain | support] = (tip & gain) | trace
            family.update(tip | t for t in row)
        # At most three four-set gains cannot create a new minimal five-set.
        # The assembly formula therefore gives the complete clause list.
        for j, trace in clauses.items():
            assert set(subsets(j)) - {x & j for x in family} == {trace}
        return family, clauses, 4095 | layers

    def descend(rows, layers, removed):
        remaining = layers ^ removed
        descended = {tip: {t for t in subsets(remaining)
                            if all(t | u in row for u in subsets(removed))}
                     for tip, row in rows.items()}
        return model(descended, remaining)

    p, q, r = 1 << 12, 1 << 13, 1 << 14
    examples = [
        ('two_layers', p | q,
         {a: set(subsets(p | q)) - {0}, b: set(subsets(p | q)) - {p | q}},
         1202, {p, q, p | q}),
        ('three_layers', p | q | r,
         {a: set(subsets(p | q | r)) - {0},
          b: set(subsets(p | q | r)) - {p | q | r}},
         2406, set(subsets(p | q | r)) - {0}),
        ('three_owners', p | q,
         {a: {0, p, q}, b: {p, p | q}, d: {q, p | q}},
         1203, {p | q}),
        ('cornered_largest_seed', p | q,
         {a: {p, p | q}, b: {0, q, p | q}, d: {0, p, q}},
         1204, {q}),
        ('terminal_endpoint', p | q,
         {a: {p, p | q}, b: {0, q, p | q}, d: {0, p}},
         1203, {p | q}),
        ('three_owners_no_resolving_addition', p | q | r,
         {a: {t for t in subsets(p | q | r) if not t & p or not t & q},
          b: {t for t in subsets(p | q | r) if t & q or not t & r},
          d: {t for t in subsets(p | q | r) if t & p or t & r}},
         2410, {p | q | r}),
    ]
    records = []
    for name, layers, rows, order, expected in examples:
        family, clauses, ambient = model(rows, layers)
        assert len(family) == order
        balanced, successes, trajectories, bad = check_family(family, ambient, clauses)
        assert successes == expected and bad is None
        maximum = 0
        for e in successes:
            maximum |= e
        starts = [t for t in trajectories if t['steps'][0]['layers'].bit_count() == 1]
        detected = sum(t['steps'][0]['layers'] for t in starts if t['status'] == 'success')
        assert detected == maximum and maximum in successes
        thin = 0
        for e in [1 << f for f in range(15) if ambient & (1 << f)]:
            zero = {x for x in family if not x & e}
            one = {x ^ e for x in family if x & e}
            if len(zero - one) == len(one - zero) == 1:
                thin |= e
        maximum_fibres = fibres(family, clauses, ambient, maximum)
        owner_count = len(maximum_fibres)
        resolving_count = sum(j & maximum == maximum for j in maximum_fibres)
        if thin:
            assert thin == maximum and owner_count == resolving_count == 2
            assert successes == set(subsets(thin)) - {0}
        else:
            assert successes == {maximum} and maximum.bit_count() >= 2
            assert owner_count >= 3 and resolving_count <= 1
        fibre_checks = 0
        for e in successes:
            full_fibres = fibres(family, clauses, ambient, e)
            for sub in subsets(e):
                if not sub or sub == e:
                    continue
                remaining = e ^ sub
                sub_fibres = fibres(family, clauses, ambient, sub)
                for j, actual in sub_fibres.items():
                    tip = next(iter(full_fibres[j]))
                    expected_fibre = {tip | t for t in subsets(remaining)
                                      if t & j == clauses[j] & remaining}
                    assert actual == expected_fibre
                    fibre_checks += 1
        inherited = []
        for e in sorted(successes):
            seed, seed_clauses, seed_ambient = descend(rows, layers, e)
            seed_balanced = balanced_coordinates(seed_clauses, seed_ambient)
            assert (balanced & ~e) & ~seed_balanced == 0
            assert all({bool(trace & (1 << f)) for j, trace in seed_clauses.items()
                        if j & (1 << f)} == {False, True}
                       for f in range(15) if seed_ambient & (1 << f))
            inherited.append({'layers': e, 'seed_order': len(seed),
                              'seed_balanced': seed_balanced})
        disjoint = overlapping = 0
        for e in successes:
            for f in successes:
                assert e | f in successes
                if (e & f) in [e, f]:
                    continue
                intersection = e & f
                y, yc, ya = descend(rows, layers, intersection)
                ep, fp = e ^ intersection, f ^ intersection
                assert successful(y, yc, ya, ep) and successful(y, yc, ya, fp)
                seed, sc, sa = descend(rows, layers, e)
                assert fp & ~balanced_coordinates(sc, sa) == 0
                assert successful(seed, sc, sa, fp)
                if intersection:
                    overlapping += 1
                else:
                    disjoint += 1
        records.append({'name': name, 'order': order, 'balanced': balanced,
                        'successful_layer_sets': sorted(successes), 'largest': maximum,
                        'thin_coordinates': thin, 'largest_owner_count': owner_count,
                        'resolving_addition_count_by_fibre_theorem': resolving_count,
                        'sublayer_fibre_checks': fibre_checks,
                        'singleton_starts': starts, 'seed_checks': inherited,
                        'ordered_disjoint_pairs': disjoint,
                        'ordered_overlapping_incomparable_pairs': overlapping})
        print(json.dumps({k: v for k, v in records[-1].items()
                          if k not in ['singleton_starts', 'seed_checks']}), flush=True)

    control = {0, 1, 2, 3, 7, 8}
    cb, cs, _, bad = check_family(control, 15, signed_clauses(control, 15))
    assert cb == 3 and cs == {1, 2} and bad is not None
    certificate_path = ROOT / 'damp_terminal_clause_recovery_verification.json'
    certificate = json.loads(certificate_path.read_text())
    assert certificate['source_mathematical_profile_sha256'] == source['mathematical_profile_sha256']
    assert {a, b, d} <= {entry['tip'] for entry in certificate['repair_proofs']}
    result = {'scope': __doc__.strip(),
              'source_file_sha256': hashlib.sha256(source_path.read_bytes()).hexdigest(),
              'source_mathematical_profile_sha256': source['mathematical_profile_sha256'],
              'repair_certificate_sha256': hashlib.sha256(certificate_path.read_bytes()).hexdigest(),
              'fixed_obstructions': records,
              'peelable_control': {'family': sorted(control), 'successes': sorted(cs),
                                   'failed_union': list(bad)},
              'all_checks_passed': True}
    (ROOT / 'damp_obstruction_layer_unions_verification.json').write_text(
        json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
