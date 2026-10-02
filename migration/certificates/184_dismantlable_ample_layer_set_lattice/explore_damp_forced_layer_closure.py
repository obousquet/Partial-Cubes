#!/usr/bin/env python3
"""Check forced layer enlargement and whether successful sets are union-closed.

The bounded domain is all nonempty ample families through Q4. The fixed
Hall control is checked separately. This is a recovery diagnostic, not a
search for forbidden pc-minors.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from verify_damp_singleton_recovery_implications import shadow, subsets


def signed_clauses(family, ambient):
    sh = shadow(family, ambient)
    clauses = {}
    for j in subsets(ambient):
        if j in sh or not all(j ^ e in sh for e in subsets(j) if e.bit_count() == 1):
            continue
        missing = set(subsets(j)) - {x & j for x in family}
        assert len(missing) == 1
        clauses[j] = missing.pop()
    return clauses


def balanced_coordinates(clauses, ambient):
    result = 0
    for e in subsets(ambient):
        if e.bit_count() != 1:
            continue
        incident = [j for j in clauses if j & e]
        if len(incident) == 2 and (clauses[incident[0]] ^ clauses[incident[1]]) & e:
            result |= e
    return result


def fibres(family, clauses, ambient, layers):
    old = ambient ^ layers
    projection = {x & old for x in family}
    return {j: {x for x in projection if x & j == clauses[j] & old}
            for j in clauses if j & layers}


def successful(family, clauses, ambient, layers):
    rows = list(fibres(family, clauses, ambient, layers).values())
    return all(len(row) == 1 for row in rows) and len({next(iter(row)) for row in rows}) == len(rows)


def closure(family, clauses, ambient, balanced, initial):
    layers = initial
    steps = []
    while True:
        rows = fibres(family, clauses, ambient, layers)
        assert all(rows.values())
        variation = 0
        for row in rows.values():
            first = next(iter(row))
            for word in row:
                variation |= first ^ word
        steps.append({'layers': layers, 'variation': variation})
        if variation & ~balanced:
            forbidden = variation & ~balanced
            bit = forbidden & -forbidden
            j, row = next((j, row) for j, row in rows.items()
                          if len({word & bit for word in row}) == 2)
            pair = [next(word for word in row if not word & bit),
                    next(word for word in row if word & bit)]
            return {'status': 'outside_balanced', 'steps': steps,
                    'witness': {'clause': j, 'words': pair, 'coordinate_mask': bit}}
        keys = list(rows)
        for k, j in enumerate(keys):
            for i in keys[:k]:
                overlap = rows[i] & rows[j]
                if overlap:
                    return {'status': 'collision', 'steps': steps,
                            'witness': {'clauses': [i, j], 'word': min(overlap)}}
        if variation:
            layers |= variation
            continue
        return {'status': 'success', 'layers': layers, 'steps': steps}


def check_family(family, ambient, clauses=None):
    if clauses is None:
        clauses = signed_clauses(family, ambient)
    balanced = balanced_coordinates(clauses, ambient)
    candidates = [e for e in subsets(balanced) if e]
    successes = {e for e in candidates if successful(family, clauses, ambient, e)}
    results = []
    for e in candidates:
        result = closure(family, clauses, ambient, balanced, e)
        containing = {s for s in successes if e & s == e}
        assert (result['status'] == 'success') == bool(containing)
        if containing:
            assert result['layers'] in containing
            assert all(result['layers'] & s == result['layers'] for s in containing)
        results.append(result)
    bad_union = next(((e, f, e | f) for e in successes for f in successes
                      if e | f not in successes), None)
    return balanced, successes, results, bad_union


def hall_control():
    root = Path(__file__).resolve().parent
    source = json.loads((root / 'computations/hall_corner_exchange/source_H.json').read_text())
    bits = int(source['source_family_hex'], 16)
    seed = {x for x in range(1 << 12) if bits >> x & 1}
    assert len(seed) == 299 and source['source_label'] == 'H'
    clauses = {}
    for support in subsets((1 << 12) - 1):
        if support.bit_count() == 4:
            missing = set(subsets(support)) - {x & support for x in seed}
            assert len(missing) == 1
            clauses[support] = missing.pop()
    p, q = 1 << 12, 1 << 13
    family = {x | t for x in seed for t in [0, p, q, p | q]}
    for tip, rows, layers, trace in [
            (183, [0, p, q], p | q, p | q),
            (315, [p, p | q], p, 0),
            (1836, [q, p | q], q, 0)]:
        support = sum(1 << e for e in range(12) if tip ^ (1 << e) in seed)
        assert clauses.pop(support) == tip & support
        clauses[support | layers] = (tip & support) | trace
        family.update(tip | t for t in rows)
    ambient = (1 << 14) - 1
    assert len(clauses) == 495 and len(family) == 1203
    balanced, successes, results, bad_union = check_family(family, ambient, clauses)
    assert balanced == p | q and successes == {p | q} and bad_union is None
    assert sum(len(r['steps']) > 1 for r in results) == 2
    return {'source_mathematical_profile_sha256': source['mathematical_profile_sha256'],
            'balanced': balanced, 'successes': sorted(successes), 'closures': results}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-dimension', type=int, default=4)
    args = parser.parse_args()
    records = []
    union_counterexample = None
    failure_examples = {}
    for n in range(1, args.max_dimension + 1):
        ambient = (1 << n) - 1
        counts = {'dimension': n, 'ample_families': 0, 'candidate_sets': 0,
                  'successes': 0, 'forced_growth': 0, 'outside_balanced': 0, 'collision': 0}
        for bits in range(1, 1 << (1 << n)):
            family = {x for x in range(1 << n) if bits >> x & 1}
            if len(shadow(family, ambient)) != len(family):
                continue
            counts['ample_families'] += 1
            balanced, successes, results, bad_union = check_family(family, ambient)
            counts['candidate_sets'] += len(results)
            counts['successes'] += len(successes)
            counts['forced_growth'] += sum(len(r['steps']) > 1 for r in results)
            for result in results:
                if result['status'] != 'success':
                    counts[result['status']] += 1
                    failure_examples.setdefault(result['status'], {
                        'dimension': n, 'family': sorted(family),
                        'balanced': balanced, 'closure': result})
            if bad_union and union_counterexample is None:
                union_counterexample = {'dimension': n, 'family': sorted(family),
                                        'balanced': balanced, 'sets_and_union': bad_union}
        records.append(counts)
        print(json.dumps(counts), flush=True)
    output = {'scope': __doc__.strip(), 'counts': records,
              'union_counterexample': union_counterexample,
              'failure_examples': failure_examples,
              'hall_control': hall_control(),
              'forced_closure_checks_passed': True}
    path = Path(__file__).with_name('damp_forced_layer_closure_verification.json')
    path.write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
