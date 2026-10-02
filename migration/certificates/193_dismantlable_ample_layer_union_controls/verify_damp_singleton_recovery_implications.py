#!/usr/bin/env python3
"""Finite checks of the automatic conclusions in clause recovery.

Exhaust all families through Q_3, then check two fixed multiple-layer
constructions. This checks the recovery lemma, not obstruction minimality
or the open assertion that every terminal obstruction admits recovery.
"""
from __future__ import annotations

import json
from pathlib import Path


def subsets(mask: int):
    sub = mask
    while True:
        yield sub
        if not sub:
            break
        sub = (sub - 1) & mask


def shadow(family: set[int], coordinates: int) -> set[int]:
    return {s for s in subsets(coordinates)
            if len({x & s for x in family}) == 1 << s.bit_count()}


def check(family: set[int], dimension: int) -> dict:
    ambient = (1 << dimension) - 1
    sh = shadow(family, ambient)
    assert family and len(sh) == len(family)
    clauses = {}
    for j in subsets(ambient):
        if j in sh or not all(j ^ (1 << e) in sh
                                 for e in range(dimension) if j >> e & 1):
            continue
        missing = set(subsets(j)) - {x & j for x in family}
        assert len(missing) == 1
        clauses[j] = missing.pop()
    balanced = 0
    for e in range(dimension):
        incident = [j for j in clauses if j >> e & 1]
        if len(incident) == 2 and ((clauses[incident[0]] ^ clauses[incident[1]]) >> e & 1):
            balanced |= 1 << e
    successes = []
    candidates = 0
    for layers in subsets(balanced):
        if not layers:
            continue
        candidates += 1
        old = ambient ^ layers
        c = {x for x in subsets(old)
             if all(x & j != trace for j, trace in clauses.items() if not j & layers)}
        assert c == {x & old for x in family}
        mixed = [j for j in clauses if j & layers]
        fibres = [{x for x in c if x & j == clauses[j] & old} for j in mixed]
        if not all(len(t) == 1 for t in fibres):
            continue
        tips = [next(iter(t)) for t in fibres]
        if len(set(tips)) != len(tips):
            continue
        seed = c - set(tips)
        strong = {x for x in subsets(old)
                  if all(x | t in family for t in subsets(layers))}
        assert seed == strong and seed
        seed_shadow = shadow(seed, old)
        assert len(seed_shadow) == len(seed)
        assert seed_shadow == {d for d in subsets(old) if d | layers in sh}
        supports = []
        for j, tip in zip(mixed, tips):
            repair = seed | {tip}
            repair_shadow = shadow(repair, old)
            assert len(repair_shadow) == len(repair)
            assert repair_shadow - seed_shadow == {j & old}
            assert seed_shadow <= repair_shadow
            supports.append(j & old)
        assert len(set(supports)) == len(supports)
        rebuilt = {x | t for x in seed for t in subsets(layers)}
        for j, tip in zip(mixed, tips):
            rebuilt.update(tip | t for t in subsets(layers)
                           if t & j != clauses[j] & layers)
        assert rebuilt == family
        successes.append({'layers_mask': layers, 'seed_order': len(seed),
                          'tips': tips, 'supports': supports})
    return {'candidates': candidates, 'successes': successes}


def main() -> None:
    exhaustive = []
    for n in range(1, 4):
        counts = {'dimension': n, 'families': (1 << (1 << n)) - 1,
                  'ample': 0, 'candidate_layer_sets': 0, 'successful_layer_sets': 0}
        for bits in range(1, 1 << (1 << n)):
            family = {x for x in range(1 << n) if bits >> x & 1}
            if len(shadow(family, (1 << n) - 1)) != len(family):
                continue
            counts['ample'] += 1
            result = check(family, n)
            counts['candidate_layer_sets'] += result['candidates']
            counts['successful_layer_sets'] += len(result['successes'])
        exhaustive.append(counts)

    # Two antipodally punctured rows over a singleton seed in Q_2.
    two_rows = {0, 4, 8, 12} | {1, 5, 9} | {6, 10, 14}
    two = check(two_rows, 4)
    assert any(r['layers_mask'] == 12 and r['seed_order'] == 1 for r in two['successes'])
    # Three facet owners over a singleton seed in Q_3.
    three_rows = {0, 8, 16, 24} | {1, 9, 17} | {10, 26} | {20, 28}
    three = check(three_rows, 5)
    assert any(r['layers_mask'] == 24 and len(r['tips']) == 3 for r in three['successes'])
    result = {'scope': __doc__.strip(), 'exhaustive': exhaustive,
              'fixed_two_row_control': two, 'fixed_three_row_control': three,
              'all_checks_passed': True}
    output = Path(__file__).with_name('damp_singleton_recovery_implications_verification.json')
    output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
