#!/usr/bin/env python3
"""A bounded repair-creation test at the fixed certified 291-concept seed.

Test whether two individually nonresolving additions can cooperate, or
whether a neighboring addition creates a new resolving tip. Stop the
initial ordered-pair test after 100 checks or the first positive result.
Separately check the twelve neighbors of the fixed one-corner extension
at 2184. These negative controls do not prove repair inheritance for
cornerless predecessors. No order, dimension, or Hall-atlas radius grows.
"""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
from analyze_damp_yang_balanced_circuits import reconstruct_291
from explore_damp_cornered_hall_repair_signatures import bounded_peeling
from verify_damp_terminal_clause_recovery import corner, subsets


def ample_tip(family, tip):
    support = sum(1 << e for e in range(12) if tip ^ (1 << e) in family)
    if not all(tip ^ u in family for u in subsets(support) if u):
        return None
    if len({x & support for x in family}) == 1 << support.bit_count():
        return None
    return support


def main():
    seed = set(reconstruct_291())
    assert len(seed) == 291 and not any(corner(seed, x, 12) for x in seed)
    extensions = []
    for tip in range(4096):
        if tip in seed:
            continue
        support = ample_tip(seed, tip)
        if support is not None:
            extensions.append({'tip': tip, 'support': support,
                               **bounded_peeling(seed | {tip}, 1000)})
    bad = [r for r in extensions if r['status'] == 'nonpeelable']
    known = {r['tip']: r for r in extensions}
    pairs = []
    for row in bad:
        c = row['tip']
        family = seed | {c}
        candidates = sorted(({c ^ (1 << e) for e in range(12)}
                             | {r['tip'] for r in bad}) - family)
        for tip in candidates:
            support = ample_tip(family, tip)
            if support is None:
                continue
            decision = bounded_peeling(family | {tip}, 1000)
            pairs.append({'c': c, 'a': tip, 'new_support': support,
                          'adjacent': (tip ^ c).bit_count() == 1,
                          'base_a_status': known.get(tip, {}).get('status', 'not_ample'),
                          **decision})
            if decision['status'] == 'peelable' or len(pairs) == 100:
                break
        if pairs and (pairs[-1]['status'] == 'peelable' or len(pairs) == 100):
            break
    c = 2184
    family = seed | {c}
    assert [x for x in sorted(family) if corner(family, x, 12)] == [c]
    nested = []
    for e in range(12):
        tip = c ^ (1 << e)
        if tip in family:
            continue
        support = ample_tip(family, tip)
        if support is not None:
            nested.append({'c': c, 'a': tip, 'new_support': support,
                           **bounded_peeling(family | {tip}, 10000)})
    result = {'scope': __doc__.strip(), 'base': sorted(seed),
              'source_model_sha256': hashlib.sha256((ROOT / 'damp_hall_carrier_defect_min8.json').read_bytes()).hexdigest(),
              'reconstruction_helper_sha256': hashlib.sha256((ROOT / 'analyze_damp_yang_balanced_circuits.py').read_bytes()).hexdigest(),
              'extensions': extensions, 'ordered_pair_checks': pairs,
              'fixed_neighbor_checks': nested,
              'initial_pair_budget': 100, 'initial_state_limit': 1000,
              'fixed_neighbor_state_limit': 10000}
    output = ROOT / 'damp_fixed_291_repair_creation.json'
    output.write_text(json.dumps(result, indent=2) + '\n')
    for name, rows in [('extensions', extensions), ('ordered_pair_checks', pairs),
                       ('fixed_neighbor_checks', nested)]:
        print(name, {status: sum(r['status'] == status for r in rows)
                     for status in ['peelable', 'nonpeelable', 'unresolved']})


if __name__ == '__main__':
    main()
