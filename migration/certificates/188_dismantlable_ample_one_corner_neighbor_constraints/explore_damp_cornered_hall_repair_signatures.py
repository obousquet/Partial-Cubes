#!/usr/bin/env python3
"""Test the four fixed cornered extensions in the original Hall atlas.

Enumerate their one-concept ample extensions by the local punctured-cube
test, obtain bounded exact peeling decisions, and track which old corners
each repair destroys. This is a fixed-adversary diagnostic, not a census
by order or dimension. State-limit exhaustion is reported as unresolved.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from verify_damp_terminal_clause_recovery import corner, subsets

ROOT = Path(__file__).resolve().parent
DIMENSION = 12


class LimitReached(Exception):
    pass


def bounded_peeling(family, limit=10000):
    failed = {}
    initial = set(family)
    visited = 0

    def visit(remaining):
        nonlocal visited
        key = frozenset(remaining)
        if key in failed:
            return None
        if not remaining:
            return []
        visited += 1
        if visited > limit:
            raise LimitReached
        choices = tuple(x for x in sorted(remaining) if corner(remaining, x, DIMENSION))
        for word in choices:
            result = visit(remaining - {word})
            if result is not None:
                return [word] + result
        failed[key] = choices
        return None

    try:
        order = visit(set(family))
    except LimitReached:
        return {'status': 'unresolved', 'visited': visited}
    if order is None:
        certificate = []
        for state, choices in failed.items():
            remaining = set(state)
            assert remaining
            assert choices == tuple(x for x in sorted(remaining) if corner(remaining, x, DIMENSION))
            assert all(frozenset(remaining - {x}) in failed for x in choices)
            certificate.append({'removed': sorted(initial - remaining), 'corners': list(choices)})
        certificate.sort(key=lambda row: (len(row['removed']), row['removed']))
        return {'status': 'nonpeelable', 'visited': visited,
                'complete_corner_deletion_certificate': certificate}
    remaining = set(family)
    for word in order:
        assert word in remaining and corner(remaining, word, DIMENSION)
        remaining.remove(word)
    assert not remaining
    return {'status': 'peelable', 'visited': visited, 'peeling': order}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--seed-tip', type=int, action='append')
    parser.add_argument('--state-limit', type=int, default=10000)
    args = parser.parse_args()
    source_path = ROOT / 'computations/hall_corner_exchange/source_H.json'
    source = json.loads(source_path.read_text())
    bits = int(source['source_family_hex'], 16)
    hall = {x for x in range(1 << DIMENSION) if bits >> x & 1}
    assert len(hall) == 299 and not any(corner(hall, x, DIMENSION) for x in hall)
    for support in subsets((1 << DIMENSION) - 1):
        if support.bit_count() == 4:
            assert len({x & support for x in hall}) == 15
    inputs = {r['added']: r for r in source['extensions'] if not r['dismantlable']}
    original_repairs = {r['added']: sum(1 << e for e in range(DIMENSION)
                        if r['added'] ^ (1 << e) in hall)
                        for r in source['extensions'] if r['dismantlable']}
    original_repair_proofs = {}
    for tip in original_repairs:
        proof = bounded_peeling(hall | {tip}, args.state_limit)
        assert proof['status'] == 'peelable'
        original_repair_proofs[tip] = proof['peeling']
    selected = args.seed_tip or sorted(inputs)
    reports = []
    for added in selected:
        stored = inputs[added]
        seed = hall | {added}
        old_support = sum(1 << e for e in range(DIMENSION) if added ^ (1 << e) in hall)
        assert old_support.bit_count() == 4
        assert all(added ^ t in hall for t in subsets(old_support) if t)
        old_corners = sorted(x for x in seed if corner(seed, x, DIMENSION))
        assert old_corners == stored['corners'] and len(old_corners) == 2
        # Both possible first corner deletions are cornerless; this is a
        # complete nonpeelability certificate for the fixed 300-word seed.
        for x in old_corners:
            child = seed - {x}
            assert not any(corner(child, y, DIMENSION) for y in child)
        candidates = []
        for tip in range(1 << DIMENSION):
            if tip in seed:
                continue
            support = sum(1 << e for e in range(DIMENSION) if tip ^ (1 << e) in seed)
            # Sh(seed) consists of all sets of size <=3 and old_support.
            # Its minimal nonfaces are exactly the other four-sets.
            if support.bit_count() != 4 or support == old_support:
                continue
            if not all(tip ^ t in seed for t in subsets(support) if t):
                continue
            extension = seed | {tip}
            decision = bounded_peeling(extension, args.state_limit)
            destroyed = [x for x in old_corners if not corner(extension, x, DIMENSION)]
            candidates.append({'tip': tip, 'support': support,
                               'destroys_old_corners': destroyed, **decision})
        summary = {'seed_added_word': added, 'corners': old_corners,
                   'ample_extensions': len(candidates),
                   'resolving_extensions': sum(r['status'] == 'peelable' for r in candidates),
                   'nonpeelable_extensions': sum(r['status'] == 'nonpeelable' for r in candidates),
                   'unresolved_extensions': sum(r['status'] == 'unresolved' for r in candidates),
                   'dependent_resolving_tips': [r['tip'] for r in candidates
                        if r['status'] == 'peelable' and r['destroys_old_corners']]}
        resolved = {row['tip']: row['support'] for row in candidates if row['status'] == 'peelable'}
        assert resolved == original_repairs
        assert old_support not in resolved.values()
        summary['same_resolving_tips_and_supports_as_Hall'] = True
        # Lift over the cornered seed with two inherited Hall repairs.
        # Delete the added corner's two layer copies to recover the usual
        # 600-concept Hall lift. Each step passes the local ample-extension
        # test in reverse; pc-minimality follows from the one-concept atlas.
        a, b, e = 183, 315, 1 << DIMENSION
        base = (hall | {b}) | {x | e for x in hall | {a}}
        current = set(base)
        additions = []
        for word in [added, added | e]:
            support = sum(1 << f for f in range(DIMENSION + 1)
                          if word ^ (1 << f) in current)
            assert all(word ^ t in current for t in subsets(support) if t)
            assert len({x & support for x in current}) < 1 << support.bit_count()
            current.add(word)
            assert corner(current, word, DIMENSION + 1)
            strong = {x for x in current if not x & e and x | e in current}
            assert strong in [hall, seed]
            additions.append({'added_word': word, 'acquired_support': support,
                              'order': len(current), 'strong_projection_order': len(strong)})
        assert len(current) == 602
        # Two-repair trace transport gives all four-set clauses except
        # the three acquired supports, plus the two mixed clauses.
        support_a, support_b = resolved[a], resolved[b]
        clause_supports = [j for j in subsets((1 << DIMENSION) - 1)
                           if j.bit_count() == 4 and j not in {old_support, support_a, support_b}]
        clause_supports += [support_a | e, support_b | e]
        degrees = [sum(bool(j & (1 << f)) for j in clause_supports) for f in range(DIMENSION + 1)]
        assert degrees == [165 - bool(old_support & (1 << f)) for f in range(DIMENSION)] + [2]
        for word in [added | e, added]:
            assert corner(current, word, DIMENSION + 1)
            current.remove(word)
        assert current == base and len(base) == 600
        reports.append({**summary, 'extensions': candidates,
                        'lift_additions': additions,
                        'lift_clause_degrees': degrees,
                        'lift_unique_successful_layer_coordinate': DIMENSION,
                        'obstruction_preserving_lift_deletions': [added | e, added]})
        print(json.dumps(summary), flush=True)
    result = {'scope': __doc__.strip(),
              'source_file_sha256': hashlib.sha256(source_path.read_bytes()).hexdigest(),
              'source_mathematical_profile_sha256': source['mathematical_profile_sha256'],
              'original_hall_repair_peelings': original_repair_proofs,
              'families': reports}
    path = ROOT / 'damp_cornered_hall_repair_signatures.json'
    path.write_text(json.dumps(result, indent=2) + '\n')


if __name__ == '__main__':
    main()
