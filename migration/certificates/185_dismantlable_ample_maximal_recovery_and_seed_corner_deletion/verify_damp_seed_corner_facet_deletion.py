#!/usr/bin/env python3
"""Check the seed-corner deletion theorem on fixed Hall assemblies.

This is a proof regression, not an obstruction census. Check the local
corner formula on the 602- and 1202-concept controls, give explicit
peelings of every elementary image after each of the two prescribed
602-to-601 deletions, and exhibit the bad minors that prevent descent
through the cornered terminal 600-concept seed of the nested control.
Peeling search is bounded; state exhaustion is unresolved verification,
never evidence of nonpeelability.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path

from verify_damp_facet_partition_repairs import build
from verify_damp_terminal_clause_recovery import corner, subsets

ROOT = Path(__file__).resolve().parent


def checked_peeling(family, dimension, limit=10000):
    failed = set()
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
            raise RuntimeError('Peeling state limit reached: verification unresolved.')
        for word in sorted(remaining):
            if corner(remaining, word, dimension):
                suffix = visit(remaining - {word})
                if suffix is not None:
                    return [word] + suffix
        failed.add(key)
        return None

    order = visit(set(family))
    assert order is not None, 'Complete search contradicts the predicted positive image.'
    remaining = set(family)
    for word in order:
        assert word in remaining and corner(remaining, word, dimension)
        remaining.remove(word)
    assert not remaining
    return order, visited


def elementary_images(family, dimension):
    for e in range(dimension):
        bit = 1 << e
        for sign in [0, 1]:
            yield ('restriction', e, sign), {x & ~bit for x in family if bool(x & bit) == bool(sign)}
        yield ('contraction', e, None), {x & ~bit for x in family}


def check_corner_formula(seed, rows, mask, dimension):
    family = build(seed, rows, mask)
    seed_corners = [s for s in sorted(seed) if corner(seed, s, dimension)]
    predicted = {a | t for a, row in rows.items() for t in row
                 if corner(row, t, dimension)}
    for s in seed_corners:
        neighbors = [a for a in rows if (s ^ a).bit_count() == 1]
        predicted.update(s | t for t in subsets(mask)
                         if all(t not in rows[a] for a in neighbors))
    actual = {x for x in family if corner(family, x, dimension)}
    assert predicted == actual
    return {'order': len(family), 'seed_corners': seed_corners,
            'assembly_corners': sorted(actual)}


def main():
    source_path = ROOT / 'computations/hall_corner_exchange/source_H.json'
    source = json.loads(source_path.read_text())
    bits = int(source['source_family_hex'], 16)
    hall = {x for x in range(4096) if bits >> x & 1}
    assert len(hall) == 299 and not any(corner(hall, x, 12) for x in hall)
    repair_path = ROOT / 'damp_cornered_hall_repair_signatures.json'
    repair_data = json.loads(repair_path.read_text())
    assert repair_data['source_mathematical_profile_sha256'] == source['mathematical_profile_sha256']
    old_orders = repair_data['original_hall_repair_peelings']
    a, b, added = 183, 315, 373
    p, q = 1 << 12, 1 << 13
    seed = hall | {added}
    rows = {a: {p}, b: {0}}
    family = build(seed, rows, p)
    assert len(family) == 602
    corner_controls = [check_corner_formula(seed, rows, p, 13)]
    for tip in rows:
        assert (tip ^ added).bit_count() > 1
        remaining = hall | {tip}
        for word in old_orders[str(tip)]:
            assert word in remaining and corner(remaining, word, 12)
            remaining.remove(word)
        assert not remaining
    # These are exactly the two layer words allowed for the added seed
    # corner: both facet owners remain resolving repairs of Hall.
    deletions = []
    for t in [0, p]:
        word = added | t
        assert corner(family, word, 13)
        deleted = family - {word}
        assert {x for x in deleted if not x & p and x | p in deleted} == hall
        images = []
        for (kind, e, sign), minor in elementary_images(deleted, 13):
            order, visited = checked_peeling(minor, 13)
            images.append({'kind': kind, 'coordinate': e, 'sign': sign,
                           'order': len(minor), 'peeling': order,
                           'search_states': visited})
        assert len(images) == 39
        deletions.append({'deleted_word': word, 'order': len(deleted),
                          'strong_projection_is_Hall': True,
                          'elementary_image_peelings': images})
        print(f'Checked all 39 elementary images after deleting {word}.', flush=True)

    # The inner seed is terminal and cornered. Each of its two corner
    # deletions exposes Hall in a layer restriction, so neither descends
    # to an obstruction. Its terminality is inherited by the nested lift.
    inner = build(hall, rows, p)
    inner_corners = {x for x in inner if corner(inner, x, 13)}
    assert inner_corners == {a | p, b}
    inner_bad_minors = []
    for word in sorted(inner_corners):
        deleted = inner - {word}
        sign = int(bool(word & p))
        restricted = {x & ~p for x in deleted if bool(x & p) == bool(sign)}
        assert restricted == hall
        inner_bad_minors.append({'deleted_word': word, 'restriction_coordinate': 12,
                                 'restriction_sign': sign, 'minor_is_Hall': True})
    outer_rows = {a: {q}, b | p: {0}}
    nested = build(inner, outer_rows, q)
    corner_controls.append(check_corner_formula(inner, outer_rows, q, 14))
    flattened = build(hall, {a: set(subsets(p | q)) - {0},
                             b: set(subsets(p | q)) - {p | q}}, p | q)
    assert nested == flattened and len(nested) == 1202

    # Directly check the coordinate-wise face test for all eligibility
    # subsets of one fixed three-owner partition. This checks the exact
    # Boolean equivalence, without claiming these subsets are realized
    # by additional obstruction seeds.
    owners = {(0, 0): 0, (0, 1): 1, (1, 0): 0, (1, 1): 2}
    face_checks = []
    for eligibility in range(8):
        good = {a for a in range(3) if eligibility >> a & 1}
        allowed = [list(t) for t in itertools.product([0, 1], repeat=2)
                   if all(owners[e, t[e]] in good for e in range(2))]
        blocked = [e for e in range(2)
                   if owners[e, 0] not in good and owners[e, 1] not in good]
        assert bool(allowed) == (not blocked)
        face_checks.append({'eligible_owners': sorted(good), 'allowed_words': allowed,
                            'blocking_coordinates': blocked})
    result = {'scope': __doc__.strip(),
              'source_file_sha256': hashlib.sha256(source_path.read_bytes()).hexdigest(),
              'source_mathematical_profile_sha256': source['mathematical_profile_sha256'],
              'repair_certificate_sha256': hashlib.sha256(repair_path.read_bytes()).hexdigest(),
              'corner_controls': corner_controls, 'positive_deletions': deletions,
              'terminal_seed_bad_minors': inner_bad_minors,
              'nested_flattening_order': len(nested), 'eligibility_face_checks': face_checks,
              'all_checks_passed': True}
    output = ROOT / 'damp_seed_corner_facet_deletion_verification.json'
    output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'corner_controls': corner_controls, 'positive_image_peelings': 78,
                      'nested_flattening_order': len(nested), 'all_checks_passed': True}, indent=2))


if __name__ == '__main__':
    main()
