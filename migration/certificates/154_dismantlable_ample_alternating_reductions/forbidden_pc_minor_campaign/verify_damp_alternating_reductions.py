#!/usr/bin/env python3
"""Verify forced alternation of graph pruning and signed-column compression.

No source search and no peeling search. Rebuild the fixed 577/579 families,
replace coordinate 0 by {0,24}, then complete one specified square. Check
actual cubes, every proper-minor order, negative certificates, both column
representatives, both pruning priorities, and intrinsic reconstruction data.
The all-size normal-form and phase-bound statements are symbolic proofs.
"""
from itertools import product
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def subsets(mask):
    t = mask
    while True:
        yield t
        if not t:
            return
        t = (t - 1) & mask


def directions(H, x, ground):
    return [q for q in ground if x ^ (1 << q) in H]


def corner(H, x, ground):
    support = sum(1 << q for q in directions(H, x, ground))
    return all(x ^ t in H for t in subsets(support))


def check_order(H, order, ground):
    assert len(order) == len(H) and set(order) == H
    left = set(H)
    for x in order:
        assert corner(left, x, ground), (len(left), x)
        left.remove(x)


def cube_shadow(H, ground):
    answer = set()
    ground = sorted(ground)
    for x in H:
        def extend(mask, vertices, start):
            answer.add(mask)
            for j in range(start, len(ground)):
                bit = 1 << ground[j]
                if x & bit:
                    continue
                other = [v | bit for v in vertices]
                if all(v in H for v in other):
                    extend(mask | bit, vertices + other, j + 1)
        extend(0, [x], 0)
    assert len(answer) == len(H)
    return answer


def columns(H, ground, shadow=None):
    if shadow is None:
        shadow = cube_shadow(H, ground)
    # Minimal nonfaces lie immediately above shadow faces. Avoid enumerating Q25.
    boundary = {s | (1 << q) for s in shadow for q in ground} - shadow
    supports = sorted(s for s in boundary
                      if all(s ^ (1 << q) in shadow for q in ground if s >> q & 1))
    signed = []
    for s in supports:
        absent = set(subsets(s)) - {x & s for x in H}
        assert len(absent) == 1
        signed.append((s, absent.pop()))
    profiles = {q: tuple(0 if not s >> q & 1 else 1 if t >> q & 1 else -1
                         for s, t in signed) for q in ground}
    groups = {}
    for q in sorted(ground):
        v = profiles[q]
        key = min(v, tuple(-t for t in v))
        groups.setdefault(key, []).append(q)
    return signed, list(groups.values()), profiles


def strong(H, mask):
    return {x for x in H if not x & mask and all(x | t in H for t in subsets(mask))}


def minor(H, q, side):
    return {x & ~(1 << q) for x in H if side is None or (x >> q) & 1 == side}


def replace(H, q, new):
    a, b = 1 << q, 1 << new
    B0, B1 = minor(H, q, 0), minor(H, q, 1)
    C = B0 | B1
    return B0 | {x | a | b for x in B1} | {x | a for x in C} | {x | b for x in C}


def lift_order(H, order, q, new):
    a, b = 1 << q, 1 << new
    left, result = set(H), []
    for x in order:
        y = x & ~a
        if x ^ a in left:
            result.append(y | (a | b if x & a else 0))
        else:
            result.extend([y | a, y | b, y | (a | b if x & a else 0)])
        left.remove(x)
    return result


def negative_states(H, states, ground):
    rows = {frozenset(r['removed']): r for r in states}
    assert frozenset() in rows
    for removed, row in rows.items():
        assert removed <= H and H - removed
        cs = sorted(x for x in H - removed if corner(H - removed, x, ground))
        assert cs == sorted(row['corners'])
        assert all(removed | {x} in rows for x in cs)


def prune(H, ground, reverse=False):
    left, removed = set(H), []
    while True:
        ready = sorted(x for x in left if len(directions(left, x, ground)) < 3)
        if not ready:
            return left, removed
        x = ready[-1] if reverse else ready[0]
        assert len(directions(left, x, ground)) == 2 and corner(left, x, ground)
        removed.append(x); left.remove(x)


def compress(H, ground, groups, profiles, reverse=False):
    flip = 0
    for group in groups:
        ref = profiles[group[0]]
        for q in group:
            if profiles[q] != ref:
                assert profiles[q] == tuple(-t for t in ref)
                flip |= 1 << q
    normalized = {x ^ flip for x in H}
    reps = [g[-1] if reverse else g[0] for g in groups]
    deleted = sum(1 << q for q in ground if q not in reps)
    small = strong(normalized, deleted)
    # Rename representatives to the first original coordinate of each group.
    result = {sum(1 << g[0] for g, q in zip(groups, reps) if x >> q & 1) for x in small}
    return result, {g[0] for g in groups}, dict(groups=groups, representatives=reps, flip_mask=flip)


def recover_squares(core, upper, removed, ground):
    names = {x: 'a' + str(i) for i, x in enumerate(sorted(upper - core))}
    def ref(x):
        return ['new', names[x]] if x in names else ['base', x]
    left, rows = set(upper), []
    for x in removed:
        p, q = sorted(directions(left, x, ground))
        u, v, w = x ^ (1 << p), x ^ (1 << q), x ^ (1 << p) ^ (1 << q)
        assert {u, v, w} <= left
        rows.append(dict(name=names[x], p=p, q=q, u=ref(u), v=ref(v), w=ref(w)))
        left.remove(x)
    assert left == core
    # The reconstruction uses abstract labels/targets, not the removed words.
    for reverse in (False, True):
        values, pending, rebuilt = {}, {r['name']: r for r in rows}, set(core)
        old_shadow = cube_shadow(core, ground)
        masks = [(1 << r['p']) | (1 << r['q']) for r in rows]
        assert len(set(masks)) == len(masks) and not set(masks) & old_shadow
        while pending:
            ready = sorted(k for k, r in pending.items()
                           if all(t[0] == 'base' or t[1] in values for t in (r['u'], r['v'], r['w'])))
            assert ready
            name = ready[-1] if reverse else ready[0]
            r = pending.pop(name)
            def value(t):
                if t[0] == 'base':
                    assert t[1] in core
                    return t[1]
                return values[t[1]]
            x = value(r['u']) ^ (1 << r['p'])
            assert x ^ (1 << r['q']) == value(r['v'])
            assert x ^ (1 << r['p']) ^ (1 << r['q']) == value(r['w'])
            assert x not in rebuilt
            values[name] = x; rebuilt.add(x)
        assert rebuilt == upper
    return rows


def reduction(H, ground, reverse=False):
    H, ground = set(H), set(ground)
    stages = []
    deficit = None
    while True:
        shadow = cube_shadow(H, ground)
        vc = max(s.bit_count() for s in shadow)
        if deficit is None:
            deficit = len(ground) - vc
        assert len(ground) - vc == deficit
        signed, groups, profiles = columns(H, ground, shadow)
        low = sorted(x for x in H if len(directions(H, x, ground)) < 3)
        repeated = [g for g in groups if len(g) > 1]
        assert not (low and repeated)
        if not low and not repeated:
            return H, ground, stages
        old_order = len(H)
        if low:
            core, removed = prune(H, ground, reverse)
            data = recover_squares(core, H, removed, ground)
            stages.append(dict(kind='square', source_order=old_order, target_order=len(core),
                               dimension=len(ground), vc_dimension=vc, removed=removed,
                               attachment_data=data))
            H = core
        else:
            small, remaining, data = compress(H, ground, groups, profiles, reverse)
            stages.append(dict(kind='columns', source_order=old_order, target_order=len(small),
                               source_dimension=len(ground), target_dimension=len(remaining),
                               source_vc_dimension=vc, column_data=data))
            H, ground = small, remaining
        assert len(H) < old_order
        assert len(stages) < 20, 'unexpected phase count'
        assert len(stages) < 2 or stages[-1]['kind'] != stages[-2]['kind']


def main():
    names = ['damp_three_core_verification.json', 'damp_noncorner_rooted_obstruction_verification.json']
    data = {name: json.loads((ROOT/name).read_text()) for name in names}
    first = data[names[0]]
    saved = first['two_stage_example']
    K, A = set(first['core_family']), set(saved['family'])
    old = next(r for r in data[names[1]]['doubles'] if r['order'] == 577)
    assert set(old['family']) == K
    assert A == K | {8194, 8198}
    ground = set(range(24))
    for H in (K, A):
        sigma = cube_shadow(H, ground)
        assert max(s.bit_count() for s in sigma) == 3
        assert all(len(g) == 1 for g in columns(H, ground, sigma)[1])
    assert not any(corner(K, x, ground) for x in K)
    negative_states(A, saved['negative_certificate'], ground)
    orders = {}
    checked = 0
    for label, H, rows in [('K', K, old['elementary_minor_peelings']), ('A', A, saved['elementary_minor_peelings'])]:
        local = {}
        for row in rows:
            side = {'contraction': None, 'root_section': 0, 'opposite_section': 1}[row['operation']]
            q, order = row['coordinate'], row['peeling']
            assert (q, side) not in local
            check_order(minor(H, q, side), order, ground)
            local[q, side] = order; checked += 1
        assert len(local) == 72
        orders[label] = local
    X = replace(A, 0, 24)
    # This word completes the fixed missing pair {1,14} at the constant 00 section.
    tip = 24578
    assert tip not in X
    Y = X | {tip}
    large_ground = set(range(25))
    assert directions(Y, tip, large_ground) == [1, 14]
    for H in (X, Y):
        sigma = cube_shadow(H, large_ground)
        assert max(s.bit_count() for s in sigma) == 4
        assert strong(H, 1 << 24) == A
    assert len(X) == 3 * len(minor(A, 0, None)) + len(strong(A, 1))
    output_orders = []
    for q in range(25):
        for side in (None, 0, 1):
            child = minor(X, q, side)
            if q not in (0, 24):
                order = lift_order(minor(A, q, side), orders['A'][q, side], 0, 24)
            else:
                keep = 24 if q == 0 else 0
                bit = 1 << keep
                if side is None:
                    order = [x | t for x in orders['A'][0, None] for t in [0, bit]]
                else:
                    # The restricted family is a peripheral expansion of C along B_side.
                    order = [x | (bit if side else 0) for x in orders['A'][0, side]]
                    order += [x | (0 if side else bit) for x in orders['A'][0, None]]
            check_order(child, order, large_ground); checked += 1
            bigger = minor(Y, q, side)
            added = sorted(bigger - child)
            assert len(added) <= 1
            final_order = added + order
            check_order(bigger, final_order, large_ground); checked += 1
            output_orders.append(dict(coordinate=q, side=side, X_order=order, Y_order=final_order))
    reduced = []
    for reverse in (False, True):
        endpoint, remaining, stages = reduction(Y, large_ground, reverse)
        assert endpoint == K and remaining == ground
        assert [s['kind'] for s in stages] == ['square', 'columns', 'square']
        assert [s['source_order'] for s in stages] == [len(Y), len(X), 579]
        assert [s['target_order'] for s in stages] == [len(X), 579, 577]
        assert len(stages) == 2 * 4 - 5
        reduced.append(dict(reverse_choices=reverse, stages=stages))
    # One application of each map is not commutative on the intermediate X.
    assert prune(X, large_ground)[0] == X
    _, groups, profiles = columns(X, large_ground)
    compressed, remaining, _ = compress(X, large_ground, groups, profiles)
    assert compressed == A and prune(compressed, remaining)[0] == K
    # Scope control: mutual exclusion is false for arbitrary peelable ample families.
    square = set(range(4))
    assert all(len(directions(square, x, {0, 1})) == 2 for x in square)
    assert columns(square, {0, 1})[1] == [[0, 1]]
    # Stopping conditions alone do not imply forbiddenness.
    positive = {x for x in range(64) if x.bit_count() <= 3}
    assert len(positive) == 42
    positive_shadow = cube_shadow(positive, set(range(6)))
    assert positive_shadow == positive
    assert min(len(directions(positive, x, set(range(6)))) for x in positive) == 3
    assert all(len(g) == 1 for g in columns(positive, set(range(6)), positive_shadow)[1])
    positive_order = sorted(positive, key=lambda x: (-x.bit_count(), x))
    check_order(positive, positive_order, set(range(6)))
    result = dict(schema='damp-alternating-reductions-v1', verifier_sha256=sha(Path(__file__)),
                  input_sha256={name: sha(ROOT/name) for name in names},
                  source_order=len(Y), intermediate_order=len(X), final_order=len(K),
                  source_dimension=25, source_vc_dimension=4, dimension_minus_vc=21,
                  proper_minor_orders_checked=checked, X=sorted(X), Y=sorted(Y),
                  added_word=tip, all_elementary_output_orders=output_orders,
                  reduction_checks=reduced, noncommuting_one_pass_orders=[579, 577],
                  peelable_fixed_point_control=dict(order=42, dimension=6, vc_dimension=3,
                                                    family=sorted(positive), peeling=positive_order),
                  bound_attained=True, all_checks_passed=True,
                  scope='Fixed forbidden examples and all proper-minor certificates. The all-size forced-alternation, uniqueness, and intrinsic reconstruction statements are proved symbolically.')
    (ROOT/'damp_alternating_reductions_verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print('All checks passed:',len(Y),'->',len(X),'-> 579 -> 577;',checked,'proper-minor orders;',
          'both reduction priorities and representatives; three phases attain the VC-four bound.',flush=True)


if __name__ == '__main__':
    main()
