#!/usr/bin/env python3
"""Verify inputs and fixed realizations of reduced-base universality.

The theorem is uniform in a nonempty peelable ample input K. Check the
fixed 299 seed, all 36 seed-minor orders, three repairs, and four K inputs.
Materialize the outputs for a singleton and an edge only. The 289 input
uses the exact clause formulas and compact row data; its output is not
materialized. No source search, peeling search, or output-minor census.
"""
from itertools import combinations
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


def corner(H, x, n):
    support = sum(1 << q for q in range(n) if x ^ (1 << q) in H)
    return all(x ^ t in H for t in subsets(support))


def check_order(H, order, n):
    assert len(order) == len(H) and set(order) == H
    left = set(H)
    for x in order:
        assert corner(left, x, n), (len(left), x)
        left.remove(x)


def shadow(H, n):
    answer = set()
    for x in H:
        def extend(mask, vertices, start):
            answer.add(mask)
            for q in range(start, n):
                bit = 1 << q
                if x & bit:
                    continue
                other = [v | bit for v in vertices]
                if all(v in H for v in other):
                    extend(mask | bit, vertices + other, q + 1)
        extend(0, [x], 0)
    assert len(answer) == len(H)
    return answer


def clauses(H, sigma, n):
    boundary = {s | (1 << q) for s in sigma for q in range(n)} - sigma
    minimal = sorted(s for s in boundary
                     if all(s ^ (1 << q) in sigma for q in range(n) if s >> q & 1))
    result = {}
    for s in minimal:
        missing = set(subsets(s)) - {v & s for v in H}
        assert len(missing) == 1
        result[s] = missing.pop()
    return result


def minor(H, e, operation):
    b = 1 << e
    if operation == 'contraction':
        return {x & ~b for x in H}
    assert operation in ('root_section', 'opposite_section')
    return {x & ~b for x in H
            if (x >> e) & 1 == int(operation == 'opposite_section')}


def main():
    names = ['computations/hall_corner_exchange/source_H.json',
             'damp_terminal_clause_recovery_verification.json',
             'computations/round810_hall_complement_peeling.json',
             'damp_rooted_clause_pivots_verification.json']
    source, repairs, seed_orders, rooted = [json.loads((ROOT / p).read_text()) for p in names]
    hashes = {p: sha(ROOT / p) for p in names}
    assert seed_orders['source_sha256'] == hashes[names[0]]
    assert repairs['source_file_sha256'] == hashes[names[0]]
    bits = int(source['source_family_hex'], 16)
    S = {x for x in range(4096) if bits >> x & 1}
    sigma = shadow(S, 12)
    assert len(S) == 299 and sigma == {s for s in range(4096) if s.bit_count() <= 3}
    assert not any(corner(S, x, 12) for x in S)
    seen = set()
    for row in seed_orders['seed_minor_peelings']:
        key = row['coordinate'], row['operation']
        assert key not in seen
        seen.add(key)
        M = minor(S, *key)
        assert 0 < len(M) < len(S)
        check_order(M, row['peeling'], 12)
    assert seen == {(q, op) for q in range(12)
                    for op in ('root_section', 'opposite_section', 'contraction')}
    tips = [183, 315, 1836]
    orders = {r['tip']: r['peeling'] for r in repairs['repair_proofs']}
    supports = []
    for a in tips:
        assert a not in S
        p = sum(1 << q for q in range(12) if a ^ (1 << q) in S)
        assert p.bit_count() == 4 and all(a ^ t in S for t in subsets(p) if t)
        assert shadow(S | {a}, 12) == sigma | {p}
        check_order(S | {a}, orders[a], 12)
        supports.append(p)
    assert len(set(supports)) == 3
    seed_clauses = clauses(S, sigma, 12)
    old_clauses = {s: t for s, t in seed_clauses.items() if s not in supports}
    assert len(old_clauses) == 492
    old_patterns = [tuple(bool(s >> q & 1) for s in old_clauses) for q in range(12)]
    assert len(set(old_patterns)) == 12 and all(any(p) for p in old_patterns)

    base = rooted['square_lift_tower'][0]
    inputs = [('singleton', {0}, 0, [0]), ('edge', {0, 1}, 1, [0, 1]),
              ('ordinarily_peelable_rooted_factor', set(base['negative_family']),
               12, base['negative_ordinary_peeling']),
              ('cornerless_negative_input', {x ^ min(S) for x in S}, 12, None)]
    reports = []
    for name, K, n, k_order in inputs:
        assert 0 in K
        if k_order is not None:
            check_order(K, k_order, n)
        else:
            assert not any(corner(K, x, n) for x in K)
        sk = shadow(K, n)
        assert all(1 << q in sk for q in range(n))
        ck = clauses(K, sk, n)
        w, p, q, z = [1 << (n + j) for j in range(4)]
        m = n + 4
        full = (1 << m) - 1
        U = K | {x | w for x in K}
        L = U | {p, q, p | q}
        assert len(L) == 2 * len(K) + 3
        if k_order is not None:
            l_order = [p | q, p, q] + [x | t for x in k_order for t in (0, w)]
            check_order(L, l_order, m)
        sl = shadow(L, m)
        c0 = clauses(L, sl, m)  # L is the row R_0, with z fixed at zero.
        expected = dict(ck)
        expected.update({(1 << e) | t: (1 << e) | t
                         for e in range(n + 1) for t in (p, q)})
        expected[z] = z
        assert c0 == expected
        signed = dict(old_clauses)
        for J, trace in c0.items():
            signed[supports[0] | (J << 12)] = (tips[0] & supports[0]) | (trace << 12)
        signed[supports[1] | (full << 12)] = tips[1] & supports[1]
        signed[supports[2] | (full << 12)] = (tips[2] & supports[2]) | (full << 12)
        assert len(signed) == 497 + 2 * n + len(ck)
        # Different support columns imply different signed columns up to sign.
        patterns = [tuple(bool(s >> e & 1) for s in signed) for e in range(12 + m)]
        assert all(any(v) for v in patterns) and len(set(patterns)) == 12 + m
        cardinality = 301 * (1 << m) + 2 * len(K) + 1
        assert max(3 + m, 4 + max(s.bit_count() for s in sl), 4 + m - 1) == m + 3
        output = dict(order=cardinality, dimension=12 + m, vc_dimension=m + 3,
                      dimension_difference=9, minimum_degree_lower_bound=4,
                      distinct_support_columns=len(patterns), clause_count=len(signed),
                      complete_crossing_graph=True, output_minor_orders_replayed=0,
                      forbidden_by_uniform_theorem=(k_order is not None),
                      nonpeelable_by_fixed_seed_derivative=True)
        if n <= 1:
            layers = set(range(1 << m))
            rows = [L, layers - {0}, layers - {full}]
            check_order(rows[1], sorted(rows[1], key=lambda x: (x.bit_count(), x)), m)
            check_order(rows[2], sorted(rows[2], key=lambda x: (-x.bit_count(), x)), m)
            for e in range(m):
                for side in (0, 1):
                    facet = {t for t in layers if (t >> e) & 1 == side}
                    assert facet <= rows[1 if side else 2]
            X = {x | (t << 12) for x in S for t in layers}
            for a, row in zip(tips, rows):
                X.update(a | (t << 12) for t in row)
            assert len(X) == cardinality
            sx = shadow(X, 12 + m)
            assert max(s.bit_count() for s in sx) == m + 3
            assert all((1 << e) | (1 << f) in sx for e, f in combinations(range(12 + m), 2))
            actual = clauses(X, sx, 12 + m)
            assert actual == signed
            degrees = [sum(x ^ (1 << e) in X for e in range(12 + m)) for x in X]
            assert min(degrees) >= 4
            # Properness of the three rows makes the full layer derivative exactly S.
            derivative = {a for a in {x & 4095 for x in X}
                          if all(a | (t << 12) in X for t in layers)}
            assert derivative == S
            recovered = {t for t in layers if tips[0] | (t << 12) in X
                         and not t & (w | p | q | z)}
            assert recovered == K
            output.update(materialized=True, actual_minimum_degree=min(degrees),
                          corner_count=sum(corner(X, x, 12 + m) for x in X),
                          cube_shadow_checked=True, full_signed_presentation_checked=True,
                          strong_projection_equals_seed=True, coordinate_face_equals_input=True)
        else:
            # The same membership formula checks the whole prescribed K face,
            # without forming the exponentially larger forbidden output.
            assert {t for t in range(1 << n) if t in L} == K
            output.update(materialized=False, coordinate_face_equals_input=True,
                          justification='Uniform theorem and compact input checks; no output enumeration.')
        reports.append(dict(name=name, input_order=len(K), input_dimension=n,
                            row_zero_order=len(L), input_peeling_replayed=(k_order is not None),
                            cornerless_input_checked=(k_order is None),
                            row_zero_peeling_replayed=(k_order is not None), output=output))
        print(json.dumps(reports[-1]), flush=True)
    result = dict(schema=1, all_checks_passed=True, verifier_sha256=sha(Path(__file__)),
                  input_sha256=hashes, seed_order=299, seed_minor_orders_replayed=36,
                  seed_repair_orders_replayed=3, models=reports,
                  scope=('Uniform universality theorem for reduced forbidden bases. Finite input '
                         'and two full-output checks only; no output-minor orders or general '
                         'intrinsic base-admissibility criterion are claimed.'))
    target = ROOT / 'damp_reduced_base_universality_verification.json'
    target.write_text(json.dumps(result, indent=2) + '\n')
    print('All checks passed:', target, flush=True)


if __name__ == '__main__':
    main()
