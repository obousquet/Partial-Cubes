#!/usr/bin/env python3
"""Check coordinate-block replacement, using actual cubes and corner orders.

Exhaust all ample families through Q3, all marked coordinates, and blocks
of sizes two and three. Check every corner lift, the shadow, clauses,
strong projections, and old-coordinate minor commutation. Independently
replay the fixed 291 input's minor certificates and construct all proper
minor orders for its 739- and 1635-concept replacements. No negative
search is used: the input and both outputs are checked cornerless.
The uniform all-size conclusions are proved in the manuscript.
"""
from itertools import product
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def submasks(mask):
    t = mask
    while True:
        yield t
        if not t:
            return
        t = (t - 1) & mask


def corner(H, x, n):
    support = sum(1 << i for i in range(n) if x ^ (1 << i) in H)
    return all(x ^ t in H for t in submasks(support))


def check_order(H, order, n):
    assert len(order) == len(H) and set(order) == H
    H = set(H)
    for x in order:
        assert corner(H, x, n), (len(H), x)
        H.remove(x)


def shadow(H, n):
    result = set()
    for x in H:
        def extend(support, vertices, start):
            result.add(support)
            for i in range(start, n):
                bit = 1 << i
                if x & bit:
                    continue
                upper = [y | bit for y in vertices]
                if all(y in H for y in upper):
                    extend(support | bit, vertices + upper, i + 1)
        extend(0, [x], 0)
    return result


def strong(H, mask):
    return {x for x in H if not x & mask
            and all(x | t in H for t in submasks(mask))}


def minor(H, q, side):
    return {x & ~(1 << q) for x in H
            if side is None or (x >> q) & 1 == side}


def replace(H, z, block):
    """Keep old bit positions; the old marked bit can belong to block."""
    B0, B1 = minor(H, z, 0), minor(H, z, 1)
    C = B0 | B1
    result = {x for x in B0} | {x | block for x in B1}
    for t in submasks(block):
        if t and t != block:
            result.update(x | t for x in C)
    return result


def corner_lift(H, x, z, new):
    """Deletion schedule for a single source corner, for block {z,new}."""
    a, b = 1 << z, 1 << new
    base = x & ~a
    if x ^ a in H:
        return [base | (a | b if x & a else 0)]
    return [base | a, base | b, base | (a | b if x & a else 0)]


def lift_order(H, order, z, new):
    H = set(H)
    result = []
    for x in order:
        result.extend(corner_lift(H, x, z, new))
        H.remove(x)
    return result


def clauses(H, n, sigma, ground=None):
    result = []
    if ground is None:
        ground = (1 << n) - 1
    for support in submasks(ground):
        if support in sigma or any(support ^ (1 << i) not in sigma
                                  for i in range(n) if support >> i & 1):
            continue
        absent = set(submasks(support)) - {x & support for x in H}
        assert len(absent) == 1
        result.append((support, absent.pop()))
    return result


def column_reduction(H, n):
    """Normalize column signs, compress, rebuild, and check all representatives."""
    signed = clauses(H, n, shadow(H, n))
    classes = {}
    flip = 0
    for q in range(n):
        v = tuple(0 if not s >> q & 1 else (1 if t >> q & 1 else -1)
                  for s, t in signed)
        key = min(v, tuple(-a for a in v))
        classes.setdefault(key, []).append(q)
        if v != key:
            flip |= 1 << q
    blocks = list(classes.values())
    normalized = {x ^ flip for x in H}
    ground = sum(1 << b[0] for b in blocks)
    deleted = ((1 << n) - 1) ^ ground
    reduced = strong(normalized, deleted)
    assert len(shadow(reduced, n)) == len(reduced)
    rebuilt = set(reduced)
    for b in blocks:
        if len(b) > 1:
            rebuilt = replace(rebuilt, b[0], sum(1 << q for q in b))
    assert rebuilt == normalized
    reduced_signed = clauses(reduced, n, shadow(reduced, n), ground)
    assert len(reduced_signed) == len(signed)
    signatures = []
    for b in blocks:
        q = b[0]
        v = tuple(0 if not s >> q & 1 else (1 if t >> q & 1 else -1)
                  for s, t in reduced_signed)
        signatures.append(min(v, tuple(-a for a in v)))
    assert len(signatures) == len(set(signatures))
    choices = 0
    for reps in product(*blocks):
        kept = sum(1 << q for q in reps)
        other = strong(normalized, ((1 << n) - 1) ^ kept)
        renamed = {sum(1 << b[0] for b, q in zip(blocks, reps) if x >> q & 1)
                   for x in other}
        assert renamed == reduced
        choices += 1
    return dict(blocks=blocks, reduced_order=len(reduced),
                reduced_dimension=len(blocks), representative_choices_checked=choices)


def toy_checks():
    counts = dict(ample_inputs=0, replacements=0, corner_lifts=0,
                  strong_identities=0, minor_identities=0, clause_identities=0,
                  column_normal_forms=0, representative_choices=0)
    for n in range(1, 4):
        for bits in range(1 << (1 << n)):
            H = {x for x in range(1 << n) if bits >> x & 1}
            sigma = shadow(H, n)
            if len(sigma) != len(H):
                continue
            counts['ample_inputs'] += 1
            reduced = column_reduction(H, n)
            counts['column_normal_forms'] += 1
            counts['representative_choices'] += reduced['representative_choices_checked']
            for z in range(n):
                old = ((1 << n) - 1) ^ (1 << z)
                C, A = minor(H, z, None), strong(H, 1 << z)
                for m in (2, 3):
                    block = (1 << z) | (((1 << (m - 1)) - 1) << n)
                    X = replace(H, z, block)
                    expected = {s | t for s in shadow(C, n) for t in submasks(block)
                                if t != block}
                    expected |= {s | block for s in shadow(A, n)}
                    assert shadow(X, n + m - 1) == expected
                    assert len(expected) == len(X)
                    wanted_clauses = {(s if not s >> z & 1 else (s ^ (1 << z)) | block,
                                       t if not s >> z & 1 else
                                       (t & ~(1 << z)) | (block if t >> z & 1 else 0))
                                      for s, t in clauses(H, n, sigma)}
                    assert set(clauses(X, n + m - 1, expected)) == wanted_clauses
                    counts['clause_identities'] += 1
                    for D in submasks(old):
                        for T in submasks(block):
                            target = (strong(H, D | (1 << z)) if T == block else
                                      replace(strong(H, D), z, block ^ T))
                            assert strong(X, D | T) == target
                            counts['strong_identities'] += 1
                    for q in range(n):
                        if q == z:
                            continue
                        for side in (None, 0, 1):
                            assert minor(X, q, side) == replace(minor(H, q, side), z, block)
                            counts['minor_identities'] += 1
                    if m == 2:
                        for x in H:
                            if not corner(H, x, n):
                                continue
                            left = set(X)
                            for y in corner_lift(H, x, z, n):
                                assert corner(left, y, n + 1)
                                left.remove(y)
                            assert left == replace(H - {x}, z, block)
                            counts['corner_lifts'] += 1
                    counts['replacements'] += 1
    return counts


def fixed_checks():
    source = json.loads((ROOT / 'damp_fixed_291_repair_creation.json').read_text())
    records = json.loads((ROOT / 'computations/round766_base_minor_peelings.json').read_text())
    H = set(source['base'])
    assert len(H) == 291 and len(shadow(H, 12)) == 291
    assert max(s.bit_count() for s in shadow(H, 12)) == 3
    assert not any(corner(H, x, 12) for x in H)
    assert column_reduction(H, 12)['reduced_dimension'] == 12
    orders = {(r['coordinate'], r['side']): r['peeling'] for r in records}
    assert len(orders) == 36
    for (q, side), order in orders.items():
        check_order(minor(H, q, side), order, 12)
    # Construct positive carrier orders; failure is unresolved, never negative.
    carriers = []
    for q in range(12):
        K = strong(H, 1 << q)
        left = set(K); order = []
        while left:
            choices = [x for x in sorted(left) if corner(left, x, 12)]
            assert choices, ('positive carrier order unresolved', q)
            order.append(choices[0]); left.remove(choices[0])
        check_order(K, order, 12)
        carriers.append(dict(coordinate=q, order=len(K), peeling=order))

    output = []
    z = 0
    for m in (2, 3):
        block = (1 << z) | (((1 << (m - 1)) - 1) << 12)
        X = replace(H, z, block)
        n = 11 + m
        sigma = shadow(X, n)
        assert len(sigma) == len(X)
        assert max(s.bit_count() for s in sigma) == m + 2
        assert not any(corner(X, x, n) for x in X)
        degrees = [sum(x ^ (1 << i) in X for i in range(n)) for x in X]
        assert min(degrees) >= 3
        reduction = column_reduction(X, n)
        assert reduction['reduced_order'] == 291 and reduction['reduced_dimension'] == 12
        assert [b for b in reduction['blocks'] if len(b) > 1] == [[i for i in range(n) if block >> i & 1]]
        for e in range(n):
            if block >> e & 1:
                assert all((1 << e) | (1 << q) in sigma for q in range(n) if q != e)
        original_clauses = clauses(H, 12, shadow(H, 12))
        expected_clauses = {(s if not s >> z & 1 else (s ^ (1 << z)) | block,
                             t if not s >> z & 1 else
                             (t & ~(1 << z)) | (block if t >> z & 1 else 0))
                            for s, t in original_clauses}
        assert set(clauses(X, n, sigma)) == expected_clauses
        minor_rows = []
        for q in range(n):
            for side in (None, 0, 1):
                child = minor(X, q, side)
                if not block >> q & 1:
                    small = minor(H, q, side)
                    order = orders[q, side]
                    for j in range(12, n):
                        order = lift_order(small, order, z, j)
                        small = replace(small, z, (1 << z) | (1 << j))
                    assert small == child
                else:
                    C = minor(H, z, None)
                    remaining = block ^ (1 << q)
                    if side is None:
                        order = [x | t for x in orders[z, None] for t in submasks(remaining)]
                    else:
                        B = minor(H, z, side)
                        # Peel the exceptional site B first, then the product
                        # of C with the cube punctured at that endpoint.
                        endpoint = remaining if side else 0
                        row = set(submasks(remaining)) - {endpoint}
                        layers = []
                        while row:
                            t = next(t for t in sorted(row) if corner(row, t, n))
                            layers.append(t); row.remove(t)
                        order = [x | endpoint for x in orders[z, side]]
                        order += [x | t for x in orders[z, None] for t in layers]
                check_order(child, order, n)
                minor_rows.append(dict(coordinate=q, side=side, order=len(child), peeling=order))
        negative = []
        for D in submasks(block):
            K = strong(X, D)
            expected = strong(H, 1 << z) if D == block else replace(H, z, block ^ D)
            assert K == expected
            if D != block:
                assert not any(corner(K, x, n) for x in K)
                negative.append(D)
        positive_rows = []
        for row in carriers:
            q = row['coordinate']
            if q == z:
                mask = block
                order = row['peeling']
            else:
                mask = 1 << q
                small = strong(H, mask)
                order = row['peeling']
                for j in range(12, n):
                    order = lift_order(small, order, z, j)
                    small = replace(small, z, (1 << z) | (1 << j))
                assert small == strong(X, mask)
            check_order(strong(X, mask), order, n)
            positive_rows.append(dict(support=mask, order=len(order), peeling=order))
        output.append(dict(block_size=m, block_coordinates=[i for i in range(n) if block >> i & 1],
                           family=sorted(X), order=len(X), dimension=n, vc_dimension=m+2,
                           minimum_degree=min(degrees), cornerless=True,
                           signed_column_reduction=reduction,
                           negative_block_supports=sorted(negative),
                           full_block_carrier_order=len(strong(X, block)),
                           positive_carrier_peelings=positive_rows,
                           elementary_minor_peelings=minor_rows))
        print('fixed replacement',m,len(X),'all',len(minor_rows),'proper-minor orders passed',flush=True)
    return carriers, output


def main():
    toy = toy_checks()
    print('toy checks', toy, flush=True)
    carriers, fixed = fixed_checks()
    inputs = ['damp_fixed_291_repair_creation.json', 'computations/round766_base_minor_peelings.json']
    result = dict(schema='damp-coordinate-block-replacement-v1', verifier_sha256=sha(Path(__file__)),
                  input_sha256={name: sha(ROOT / name) for name in inputs},
                  toy_checks=toy, base_carrier_peelings=carriers, fixed_replacements=fixed,
                  proper_minor_orders_checked=36 + sum(len(x['elementary_minor_peelings']) for x in fixed),
                  positive_carrier_orders_checked=len(carriers) + sum(len(x['positive_carrier_peelings']) for x in fixed),
                  all_checks_passed=True,
                  scope='Finite checks support the uniform symbolic proof. They do not prove exhaustive generation or confluence up to isomorphism.')
    (ROOT / 'damp_coordinate_block_replacement_verification.json').write_text(json.dumps(result, indent=2)+'\n')
    print('All checks passed.', flush=True)


if __name__ == '__main__':
    main()
