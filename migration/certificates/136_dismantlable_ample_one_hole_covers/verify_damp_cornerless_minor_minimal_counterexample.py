#!/usr/bin/env python3
"""Check a cornerless ample family whose sole forbidden pc-minor has a corner.

No peeling search. Replay the fixed 292-concept obstruction, 24 new
root-retaining orders, and the known 289-concept rooted factor. Their
vertex union has 580 concepts and no corners. Directly check all 72
elementary images; the all-proper-minor conclusion uses the uniform
coordinate-block theorem proved in the manuscripts.
"""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def directions(H, x, dimension):
    return sum(1 << e for e in range(dimension) if x ^ (1 << e) in H)


def corner(H, x, dimension):
    support = directions(H, x, dimension)
    sub = support
    while sub:
        if x ^ sub not in H:
            return False
        sub = (sub - 1) & support
    return True


def corners(H, dimension):
    return [x for x in sorted(H) if corner(H, x, dimension)]


def shadow(H):
    return {mask for mask in range(1 << 12)
            if len({x & mask for x in H}) == 1 << mask.bit_count()}


def check_order(H, order, dimension, root=None):
    assert len(H) == len(order) and set(order) == H
    if root is not None:
        assert order[-1] == root
    left = set(H)
    for x in order:
        assert corner(left, x, dimension)
        left.remove(x)
    assert not left


def child(H, e, kind):
    selected = H if kind == 'contraction' else {
        x for x in H if bool(x & (1 << e)) == (kind == 'opposite_section')}
    return {x & ~(1 << e) for x in selected}


def expand_dense(order, H):
    active = 0
    for x in H:
        active |= x
    coordinates = [e for e in range(active.bit_length()) if active & (1 << e)]
    return [sum(((x >> j) & 1) << e for j, e in enumerate(coordinates))
            for x in order]


def main():
    names = ['damp_fixed_291_repair_creation.json',
             'damp_hall_carrier_defect_min8.json',
             'damp_empty_core_counterexample.json',
             'computations/round776_rooted_corner_descent.json',
             'computations/round818_one_corner_rooted_images.json']
    data = {name: json.loads((ROOT / name).read_text()) for name in names}
    source = data[names[0]]
    one = data[names[2]]
    assert one['source_file_sha256'] == sha(ROOT / names[0])
    assert one['source_model_sha256'] == sha(ROOT / names[1])
    predecessor = ((set(source['base']) - set(one['removed_from_original']))
                   | set(one['added_to_original']))
    assert predecessor == set(one['predecessor']) and len(predecessor) == 291
    original = predecessor | {3959}
    assert original == set(one['counterexample']) and one['corner'] == 3959
    assert len(shadow(predecessor)) == 291 and not corners(predecessor, 12)
    assert len(shadow(original)) == 292 and corners(original, 12) == [3959]
    S = {x ^ 3959 for x in original}
    ordinary_images = {}
    record = next(r for r in one['elementary_minor_peelings']
                  if r['label'] == 'counterexample')
    seen = set()
    for row in record['minor_checks']:
        e, side = row['coordinate'], row['side']
        assert (e, side) not in seen
        seen.add((e, side))
        image = {x & ~(1 << e) for x in original
                 if side is None or (x >> e & 1) == side}
        check_order(image, row['peeling'], 12)
        kind = ('contraction' if side is None else
                'root_section' if side == (3959 >> e & 1) else 'opposite_section')
        translated = [x ^ (3959 & ~(1 << e)) for x in row['peeling']]
        check_order(child(S, e, kind), translated, 12)
        ordinary_images[e, kind] = translated
    assert seen == {(e, side) for e in range(12) for side in [None, 0, 1]}

    rooted = data[names[4]]
    assert rooted['source_sha256'] == sha(ROOT / names[2])
    assert rooted['generator_sha256'] == sha(
        ROOT / 'computations/probe_round818_one_corner_rooted_images.py')
    assert rooted['corner'] == 3959 and set(rooted['family']) == S
    S_orders = {}
    for row in rooted['rows']:
        e, kind = row['coordinate'], row['operation']
        assert (e, kind) not in S_orders
        image = child(S, e, kind)
        assert len(image) == row['order']
        check_order(image, row['rooted_peeling'], 12, root=0)
        S_orders[e, kind] = row['rooted_peeling']
    assert set(S_orders) == {(e, k) for e in range(12)
                            for k in ['root_section', 'contraction']}

    P = {x ^ 29 for x in (set(source['base']) | {91}) - {61, 125, 93}}
    ext = next(r for r in source['extensions'] if r['tip'] == 91)
    P_ordinary = [x ^ 29 for x in ext['peeling'] if x not in {61, 125, 93}]
    check_order(P, P_ordinary, 12)
    assert len(P) == 289 and corners(P, 12) == [0]
    final = next(r for r in data[names[3]]['tests']
                 if isinstance(r['deleted'], list) and set(r['deleted']) == {61, 125, 93})
    P_orders = {}
    for row in final['elementary_checks']:
        e, kind = row['coordinate'], row['operation']
        assert (e, kind) not in P_orders
        image = child(P, e, kind)
        order = expand_dense(row['peeling'], image)
        check_order(image, order, 12, root=0)
        P_orders[e, kind] = order
    assert set(P_orders) == set(S_orders)
    S_shadow, P_shadow = shadow(S), shadow(P)
    assert len(S_shadow) == 292 and len(P_shadow) == 289
    assert max(map(int.bit_count, S_shadow | P_shadow)) == 3
    assert corners(S, 12) == [0]

    X = S | {x << 12 for x in P}
    assert len(X) == 580 and not corners(X, 24)
    for e in range(12):
        for f in range(12, 24):
            assert {x & ((1 << e) | (1 << f)) for x in X} == {0, 1 << e, 1 << f}
    assert {x & 4095 for x in X} == S
    assert {x >> 12 for x in X} == P
    assert len(S_shadow) + len(P_shadow) - 1 == len(X)

    images = []
    for e in range(24):
        for kind in ['root_section', 'opposite_section', 'contraction']:
            image = child(X, e, kind)
            cs = corners(image, 24)
            assert cs and len(image) < len(X)
            if e < 12:
                if kind == 'opposite_section':
                    order = ordinary_images[e, kind]
                else:
                    order = S_orders[e, kind][:-1] + [x << 12 for x in P_ordinary]
            elif kind == 'opposite_section':
                f = e - 12
                order = [(x & ~(1 << f)) << 12 for x in P_ordinary if x & (1 << f)]
            else:
                order = None
                assert {x & 4095 for x in image} == S
                assert len(image) > len(S)
                assert P_orders[e - 12, kind][0] << 12 in cs
            if order is not None:
                check_order(image, order, 24)
            images.append({'coordinate': e, 'operation': kind, 'order': len(image),
                           'corners': cs, 'peeling': order,
                           'nonpeelable_by_S_projection': order is None})
    assert sum(r['peeling'] is not None for r in images) == 48
    assert sum(r['nonpeelable_by_S_projection'] for r in images) == 24
    report = {
        'schema': 'damp-cornerless-minor-minimal-counterexample-v1',
        'verifier_sha256': sha(Path(__file__)),
        'input_sha256': {name: sha(ROOT / name) for name in names},
        'family': sorted(X), 'order': 580, 'dimension': 24, 'vc_dimension': 3,
        'corner_count': 0,
        'minimum_degree': min(directions(X, x, 24).bit_count() for x in X),
        'S': sorted(S), 'P': sorted(P),
        'S_rooted_minor_orders': [dict(coordinate=e, operation=k, peeling=o)
                                 for (e, k), o in S_orders.items()],
        'P_rooted_minor_orders': [dict(coordinate=e, operation=k, peeling=o)
                                 for (e, k), o in P_orders.items()],
        'P_ordinary_order': P_ordinary,
        'elementary_images': images,
        'all_checks_passed': True,
        'scope': 'Fixed hypotheses and all 72 elementary images checked directly. '
                 'All proper nonempty pc-minors having a corner and the unique '
                 'forbidden pc-minor S follow from the uniform block-minor proof, '
                 'not from extrapolating the elementary corner checks.'}
    (ROOT / 'damp_cornerless_minor_minimal_counterexample_verification.json').write_text(
        json.dumps(report, indent=2) + '\n')
    print({k: report[k] for k in ['all_checks_passed', 'order', 'dimension',
                                'vc_dimension', 'corner_count', 'minimum_degree']})
    print('Rooted orders: 48; old ordinary minor orders: 36; union images: 48 peelable, 24 negative.')


if __name__ == '__main__':
    main()
