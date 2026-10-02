#!/usr/bin/env python3
"""Rooted clause pivots, an antipodal nonrepair, and remote endpoint changes.

Replay a fixed 290-concept negative extension and its 24 rooted minors.
Its two nonroot corners are antipodes of the completed four-cube. Check
the exterior-neighbor witnesses in the uniform proof, and compare with
the restoring tip 604. Check two iterations of the uniform square lift,
including all rooted minors and the single changed clause. No search.
"""
import hashlib
import itertools
import json
from pathlib import Path

from verify_damp_cut_vertex_obstruction import check_order, expand_dense
from verify_damp_layer_midpoint_cores import shadow, subsets
from verify_damp_paired_root_transfer import check_negative, square
from verify_damp_rank_three_corner_neighbors import corners, directions
from verify_damp_terminal_cut_vertex import elementary

ROOT = Path(__file__).resolve().parent


def lift(C, ordinary, alpha, delta, a_order, d_order, dimension,
         rooted_orders=None, endpoint_order=None):
    """The two corner-tip repairs propagate to their missing 11 copies."""
    A, D = C | {alpha}, C | {delta}
    sc, sa, sd, su = [shadow(F, dimension) for F in [C, A, D, A | D]]
    pa, pd = directions(C, alpha, dimension), directions(C, delta, dimension)
    assert sa == sc | {pa} and sd == sc | {pd} and pa != pd
    assert su == sc | {pa, pd} and all(len(s) == len(F) for s, F in
                                      [(sc, C), (sa, A), (sd, D), (su, A | D)])
    check_order(A, a_order, dimension, last=0)
    check_order(D, d_order, dimension, last=0)
    Y = square(C, A, D)
    sy = {s << 2 for s in su} | {s << 2 | t for s in sc for t in [1, 2]} | {3}
    assert len(sy) == len(Y)
    y_ordinary = [0, alpha << 2 | 2, delta << 2 | 1] + [v << 2 | t for v in ordinary for t in [1, 3, 2]]
    check_order(Y, y_ordinary, dimension + 2)
    new_alpha, new_delta = alpha << 2 | 3, delta << 2 | 3
    assert new_alpha not in Y and new_delta not in Y
    assert directions(Y, new_alpha, dimension + 2) == pa << 2 | 1
    assert directions(Y, new_delta, dimension + 2) == pd << 2 | 2
    for tip, support in [(new_alpha, pa << 2 | 1), (new_delta, pd << 2 | 2)]:
        assert support not in sy
        assert all(tip ^ sub in Y for sub in subsets(support) if sub)
    new_a = ([v << 2 | 2 for v in a_order[:-1]] + [new_alpha, 2]
             + [v << 2 | 3 for v in ordinary] + [v << 2 | 1 for v in d_order] + [0])
    new_d = ([v << 2 | 1 for v in d_order[:-1]] + [new_delta, 1]
             + [v << 2 | 3 for v in ordinary] + [v << 2 | 2 for v in a_order] + [0])
    check_order(Y | {new_alpha}, new_a, dimension + 2, last=0)
    check_order(Y | {new_delta}, new_d, dimension + 2, last=0)
    end = None
    if endpoint_order is not None:
        end = ([alpha << 2 | 2, delta << 2 | 1]
               + [v << 2 | t for v in endpoint_order[:-1] for t in [1, 3, 2]] + [3, 2, 1, 0])
        check_order(Y, end, dimension + 2, last=0)
    minors = None
    if rooted_orders is not None:
        minors = {}
        for e, kind in itertools.product(range(dimension + 2), ['root_section', 'contraction']):
            child = elementary(Y, e, kind)
            if e == 0:
                order = ([v << 2 | 2 for v in a_order] + [0] if kind == 'root_section' else
                         [alpha << 2 | 2] + [v << 2 | 2 for v in ordinary] + [v << 2 for v in d_order])
            elif e == 1:
                order = ([v << 2 | 1 for v in d_order] + [0] if kind == 'root_section' else
                         [delta << 2 | 1] + [v << 2 | 1 for v in ordinary] + [v << 2 for v in a_order])
            else:
                base = elementary(C, e - 2, kind)
                template = square(base, base, base)
                extra = sorted(child - template)
                assert template <= child and len(extra) <= 2
                order = extra + [v << 2 | t for v in rooted_orders[e - 2, kind][:-1]
                                 for t in [1, 3, 2]] + [3, 2, 1, 0]
            check_order(child, order, dimension + 2, last=0)
            minors[e, kind] = order
    return {'family': Y, 'ordinary': y_ordinary, 'alpha': new_alpha, 'delta': new_delta,
            'a_order': new_a, 'd_order': new_d, 'shadow': sy, 'dimension': dimension + 2,
            'minors': minors, 'endpoint': end}


def main():
    names = ['damp_fixed_291_repair_creation.json',
             'computations/round776_rooted_corner_descent.json',
             'computations/round793_third_endpoint_repair.json',
             'computations/round790_paired_endpoint_repairs.json']
    raw, descent, positive, repairs = [json.loads((ROOT / name).read_text()) for name in names]
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    assert positive['source_sha256'] == sha(ROOT / names[0])
    assert repairs['source_sha256'] == sha(ROOT / names[0])
    B = {v ^ 29 for v in (set(raw['base']) | {91}) - {61, 125, 93}}
    core = next(row for row in descent['tests'] if row['deleted'] == [61, 125, 93])
    assert len(B) == 289 and corners(B, 12) == [0]
    assert all(any(v ^ (1 << e) in B for v in B) for e in range(12))
    b_shadow = shadow(B, 12)
    assert len(b_shadow) == len(B)
    assert directions(B, 0, 12) == 262
    check_order(B, core['ordinary_peeling'], 12)
    check_negative(B, 12, core['nonrootability']['negative_certificate'])
    c, P = 186 ^ 29, 323
    d = c ^ P
    H = B | {c}
    cube = {c ^ u for u in subsets(P)}
    assert c == 167 and d == 484 and len(cube) == 16
    assert cube - B == {c} and 0 not in cube
    assert directions(B, c, 12) == P
    assert shadow(H, 12) == b_shadow | {P}
    assert max(map(int.bit_count, b_shadow | {P})) == 4
    assert corners(H, 12) == [0, c, d]
    assert directions(H, c, 12) == directions(H, d, 12) == P
    available = sorted(v for v in cube - {0, c} if not directions(H, v, 12) & ~P)
    assert available == [d] and d not in corners(B, 12)
    assert [directions(H, v, 12).bit_count() for v in corners(H, 12)] == [3, 4, 4]
    exterior = []
    for v in sorted(cube - {c, d}):
        E = directions(H, v, 12) & ~P
        assert E
        e_bit = E & -E
        assert v ^ e_bit in H and v ^ e_bit not in cube
        assert c ^ e_bit not in H and d ^ e_bit not in H
        exterior.append({'cube_vertex': v, 'outside_direction': e_bit.bit_length() - 1,
                         'outside_neighbor': v ^ e_bit})
    # These are all root-protected states, not a sampled negative path.
    neg = [{'removed': [], 'corners_except_root': [c, d]},
           {'removed': [c], 'corners_except_root': []},
           {'removed': [d], 'corners_except_root': []}]
    check_negative(H, 12, neg)
    assert corners(H - {c}, 12) == corners(H - {d}, 12) == [0]
    ordinary = [c] + core['ordinary_peeling']
    check_order(H, ordinary, 12)
    minor_orders, b_minor_orders = [], {}
    checked = set()
    for row in core['elementary_checks']:
        key = row['coordinate'], row['operation']
        child_b, child_h = elementary(B, *key), elementary(H, *key)
        order_b = expand_dense(row['peeling'], child_b)
        check_order(child_b, order_b, 12, last=0)
        b_minor_orders[key] = order_b
        extras = sorted(child_h - child_b)
        assert len(extras) <= 1
        order_h = extras + order_b
        check_order(child_h, order_h, 12, last=0)
        minor_orders.append({'coordinate': key[0], 'operation': key[1],
                             'new_concept_prefix': extras, 'peeling': order_h})
        checked.add(key)
    assert checked == set(itertools.product(range(12), ['root_section', 'contraction']))
    gamma, Q = positive['tip'], positive['acquired_support']
    assert gamma == 604 and Q == 1632
    K = B | {gamma}
    cube_positive = {gamma ^ u for u in subsets(Q)}
    assert shadow(K, 12) == b_shadow | {Q}
    check_order(K, positive['peeling'], 12, last=0)
    available_positive = sorted(v for v in cube_positive - {0, gamma}
                                if not directions(K, v, 12) & ~Q)
    assert available_positive == [1084, 1116, 1148]
    assert P.bit_count() == Q.bit_count() == 4
    assert c.bit_count() == gamma.bit_count() == 5
    assert (c & P).bit_count() == (gamma & Q).bit_count() == 2
    # Each first old-corner deletion removes the newly acquired support.
    # Check the full signed minimal-nonface presentations of both pivots.
    minimal_nonfaces = [mask for mask in range(1, 1 << 12) if mask not in b_shadow
                        and all(mask ^ (1 << e) in b_shadow
                                for e in range(12) if mask >> e & 1)]
    def signing(family):
        result = {}
        for mask in minimal_nonfaces:
            missing = set(subsets(mask)) - {v & mask for v in family}
            assert len(missing) == 1
            result[mask] = missing.pop()
        return result
    old_signing = signing(B)
    pivots = []
    for tip, support, deleted, succeeds in [(c, P, d, False),
                                           (gamma, Q, positive['peeling'][0], True)]:
        assert deleted == tip ^ support
        pivot = (B | {tip}) - {deleted}
        assert shadow(pivot, 12) == b_shadow
        new_signing = signing(pivot)
        changed = [{'support': mask, 'before': old_signing[mask], 'after': new_signing[mask]}
                   for mask in minimal_nonfaces if old_signing[mask] != new_signing[mask]]
        assert changed == [{'support': support, 'before': tip & support,
                            'after': deleted & support}]
        cube_here = {tip ^ u for u in subsets(support)}
        E = {v for v in cube_here - {0, tip} if not directions(B | {tip}, v, 12) & ~support}
        expected = sorted(v for v in E if (v ^ deleted).bit_count() == 1)
        assert [v for v in corners(pivot, 12) if v] == expected
        row = {'family': sorted(pivot), 'tip': tip, 'deleted_old_corner': deleted,
               'changed_clauses': changed, 'corners': corners(pivot, 12),
               'peels_to_root': succeeds}
        if succeeds:
            order = positive['peeling'][1:]
            check_order(pivot, order, 12, last=0)
            assert tip.bit_count() == deleted.bit_count() == 5
            assert {v for v in pivot if v.bit_count() <= 4} == {v for v in B if v.bit_count() <= 4}
            row['endpoint_peeling'] = order
            row['unchanged_root_ball_radius'] = 4
        else:
            check_negative(pivot, 12, [{'removed': [], 'corners_except_root': []}])
        pivots.append(row)
    # Match the marked root, added tip, and entire punctured cube. This
    # comparison deliberately does not identify the exterior attachments.
    permutation = [-1] * 12
    for in_support, tip_bit in itertools.product([0, 1], repeat=2):
        first = [e for e in range(12) if (P >> e & 1, c >> e & 1) == (in_support, tip_bit)]
        second = [e for e in range(12) if (Q >> e & 1, gamma >> e & 1) == (in_support, tip_bit)]
        assert len(first) == len(second)
        for e, f in zip(first, second):
            permutation[e] = f
    assert sorted(permutation) == list(range(12))
    permute = lambda x: sum((x >> e & 1) << f for e, f in enumerate(permutation))
    assert permute(0) == 0 and permute(c) == gamma and permute(P) == Q
    assert {permute(v) for v in cube - {c}} == cube_positive - {gamma}
    # The negative and positive families have equal shadows and differ
    # in just one signed minimal nonface, increasingly far from the root.
    K0 = (B | {gamma}) - {positive['peeling'][0]}
    k_order = positive['peeling'][1:]
    base_a, base_d = [r['peeling'] for r in repairs['repairs']]
    neg_state = {'family': B, 'ordinary': core['ordinary_peeling'], 'alpha': 128,
                 'delta': 290, 'a_order': base_a, 'd_order': base_d, 'dimension': 12,
                 'minors': b_minor_orders, 'endpoint': None, 'shadow': b_shadow}
    pos_state = {'family': K0, 'ordinary': k_order, 'alpha': 128, 'delta': 290,
                 'a_order': [128] + k_order, 'd_order': [290] + k_order, 'dimension': 12,
                 'minors': None, 'endpoint': k_order, 'shadow': b_shadow}
    tower = []
    for level in range(3):
        F, G = neg_state['family'], pos_state['family']
        sigma, dim = neg_state['shadow'], neg_state['dimension']
        assert sigma == pos_state['shadow'] and len(F) == len(G) == (581 * 3**level - 3) // 2
        assert dim == 12 + 2 * level and max(map(int.bit_count, sigma)) == 3 + level
        assert directions(F, 0, dim) == directions(G, 0, dim) == (262 if level == 0 else 3)
        difference = F ^ G
        assert min(v.bit_count() for v in difference) == 5 + level
        assert {v for v in F if v.bit_count() <= 4 + level} == {v for v in G if v.bit_count() <= 4 + level}
        I = 1632 << (2 * level)
        minimal = [mask for mask in range(1, 1 << dim) if mask not in sigma and
                   all(mask ^ (1 << e) in sigma for e in range(dim) if mask >> e & 1)]
        changed = []
        for mask in minimal:
            traces_f, traces_g = {v & mask for v in F}, {v & mask for v in G}
            if traces_f != traces_g:
                missing_f, missing_g = set(subsets(mask)) - traces_f, set(subsets(mask)) - traces_g
                assert len(missing_f) == len(missing_g) == 1
                changed.append({'support': mask, 'before': missing_f.pop(), 'after': missing_g.pop()})
        assert changed == [{'support': I, 'before': 576 << (2 * level), 'after': 1056 << (2 * level)}]
        row = {'level': level, 'order': len(F), 'dimension': dim, 'vc_dimension': 3 + level,
               'negative_family': sorted(F), 'positive_family': sorted(G),
               'equal_shadow': sorted(sigma), 'changed_clause': changed[0],
               'unchanged_root_ball_radius': 4 + level, 'first_difference_distance': 5 + level,
               'negative_ordinary_peeling': neg_state['ordinary'], 'positive_endpoint_peeling': pos_state['endpoint'],
               'repair_tips': [neg_state['alpha'], neg_state['delta']],
               'negative_repair_endpoint_orders': [neg_state['a_order'], neg_state['d_order']],
               'positive_repair_endpoint_orders': [pos_state['a_order'], pos_state['d_order']],
               'rooted_minor_peelings': [{'coordinate': e, 'operation': kind, 'peeling': order}
                                        for (e, kind), order in neg_state['minors'].items()]}
        tower.append(row)
        print('Tower level', level, 'order', len(F), 'equal root ball', 4 + level,
              'one changed clause', I, flush=True)
        if level < 2:
            def next_state(s):
                return lift(s['family'], s['ordinary'], s['alpha'], s['delta'], s['a_order'],
                            s['d_order'], s['dimension'], s['minors'], s['endpoint'])
            neg_state, pos_state = next_state(neg_state), next_state(pos_state)
    report = {'schema': 1, 'input_sha256': {n: sha(ROOT / n) for n in names},
              'verifier_sha256': sha(Path(__file__)), 'all_checks_passed': True,
              'core': sorted(B), 'family': sorted(H), 'root': 0, 'tip': c,
              'acquired_support': P, 'cube': sorted(cube), 'antipodal_old_corner': d,
              'corners': corners(H, 12), 'available_old_nonroot_corners': available,
              'exterior_neighbor_witnesses': exterior, 'negative_certificate': neg,
              'ordinary_peeling': ordinary, 'rooted_elementary_peelings': minor_orders,
              'positive_control': {'tip': gamma, 'acquired_support': Q,
                                   'available_old_nonroot_corners': available_positive,
                                   'endpoint_peeling': positive['peeling']},
              'marked_punctured_cube_coordinate_permutation': permutation,
              'signed_clause_pivots': pivots,
              'square_lift_tower': tower,
              'counts': {'order': 290, 'dimension': 12, 'vc_dimension': 4,
                         'negative_states': 3, 'exterior_neighbors': 14,
                         'core_rooted_minor_orders': 24, 'new_rooted_minor_orders': 24,
                         'tower_additional_rooted_minor_orders': 60},
              'scope': 'An exposed old nonroot corner is necessary but not sufficient for repair '
                       'when the core root is its only corner. The antipodal singleton case is '
                       'intrinsically nonrestoring. Marked punctured cubes without exterior '
                       'attachments do not determine the endpoint. The uniform tower defeats every '
                       'fixed root-ball radius even with the entire shadow retained; levels 0,1,2 '
                       'are checked here. Full classification remains open.'}
    target = ROOT / 'damp_rooted_clause_pivots_verification.json'
    target.write_text(json.dumps(report, indent=2) + '\n')
    print('PASS: nonrestoring tip', c, 'new old corner', d, 'negative states', len(neg))
    print('PASS: 14 exterior witnesses; 24 core and 24 new rooted minor orders')
    print('Positive comparison:', gamma, 'old corners', available_positive)
    print('Single-clause pivots:', [(r['changed_clauses'], r['corners'], r['peels_to_root']) for r in pivots])
    print('Marked punctured-cube permutation:', permutation)


if __name__ == '__main__':
    main()
