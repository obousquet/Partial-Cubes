#!/usr/bin/env python3
"""Replay a terminal cut-vertex obstruction; no endpoint search is performed.

The 870-concept rooted factor has two nonroot corners, neither permitting
a rooted-minimal deletion. Its wedge with the 289-concept unique-corner
factor is a 1158-concept terminal ample forbidden pc-minor.
"""
import hashlib
import itertools
import json
from pathlib import Path

from verify_damp_cut_vertex_obstruction import check_order, expand_dense
from verify_damp_layer_midpoint_cores import shadow
from verify_damp_rank_three_corner_neighbors import corners, directions

ROOT = Path(__file__).resolve().parent


def elementary(family, coordinate, operation):
    bit = 1 << coordinate
    return {x & ~bit for x in family
            if operation == 'contraction'
            or bool(x & bit) == (operation == 'opposite_section')}


def components(family, dimension):
    todo, result = set(family), []
    while todo:
        pending = [todo.pop()]
        found = set(pending)
        while pending:
            x = pending.pop()
            for e in range(dimension):
                y = x ^ (1 << e)
                if y in todo:
                    todo.remove(y)
                    found.add(y)
                    pending.append(y)
        result.append(found)
    return result


def main():
    inputs = ['damp_fixed_291_repair_creation.json',
              'computations/round776_rooted_corner_descent.json',
              'computations/round790_paired_endpoint_repairs.json']
    raw, rooted, repairs = [json.loads((ROOT / p).read_text()) for p in inputs]
    hashes = {p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest() for p in inputs}
    assert repairs['source_sha256'] == hashes[inputs[0]]
    C = {x ^ 29 for x in (set(raw['base']) | {91}) - {61, 125, 93}}
    core_row = next(r for r in rooted['tests'] if r['deleted'] == [61, 125, 93])
    c_order = core_row['ordinary_peeling']
    assert len(C) == 289 and corners(C, 12) == [0]
    check_order(C, c_order, 12)
    sh_c = shadow(C, 12)
    assert len(sh_c) == len(C) and max(map(int.bit_count, sh_c)) == 3
    core_orders = {}
    for row in core_row['elementary_checks']:
        key = row['coordinate'], row['operation']
        child = elementary(C, *key)
        order = expand_dense(row['peeling'], child)
        check_order(child, order, 12, last=0)
        core_orders[key] = order
    assert set(core_orders) == set(itertools.product(range(12), ['root_section', 'contraction']))

    alpha, delta = 128, 290
    A, D = C | {alpha}, C | {delta}
    a_order, d_order = [r['peeling'] for r in repairs['repairs']]
    assert [r['tip'] for r in repairs['repairs']] == [alpha, delta]
    check_order(A, a_order, 12, last=0)
    check_order(D, d_order, 12, last=0)
    sh_a, sh_d, sh_union = [shadow(F, 12) for F in [A, D, A | D]]
    assert sh_a == sh_c | {390} and sh_d == sh_c | {480}
    assert sh_union == sh_c | {390, 480} and len(sh_union) == len(A | D)
    Y = {0} | {x << 2 | 2 for x in A} | {x << 2 | 1 for x in D} | {x << 2 | 3 for x in C}
    tips = [alpha << 2 | 2, delta << 2 | 1]
    sh_y = shadow(Y, 14)
    expected_shadow = {s << 2 for s in sh_union} | {s << 2 | t for s in sh_c for t in [1, 2]} | {3}
    assert sh_y == expected_shadow and len(sh_y) == len(Y) == 870
    assert directions(Y, 0, 14) == 3
    y_order = [0] + tips + [x << 2 | t for x in c_order for t in [1, 3, 2]]
    check_order(Y, y_order, 14)
    negative_y, negative_x = [], []
    X = Y | {x << 14 for x in C}
    for size in range(3):
        for removed in itertools.combinations(tips, size):
            remaining_tips = sorted(set(tips) - set(removed))
            assert corners(Y - set(removed), 14) == [0] + remaining_tips
            assert corners(X - set(removed), 26) == remaining_tips
            assert len(Y - set(removed)) > 1
            negative_y.append({'removed': list(removed), 'corners_except_root': remaining_tips})
            negative_x.append({'removed': list(removed), 'corners': remaining_tips})
    assert len(negative_y) == len(negative_x) == 4

    y_rows, y_orders = [], {}
    for e, kind in itertools.product(range(14), ['root_section', 'contraction']):
        child = elementary(Y, e, kind)
        if e == 0:
            order = ([x << 2 | 2 for x in a_order] + [0] if kind == 'root_section'
                     else [tips[0]] + [x << 2 | 2 for x in c_order] + [x << 2 for x in d_order])
        elif e == 1:
            order = ([x << 2 | 1 for x in d_order] + [0] if kind == 'root_section'
                     else [tips[1]] + [x << 2 | 1 for x in c_order] + [x << 2 for x in a_order])
        else:
            core = elementary(C, e - 2, kind)
            core_order = core_orders[e - 2, kind]
            product = {0} | {x << 2 | t for x in core for t in [1, 2, 3]}
            assert product <= child and len(child - product) <= 2
            order = sorted(child - product)
            order += [x << 2 | t for x in core_order[:-1] for t in [1, 3, 2]]
            order += [3, 2, 1, 0]
        check_order(child, order, 14, last=0)
        assert len(child) < len(Y)
        y_orders[e, kind] = order
        y_rows.append({'coordinate': e, 'operation': kind, 'order': len(child), 'peeling': order})

    # No mixed support can be shattered: each cross-block pair omits 11.
    for e, f in itertools.product(range(14), range(14, 26)):
        mask = (1 << e) | (1 << f)
        assert {x & mask for x in X} == {0, 1 << e, 1 << f}
    sh_x = sh_y | {s << 14 for s in sh_c}
    assert len(X) == len(sh_x) == 1158 and max(map(int.bit_count, sh_x)) == 4
    split = components(X - {0}, 26)
    assert {frozenset(F) for F in split} == {frozenset(Y - {0}), frozenset(x << 14 for x in C - {0})}
    assert [directions(X, x, 26).bit_count() for x in tips] == [4, 4]
    assert directions(X, 0, 26).bit_count() == 5
    min_degree = min(directions(X, x, 26).bit_count() for x in X)
    assert min_degree == 4
    x_rows = []
    for e, kind in itertools.product(range(26), ['root_section', 'opposite_section', 'contraction']):
        child = elementary(X, e, kind)
        if e < 14:
            if kind == 'opposite_section':
                order = [x & ~(1 << e) for x in y_order if x >> e & 1]
            else:
                order = y_orders[e, kind][:-1] + [x << 14 for x in c_order]
        else:
            old_e = e - 14
            if kind == 'opposite_section':
                order = [(x & ~(1 << old_e)) << 14 for x in c_order if x >> old_e & 1]
            else:
                order = [x << 14 for x in core_orders[old_e, kind][:-1]] + y_order
        check_order(child, order, 26)
        assert len(child) < len(X)
        x_rows.append({'coordinate': e, 'operation': kind, 'order': len(child), 'peeling': order})
    assert len(x_rows) == 78
    double = {x << 2 for x in C} | {x << 14 for x in C}
    assert len(double) == 577 and corners(double, 26) == []
    terminal = []
    for tip, restrict, contract in [(tips[0], 0, 1), (tips[1], 1, 0)]:
        intermediate = elementary(X - {tip}, restrict, 'root_section')
        witness = elementary(intermediate, contract, 'contraction')
        assert witness == double and len(witness) < len(X) - 1
        rooted_witness = elementary(elementary(Y - {tip}, restrict, 'root_section'), contract, 'contraction')
        assert rooted_witness == {x << 2 for x in C}
        terminal.append({'deleted_corner': tip, 'root_section': restrict, 'contract': contract,
                         'bad_proper_minor_order': len(witness)})
    report = {'schema': 'damp-terminal-cut-vertex-verification-v1', 'input_sha256': hashes,
              'verifier_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'core': sorted(C), 'core_ordinary_peeling': c_order,
              'rooted_factor': sorted(Y), 'rooted_ordinary_peeling': y_order,
              'rooted_elementary_peelings': y_rows, 'rooted_negative_certificate': negative_y,
              'family': sorted(X), 'order': len(X), 'dimension': 26, 'vc_dimension': 4,
              'corners': tips, 'corner_degrees': [4, 4], 'minimum_degree': min_degree,
              'cut_vertex': 0, 'components_after_cut_deletion': sorted(map(len, split)),
              'negative_certificate': negative_x, 'elementary_minor_peelings': x_rows,
              'terminal_minor_witnesses': terminal, 'cornerless_bad_minor': sorted(double),
              'checks': {'core_rooted_minors': 24, 'new_rooted_minors': 28,
                         'ordinary_elementary_minors': 78, 'cross_pairs': 168,
                         'negative_states_each': 4}, 'all_checks_passed': True,
              'scope': 'Refutes universal preserving nonroot-corner deletion and universal '
                       'nontrivial common-seed completeness for cornered terminal inputs. '
                       'Intrinsic core/factor admissibility and cornerlessness of the seed '
                       'conditional on a successful maximal layer set remain open.'}
    (ROOT / 'damp_terminal_cut_vertex_verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print({k: report[k] for k in ['order', 'dimension', 'vc_dimension', 'corners', 'minimum_degree',
                                 'components_after_cut_deletion', 'checks', 'all_checks_passed']})


if __name__ == '__main__':
    main()
