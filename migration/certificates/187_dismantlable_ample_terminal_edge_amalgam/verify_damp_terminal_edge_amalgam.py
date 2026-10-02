#!/usr/bin/env python3
"""Replay a terminal two-connected obstruction with no layer presentation.

Only stored positive orders are used. All 75 elementary minors receive
explicit constructed peelings. Complete corner-state and signed-clause
checks certify nonpeelability, terminality, and failure of layer recovery.
"""
import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path

from verify_damp_cornered_maximal_seed import cut_vertices, exact_shadow
from verify_damp_cut_vertex_obstruction import check_order
from verify_damp_layer_midpoint_cores import corner, shadow, subsets
from verify_damp_rank_three_corner_neighbors import corners, directions
from verify_damp_terminal_cut_vertex import elementary

ROOT = Path(__file__).resolve().parent
KINDS = ['root_section', 'opposite_section', 'contraction']


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def relative_check(family, order, target, dimension):
    assert target <= family and set(order) == family - target
    assert len(order) == len(family) - len(target)
    remaining = set(family)
    for x in order:
        assert corner(remaining, x, dimension)
        remaining.remove(x)
    assert remaining == target


def projected_order(order, coordinate, operation):
    bit = 1 << coordinate
    selected = [x for x in order if operation == 'contraction'
                or bool(x & bit) == (operation == 'opposite_section')]
    counts = Counter(x & ~bit for x in selected)
    result = []
    for x in selected:
        y = x & ~bit
        counts[y] -= 1
        if not counts[y]:
            result.append(y)
    return result


def main():
    names = ['damp_terminal_cut_vertex_verification.json',
             'computations/round806_edge_amalgam_orders.json']
    source, saved = [json.loads((ROOT / name).read_text()) for name in names]
    assert saved['source_sha256'] == sha(ROOT / names[0])
    assert saved['probe_sha256'] == sha(ROOT / 'computations/probe_round806_edge_amalgam.py')
    record, = saved['records']
    assert record['edge'] == 1
    C, D = set(source['core']), {0, 2}
    c_order = source['core_ordinary_peeling']
    check_order(C, c_order, 12)
    assert len(C) == 289 and corners(C, 12) == [0]
    assert not cut_vertices(C, 12)
    sc = shadow(C, 12)
    assert len(sc) == len(C)
    for support, trace in [(38, 32), (70, 64), (15, 8)]:
        assert support not in sc and support & 2
        assert all(support ^ (1 << f) in sc for f in range(12) if support >> f & 1)
        assert set(subsets(support)) - {v & support for v in C} == {trace}
    relative = {}
    for row in record['tests']:
        f, k = row['coordinate'], row['operation']
        relative_check(elementary(C, f, k), row['order'], elementary(D, f, k), 12)
        relative[f, k] = row['order']
    expected = {(f, k) for f in range(12) if f != 1
                for k in ['root_section', 'contraction']}
    expected |= {(1, k) for k in KINDS}
    assert set(relative) == expected and len(relative) == 25
    a, b = 128, 290
    repair = {r['tip']: r['order'] for r in record['repairs']}
    assert set(repair) == {a, b}
    for tip in [a, b]:
        relative_check(C | {tip}, repair[tip], D, 12)
    pa, pb = directions(C, a, 12), directions(C, b, 12)
    assert (pa, pb) == (390, 480)
    assert shadow(C | {a}, 12) == sc | {pa}
    assert shadow(C | {b}, 12) == sc | {pb}
    assert shadow(C | {a, b}, 12) == sc | {pa, pb}

    Y = {v << 2 for v in D}
    Y |= {v << 2 | 2 for v in C | {a}}
    Y |= {v << 2 | 1 for v in C | {b}}
    Y |= {v << 2 | 3 for v in C}
    tips = [a << 2 | 2, b << 2 | 1]
    y_order = [v << 2 for v in sorted(D)] + tips
    y_order += [v << 2 | t for v in c_order for t in [1, 3, 2]]
    check_order(Y, y_order, 14)
    sy = {s << 2 for s in sc | {pa, pb}}
    sy |= {s << 2 | t for s in sc for t in [1, 2]}
    sy |= {s << 2 | 3 for s in [0, 2]}
    exact_shadow(Y, sy, 14)
    assert len(Y) == 871

    mapping = {1: 3}
    mapping.update({f: 14 + i for i, f in enumerate(f for f in range(12) if f != 1)})
    inverse = {v: k for k, v in mapping.items()}
    phi = lambda v: sum(1 << mapping[f] for f in range(12) if v >> f & 1)
    B = {phi(v) for v in C}
    X = Y | B
    interface = {v << 2 for v in D}
    assert Y & B == interface == {0, 8}
    sx = sy | {phi(s) for s in sc}
    clause_count = exact_shadow(X, sx, 25)
    assert len(X) == 1158 and max(map(int.bit_count, sx)) == 4
    assert not cut_vertices(X, 25)
    assert corners(X, 25) == tips
    degrees = [directions(X, v, 25).bit_count() for v in X]
    assert min(degrees) >= 3

    negative = []
    for size in range(3):
        for removed in itertools.combinations(tips, size):
            cs = corners(X - set(removed), 25)
            assert cs == sorted(set(tips) - set(removed))
            negative.append({'removed': list(removed), 'corners': cs})

    def y_relative(f, k):
        child = elementary(Y, f, k)
        if f < 2:
            t, tip = (2, a) if f == 0 else (1, b)
            if k == 'root_section':
                order = [v << 2 | t for v in repair[tip]]
                order += [v << 2 | t for v in sorted(D)]
            else:
                other = b if f == 0 else a
                order = [tip << 2 | t] + [v << 2 | t for v in c_order]
                order += [v << 2 for v in repair[other]]
        else:
            base = elementary(C, f - 2, k)
            product = interface | {v << 2 | t for v in base for t in [1, 2, 3]}
            assert product <= child and len(child - product) <= 2
            order = sorted(child - product)
            order += [v << 2 | t for v in relative[f - 2, k] for t in [1, 3, 2]]
            order += [v << 2 | t for t in [3, 2, 1] for v in sorted(D)]
        relative_check(child, order, interface, 25)
        return order

    minor_rows = []
    for f, k in itertools.product(range(25), KINDS):
        child = elementary(X, f, k)
        if f == 3:
            order = [phi(v) for v in relative[1, k]]
            order += projected_order(y_order, f, k)
        elif f < 14:
            if k == 'opposite_section':
                order = projected_order(y_order, f, k)
            else:
                order = y_relative(f, k) + [phi(v) for v in c_order]
        else:
            g = inverse[f]
            if k == 'opposite_section':
                order = [phi(v) for v in projected_order(c_order, g, k)]
            else:
                order = [phi(v) for v in relative[g, k]] + y_order
        assert len(child) < len(X)
        check_order(child, order, 25)
        minor_rows.append({'coordinate': f, 'operation': k,
                           'order': len(child), 'peeling': order})

    double = {v << 2 for v in C} | B
    assert len(double) == 576 and not corners(double, 25)
    exact_shadow(double, {s << 2 for s in sc} | {phi(s) for s in sc}, 25)
    terminal = []
    for tip, f, g in [(tips[0], 0, 1), (tips[1], 1, 0)]:
        witness = elementary(elementary(X - {tip}, f, 'root_section'), g, 'contraction')
        assert witness == double and len(witness) < len(X) - 1
        terminal.append({'deleted_corner': tip, 'root_section': f,
                         'contraction': g, 'proper_cornerless_minor_order': len(double)})

    boundary = {s | (1 << f) for s in sx for f in range(25)} - sx
    minimal = sorted(s for s in boundary if all(s ^ (1 << f) in sx
                     for f in range(25) if s >> f & 1))
    clauses = []
    for s in minimal:
        absent = set(subsets(s)) - {v & s for v in X}
        assert len(absent) == 1
        clauses.append({'support': s, 'missing_trace': absent.pop()})
    assert len(clauses) == clause_count
    clause_degrees = [sum(bool(s & (1 << f)) for s in minimal) for f in range(25)]
    assert min(clause_degrees) >= 3
    thin = []
    for f in range(25):
        lower, upper = [elementary(X, f, k) for k in KINDS[:2]]
        if len(lower - upper) == len(upper - lower) == 1:
            thin.append(f)
    assert not thin

    imports = ['verify_damp_cornered_maximal_seed.py', 'verify_damp_cut_vertex_obstruction.py',
               'verify_damp_layer_midpoint_cores.py', 'verify_damp_rank_three_corner_neighbors.py',
               'verify_damp_terminal_cut_vertex.py']
    report = {'schema': 'damp-terminal-edge-amalgam-verification-v1',
              'input_sha256': {name: sha(ROOT / name) for name in names},
              'import_sha256': {name: sha(ROOT / name) for name in imports},
              'verifier_sha256': sha(Path(__file__)),
              'family': sorted(X), 'order': len(X), 'dimension': 25, 'vc_dimension': 4,
              'minimum_degree': min(degrees), 'corners': tips, 'interface': sorted(interface),
              'two_connected': True, 'negative_certificate': negative,
              'elementary_minor_peelings': minor_rows, 'terminal_minor_witnesses': terminal,
              'cornerless_bad_minor': sorted(double), 'signed_minimal_nonfaces': clauses,
              'clause_degrees': clause_degrees, 'balanced_coordinates': [], 'thin_coordinates': [],
              'checks': {'input_relative_orders': 25, 'input_repair_orders': 2,
                         'elementary_minor_orders': len(minor_rows), 'negative_states': 4,
                         'terminal_minor_witnesses': 2}, 'all_checks_passed': True,
              'scope': 'Refutes successful common-seed layer recovery for all two-connected cornered terminal ample obstructions. Primitive admissibility and the full intrinsic classification remain open.'}
    (ROOT / 'damp_terminal_edge_amalgam_verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: report[k] for k in ['order', 'dimension', 'minimum_degree', 'corners',
                     'clause_degrees', 'checks', 'all_checks_passed']}, indent=2))


if __name__ == '__main__':
    main()
