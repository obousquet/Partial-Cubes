#!/usr/bin/env python3
"""Fixed two-connected terminal towers with a cornered largest recovered seed.

Use the saved 1158-concept terminal cut-vertex obstruction, transfer its
two endpoint repairs to ordinary repairs, and verify layer dimensions 1,2.
All orders are assembled from stored orders; no peeling search is used.
"""
import hashlib
import itertools
import json
from pathlib import Path

from verify_damp_cut_vertex_obstruction import check_order
from verify_damp_layer_midpoint_cores import shadow, subsets
from verify_damp_rank_three_corner_neighbors import corners, directions
from verify_damp_terminal_cut_vertex import elementary

ROOT = Path(__file__).resolve().parent
OLD = 26
OLD_MASK = (1 << OLD) - 1


def digest(value):
    return hashlib.sha256(json.dumps(value, separators=(',', ':')).encode()).hexdigest()


def exact_shadow(family, expected, dimension):
    """Check all asserted faces and exclude every minimal nonface.

    This certifies the entire shadow without enumerating Q_dimension.
    """
    assert 0 in expected
    assert all(s ^ (1 << e) in expected for s in expected
               for e in range(dimension) if s >> e & 1)
    for s in expected:
        assert len({v & s for v in family}) == 1 << s.bit_count()
    boundary = {s | (1 << e) for s in expected for e in range(dimension)} - expected
    minimal = {s for s in boundary if all(s ^ (1 << e) in expected
               for e in range(dimension) if s >> e & 1)}
    for s in minimal:
        assert len({v & s for v in family}) < 1 << s.bit_count()
    assert len(expected) == len(family)
    return len(minimal)


def cut_vertices(family, dimension):
    adjacency = {v: [v ^ (1 << e) for e in range(dimension)
                     if v ^ (1 << e) in family] for v in family}
    root = min(family)
    depth, low, parent, children = {root: 0}, {root: 0}, {root: None}, {root: 0}
    stack, cuts = [(root, iter(adjacency[root]))], set()
    while stack:
        v, edges = stack[-1]
        w = next(edges, None)
        if w is None:
            stack.pop()
            p = parent[v]
            if p is not None:
                low[p] = min(low[p], low[v])
                if parent[p] is not None and low[v] >= depth[p]:
                    cuts.add(p)
            elif children[v] > 1:
                cuts.add(v)
        elif w not in depth:
            parent[w], depth[w], low[w], children[w] = v, len(depth), len(depth), 0
            children[v] += 1
            stack.append((w, iter(adjacency[w])))
        elif w != parent[v]:
            low[v] = min(low[v], depth[w])
    assert len(depth) == len(family)
    return cuts


def tower(seed, a, b, r):
    full = (1 << r) - 1
    return ({v | (t << OLD) for v in seed for t in range(full + 1)}
            | {a | (t << OLD) for t in range(1, full + 1)}
            | {b | (t << OLD) for t in range(full)})


def strong(family, layer_bits):
    return {v & ~layer_bits for v in family
            if all((v & ~layer_bits) | t in family for t in subsets(layer_bits))}


def assembled_order(family, base, base_order, layer_bits, tip_kinds):
    layer_words = sorted(subsets(layer_bits), reverse=True)
    product = {v | t for v in base for t in layer_words}
    assert product <= family
    extras = family - product
    order = []
    for v in sorted({x & OLD_MASK for x in extras}):
        row = sorted((x for x in extras if x & OLD_MASK == v),
                     reverse=tip_kinds[v] == 'down')
        order.extend(row)
    order.extend(v | t for v in base_order for t in layer_words)
    return order


def main():
    names = ['damp_terminal_cut_vertex_verification.json',
             'computations/round790_paired_endpoint_repairs.json']
    source, repair_source = [json.loads((ROOT / name).read_text()) for name in names]
    assert source['all_checks_passed']
    for name, h in source['input_sha256'].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == h
    assert hashlib.sha256((ROOT / 'verify_damp_terminal_cut_vertex.py').read_bytes()).hexdigest() == source['verifier_sha256']
    C, Y, S = (set(source[k]) for k in ['core', 'rooted_factor', 'family'])
    assert len(S) == 1158 and cut_vertices(S, OLD) == {0}
    c_shadow, y_shadow = shadow(C, 12), shadow(Y, 14)
    sigma = y_shadow | {s << 14 for s in c_shadow}
    exact_shadow(S, sigma, OLD)
    seed_corners = corners(S, OLD)
    assert seed_corners == sorted(source['corners']) and len(seed_corners) == 2
    for removed in [set(t) for k in range(3) for t in itertools.combinations(seed_corners, k)]:
        assert corners(S - removed, OLD) == sorted(set(seed_corners) - removed)
    source_orders = {}
    for row in source['elementary_minor_peelings']:
        e, operation = row['coordinate'], row['operation']
        child = elementary(S, e, operation)
        check_order(child, row['peeling'], OLD)
        source_orders[e, operation] = row['peeling']
    assert len(source_orders) == 78

    a, b = 128 << 14, 290 << 14
    supports = {a: 390 << 14, b: 480 << 14}
    repaired_orders = {}
    for tip, row in zip([a, b], repair_source['repairs']):
        assert row['tip'] << 14 == tip
        check_order(C | {row['tip']}, row['peeling'], 12, last=0)
        order = [v << 14 for v in row['peeling'][:-1]] + source['rooted_ordinary_peeling']
        check_order(S | {tip}, order, OLD)
        assert directions(S, tip, OLD) == supports[tip]
        exact_shadow(S | {tip}, sigma | {supports[tip]}, OLD)
        repaired_orders[tip] = order
    union_order = [b] + repaired_orders[a]
    check_order(S | {a, b}, union_order, OLD)
    exact_shadow(S | {a, b}, sigma | set(supports.values()), OLD)

    double = {v << 2 for v in C} | {v << 14 for v in C}
    assert len(double) == 577 and not corners(double, OLD)
    fixed = []
    for r in [1, 2]:
        dimension, full = OLD + r, (1 << r) - 1
        layer_bits = full << OLD
        X = tower(S, a, b, r)
        expected = {s | (t << OLD) for s in sigma for t in range(full + 1)}
        expected |= {p | (t << OLD) for p in supports.values() for t in range(full)}
        minimal_nonfaces = exact_shadow(X, expected, dimension)
        assert len(X) == 1160 * (1 << r) - 2
        assert max(map(int.bit_count, expected)) == r + 4
        assert strong(X, layer_bits) == S
        assert cut_vertices(X, dimension) == set()
        cs = corners(X, dimension)
        expected_corners = {v | (t << OLD) for v in seed_corners for t in range(full + 1)}
        expected_corners |= {a | ((1 << j) << OLD) for j in range(r)}
        expected_corners |= {b | ((full ^ (1 << j)) << OLD) for j in range(r)}
        assert set(cs) == expected_corners and len(cs) == 2 ** (r + 1) + 2 * r
        minimum_degree = min(directions(X, v, dimension).bit_count() for v in X)
        assert minimum_degree == r + 3
        thin = []
        for e in range(dimension):
            left, right = [elementary(X, e, k) for k in ['root_section', 'opposite_section']]
            if len(left - right) == len(right - left) == 1:
                thin.append(e)
        assert thin == list(range(OLD, dimension))

        rows = []
        for e, operation in itertools.product(range(dimension),
                ['root_section', 'opposite_section', 'contraction']):
            child = elementary(X, e, operation)
            if e < OLD:
                base = elementary(S, e, operation)
                base_order = source_orders[e, operation]
                kinds = {a & ~(1 << e): 'up', b & ~(1 << e): 'down'}
                remaining_layers = layer_bits
            else:
                remaining_layers = layer_bits & ~(1 << e)
                kinds = {a: 'up', b: 'down'}
                if operation == 'contraction':
                    base, base_order = S | {a, b}, union_order
                else:
                    full_tip = a if operation == 'opposite_section' else b
                    base, base_order = S | {full_tip}, repaired_orders[full_tip]
            order = assembled_order(child, base, base_order, remaining_layers, kinds)
            check_order(child, order, dimension)
            rows.append({'coordinate': e, 'operation': operation,
                         'order': len(child), 'peeling_sha256': digest(order)})

        extension_rows = []
        for tip, missing_word in [(a, 0), (b, full)]:
            new = tip | (missing_word << OLD)
            upper = X | {new}
            order = assembled_order(upper, S | {tip}, repaired_orders[tip],
                                    layer_bits, {a: 'up', b: 'down'})
            check_order(upper, order, dimension)
            extension_rows.append({'tip': new, 'peeling_sha256': digest(order)})

        terminal_rows = []
        for deleted in cs:
            old, t = deleted & OLD_MASK, deleted >> OLD
            residual = X - {deleted}
            if old in [a, b]:
                difference = t if old == a else t ^ full
                assert difference.bit_count() == 1
                j = difference.bit_length() - 1
                operation = 'opposite_section' if t >> j & 1 else 'root_section'
                witness = elementary(residual, OLD + j, operation)
                remaining_layers = layer_bits & ~(1 << (OLD + j))
                assert strong(witness, remaining_layers) == S
                target = 'terminal_1158'
                steps = [[OLD + j, operation]]
            else:
                entry = next(q for q in source['terminal_minor_witnesses'] if q['deleted_corner'] == old)
                steps = [[entry['root_section'], 'root_section'], [entry['contract'], 'contraction']]
                witness = residual
                for e, operation in steps:
                    witness = elementary(witness, e, operation)
                assert strong(witness, layer_bits) == double
                # The missing seed-corner occurrence becomes a third proper row.
                assert len(witness) == 580 * (1 << r) - 3
                target = 'cornerless_577'
            assert len(witness) < len(residual)
            terminal_rows.append({'deleted_corner': deleted, 'steps': steps,
                                  'proper_minor_order': len(witness),
                                  'nonpeelable_strong_projection': target})
        row = {'layers': r, 'order': len(X), 'dimension': dimension,
               'vc_dimension': r + 4, 'minimum_degree': minimum_degree,
               'corners': cs, 'thin_coordinates': thin,
               'largest_recovered_seed_order': len(S), 'largest_recovered_seed_corners': seed_corners,
               'cut_vertices': [], 'minimal_nonface_count': minimal_nonfaces,
               'elementary_minor_peelings': rows, 'resolving_extensions': extension_rows,
               'terminal_minor_witnesses': terminal_rows, 'family_sha256': digest(sorted(X))}
        fixed.append(row)
        print({k: row[k] for k in ['layers', 'order', 'dimension', 'minimum_degree',
                                   'thin_coordinates', 'largest_recovered_seed_order']}, flush=True)

    report = {'schema': 'damp-cornered-maximal-seed-v1', 'all_checks_passed': True,
              'input_sha256': {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in names},
              'verifier_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'helper_sha256': {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in
                               ['verify_damp_cut_vertex_obstruction.py', 'verify_damp_layer_midpoint_cores.py',
                                'verify_damp_rank_three_corner_neighbors.py', 'verify_damp_terminal_cut_vertex.py']},
              'seed': {'order': len(S), 'corners': seed_corners, 'cut_vertex': 0,
                       'resolving_tips': [a, b], 'acquired_supports': list(supports.values()),
                       'ordinary_repair_peelings': {str(t): o for t, o in repaired_orders.items()}},
              'fixed_towers': fixed,
              'counts': {'source_elementary_orders': 78, 'new_elementary_orders': sum(len(x['elementary_minor_peelings']) for x in fixed),
                         'terminal_minor_witnesses': sum(len(x['terminal_minor_witnesses']) for x in fixed),
                         'tower_resolving_orders': 4},
              'scope': 'Refutes cornerlessness of the largest recovered seed for terminal inputs, even for two-connected graph three-cores. It does not refute common-seed recovery with cornered seeds or complete intrinsic core admissibility.'}
    (ROOT / 'damp_cornered_maximal_seed_verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(report['counts'], flush=True)
    print('All cornered largest-seed checks passed.', flush=True)


if __name__ == '__main__':
    main()
