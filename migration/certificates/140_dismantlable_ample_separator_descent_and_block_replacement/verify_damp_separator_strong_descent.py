#!/usr/bin/env python3
"""Verify the fixed hypotheses and transports for separator strong descent.

No peeling search. Check ampleness from actual cubes, replay all proper
elementary-minor orders and complete negative state certificates, and
transport positive contraction orders to carriers. One saved 261-concept
order checks that a universal crossing coordinate need not give descent.
The 2318-concept tower checks the nontrivial, one-step conclusion. One
square attachment checks invariance and removes the full graph's separators.
The theorem for arbitrary inputs and all strong projections is symbolic.
"""
from itertools import combinations
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent
KINDS = ('root_section', 'opposite_section', 'contraction')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def directions(H, x, n):
    return sum(1 << e for e in range(n) if x ^ (1 << e) in H)


def corners(H, n):
    answer = []
    for x in sorted(H):
        support = directions(H, x, n)
        sub = support
        while sub and x ^ sub in H:
            sub = (sub - 1) & support
        if not sub:
            answer.append(x)
    return answer


def check_order(H, order, n):
    assert len(order) == len(H) and set(order) == H
    left = set(H)
    for x in order:
        support = directions(left, x, n)
        sub = support
        while sub:
            assert x ^ sub in left, (x, support, sub)
            sub = (sub - 1) & support
        left.remove(x)
    assert not left


def cube_shadow(H, n):
    """Enumerate actual cubes from their coordinatewise smallest vertex.

    Equal strongly shattered support count and vertex count is the lower
    Sandwich equality. No claim depends on an assumed VC-rank cutoff.
    """
    supports = set()
    for x in H:
        def extend(mask, vertices, start):
            supports.add(mask)
            for e in range(start, n):
                bit = 1 << e
                if x & bit:
                    continue
                opposite = [v | bit for v in vertices]
                if all(v in H for v in opposite):
                    extend(mask | bit, vertices + opposite, e + 1)
        extend(0, [x], 0)
    assert len(supports) == len(H)
    return supports


def strong(H, e):
    return {x for x in H if not x & (1 << e) and x ^ (1 << e) in H}


def image(H, e, kind):
    selected = H if kind == 'contraction' else {
        x for x in H if bool(x & (1 << e)) == (kind == 'opposite_section')}
    return {x & ~(1 << e) for x in selected}


def crossing(H, n):
    graph = {e: set() for e in range(n)}
    for e, f in combinations(range(n), 2):
        mask = (1 << e) | (1 << f)
        if len({x & mask for x in H}) == 4:
            graph[e].add(f)
            graph[f].add(e)
    return graph


def components(vertices, graph):
    todo, answer = set(vertices), []
    while todo:
        x = todo.pop()
        found, queue = {x}, [x]
        for v in queue:
            for w in graph[v] & todo:
                todo.remove(w)
                found.add(w)
                queue.append(w)
        answer.append(found)
    return answer


def check_negative(H, rows, n):
    states = {frozenset(r['removed']): r for r in rows}
    assert len(states) == len(rows) and frozenset() in states
    for removed, row in states.items():
        assert removed <= H and H - removed
        cs = corners(H - removed, n)
        assert cs == sorted(row['corners'])
        assert all(removed | {x} in states for x in cs)
    # All corner-deletion successors stay in the supplied nonempty states.
    return len(states)


def strong_order(H, order, e):
    """Order full e-fibres by their first disappearing endpoint."""
    carrier = strong(H, e)
    seen, result = set(), []
    for x in order:
        y = x & ~(1 << e)
        if y in carrier and y not in seen:
            seen.add(y)
            result.append(y)
    assert seen == carrier
    return result


def main():
    names = [
        'damp_terminal_cut_vertex_verification.json',
        'damp_terminal_edge_amalgam_verification.json',
        'damp_coupled_edge_minors_verification.json',
        'damp_cornered_maximal_seed_verification.json',
        'computations/round819_separator_carrier_order.json',
    ]
    data = [json.loads((ROOT / name).read_text()) for name in names]
    assert all(d['all_checks_passed'] for d in data[:4])
    seed = set(data[0]['family'])
    positive = data[4]
    assert positive['source_sha256'] == sha(ROOT / names[1])
    assert positive['generator_sha256'] == sha(
        ROOT / 'computations/probe_round819_separator_carrier_order.py')

    a, b, bit = 128 << 14, 290 << 14, 1 << 26
    tower = seed | {x | bit for x in seed} | {a | bit, b}
    tower_record = next(r for r in data[3]['fixed_towers'] if r['layers'] == 1)
    assert hashlib.sha256(json.dumps(sorted(tower), separators=(',', ':')).encode()
                          ).hexdigest() == tower_record['family_sha256']
    assert len(tower) == 2318 and strong(tower, 26) == seed
    seed_orders = {(r['coordinate'], r['operation']): r['peeling']
                   for r in data[0]['elementary_minor_peelings']}
    repair_orders = {int(t): order for t, order in
                     data[3]['seed']['ordinary_repair_peelings'].items()}
    assert set(repair_orders) == {a, b}
    for tip, order in repair_orders.items():
        check_order(seed | {tip}, order, 26)

    def tower_order(e, kind):
        child = image(tower, e, kind)
        if e < 26:
            base, base_order = image(seed, e, kind), seed_orders[e, kind]
            layers = [bit, 0]
        else:
            layers = [0]
            if kind == 'contraction':
                base, base_order = seed | {a, b}, [b] + repair_orders[a]
            else:
                tip = a if kind == 'opposite_section' else b
                base, base_order = seed | {tip}, repair_orders[tip]
        product = {x | t for x in base for t in layers}
        assert product <= child
        extras = child - product
        old_mask = bit - 1
        kinds = {a & ~(1 << e): 'up', b & ~(1 << e): 'down'}
        order = []
        for v in sorted({x & old_mask for x in extras}):
            order.extend(sorted((x for x in extras if x & old_mask == v),
                                reverse=kinds[v] == 'down'))
        order.extend(x | t for x in base_order for t in layers)
        return order

    families = [(names[i], set(d['family']), d['dimension'], d)
                for i, d in enumerate(data[:3])]
    families.append(('one_layer_cornered_tower', tower, 27, tower_record))
    square_tip = 32769
    enlarged = tower | {square_tip}
    assert square_tip not in tower
    assert directions(enlarged, square_tip, 27) == (1 << 0) | (1 << 15)
    assert all(square_tip ^ t in tower for t in [1, 1 << 15, (1 << 0) | (1 << 15)])
    assert len({x & square_tip for x in tower}) == 3
    assert all(directions(tower, x, 27).bit_count() >= 3 for x in tower)
    assert {x for x in enlarged if directions(enlarged, x, 27).bit_count() < 3} == {square_tip}
    assert strong(enlarged, 26) == seed
    families.append(('square_extended_tower', enlarged, 27, tower_record))

    output, parent_order_count, transported_count = [], 0, 0
    for label, H, n, record in families:
        sigma = cube_shadow(H, n)
        assert all(any(x ^ (1 << e) in H for x in H) for e in range(n))
        orders = {}
        for row in record['elementary_minor_peelings']:
            e, kind = row['coordinate'], row['operation']
            assert (e, kind) not in orders
            child = image(H, e, kind)
            assert len(child) < len(H)
            order = row.get('peeling')
            if label == 'square_extended_tower':
                base = image(tower, e, kind)
                assert base <= child and len(child - base) <= 1
                order = sorted(child - base) + tower_order(e, kind)
            elif order is None:
                assert label == 'one_layer_cornered_tower'
                order = tower_order(e, kind)
                assert hashlib.sha256(json.dumps(order, separators=(',', ':')).encode()
                                      ).hexdigest() == row['peeling_sha256']
            check_order(child, order, n)
            orders[e, kind] = order
            parent_order_count += 1
        assert set(orders) == {(e, k) for e in range(n) for k in KINDS}

        negative_states = None
        if label not in ['one_layer_cornered_tower', 'square_extended_tower']:
            rows = (record['negative_certificate'] if 'negative_certificate' in record
                    else record['negative_certificates']['X'])
            negative_states = check_negative(H, rows, n)
        else:
            # Its strong projection is the independently checked negative seed.
            assert output[0]['negative_state_count'] == 4

        gamma = crossing(H, n)
        assert {sum(1 << e for e in pair) for pair in combinations(range(n), 2)
                if pair[1] in gamma[pair[0]]} == {s for s in sigma if s.bit_count() == 2}
        parts = components(range(n), gamma)
        cuts = [e for e in range(n)
                if len(components(set(range(n)) - {e}, gamma)) > len(parts)]
        universal = [e for e in range(n) if len(gamma[e]) == n - 1]
        if label == 'square_extended_tower':
            assert len(parts) == 1 and not cuts
        else:
            assert len(parts) > 1 or cuts
        assert len(universal) <= 1
        carrier_rows = []
        for e in range(n):
            K = strong(H, e)
            lost = [f for f in range(n) if f != e and
                    not any(x ^ (1 << f) in K for x in K)]
            assert set(lost) == set(range(n)) - {e} - gamma[e]
            if lost:
                q = lost[0]
                projected = image(H, q, 'contraction')
                projected_carrier = strong(projected, e)
                assert projected_carrier == image(K, q, 'contraction')
                constants = {x & (1 << q) for x in K}
                assert len(constants) <= 1
                constant = next(iter(constants), 0)
                order = [x | constant for x in
                         strong_order(projected, orders[q, 'contraction'], e)]
                check_order(K, order, n)
                transported_count += 1
                carrier_rows.append(dict(coordinate=e, order=len(K), status='peelable',
                                         lost_coordinate=q, peeling=order))
            elif label == names[1]:
                assert e == positive['coordinate'] == 3
                assert K == set(positive['family']) and len(K) == 261
                check_order(K, positive['peeling'], n)
                carrier_rows.append(dict(coordinate=e, order=261, status='peelable',
                                         peeling=positive['peeling']))
            else:
                assert label in ['one_layer_cornered_tower', 'square_extended_tower'] and e == 26
                assert K == seed
                assert all(any(x ^ (1 << f) in K for x in K) for f in range(26))
                kgraph = {x: {x ^ (1 << f) for f in range(26)
                              if x ^ (1 << f) in K} for x in K}
                assert len(components(K - {0}, kgraph)) == 2
                if label == 'square_extended_tower':
                    assert len(components(set(range(n)) - {e}, gamma)) == 1
                    carrier_rows.append(dict(coordinate=e, order=len(K), status='nonpeelable',
                                             negative_source=names[0],
                                             unchanged_from_three_core=True))
                    continue
                assert len(components(set(range(n)) - {e}, gamma)) == 2
                hgraph = {x: {x ^ (1 << f) for f in range(n)
                              if x ^ (1 << f) in H} for x in H}
                edge = {0, 1 << e}
                hparts = components(H - edge, hgraph)
                assert len(hparts) == 2
                recovered = [strong(part | edge, e) for part in hparts]
                expected = [set(data[0]['rooted_factor']),
                            {x << 14 for x in data[0]['core']}]
                assert {frozenset(x) for x in recovered} == {frozenset(x) for x in expected}
                assert recovered[0] & recovered[1] == {0}
                assert recovered[0] | recovered[1] == K
                carrier_rows.append(dict(coordinate=e, order=len(K), status='nonpeelable',
                                         negative_source=names[0],
                                         recovered_rooted_block_orders=sorted(map(len, recovered)),
                                         separating_edge=sorted(edge)))
        negatives = [r['coordinate'] for r in carrier_rows if r['status'] == 'nonpeelable']
        assert negatives == ([26] if label in ['one_layer_cornered_tower', 'square_extended_tower'] else [])
        output.append(dict(source=label, order=len(H), dimension=n,
                           vc_dimension=max(map(int.bit_count, sigma)),
                           negative_state_count=negative_states,
                           crossing_components=sorted(sorted(p) for p in parts),
                           crossing_cut_vertices=cuts, universal_coordinates=universal,
                           carriers=carrier_rows, negative_carrier_coordinates=negatives))
        print(label, len(H), 'universal', universal, 'negative carriers', negatives, flush=True)

    report = {
        'schema': 'damp-separator-strong-descent-v1',
        'verifier_sha256': sha(Path(__file__)),
        'input_sha256': {name: sha(ROOT / name) for name in names},
        'families': output,
        'proper_minor_orders_checked': parent_order_count,
        'transported_carrier_orders_checked': transported_count,
        'direct_carrier_orders_checked': 1,
        'all_checks_passed': True,
        'scope': 'Five fixed ample forbidden inputs checked independently from cubes, '
                 'all proper elementary orders, and complete negative certificates or '
                 'the checked negative strong projection. All one-coordinate carriers '
                 'classified by explicit positive orders or the checked forbidden seed. '
                 'The square extension retains exactly its original three-core and '
                 'negative carrier while making the crossing graph two-connected. '
                 'All-size and all-subset conclusions use the uniform theorem, not a census.',
    }
    (ROOT / 'damp_separator_strong_descent_verification.json').write_text(
        json.dumps(report, indent=2) + '\n')
    print('All checks passed;', parent_order_count, 'proper minor orders;',
          transported_count, 'transported carrier orders; one direct carrier order.')


if __name__ == '__main__':
    main()
