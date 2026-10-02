"""Replay a fixed Q4 pair whose feasible layers require two partitions.

No search for forbidden supports: every support here is Q5. The supplied
orders were found by a bounded, seeded random-peeling diagnostic. This file
checks the certificate, including all 128 assignments after forced equalities.
"""
from pathlib import Path
from itertools import product
import json
from joint_switch_primitives import orientation, descendants, topsort, constraints

PA = [9, 1, 11, 10, 14, 15, 13, 3, 12, 7, 2, 6, 8, 5, 0, 4]
PC = [11, 3, 15, 2, 13, 10, 7, 5, 14, 9, 8, 12, 6, 1, 0, 4]
BLOCKS = [{0}, {1, 3, 5, 7}, {4}, {6}, {2, 10},
          {9, 11, 13, 14, 15}, {8, 12}]


def check_order(vertices, edges, order):
    assert len(order) == len(vertices) and set(order) == set(vertices)
    rank = {x: i for i, x in enumerate(order)}
    assert all(rank[x] < rank[y] for x, y in edges)


def check_peeling(order, n):
    """Independent local cube check at every deletion, including the last."""
    remaining = set(order)
    assert len(order) == len(remaining)
    for x in order:
        directions = [1 << e for e in range(n) if x ^ (1 << e) in remaining]
        cube = {x}
        for bit in directions:
            cube |= {y ^ bit for y in tuple(cube)}
        assert cube <= remaining
        remaining.remove(x)


def main():
    A = set(range(16))
    check_peeling(PA, 4)
    check_peeling(PC, 4)
    da, dc = orientation(PA, 4), orientation(PC, 4)
    reach = descendants(A, dc)
    safe, individual = set(), []
    for x, y in da:
        edges = set(da) | {(x, z) for z in reach[y]}
        order = topsort(A, edges)
        if order is not None:
            check_order(A, edges, order)
            safe.add((x, y))
            individual.append(dict(edge=[x, y], order=order))
    groups = [{x} for x in sorted(A)]
    for x, y in set(da) - safe:
        bx = next(B for B in groups if x in B)
        by = next(B for B in groups if y in B)
        if bx != by:
            bx.update(by)
            groups.remove(by)
    assert {frozenset(B) for B in groups} == {frozenset(B) for B in BLOCKS}

    # Each branch allows ALL bit assignments on its merged equality blocks.
    branch_orders = []
    for j in (4, 5):
        groups = [B for i, B in enumerate(BLOCKS) if i not in (1, j)]
        groups.append(BLOCKS[1] | BLOCKS[j])
        label = {x: i for i, B in enumerate(groups) for x in B}
        edges = set(da) | {(x, z) for x, y in da
                          if label[x] != label[y] for z in reach[y]}
        order = topsort(A, edges)
        assert order is not None
        check_order(A, edges, order)
        branch_orders.append(dict(merge=[1, j], order=order))

    good = 0
    for bits in product((0, 1), repeat=len(BLOCKS)):
        s = {x: b for B, b in zip(BLOCKS, bits) for x in B}
        actual = constraints(A, A, da, dc, s)[0] is not None
        predicted = bits[4] == bits[1] or bits[5] == bits[1]
        assert actual == predicted
        good += actual
    assert good == 96
    S, T = BLOCKS[5], BLOCKS[4]
    for support in (S, T):
        s = {x: int(x in support) for x in A}
        assert constraints(A, A, da, dc, s)[0] is not None
    s = {x: int(x in S | T) for x in A}
    _, edges, _, invalid = constraints(A, A, da, dc, s)
    crossing = [(x, y) for x, y in da if s[x] != s[y]]
    assert not invalid and set(crossing) <= safe
    cycle = [1, 3, 14, 15, 1]
    assert all(e in edges for e in zip(cycle, cycle[1:]))
    assert {(2, 10), (10, 14), (7, 5), (5, 1)} <= set(dc)
    report = dict(status='PASS', scope='Fixed acyclic input maps, not forbidden supports',
                  contraction_order=PA, reduction_order=PC,
                  equality_blocks=[sorted(B) for B in BLOCKS],
                  individually_safe_edges=individual, branch_orders=branch_orders,
                  valid_assignments=good, assignments_after_equalities=128,
                  good_one_sets=[sorted(S), sorted(T)], bad_one_set=sorted(S | T),
                  bad_crossing_edges=crossing, augmented_cycle=cycle,
                  discovery=dict(seed=20261001, zero_based_pair_index=1367,
                                 sampling='Independent uniform choices among current corners of Q4'))
    Path(__file__).with_suffix('.json').write_text(json.dumps(report, indent=2) + '\n')
    print('PASS: 96/128 assignments; two equality partitions required; joint four-cycle.')


if __name__ == '__main__':
    main()
