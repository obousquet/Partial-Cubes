"""Verify component-label corner rules and the 267-seed equality obstruction.

Exact partition recursion is run on small positive/negative R_per controls;
the large seed uses a closed first-corner obstruction plus explicit positive
peeling orders for its contraction and overlap, never a large-state search.
"""
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'papers/cover-recursive/research'))
from component_label_recursion import components, lift, shadow
from explore_rper_cyclic_families import contraction, section


def corner(H, x, n):
    directions = [1 << e for e in range(n) if x ^ (1 << e) in H]
    masks = [0]
    for d in directions:
        masks += [s | d for s in masks]
    return all(x ^ s in H for s in masks)


def replay(H, order, n):
    left = set(H)
    assert len(order) == len(left) - 1
    for x in order:
        assert x in left and corner(left, x, n)
        left.remove(x)
    assert len(left) == 1


def greedy_order(H, n):
    left, out = set(H), []
    while len(left) > 1:
        candidates = [x for x in sorted(left) if corner(left, x, n)]
        if not candidates:
            return None  # Unknown: failure is not a nonpeelability certificate.
        x = candidates[0]
        out.append(x)
        left.remove(x)
    replay(H, out, n)
    return out


def equalities_partition(k, groups):
    parent = list(range(k))
    def root(i):
        while parent[i] != i:
            i = parent[i]
        return i
    for group in groups:
        group = list(group)
        for i in group[1:]:
            parent[root(i)] = root(group[0])
    blocks = {}
    for i in range(k):
        blocks.setdefault(root(i), []).append(i)
    return tuple(sorted(tuple(b) for b in blocks.values()))


def refines(p, q):
    """The equalities required by p are a subset of those required by q."""
    return all(any(set(a) <= set(b) for b in q) for a in p)


def partition_recursion(n):
    @lru_cache(None)
    def solve(A, C):
        blocks = components(A, C, n)
        k = len(blocks)
        discrete = tuple((i,) for i in range(k))
        if len(A) + len(C) <= 1:
            return (discrete,)
        parent_ids = {x: i for i, block in enumerate(blocks) for x in block}
        moves = [(A - {x}, C) for x in A - C if corner(A, x, n)]
        moves += [(A, C - {x}) for x in C if corner(C, x, n)]
        answers = []
        for a, c in moves:
            child_blocks = components(a, c, n)
            # A newly exclusive core word has no parent variable: existential
            # elimination is encoded by dropping that word from this map.
            maps = [{parent_ids[x] for x in block if x in parent_ids}
                    for block in child_blocks]
            for p in solve(frozenset(a), frozenset(c)):
                groups = maps + [set().union(*(maps[i] for i in b)) for b in p]
                q = equalities_partition(k, groups)
                if q == discrete:
                    return (discrete,)  # All labels accepted; no further search needed.
                if any(refines(old, q) for old in answers):
                    continue
                answers = [old for old in answers if not refines(q, old)] + [q]
        return tuple(sorted(answers))
    return solve


def accepted(partitions, labels):
    return any(all(len({labels[i] for i in block}) <= 1 for block in p)
               for p in partitions)


def check_local_rules(A, C, labels, n):
    blocks = components(A, C, n)
    ids = {x: i for i, b in enumerate(blocks) for x in b}
    X = lift(A, C, blocks, labels, n)
    for x in A:
        for b in range(2):
            word = x | (b << n)
            if word not in X:
                continue
            if x not in C:
                expected = corner(A, x, n)
            else:
                neighbor_ids = {ids[x ^ (1 << e)] for e in range(n)
                                if x ^ (1 << e) in A - C}
                expected = corner(C, x, n) and all(labels[i] != b for i in neighbor_ids)
            assert corner(X, word, n + 1) == expected
    return X


def main():
    source_path = ROOT / 'papers/minor-atlas/research/coordinate_block_amplification_verification.json'
    O = frozenset(json.loads(source_path.read_text())['source_family'])
    small = []
    for e in range(5):
        A, C = contraction(O, e), section(O, e, 0) & section(O, e, 1)
        solve = partition_recursion(4)
        partitions = solve(A, C)
        blocks = components(A, C, 4)
        checks = []
        for t in range(1 << len(blocks)):
            labels = [(t >> i) & 1 for i in range(len(blocks))]
            X = check_local_rules(A, C, labels, 4)
            order = greedy_order(X, 5)
            assert order is not None and accepted(partitions, labels)
            checks.append(dict(labels=labels, order=order))
        assert partitions == (((0,), (1,)),)
        small.append(dict(coordinate=e, partitions=partitions, states=solve.cache_info().currsize,
                          checks=checks))
    seed_path = ROOT / 'papers/dismantlable-ample/evidence/order267_obstruction/record.json'
    seed = json.loads(seed_path.read_text())
    K, n = frozenset(seed['family']), seed['n']
    seed_checks = []
    for e in range(n):
        A = contraction(K, e)
        C = section(K, e, 0) & section(K, e, 1)
        assert len(A) == len(shadow(A, 11)) == 208
        assert len(C) == len(shadow(C, 11)) == 59
        blocks = components(A, C, 11)
        assert len(blocks) == 2
        ids = {x: i for i, block in enumerate(blocks) for x in block}
        outside_corners = [x for x in A - C if corner(A, x, 11)]
        core_corners = [x for x in C if corner(C, x, 11)]
        requirements = [{ids[x ^ (1 << f)] for f in range(11)
                         if x ^ (1 << f) in A - C} for x in core_corners]
        assert not outside_corners and len(core_corners) == 3
        assert all(s == {0, 1} for s in requirements)
        a_order, c_order = greedy_order(A, 11), greedy_order(C, 11)
        assert a_order is not None and c_order is not None
        labels = []
        for block in blocks:
            colors = {b for x in block for b in (0, 1)
                      if (x | (b << 11)) in
                      {y & ((1 << e) - 1) | ((y >> (e + 1)) << e)
                       | (((y >> e) & 1) << 11) for y in K}}
            assert len(colors) == 1
            labels.append(colors.pop())
        assert labels[0] != labels[1]
        for word in range(4):
            colors = [word & 1, (word >> 1) & 1]
            X = check_local_rules(A, C, colors, 11)
            count = sum(corner(X, x, 12) for x in X)
            assert count == (3 if word in (0, 3) else 0)
        seed_checks.append(dict(coordinate=e, contraction=sorted(A), overlap=sorted(C),
                                components=[sorted(b) for b in blocks], actual_labels=labels,
                                core_corners=core_corners,
                                first_move_requirements=[sorted(s) for s in requirements],
                                exact_positive_partition=[[0, 1]],
                                contraction_order=a_order, overlap_order=c_order))
    result = dict(status='passed', small_controls=small, seed_checks=seed_checks,
                  scope='Exact small partition recursions; all 12 seed pairs classified by first-corner equality obstruction and positive input orders. No exponential partition search on 208-concept inputs.',
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                                for p in (source_path, seed_path,
                                          ROOT / 'papers/cover-recursive/research/component_label_recursion.py')})
    Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2) + '\n')
    print('PASS: exact small partition recursion and all local corner rules.', flush=True)
    print('PASS: all 12 seed pairs have exactly the equality t0=t1;', flush=True)
    print('positive contraction/overlap orders replay; opposite labels are cornerless.', flush=True)


if __name__ == '__main__':
    main()
