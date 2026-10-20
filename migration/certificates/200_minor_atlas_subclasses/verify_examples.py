#!/usr/bin/env python3
"""Check the fixed atlas witnesses directly from finite cube definitions."""
import itertools
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def submasks(mask):
    sub = mask
    while True:
        yield sub
        if not sub:
            return
        sub = (sub - 1) & mask


def maximal_cubes(words, n):
    cubes = set()
    strong = set()
    for free in range(1 << n):
        for fixed in submasks(((1 << n) - 1) ^ free):
            vertices = frozenset(fixed | s for s in submasks(free))
            if vertices <= words:
                cubes.add(vertices)
                strong.add(free)
    maximal = [q for q in cubes if not any(q < r for r in cubes)]
    return sorted(maximal, key=lambda q: sorted(q)), strong


def chordal(cubes):
    neighbors = {i: {j for j, r in enumerate(cubes) if i != j and q & r}
                 for i, q in enumerate(cubes)}
    remaining = set(neighbors)
    while remaining:
        for i in sorted(remaining):
            ns = neighbors[i] & remaining
            if all(b in neighbors[a] for a, b in itertools.combinations(ns, 2)):
                remaining.remove(i)
                break
        else:
            return False
    return True


def join_tree(cubes, words):
    # A maximum-weight spanning tree of hyperedges is a join tree iff
    # each word's occurrence nodes are connected (the running-intersection test).
    parent = list(range(len(cubes)))
    def root(i):
        while parent[i] != i:
            i = parent[i]
        return i
    edges = sorted(((len(a & b), i, j) for i, a in enumerate(cubes)
                    for j, b in enumerate(cubes) if i < j), reverse=True)
    tree = {i: set() for i in range(len(cubes))}
    for _, i, j in edges:
        a, b = root(i), root(j)
        if a != b:
            parent[a] = b
            tree[i].add(j)
            tree[j].add(i)
    for v in words:
        nodes = {i for i, q in enumerate(cubes) if v in q}
        reached = {next(iter(nodes))}
        pending = list(reached)
        while pending:
            for j in tree[pending.pop()] & nodes - reached:
                reached.add(j)
                pending.append(j)
        if reached != nodes:
            return False
    return True


def properties(words, n):
    words = set(words)
    cubes, strong = maximal_cubes(words, n)
    shattered = {s for s in range(1 << n)
                 if {v & s for v in words} == set(submasks(s))}
    daisy = any(all(all((r ^ s) in words for s in submasks(r ^ v))
                    for v in words) for r in words)
    remainder = set(words)
    while remainder:
        for v in sorted(remainder):
            directions = sum(1 << i for i in range(n) if (v ^ (1 << i)) in remainder)
            if all((v ^ s) in remainder for s in submasks(directions)):
                remainder.remove(v)
                break
        else:
            break
    return dict(order=len(words), maximal_cubes=len(cubes), ample=shattered == strong,
                cube_chordal=chordal(cubes), daisy=daisy,
                corner_peelable=not remainder,
                median=all(((a & b) | (a & c) | (b & c)) in words
                           for a, b, c in itertools.combinations(words, 3)),
                maximal_cube_join_tree=join_tree(cubes, words))


def verify():
    inputs = json.loads((ROOT / 'inputs.json').read_text())
    report = {}
    for name, case in inputs.items():
        result = properties(case['words'], case['dimension'])
        for key, expected in case['expected'].items():
            assert result[key] == expected, (name, key, result[key], expected)
        report[name] = result
    projected = {v & 15 for v in inputs['cube_chordal_projection_counterexample']['words']}
    assert projected == set(inputs['non_cube_chordal_dismantlable_ten_words']['words'])
    return dict(status='passed', scope='Fixed finite witnesses only; no universal theorem is inferred.',
                projection_to_ten_word_class=True, cases=report)


if __name__ == '__main__':
    print(json.dumps(verify(), indent=2, sort_keys=True))
