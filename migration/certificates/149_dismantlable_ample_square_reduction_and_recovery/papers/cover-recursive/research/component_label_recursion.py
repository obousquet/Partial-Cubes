"""Check the uncoloured face circuit for peripheral component expansions.

Bounded scope: all 64 clause-sign splits at coordinate 0 of P5, all
18 elementary minors of each ample output, and explicit positive/negative
component-expansion controls. No claim about higher-dimensional coverage.
Run from the repository root with the memory/CPU limits in AGENTS.md.
"""
from collections import Counter
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
from explore_rper_cyclic_families import contraction, normalize, rper, section


def components(A, C, n):
    left = set(A - C)
    out = []
    while left:
        todo = [min(left)]
        left.remove(todo[0])
        found = set(todo)
        while todo:
            x = todo.pop()
            for e in range(n):
                y = x ^ (1 << e)
                if y in left:
                    left.remove(y)
                    found.add(y)
                    todo.append(y)
        out.append(frozenset(found))
    return out


def circuit(A, C, n):
    """Build without evaluating R_per membership or constructing coloured lifts.

    Correctness requires A,C in R_per (empty allowed). The caller supplies
    that hypothesis separately. Fixed-face vertices retain original bit indices.
    """
    assert C <= A
    blocks = components(A, C, n)
    component = {x: i for i, block in enumerate(blocks) for x in block}
    nodes = []

    @lru_cache(None)
    def face(remaining, fixed):
        mask = ((1 << n) - 1) ^ remaining
        a = frozenset(x for x in A if x & mask == fixed)
        c = frozenset(x for x in C if x & mask == fixed)
        mono = sorted({component[x] for x in a - c})
        cuts = []
        for e in range(n):
            if not remaining >> e & 1:
                continue
            bit = 1 << e
            orientations = [b for b in (0, 1)
                            if all(x ^ bit in a for x in a if (x >> e) & 1 == b)
                            and all(x ^ bit in c for x in c if (x >> e) & 1 == b)]
            if orientations:
                cuts.append([e, orientations[0], face(remaining ^ bit, fixed),
                             face(remaining ^ bit, fixed | bit)])
        index = len(nodes)
        nodes.append(dict(remaining=remaining, fixed=fixed, mono=mono, cuts=cuts))
        return index

    root = face((1 << n) - 1, 0)
    return dict(dimension=n, components=[sorted(x) for x in blocks],
                nodes=nodes, root=root)


def evaluate(graph, labels):
    values = []
    for node in graph['nodes']:
        mono = len({labels[i] for i in node['mono']}) <= 1
        values.append(mono or any(values[a] and values[b]
                                 for _, _, a, b in node['cuts']))
    return values[graph['root']]


def lift(A, C, blocks, labels, n):
    return frozenset(set(C) | {x | (1 << n) for x in C}
                     | {x | (labels[i] << n) for i, block in enumerate(blocks)
                        for x in block})


def shadow(H, n):
    return {s for s in range(1 << n)
            if len({x & s for x in H}) == 1 << s.bit_count()}


def certificate(H, n):
    """Produce a peripheral construction tree; never accepts a negative family."""
    H, n = normalize(H, n)
    if len(H) <= 1:
        return dict(dimension=n, family=sorted(H), leaf=True)
    for e in range(n):
        a, b = section(H, e, 0), section(H, e, 1)
        for small, big, side in ((a, b, 0), (b, a, 1)):
            if small <= big and rper(*normalize(small, n - 1)) and rper(*normalize(big, n - 1)):
                return dict(dimension=n, family=sorted(H), coordinate=e, small_side=side,
                            site=certificate(small, n - 1), base=certificate(big, n - 1))
    raise AssertionError('No positive construction certificate')


def compare_all_labels(A, C, n):
    assert rper(*normalize(A, n)) and rper(*normalize(C, n))
    graph = circuit(A, C, n)
    blocks = [frozenset(x) for x in graph['components']]
    accepted = []
    for word in range(1 << len(blocks)):
        labels = [(word >> i) & 1 for i in range(len(blocks))]
        X = lift(A, C, blocks, labels, n)
        assert len(X) == len(shadow(X, n + 1))
        predicted = evaluate(graph, labels)
        assert predicted == rper(*normalize(X, n + 1))
        if predicted:
            accepted.append(word)
    return dict(graph=graph, accepted=accepted, labelings=1 << len(blocks))


def main():
    source = ROOT / 'papers/minor-atlas/research/coordinate_block_amplification_verification.json'
    record = json.loads(source.read_text())
    O = frozenset(record['source_family'])
    rows = record['source_clauses']
    relevant = [i for i, (s, t) in enumerate(rows) if s & 1]
    counts = Counter()
    cases = []
    comparisons = 0
    example = None
    for twist in range(64):
        clauses = [(s | 32, t | (((t & 1) ^ ((twist >> relevant.index(i)) & 1)) << 5))
                   if s & 1 else (s, t) for i, (s, t) in enumerate(rows)]
        X = frozenset(x for x in range(64) if all(x & s != t for s, t in clauses))
        XS = shadow(X, 6)
        if len(X) != len(XS):
            counts['nonample'] += 1
            continue
        counts['ample'] += 1
        A = contraction(X, 5)
        C = section(X, 5, 0) & section(X, 5, 1)
        assert C == O and not rper(*normalize(X, 6))
        positive_inputs = all(rper(*normalize(M, 5)) for M in
                              (A, section(X, 5, 0), section(X, 5, 1)))
        failures = []
        checks = []
        for e in range(5):
            for op in (0, 1, 2):
                mu = lambda H: contraction(H, e) if op == 2 else section(H, e, op)
                Y, a, c = mu(X), mu(A), mu(C)
                # Commutation is checked from actual families, independently.
                assert contraction(Y, 4) == a
                assert section(Y, 4, 0) & section(Y, 4, 1) == c
                if not rper(*normalize(a, 4)):
                    failures.append([e, op, 'negative contraction input'])
                    continue
                comparison = compare_all_labels(a, c, 4)
                comparisons += comparison['labelings']
                labels = []
                for block in comparison['graph']['components']:
                    colors = {b for x in block for b in (0, 1) if x | (b << 4) in Y}
                    assert len(colors) == 1
                    labels.append(colors.pop())
                ok = evaluate(comparison['graph'], labels)
                assert ok == rper(*normalize(Y, 5))
                checks.append(dict(coordinate=e, operation=op, labels=labels,
                                   circuit=ok, face_nodes=len(comparison['graph']['nodes'])))
                if not ok:
                    failures.append([e, op, 'label circuit'])
        generated = positive_inputs and not failures
        direct = all(rper(*normalize(M, 5)) for e in range(6) for M in
                     (section(X, e, 0), section(X, e, 1), contraction(X, e)))
        assert generated == direct
        counts['forbidden' if generated else 'nonminimal'] += 1
        nf = [s for s in range(64) if s not in XS and
              all(s ^ (1 << e) in XS for e in range(6) if s >> e & 1)]
        canonical = [(s, next(t for t in range(64) if t & ~s == 0
                             and t not in {x & s for x in X})) for s in nf]
        columns = [tuple(0 if not s >> e & 1 else 2 * ((t >> e) & 1) - 1
                         for s, t in canonical) for e in range(6)]
        distinct = len({min(v, tuple(-a for a in v)) for v in columns})
        if generated and distinct == 6:
            counts['column_distinct_forbidden'] += 1
        cases.append(dict(twist=twist, order=len(X), forbidden=generated,
                          distinct_columns=distinct, positive_inputs=positive_inputs,
                          failures=failures, checks=checks))
        if twist == 2:
            assert generated and len(X) == 37 and distinct == 6
            example = dict(family=sorted(X), dimension=6, order=37,
                           vc_dimension=max(s.bit_count() for s in XS), clauses=canonical,
                           overlap=sorted(C), contraction=sorted(A),
                           positive_input_certificates=[certificate(M, 5) for M in
                                                       (A, section(X, 5, 0), section(X, 5, 1))])
    controls = []
    for e in range(5):
        A = contraction(O, e)
        C = section(O, e, 0) & section(O, e, 1)
        check = compare_all_labels(A, C, 4)
        assert check['accepted'] == [0, 3]
        controls.append(check)
        comparisons += check['labelings']
    star = frozenset([0, 1, 2, 4])
    check = compare_all_labels(star, frozenset([0]), 3)
    assert check['accepted'] == list(range(8))
    controls.append(check)
    comparisons += check['labelings']
    assert counts == dict(nonample=16, ample=48, forbidden=32, nonminimal=16,
                          column_distinct_forbidden=30)
    result = dict(status='passed', counts=dict(counts), comparisons=comparisons,
                  cases=cases, example=example, controls=controls,
                  scope='All 64 P5 sign splits; direct elementary minors and all labels for their positive minor pairs; no general finite census.',
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  input_sha256={str(source.relative_to(ROOT)):hashlib.sha256(source.read_bytes()).hexdigest()})
    Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2) + '\n')
    print('PASS:', dict(counts), flush=True)
    print('Independent label comparisons:', comparisons, flush=True)
    print('Column-distinct 37-concept obstruction strongly projects to P5.', flush=True)


if __name__ == '__main__':
    main()
