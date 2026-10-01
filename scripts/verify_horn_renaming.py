"""Bounded exact diagnostic for the Horn-renaming pair defect.

Nonempty classes through Q3; elementary minors; products on Q2 x Q2;
rooted wedges on Q2 x Q3; the inherited 13-concept Q6 example.
No symmetry quotient, randomized sampling, or larger search.
"""
from collections import deque
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1] / 'data' / 'examples' / 'assets'


@lru_cache(None)
def patterns(n, C):
    """(support mask, prescribed one mask); minimality checked by deletions."""
    C = tuple(C)
    answer = []
    for support in range(1, 1 << n):
        present = {c & support for c in C}
        bits = [i for i in range(n) if support >> i & 1]
        for choices in product((0, 1), repeat=len(bits)):
            word = sum(b << i for i, b in zip(bits, choices))
            if word in present:
                continue
            if all(any((c & (support ^ (1 << i))) ==
                       (word & (support ^ (1 << i))) for c in C) for i in bits):
                answer.append((support, word))
    return tuple(answer)


def cost(rows, orientation):
    counts = [((word ^ orientation) & support).bit_count()
              for support, word in rows]
    return sum(k * (k - 1) // 2 for k in counts)


@lru_cache(None)
def values(n, C):
    rows = patterns(n, C)
    return tuple(cost(rows, v) for v in range(1 << n))


def defect(n, C):
    return min(values(n, C))


def union_closed(C):
    C = set(C)
    return all(a | b in C for a in C for b in C)


def implication_test(n, rows):
    """Reachability version of the classical 2-SAT criterion, with paths."""
    edges = [[] for _ in range(2*n)]
    for p, (support, word) in enumerate(rows):
        for i, j in combinations([i for i in range(n) if support >> i & 1], 2):
            a, b = 2*i + (word >> i & 1), 2*j + (word >> j & 1)
            edges[a ^ 1].append((b, p, i, j))
            edges[b ^ 1].append((a, p, i, j))

    def path(start, end):
        queue = deque([start]); previous = {start: None}
        while queue:
            x = queue.popleft()
            if x == end:
                result = []
                while previous[x] is not None:
                    old, p, i, j = previous[x]
                    result.append([old, x, p, i, j]); x = old
                return result[::-1]
            for y, p, i, j in edges[x]:
                if y not in previous:
                    previous[y] = (x, p, i, j); queue.append(y)
        return None

    for i in range(n):
        a, b = path(2*i, 2*i+1), path(2*i+1, 2*i)
        if a is not None and b is not None:
            return {'satisfiable': False, 'coordinate': i, 'paths': [a, b]}
    return {'satisfiable': True}


def classes(n, rooted=False):
    for mask in range(1, 1 << (1 << n)):
        if rooted and not mask & 1:
            continue
        yield tuple(c for c in range(1 << n) if mask >> c & 1)


def drop(word, i):
    return (word & ((1 << i)-1)) | ((word >> (i+1)) << i)


def ample(n, C):
    C = set(C)
    for s in range(1 << n):
        if len({c & s for c in C}) != 1 << s.bit_count():
            continue
        if not any(all(c ^ t in C for t in range(1 << n) if t & ~s == 0)
                   for c in C):
            return False
    return True


def run():
    counts = dict(classes=0, orientations=0, elementary_minors=0,
                  products=0, rooted_wedges=0)
    zero_counts = {}
    for n in range(4):
        zero_counts[n] = 0
        for C in classes(n):
            counts['classes'] += 1
            scores = values(n, C); h = min(scores)
            for v, score in enumerate(scores):
                assert (score == 0) == union_closed([c ^ v for c in C])
                counts['orientations'] += 1
            assert implication_test(n, patterns(n, C))['satisfiable'] == (h == 0)
            zero_counts[n] += h == 0
            for i in range(n):
                for fixed in (None, 0, 1):
                    D = tuple(sorted({drop(c, i) for c in C
                                      if fixed is None or c >> i & 1 == fixed}))
                    if D:
                        assert defect(n-1, D) <= h
                        counts['elementary_minors'] += 1
    for A, B in product(classes(2), repeat=2):
        C = tuple(sorted(a | b << 2 for a in A for b in B))
        assert defect(4, C) == defect(2, A) + defect(2, B)
        counts['products'] += 1
    for A, B in product(classes(2, True), classes(3, True)):
        C = tuple(sorted(set(A) | {b << 2 for b in B}))
        activeA = 0; activeB = 0
        for a in A: activeA |= a
        for b in B: activeB |= b
        expected = tuple(values(2, A)[v & 3] + values(3, B)[v >> 2]
                         + (activeA & ~v).bit_count()
                         * (activeB & ~(v >> 2)).bit_count() for v in range(32))
        assert values(5, C) == expected
        counts['rooted_wedges'] += 1
    # Source writes s=110 in named coordinate order; low-bit encoding is 3.
    D = tuple(range(7)); s = 3
    C = tuple(sorted({d | s << 3 for d in D} | {s | d << 3 for d in D}))
    scores = values(6, C); rows = patterns(6, C)
    assert len(C) == 13 and ample(6, C)
    assert len(rows) == 11
    witness = scores.index(min(scores))
    obstruction = implication_test(6, rows)
    assert min(scores) > 0 and not obstruction['satisfiable']
    result = {
        'status': 'passed', 'date': '2026-09-13', 'counts': counts,
        'zero_counts_by_dimension': zero_counts,
        'calibration': {'words': C, 'minimal_missing_patterns': rows,
                        'defect': min(scores), 'optimal_orientation': witness,
                        'all_orientation_costs': scores,
                        'contradiction': obstruction},
        'scope': 'Exact finite diagnostics. General laws use the written proof; '
                 'MaxMin membership is imported from the named source, not tested here.',
        'script_sha256': sha256(Path(__file__).read_bytes()).hexdigest()
    }
    out = ROOT/'horn_renaming_verification.json'
    out.write_text(json.dumps(result, indent=2)+'\n')
    saved = json.loads(out.read_text())
    # Replay every directed edge with its signed-pattern provenance.
    for path in saved['calibration']['contradiction']['paths']:
        for x, y, p, i, j in path:
            support, word = saved['calibration']['minimal_missing_patterns'][p]
            assert support >> i & 1 and support >> j & 1
            a, b = 2*i + (word >> i & 1), 2*j + (word >> j & 1)
            assert (x, y) in ((a ^ 1, b), (b ^ 1, a))
        assert all(path[k][1] == path[k+1][0] for k in range(len(path)-1))
    a,b = saved['calibration']['contradiction']['paths']
    assert a[0][0] == b[-1][1] and a[-1][1] == b[0][0]
    assert a[0][0] ^ a[-1][1] == 1
    print(json.dumps({'status': 'passed', 'counts': counts,
                      'calibration_defect': min(scores),
                      'optimal_orientation': witness}, indent=2))


if __name__ == '__main__':
    run()
