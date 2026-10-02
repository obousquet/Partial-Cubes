"""Self-contained orientation primitives used by the joint-switch replay.

These functions are copied without mathematical change from
``orientation_gluing.py``.  Keeping this narrow module avoids importing the
unrelated large-control replay dependencies of that source-environment file.
"""

import heapq


def orientation(order, n):
    rank = {x: i for i, x in enumerate(order)}
    return tuple(sorted((x, x ^ (1 << e)) for x in order for e in range(n)
                        if x ^ (1 << e) in rank and rank[x] < rank[x ^ (1 << e)]))


def topsort(vertices, edges):
    out = {x: set() for x in vertices}
    indeg = dict.fromkeys(vertices, 0)
    for x, y in edges:
        if y not in out[x]:
            out[x].add(y)
            indeg[y] += 1
    ready = [x for x in vertices if not indeg[x]]
    heapq.heapify(ready)
    order = []
    while ready:
        x = heapq.heappop(ready)
        order.append(x)
        for y in sorted(out[x]):
            indeg[y] -= 1
            if not indeg[y]:
                heapq.heappush(ready, y)
    return order if len(order) == len(vertices) else None


def descendants(vertices, edges):
    order = topsort(vertices, edges)
    assert order is not None
    successors = {y: set() for y in vertices}
    for x, y in edges:
        successors[x].add(y)
    reach = {}
    for y in reversed(order):
        reach[y] = {y}.union(*(reach[z] for z in successors[y]))
    return reach


def constraints(A, C, da, dc, signs, endpoint=None):
    reach = descendants(C, dc)
    added = []
    invalid = []
    for x, y in da:
        if signs[x] == signs[y]:
            continue
        if y not in C:
            invalid.append([x, y])
        else:
            added.extend((x, z, x, y) for z in sorted(reach[y]))
    edges = set(da) | {(x, z) for x, z, _, _ in added}
    if endpoint is not None:
        r, b = endpoint
        if signs[r] != b:
            invalid.append(['endpoint', r, b])
        edges.update((x, r) for x in A if x != r)
    alpha = None if invalid else topsort(A, edges)
    return alpha, edges, added, invalid
