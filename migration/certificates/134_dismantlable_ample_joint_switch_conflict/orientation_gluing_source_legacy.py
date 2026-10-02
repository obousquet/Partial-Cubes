"""Eliminate both input total orders for fixed peeling orientations and bits.

Compare every orientation/bit choice with independently enumerated parent
peelings on Q3. Also reconstruct one archived 577-word positive witness.
No enumeration of orientations or bits is attempted on the large input.
"""
from pathlib import Path
from itertools import product
import json
import heapq
from latest_overlap_schedule import peelings, latest
from orientation_deadlines import normalize
from filtration_lifting import lift_heights, replay_heights
from ample_filtrations import transport


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


def strong_components(vertices, edges):
    """Kosaraju, used only to verify the color-independent decomposition."""
    out, rev = {x: [] for x in vertices}, {x: [] for x in vertices}
    for x, y in edges:
        out[x].append(y)
        rev[y].append(x)
    seen, finish = set(), []
    def visit(x):
        seen.add(x)
        for y in out[x]:
            if y not in seen:
                visit(y)
        finish.append(x)
    for x in sorted(vertices):
        if x not in seen:
            visit(x)
    seen = set()
    blocks = []
    for x in reversed(finish):
        if x in seen:
            continue
        todo, block = [x], {x}
        seen.add(x)
        while todo:
            x = todo.pop()
            for y in rev[x]:
                if y not in seen:
                    seen.add(y)
                    block.add(y)
                    todo.append(y)
        blocks.append(sorted(block))
    return sorted(blocks)


def descendants(C, dc):
    gamma = topsort(C, dc)
    assert gamma is not None
    successors = {y: set() for y in C}
    for x, y in dc:
        successors[x].add(y)
    reach = {}
    for y in reversed(gamma):
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


def generate_assignments(A, C, da, dc):
    """Direct edge-selection/partition generator; no parent recognition call."""
    reach = descendants(C, dc)
    eligible = [(x, y) for x, y in da if y in C]
    generated = set()
    for bits in product((0, 1), repeat=len(eligible)):
        chosen = {e for e, b in zip(eligible, bits) if b}
        edges = set(da) | {(x, z) for x, y in chosen for z in reach[y]}
        if topsort(A, edges) is None:
            continue
        # Equality components of the remaining undirected edges.
        adj = {x: set() for x in A}
        for x, y in set(da)-chosen:
            adj[x].add(y)
            adj[y].add(x)
        left = set(A)
        blocks = []
        while left:
            todo = [min(left)]
            block = set(todo)
            left -= block
            while todo:
                x = todo.pop()
                new = adj[x] & left
                block |= new
                left -= new
                todo.extend(new)
            blocks.append(block)
        for colors in product((0, 1), repeat=len(blocks)):
            signs = {x: b for B, b in zip(blocks, colors) for x in B}
            generated.add(tuple(signs[x] for x in sorted(A)))
    return generated


def reconstruct(A, C, da, dc, signs, n):
    alpha, edges, added, invalid = constraints(A, C, da, dc, signs)
    if alpha is None:
        return None
    gamma = topsort(C, dc)
    gamma, _, _ = normalize(alpha, gamma, n)
    a, c, _, _ = latest(alpha, gamma)
    h = lift_heights(a, c, signs, n)
    parent = replay_heights(h, n+1)
    assert orientation(alpha, n) == da and orientation(gamma, n) == dc
    return parent


def cube_uso(da, dc, signs, n):
    """Check the independently assembled outmap on a full cube prism."""
    A = range(1 << n)
    oa, oc = dict.fromkeys(A, 0), dict.fromkeys(A, 0)
    for x, y in da:
        oa[x] |= x ^ y
    for x, y in dc:
        oc[x] |= x ^ y
    out = {}
    for x in A:
        out[x | (signs[x] << n)] = oa[x]
        out[x | ((1-signs[x]) << n)] = oc[x] | (1 << n)
    for x in out:
        for e in range(n+1):
            y = x ^ (1 << e)
            if bool(out[x] & (1 << e)) == bool(out[y] & (1 << e)):
                return None
    # Every cube face has exactly one sink. Outgoing-cube containment is
    # automatic here because the support is the full cube.
    for face in product((-1, 0, 1), repeat=n+1):
        free = sum(1 << e for e, b in enumerate(face) if b == -1)
        words = [x for x in out if all(b == -1 or ((x >> e) & 1) == b
                                      for e, b in enumerate(face))]
        if sum(not (out[x] & free) for x in words) != 1:
            return None
    return out


def main():
    A = C = frozenset(range(4))
    n = 2
    # Independent truth set: enumerate actual parent corner orders, then
    # project max/min times and record edge orientations and upper layers.
    truth = set()
    endings = {}
    parent_count = 0
    for order in peelings(frozenset(range(8)), 3):
        parent_count += 1
        h = {x: i for i, x in enumerate(order)}
        a, c = transport(h, n, 'max'), transport(h, n, 'min')
        da = orientation(sorted(a, key=a.get), n)
        dc = orientation(sorted(c, key=c.get), n)
        signs = tuple(int(h[x | (1 << n)] > h[x]) for x in sorted(A))
        truth.add((da, dc, signs))
        endings.setdefault((da, dc, signs), set()).add(order[-1])
    os = sorted({orientation(order, n) for order in peelings(A, n)})
    checks = accepted = 0
    cyclic_uso = None
    for da, dc, bits in product(os, os, product((0, 1), repeat=4)):
        signs = dict(zip(sorted(A), bits))
        alpha, edges, added, invalid = constraints(A, C, da, dc, signs)
        expected = (da, dc, bits) in truth
        assert (alpha is not None) == expected
        checks += 1
        accepted += expected
        for r, b in product(A, (0, 1)):
            rooted, _, _, _ = constraints(A, C, da, dc, signs, (r, b))
            assert (rooted is not None) == (r | (b << n) in endings.get((da, dc, bits), set()))
        if expected:
            reconstruct(A, C, da, dc, signs, n)
        elif cyclic_uso is None:
            out = cube_uso(da, dc, signs, n)
            if out is not None:
                parent_edges = [(x, x ^ (1 << e)) for x in out for e in range(n+1)
                                if out[x] & (1 << e)]
                assert topsort(out, parent_edges) is None
                cyclic_uso = dict(contraction_orientation=da, overlap_orientation=dc,
                    upper_layers=signs, parent_outmap=out, parent_edges=parent_edges,
                    augmented_edges=sorted(edges), added_edges_with_cut_witness=added)
    assert cyclic_uso is not None
    for da, dc in product(os, repeat=2):
        generated = generate_assignments(A, C, da, dc)
        assert generated == {bits for aa, cc, bits in truth if aa == da and cc == dc}
    # A proper overlap tests forced exclusive labels and outward cut edges.
    # Every assignment here is ample: the two exclusive vertices are separate
    # components of A\C. Truth again comes from parent corner orders.
    P, K, pn = frozenset([0, 1, 3]), frozenset([1]), 2
    path_checks = path_accepted = outward_rejections = 0
    pos = {orientation(order, pn) for order in peelings(P, pn)}
    for bits in product((0, 1), repeat=3):
        signs = dict(zip(sorted(P), bits))
        X = frozenset([1, 1 | (1 << pn)] +
                      [x | (signs[x] << pn) for x in P-K])
        actual = set()
        for order in peelings(X, pn+1):
            h = {x: i for i, x in enumerate(order)}
            a, c = transport(h, pn, 'max'), transport(h, pn, 'min')
            upper = tuple(int(h.get(x | (1 << pn), -1) == a[x]) for x in sorted(P))
            actual.add((orientation(sorted(a, key=a.get), pn), upper))
        for da in pos:
            alpha, _, _, invalid = constraints(P, K, da, (), signs)
            assert (alpha is not None) == ((da, bits) in actual)
            path_checks += 1
            path_accepted += alpha is not None
            outward_rejections += bool(invalid)
    for da in pos:
        generated = generate_assignments(P, K, da, ())
        expected = {bits for bits in product((0, 1), repeat=3)
                    if constraints(P, K, da, (), dict(zip(sorted(P), bits)))[0] is not None}
        assert generated == expected
    path = Path(__file__).with_name('orientation_deadlines.json')
    record = json.loads(path.read_text())['trap_records'][0]
    order = record['successful_lift_order']
    data = json.loads(Path(__file__).with_name('cube_corner_compatibility.json').read_text())
    A, C, n = frozenset(data['A']), frozenset(data['C']), data['dimension']
    h = {x: i for i, x in enumerate(order)}
    a, c = transport(h, n, 'max'), transport(h, n, 'min')
    da, dc = orientation(sorted(a, key=a.get), n), orientation(sorted(c, key=c.get), n)
    signs = {x: (1 if h.get(x | (1 << n), -1) == a[x] else 0) for x in A}
    lifted = reconstruct(A, C, da, dc, signs, n)
    assert lifted is not None and set(lifted) == set(order)
    alpha, edges, added, invalid = constraints(A, C, da, dc, signs)
    blocks = strong_components(A, set(da) | set(dc))
    assert all(topsort(B, [(x, y) for x, y in edges if x in B and y in B])
               is not None for B in blocks)
    result = dict(status='PASS: exact fixed-orientation gluing; full obstruction generation remains open',
        cube_parent_peelings=parent_count, cube_orientation_count=len(os),
        cube_choices_checked=checks, cube_choices_accepted=accepted,
        cube_endpoint_checks=8*checks,
        partition_generator_orientation_pairs=len(os)**2 + len(pos),
        path_control=dict(choices_checked=path_checks, accepted=path_accepted,
                          outward_edge_rejections=outward_rejections),
        cyclic_uso_control=cyclic_uso,
        large_control=dict(contraction_orientation=da, overlap_orientation=dc,
            upper_layers=signs, augmented_edges=sorted(edges),
            added_edges_with_cut_witness=added, reconstructed_order=lifted,
            union_strong_components=blocks,
            changed_contraction_order=alpha != sorted(a, key=a.get)))
    Path(__file__).with_suffix('.json').write_text(json.dumps(result, indent=2)+'\n')
    print('PASS Q3:', parent_count, 'parent peelings;', checks, 'choices;', accepted, 'accepted')
    print('PASS: cyclic USO with acyclic inputs; large lift', len(lifted),
          'words, contraction order changed:', result['large_control']['changed_contraction_order'])


if __name__ == '__main__':
    main()
