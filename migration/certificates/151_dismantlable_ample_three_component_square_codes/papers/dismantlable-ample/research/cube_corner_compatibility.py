"""Fixed cube-fiber incompatibility witness, not a minimal-obstruction claim.

Reconstruct the ample pair from the archived 580-word construction's
cornerless deletion leaf. Replay all positive orders; certify negatives by
absence of corners, never by a failed greedy search. Run with AGENTS limits.
"""
from pathlib import Path
import hashlib, json, sys
ROOT = Path(__file__).resolve().parents[3]
sys.path[:0] = [str(ROOT), str(Path(__file__).parent),
                str(ROOT / 'papers/cover-recursive/research')]
from component_label_recursion import components, lift
from component_label_partitions import corner, greedy_order, replay
from cube_code_expansion import cube_code
from explore_rper_cyclic_families import section, contraction



def expand_binary_trace(A, C, labels, order, n, d, codes):
    """Extract equalities from a binary order, then realize its full-word trace.

    Reject incompatible codes rather than pretending independently successful
    bit orders suffice. Fresh variables are eliminated after recording the trace.
    """
    ks = components(A,C,n)
    ids = {x:i for i,K in enumerate(ks) for x in K}
    parent = list(range(len(ks)))
    def root(i):
        while parent[i] != i:
            i = parent[i]
        return i
    X = lift(A,C,ks,labels,n)
    replay(X,order,n+1)
    a,c = set(A),set(C)
    moves = []
    for word in order:
        x = word & ((1<<n)-1)
        if x not in c:
            assert corner(a,x,n)
            moves.append(('exclusive',x,ids.pop(x)))
            a.remove(x)
        else:
            assert corner(c,x,n)
            adjacent = [ids[x^(1<<f)] for f in range(n) if x^(1<<f) in a-c]
            if adjacent:
                variable = root(adjacent[0])
                for i in adjacent[1:]:
                    parent[root(i)] = variable
            else:
                variable = len(parent); parent.append(variable)
            ids[x] = variable
            moves.append(('core',x,variable))
            c.remove(x)
    assert len(a) == 1 and not c
    values = {}
    for i,q in enumerate(codes):
        r = root(i)
        if r in values and values[r] != q:
            return None
        values[r] = q
    value = lambda i: values.get(root(i),0)
    output = []
    def punctured(q):
        return sorted((y for y in range(1<<d) if y != q),
                      key=lambda y: ((y^q).bit_count(),y))
    for kind,x,i in moves:
        q = value(i)
        if kind == 'core':
            output.append(x | (q<<n))
        else:
            output.extend(x | (y<<n) for y in punctured(q))
    last = next(iter(a)); q = value(ids[last])
    output.extend(last | (y<<n) for y in punctured(q)[:-1])
    Z = cube_code(A,C,n,d,codes)
    replay(Z,output,n+d)
    return dict(dimension=d, codes=codes, peeling=output,
                trace=[list(m) for m in moves],
                partition=[[i for i in range(len(ks)) if root(i)==r]
                           for r in sorted({root(i) for i in range(len(ks))})])


def main():
    source = Path(__file__).with_name('root_cube_completion.json')
    data = json.loads(source.read_text())['cube_completion']
    leaves = [frozenset(r['family']) for r in data['negative_states'] if not r['corners']]
    assert len(leaves) == 1
    H = leaves[0]; assert len(H) == 577
    n, e = 23, 1
    A = contraction(H, e); C = section(H,e,0) & section(H,e,1)
    ks = components(A,C,n); assert len(ks) == 3
    exclusive = [x for x in A-C if corner(A,x,n)]
    hyper = [sorted({i for i,K in enumerate(ks)
                     if any(x^(1<<f) in K for f in range(n))})
             for x in sorted(C) if corner(C,x,n)]
    assert not exclusive and sorted(hyper) == [[0,1],[0,1],[0,1],[1,2]]
    rows = []
    for labels in [(0,0,0),(0,0,1),(0,1,0),(0,1,1)]:
        X = lift(A,C,ks,labels,n)
        cs = sorted(x for x in X if corner(X,x,n+1))
        if labels == (0,1,0):
            assert not cs; order = None
        else:
            order = greedy_order(X,n+1); assert order is not None
            replay(X,order,n+1)
        rows.append(dict(labels=labels, corners=cs, peeling=order))
    codes = [0,2,3]; Z = cube_code(A,C,n,2,codes)
    assert len(Z) == 1601
    predicted = set()
    # Verify the stated uniform corner formula against the definition.
    for x in C:
        adjacent = {i for i,K in enumerate(ks) if any(x^(1<<f) in K for f in range(n))}
        for y in range(4):
            if corner(C,x,n) and all(codes[i] == y for i in adjacent):
                predicted.add(x | (y<<n))
    for K,q in zip(ks,codes):
        for x in K:
            for y in range(4):
                if (y^q).bit_count() == 1 and corner(A,x,n):
                    predicted.add(x | (y<<n))
    actual = {x for x in Z if corner(Z,x,n+2)}
    assert actual == predicted == set()
    projections = []
    for f, row in [(n,rows[3]),(n+1,rows[1])]:
        M = section(Z,f,0) & section(Z,f,1)
        X = lift(A,C,ks,row['labels'],n)
        assert M == frozenset(x^(1<<n) for x in X)
        order = [x^(1<<n) for x in row['peeling']]
        replay(M,order,n+1)
        projections.append(dict(coordinate=f, family=sorted(M), peeling=order))
    bad = section(Z,n,0)
    assert len(bad) == 865 and not any(corner(bad,x,n+1) for x in bad)
    # Whole-word compatibility constructs an actual order, including fresh
    # variables, full core fibers, exclusive batches, and the final fiber.
    coherent = []
    for row in (rows[1],rows[3]):
        codes2 = [0 if t == 0 else 5 for t in row['labels']]
        result = expand_binary_trace(A,C,row['labels'],row['peeling'],n,3,codes2)
        assert result is not None
        coherent.append(result)
        assert expand_binary_trace(A,C,row['labels'],row['peeling'],n,2,[0,2,3]) is None
    starA,starC = frozenset([0,1,2,4]),frozenset([0])
    starlabels = [0,0,0]; starX = lift(starA,starC,components(starA,starC,3),starlabels,3)
    # Remove all three exclusive leaves before breaking the core fiber.
    starorder = [1,2,4,0]; replay(starX,starorder,4)
    result = expand_binary_trace(starA,starC,starlabels,starorder,3,3,[0,3,5])
    assert result is not None and result['partition'] == [[0],[1],[2]]
    coherent.append(dict(input_A=sorted(starA),input_C=sorted(starC),**result))
    out = dict(status='PASS: bitwise positivity fails; square is not pc-minimal',
        source=str(source.relative_to(ROOT)), source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
        source_leaf=sorted(H), marked_coordinate=e, dimension=n,
        A=sorted(A), C=sorted(C), components=[sorted(K) for K in ks],
        core_corner_hyperedges=hyper, binary_cases_up_to_complement=rows,
        square_codes=codes, square_family=sorted(Z), positive_strong_projections=projections,
        coherent_trace_controls=coherent,
        cornerless_proper_halfspace=dict(coordinate=n,side=0,family=sorted(bad)),
        scope='Ample inputs inherited from the archived construction and its corner deletions; no terminal-seed or minimality theorem.')
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print('PASS: six positive binary labels; cornerless square 1601; negative halfspace 865; three full-word trace orders replayed')

if __name__ == '__main__':
    main()
