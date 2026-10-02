"""Replay the inputs and finite controls for unbounded protected-leaf depth.

Uniform theorem is symbolic. Controls m=1,2 check a fixed-envelope pruning
order, actual tip-fibre deletions, and independent canonical-core emptiness.
No exhaustive corner-choice tree or claim of exact depth is made.
"""
from pathlib import Path
import json, hashlib
from first_deletion_constraints import canonical_core
from component_label_partitions import corner, replay

HERE = Path(__file__).resolve().parent


def row_masks(m):
    bits = [1] + [1 << (24+j) for j in range(m-1)]
    return [sum(b for j,b in enumerate(bits) if t >> j & 1) for t in range(1 << m)]


def block(X, m):
    masks = row_masks(m)
    zero = {x for x in X if not x & 1}
    one = {x ^ 1 for x in X if x & 1}
    return {x | t for i,t in enumerate(masks)
            for x in (zero if i == 0 else one if i == len(masks)-1 else zero | one)}


def complete_order(X, prefix):
    return prefix + sorted(X-set(prefix))


def lifted_order(X, base_order, m):
    left = set(X)
    masks = row_masks(m)
    result = []
    for x in base_order:
        y, side = x & ~1, bool(x & 1)
        pole = masks[-1] if side else 0
        if x ^ 1 in left:
            result.append(y | pole)
        else:
            omitted = 0 if side else masks[-1]
            row = sorted(set(masks)-{omitted}, key=lambda t: (-(t ^ pole).bit_count(), t))
            result.extend(y | t for t in row)
        left.remove(x)
    return result


def weak_check(X, order, n):
    """Later incident directions must span a cube in the FIXED original X."""
    assert len(order) == len(X) and set(order) == X
    remaining = set(X)
    for x in order:
        J = sum(1 << e for e in range(n) if x ^ (1 << e) in remaining)
        sub = J
        while sub:
            assert x ^ sub in X, (x, sub)
            sub = (sub-1) & J
        remaining.remove(x)


def main():
    source = HERE/'root_cube_completion.json'
    delay = HERE/'protected_core_delay.json'
    d, cert = json.loads(source.read_text()), json.loads(delay.read_text())
    A = set(d['cube_completion']['factor'])
    ao = complete_order(A, d['cube_completion']['factor_peeling'])
    replay(A, ao[:-1], 12)
    assert {x for x in A if corner(A, x, 12)} == {0}
    B = {x << 12 for x in A}; bo = [x << 12 for x in ao]
    D = sum(1 << e for e in d['cube_completion']['root_directions'])
    cube = {0};sub = D
    while sub:
        cube.add(sub);sub = (sub-1) & D
    S = A | B;H = S | cube;K = H | {32}
    assert len(H)==626 and H==set(cert['base_family']) and len(K)==627
    assert 1 not in A and not D & 1
    assert {e for e in range(12) if 1 << e in A} == {1,2,8}
    assert sum(1<<e for e in range(24) if 32 ^ (1<<e) in H)==38
    assert 33 not in H and 32 not in H
    assert len({x & ~1 for x in H})==559
    assert len({x for x in H if not x & 1 and x ^ 1 in H})==67
    controls=[]
    for m in (1,2):
        n=23+m;M=row_masks(m);P=set(M)-{M[-1]}
        po=sorted(P,key=lambda t:(-t.bit_count(),t))
        HM,KM=block(H,m),block(K,m)
        WA=block(A,m);BP={x|t for x in B for t in P}
        assert HM==WA|BP|{x|t for x in cube for t in P}
        waorder=lifted_order(A,ao,m);replay(WA,waorder[:-1],n)
        bporder=[x|t for x in bo for t in po];replay(BP,bporder[:-1],n)
        # First remove all mixed old words and the old root, fibre by fibre.
        interface=[x|t for x in sorted((cube-S)|{0}) for t in po]
        used=set(interface)
        weak=interface+[x for x in waorder if x not in used]
        used=set(weak);weak += [x for x in bporder if x not in used]
        weak_check(HM,weak,n)
        assert not canonical_core(HM,n)['core']
        fibre=[32|t for t in po]
        assert set(fibre)==KM-HM and len(fibre)==(1<<m)-1
        remaining=set(KM)
        for x in fibre:
            assert corner(remaining,x,n)
            remaining.remove(x)
        assert remaining==HM
        assert not canonical_core(KM,n)['core']
        assert len(HM)==559*(1<<m)-492 and len(KM)==560*(1<<m)-493
        controls.append(dict(block_size=m,base_size=len(HM),output_size=len(KM),
                             base_fixed_envelope_order=weak,tip_fibre_deletion_order=fibre,
                             protected_leaf_depth_lower_bound=1<<m))
    output=dict(status='PASS',controls=controls,
                source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (source,delay,Path(__file__))},
                uniform_claim='For all m>=1, forbidden K_m has order560*2^m-493, dimension m+23, VC dimension m+5 and protected-leaf depth at least2^m.',
                scope='Uniform symbolic proof; finite replay at m=1,2. These outputs are not reduction-terminal for m>1 and compress to the same627 seed. No exact-depth, primitive-coverage or terminal-depth claim.')
    Path(__file__).with_suffix('.json').write_text(json.dumps(output,indent=2)+'\n')
    print('PASS: fixed-envelope pruning, corner-fibre deletion, independent empty cores; m=1,2.',flush=True)


if __name__=='__main__': main()
