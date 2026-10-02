"""Finite replay for the uniform reduced, reduction-terminal depth theorem.

The proof, not a census, establishes every m. Check m=1,2: factor orders,
lifted frozen-root witnesses, all possible old first-move certificates,
actual tip-fibre deletion, empty core, degrees and column-separating pairs.
"""
from pathlib import Path
import json,hashlib
from component_label_partitions import corner,replay
from first_deletion_constraints import canonical_core
from protected_block_delay import complete_order,weak_check

HERE=Path(__file__).resolve().parent
EOLD=1<<8


def masks(m):
    bits=[EOLD]+[1<<(24+j) for j in range(m-1)]
    return [sum(b for j,b in enumerate(bits) if t>>j&1) for t in range(1<<m)]


def block(X,m):
    ms=masks(m);zero={x for x in X if not x&EOLD};one={x^EOLD for x in X if x&EOLD}
    return {x|t for i,t in enumerate(ms) for x in
            (zero if i==0 else one if i==len(ms)-1 else zero|one)}


def lift_order(X,order,m):
    left=set(X);out=[];ms=masks(m)
    for x in order:
        y=x&~EOLD;s=bool(x&EOLD);pole=ms[-1] if s else 0
        if x^EOLD in left:out.append(y|pole)
        else:
            row=set(ms)-{0 if s else ms[-1]}
            out += [y|t for t in sorted(row,key=lambda t:(-(t^pole).bit_count(),t))]
        left.remove(x)
    return out


def cube_at_zero(D):
    out={0};s=D
    while s:out.add(s);s=(s-1)&D
    return out


def check_witness(S,X,row,n):
    x,J=row['word'],row['directions']
    assert J.bit_count()>=2 and row['missing']==x^J and x^J not in X
    assert all(x^(1<<e) in S for e in range(n) if J>>e&1)


def main():
    source=HERE/'root_cube_completion.json';delay=HERE/'protected_core_delay.json'
    d,c=json.loads(source.read_text()),json.loads(delay.read_text())
    A=set(d['cube_completion']['factor']);ao=complete_order(A,d['cube_completion']['factor_peeling'])
    replay(A,ao[:-1],12)
    B={x<<12 for x in A};bo=[x<<12 for x in ao]
    assert {i for i in range(12) if 1<<i in A}=={1,2,8}
    witnesses={r['word']:r for r in c['rooted_factor_growth']['witnesses']}
    assert set(witnesses)==A-{0}
    for r in witnesses.values():check_witness(A,A|{32},r,12)
    assert len({x&~EOLD for x in A})==224
    assert len({x for x in A if not x&EOLD and x^EOLD in A})==65
    records=[]
    for m in (1,2):
        n=m+23;ms=masks(m);eb=ms[-1];P=set(ms)-{eb}
        po=sorted(P,key=lambda t:(-t.bit_count(),t))
        Am,Ap=block(A,m),block(A|{32},m)
        amorder=lift_order(A,ao,m);replay(Am,amorder[:-1],n)
        assert {x for x in Am if corner(Am,x,n)}=={0}
        DA=(1<<1)|(1<<2)|eb;DB=sum(1<<i for i in (13,14,20));D=DA|DB
        Q=cube_at_zero(D);S=Am|B;H=S|Q;tips={32|t for t in P};K=H|tips
        assert Ap==Am|tips and len(H)==252*(1<<m)+122 and len(K)==253*(1<<m)+121
        frozen=[]
        for v in sorted(Am-{0}):
            t=v&eb;y=v&~eb
            if t==0:x=y
            elif t==eb:x=y|EOLD
            elif y==0:x=EOLD
            else:x=next(z for z in (y,y|EOLD) if z in A)
            assert x in witnesses
            J=witnesses[x]['directions']
            target_side=bool(x&EOLD)^bool(J&EOLD)
            target=eb if target_side else 0
            L=(J&~EOLD)|(t^target)
            r=dict(word=v,directions=L,missing=v^L)
            check_witness(Am,Ap,r,n);check_witness(S,K,r,n);frozen.append(r)
        for x in sorted(B-{0}):
            r=witnesses[x>>12]
            row=dict(word=x,directions=r['directions']<<12,missing=r['missing']<<12)
            check_witness(S,K,row,n);frozen.append(row)
        mixed=Q-S
        assert {x for x in K if corner(K,x,n)}==mixed|{32|t for t in P if corner(P,t,n)}
        # Every mixed first deletion protects S; all old nonroot rows unchanged.
        for x in mixed:
            check_witness(S,K-{x},dict(word=0,directions=x,missing=x),n)
        # If the root tip is absent, root deletion protects the entire old base.
        root_child=H-{0};root_rows=[]
        frozen_by_word={r['word']:r for r in frozen}
        for x in sorted(root_child):
            if x in mixed:J=x
            elif x.bit_count()==1:
                factor,own,other=(Am,DA,DB) if x in Am else (B,DB,DA)
                g=next(1<<i for i in range(n) if not own>>i&1 and x^(1<<i) in factor)
                f=other&-other;J=g|f
            else:J=frozen_by_word[x]['directions']
            row=dict(word=x,directions=J,missing=x^J)
            check_witness(root_child,(K-{32})-{0},row,n);root_rows.append(row)
        interface=sorted(mixed|{0});used=set(interface)
        weak=interface+[x for x in amorder if x not in used];used=set(weak)
        weak += [x for x in bo if x not in used]
        weak_check(H,weak,n);assert not canonical_core(H,n)['core']
        remaining=set(K);tip_order=[32|t for t in po]
        for x in tip_order:
            assert corner(remaining,x,n);remaining.remove(x)
        assert remaining==H and not canonical_core(K,n)['core']
        assert min(sum(x^(1<<i) in K for i in range(n)) for x in K)>=3
        left=set(range(12))|set(range(24,n));right=set(range(12,24))
        offleft={i for i in left if not DA>>i&1};offright={i for i in right if not DB>>i&1}
        pairs={(f,q) for f in left for q in offright}|{(f,q) for f in offleft for q in right}
        for f,q in pairs:
            mask=(1<<f)|(1<<q)
            assert {x&mask for x in K}=={0,1<<f,1<<q}
        patterns={i:frozenset(pair for pair in pairs if i in pair) for i in range(n)}
        assert len(set(patterns.values()))==n
        records.append(dict(block_size=m,order=len(K),base_order=len(H),
            frozen_witnesses=frozen,root_deletion_witnesses=root_rows,
            fixed_envelope_pruning_order=weak,tip_fibre_order=tip_order,
            column_separating_pairs=sorted(pairs),protected_depth_lower_bound=1<<m))
    report=dict(status='PASS',controls=records,
        source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (source,delay,Path(__file__))},
        scope='Uniform symbolic theorem: reduced and reduction-terminal forbidden classes of order253*2^m+121 and depth>=2^m. Finite replay m=1,2; no exact-depth or exhaustive-generator claim.')
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS: terminal construction controls m=1,2; frozen witnesses, first moves, empty cores, fibre orders and reduced-column pairs.',flush=True)


if __name__=='__main__':main()
