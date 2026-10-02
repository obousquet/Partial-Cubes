"""Check the cube-code construction, positivity partition, and unsigned module.

This reuses a formerly exploratory 166-word square example from the fixed
80-word seed, rather than searching a new census. The uniform theorem is in
component_label_recursion.tex; these exact controls verify its encoding.
"""
from pathlib import Path
import hashlib,itertools,json,sys
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT),str(Path(__file__).parent)]
from component_label_recursion import components,shadow,certificate
from canonical_label_partition import normal_partition,monochrome,proper_partition,partition
from positive_overlap_generator import check_tree
from explore_rper_cyclic_families import section,contraction,normalize,rper,peripheral_coords

def cube_code(A,C,n,d,codes):
    ks=components(A,C,n);assert len(ks)==len(codes)
    return frozenset({x|(b<<n) for x in C for b in range(1<<d)}|
                     {x|(b<<n) for K,q in zip(ks,codes) for x in K for b in range(1<<d) if b!=q})

def predicted_shadow(A,C,n,d):
    SA,SC=shadow(A,n),shadow(C,n)
    return {s|(t<<n) for s in SA for t in range((1<<d)-1)}|{s|(((1<<d)-1)<<n) for s in SC}

def main():
    source=Path(__file__).with_name('peripheral_exchange_audit.json');r=json.loads(source.read_text());H=frozenset(r['family'])
    A=contraction(H,6);C=section(H,6,0)&section(H,6,1);ks=components(A,C,6);codes=[0,1,2]
    X=cube_code(A,C,6,2,codes);assert len(X)==166
    SX=shadow(X,8);assert SX==predicted_shadow(A,C,6,2) and len(SX)==len(X)
    assert not peripheral_coords(X,8)
    trees=[]
    for e in range(8):
        for op in range(3):
            M=contraction(X,e) if op==2 else section(X,e,op);t=certificate(M,7);check_tree(t)
            trees.append(dict(coordinate=e,operation=op,tree=t))
    nf=[s for s in range(256) if s not in SX and all(s^(1<<i) in SX for i in range(8) if s>>i&1)]
    clauses=[(s,next(t for t in range(256) if not t&~s and t not in {x&s for x in X})) for s in nf]
    columns=[tuple(0 if not s>>i&1 else 2*((t>>i)&1)-1 for s,t in clauses) for i in range(8)]
    assert len({min(c,tuple(-x for x in c)) for c in columns})==8
    assert all(bool(s&64)==bool(s&128) for s,t in clauses)
    positive=[e for e in range(8) if rper(*normalize(section(X,e,0)&section(X,e,1),7))]
    assert positive==[1,2,3,4]
    # The four positive old overlaps are forced for every cube dimension/code.
    prop,_=proper_partition(A,C,6);assert prop==[[0],[1],[2]]
    forced=[]
    for f in (1,2,3,4):
        a=section(A,f,0)&section(A,f,1);c=section(C,f,0)&section(C,f,1)
        children=components(a,c,5);images=[contraction(K,f)&(a-c) for K in ks]
        parents=[]
        for K in children:
            ids=[i for i,M in enumerate(images) if M&K];assert len(ids)==1;parents.append(ids[0])
        p=normal_partition(a,c,5)
        pulled=partition(len(ks),[[parents[i] for i in B] for B in p]);assert pulled==[[0],[1],[2]]
        ta,tc=certificate(a,5),certificate(c,5);check_tree(ta);check_tree(tc)
        forced.append(dict(coordinate=f,partition=p,component_parents=parents,pulled_partition=pulled,
                           positive_contraction_tree=ta,positive_core_tree=tc))
    # Positive-input controls exercise all square codes, including nonpositive outputs.
    P=frozenset([0,1,2,3,4,6,7,8,12,14,16,17,19,24,25,28])
    checks=0
    for e in range(5):
        a=contraction(P,e);c=section(P,e,0)&section(P,e,1);p=normal_partition(a,c,4);k=len(components(a,c,4))
        for qs in itertools.product(range(4),repeat=k):
            Y=cube_code(a,c,4,2,qs)
            assert shadow(Y,6)==predicted_shadow(a,c,4,2)
            assert monochrome(qs,p)==rper(*normalize(Y,6))
            checks+=1
    out=dict(status='PASS; global reduced-core conjecture remains open',family=sorted(X),dimension=8,
             contraction_input=sorted(A),strong_input=sorted(C),components=[sorted(K) for K in ks],corner_codes=codes,
             clauses=clauses,proper_minor_trees=trees,positive_projection_coordinates=positive,
             equal_unsigned_columns=[6,7],all_signed_columns_distinct=True,positive_pair_code_checks=checks,proper_partition=prop,uniformly_forced_positive_overlaps=forced,
             source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print('PASS: ample reduced 166-word cube-code obstruction; 24 minor trees; four positive overlaps;',checks,'square-code partition checks; four old overlaps uniformly forced positive')
if __name__=='__main__':main()
