"""Uncoloured equality normal form and a certificate-free minor selection rule.

Bounded validation: archived positive-pair controls, all 48 ample sign splits
above P5, and two positive-overlap presentations. No output recognition is used
by structural_select; rper comparisons below are independent checks only.
"""
from pathlib import Path
import sys,json,hashlib
from itertools import product
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT),str(Path(__file__).parent)]
from component_label_recursion import components,lift,circuit,evaluate
from explore_rper_cyclic_families import section,contraction,normalize,rper

def partition(k,groups):
    parent=list(range(k))
    def find(i):
        while parent[i]!=i:
            parent[i]=parent[parent[i]];i=parent[i]
        return i
    for group in groups:
        group=list(group)
        for i in group[1:]:parent[find(i)]=find(group[0])
    return sorted(sorted(i for i in range(k) if find(i)==r) for r in {find(i) for i in range(k)})

def common(A,C,remaining):
    return [(e,b) for e in range(remaining.bit_length()) if remaining>>e&1 for b in (0,1)
            if all(x^(1<<e) in A for x in A if x>>e&1==b)
            and all(x^(1<<e) in C for x in C if x>>e&1==b)]

def normal_partition(A,C,n,reverse=False):
    """Strip any common coordinate, irrespective of colours; unite leaf indices."""
    blocks=components(A,C,n);owner={x:i for i,K in enumerate(blocks) for x in K}
    leaves=[]
    def split(a,c,remaining):
        cuts=common(a,c,remaining)
        if not cuts:
            leaves.append(sorted({owner[x] for x in a-c}));return
        e,_=cuts[-1] if reverse else cuts[0];bit=1<<e
        for b in (0,1):split(frozenset(x for x in a if x>>e&1==b),
                             frozenset(x for x in c if x>>e&1==b),remaining^bit)
    split(A,C,(1<<n)-1)
    return partition(len(blocks),leaves)

def graph_partition(A,C,n):
    graph=circuit(A,C,n)
    return partition(len(graph['components']),[q['mono'] for q in graph['nodes'] if not q['cuts']])

def monochrome(t,p):return all(len({t[i] for i in B})<=1 for B in p)

def proper_partition(A,C,n):
    blocks=components(A,C,n);groups=[];records=[]
    for e in range(n):
        for kind in (0,1,2):
            mu=lambda F:contraction(F,e) if kind==2 else section(F,e,kind)
            a,c=mu(A),mu(C);ks=components(a,c,n-1)
            images=[mu(K)-c for K in blocks]
            parents=[]
            for K in ks:
                ids=[i for i,M in enumerate(images) if M&K]
                assert len(ids)==1;parents.append(ids[0])
            p=normal_partition(a,c,n-1)
            q=normal_partition(a,c,n-1,True);assert p==q
            groups += [[parents[j] for j in B] for B in p]
            records.append(dict(coordinate=e,operation=kind,partition=p,component_parents=parents))
    return partition(len(blocks),groups),records

def peripheral(H,n):
    return [(e,b) for e in range(n) if len({x>>e&1 for x in H})==2 for b in (0,1)
            if all(x^(1<<e) in H for x in H if x>>e&1==b)]

def orientation_forces(A,C,n,wing,e,side):
    owner={x:i for i,K in enumerate(components(A,C,n)) for x in K};forced=[]
    for x in A:
        if x>>e&1!=side:continue
        y=x^(1<<e)
        if x in C:
            if y not in A:return None
            if y not in C:forced.append((owner[y],wing))
        elif y not in A:forced.append((owner[x],1-wing))
    return forced

def check_forces(A,C,n,t):
    ks=components(A,C,n)
    for wing in (0,1):
        B=frozenset(set(C)|set().union(*(K for i,K in enumerate(ks) if t[i]==wing)))
        for e in range(n):
            for side in (0,1):
                forced=orientation_forces(A,C,n,wing,e,side)
                predicted=forced is not None and all(t[i]==b for i,b in forced)
                assert predicted==all(x^(1<<e) in B for x in B if x>>e&1==side)

def structural_select(A,C,n,t,p,positive_overlap):
    """Caller generates positive A and positive C, or already-generated negative C."""
    blocks=components(A,C,n)
    if len(set(t))!=2 or not monochrome(t,p):return False
    if positive_overlap and common(A,C,(1<<n)-1):return False
    wings=[frozenset(set(C)|set().union(*(K for i,K in enumerate(blocks) if t[i]==b))) for b in (0,1)]
    return all(peripheral(B,n) for B in wings)

def main():
    source=Path(__file__).with_name('component_label_recursion.json');old=json.loads(source.read_text())
    comparisons=0;controls=[]
    P5=frozenset([0,1,2,3,4,6,7,8,12,14,16,17,19,24,25,28])
    pairs=[(contraction(P5,e),section(P5,e,0)&section(P5,e,1),4) for e in range(5)]
    pairs.append((frozenset([0,1,2,4]),frozenset([0]),3))
    for A,C,n in pairs:
        p=normal_partition(A,C,n);assert p==normal_partition(A,C,n,True)==graph_partition(A,C,n)
        ks=components(A,C,n);g=circuit(A,C,n);accepted=[]
        for t in product((0,1),repeat=len(ks)):
            ok=monochrome(t,p);assert ok==evaluate(g,t)==rper(*normalize(lift(A,C,ks,t,n),n+1))
            comparisons+=1
            if ok:accepted.append(list(t))
        controls.append(dict(partition=p,accepted=accepted))
    # All ample twists retained by the previous experiment, including nonminimal outputs.
    amp=json.loads((ROOT/'papers/minor-atlas/research/coordinate_block_amplification_verification.json').read_text())
    rows=amp['source_clauses'];relevant=[i for i,(s,t) in enumerate(rows) if s&1];twists=[]
    for row in old['cases']:
        twist=row['twist']
        clauses=[(s|32,t|(((t&1)^((twist>>relevant.index(i))&1))<<5)) if s&1 else (s,t) for i,(s,t) in enumerate(rows)]
        X=frozenset(x for x in range(64) if all(x&s!=t for s,t in clauses))
        A=contraction(X,5);C=section(X,5,0)&section(X,5,1)
        if not rper(*normalize(A,5)):
            twists.append(dict(twist=twist,ineligible_negative_base=True));continue
        p,records=proper_partition(A,C,5);ks=components(A,C,5)
        t=[0 if K<=section(X,5,0) else 1 for K in ks]
        # Challenge the single-subspace law on every label of every old-minor pair.
        for e in range(5):
            for kind in (0,1,2):
                mu=lambda F:contraction(F,e) if kind==2 else section(F,e,kind)
                a,c=mu(A),mu(C);parts=normal_partition(a,c,4);g=circuit(a,c,4);blocks=components(a,c,4)
                assert parts==graph_partition(a,c,4)
                for labels in product((0,1),repeat=len(blocks)):
                    assert monochrome(labels,parts)==evaluate(g,labels)==rper(*normalize(lift(a,c,blocks,labels,4),5))
                    comparisons+=1
        check_forces(A,C,5,t)
        selected=structural_select(A,C,5,t,p,False)
        assert selected==row['forbidden']
        twists.append(dict(twist=twist,partition=p,labels=t,selected=selected))
        comparisons+=1
    positive=json.loads(Path(__file__).with_name('positive_overlap_generator.json').read_text());examples=[]
    for r in positive['examples']:
        A=frozenset(r['positive_inputs']['A']['original_family']);C=frozenset(r['positive_inputs']['C']['original_family']);n=r['dimension']-1
        check_forces(A,C,n,r['labels'])
        p,records=proper_partition(A,C,n);assert structural_select(A,C,n,r['labels'],p,True)
        wings=[frozenset(r['positive_inputs'][key]['original_family']) for key in ('B0','B1')]
        examples.append(dict(order=r['order'],proper_partition=p,labels=r['labels'],wing_orientations=[peripheral(B,n) for B in wings],old_minor_partitions=records))
    out=dict(status='PASS',controls=controls,twists=twists,positive_overlap_examples=examples,comparisons=comparisons,
             verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,ROOT/'papers/minor-atlas/research/coordinate_block_amplification_verification.json',Path(__file__).with_name('positive_overlap_generator.json')]})
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print('PASS: choice-independent equality partitions and structural selection; 48 ample twists and 2 positive-overlap controls;',comparisons,'comparisons')
if __name__=='__main__':main()
