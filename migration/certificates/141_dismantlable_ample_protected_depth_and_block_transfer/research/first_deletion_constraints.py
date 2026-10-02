"""Uncoloured first-deletion clauses and protected-cube terminal seeds.

Exact tests: all nested ample pairs in Q_2 and Q_3; all labelings of
267/292/positive-control coordinate pairs; deeper negative 580 control.
Protected witnesses prove the historical 14-tip batch without SAT or subset
enumeration. A failed greedy order is UNKNOWN, never a negative certificate.
"""
from pathlib import Path
import hashlib,json,sys
ROOT=Path(__file__).resolve().parents[3]
HERE=Path(__file__).parent
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(HERE))
from component_label_partitions import corner,greedy_order,replay
from component_label_recursion import components,lift,shadow
from explore_rper_cyclic_families import section,contraction

def bc(S,t):return len({t[i] for i in S})==2

def child_constraints(A,C,ids,n,drop_fresh=False):
    outside=[x for x in sorted(A-C) if corner(A,x,n)]
    if outside:return dict(possible=False,exclusive_corners=outside,edges=[])
    blocks=components(A,C,n);mapping={}
    for block in blocks:
        old={ids[x] for x in block if x in ids}
        # Under a core-deletion guard, merged old labels are equal.
        var=min(old) if old else None
        assert var is not None or drop_fresh
        for x in block:mapping[x]=var
    edges=[]
    for x in sorted(C):
        if not corner(C,x,n):continue
        variables={mapping[x^(1<<e)] for e in range(n) if x^(1<<e) in A-C}
        # A free new component is universally quantified. Bichromaticity
        # for both its colors is equivalent to bichromaticity without it.
        variables.discard(None)
        edges.append(sorted(variables))
    return dict(possible=True,exclusive_corners=[],edges=edges)

def generate(A,C,n):
    blocks=components(A,C,n);ids={x:i for i,b in enumerate(blocks) for x in b}
    exclusive=[];core=[]
    for x in sorted(A-C):
        if corner(A,x,n):
            exclusive.append(dict(word=x,child=child_constraints(A-{x},C,ids,n)))
    for x in sorted(C):
        if corner(C,x,n):
            N=sorted({ids[x^(1<<e)] for e in range(n) if x^(1<<e) in A-C})
            core.append(dict(word=x,neighbors=N,
                             child=child_constraints(A,C-{x},ids,n,drop_fresh=not N)))
    return dict(components=[sorted(b) for b in blocks],exclusive_moves=exclusive,core_moves=core,
                pair_states=1+len(exclusive)+len(core))

def child_accepts(child,t):return child['possible'] and all(bc(S,t) for S in child['edges'])

def evaluate(model,t):
    return all(child_accepts(m['child'],t) for m in model['exclusive_moves']) and all(
        bc(m['neighbors'],t) or child_accepts(m['child'],t) for m in model['core_moves'])

def cnf(model):
    """Two signed clauses express BC; distribute each core implication."""
    clauses=[]
    def bc_clauses(S):return [list(i+1 for i in S),list(-i-1 for i in S)]
    for m in model['exclusive_moves']:
        child=m['child']
        if not child['possible']:clauses.append([])
        else:
            for S in child['edges']:clauses+=bc_clauses(S)
    for m in model['core_moves']:
        guard=bc_clauses(m['neighbors']);child=m['child']
        if not child['possible']:clauses+=guard
        else:
            for S in child['edges']:
                for a in guard:
                    for b in bc_clauses(S):clauses.append(sorted(set(a+b)))
    # Omit tautologies and duplicate clauses; an empty clause is false.
    return sorted({tuple(c) for c in clauses if not any(-x in c for x in c)})

def eval_cnf(clauses,t):return all(any(t[abs(x)-1]==(x>0) for x in c) for c in clauses)

def direct(X,n):
    assert len(X)>1
    for x in X:
        if corner(X,x,n) and any(corner(X-{x},y,n) for y in X-{x}):return False
    return True

def test_pair(A,C,n):
    model=generate(A,C,n);blocks=[frozenset(b) for b in model['components']];clauses=cnf(model);tests=[]
    for mask in range(1<<len(blocks)):
        t=[mask>>i&1 for i in range(len(blocks))];X=lift(A,C,blocks,t,n)
        if len(X)<=1:continue
        answer=evaluate(model,t)
        assert answer==direct(X,n+1)==eval_cnf(clauses,t)
        tests.append(dict(labels=t,accepted=answer,corners=[x for x in sorted(X) if corner(X,x,n+1)]))
    return model,clauses,tests

def canonical_core(H,n,max_rank=None):
    S=set(H);rounds=[]
    while S:
        deleted=[];witnesses=[]
        for x in sorted(S):
            ds=[1<<e for e in range(n) if x^(1<<e) in S]
            masks=[0]
            for d in ds:masks += [j|d for j in masks]
            J=next((j for j in sorted(masks,key=lambda j:(j.bit_count(),j)) if (max_rank is None or j.bit_count()<=max_rank) and x^j not in H),None)
            if J is None:deleted.append(x)
            else:witnesses.append(dict(word=x,directions=J,missing=x^J))
        if not deleted:return dict(core=sorted(S),rounds=rounds,witnesses=witnesses)
        rounds.append(deleted);S.difference_update(deleted)
    return dict(core=[],rounds=rounds,witnesses=[])

def retained_certificate(H,n):
    tips={x for x in H if corner(H,x,n)}
    S=H-tips;assert len(S)>=2
    rows=[]
    for x in sorted(S):
        ds=[1<<e for e in range(n) if x^(1<<e) in S]
        masks=[0]
        for d in ds:masks += [j|d for j in masks]
        J=next((j for j in sorted(masks,key=lambda j:(j.bit_count(),j)) if x^j not in H),None)
        assert J is not None,(x,'no protected witness')
        assert all(x^(1<<e) in S for e in range(n) if J>>e&1)
        rows.append(dict(word=x,directions=J,missing=x^J))
    return dict(core=sorted(S),core_order=len(S),corners=sorted(tips),witnesses=rows)

def main():
    small=[]
    for n in (2,3):
        families=[]
        for mask in range(1,1<<(1<<n)):
            A=frozenset(x for x in range(1<<n) if mask>>x&1)
            if len(shadow(A,n))==len(A):families.append(A)
        pairs=tests=0
        for A in families:
            for C in [frozenset()]+families:
                if not C<=A or len(A)+len(C)<=1:continue
                model,clauses,checks=test_pair(A,C,n);pairs+=1;tests+=len(checks)
        small.append(dict(dimension=n,ample_families=len(families),nested_pairs=pairs,label_checks=tests))
        print('small',small[-1],flush=True)
    paths=['papers/dismantlable-ample/evidence/order267_obstruction/record.json',
           'papers/dismantlable-ample/research/root_cube_completion.json',
           'papers/dismantlable-ample/research/component_label_partitions.json',
           'forbidden_pc_minor_campaign/damp_fixed_305_completion_verification.json']
    src=[json.loads((ROOT/p).read_text()) for p in paths]
    K=frozenset(src[0]['family']);H=frozenset(src[1]['cornered_terminal_control']['family'])
    P=frozenset(src[1]['cornered_terminal_control']['positive_control'])
    controls=[]
    for name,X in [('267',K),('292',H),('292_positive',P)]:
        rows=[]
        for e in range(12):
            A=contraction(X,e);C=section(X,e,0)&section(X,e,1)
            model,clauses,checks=test_pair(A,C,11)
            assert len(model['components'])==2
            assert [r['labels'] for r in checks if r['accepted']]==([[1,0],[0,1]] if name!='292_positive' else [])
            rows.append(dict(coordinate=e,model=model,cnf=clauses,label_checks=checks))
        controls.append(dict(name=name,order=len(X),coordinate_models=rows))
        print('large',name,'all 12 pairs pass',flush=True)
    # Replay the existing 292 proper-minor and strong-projection certificates.
    for row in src[1]['cornered_terminal_control']['minor_checks']:
        for c in row['certificates']:replay(frozenset(c['family']),c['order'],11)
    # Supply all 36 direct 267 proper-minor orders; strong orders already exist.
    minors=[]
    for e in range(12):
        row=dict(coordinate=e,certificates=[])
        for name,M in [('zero',section(K,e,0)),('one',section(K,e,1)),('contract',contraction(K,e))]:
            p=greedy_order(M,11);assert p is not None
            row['certificates'].append(dict(operation=name,family=sorted(M),order=p))
        replay(frozenset(src[2]['seed_checks'][e]['overlap']),src[2]['seed_checks'][e]['overlap_order'],11)
        minors.append(row)
    # A deeper connected-support seed, produced by one ample addition to H292.
    J=H|{553};jshadow=shadow(J,12)
    assert len(jshadow)==293 and jshadow-shadow(H,12)=={585}
    assert not direct(J,12)
    from root_face_recovery import support_components
    assert all(len(support_components(J,12,r)[2])==1 for r in J)
    retained=[]
    for name,X in [('267',K),('292',H),('293',J)]:
        cert=retained_certificate(X,12)
        core=frozenset(cert['core']);cert['core_shadow_count']=len(shadow(core,12))
        retained.append(dict(name=name,certificate=cert))
    assert retained[-1]['certificate']['core_order']==289
    subcore=retained[-1]['certificate']
    assert set(subcore['core'])<=H
    assert all(row['missing'] not in H and row['missing']!=553 for row in subcore['witnesses'])
    from root_cube_completion import certificate
    states=certificate(J,12);assert len(states)==9
    print('293: ample addition acquires support 585; connected at all roots; protected core 289; nine negative states',flush=True)
    canonical=[]
    for name,X,n in [('267',K,12),('292',H,12),('293',J,12),('580',frozenset(src[1]['cube_completion']['family']),24),('626',frozenset(src[1]['full_root_cube_control']['family']),24)]:
        row=canonical_core(X,n)
        expected={'267':267,'292':290,'293':289,'580':577,'626':0}[name]
        assert len(row['core'])==expected
        canonical.append(dict(name=name,certificate=row))
        print('canonical protected core',name,len(row['core']),'pruning rounds',len(row['rounds']),flush=True)
    # A cube-witness proof for the historical 14-tip batch, independent of SAT.
    G=H-{186};T=frozenset(src[3]['tips']);batch=G|T;gshadow=shadow(G,12)
    assert len(G)==len(gshadow)==291 and not any(corner(G,x,12) for x in G)
    assert len(batch)==len(shadow(batch,12))==305
    acquired={}
    for a in sorted(T):
        sh=shadow(G|{a},12);I=int(src[3]['acquired_supports'][str(a)])
        assert len(sh)==292 and sh-gshadow=={I}
        acquired[a]=I
    assert len(set(acquired.values()))==14
    rank_cores=[]
    for cap,expected in [(2,0),(3,248),(4,281),(None,281)]:
        cert=canonical_core(batch,12,cap);assert len(cert['core'])==expected
        assert set(cert['core'])<=G
        assert all(r['missing'] not in batch for r in cert['witnesses'])
        rank_cores.append(dict(rank=cap,certificate=cert))
    base_minors=[];base_strong=[]
    for e in range(12):
        row=dict(coordinate=e,certificates=[])
        for name,M in [('zero',section(G,e,0)),('one',section(G,e,1)),('contract',contraction(G,e))]:
            order=greedy_order(M,11);assert order is not None
            row['certificates'].append(dict(operation=name,family=sorted(M),order=order))
        base_minors.append(row)
        C=section(G,e,0)&section(G,e,1);order=greedy_order(C,11);assert order is not None
        base_strong.append(dict(coordinate=e,family=sorted(C),order=order))
    batch_control=dict(order=305,base=sorted(G),tips=sorted(T),family=sorted(batch),acquired_supports=acquired,
                       rank_cores=rank_cores,base_proper_minor_orders=base_minors,base_strong_orders=base_strong,
                       conclusion='All 2^14 embedded subsets are terminal forbidden families by uniform compatible-batch and protected-witness proofs. No subset census or SAT run.')
    print('305 batch: rank-2 core 0, rank-3 core 248, rank-4/full core 281; base proper and strong certificates replay',flush=True)
    deep=frozenset(src[1]['cube_completion']['family']);assert not direct(deep,24)
    out=dict(schema='first-deletion-constraints-v1',all_checks_passed=True,
             source_sha256={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
             verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),small_checks=small,
             controls=controls,order267_proper_minor_orders=minors,retained_core_controls=retained,canonical_core_controls=canonical,compatible_batch_control=batch_control,
             deeper_connected_seed=dict(family=sorted(J),order=293,acquired_support=585,negative_states=states,
                                        terminal_and_minor_scope="Uniform ample-addition stability and atlas applied to certified H292; no new full minor scan."),
             deeper_negative_control=dict(order=580,accepted=False),
             scope='Exact one-deletion and protected-cube branches; all 2^14 compatible batch outputs follow from a uniform proof, not subset enumeration. Global seed exhaustion remains open.')
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print('PASS first-deletion clauses, calibration, and positive minor replay',flush=True)
if __name__=='__main__':main()
