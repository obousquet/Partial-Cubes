"""Replay a positive-input generator for the positive-overlap R_per branch.

All output negative decisions use common-periphery absence and a nonconstant
component split. All old minor certificates use uncoloured face forests.
Direct output recognition is a verification comparison only.
"""
from pathlib import Path
import sys,json,hashlib
ROOT=Path(__file__).resolve().parents[3]
sys.path[:0]=[str(ROOT),str(Path(__file__).parent)]
from component_label_recursion import components,circuit,evaluate,lift,shadow,certificate
from explore_rper_cyclic_families import section,contraction,rper,normalize

def check_tree(t):
    H=set(t['family']);n=t['dimension']
    if t.get('leaf'):assert len(H)<=1;return
    e,b=t['coordinate'],t['small_side'];small=section(H,e,b);large=section(H,e,1-b)
    assert small<=large
    for child,F in [('site',small),('base',large)]:
        h,m=normalize(F,n-1);assert t[child]['dimension']==m and set(t[child]['family'])==set(h)
        check_tree(t[child])

def positive(H,n):
    t=certificate(H,n);check_tree(t)
    return dict(original_family=sorted(H),original_dimension=n,tree=t)

def face_forest(g,labels):
    values=[]
    for node in g['nodes']:
        values.append(len({labels[i] for i in node['mono']})<=1 or any(values[a] and values[b] for _,_,a,b in node['cuts']))
    def build(i):
        q=g['nodes'][i]
        if len({labels[j] for j in q['mono']})<=1:
            return dict(remaining=q['remaining'],fixed=q['fixed'],leaf_components=q['mono'])
        e,b,a,c=next(cut for cut in q['cuts'] if values[cut[2]] and values[cut[3]])
        return dict(remaining=q['remaining'],fixed=q['fixed'],coordinate=e,orientation=b,children=[build(a),build(c)])
    assert values[g['root']]
    return build(g['root'])

def check_forest(A,C,ks,n,t,q):
    mask=((1<<n)-1)^q['remaining'];a={x for x in A if x&mask==q['fixed']};c={x for x in C if x&mask==q['fixed']}
    if 'leaf_components' in q:
        ids=[i for i,K in enumerate(ks) if K & (a-c)]
        assert ids==q['leaf_components'] and len({t[i] for i in ids})<=1
        return
    e,b=q['coordinate'],q['orientation'];bit=1<<e;assert q['remaining']&bit
    assert all(x^bit in a for x in a if x>>e&1==b)
    assert all(x^bit in c for x in c if x>>e&1==b)
    for side,child in enumerate(q['children']):
        assert child['remaining']==q['remaining']^bit and child['fixed']==q['fixed']|(side<<e)
        check_forest(A,C,ks,n,t,child)

def common(A,C,n):
    return [(e,b) for e in range(n) for b in (0,1)
            if all(x^(1<<e) in A for x in A if (x>>e)&1==b)
            and all(x^(1<<e) in C for x in C if (x>>e)&1==b)]

def run(H,n,w):
    B=[set(section(H,w,b)) for b in (0,1)];A=B[0]|B[1];C=B[0]&B[1];m=n-1
    assert C and not common(A,C,m)
    trees={k:positive(F,m) for k,F in [('A',A),('C',C),('B0',B[0]),('B1',B[1])]}
    blocks=components(A,C,m);labels=[]
    for K in blocks:
        assert K<=B[0] or K<=B[1];labels.append(0 if K<=B[0] else 1)
    assert len(set(labels))==2
    X=lift(A,C,blocks,labels,m);assert len(shadow(X,n))==len(X)
    # Marked coordinate placed last; this is the recovered output embedding.
    assert contraction(X,m)==frozenset(A)
    tests=[]
    for e in range(m):
        for kind in ['restriction0','restriction1','contraction']:
            mu=lambda F:contraction(F,e) if kind=='contraction' else section(F,e,int(kind[-1]))
            a,c=set(mu(A)),set(mu(C));Y=mu(X)
            assert set(contraction(Y,m-1))==a and set(section(Y,m-1,0))&set(section(Y,m-1,1))==c
            ks=components(a,c,m-1);ys=[set(section(Y,m-1,b)) for b in (0,1)]
            t=[0 if K<=ys[0] else 1 for K in ks];assert all(K<=ys[t[i]] for i,K in enumerate(ks))
            g=circuit(a,c,m-1);assert evaluate(g,t)
            f=face_forest(g,t);check_forest(a,c,ks,m-1,t,f)
            # Independent verification only; not the construction's admissibility test.
            assert rper(*normalize(Y,m))
            tests.append(dict(coordinate=e,operation=kind,labels=t,face_forest=f))
    assert not rper(*normalize(X,n))
    terminal=all(rper(*normalize(set(section(X,e,0))&set(section(X,e,1)),n-1)) for e in range(n))
    return dict(family=sorted(X),dimension=n,order=len(X),marked_coordinate=w,
                positive_inputs=trees,components=[sorted(K) for K in blocks],labels=labels,
                common_peripheral_orientations=[],old_minor_face_forests=tests,
                terminal_strong_projection_status=terminal)

def main():
    source=ROOT/'papers/cover-recursive/research/component_label_recursion.json';old=json.loads(source.read_text())
    P5={0,1,2,3,4,6,7,8,12,14,16,17,19,24,25,28}
    examples=[]
    for H,n in [(P5,5),(set(old['example']['family']),6)]:
        w=next(e for e in range(n) if rper(*normalize(set(section(H,e,0))&set(section(H,e,1)),n-1)))
        r=run(frozenset(H),n,w);examples.append(r)
        print('output',len(H),'marked',w,'positive overlap',len(r['positive_inputs']['C']['original_family']),
              'all',len(r['old_minor_face_forests']),'old-minor face forests PASS; terminal',r['terminal_strong_projection_status'],flush=True)
    result=dict(all_checks_passed=True,examples=examples,
                source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                scope='Fixed terminal P5 and nonterminal 37 controls. General coverage is the positive-overlap theorem and lower-dimensional recursion. Fixed tests do not themselves establish coverage or publication novelty.')
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
