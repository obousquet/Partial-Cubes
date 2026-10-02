"""Fixed diagnostic for the cornerless-minor descent claim; failures are inconclusive."""
from pathlib import Path
import json
import hashlib
from time import monotonic

ROOT=Path(__file__).resolve().parents[1]
d=json.loads((ROOT/'damp_empty_core_counterexample.json').read_text())
c=d['corner'];H={x^c for x in d['counterexample']}

def corner(X,x):
    mask=sum(1<<e for e in range(12) if x^(1<<e) in X)
    sub=mask
    while sub:
        if x^sub not in X:return False
        sub=(sub-1)&mask
    return True

def greedy(X,reverse):
    left=set(X);order=[]
    while len(left)>1:
        choices=[x for x in sorted(left,reverse=reverse) if x and corner(left,x)]
        if not choices:return None,sorted(left)
        x=choices[0];left.remove(x);order.append(x)
    return order+[0],[]

assert len(H)==292 and [x for x in H if corner(H,x)]==[0]
rows=[]
for e in range(12):
    bit=1<<e
    for kind in ['root_section','contraction']:
        selected=H if kind=='contraction' else {x for x in H if not x&bit}
        child={x&~bit for x in selected}
        tried=[];order=None
        for reverse in [False,True]:
            order,remaining=greedy(child,reverse)
            tried.append({'reverse':reverse,'remaining':remaining})
            if order:break
        row={'coordinate':e,'operation':kind,'order':len(child),'rooted_peeling':order,'greedy_trials':tried}
        rows.append(row)
        print(e,kind,len(child),'positive' if order else 'unresolved',flush=True)
(ROOT/'computations/round818_one_corner_rooted_images.json').write_text(json.dumps({'source_sha256':hashlib.sha256((ROOT/'damp_empty_core_counterexample.json').read_bytes()).hexdigest(),'generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'corner':c,'family':sorted(H),'rows':rows,'scope':'Fixed root-retaining greedy diagnostic; failures are unresolved, not negative proofs.'},indent=2)+'\n')
print('positive',sum(r['rooted_peeling'] is not None for r in rows),'of',len(rows),flush=True)
