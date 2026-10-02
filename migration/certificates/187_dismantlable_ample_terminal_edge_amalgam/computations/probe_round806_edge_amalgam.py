import hashlib, json, sys
from pathlib import Path
sys.path.insert(0, str(Path.cwd()/'forbidden_pc_minor_campaign'))
from verify_damp_layer_midpoint_cores import corner
from verify_damp_terminal_cut_vertex import elementary
from verify_damp_rank_three_corner_neighbors import corners
root=Path('forbidden_pc_minor_campaign')
source=json.loads((root/'damp_terminal_cut_vertex_verification.json').read_text())
C=set(source['core'])

def peel(F, keep, mode=0):
    R=set(F); out=[]
    priority=sorted(R-keep, key=lambda x: ((-min((x^v).bit_count() for v in keep)) if mode==0 else (x if mode==1 else -x), x))
    while len(R)>len(keep):
        v=next((x for x in priority if x in R and corner(R,x,12)),None)
        if v is None: return None, sorted(R)
        R.remove(v);out.append(v)
    return out, sorted(R)

records=[]
for e in [1]:
    keep={0,1<<e}; tests=[]
    cases=[(f,k) for f in range(12) if f!=e for k in ['root_section','contraction']]
    cases += [(e,'root_section'),(e,'opposite_section'),(e,'contraction')]
    for f,k in cases:
        F=elementary(C,f,k)
        K=elementary(keep,f,k)
        for mode in range(3):
            order, rem=peel(F,K,mode)
            if order is not None:break
        tests.append({'coordinate':f,'operation':k,'order':order,'remaining':rem if order is None else sorted(K)})
    failed=[(r['coordinate'],r['operation'],len(r['remaining'])) for r in tests if r['order'] is None]
    print({'edge':e,'cases':len(tests),'failed':failed},flush=True)
    repairs=[]
    for a in [128,290]:
        for mode in range(3):
            order, rem=peel(C|{a},keep,mode)
            if order is not None:break
        repairs.append({'tip':a,'order':order,'remaining':rem})
    print({'edge':e,'repair_failures':[(r['tip'],len(r['remaining'])) for r in repairs if r['order'] is None]},flush=True)
    records.append({'edge':e,'tests':tests,'repairs':repairs})
report={'schema':'damp-edge-amalgam-orders-v1', 'source_sha256':hashlib.sha256((root/'damp_terminal_cut_vertex_verification.json').read_bytes()).hexdigest(), 'probe_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), 'records':records, 'scope':'Bounded positive-order search for one fixed edge; failed greedy searches would not prove nonpeelability.'}
(root/'computations/round806_edge_amalgam_orders.json').write_text(json.dumps(report,indent=2)+'\n')
