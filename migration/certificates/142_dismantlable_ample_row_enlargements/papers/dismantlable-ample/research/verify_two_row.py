"""Independent two-row identity, USO input and fiber-exclusion certificates."""
from pathlib import Path
import sys,json,hashlib,itertools,time
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'papers/ample-corners/scripts'))
from verify_full_base_defect import shadow,submasks
from verify_rectangular_obstruction import peel
from explore_seed267_endpoints import verify,iscorner

def main():
 start=time.monotonic();toy=0
 for hm in range(256):
  H={x for x in range(8) if hm>>x&1}
  if len(shadow(H,3))!=len(H):continue
  B0={x&3 for x in H if not x&4};B1={x&3 for x in H if x&4};C=B0|B1
  supers=[{x for x in range(4) if d>>x&1} for d in range(16) if all(d>>x&1 for x in C)]
  for D0 in supers:
   for D1 in supers:
    K=D1|{x|4 for x in D0};V=B0|{x|4 for x in D0}|{x|8 for x in D1}|{x|12 for x in B1}
    assert (len(shadow(V,4))==len(V)) == (D0&D1==C and len(shadow(K,3))==len(K));toy+=1
 src=ROOT/'papers/dismantlable-ample/evidence/order267_obstruction/record.json';X=set(json.loads(src.read_text())['family']);assert len(shadow(X,12))==267 and not any(iscorner(X,x,12) for x in X)
 rows=[{x&~3 for x in X if 2*(x&1)+((x>>1)&1)==(t^1)} for t in range(4)]
 B0,D0,D1,B1=rows;C=B0|B1;H=B0|{x|2 for x in B1};K=D1|{x|2 for x in D0}
 assert list(map(len,rows))==[13,136,72,46] and len(C)==48 and D0&D1==C
 orders={}
 for name,Y in [('H',H),('K',K)]:
  assert len(shadow(Y,12))==len(Y);seq=peel(Y,12);verify(Y,seq,12);orders[name]=seq
 # Independently enumerate fibers by support/trace, rather than iterative intersections.
 probe=json.loads((HERE/'two_row_probe.json').read_text());pats=probe['results'][0]['patterns'];deletions=set()
 for pat in pats:
  e,f=pat['pair'];flip=pat['flip'];mask=(1<<e)|(1<<f);coords=[j for j in range(12) if j not in (e,f)]
  rs=[{x&~mask for x in X if 2*((x>>e)&1)+((x>>f)&1)==(t^flip)} for t in range(4)];C=rs[0]|rs[3]
  for i in (1,2):
   for sm in range(1<<len(coords)):
    support=sum(1<<j for q,j in enumerate(coords) if sm>>q&1);forbidden={x&support for x in C};groups={}
    for x in rs[i]:groups.setdefault(x&support,set()).add(x)
    for trace,G in groups.items():
     if trace in forbidden:continue
     t=i^flip;deletions.add(tuple(sorted(x|((t>>1)<<e)|((t&1)<<f) for x in G)))
 r=json.loads((HERE/'two_row_fiber_shrink.json').read_text());assert len(deletions)==r['unique_outputs'] and not r['hits']
 # Precompute seed cubes only to find short witnesses. Check each witness on Y directly.
 supports=sorted(shadow(X,12),key=lambda s:(s.bit_count(),s));cube_sets={s:[{x^t for t in submasks(s)} for x in X if x&s==0 and all(x^t in X for t in submasks(s))] for s in supports}
 certs=[];counts={'corner':0,'nonample':0}
 for rem in sorted(deletions):
  Y=X-set(rem);v=next((x for x in Y if iscorner(Y,x,12)),None)
  if v is not None:
   counts['corner']+=1
  s=next(s for s in supports if len({x&s for x in Y})==1<<s.bit_count() and not any(Q<=Y for Q in cube_sets[s]))
  # A shattered but not strongly shattered support refutes ampleness.
  assert len({x&s for x in Y})==1<<s.bit_count()
  assert all(any(x^t not in Y for t in submasks(s)) for x in Y if x&s==0)
  certs.append(dict(removed=rem,shattered_not_strong=s));counts['nonample']+=1
 report=dict(status='PASS',toy_pairs=toy,seed_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),positive_input_orders=orders,row_sizes=[13,136,72,46],C_order=48,H_order=len(H),K_order=len(K),unique_fiber_deletions=len(deletions),exclusion_counts=counts,certificates=certs,seconds=time.monotonic()-start)
 (HERE/'two_row_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k not in ['certificates','positive_input_orders']}),flush=True)
if __name__=='__main__':main()
