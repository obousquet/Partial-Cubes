"""Bounded exact closure of the <=2-corner mutation component of one Hall projection.

No isomorphism quotient. All legal maximum VC3 single-word exchanges are
examined at each processed state. A completed queue certifies that leaving
this component requires at least three corners, unless a <=1-corner state
is found. A time/node limit is PARTIAL, not a barrier certificate.
"""
import collections,hashlib,json,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'forbidden_pc_minor_campaign'))
from explore_damp_maximum_vc3_corner_mutations import hall_projection,build_state,legal_swaps,swap_corner_count,subsets

def main():
 start=time.monotonic();seed=hall_projection(0);families=[seed];ids={seed:0};queue=collections.deque([0]);rows=[];found=None
 while queue and len(rows)<128 and time.monotonic()-start<50:
  i=queue.popleft();X=families[i];state=build_state(X,11);edges=[]
  # Independent local corner definition, not maximal-cube incidence.
  direct={v for v in X if all(v^t in X for t in subsets(sum(1<<j for j in range(11) if v^(1<<j) in X)))}
  assert direct==set(state.corners)
  for removed,added in legal_swaps(X,11):
   count=swap_corner_count(state,removed,added,11);target=None
   if count<=2:
    Y=(X-{removed})|{added}
    if Y not in ids:ids[Y]=len(families);families.append(Y);queue.append(ids[Y])
    target=ids[Y]
    assert len(build_state(Y,11).corners)==count
    if count<2:found=target
   edges.append([removed,added,count,target])
  rows.append(dict(id=i,corners=sorted(state.corners),exchanges=edges))
  if found is not None:break
 out=ROOT/'papers/ample-corners/evidence/maximum_two_corner_component';out.mkdir(parents=True,exist_ok=True)
 checks=dict(families=[sorted(X) for X in families],processed=rows,pending=list(queue),found=found)
 (out/'checks.json').write_text(json.dumps(checks,separators=(',',':'))+'\n')
 boundary=collections.Counter(e[2] for r in rows for e in r['exchanges'] if e[2]>2)
 report=dict(status='IMPROVEMENT' if found is not None else ('COMPLETE' if not queue else 'PARTIAL'),seconds=time.monotonic()-start,seed='Hall projection at coordinate0',dimension=11,order=232,states=len(families),processed=len(rows),boundary_corner_histogram=dict(boundary),checks_sha256=hashlib.sha256((out/'checks.json').read_bytes()).hexdigest(),sources={})
 for source in [Path(__file__),ROOT/'forbidden_pc_minor_campaign/explore_damp_maximum_vc3_corner_mutations.py',ROOT/'explore_smallest_additional_peripheral_obstruction.py']:
  data=source.read_bytes();(out/source.name).write_bytes(data);report['sources'][str(source.relative_to(ROOT))]=hashlib.sha256(data).hexdigest()
 (out/'manifest.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)

if __name__=='__main__':main()
