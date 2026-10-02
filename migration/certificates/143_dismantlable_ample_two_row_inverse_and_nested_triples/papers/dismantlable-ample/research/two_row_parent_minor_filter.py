"""Necessary cardinality test for two-row inverse failures over fixed H3480.

The marked-reduction theorem forces every12-dimensional forbidden proper
minor to have a coordinate reduction isomorphic to an at-most-two-operation
old-coordinate minor of H3480 (ambient dimension13). Enumerate all such
labelled operations, retaining constants; cardinalities alone are necessary,
never sufficient. Targets: recorded267 (all reductions59) and rectangular
1364 family (all reductions472). No statement about unknown obstructions.
"""
import json,hashlib,itertools,time
from pathlib import Path
from one_row_enlargement import rows,image
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent

def main():
 start=time.monotonic()
 ps=ROOT/'papers/ample-corners/evidence/relative_pair_q6/manifest.json'
 hs=ROOT/'papers/ample-corners/evidence/rectangular_boundary_1364/certificate.json'
 pair=json.loads(ps.read_text());h=json.loads(hs.read_text())
 P=set(pair['outer']);O=set(pair['inner']);P0={x>>1 for x in P}
 H={(x>>1)|((x&1)<<11) for x in h['family']}
 D={x|(y<<5) for x in P0 for y in P}|{x|(y<<5) for x in range(32) for y in O}
 X=rows(H,D,11);assert len(X)==3480
 records=[]
 for length in (1,2):
  for coords in itertools.combinations(range(13),length):
   for ops in itertools.product((0,1,None),repeat=length):
    Y=set(X)
    for q,b in zip(coords,ops):Y=image(Y,q,b)
    mask=sum(1<<q for q in coords)
    direct={x&~mask for x in X if all(b is None or (x>>q)&1==b for q,b in zip(coords,ops))}
    assert Y==direct
    active=sum(any(x^(1<<q) in Y for x in Y) for q in range(13))
    records.append(dict(coordinates=coords,operations=ops,order=len(Y),active_coordinates=active))
 assert len(records)==39+702
 target_results={}
 for order,reduction in [(267,59),(1364,472)]:
  hits=[r for r in records if r['order']==reduction]
  target_results[str(order)]=dict(required_reduction_order=reduction,cardinality_hits=hits,excluded_by_cardinality=not hits)
 report=dict(status='PASS',direct_trace_replay=True,parent_order=len(X),parent_dimension=13,records=records,targets=target_results,minimum_nonempty_order=min(r['order'] for r in records if r['order']),source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (ps,hs)},seconds=time.monotonic()-start,scope=__doc__)
 (HERE/'two_row_parent_minor_filter.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps({k:v for k,v in report.items() if k not in ('records','scope')}),flush=True)
if __name__=='__main__':main()
