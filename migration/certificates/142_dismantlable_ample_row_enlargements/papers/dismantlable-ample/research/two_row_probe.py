"""Test opposite nested row patterns in fixed certified obstructions.
No absent pattern is interpreted as a general nonexistence result.
"""
import json,itertools,time,hashlib,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'papers/ample-corners/scripts'))
from verify_rectangular_obstruction import peel
from verify_shrinking_rooted_inputs import replay
HERE=Path(__file__).resolve().parent

def patterns(X,n):
 out=[]
 for e,f in itertools.combinations(range(n),2):
  mask=(1<<e)|(1<<f);rows=[set() for _ in range(4)]
  for x in X:rows[2*((x>>e)&1)+((x>>f)&1)].add(x&~mask)
  for flip in (0,1):
   B0,D0,D1,B1=[rows[i^flip] for i in range(4)];C=B0|B1
   if not C or not C<=D0&D1:continue
   assert D0&D1==C, 'Ample family violates intersection identity'
   H=B0|{x|(1<<f) for x in B1}
   try:seq=peel(H,n)
   except AssertionError as err:
    if not err.args or not isinstance(err.args[0],tuple) or err.args[0][0]!='stalled':raise
    seq=None
   if seq is not None:replay(H,seq,n)
   out.append(dict(pair=[e,f],flip=flip,row_orders=list(map(len,[B0,D0,D1,B1])),H_order=len(H),H_positive=seq is not None,H_order_certificate=seq))
 return out

def main():
 start=time.monotonic();res=[]
 sources=[('record267','papers/dismantlable-ample/evidence/order267_obstruction/record.json','family',12),('rectangle1364','papers/ample-corners/evidence/rectangular_boundary_1364/certificate.json','family',12),('row3348','papers/dismantlable-ample/research/one_row_enlargement.json','family',13)]
 for name,path,key,n in sources:
  p=ROOT/path;d=json.loads(p.read_text());X=set(d[key]);ps=patterns(X,n)
  res.append(dict(name=name,path=path,source_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),order=len(X),dimension=n,patterns=ps));print(name,len(ps),'patterns',sum(p['H_positive'] for p in ps),'positive reduced inputs',flush=True)
 report=dict(status='PASS',results=res,seconds=time.monotonic()-start,scope='All unordered coordinate pairs and the two diagonal choices of each fixed input; no family-wide coverage. A failed greedy order is UNKNOWN, not a nonpeelability proof.')
 (HERE/'two_row_probe.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
