import sys,json,heapq,hashlib
from pathlib import Path
r=Path(__file__).resolve().parents[3];sys.path.insert(0,str(r/'papers/dismantlable-ample/research'));sys.path.insert(0,str(r/'forbidden_pc_minor_campaign'))
from cube_corner_compatibility import cube_code
from verify_damp_alternating_reductions import corner,columns,cube_shadow
p=r/'papers/dismantlable-ample/research/cube_corner_compatibility.json';d=json.loads(p.read_text());A=set(d['A']);C=set(d['C']);n=25;H=set(cube_code(A,C,23,2,[0,3,1]));ground=range(n)
print('family',len(H),flush=True);sh=cube_shadow(H,ground);cl,groups,_=columns(H,ground,sh);print('columns',groups,flush=True)
def greedy(F):
 L=set(F);heap=[v for v in L if corner(L,v,ground)];heapq.heapify(heap);queued=set(heap);order=[]
 while heap:
  v=heapq.heappop(heap);queued.remove(v)
  if not corner(L,v,ground):continue
  L.remove(v);order.append(v)
  for e in ground:
   u=v^(1<<e)
   if u in L and u not in queued and corner(L,u,ground):heapq.heappush(heap,u);queued.add(u)
 return order,L
rows=[]
for e in ground:
 for b in [0,1,None]:
  M={x&~(1<<e) for x in H if b is None or (x>>e&1)==b};o,left=greedy(M)
  rows.append(dict(coordinate=e,side=b,order=o,residual=sorted(left)))
  print(e,b,'positive' if not left else 'unknown',len(M),len(left),flush=True)
out=dict(codes=[0,3,1],family=sorted(H),dimension=n,columns=groups,signed_clauses=cl,minor_rows=rows,source_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),scope='Greedy orders are positive certificates; nonempty residuals are UNKNOWN, not negative claims.')
(r/'papers/dismantlable-ample/research/three_component_code_minors.json').write_text(json.dumps(out,indent=2)+'\n')
