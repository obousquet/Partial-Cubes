"""Replay a 909-word tree-reduction obstruction from two existing inputs.
No search or negative DFS. Negativity uses the first intact-path deletion,
the unique core/root corner, and the established asymmetric edge profile.
"""
from pathlib import Path
import json,hashlib
from component_label_partitions import corner,shadow
HERE=Path(__file__).resolve().parent

def outmap(A,order,n):
 assert len(order)==len(A) and set(order)==set(A)
 left=set(A);r={}
 for u in order:
  assert corner(left,u,n),('not corner',u)
  r[u]=sum(1<<j for j in range(n) if u^(1<<j) in left)
  left.remove(u)
 return r

def topological(r,n):
 deg=dict.fromkeys(r,0)
 for u,S in r.items():
  for j in range(n):
   if S>>j&1:
    assert u^(1<<j) in r
    deg[u^(1<<j)]+=1
 todo=sorted(u for u,d in deg.items() if d==0);order=[]
 while todo:
  u=todo.pop();order.append(u)
  for j in range(n):
   if r[u]>>j&1:
    v=u^(1<<j);deg[v]-=1
    if not deg[v]:todo.append(v)
 assert len(order)==len(r)
 return order

def main():
 src=HERE/'relative_wing_square.json';psrc=HERE/'three_piece_edge_witness.json'
 data=json.loads(src.read_text())['controls'][0];pd=json.loads(psrc.read_text())
 B=set(data['core']);P=set(pd['family']);H={2,0,4,5};a,b,c,d=2,0,4,5
 assert B==set(pd['core']) and len(B)==289 and len(P)==330
 assert [u for u in sorted(B) if corner(B,u,12)]==[0]
 assert pd['edge']==[0,4096] and pd['profile']==[0,1,0,1,0]
 # Deterministic certificate discovery replay: retain a, and delay d until leaf exit.
 left=set(B);bo=[]
 while len(left)>1:
  choices=[u for u in sorted(left) if u!=a and corner(left,u,12)
           and (u!=d or {u^(1<<j) for j in range(12) if u^(1<<j) in left}=={c})]
  assert choices
  u=choices[0];bo.append(u);left.remove(u)
 assert left=={a};bo.append(a)
 br=outmap(B,bo,12);assert bo[0]==b and br[a]==0 and br[d]==1
 ao=data['ordinary_core_order'];ar=outmap(B,ao,12)
 po=pd['sink_order'];pr=outmap(P,po,13);assert po[-1]==0
 # The upper shore at P's second endpoint is B, whose unique corner is its root.
 assert {v&4095 for v in P if v&4096}==B
 # Thus the established shore argument excludes a leaf exit at that endpoint.
 embA=lambda v:a^(v<<12)
 embP=lambda v:c^((v&4095)<<24)^(1 if v&4096 else 0)
 supP=lambda S:((S&4095)<<24)|(1 if S&4096 else 0)
 Awords={embA(v) for v in B};Pwords={embP(v) for v in P}
 U=B|Awords|Pwords
 assert B&Awords=={a} and B&Pwords=={c,d} and Awords&Pwords==set()
 r=dict(br)
 for v,S in ar.items():r[embA(v)]=r.get(embA(v),0)|(S<<12)
 for v,S in pr.items():r[embP(v)]=r.get(embP(v),0)|supP(S)
 uo=topological(r,36);assert outmap(U,uo,36)==r
 # Independent representation-map check, with exact shadows from block gluing.
 shB=shadow(B,12);shP=shadow(P,13)
 shA={S<<12 for S in shB};shPe={supP(S) for S in shP}
 shU=shB|shA|shPe
 assert len(shB)==289 and len(shP)==330 and len(shU)==len(U)==905
 assert len(set(r.values()))==905 and set(r.values())==shU
 assert all((u^v)&(r[u]|r[v]) for u in U for v in U if u<v)
 E=1<<36
 upper=H|Awords|Pwords
 X=B|{v|E for v in upper}
 assert len(upper)==620 and len(X)==909 and (B&upper)==H and (B|upper)==U
 shUpper=shA|shPe|{2,4};shH={0,1,2,4}
 assert len(shUpper)==620 and shB&shUpper==shH
 assert len(shU|{S|E for S in shH})==len(X) # exact section shadow identity
 upper_order=[embP(v) for v in po[:-1]]+[c,b]+[embA(v) for v in ao]
 outmap(upper,upper_order,36)
 # Check the first-path-deletion geometry at all eight lifts.
 assert all(sum(v^(1<<j) in H for j in range(12))==(1 if v in {a,d} else 2) for v in H)
 assert any(a^(1<<j) in B-H for j in range(12))
 assert any(d^(1<<j) in B-H for j in range(12))
 # The certified Round 736 root is a; there are no free vertical bits.
 # Outgoing cube labels on the path are A, B, B, P respectively.
 bits={a:0,b:1,c:1,d:0};tree_support={a:0,b:2,c:4,d:1}
 R={v:r[v] for v in B-H}
 R.update({v|E:r[v] for v in upper-H})
 for v in H:
  R[v|(bits[v]*E)]=tree_support[v]|E
  R[v|((1-bits[v])*E)]=r[v]
 shX=shU|{S|E for S in shH}
 assert set(R)==X and set(R.values())==shX and len(set(R.values()))==len(X)
 assert all((u^v)&(R[u]|R[v]) for u in X for v in X if u<v)
 def path(rr,start,end,n):
  queue=[start];prev={start:None}
  for u in queue:
   if u==end:break
   for j in range(n):
    if rr[u]>>j&1:
     v=u^(1<<j)
     if v not in prev:prev[v]=u;queue.append(v)
  assert end in prev
  answer=[];u=end
  while u is not None:answer.append(u);u=prev[u]
  return answer[::-1]
 bp=path(br,b,d,12);pp=path(pr,4096,0,13)
 cycle=bp+[embP(v)|E for v in pp]+[b|E,b]
 assert all((u^v).bit_count()==1 and R[u]&(u^v) for u,v in zip(cycle,cycle[1:]))
 assert len(set(cycle[:-1]))==len(cycle)-1
 report=dict(all_checks_passed=True,source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [src,psrc]},
   canonical_lift_cycle=cycle,canonical_root=a,free_vertical_vertices=[],canonical_lift_outmap=R,
   path=[a,b,c,d],core_order=bo,core_outmap=br,upper_section_order=upper_order,contraction_order=uo,
   family=sorted(X),contraction=sorted(U),order=909,dimension=37,new_coordinate=36,
   lower_section_order=bo,negative_proof='Before the first doubled-path deletion, no lower core vertex can be deleted. Internal path lifts are not corners; lower endpoints have unmatched core neighbours. The upper near endpoint would require the rooted 289-piece to peel to its root. The upper far endpoint would give leaf-exit value 1 at the forbidden endpoint of the 330-piece. Both contradict established input properties.',
   scope='Explicit ample nonpeelable expansion with positive contraction and both positive sections. Minor-minimality is separately certified by four_vertex_path_minors.py/json, proving the four-edge bound sharp with this explicit class.')
 (HERE/'four_vertex_path_obstruction.json').write_text(json.dumps(report,indent=2)+'\n')
 print('PASS: 909 words, dimension 37; 905-word contraction and both sections peel; exact shadow and independent representation-map checks pass.')
 print('Negative proof is structural; proper-minor positivity is certified separately by four_vertex_path_minors.py.')
if __name__=='__main__':main()
