"""An eligible both-nonrelative two-row presentation from a nested ample triple.
Soundness uses two established one-row steps. All42 final minor orders replay.
Relative failures use coordinate faces equal to the certified42/22 pair.
"""
import json,sys,hashlib,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'papers/ample-corners/scripts'))
from verify_rectangular_obstruction import peel
from verify_shrinking_rooted_inputs import replay
from verify_full_base_defect import shadow,corners
from one_row_enlargement import rows,lift,image,transport,relative_replay
HERE=Path(__file__).resolve().parent

def expand(H,orders,D,dorder,n):
 z=n-1;f=1<<z;e=1<<n;V=rows(H,D,z);out=[]
 def get(q,b):return next(r['order'] for r in orders if r['coordinate']==q and r['operation']==b)
 for q in range(z):
  for b in (0,1,None):
   MH=image(H,q,b);MD=image(D,q,b)
   seq=lift(MH,get(q,b),z)+[x|f for x in transport(dorder,q,b)]
   replay(image(V,q,b),seq,n+1)
   out.append(dict(coordinate=q,operation=b,order=seq))
 B0=get(z,0);B1=get(z,1);co=get(z,None)
 special={(n,0):B0+[x|f for x in dorder],(n,1):[x|f for x in B1]+co,(n,None):co+[x|f for x in dorder],(z,0):B0+[x|e for x in co],(z,1):[x|e for x in B1]+dorder,(z,None):[x|e for x in co]+dorder}
 for (q,b),seq in special.items():
  replay(image(V,q,b),seq,n+1);out.append(dict(coordinate=q,operation=b,order=seq))
 return V,out

def main():
 start=time.monotonic()
 ps=ROOT/'papers/ample-corners/evidence/relative_pair_q6/manifest.json';hs=ROOT/'papers/ample-corners/evidence/rectangular_boundary_1364/certificate.json'
 pair=json.loads(ps.read_text());cert=json.loads(hs.read_text())
 P=set(pair['outer']);O=set(pair['inner']);P0={x>>1 for x in P};O0={x>>1 for x in O};T=set(range(32))
 normalize=lambda x:(x>>1)|((x&1)<<11)
 H={normalize(x) for x in cert['family']}
 orders=[dict(coordinate=(r['coordinate']-1)%12,operation=None if r['operation']=='contract' else int(r['operation'][-1]),order=[normalize(x) for x in r['peeling']]) for r in cert['minor_certificates']]
 C={x&2047 for x in H}
 D={x|(y<<5) for x in P0 for y in P}|{x|(y<<5) for x in T for y in O}
 E={x|(y<<5) for x in T for y in P}
 assert C<D<E and (len(C),len(D),len(E))==(892,1224,1344)
 po=peel(P,6);oo=peel(O,6);p0o=peel(P0,5);to=peel(T,5)
 prefix=[x|(y<<5) for x in peel(T,5,P0) for y in oo]
 product={x|(y<<5) for x in P0 for y in P}
 relative_replay(D,prefix,product,11)
 do=prefix+[x|(y<<5) for x in p0o for y in po];eo=[x|(y<<5) for x in to for y in po]
 replay(D,do,11);replay(E,eo,11)
 # Exact face witnesses: neither relative peeling is possible.
 u=min(P0-O0);v=min(T-P0)
 section=lambda A,x:{w>>5 for w in A if w&31==x}
 assert section(D,u)==P and section(C,u)==O
 assert section(E,v)==P and section(D,v)==O
 assert set(corners(P,6))<=O and P!=O
 X,xorders=expand(H,orders,D,do,12)
 L=E|{x|2048 for x in D};lo=[x|2048 for x in do]+eo;replay(L,lo,12)
 Z,zorders=expand(X,xorders,L,lo,13)
 assert len(X)==3480 and len(Z)==len(shadow(Z,14))==8164
 # Read rows on (b,f)=(12,11), retaining a=13 as an old coordinate.
 rr=[{x&~(2048|4096) for x in Z if 2*((x>>12)&1)+((x>>11)&1)==t} for t in range(4)]
 U0,R0,R1,U1=rr;core=U0|U1
 assert R0&R1==core
 assert R0==D|{x|8192 for x in D}
 assert R1==E|{x|8192 for x in C}
 assert core==D|{x|8192 for x in C}
 Hstar=U0|{x|2048 for x in U1}
 assert Hstar=={(x&~4096)|((x&4096)<<1) for x in X}
 # Reduction in a proves first relative failure; a=0 face proves the second.
 reduce_a=lambda A:{x for x in A if not x&8192 and x|8192 in A}
 assert reduce_a(R0)==D and reduce_a(core)==C
 assert {x for x in R1 if not x&8192}==E
 assert {x for x in core if not x&8192}==D
 result=dict(status='PASS',order=len(Z),dimension=14,reduced_input_order=len(Hstar),reduced_input_corners=len(corners(X,13)),large_row_orders=[len(R0),len(R1)],intersection_order=len(core),triple_orders=[len(C),len(D),len(E)],relative_failure_face_words=[u,v],family=sorted(Z),proper_minor_orders=zorders,first_step_orders=xorders,positive_D_order=do,positive_E_order=eo,source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [ps,hs]},seconds=time.monotonic()-start,scope='Two one-row construction identities; all42 final elementary-minor orders and ampleness replayed. Both relative failures follow from recorded42/22 face witnesses, not greedy stalls. No new primitive input or global exhaustiveness.')
 (HERE/'two_row_both_nonrelative.json').write_text(json.dumps(result,separators=(',',':'))+'\n')
 print(json.dumps({k:v for k,v in result.items() if k not in ['family','proper_minor_orders','first_step_orders','positive_D_order','positive_E_order']}),flush=True)
if __name__=='__main__':main()
