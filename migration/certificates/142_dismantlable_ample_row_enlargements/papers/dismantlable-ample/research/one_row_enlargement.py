"""Test arbitrary ample enlargement of one nonconstant binary block row.

Rows00,01,10,11 are B0,D,C,B1, where C=B0 union B1 is the
contraction of H and D contains C. A peeling of H lifts relative to D,
without deleting any word of D. Verify all small ample inputs and a
3348-word forbidden example using the certified1364 rectangle.
"""
from pathlib import Path
import sys,json,hashlib,time
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'papers/ample-corners/scripts'))
from verify_full_base_defect import shadow,corners,submasks
from verify_rectangular_obstruction import peel
from verify_shrinking_rooted_inputs import replay
HERE=Path(__file__).resolve().parent


def rows(H,D,n):
 mask=(1<<n)-1;B0={v&mask for v in H if not v>>n&1};B1={v&mask for v in H if v>>n&1};C=B0|B1
 assert C<=D
 return B0|{x|(1<<n) for x in D}|{x|(2<<n) for x in C}|{x|(3<<n) for x in B1}


def lift(H,seq,n):
 H=set(H);out=[];mask=(1<<n)-1
 assert set(seq)==H and len(seq)==len(H)
 for v in seq:
  x=v&mask;b=(v>>n)&1
  if v^(1<<n) not in H:out.append(x|(2<<n))
  out.append(x|((3 if b else 0)<<n));H.remove(v)
 return out


def relative_replay(X,seq,D,n):
 L=set(X)
 for v in seq:
  assert v in L-D
  directions=sum(1<<i for i in range(n) if v^(1<<i) in L)
  assert all(v^t in L for t in submasks(directions)),(v,len(L))
  L.remove(v)
 assert L==D


def image(X,j,b):return {x&~(1<<j) for x in X if b is None or (x>>j)&1==b}

def transport(seq,j,b):
 vals=[x&~(1<<j) for x in seq if b is None or (x>>j)&1==b]
 last={x:i for i,x in enumerate(vals)}
 return [x for i,x in enumerate(vals) if last[x]==i]


def main():
 start=time.monotonic();tests=0
 for mask in range(256):
  H={x for x in range(8) if mask>>x&1}
  if len(shadow(H,3))!=len(H):continue
  seq=peel(H,3);C={x&3 for x in H}
  for dm in range(16):
   D={x for x in range(4) if dm>>x&1}
   if not C<=D or len(shadow(D,2))!=len(D):continue
   V=rows(H,D,2);assert len(shadow(V,4))==len(V)
   relative_replay(V,lift(H,seq,2),{x|4 for x in D},4);tests+=1
 source=ROOT/'papers/ample-corners/evidence/rectangular_boundary_1364/certificate.json'
 inputs=ROOT/'papers/ample-corners/evidence/relative_pair_q6/manifest.json'
 cert=json.loads(source.read_text());inp=json.loads(inputs.read_text());oldH=set(cert['family']);P=set(inp['outer']);O=set(inp['inner'])
 normalize=lambda v:(v>>1)|((v&1)<<11)
 H={normalize(v) for v in oldH};C={v&2047 for v in H};P0={v>>1 for v in P}
 D={x|(y<<5) for x in P0 for y in P};V=rows(H,D,11)
 assert len(H)==1364 and len(C)==892 and len(D)==1092 and len(V)==3348
 assert C<D and set(corners(D,11))<=C and not corners(V,13)
 assert len(shadow(V,13))==len(V)
 p0=peel(P0,5);po=peel(P,6);dorder=[x|(y<<5) for x in p0 for y in po];replay(D,dorder,11)
 minor_records=[]
 for q in range(1,12):
  for b in [0,1,None]:
   old=next(r for r in cert['minor_certificates'] if r['coordinate']==q and r['operation']==('contract' if b is None else f'restrict-{b}'))
   MH={normalize(v) for v in image(oldH,q,b)};ho=[normalize(v) for v in old['peeling']];MD=image(D,q-1,b);do=transport(dorder,q-1,b)
   MV=image(V,q-1,b);assert MV==rows(MH,MD,11)
   seq=lift(MH,ho,11)+[v|(1<<11) for v in do];replay(MV,seq,13)
   minor_records.append(dict(coordinate=q-1,operation=b,order=seq))
 baseorders={}
 for b in [0,1,None]:
  old=next(r for r in cert['minor_certificates'] if r['coordinate']==0 and r['operation']==('contract' if b is None else f'restrict-{b}'))
  baseorders[b]=[v>>1 for v in old['peeling']]
 B0=baseorders[0];B1=baseorders[1];co=baseorders[None]
 special={(12,0):B0+[x|2048 for x in dorder],(12,1):[x|2048 for x in B1]+co,(12,None):co+[x|2048 for x in dorder],(11,0):B0+[x|4096 for x in co],(11,1):[x|4096 for x in B1]+dorder,(11,None):[x|4096 for x in co]+dorder}
 for (q,b),seq in special.items():replay(image(V,q,b),seq,13);minor_records.append(dict(coordinate=q,operation=b,order=seq))
 assert len(minor_records)==39
 report=dict(status='PASS',small_input_pairs=tests,H_order=len(H),contraction_order=len(C),enlarged_row_order=len(D),output_order=len(V),output_dimension=13,shattered_supports=len(V),corners=0,enlarged_row_corners=sorted(corners(D,11)),enlarged_row_corners_inside_C=True,family=sorted(V),H=sorted(H),C=sorted(C),D=sorted(D),D_order=dorder,proper_minor_orders=minor_records,seconds=time.monotonic()-start,source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,inputs]},script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='All ample H in Q3 and ample D in Q2 containing its contraction; constructive39 minor orders and direct ampleness/cornerlessness for one3348 output. D cannot peel to C because all its corners lie in C. General claims require the symbolic proof.')
 (HERE/'one_row_enlargement.json').write_text(json.dumps(report,separators=(',',':'))+'\n');print(json.dumps({k:v for k,v in report.items() if k not in ['family','H','C','D','D_order','proper_minor_orders','enlarged_row_corners']}),flush=True)
if __name__=='__main__':main()
