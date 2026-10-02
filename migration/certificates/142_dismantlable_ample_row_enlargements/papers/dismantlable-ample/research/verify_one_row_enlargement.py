"""Independent USO/minor replay for the3348-word row enlargement."""
from pathlib import Path
import json,hashlib,time,sys
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'papers/ample-corners/scripts'))
from verify_full_base_defect import corners
from explore_seed267_endpoints import verify
HERE=Path(__file__).resolve().parent

def main():
 start=time.monotonic();src=HERE/'one_row_enlargement.json';r=json.loads(src.read_text())
 for fn,h in r['source_sha256'].items():assert hashlib.sha256((ROOT/fn).read_bytes()).hexdigest()==h
 V=set(r['family']);H=set(r['H']);C=set(r['C']);D=set(r['D']);assert C<D
 actual=[{x&2047 for x in V if x>>11==t} for t in range(4)]
 assert actual==[{x&2047 for x in H if not x&2048},D,C,{x&2047 for x in H if x&2048}]
 # Projection bitsets count all shattered supports, independently of set traces.
 n=13;N=1<<n;P=[0]*N;P[N-1]=sum(1<<v for v in V)
 zero=[sum(1<<v for v in range(N) if not v>>i&1) for i in range(n)]
 for s in range(N-2,-1,-1):
  b=((N-1)^s)&-((N-1)^s);i=b.bit_length()-1;v=P[s|b];P[s]=(v|(v>>b))&zero[i]
 assert sum(P[s].bit_count()==1<<s.bit_count() for s in range(N))==len(V)==3348
 assert not corners(V,13) and set(corners(D,11))<=C
 verify(D,r['D_order'],11)
 keys=set()
 for m in r['proper_minor_orders']:
  q=m['coordinate'];b=m['operation'];keys.add((q,b));M={x&~(1<<q) for x in V if b is None or (x>>q)&1==b}
  verify(M,m['order'],13)
 assert keys=={(q,b) for q in range(13) for b in [0,1,None]}
 report=dict(status='PASS',order=len(V),proper_minor_orders=39,cornerless=True,ample=True,positive_enlarged_row=True,no_relative_first_move=True,seconds=time.monotonic()-start,source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Direct independent ampleness/cornerlessness and every39 elementary-minor USO order. Symbolic theorem supplies sound generation and marked recovery for the stated four-row subfamily.')
 (HERE/'one_row_enlargement_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)
if __name__=='__main__':main()
