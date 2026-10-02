"""Exact first-corner compatibility for two-row assemblies; not minimality."""
from pathlib import Path
import json,sys,time,hashlib
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'papers/ample-corners/scripts'))
from verify_full_base_defect import shadow,corners

def predict(B0,B1,D0,D1,n):
 C=B0|B1;assert D0&D1==C
 H=B0|{x|(1<<n) for x in B1};K=D1|{x|(1<<n) for x in D0};mask=(1<<n)-1
 out=set();blocked=[];hpaired=[]
 for v in corners(H,n+1):
  if v^(1<<n) in H:
   x=v&mask;b=v>>n;out.add(x|((3*b)<<n));hpaired.append(v)
 kc=corners(K,n+1)
 for v,Q in kc.items():
  x=v&mask;b=v>>n;M={q&mask for q in Q}
  if x not in C:out.add(x|((2-b)<<n));continue
  targets=[B for B in (B0,B1) if x in B];T=set.intersection(*targets)
  if M<=T:out.add(x|((2-b)<<n))
  else:blocked.append(dict(K_corner=v,old_word=x,old_cube=sorted(M),missing_from_required_section=min(M-T)))
 return out,dict(H_corners=sorted(corners(H,n+1)),H_corners_incident_with_z=hpaired,K_corners=sorted(kc),blocked_K_corners=blocked)

def main():
 start=time.monotonic();tests=0
 for hm in range(256):
  H={x for x in range(8) if hm>>x&1}
  if len(shadow(H,3))!=len(H):continue
  B0={x&3 for x in H if not x&4};B1={x&3 for x in H if x&4};C=B0|B1
  ds=[{x for x in range(4) if d>>x&1} for d in range(16) if all(d>>x&1 for x in C)]
  for D0 in ds:
   for D1 in ds:
    if D0&D1!=C:continue
    K=D1|{x|4 for x in D0}
    if len(shadow(K,3))!=len(K):continue
    V=B0|{x|4 for x in D0}|{x|8 for x in D1}|{x|12 for x in B1}
    predicted,_=predict(B0,B1,D0,D1,2);assert predicted==set(corners(V,4));tests+=1
 src=ROOT/'papers/dismantlable-ample/evidence/order267_obstruction/record.json';X=set(json.loads(src.read_text())['family'])
 # Delete the two marked bits, so the base uses ten consecutive bits.
 rows=[{x>>2 for x in X if 2*(x&1)+((x>>1)&1)==(t^1)} for t in range(4)]
 B0,D0,D1,B1=rows;predicted,detail=predict(B0,B1,D0,D1,10)
 V=B0|{x|1024 for x in D0}|{x|2048 for x in D1}|{x|3072 for x in B1}
 assert predicted==set(corners(V,12))==set()
 report=dict(status='PASS',small_ample_assemblies=tests,seed_detail=detail,source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),seconds=time.monotonic()-start,scope='Exact corner formula, not an acyclic-USO gluing theorem or forbidden-minor minimality criterion.')
 (HERE/'two_row_corner_compatibility.json').write_text(json.dumps(report,indent=2)+'\n')
 print(json.dumps(dict(status='PASS',small_ample_assemblies=tests,H_corners=len(detail['H_corners']),H_paired=len(detail['H_corners_incident_with_z']),K_corners=len(detail['K_corners']),blocked=len(detail['blocked_K_corners']))),flush=True)
if __name__=='__main__':main()
