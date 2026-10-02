"""Check every ample output of the complete fixed two-row fiber deletion scan."""
from pathlib import Path
import json,sys,time,hashlib
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'papers/ample-corners/scripts'))
from verify_full_base_defect import shadow
from verify_rectangular_obstruction import peel
from explore_seed267_endpoints import verify

def main():
 start=time.monotonic();src=ROOT/'papers/dismantlable-ample/evidence/order267_obstruction/record.json';X=set(json.loads(src.read_text())['family']);r=json.loads((HERE/'two_row_verification.json').read_text());vals=sorted(X);bits={x:1<<i for i,x in enumerate(vals)};full=(1<<len(vals))-1;groups=[]
 for s in shadow(X,12):
  g={}
  for x in X:g[x&s]=g.get(x&s,0)|bits[x]
  groups.append(list(g.values()))
 ample=[];unknown=[];nonample=0
 for c in r['certificates']:
  rem=c['removed'];Y=X-set(rem);remain=full^sum(bits[x] for x in rem)
  if sum(all(t&remain for t in g) for g in groups)!=len(Y):nonample+=1;continue
  try:seq=peel(Y,12)
  except AssertionError as err:
   if not err.args or not isinstance(err.args[0],tuple) or err.args[0][0]!='stalled':raise
   unknown.append(rem);continue
  verify(Y,seq,12);ample.append(dict(removed=rem,order=seq))
 report=dict(status='COMPLETE' if not unknown else 'UNKNOWN',checked=len(r['certificates']),nonample=nonample,positive=len(ample),positive_certificates=ample,unknown=unknown,seconds=time.monotonic()-start,source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),scope='Every ample output of the7148 fixed single-row cylinder deletions; cornered outputs included.')
 (HERE/'two_row_fiber_positivity.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k not in ['positive_certificates','unknown']}),flush=True)
if __name__=='__main__':main()
