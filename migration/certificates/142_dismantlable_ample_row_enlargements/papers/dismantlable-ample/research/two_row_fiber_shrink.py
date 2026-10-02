"""Exhaust all single coordinate-trace-fiber deletions in the large rows
of the six opposite-nested presentations of the267 seed. Not all subsets.
"""
from pathlib import Path
import json,sys,time,hashlib
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'papers/ample-corners/scripts'))
from verify_full_base_defect import corners,shadow
from two_row_probe import patterns

def fibers(D,C,coords):
 vals=sorted(D);full=(1<<len(vals))-1
 keep=sum(1<<i for i,x in enumerate(vals) if x in C)
 states={full}
 for j in coords:
  m=sum(1<<i for i,x in enumerate(vals) if x>>j&1)
  states|={t&m for t in states}|{t&(full^m) for t in states}
 # Includes every cylinder, deduplicated by its intersection with D.
 return [set(vals[i] for i in range(len(vals)) if t>>i&1) for t in sorted(states) if t and not t&keep]

def main():
 start=time.monotonic();src=ROOT/'papers/dismantlable-ample/evidence/order267_obstruction/record.json';X=set(json.loads(src.read_text())['family']);out=[];seen=set();hits=[]
 # A subset cannot shatter a support absent from X. Cache all original
 # trace fibers; a support survives exactly when every fiber survives.
 vals=sorted(X);bits={x:1<<i for i,x in enumerate(vals)};full=(1<<len(vals))-1
 trace_fibers=[]
 for support in shadow(X,12):
  groups={}
  for x in vals:groups[x&support]=groups.get(x&support,0)|bits[x]
  trace_fibers.append(list(groups.values()))
 assert len(trace_fibers)==len(X)
 for pat in patterns(X,12):
  e,f=pat['pair'];flip=pat['flip'];mask=(1<<e)|(1<<f)
  rows=[{x&~mask for x in X if 2*((x>>e)&1)+((x>>f)&1)==(t^flip)} for t in range(4)];C=rows[0]|rows[3]
  for t in (1,2):
   fs=fibers(rows[t],C,[j for j in range(12) if j not in (e,f)])
   row=dict(pair=[e,f],flip=flip,row=t,fibers=len(fs),unique_outputs=0,cornerless=0,ample_cornerless=0,deletion_sizes={})
   for removed in fs:
    actual=t^flip;removedV={x|((actual>>1)<<e)|((actual&1)<<f) for x in removed};key=tuple(sorted(removedV))
    row['deletion_sizes'][len(removed)]=row['deletion_sizes'].get(len(removed),0)+1
    if key in seen:continue
    seen.add(key);row['unique_outputs']+=1;Y=X-removedV
    if corners(Y,12):continue
    row['cornerless']+=1
    # Exact shattering equality, no greedy failure as negative evidence.
    remain=full^sum(bits[x] for x in removedV)
    sh_count=sum(all(t&remain for t in group) for group in trace_fibers)
    if sh_count!=len(Y):continue
    row['ample_cornerless']+=1;hits.append(dict(pair=[e,f],flip=flip,row=t,removed=sorted(removedV),family=sorted(Y),order=len(Y)))
   out.append(row);print(row,flush=True)
 report=dict(status='PASS',source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),rows=out,unique_outputs=len(seen),hits=hits,seconds=time.monotonic()-start,scope='Exhaustive single-cylinder deletions disjoint from C in either large row of the six recorded presentations. No simultaneous two-row deletion or arbitrary subset coverage.')
 (HERE/'two_row_fiber_shrink.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
