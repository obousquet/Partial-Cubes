"""Exact cornerless-ample search over ALL proper nonempty subsets of Y267.
No fixed rows, prescribed shadow, degree proxy, or deletion-radius bound.
"""
from pathlib import Path
from itertools import combinations
import json,time,threading,subprocess,hashlib,sys
from pysat.formula import IDPool,CNF
from pysat.solvers import Solver
from verify_full_base_defect import shadow,corners,submasks
ROOT=Path(__file__).resolve().parents[3]

def build(X,n,require_ample=True,require_cornerless=True,protected=()):
 pool=IDPool();clauses=[];p={x:pool.id(('member',x)) for x in sorted(X)}
 def AND(xs):
  xs=list(dict.fromkeys(xs));g=pool.id();clauses.extend([[-g,x] for x in xs]);clauses.append([g]+[-x for x in xs]);return g
 def OR(xs):return -AND([-x for x in xs])
 for S in (sorted(shadow(X,n)) if require_ample else []):
  sh=AND(OR(p[x] for x in sorted(X) if x&S==t) for t in submasks(S))
  cubes=[{x^t for t in submasks(S)} for x in sorted(X) if x&S==0 and all(x^t in X for t in submasks(S))]
  strong=OR(AND(p[x] for x in sorted(Q)) for Q in cubes)
  clauses.append([-sh,strong])
 for v in ([x for x in sorted(X) if x not in protected] if require_cornerless else []):
  dirs=[j for j in range(n) if v^(1<<j) in X];bad=[]
  for k in range(2,len(dirs)+1):
   for T in combinations(dirs,k):
    w=v^sum(1<<j for j in T);terms=[p[v^(1<<j)] for j in T]
    if w in p:terms.append(-p[w])
    bad.append(AND(terms))
  clauses.append([-p[v]]+bad)
 clauses.append(list(p.values())) # nonempty
 return pool,clauses,p

def main():
 start=time.monotonic();source=ROOT/'papers/dismantlable-ample/evidence/order267_obstruction/record.json';X=set(json.loads(source.read_text())['family']);assert len(X)==267 and len(shadow(X,12))==267 and not corners(X,12)
 out=ROOT/'papers/ample-corners/evidence/seed267_all_subfamilies';out.mkdir(exist_ok=True);pool,base,p=build(X,12)
 with Solver(name='g4',bootstrap_with=base) as s:assert s.solve(assumptions=list(p.values()))
 clauses=base+[[-v for v in p.values()]];cnf=out/'proper_subfamily.cnf';CNF(from_clauses=clauses).to_file(str(cnf));print('formula',pool.top,len(clauses),flush=True)
 report=dict(source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),variables=pool.top,clauses=len(clauses),cnf_sha256=hashlib.sha256(cnf.read_bytes()).hexdigest(),scope='All nonempty proper subsets of the recorded267 seed, with exact ampleness and cornerlessness. No statement about classes using other vertices.')
 with Solver(name='g4',bootstrap_with=clauses,with_proof=True) as s:
  timer=threading.Timer(45,s.interrupt);timer.start()
  try:ans=s.solve_limited(expect_interrupt=True)
  finally:timer.cancel()
  if ans is None:report['status']='UNKNOWN'
  elif ans:
   model=set(s.get_model());Y={x for x,v in p.items() if v in model};assert 0<len(Y)<267 and len(shadow(Y,12))==len(Y) and not corners(Y,12);report.update(status='SAT-verified',order=len(Y),family=sorted(Y))
  else:
   proof=out/'proper_subfamily.drat';proof.write_text('\n'.join(s.get_proof())+'\n');log=out/'proper_subfamily.check.log'
   with log.open('w') as f:
    try:cp=subprocess.run(['/home/ec2-user/damp-proof-archive/checkers/drat-trim',str(cnf),str(proof),'-t','150'],stdout=f,stderr=subprocess.STDOUT,timeout=180);rc=cp.returncode
    except subprocess.TimeoutExpired:rc=None
   report.update(status='UNSAT-DRAT-verified' if rc==0 and 'VERIFIED' in log.read_text() else 'UNSAT-unverified',checker_returncode=rc,proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest())
 report['seconds']=time.monotonic()-start;(out/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if k!='family'}),flush=True)
if __name__=='__main__':main()
