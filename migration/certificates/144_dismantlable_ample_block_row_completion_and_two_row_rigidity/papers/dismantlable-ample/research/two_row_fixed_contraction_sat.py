"""Exact all-diagonal-assignment search with the208 contraction fixed.
H ampleness is encoded by Sh=strong-Sh. Cornerlessness uses full incident
cube witnesses, not a degree surrogate. UNSAT requires external DRAT replay.
"""
from pathlib import Path
import sys,json,time,threading,subprocess,hashlib
from itertools import combinations
from pysat.formula import IDPool,CNF
from pysat.card import CardEnc
from pysat.solvers import Solver
ROOT=Path(__file__).resolve().parents[3];HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'papers/ample-corners/scripts'))
from verify_full_base_defect import shadow,corners,submasks

def build(D0,D1,bound):
 C=D0&D1;pool=IDPool();cs=[]
 def AND(xs):
  xs=list(xs);g=pool.id();cs.extend([[-g,x] for x in xs]);cs.append([g]+[-x for x in xs]);return g
 def OR(xs):return -AND([-x for x in xs])
 a={(x,b):pool.id(('member',x,b)) for x in sorted(C) for b in (0,1)}
 overlap={x:AND([a[x,0],a[x,1]]) for x in sorted(C)}
 for x in sorted(C):cs.append([a[x,0],a[x,1]])
 # Supports avoiding z are exactly Sh(C), regardless of assignment.
 for s in sorted(shadow(C,10)):
  cubes=[{x^t for t in submasks(s)} for x in sorted(C) if x&s==0 and all(x^t in C for t in submasks(s))]
  strong_old=[AND(a[x,b] for x in sorted(Q)) for b in (0,1) for Q in cubes];cs.append(strong_old)
  sh=[]
  for b in (0,1):
   traces=[OR(a[x,b] for x in sorted(C) if x&s==t) for t in submasks(s)]
   sh.append(AND(traces))
  strong_new=OR(AND(overlap[x] for x in sorted(Q)) for Q in cubes)
  cs.append([-sh[0],-sh[1],strong_new])
 # All vertices of the maximal four-row envelope; constants have variables.
 truth=pool.id();cs.append([truth])
 member={x|1024:truth for x in sorted(D0)};member.update({x|2048:truth for x in sorted(D1)})
 member.update({x:a[x,0] for x in sorted(C)});member.update({x|3072:a[x,1] for x in sorted(C)})
 for v,present in member.items():
  dirs=[i for i in range(12) if v^(1<<i) in member];bad=[]
  for k in range(2,len(dirs)+1):
   for T in combinations(dirs,k):
    target=v^sum(1<<j for j in T);q=member.get(target)
    if q==truth:continue
    terms=[member[v^(1<<j)] for j in T]
    if q is not None:terms.append(-q)
    bad.append(AND(terms))
  cs.append([-present]+bad)
 cs.extend(CardEnc.atmost(list(overlap.values()),bound=bound,vpool=pool).clauses)
 # Normalize only repeated literals and tautologies; preserve satisfiability.
 clean=[]
 for clause in cs:
  clause=list(dict.fromkeys(clause))
  if not any(-lit in clause for lit in clause):clean.append(clause)
 return pool,clean,a

def solve_case(name,D0,D1,bound,out,seconds,forced=None):
 pool,cs,a=build(D0,D1,bound)
 if forced is not None:cs.extend([[lit if x in forced[b] else -lit] for (x,b),lit in a.items()])
 cnf=out/(name+'.cnf');CNF(from_clauses=cs).to_file(str(cnf))
 result=dict(bound=bound,variables=pool.top,clauses=len(cs))
 start=time.monotonic()
 with Solver(name='g4',bootstrap_with=cs,with_proof=True) as solver:
  timer=threading.Timer(seconds,solver.interrupt);timer.start()
  try:ans=solver.solve_limited(expect_interrupt=True)
  finally:timer.cancel()
  if ans is None:result['status']='UNKNOWN'
  elif ans:
   model=set(solver.get_model());B=[{x for (x,b),lit in a.items() if b==j and lit in model} for j in (0,1)]
   V=B[0]|{x|1024 for x in sorted(D0)}|{x|2048 for x in sorted(D1)}|{x|3072 for x in B[1]}
   assert len(shadow(V,12))==len(V) and not corners(V,12)
   assert len(V)<=len(D0)+len(D1)+len(D0&D1)+bound
   result.update(status='SAT-verified',order=len(V),family=sorted(V),B0=sorted(B[0]),B1=sorted(B[1]))
  else:
   proof=out/(name+'.drat');proof.write_text('\n'.join(solver.get_proof())+'\n');log=out/(name+'.check.log')
   with log.open('w') as f:
    checked=subprocess.run(['/home/ec2-user/damp-proof-archive/checkers/drat-trim',str(cnf),str(proof),'-t','30'],stdout=f,stderr=subprocess.STDOUT,timeout=40)
   result.update(status='UNSAT-DRAT-verified' if checked.returncode==0 and 'VERIFIED' in log.read_text() else 'UNSAT-unverified',checker_returncode=checked.returncode,proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest())
 result.update(seconds=time.monotonic()-start,cnf_sha256=hashlib.sha256(cnf.read_bytes()).hexdigest());return result

def main():
 src=ROOT/'papers/dismantlable-ample/evidence/order267_obstruction/record.json';X=set(json.loads(src.read_text())['family'])
 out=HERE/'two_row_fixed_contraction';out.mkdir(exist_ok=True)
 report=dict(source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),pairs=[])
 for e,f in combinations(range(12),2):
  old=[j for j in range(12) if j not in (e,f)]
  norm=lambda x:sum(((x>>j)&1)<<q for q,j in enumerate(old))
  for flip in (0,1):
   rs=[{norm(x) for x in X if 2*((x>>e)&1)+((x>>f)&1)==(t^flip)} for t in range(4)];B0,D0,D1,B1=rs;C=B0|B1
   if not C or not C<=D0&D1:continue
   assert D0&D1==C
   row=dict(pair=[e,f],flip=flip,D0=sorted(D0),D1=sorted(D1),C_order=len(C),cases={})
   allowance=266-len(D0)-len(D1)-len(C)
   for suffix,bound,forced in [('baseline267',allowance+1,[B0,B1]),('at_most266',allowance,None)]:
    name=f'pair_{e}_{f}_{suffix}';result=solve_case(name,D0,D1,bound,out,45,forced);row['cases'][suffix]=result
    print(name,{k:v for k,v in result.items() if k not in ['family','B0','B1']},flush=True)
   report['pairs'].append(row);(out/'report.json').write_text(json.dumps(report,indent=2)+'\n')
 assert len(report['pairs'])==6
if __name__=='__main__':main()
