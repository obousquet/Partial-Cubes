"""Bounded complete-model test for each of the six fixed contractions.
Stops at64 models; coverage only when the blocking CNF has a checked proof.
"""
from pathlib import Path
import json,time,threading,subprocess,hashlib
from pysat.solvers import Solver
from pysat.formula import CNF
from two_row_fixed_contraction_sat import build,ROOT,HERE,shadow,corners

def solve_pair(p,out):
 D0=set(p['D0']);D1=set(p['D1']);pool,cs,a=build(D0,D1,len(D0&D1));models=[];blocks=[];start=time.monotonic();status='UNKNOWN'
 prefix='pair_'+str(p['pair'][0])+'_'+str(p['pair'][1])+'_rigidity'
 with Solver(name='g4',bootstrap_with=cs,with_proof=True) as solver:
  for attempt in range(65):
   timer=threading.Timer(20,solver.interrupt);timer.start()
   try:ans=solver.solve_limited(expect_interrupt=True)
   finally:timer.cancel()
   if ans is None:break
   if not ans:
    cnf=out/(prefix+'.cnf');CNF(from_clauses=cs+blocks).to_file(str(cnf));proof=out/(prefix+'.drat');proof.write_text('\n'.join(solver.get_proof())+'\n');log=out/(prefix+'.check.log')
    with log.open('w') as f:cp=subprocess.run(['/home/ec2-user/damp-proof-archive/checkers/drat-trim',str(cnf),str(proof),'-t','30'],stdout=f,stderr=subprocess.STDOUT,timeout=40)
    assert cp.returncode==0 and 'VERIFIED' in log.read_text();status='COMPLETE-DRAT-verified';break
   if attempt==64:status='MODEL-CAP';break
   m=set(solver.get_model());B=[{x for (x,b),lit in a.items() if b==j and lit in m} for j in (0,1)];V=B[0]|{x|1024 for x in D0}|{x|2048 for x in D1}|{x|3072 for x in B[1]}
   assert len(shadow(V,12))==len(V) and not corners(V,12)
   models.append(dict(B0=sorted(B[0]),B1=sorted(B[1]),order=len(V)))
   clause=[-lit if lit in m else lit for lit in a.values()];blocks.append(clause);solver.add_clause(clause)
 report=dict(status=status,pair=p['pair'],models=models,model_count=len(models),seconds=time.monotonic()-start)
 # Both models must be the original diagonal rows or their exchange.
 baseline=p['cases']['baseline267'];expected={(tuple(baseline['B0']),tuple(baseline['B1'])),(tuple(baseline['B1']),tuple(baseline['B0']))}
 actual={(tuple(m['B0']),tuple(m['B1'])) for m in models}
 report['only_original_and_exchange']=(status=='COMPLETE-DRAT-verified' and actual==expected)
 print(p['pair'],status,len(models),sorted({m['order'] for m in models}),report['only_original_and_exchange'],flush=True)
 return report

def main():
 src=HERE/'two_row_fixed_contraction/report.json';r=json.loads(src.read_text());out=src.parent
 report=dict(source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),pairs=[],scope='All cornerless ample diagonal assignments for each listed fixed contraction only when COMPLETE. Fixed input rigidity, not a global classification.')
 for p in r['pairs']:
  report['pairs'].append(solve_pair(p,out));(out/'rigidity.json').write_text(json.dumps(report,indent=2)+'\n')
if __name__=='__main__':main()
