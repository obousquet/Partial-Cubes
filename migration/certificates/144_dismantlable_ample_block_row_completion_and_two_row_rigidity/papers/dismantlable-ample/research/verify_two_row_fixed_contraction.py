"""Reconstruct inputs and blocked-model CNFs, then independently replay DRAT.
The manuscript proves encoder coverage; DRAT checks propositional refutations.
"""
from pathlib import Path
import json,hashlib,subprocess
from pysat.formula import CNF
from two_row_fixed_contraction_sat import build,ROOT,HERE,shadow,corners

def main():
 out=HERE/'two_row_fixed_contraction';base=json.loads((out/'report.json').read_text());rig=json.loads((out/'rigidity.json').read_text());src=ROOT/'papers/dismantlable-ample/evidence/order267_obstruction/record.json';assert hashlib.sha256(src.read_bytes()).hexdigest()==base['source_sha256'];assert hashlib.sha256((out/'report.json').read_bytes()).hexdigest()==rig['source_sha256'];X=set(json.loads(src.read_text())['family']);checks=[]
 assert {tuple(p['pair']) for p in base['pairs']}=={(j,j+1) for j in range(0,12,2)}
 for p,r in zip(base['pairs'],rig['pairs']):
  assert p['pair']==r['pair'];e,f=p['pair'];old=[j for j in range(12) if j not in (e,f)];norm=lambda x:sum(((x>>j)&1)<<q for q,j in enumerate(old))
  rs=[{norm(x) for x in X if 2*((x>>e)&1)+((x>>f)&1)==(t^p['flip'])} for t in range(4)];B0,D0,D1,B1=rs;assert D0==set(p['D0']) and D1==set(p['D1']);assert B0|B1==D0&D1
  _,clauses,a=build(D0,D1,len(D0&D1));expected={(tuple(sorted(B0)),tuple(sorted(B1))),(tuple(sorted(B1)),tuple(sorted(B0)))};assert r['model_count']==2 and {(tuple(m['B0']),tuple(m['B1'])) for m in r['models']}==expected
  for m in r['models']:
   Bs=[set(m['B0']),set(m['B1'])];V=Bs[0]|{x|1024 for x in D0}|{x|2048 for x in D1}|{x|3072 for x in Bs[1]};assert len(V)==267 and len(shadow(V,12))==267 and not corners(V,12)
   clauses.append([-lit if x in Bs[b] else lit for (x,b),lit in a.items()])
  prefix='pair_'+str(e)+'_'+str(f)+'_rigidity';cnf=out/(prefix+'.cnf');proof=out/(prefix+'.drat');assert CNF(from_file=str(cnf)).clauses==clauses
  log=out/(prefix+'.replay.log')
  with log.open('w') as stream:cp=subprocess.run(['/home/ec2-user/damp-proof-archive/checkers/drat-trim',str(cnf),str(proof),'-t','30'],stdout=stream,stderr=subprocess.STDOUT,timeout=40)
  assert cp.returncode==0 and 'VERIFIED' in log.read_text()
  checks.append(dict(pair=[e,f],models=2,cnf_sha256=hashlib.sha256(cnf.read_bytes()).hexdigest(),proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest()))
 report=dict(status='PASS',pairs=checks,checker_sha256=hashlib.sha256(Path('/home/ec2-user/damp-proof-archive/checkers/drat-trim').read_bytes()).hexdigest(),scope='Six fixed contractions, all diagonal assignments. Independent DRAT replay plus direct surviving-model checks; symbolic coverage proof in manuscript.')
 (out/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print('PASS: six full rigidity refutations, twelve directly checked models',flush=True)
if __name__=='__main__':main()
