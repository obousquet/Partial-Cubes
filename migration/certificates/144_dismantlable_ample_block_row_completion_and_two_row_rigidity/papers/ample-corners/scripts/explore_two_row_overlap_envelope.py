"""All ample subfamilies of the 96-word envelope above one recorded overlap.

Decision: can the marked 267 two-row contraction supply a forbidden reduced
input? Search every nonempty ample cornerless subfamily of C x Q1, not just
lifts projecting onto C. UNSAT implies all ample subfamilies peel, by induction.
Only the coordinate pair (0,1), flip1 is covered; no global size lower bound.
"""
from pathlib import Path
import json, hashlib, time, threading, subprocess, gzip
from pysat.formula import CNF
from pysat.solvers import Solver
from explore_seed267_subfamilies import build, ROOT
from verify_full_base_defect import shadow, corners

OUT=ROOT/'papers/ample-corners/evidence/two_row_overlap_envelope'
SOURCE=ROOT/'papers/dismantlable-ample/evidence/order267_obstruction/record.json'
def inputs():
    X=set(json.loads(SOURCE.read_text())['family'])
    rows=[{x>>2 for x in X if (2*(x&1)+((x>>1)&1))==t} for t in range(4)]
    B0,D0,D1,B1=[rows[t^1] for t in range(4)]
    C=B0|B1
    assert C==D0&D1 and len(C)==48
    return C, C|{x|1024 for x in C}
def main():
    start=time.monotonic();OUT.mkdir(exist_ok=True)
    C,E=inputs();assert len(shadow(C,10))==48 and len(E)==96
    pool,clauses,p=build(E,11)
    cnf=OUT/'search.cnf';CNF(from_clauses=clauses).to_file(str(cnf))
    report=dict(base=sorted(C),envelope=sorted(E),variables=pool.top,clauses=len(clauses),source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),cnf_sha256=hashlib.sha256(cnf.read_bytes()).hexdigest(),scope=__doc__)
    print('formula',pool.top,len(clauses),flush=True)
    with Solver(name='g4',bootstrap_with=clauses,with_proof=True) as solver:
        timer=threading.Timer(30,solver.interrupt);timer.start()
        try:ans=solver.solve_limited(expect_interrupt=True)
        finally:timer.cancel()
        if ans is None:report['status']='UNKNOWN'
        elif ans:
            model=set(solver.get_model());Y={x for x,v in p.items() if v in model}
            assert Y and len(shadow(Y,11))==len(Y) and not corners(Y,11)
            report.update(status='SAT-directly-verified',family=sorted(Y),order=len(Y))
        else:
            proof=OUT/'search.drat';proof.write_text('\n'.join(solver.get_proof())+'\n')
            report['proof_sha256']=hashlib.sha256(proof.read_bytes()).hexdigest()
            log=OUT/'check.log'
            with log.open('w') as f:
                cp=subprocess.run(['/home/ec2-user/damp-proof-archive/checkers/drat-trim',str(cnf),str(proof),'-t','45'],stdout=f,stderr=subprocess.STDOUT,timeout=55)
            report['status']='UNSAT-DRAT-verified' if cp.returncode==0 and 's VERIFIED' in log.read_text() else 'UNSAT-unverified'
            with gzip.open(OUT/'search.drat.gz','wb') as f:f.write(proof.read_bytes())
            proof.unlink()
    report['seconds']=time.monotonic()-start
    (OUT/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('base','envelope','family','scope')}),flush=True)
if __name__=='__main__':main()
