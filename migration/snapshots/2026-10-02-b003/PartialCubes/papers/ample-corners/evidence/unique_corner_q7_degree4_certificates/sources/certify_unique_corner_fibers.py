"""Regenerate each Q7 degree-four fiber case and independently check its DRAT proof.

No solver assumptions or inter-case learned clauses are used. Files and checker
fingerprints are retained. Interrupted cases remain unresolved.
"""
import argparse
import hashlib
import json
import subprocess
import threading
import time
from pathlib import Path
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from pysat.solvers import Solver
from explore_damp_cornerless_ample_sat import Sink
from explore_unique_corner_asymmetry import build, add_unique_corner
from explore_unique_corner_fibers import catalog


def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser();p.add_argument('--seconds',type=float,default=60)
    p.add_argument('--checker',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    report=dict(cases=[],checker_sha256=digest(args.checker),script_sha256=digest(Path(__file__)),
        dependencies={n:digest(Path(__file__).with_name(n)) for n in ['explore_unique_corner_fibers.py','explore_unique_corner_asymmetry.py','explore_damp_cornerless_ample_sat.py']},
        scope='Complete Q7 degree-four fiber orbit cover. Each case solved fresh with unit constraints and checked by drat-trim.')
    for case in catalog():
        rep=case['representative'];prefix=args.output/str(rep);pool=IDPool();clauses=[]
        class Collect:
            def append(self,cl):clauses.append(cl)
            def extend(self,cls):clauses.extend(cls)
        cnf=Collect();cov=build(7,cnf,pool);add_unique_corner(7,4,cnf,cov)
        cnf.extend(CardEnc.atleast(cov,bound=33,vpool=pool,encoding=EncType.seqcounter).clauses)
        for x in range(1,16):cnf.append([cov[x|(y<<4)] for y in range(1,8)])
        for y in range(8):cnf.append([cov[1|(y<<4)] if rep & (1<<y) else -cov[1|(y<<4)]])
        formula=prefix.with_suffix('.cnf');proof=prefix.with_suffix('.drat');log=prefix.with_suffix('.check.log')
        with formula.open('w') as f:
            f.write(f'p cnf {pool.top} {len(clauses)}\n')
            for cl in clauses:f.write(' '.join(map(str,cl))+' 0\n')
        start=time.monotonic()
        with Solver(name='g4',bootstrap_with=clauses,with_proof=True) as solver:
            timer=threading.Timer(args.seconds,solver.interrupt);timer.start()
            try: answer=solver.solve_limited(expect_interrupt=True)
            finally:timer.cancel()
            row=dict(representative=rep,seconds=args.seconds,solve_seconds=time.monotonic()-start,
                     status={True:'SAT-unexpected',False:'UNSAT-answer-only',None:'UNKNOWN'}[answer],
                     cnf_sha256=digest(formula),variables=pool.top,clauses=len(clauses))
            if answer is False:
                proof.write_text('\n'.join(solver.get_proof())+'\n');row['proof_sha256']=digest(proof)
                start=time.monotonic()
                with log.open('w') as f:
                    checked=subprocess.run([str(args.checker),str(formula),str(proof),'-t','60'],stdout=f,stderr=subprocess.STDOUT,timeout=70)
                row.update(check_seconds=time.monotonic()-start,checker_returncode=checked.returncode,check_log_sha256=digest(log))
                if checked.returncode==0 and 's VERIFIED' in log.read_text():row['status']='UNSAT-DRAT-checked'
        report['cases'].append(row);(args.output/'manifest.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(row),flush=True)
        if answer is True:break

if __name__=='__main__':main()
