"""Exhaustive symmetry split of the Q7 degree-four unique-corner target.

The root cube uses bits 0..3. Its fiber at x=1 is a nontrivial ample
Q3 family containing zero. Canonicalize this fiber under permutations of
bits 4..6 only. This changes no root data and assumes no invariant family.
All 63 eligible labelled fibers are partitioned into orbits below.
Bounded case answers are UNKNOWN or uncertified UNSAT, never proofs.
"""
import argparse
import hashlib
from itertools import permutations
import json
from pathlib import Path
import threading
import time
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from pysat.solvers import Solver
from explore_damp_cornerless_ample_sat import Sink, submasks, verify
from explore_unique_corner_asymmetry import build, add_unique_corner


def catalog():
    perms=list(permutations(range(3)))
    def moved(mask,p):
        return sum(1 << sum(1 << p[j] for j in range(3) if v & (1 << j))
                   for v in range(8) if mask & (1 << v))
    eligible=[m for m in range(3,256,2) if verify([v for v in range(8) if m & (1<<v)],3)[0]]
    orbits={}
    for m in eligible:
        rep=min(moved(m,p) for p in perms)
        orbits.setdefault(rep,[]).append(m)
    assert sum(map(len,orbits.values())) == len(eligible) == 63
    assert {x for o in orbits.values() for x in o} == set(eligible)
    return [dict(representative=m, words=[v for v in range(8) if m & (1<<v)],
                 orbit_size=len(o), orbit=o) for m,o in sorted(orbits.items())]


def main():
    p=argparse.ArgumentParser();p.add_argument('--seconds',type=float,default=8)
    p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    cases=catalog();report=dict(catalog=cases, searches=[],
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        encoder_sha256=hashlib.sha256(Path(__file__).with_name('explore_unique_corner_asymmetry.py').read_bytes()).hexdigest(),
        scope='Exhaustive fiber-orbit split of degree-four Q7 unique-corner existence. Learned clauses reused under assumptions. UNSAT answers uncertified; incomplete cases remain UNKNOWN.')
    with Solver(name='g4') as solver:
        pool=IDPool();cnf=Sink(solver);cov=build(7,cnf,pool)
        add_unique_corner(7,4,cnf,cov)
        cnf.extend(CardEnc.atleast(cov,bound=33,vpool=pool,encoding=EncType.seqcounter).clauses)
        for x in range(1,16):
            cnf.append([cov[x|(y<<4)] for y in range(1,8)])
        report.update(variables=pool.top,clauses=cnf.count)
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        for case in cases:
            mask=case['representative']
            assumptions=[cov[1|(y<<4)] if mask & (1<<y) else -cov[1|(y<<4)] for y in range(8)]
            solver.clear_interrupt();start=time.monotonic()
            timer=threading.Timer(args.seconds,solver.interrupt);timer.start()
            try: answer=solver.solve_limited(assumptions=assumptions,expect_interrupt=True)
            finally: timer.cancel()
            row=dict(representative=mask,seconds=args.seconds,elapsed_seconds=time.monotonic()-start,
                     status={True:'SAT',False:'UNSAT-answer-only',None:'UNKNOWN'}[answer])
            if answer:
                model=set(solver.get_model());family=[v for v,x in enumerate(cov) if x in model]
                ample,nsh,nc=verify(family,7);words=set(family)
                corners=[v for v in family if all(v^m in words for m in submasks(sum(1<<j for j in range(7) if v^(1<<j) in words)))]
                assert ample and nc==1 and corners==[0]
                row.update(family=family,order=len(family),shattered_count=nsh,corners=corners,independent_check='PASS')
            report['searches'].append(row)
            args.output.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(row),flush=True)
            if answer:break

if __name__=='__main__':main()
