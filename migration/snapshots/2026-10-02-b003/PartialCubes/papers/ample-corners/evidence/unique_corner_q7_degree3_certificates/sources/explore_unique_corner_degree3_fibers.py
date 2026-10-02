"""Complete proper-unit-fiber split of the Q7 degree-three unique-corner target.

The proper-fiber lemma covers existence by nontrivial proper ample Q4 fibers
containing zero, modulo coordinate permutations. Fix one such fiber over root
unit 1. Every other nonzero root-cube fiber must also be nontrivial.
SAT models get independent Sh/sSh and corner checks. UNSAT answers are only
answers until independently certified. Timeouts remain UNKNOWN.
"""
import argparse
from collections import Counter
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
    # Projection cells: a support T is shattered iff the family meets every cell.
    cells=[[sum(1<<v for v in range(16) if v & T == b) for b in submasks(T)]
           for T in range(16)]
    eligible=[mask for mask in range(3,65535,2)
              if sum(all(mask & cell for cell in part) for part in cells)==mask.bit_count()]
    maps=[[sum(1<<p[j] for j in range(4) if v & (1<<j)) for v in range(16)]
          for p in permutations(range(4))]
    orbits={}
    for mask in eligible:
        words=[v for v in range(16) if mask & (1<<v)]
        rep=min(sum(1<<mp[v] for v in words) for mp in maps)
        orbits.setdefault(rep,[]).append(mask)
    rows=[dict(representative=m,words=[v for v in range(16) if m & (1<<v)],
               orbit_size=len(o),orbit=o) for m,o in sorted(orbits.items())]
    assert {v for x in rows for v in x['orbit']}==set(eligible)
    # Independent direct Sh/sSh check on every orbit representative.
    assert all(verify(x['words'],4)[0] for x in rows)
    return rows


def main():
    p=argparse.ArgumentParser();p.add_argument('--seconds',type=float,default=1)
    p.add_argument('--catalog-only',action='store_true');p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();cases=catalog()
    report=dict(catalog=cases,searches=[],script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        encoder_sha256=hashlib.sha256(Path(__file__).with_name('explore_unique_corner_asymmetry.py').read_bytes()).hexdigest(),
        labelled_fibers=sum(x['orbit_size'] for x in cases),fiber_orbits=len(cases),
        scope='Complete proper-unit Q4 fiber cover for the Q7 degree-three unique-corner target. No invariant-family or VC restriction. Uncertified UNSAT and UNKNOWN remain distinct.')
    args.output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(labelled_fibers=report['labelled_fibers'],fiber_orbits=len(cases))),flush=True)
    if args.catalog_only:return
    with Solver(name='g4') as solver:
        pool=IDPool();cnf=Sink(solver);cov=build(7,cnf,pool);add_unique_corner(7,3,cnf,cov)
        cnf.extend(CardEnc.atleast(cov,bound=33,vpool=pool,encoding=EncType.seqcounter).clauses)
        for x in range(1,8):cnf.append([cov[x|(y<<3)] for y in range(1,16)])
        report.update(variables=pool.top,clauses=cnf.count)
        for index,case in enumerate(cases):
            mask=case['representative'];assumptions=[cov[1|(y<<3)] if mask & (1<<y) else -cov[1|(y<<3)] for y in range(16)]
            solver.clear_interrupt();start=time.monotonic();timer=threading.Timer(args.seconds,solver.interrupt);timer.start()
            try:answer=solver.solve_limited(assumptions=assumptions,expect_interrupt=True)
            finally:timer.cancel()
            row=dict(representative=mask,seconds=args.seconds,elapsed_seconds=time.monotonic()-start,
                     status={True:'SAT',False:'UNSAT-answer-only',None:'UNKNOWN'}[answer])
            if answer:
                model=set(solver.get_model());family=[v for v,x in enumerate(cov) if x in model]
                ample,nsh,nc=verify(family,7);words=set(family)
                corners=[v for v in family if all(v^m in words for m in submasks(sum(1<<j for j in range(7) if v^(1<<j) in words)))]
                assert ample and nc==1 and corners==[0]
                row.update(family=family,order=len(family),shattered_count=nsh,corners=corners,independent_check='PASS')
            report['searches'].append(row)
            args.output.write_text(json.dumps(report,indent=2)+'\n')
            if index%10==9 or answer is True or index==len(cases)-1:
                print(json.dumps(dict(completed=index+1,total=len(cases),statuses=dict(Counter(x['status'] for x in report['searches'])))),flush=True)
            if answer is True:break

if __name__=='__main__':main()
