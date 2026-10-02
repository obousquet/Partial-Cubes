"""Exact vertex encoding via Lawrence (1983), Theorems 3--4.

A family is ample iff each antipodally invariant face section is empty or full.
Auxiliary variables record XOR of membership at opposite face vertices.
For each face, if all its XORs vanish, all pair representatives must agree.
No restriction on VC dimension, symmetry, or shattered complex is imposed.
A Q7 unique-corner model is checked independently by Sh/sSh counts.
UNSAT answers are NOT certificates. UNKNOWN is not an exclusion.
"""
import argparse
import hashlib
import json
from pathlib import Path
import threading
import time
from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from pysat.solvers import Solver
from explore_damp_cornerless_ample_sat import Sink, submasks, verify


def build(n, cnf, pool):
    cov = [pool.id(('member', v)) for v in range(1 << n)]
    dif = {}
    for a in range(1 << n):
        for b in range(a + 1, 1 << n):
            if (a ^ b).bit_count() < 2:
                continue
            z = pool.id(('xor', a, b)); dif[a, b] = z
            x, y = cov[a], cov[b]
            cnf.extend([[-z, x, y], [-z, -x, -y], [z, -x, y], [z, x, -y]])
    full = (1 << n) - 1
    for T in range(1 << n):
        if T.bit_count() < 2:
            continue  # All sections of a zero/one-dimensional face are ample.
        for b in submasks(full ^ T):
            pairs = [(b | m, (b | m) ^ T) for m in submasks(T)
                     if (b | m) < ((b | m) ^ T)]
            unequal = [dif[a, c] for a, c in pairs]
            ref = cov[pairs[0][0]]
            for a, _ in pairs[1:]:
                cnf.append(unequal + [-ref, cov[a]])
                cnf.append(unequal + [ref, -cov[a]])
    return cov


def add_unique_corner(n, d, cnf, cov):
    support = (1 << d) - 1
    for v in submasks(support):
        cnf.append([cov[v]])
    for j in range(d, n):
        cnf.append([-cov[1 << j]])
    for v in range(1, 1 << n):
        if not v & support:
            cnf.append([-cov[v]])  # Connected root-zero section is a singleton.
        for T in range(1 << n):
            cnf.append([-cov[v]] + [cov[v ^ (1 << j)] for j in range(n)
                                   if not T & (1 << j)] +
                       [-cov[v ^ m] for m in submasks(T) if m])
    # Every coordinate is active. In an ample family this also gives an edge.
    for j in range(n):
        cnf.append([cov[v] for v in range(1 << n) if v & (1 << j)])


def controls():
    checked = 0; ample_count = 0
    with Solver(name='g4') as solver:
        pool = IDPool(); cnf = Sink(solver); cov = build(3, cnf, pool)
        for bits in range(256):
            family = [v for v in range(8) if bits & (1 << v)]
            actual = solver.solve(assumptions=[x if bits & (1 << v) else -x
                                               for v, x in enumerate(cov)])
            expected = verify(family, 3)[0]
            assert actual == expected, (bits, actual, expected)
            checked += 1; ample_count += expected
    with Solver(name='g4') as solver:
        pool = IDPool(); cnf = Sink(solver); cov = build(3, cnf, pool)
        add_unique_corner(3, 1, cnf, cov)
        assert not solver.solve()
    return dict(all_Q3_families_checked=checked, ample_count=ample_count,
                unique_corner_Q3_degree_one='UNSAT-answer-only', status='PASS')


def run(d, seconds):
    start = time.monotonic(); n = 7; pool = IDPool()
    with Solver(name='g4') as solver:
        cnf = Sink(solver); cov = build(n, cnf, pool)
        add_unique_corner(n, d, cnf, cov)
        cnf.extend(CardEnc.atleast(cov, bound=33, vpool=pool,
                                  encoding=EncType.seqcounter).clauses)
        timer = threading.Timer(seconds, solver.interrupt); timer.start()
        try:
            answer = solver.solve_limited(expect_interrupt=True)
        finally:
            timer.cancel()
        result = dict(n=n, root_degree=d, seconds=seconds,
                      status={True:'SAT', False:'UNSAT-answer-only', None:'UNKNOWN'}[answer],
                      elapsed_seconds=time.monotonic()-start, variables=pool.top,
                      clauses=cnf.count, statistics=solver.accum_stats())
        if answer:
            model = set(solver.get_model()); family = [v for v, x in enumerate(cov) if x in model]
            ample, nsh, nc = verify(family, n)
            words = set(family)
            corners = [v for v in family if all(v ^ m in words for m in submasks(
                sum(1 << j for j in range(n) if v ^ (1 << j) in words)))]
            assert ample and nc == 1 and corners == [0]
            result.update(family=family, order=len(family), shattered_count=nsh,
                          corners=corners, independent_check='PASS')
        return result


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--seconds', type=float, default=60)
    parser.add_argument('--degrees', default='3,4')
    parser.add_argument('--output', type=Path, required=True)
    args=parser.parse_args()
    report=dict(controls=controls(), searches=[],
                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                verifier_sha256=hashlib.sha256(Path(__file__).with_name(
                    'explore_damp_cornerless_ample_sat.py').read_bytes()).hexdigest(),
                source='Lawrence 1983, Theorems 3--4, shared library math/lawrence1983lopsided',
                scope='Exact Q7 unique-corner degree cases; no VC or family-symmetry restriction. UNSAT answers uncertified.')
    for d in map(int,args.degrees.split(',')):
        result=run(d,args.seconds); report['searches'].append(result)
        args.output.write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(result),flush=True)
        if result['status']=='SAT': break

if __name__=='__main__': main()
