"""SAT search for cornerless ample families on n coordinates (optionally of order <= N).

Structure used (ample C, CCMW Lemma 4.1):
  * for every maximal shattered set S of C the reduction C^S is ample with only the empty set
    shattered, so C contains exactly one S-cube;
  * every maximal cube of C has a maximal shattered direction set (its placement is an isolated
    point of the connected ample reduction C^T, so C^T is a single point);
  * hence C is the union of one cube per facet of the shattered complex Sh(C), and a corner
    (a vertex in a unique maximal cube) is a vertex covered by exactly one of these cubes.
Encoding: choose a down-closed complex S on [n] (all singletons included: irredundant), one cube
per facet of S, and require that the union U shatters no set outside S.  Then
S <= sSh(U) and Sh(U) <= S, and the Sandwich inequalities |sSh| <= |U| <= |Sh| force
|S| = |sSh(U)| = |U| = |Sh(U)|: U is ample with Sh(U) = S.  Conversely every ample family arises
with S = Sh(C).  Cornerless: every covered vertex lies in at least two chosen cubes.

Options: n [--max-order=N] [--no-cornerless] [--lex] [--fix-max-vc=d] [--dimacs=F --cnf-only]
         [--solver=cadical195]
Models are verified independently (ampleness by counting, corners by maximal cubes).
"""
import sys
from itertools import combinations

from pysat.card import CardEnc, EncType
from pysat.formula import IDPool
from pysat.solvers import Solver


class Sink:
    def __init__(self, solver, log=None):
        self.solver, self.count, self.log = solver, 0, log

    def append(self, cl):
        if self.solver is not None:
            self.solver.add_clause(cl)
        self.count += 1
        if self.log is not None:
            self.log.write(" ".join(map(str, cl)) + " 0\n")

    def extend(self, cls):
        for cl in cls:
            self.append(cl)


def submasks(m):
    sub = m
    while True:
        yield sub
        if sub == 0:
            return
        sub = (sub - 1) & m


def build(n, cnf, pool, cornerless=True, max_order=None, lex=False, fix_max_vc=None, vc_exact=None, sh_split=None):
    full = (1 << n) - 1
    subsets = sorted(range(1 << n), key=lambda T: (bin(T).count("1"), T))
    sh = {T: pool.id(("sh", T)) for T in subsets}
    fac = {T: pool.id(("fac", T)) for T in subsets}
    s = lambda T, b: pool.id(("s", T, b))
    # complex: contains empty set and all singletons; down-closed
    cnf.append([sh[0]])
    for i in range(n):
        cnf.append([sh[1 << i]])
    for T in subsets:
        for i in range(n):
            if T >> i & 1:
                cnf.append([-sh[T], sh[T ^ (1 << i)]])
    if fix_max_vc is not None:
        for T in subsets:
            cnf.append([sh[T]] if bin(T).count("1") <= fix_max_vc else [-sh[T]])
    if vc_exact is not None:
        # case split: VC dimension exactly vc_exact (some face of that size, none larger)
        cnf.append([sh[T] for T in subsets if bin(T).count("1") == vc_exact])
        for T in subsets:
            if bin(T).count("1") == vc_exact + 1:
                cnf.append([-sh[T]])
    if sh_split is not None:
        # case split: signs of sh[T] for the last k sets of size vc_exact in the order (bits of idx)
        k, idx = sh_split
        level = [T for T in subsets if bin(T).count("1") == vc_exact][-k:]
        for j, T in enumerate(level):
            cnf.append([sh[T]] if idx >> j & 1 else [-sh[T]])
    # facets
    for T in subsets:
        ups = [T | (1 << i) for i in range(n) if not T >> i & 1]
        cnf.append([-fac[T], sh[T]])
        for U in ups:
            cnf.append([-fac[T], -sh[U]])
        cnf.append([-sh[T]] + [sh[U] for U in ups] + [fac[T]])
    # one cube per facet, none otherwise
    for T in subsets:
        comp = full ^ T
        lits = [s(T, b) for b in submasks(comp)]
        cnf.append([-fac[T]] + lits)
        for x in lits:
            cnf.append([-x, fac[T]])
        cnf.extend(CardEnc.atmost(lits=lits, bound=1, vpool=pool, encoding=EncType.seqcounter).clauses)
    # coverage (and cornerless double coverage)
    cov = []
    for v in range(1 << n):
        cands = [s(T, v & (full ^ T)) for T in subsets]
        c = pool.id(("c", v)); cov.append(c)
        for x in cands:
            cnf.append([-x, c])
        cnf.append([-c] + cands)
        if cornerless:
            prev1 = prev2 = None
            for x in cands:
                r1, r2 = pool._next(), pool._next()
                cnf.append([-r1, x] + ([prev1] if prev1 else []))
                if prev1 is None:
                    cnf.append([-r2])
                else:
                    cnf.append([-r2, prev2, prev1] if prev2 else [-r2, prev1])
                    cnf.append([-r2, prev2, x] if prev2 else [-r2, x])
                prev1, prev2 = r1, r2
            cnf.append([-c, prev2])
    # U shatters no set outside the complex
    for T in subsets:
        if T == 0:
            continue
        ys = [pool.id(("y", T, p)) for p in submasks(T)]
        cnf.append([sh[T]] + ys)
        for v in range(1 << n):
            cnf.append([-pool.id(("y", T, v & T)), -cov[v]])
    if max_order is not None:
        cnf.extend(CardEnc.atmost(lits=cov, bound=max_order, vpool=pool, encoding=EncType.totalizer).clauses)
    if lex:
        order = [(T, b) for T in subsets for b in submasks(full ^ T)]
        def tr(i):
            p_ = list(range(n)); p_[i], p_[i + 1] = p_[i + 1], p_[i]
            pm = lambda m: sum(1 << p_[k] for k in range(n) if m >> k & 1)
            return lambda T, b: (pm(T), pm(b))
        def fl(j):
            return lambda T, b: (T, (b ^ (1 << j)) if not (T >> j & 1) else b)
        for f in [tr(i) for i in range(n - 1)] + [fl(j) for j in range(n)]:
            prev = None
            for (T, b) in order:
                x = s(T, b); gx = s(*f(T, b))
                if x == gx:
                    continue
                cnf.append(([-prev] if prev else []) + [-x, gx])
                e = pool._next()
                cnf.append(([-prev] if prev else []) + [-x, -gx, e])
                cnf.append(([-prev] if prev else []) + [x, gx, e])
                prev = e
    return cov


def verify(U, n):
    """Independent check: ampleness by counting (|Sh| = |sSh| = |U|) and corner count."""
    S = set(U)
    sh = [T for T in range(1 << n) if len({x & T for x in U}) == 1 << bin(T).count("1")]
    def strongly(T):
        for x in U:
            if x & T:
                continue
            if all((x | sub) in S for sub in submasks(T)):
                return True
        return False
    ssh = [T for T in sh if strongly(T)]
    ample = len(sh) == len(ssh) == len(U)
    corners = 0
    for v in U:
        dirs = [i for i in range(n) if (v ^ (1 << i)) in S]
        maxes = []
        for k in range(len(dirs), -1, -1):
            for D in combinations(dirs, k):
                m = sum(1 << i for i in D)
                if any((m | E) == E for E in maxes):
                    continue
                if all((v ^ sub) in S for sub in submasks(m)):
                    maxes.append(m)
        corners += len(maxes) == 1
    return ample, len(sh), corners


def decode_and_verify(n, map_file, model_file):
    """Decode a solver model (v-lines) into the family U and verify it independently."""
    import json
    m = json.load(open(map_file))
    true = set()
    for line in open(model_file):
        if line.startswith("v"):
            true.update(int(t) for t in line.split()[1:] if int(t) > 0)
    U = [v for v in range(1 << n) if m["cov"][v] in true]
    ample, nsh, corners = verify(U, n)
    return U, ample, nsh, corners


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "decode":
    n, map_file, model_file = int(sys.argv[2]), sys.argv[3], sys.argv[4]
    U, ample, nsh, corners = decode_and_verify(n, map_file, model_file)
    print(f"decoded family: order {len(U)}, ample {ample}, |Sh| {nsh}, corners {corners}")
    print("family:", U)
    sys.exit(0)

if __name__ == "__main__":
    n = int(sys.argv[1])
    arg = lambda k, d=None: next((a.split("=", 1)[1] for a in sys.argv if a.startswith(k + "=")), d)
    max_order = int(arg("--max-order")) if arg("--max-order") else None
    fix = int(arg("--fix-max-vc")) if arg("--fix-max-vc") else None
    vce = int(arg("--vc-exact")) if arg("--vc-exact") else None
    shs = tuple(map(int, arg("--sh-split").split(":"))) if arg("--sh-split") else None
    dimacs = arg("--dimacs")
    cornerless = "--no-cornerless" not in sys.argv
    lex = "--lex" in sys.argv
    pool = IDPool()
    if dimacs and "--cnf-only" in sys.argv:
        body = open(dimacs + ".body", "w")
        sink = Sink(None, body)
        cov = build(n, sink, pool, cornerless, max_order, lex, fix, vce, shs)
        body.close()
        if arg("--map"):
            import json
            json.dump({"n": n, "cov": cov}, open(arg("--map"), "w"))
        with open(dimacs, "w") as out:
            out.write(f"p cnf {pool.top} {sink.count}\n")
            for line in open(dimacs + ".body"):
                out.write(line)
        import os; os.remove(dimacs + ".body")
        print(f"n={n} max_order={max_order} cornerless={cornerless} lex={lex} fix_max_vc={fix}: "
              f"CNF written, vars {pool.top}, clauses {sink.count}", flush=True)
        sys.exit(0)
    with Solver(name=arg("--solver", "cadical195")) as S:
        sink = Sink(S)
        cov = build(n, sink, pool, cornerless, max_order, lex, fix, vce, shs)
        print(f"n={n} max_order={max_order} cornerless={cornerless} lex={lex} fix_max_vc={fix}: "
              f"vars {pool.top}, clauses {sink.count}", flush=True)
        sat = S.solve()
        print("SAT" if sat else "UNSAT", flush=True)
        if sat:
            model = set(l for l in S.get_model() if l > 0)
            U = [v for v in range(1 << n) if cov[v] in model]
            ample, nsh, corners = verify(U, n)
            print(f"verified: order {len(U)}, ample {ample}, |Sh| {nsh}, corners {corners}", flush=True)
            print("family:", U)
