"""Test explicit circulant families as forbidden pc-minors of R_per.

R_per = convex-treelike (cover-recursively peripheral) partial cubes.
O is a forbidden pc-minor iff O is a partial cube, periphery-free, and every
elementary pc-minor (restriction to a semicube, contraction) lies in R_per.

Families are sets of ints (bitmasks) on n coordinates.
"""
import sys
from functools import lru_cache
from itertools import combinations
from collections import deque


def normalize(fam, n):
    """Drop constant coordinates, translate so that 0 is in, return (frozenset, dim)."""
    fam = list(fam)
    active = [i for i in range(n) if len({(x >> i) & 1 for x in fam}) == 2]
    out = set()
    for x in fam:
        y = 0
        for j, i in enumerate(active):
            if (x >> i) & 1:
                y |= 1 << j
        out.add(y)
    return frozenset(out), len(active)


def section(fam, e, b):
    """Projected section: members with bit e == b, with bit e removed."""
    lo = (1 << e) - 1
    out = set()
    for x in fam:
        if ((x >> e) & 1) == b:
            out.add((x & lo) | ((x >> (e + 1)) << e))
    return frozenset(out)


def contraction(fam, e):
    lo = (1 << e) - 1
    return frozenset((x & lo) | ((x >> (e + 1)) << e) for x in fam)


@lru_cache(maxsize=None)
def rper(fam, n):
    """fam normalized (all n coords active)."""
    if len(fam) <= 1:
        return True
    for e in range(n):
        s0, s1 = section(fam, e, 0), section(fam, e, 1)
        for small, big in ((s0, s1), (s1, s0)):
            if small <= big:
                a = normalize(big, n - 1)
                b = normalize(small, n - 1)
                if rper(*a) and rper(*b):
                    return True
    return False


def peripheral_coords(fam, n):
    res = []
    for e in range(n):
        s0, s1 = section(fam, e, 0), section(fam, e, 1)
        if s0 <= s1:
            res.append((e, 0))
        if s1 <= s0:
            res.append((e, 1))
    return res


def is_partial_cube(fam, n):
    """Isometric in Q_n: BFS distance == Hamming distance for all pairs."""
    fam = list(fam)
    S = set(fam)
    for s in fam:
        dist = {s: 0}
        dq = deque([s])
        while dq:
            x = dq.popleft()
            for i in range(n):
                y = x ^ (1 << i)
                if y in S and y not in dist:
                    dist[y] = dist[x] + 1
                    dq.append(y)
        if len(dist) != len(S):
            return False
        for t, d in dist.items():
            if bin(s ^ t).count("1") != d:
                return False
    return True


def shattered(fam, n, maxk=None):
    """Return list of shattered coordinate sets (as frozensets)."""
    res = []
    maxk = n if maxk is None else maxk
    for k in range(0, maxk + 1):
        found = False
        for T in combinations(range(n), k):
            pats = {tuple((x >> i) & 1 for i in T) for x in fam}
            if len(pats) == 1 << k:
                res.append(frozenset(T))
                found = True
        if not found:
            break
    return res


def strongly_shattered_count(fam, n):
    """Number of X strongly shattered: exists subcube on X in fam."""
    S = set(fam)
    cnt = 0
    for k in range(0, n + 1):
        found = False
        for T in combinations(range(n), k):
            mask = sum(1 << i for i in T)
            # base points: members with zeros on T, check full face
            ok = False
            for x in fam:
                if x & mask:
                    continue
                sub = mask
                good = True
                while True:
                    if (x | sub) not in S:
                        good = False
                        break
                    if sub == 0:
                        break
                    sub = (sub - 1) & mask
                if good:
                    ok = True
                    break
            if ok:
                cnt += 1
                found = True
        if not found:
            break
    return cnt


def is_ample(fam, n):
    return len(shattered(fam, n)) == strongly_shattered_count(fam, n)


def vc(fam, n):
    return max(len(T) for T in shattered(fam, n))


def obstruction_report(fam, n):
    fam, n = normalize(fam, n)
    rep = {"order": len(fam), "dim": n}
    rep["pc"] = is_partial_cube(fam, n)
    if not rep["pc"]:
        return rep
    per = peripheral_coords(fam, n)
    rep["periph"] = per
    if per:
        rep["obs"] = False
        return rep
    bad = []
    for e in range(n):
        for b in (0, 1):
            m = normalize(section(fam, e, b), n - 1)
            if not rper(*m):
                bad.append(("r", e, b))
        m = normalize(contraction(fam, e), n - 1)
        if not rper(*m):
            bad.append(("c", e))
    rep["bad_minors"] = bad
    rep["obs"] = not bad
    return rep


def from_sets(sets):
    return frozenset(sum(1 << i for i in s) for s in sets)


def cyclic_intervals(n, k, include_empty=True, include_full=False):
    fam = set()
    if include_empty:
        fam.add(frozenset())
    for L in range(1, k + 1):
        for a in range(n):
            fam.add(frozenset((a + t) % n for t in range(L)))
    if include_full:
        fam.add(frozenset(range(n)))
    return from_sets(fam)


if __name__ == "__main__":
    # sanity: P5 = cyclic intervals of Z5 of length <= 3
    P5 = cyclic_intervals(5, 3)
    print("P5", obstruction_report(P5, 5), "ample", is_ample(P5, 5))
    Qm = frozenset(x for x in range(8) if x not in (0, 7))
    print("Q3--", obstruction_report(Qm, 3))
    sys.setrecursionlimit(10000)
    print("\ncyclic intervals I(n,k) = {empty} + arcs of length <= k in Z_n")
    print(f"{'n':>3}{'k':>3}{'ord':>6}{'pc':>4}{'amp':>5}{'vc':>4}{'obs':>5}  info")
    for n in range(4, int(sys.argv[1]) if len(sys.argv) > 1 else 10):
        for k in range(1, n):
            fam = cyclic_intervals(n, k)
            rep = obstruction_report(fam, n)
            amp = is_ample(fam, n) if rep["pc"] else None
            v = vc(fam, n)
            info = ""
            if rep["pc"]:
                info = "periph=%d" % len(rep["periph"]) if rep["periph"] else "bad=%s" % rep["bad_minors"][:3]
            print(f"{n:>3}{k:>3}{rep['order']:>6}{str(rep['pc'])[0]:>4}{str(amp)[0]:>5}{v:>4}{str(rep.get('obs'))[0]:>5}  {info}")
            sys.stdout.flush()
