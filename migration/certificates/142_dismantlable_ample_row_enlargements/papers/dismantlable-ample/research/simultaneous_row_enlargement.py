"""Audit the common-completion closed form for iterated row enlargement.

Old bit i becomes f_i; fresh bit n+position(i) is e_i. Constant pairs
encode the old bit. Nonconstant pairs forget it; if any pair is01,
read the partial old word in D, otherwise in H.
"""
from pathlib import Path
from itertools import combinations, permutations
import hashlib, json, sys, time
ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / 'papers/ample-corners/scripts'))
from verify_full_base_defect import shadow, submasks


def project(X, S):
    return {x & ~S for x in X}


def closed(H, D, n, K, positions=None):
    positions = positions or {i:n+j for j,i in enumerate(K)}
    cache = {}
    for S in submasks(sum(1 << i for i in K)):
        cache[S] = (project(H,S), project(D,S))
    out = set()
    for old in range(1 << n):
        for code in range(1 << len(K)):
            word = old; S = 0; enlarged = False
            for j,i in enumerate(K):
                e = (code >> j) & 1; f = (old >> i) & 1
                word |= e << positions[i]
                if e != f:
                    S |= 1 << i
                    enlarged |= e == 0 and f == 1
            if old & ~S in cache[S][int(enlarged)]:
                out.add(word)
    return out


def iterate(H, D, n, K, order):
    positions = {i:n+j for j,i in enumerate(K)}
    current = set(H); done = []
    for i in order:
        # The supplied side is W_done(pi_i D), not a new free choice.
        P = project(D, 1 << i)
        side = closed(P, P, n, done, positions)
        lo = {x for x in current if not x >> i & 1}
        hi = {x & ~(1 << i) for x in current if x >> i & 1}
        union = lo | hi
        assert union <= side
        current = lo | {x | (1 << i) for x in side} | {
            x | (1 << positions[i]) for x in union} | {
            x | (1 << i) | (1 << positions[i]) for x in hi}
        done.append(i)
    return current


def reduce_bits(X, bits):
    out = set(X)
    for i in bits:
        out = {x for x in out if not x >> i & 1 and x | (1 << i) in out}
    return out


def main():
    start = time.monotonic(); checked = order_checks = 0
    for n in (1,2,3):
        ample = []
        for bits in range(1, 1 << (1 << n)):
            X = {x for x in range(1 << n) if bits >> x & 1}
            if len(shadow(X,n)) == len(X):
                ample.append(X)
        for D in ample:
            for H in ample:
                if not H <= D:
                    continue
                for km in range(1, 1 << n):
                    K = [i for i in range(n) if km >> i & 1]
                    R = closed(H,D,n,K)
                    assert len(shadow(R,n+len(K))) == len(R)
                    expected = sum(len(project(H,S)) + ((1 << S.bit_count())-1)*len(project(D,S))
                                   for S in submasks(km))
                    assert len(R) == expected
                    assert reduce_bits(R, range(n,n+len(K))) == H
                    for j,i in enumerate(K):
                        row = {x & ~(1 << i) for x in R if not x >> (n+j) & 1 and x >> i & 1}
                        recovered = reduce_bits(row, [n+k for k in range(len(K)) if k != j])
                        assert recovered == project(D,1 << i)
                    for order in permutations(K):
                        assert iterate(H,D,n,K,order) == R
                        order_checks += 1
                    checked += 1
    # The maximal family with the same selected projections need not be ample.
    path = {0,1,3,7,15}; D = set(range(16)) - path; K = [0,3]
    hull = {x for x in range(16) if all(x & ~(1 << i) in project(D,1 << i) for i in K)}
    assert len(D) == len(shadow(D,4)) == 11
    assert hull == D | {3} and len(shadow(hull,4)) != len(hull)
    hull_example = dict(completion=sorted(D), coordinates=K,
                        hull=sorted(hull), hull_shadow=len(shadow(hull,4)))
    source = ROOT / 'papers/dismantlable-ample/evidence/order267_obstruction/record.json'
    H = set(json.loads(source.read_text())['family']); D = set(range(4096)); K = [0,1]
    R = closed(H,D,12,K)
    assert len(R) == 8011 and reduce_bits(R,[12,13]) == H
    assert iterate(H,D,12,K,[0,1]) == R == iterate(H,D,12,K,[1,0])
    for j,i in enumerate(K):
        row = {x & ~(1 << i) for x in R if not x >> (12+j) & 1 and x >> i & 1}
        assert reduce_bits(row,[12+(1-j)]) == project(D,1 << i)
    report = dict(status='PASS', nested_ample_inputs_and_block_choices=checked,
                  iteration_orders_checked=order_checks, seconds=time.monotonic()-start,
                  nonample_projection_hull=hull_example,
                  forbidden_seed_example=dict(seed_order=267, output_order=8011, dimension=14,
                                              coordinates=K, completion='Q12',
                                              source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                                              output_sha256=hashlib.sha256(json.dumps(sorted(R),separators=(',',':')).encode()).hexdigest(),
                                              scope='Exact formula, both iteration orders and marked recovery checked; forbiddenness follows from the symbolic theorem.'),
                  scope='All nonempty nested ample H subset D on at most3 coordinates, '
                        'all nonempty selected sets, all iteration orders; symbolic theorem supplies general forbiddenness.')
    target = Path(__file__).with_suffix('.json')
    target.write_text(json.dumps(report,indent=2)+'\n'); print(json.dumps(report),flush=True)


if __name__ == '__main__':
    main()
