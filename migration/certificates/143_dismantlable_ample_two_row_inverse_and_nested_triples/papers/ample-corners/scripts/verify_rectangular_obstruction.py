"""Verify the rectangular-boundary obstruction from the explicit 42/22 pair.

Checks direct shattering, complete corner absence, and explicit corner orders
for all 36 elementary pc-minors. Coordinates are kept as 12-bit words with
the removed bit fixed to zero; this is isomorphic to deleting that coordinate.
"""
import hashlib
import json
from pathlib import Path
from verify_full_base_defect import corners, shadow, submasks


def corner(v,F,n):
    T=sum(1<<j for j in range(n) if v^(1<<j) in F)
    return all(v^m in F for m in submasks(T))


def peel(F,n,retain=()):
    F=set(F);retain=set(retain);order=[]
    while F!=retain:
        v=next((v for v in sorted(F-retain) if corner(v,F,n)),None)
        assert v is not None,('stalled',len(F),len(retain))
        order.append(v);F.remove(v)
    return order


def image(F,j,b):
    return {v&~(1<<j) for v in F if b is None or (v>>j)&1==b}


def main():
    root=Path(__file__).parent.parent
    source=root/'evidence/relative_pair_q6/manifest.json';m=json.loads(source.read_text())
    P=set(m['outer']);O=set(m['inner']);H={x|(y<<6) for x in P for y in P if x in O or y in O}
    assert len(H)==1364
    sh=shadow(H,12);assert len(sh)==len(H) and max(s.bit_count() for s in sh)==5
    assert not any(corner(v,H,12) for v in H)
    porder=peel(P,6);oorder=peel(O,6);certs=[]
    for j in range(12):
      for b in [0,1,None]:
        k=j%6;A=image(P,k,b);K=image(O,k,b);rel=peel(A,6,K);ko=peel(K,6)
        combine=(lambda x,y:x|(y<<6)) if j<6 else (lambda x,y:y|(x<<6))
        order=[combine(x,y) for x in rel for y in oorder]+[combine(x,y) for x in ko for y in porder]
        minor=image(H,j,b);remaining=set(minor)
        assert len(order)==len(set(order))==len(minor) and set(order)==minor
        for v in order:
            assert v in remaining and corner(v,remaining,12)
            remaining.remove(v)
        assert not remaining
        certs.append(dict(coordinate=j,operation='contract' if b is None else f'restrict-{b}',order=len(minor),peeling=order,relative_input_order=rel))
    out=root/'evidence/rectangular_boundary_1364';out.mkdir(exist_ok=True)
    report=dict(status='PASS',order=len(H),dimension=12,vc_dimension=5,shattered_supports=len(sh),corners=0,
                elementary_minors_checked=len(certs),family=sorted(H),minor_certificates=certs,
                input_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (out/'certificate.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ['family','minor_certificates']}))

if __name__=='__main__':main()
