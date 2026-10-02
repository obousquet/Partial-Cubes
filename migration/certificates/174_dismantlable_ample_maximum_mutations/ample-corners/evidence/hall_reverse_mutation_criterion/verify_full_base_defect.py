"""Independent finite check of the full-base shadow defect and corner formula.

Exhausts every labelled Q4 family containing the fixed base Q2 x {00}.
Arbitrary (including empty and nonample) sections are included.
The symbolic theorem, not this finite check, proves the general identity.
"""
import json
from pathlib import Path


def submasks(t):
    s=t
    while True:
        yield s
        if s==0:break
        s=(s-1)&t


def shadow(a,n):
    return {s for s in range(1<<n) if len({x&s for x in a})==1<<s.bit_count()}


def corners(a,n):
    out={}
    for v in a:
        t=sum(1<<j for j in range(n) if v^(1<<j) in a)
        k={v^s for s in submasks(t)}
        if k<=a:out[v]=k
    return out


def main():
    n=2;N=set(range(1<<n));families=[{v for v in N if b&(1<<v)} for b in range(1<<(1<<n))]
    checked=ample=0
    for P in families:
      for Q in families:
       for C in families:
        X=N|{v|(1<<n) for v in P}|{v|(2<<n) for v in Q}|{v|(3<<n) for v in C}
        A,B,D=shadow(P,n),shadow(Q,n),shadow(C,n)
        pc=shadow(P|C,n)-(A|D);qc=shadow(Q|C,n)-(B|D);g=D-(A|B)
        defect=len(shadow(X,n+2))-len(X)
        assert defect==(len(A)-len(P))+(len(B)-len(Q))+(len(D)-len(C))+len(pc)+len(qc)+len(g)
        predicted=set()
        # Base corners: old incident cube is full N.
        for v in N:
            neighbours=[H for H in (P,Q) if v in H]
            if not neighbours:predicted.add(v)
            elif len(neighbours)==1 and neighbours[0]==N:predicted.add(v)
            elif len(neighbours)==2 and P==N and Q==N and C==N:predicted.add(v)
        for v,K in corners(P,n).items():
            if v not in C or K<=C&Q:predicted.add(v|(1<<n))
        for v,K in corners(Q,n).items():
            if v not in C or K<=C&P:predicted.add(v|(2<<n))
        for v,K in corners(C,n).items():
            if (v not in P or K<=P) and (v not in Q or K<=Q):predicted.add(v|(3<<n))
        assert predicted==set(corners(X,n+2))
        if defect==0:
            ample+=1
            # For ample X the base condition simplifies to P,Q proper and cover.
            cornerless=(P|Q==N and P!=N and Q!=N)
            cornerless &= all(v in C and not K<=C&Q for v,K in corners(P,n).items())
            cornerless &= all(v in C and not K<=C&P for v,K in corners(Q,n).items())
            R=C&P&Q
            cornerless &= all(v in R and not K<=R for v,K in corners(C,n).items())
            assert cornerless==(not predicted)
        checked+=1
    report=dict(status='PASS',labelled_full_base_Q4_families=checked,ample_families=ample,
                scope='Every triple of arbitrary Q2 sections; defect identity and exact local corner formula. Simplified cornerless criterion checked on every ample triple.')
    out=Path(__file__).parent.parent/'evidence/full_base_defect';out.mkdir(exist_ok=True)
    (out/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))

if __name__=='__main__':main()
