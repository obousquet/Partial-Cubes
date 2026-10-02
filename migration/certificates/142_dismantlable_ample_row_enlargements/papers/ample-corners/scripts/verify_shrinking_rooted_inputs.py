"""Replay rooted inputs and all72 elementary minors of one527-word double."""
from pathlib import Path
import json,hashlib,time
from verify_full_base_defect import shadow,corners,submasks


def replay(H,seq,n,root=None):
    Y=set(H);assert len(seq)==len(Y) and set(seq)==Y
    if root is not None:assert seq[-1]==root
    for v in seq:
        d=sum(1<<i for i in range(n) if v^(1<<i) in Y)
        assert all(v^s in Y for s in submasks(d))
        Y.remove(v)
    assert not Y


def main():
    start=time.monotonic();root=Path(__file__).resolve().parents[1];out=root/'evidence/shrinking_rooted_inputs'
    m=json.loads((out/'manifest.json').read_text());rows=json.loads((out/'checks.json').read_text())
    repo=root.parent.parent
    for path,h in m['sources'].items():assert hashlib.sha256((repo/path).read_bytes()).hexdigest()==h
    assert hashlib.sha256((out/'checks.json').read_bytes()).hexdigest()==m['checks_sha256']
    seed=set(json.loads((repo/'papers/dismantlable-ample/evidence/order267_obstruction/record.json').read_text())['family'])
    original=json.loads((root/'evidence/order267_two_for_one/checks.json').read_text())['candidates']
    checks=0
    for row in rows:
        source_row=original[row['index']]
        reconstructed=(seed-set(source_row['removed']))|{source_row['added']}
        for v in row['prefix']:
            assert set(corners(reconstructed,12))=={v}
            reconstructed.remove(v)
        assert reconstructed==set(row['family'])
        H=set(row['family']);r=row['root'];assert len(H) in [264,265,266] and len(shadow(H,12))==len(H)
        assert set(corners(H,12))=={r};replay(H,row['ordinary_order'],12)
        assert {(v['coordinate'],v['operation']) for v in row['minor_checks']}=={(e,op) for e in range(12) for op in ['contraction','root_section']}
        assert all(any((x^r)&(1<<e) for x in H) for e in range(12))
        for v in row['minor_checks']:
            bit=1<<v['coordinate'];M={x&~bit for x in H if v['operation']=='contraction' or bool(x&bit)==bool(r&bit)}
            assert M==set(v['family']);replay(M,v['order'],12,r&~bit);checks+=1
    row=min(rows,key=lambda v:len(v['family']));assert set(row['family'])==(seed-{236,492,1004,748})|{220} and row['root']==750;r=row['root'];A={x^r for x in row['family']};X=A|{x<<12 for x in A}
    assert len(X)==527 and not corners(X,24)
    records=[]
    for e in range(24):
        j=e%12;other=12 if e<12 else 0;this=0 if e<12 else 12;bit=1<<e
        for op in ['contraction','zero','one']:
            M={x&~bit for x in X if op=='contraction' or bool(x&bit)==(op=='one')}
            if op=='one':
                seq=[((x^r)<<this)&~bit for x in row['ordinary_order'] if ((x^r)>>j)&1]
            else:
                v=next(v for v in row['minor_checks'] if v['coordinate']==j and v['operation']==('contraction' if op=='contraction' else 'root_section'))
                normalized=[(x^(r&~(1<<j)))<<this for x in v['order']]
                seq=normalized[:-1]+[(x^r)<<other for x in row['ordinary_order']]
            replay(M,seq,24)
            records.append(dict(coordinate=e,operation=op,order=seq))
    (out/'double527.json').write_text(json.dumps(dict(family=sorted(X),dimension=24,seed_index=row['index'],prefix=row['prefix'],root=0,elementary_minor_orders=records),separators=(',',':'))+'\n')
    report=dict(status='PASS',rooted_inputs=len(rows),proper_rooted_minor_orders=checks,double_order=527,double_minors_replayed=len(records),seconds=time.monotonic()-start,verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='All distinct labelled forced-prefix rooted-critical inputs; direct cornerlessness and all72 proper-minor orders for one double. Ampleness and all pairwise doubles follow from the established coordinate-block gluing theorem, not an isomorphism census.')
    (out/'verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
