"""Verify the small rows for simultaneous forbidden-parent corner deletion.
The large parents are certified symbolically by the existing row theorem.
"""
from pathlib import Path
import json,hashlib
from explore_seed267_endpoints import iscorner,verify
HERE=Path(__file__).resolve().parent

def main():
    source=HERE/'all_sinks_corner_failure.json';d=json.loads(source.read_text())
    K=set(d['family']);A=set(d['deleted_family']);c=d['added_corner']
    w,p,q=1<<12,1<<13,1<<14
    square={0,p,q,p|q};L=K|{x|w for x in K}|square
    z=c|w;T=L-{z};assert 0 in A and len(K)==289
    ko=d['all_endpoint_orders']['0'];verify(K,ko,12)
    ao=next(o[1:] for r,o in d['all_endpoint_orders'].items() if int(r) not in [0,c])
    verify(A,ao,12)
    lo=[p|q,p,q]+[x|w for x in ko]+ko
    to=[p|q,p,q]+[x|w for x in ao]+ko
    verify(L,lo,15);verify(T,to,15)
    assert iscorner(L,z,15)
    assert {x&~w for x in L if x&w and not x&(p|q)}==K
    assert {x&~w for x in T if x&w and not x&(p|q)}==A
    result={'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
            'row':sorted(L),'deleted_row':sorted(T),'row_order':lo,'deleted_row_order':to,
            'deleted_corner':z,'input_order':len(K),'parent_order':301*(1<<16)+2*len(K)+1,
            'parent_dimension':28,'parent_vc':19,
            'scope':'Small581/580 occurrence rows and marked289/288 faces are explicit. Both large parents are forbidden by common-seed row classification; minimum degree>=4, distinct signed columns, and complete crossing graph follow symbolically. No assertion that every proper minor of either parent has all sinks attainable.'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print('PASS:581/580 positive rows; corner deletion and289/288 face recovery; symbolic parent order',result['parent_order'])
if __name__=='__main__':main()
