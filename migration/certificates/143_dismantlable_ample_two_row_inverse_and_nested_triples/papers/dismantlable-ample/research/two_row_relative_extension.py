"""Replay a genuine two-row extension via one relative enlarged row.

The 3480-word output extends the certified3348 one-row example by132 words.
General forbiddenness follows from relative-chain transport, not this example.
"""
import json,sys,hashlib,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
sys.path.insert(0,str(ROOT/'papers/ample-corners/scripts'))
from verify_rectangular_obstruction import peel
from verify_full_base_defect import corners,shadow
from verify_shrinking_rooted_inputs import replay
from one_row_enlargement import relative_replay,image,transport
HERE=Path(__file__).resolve().parent

def main():
 start=time.monotonic();source=HERE/'one_row_enlargement.json';d=json.loads(source.read_text())
 pair=ROOT/'papers/ample-corners/evidence/relative_pair_q6/manifest.json';inp=json.loads(pair.read_text())
 P=set(inp['outer']);O=set(inp['inner']);P0={x>>1 for x in P};O0={x>>1 for x in O}
 H=set(d['H']);C=set(d['C']);D0=set(d['D']);T=set(range(32))
 D1={x|(y<<5) for x in T for y in O}|{x|(y<<5) for x in O0 for y in P}
 assert D0&D1==C and C<D0 and C<D1
 rel=peel(T,5,P0);seq=[x|(y<<5) for x in rel for y in peel(O,6)]
 relative_replay(D1,seq,C,11)
 X=set(d['family']);prefix=[x|4096 for x in seq];V=X|set(prefix)
 relative_replay(V,prefix,X,13)
 assert len(V)==len(shadow(V,13))==3480
 assert not corners(H,12) and not corners(X,13)
 expected={x|2048 for x in corners(D0,11) if x not in C}|{x|4096 for x in corners(D1,11) if x not in C}
 assert set(corners(V,13))==expected
 # Every projected relative chain followed by the old certified minor order.
 records=[]
 for r in d['proper_minor_orders']:
  q,b=r['coordinate'],r['operation'];base=image(X,q,b);target=image(V,q,b)
  vals=[x for x in transport(prefix,q,b) if x not in base]
  order=vals+r['order'];replay(target,order,13)
  records.append(dict(coordinate=q,operation=b,relative_prefix=vals,order=order))
 assert len(records)==39
 result=dict(status='PASS',order=len(V),relative_deletions=len(prefix),core_order=len(X),row_orders=[len(D0),len(D1)],intersection_order=len(C),corners=sorted(expected),family=sorted(V),relative_prefix=prefix,proper_minor_orders=records,source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,pair]},seconds=time.monotonic()-start,scope='One explicit two-row instance; ampleness, all corners, relative deletion and39 elementary minor orders checked. General criterion is symbolic.')
 (HERE/'two_row_relative_extension.json').write_text(json.dumps(result,separators=(',',':'))+'\n')
 print(json.dumps({k:v for k,v in result.items() if k not in ['family','relative_prefix','proper_minor_orders']}))
if __name__=='__main__':main()
