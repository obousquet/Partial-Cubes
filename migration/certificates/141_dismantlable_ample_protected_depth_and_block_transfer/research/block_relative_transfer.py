"""Fixed 683-word witness outside the one-relative-wing inverse branch.
Forbiddenness uses the proved block-replacement theorem and the certified
267-word seed. This replay checks the seed's 36 proper minors and new section
orders; it certifies relative failure by absence of any first deletion.
"""
from pathlib import Path
import json,hashlib
from component_label_partitions import corner,shadow,greedy_order,replay
HERE=Path(__file__).resolve().parent

def main():
 source=HERE.parent/'evidence/order267_obstruction/record.json'
 H=set(json.loads(source.read_text())['family']);assert len(H)==267
 assert len(shadow(H,12))==267 and not any(corner(H,x,12) for x in H)
 seed_orders=[]
 for j in range(12):
  for op in ('zero','one','contract'):
   F={x&~(1<<j) for x in H if op=='contract' or (x>>j&1)==(op=='one')}
   o=greedy_order(F,12);assert o is not None;seed_orders.append(dict(coordinate=j,operation=op,order=o))
 B0={x>>1 for x in H if not x&1};B1={x>>1 for x in H if x&1};A=B0&B1;C=B0|B1
 # New coordinates are 11 and 12, with 12 marking the expansion.
 W=B0|{x|(3<<11) for x in B1}|{x|(t<<11) for x in C for t in (1,2)}
 O=B0|{x|(1<<11) for x in B1}
 L0={x for x in W if not x&(1<<12)};L1={x&~(1<<12) for x in W if x&(1<<12)}
 assert L0&L1==O and len(W)==683 and len(O)==267
 assert len(shadow(W,13))==len(W) and not any(corner(W,x,13) for x in W)
 orders=[]
 for b,L,B in [(0,L0,B0),(1,L1,B1)]:
  bo=greedy_order(B,11);co=greedy_order(C,11);assert bo is not None and co is not None
  bo=bo+list(B-set(bo));co=co+list(C-set(co))
  order=([x for x in bo]+[x|(1<<11) for x in co]) if b==0 else ([x|(1<<11) for x in bo]+co)
  replay(L,order[:-1],12);orders.append(order)
  assert not any(corner(L,x,12) for x in L-O)
 report=dict(source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),family=sorted(W),dimension=13,order=683,negative_reduction=sorted(O),exclusive_wing_orders=[len(L0-O),len(L1-O)],section_orders=orders,seed_minor_orders=seed_orders,relative_failure='Neither section has a corner outside the retained overlap. Hence neither relative deletion can start.',forbiddenness='Coordinate-block replacement of the certified 267-word forbidden seed. Output proper-minor positivity follows uniformly from that theorem.',scope='New relative-boundary application of an existing construction, not a new exhaustive generator or a claim of a previously unknown family.')
 Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
 print('PASS: 683-word forbidden block expansion; overlap 267; wings 90 and 59; both relative peelings have no first move; all 36 seed minor orders replayed.')
if __name__=='__main__':main()
