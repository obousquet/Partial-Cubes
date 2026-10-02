"""Bounded falsification: flip one block bit in one mixed clause of W683.
This is not an exhaustive seed search. Stop at this domain; no greedy failure
is used as a negative certificate. Also classify all two-component square codes.
"""
from pathlib import Path
import sys,json,hashlib
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[2]/'forbidden_pc_minor_campaign'))
from verify_damp_alternating_reductions import columns,cube_shadow,corner

def main():
 src=HERE/'block_relative_transfer.json';W=set(json.loads(src.read_text())['family']);cl,_,_=columns(W,range(13))
 rows=[]
 for j,(support,trace) in enumerate(cl):
  if not support&(1<<12):continue
  changed=cl[:j]+[(support,trace^(1<<12))]+cl[j+1:]
  F={x for x in range(8192) if all(x&s!=t for s,t in changed)}
  try:
   sh=cube_shadow(F,range(13));ample=True
  except AssertionError:ample=False
  row=dict(clause=j,support=support,old_trace=trace,order=len(F),ample=ample)
  if ample:
   row.update(corners=sum(corner(F,x,range(13)) for x in F),column_classes=columns(F,range(13),sh)[1])
  rows.append(row)
 codes=[dict(q0=a,q1=b,status='positive' if a==b else 'forbidden_block' if a^b==3 else 'negative_proper_face') for a in range(4) for b in range(4)]
 report=dict(source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),single_mixed_clause_bit_flips=rows,two_component_square_codes=codes,scope='Only one sign flip on the fixed W683 presentation. Code classifications follow the uniform two-component theorem; they are not extrapolation from this search.')
 Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
 print('Single sign changes:',len(rows),'ample:',sum(r['ample'] for r in rows))
 print('Two-component square codes: 4 positive, 8 negative with a negative proper face, 4 forbidden block replacements.')
if __name__=='__main__':main()
