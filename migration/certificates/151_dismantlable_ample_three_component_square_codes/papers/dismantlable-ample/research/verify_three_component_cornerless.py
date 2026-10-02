"""Replay all 75 stored minor orders for a cornerless reduced obstruction.
Also certify all 64 square-code assignments of the fixed three-component pair.
This verifies a fixed construction family, not arbitrary primitive generation.
"""
from pathlib import Path
import sys,json,hashlib,itertools
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parents[2]/'forbidden_pc_minor_campaign'))
from verify_damp_alternating_reductions import cube_shadow,columns
from cube_corner_compatibility import cube_code,expand_binary_trace,section
from component_label_partitions import corner,replay

def main():
 src=HERE/'cube_corner_compatibility.json';cert=HERE/'three_component_code_minors.json'
 d=json.loads(src.read_text());c=json.loads(cert.read_text());A=set(d['A']);C=set(d['C']);H=set(c['family']);g=range(25)
 assert H==set(cube_code(A,C,23,2,[0,3,1])) and len(H)==1601
 assert len(cube_shadow(H,g))==1601 and not any(corner(H,x,25) for x in H)
 cl,groups,_=columns(H,g);assert groups==[[e] for e in g]
 degrees=[sum(x^(1<<e) in H for e in g) for x in H];assert min(degrees)>=3
 assert len(c['minor_rows'])==75
 seen=set()
 for row in c['minor_rows']:
  e,b=row['coordinate'],row['side'];assert (e,b) not in seen;seen.add((e,b))
  M={x&~(1<<e) for x in H if b is None or (x>>e&1)==b}
  o=row['order'];assert not row['residual'] and len(o)==len(M) and set(o)==M
  replay(M,o[:-1],25)
 assert seen=={(e,b) for e in range(25) for b in (0,1,None)}
 O={x for x in H if not x&(1<<23) and x|(1<<23) in H}
 original=set(d['source_leaf']);move=lambda x:(x&1)|((x>>2)<<1)|(((x>>1)&1)<<24)
 target={move(x) for x in original}
 assert O==target or O=={x^(1<<24) for x in target}
 corepath=HERE.parents[2]/'forbidden_pc_minor_campaign/damp_three_core_verification.json'
 assert original==set(json.loads(corepath.read_text())['core_family'])
 rows=[];cases={tuple(z['labels']):z for z in d['binary_cases_up_to_complement']}
 transforms=[]
 for flip in range(4):
  for swap in (False,True):
   f=lambda q:(((q&1)<<1)|((q>>1)&1))^flip if swap else q^flip
   transforms.append((tuple(f(q) for q in (0,3,1)),flip,swap))
 for q in itertools.product(range(4),repeat=3):
  Z=set(cube_code(A,C,23,2,list(q)))
  if q[0]==q[1] or q[1]==q[2]:
   labels=(0,0,1) if q[0]==q[1] else (0,1,1)
   trace=expand_binary_trace(A,C,labels,cases[labels]['peeling'],23,2,list(q));assert trace is not None
   rows.append(dict(codes=q,status='positive',binary_trace=labels))
  elif q[0]^q[1]==3:
   if q[2]==q[0]:
    rows.append(dict(codes=q,status='forbidden',reason='block replacement of certified source'))
   else:
    orbit=[(flip,swap) for values,flip,swap in transforms if values==q];assert len(orbit)==1
    flip,swap=orbit[0]
    def iso(x):
     old=x&((1<<23)-1);t=x>>23;t=(((t&1)<<1)|(t>>1)) if swap else t
     return old|((t^flip)<<23)
    assert {iso(x) for x in H}==Z
    rows.append(dict(codes=q,status='forbidden',reason='isometry of verified 1601 source',flip=flip,swap=swap))
  else:
   found=[]
   for e in (23,24):
    for b in (0,1):
     F=section(Z,e,b)
     if not any(corner(F,x,24) for x in F):found.append((e,b,len(F)))
   assert found;rows.append(dict(codes=q,status='not forbidden',cornerless_proper_face=found[0]))
 counts={s:sum(row['status']==s for row in rows) for s in ['positive','not forbidden','forbidden']};assert counts==dict(positive=28,**{'not forbidden':24},forbidden=12)
 out=dict(order=1601,dimension=25,minimum_degree=min(degrees),negative_overlap_order=len(O),negative_overlap=sorted(O),column_classes=groups,replayed_minor_orders=75,selector='q0 xor q1 == 3 and q2 != q1',code_counts=counts,code_cases=rows,source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [src,cert,corepath]},scope='Fixed 3-component input: all 64 square codes classified. New representative has no corners and no repeated signed columns, but a negative reduction. General generation and input selection remain open.')
 Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
 print('PASS: 1601 words; cornerless; 25 distinct columns; negative overlap577; all75 minor orders replayed; square-code counts',counts)
if __name__=='__main__':main()
