"""Certify every elementary pc-minor of the 909-word path construction.

All 111 operations are covered: 109 complete corner orders, and two path
contractions by the proved star-reduction theorem. Eighteen relative core
orders expose the structural reason that the core-exclusive minors peel.
No failure of greedy peeling is treated as a negative certificate.
"""
from pathlib import Path
import hashlib,json
from component_label_partitions import corner
from four_vertex_path_obstruction import outmap
HERE=Path(__file__).resolve().parent

def greedy_to(F,keep,n):
 left=set(F);order=[]
 while left!=keep:
  progress=False
  for v in sorted(left-keep):
   if corner(left,v,n):
    left.remove(v);order.append(v);progress=True
  if not progress:raise AssertionError(('unknown greedy outcome',len(left)))
 return order

def image(F,j,op):
 bit=1<<j
 return {v&~bit for v in F} if op=='contract' else {v for v in F if bool(v&bit)==(op=='one')}

def main():
 src=HERE/'four_vertex_path_obstruction.json';data=json.loads(src.read_text())
 X=set(data['family']);E=1<<36;B={v for v in X if not v&E};H={0,2,4,5}
 upper={v&~E for v in X if v&E};U=set(data['contraction'])
 assert B&upper==H and len(X)==909
 active=0
 for v in X:active|=v^min(X)
 assert active==(1<<37)-1
 rows=[];relative=[]
 for j in range(37):
  for op in ['contract','zero','one']:
   F=image(X,j,op);row=dict(coordinate=j,operation=op,size=len(F))
   if 3<=j<12 and op!='one':
    bp=image(B,j,op);order=greedy_to(bp,H,12)
    # In the remaining lower path each endpoint is a corner of its prism.
    full=order+[2,0,4,5]+[v|E for v in data['upper_section_order']]
    rr=outmap(bp,order+[2,0,4,5],12)
    assert len(set(rr.values()))==len(bp)
    assert all((v^w)&(rr[v]|rr[w]) for v in bp for w in bp if v<w)
    relative.append(dict(coordinate=j,operation=op,order=order,retained=sorted(H)))
   elif j==36:
    full=(data['contraction_order'] if op=='contract' else
          data['lower_section_order'] if op=='zero' else
          [v|E for v in data['upper_section_order']])
   elif j in {1,2} and op=='contract':
    # Contraction positivity follows from the stored positive order on U.
    # Reduction commutes with contraction in ample classes; the reduced
    # three-vertex path is a star. These are theorem applications, not DFS.
    D=image(H,j,op);assert len(D)==3
    assert {v&~E for v in F if not v&E}&{v&~E for v in F if v&E}==D
    assert {v&~E for v in F}==image(U,j,op)
    row['certificate']='thm:damp-star-reduction-positive; contraction is a pc-minor of the positive U'
    rows.append(row);continue
   else:full=greedy_to(F,set(),37)
   outmap(F,full,37)
   row['order']=full;rows.append(row)
 assert len(rows)==111 and len(relative)==18
 report=dict(all_checks_passed=True,source_sha256={src.name:hashlib.sha256(src.read_bytes()).hexdigest()},
  family_order=909,dimension=37,relative_core_orders=relative,elementary_minors=rows,
  complete_orders=109,star_theorem_applications=2,
  scope='All elementary minors are positive; minor closure proves every proper pc-minor positive. Together with the structural negative proof, the explicit 909-word class is forbidden. No exhaustive generator claimed.')
 Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
 print('PASS: 111 elementary minors; 109 complete orders, 2 star-theorem applications; 18 relative core orders with independent nonclashing checks. The 909-word class is minor-minimal.')
if __name__=='__main__':main()
