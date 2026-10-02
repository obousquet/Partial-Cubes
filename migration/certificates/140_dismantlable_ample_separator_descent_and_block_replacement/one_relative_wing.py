"""Replay a 1157-vertex application of the one-relative-wing inverse rule.
Negative core and its minimality reuse the established rooted-factor gluing.
The new theorem supplies all proper-minor positivity uniformly, not by search.
"""
from pathlib import Path
import hashlib,json
from component_label_partitions import corner,shadow
from four_vertex_path_obstruction import outmap
HERE=Path(__file__).resolve().parent

def main():
 src=HERE/'relative_wing_square.json';z=json.loads(src.read_text())['controls'][0]
 C=set(z['core']);A=set(z['A']);B=set(z['D']);co=z['ordinary_core_order']
 assert len(C)==289 and [x for x in C if corner(C,x,12)]==[0]
 outmap(C,co,12);outmap(A,z['A_rooted_order'],12);outmap(B,z['D_rooted_order'],12)
 msrc=HERE/'root_face_recovery.json'
 minor_rows=json.loads(msrc.read_text())['factor_rooted_minor_orders'];assert len(minor_rows)==12
 for row in minor_rows:
  j=row['coordinate']
  for op,key in [('zero','root_restriction'),('contract','contraction')]:
   image={x&~(1<<j) for x in C} if op=='contract' else {x for x in C if not x&(1<<j)}
   lift=lambda v:(v&((1<<j)-1))|((v>>j)<<(j+1))
   order=[lift(v) for v in row[key]['order']+[0]]
   assert order[-1]==0;outmap(image,order,12)
 O=C|{x<<12 for x in C};B0=A|{x<<12 for x in C};B1=C|{x<<12 for x in B};U=B0|B1
 R0=B0-O;R1=B1-O
 assert R0=={128,160} and R1=={290<<12} and B0&B1==O
 assert all((x^y).bit_count()!=1 for x in R0 for y in R1)
 o0=z['A_rooted_order'][:-1]+[x<<12 for x in co]
 o1=[x<<12 for x in z['D_rooted_order'][:-1]]+co
 outmap(B0,o0,24);outmap(B1,o1,24)
 lower_prefix=[160,128];left=set(B0)
 for x in lower_prefix:assert corner(left,x,24);left.remove(x)
 assert left==O
 outmap(U,lower_prefix+o1,24)
 E=1<<24;X=B0|{x|E for x in B1};left=set(X)
 for x in lower_prefix:assert corner(left,x,25);left.remove(x)
 assert left==O|{x|E for x in B1}
 assert len(X)==1157 and len(O)==577
 sc,sa,sb=[shadow(F,12) for F in [C,A,B]]
 so=sc|{s<<12 for s in sc};s0=sa|{s<<12 for s in sc};s1=sc|{s<<12 for s in sb}
 assert len(so)==len(O) and len(s0)==len(B0) and len(s1)==len(B1)
 assert s0&s1==so and len(s0|s1)==len(U)
 expected=(s0|s1)|{s|E for s in so};assert len(expected)==len(X)
 assert not any(corner(O,x,24) for x in O)
 report=dict(family=sorted(X),overlap=sorted(O),dimension=25,order=1157,
             wing_orders=[len(B0),len(B1)],exclusive_wings=[sorted(R0),sorted(R1)],
             positive_section_orders=[o0,o1],positive_contraction_order=lower_prefix+o1,
             one_wing_relative_order=lower_prefix,
             source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [src,msrc]},
             replayed_rooted_core_minors=24,
             negative_and_minimal_overlap='O is the established vertex gluing of two root-negative minimal rooted C289 factors; use the cut-vertex characterization.',
             proper_minor_proof='The one-relative-wing inverse theorem transports the two-deletion prefix under old operations, leaving a positive peripheral expansion. New-coordinate images have explicit positive orders above.',
             scope='Application of a uniform sufficient inverse rule and known rooted inputs. All negative-overlap expansions with one forest wing are covered by that rule; arbitrary cyclic wings and terminal inputs remain unclassified.')
 Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
 print('PASS: 1157 vertices, 25 coordinates; positive sections/contraction; one forest wing clears in two steps; all minors covered by theorem')
if __name__=='__main__':main()
