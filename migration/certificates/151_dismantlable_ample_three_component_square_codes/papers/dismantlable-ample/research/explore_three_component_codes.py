from pathlib import Path
import sys,json,itertools
r=Path(__file__).resolve().parents[3];sys.path.insert(0,str(r/'papers/dismantlable-ample/research'))
from cube_corner_compatibility import cube_code,corner,section
p=r/'papers/dismantlable-ample/research/cube_corner_compatibility.json';d=json.loads(p.read_text());A=set(d['A']);C=set(d['C']);n=23
rows=[]
for q in itertools.product(range(4),repeat=3):
 if q[0]==q[1] or q[1]==q[2]:continue
 Z=cube_code(A,C,n,2,list(q));assert not any(corner(Z,x,n+2) for x in Z)
 bad=[]
 for e in [n,n+1]:
  for b in (0,1):
   F=section(Z,e,b)
   if not any(corner(F,x,n+1) for x in F):bad.append((e,b,len(F)))
 rows.append(dict(codes=q,negative_faces=bad,block=(q[0]==q[2] and q[0]^q[1]==3)))
print('negative',len(rows),'blocks',sum(x['block'] for x in rows),'negative face',sum(bool(x['negative_faces']) for x in rows))
print('unresolved',[x for x in rows if not x['block'] and not x['negative_faces']])
(r/'papers/dismantlable-ample/research/three_component_code_probe.json').write_text(json.dumps(rows,indent=2)+'\n')
