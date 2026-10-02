"""Independent checks of all discovered classes; no completeness claim."""
import hashlib,itertools,json,time
from pathlib import Path

def main():
 out=Path(__file__).resolve().parents[1]/'evidence/maximum_two_corner_component'
 r=json.loads((out/'checks.json').read_text());m=json.loads((out/'manifest.json').read_text());start=time.monotonic()
 assert hashlib.sha256((out/'checks.json').read_bytes()).hexdigest()==m['checks_sha256']
 supports={k:[sum(1<<i for i in c) for c in itertools.combinations(range(11),k)] for k in [3,4]}
 families=[set(x) for x in r['families']];corners=[]
 assert len({tuple(sorted(x)) for x in families})==len(families)
 for X in families:
  assert len(X)==232
  assert all(len({x&s for x in X})==8 for s in supports[3])
  assert all(len({x&s for x in X})<16 for s in supports[4])
  # Exact maximum cardinality and trace tests imply ampleness.
  C=set()
  for v in X:
   cube={v}
   for i in range(11):
    if v^(1<<i) in X:cube|={w^(1<<i) for w in tuple(cube)}
   if cube<=X:C.add(v)
  assert len(C)==2;corners.append(C)
 reach={0};changed=True;edges=0
 while changed:
  changed=False
  for row in r['processed']:
   assert set(row['corners'])==corners[row['id']]
   for a,b,count,t in row['exchanges']:
    if t is None:continue
    assert count==2 and a not in corners[row['id']]
    assert families[t]==(families[row['id']]-{a})|{b}
    if row['id'] in reach and t not in reach:reach.add(t);changed=True
 assert reach==set(range(len(families)))
 combined=json.loads((out/'combined_corner_moves.json').read_text());X=set(families[0])
 for row in combined['steps']:
  a,b=row['exchange'];assert a in X and b not in X;X.remove(a);X.add(b)
  assert len(X)==232 and all(len({x&s for x in X})==8 for s in supports[3])
  assert all(len({x&s for x in X})<16 for s in supports[4])
  C=set()
  for v in X:
   cube={v}
   for i in range(11):
    if v^(1<<i) in X:cube|={w^(1<<i) for w in tuple(cube)}
   if cube<=X:C.add(v)
  assert C==set(row['corners'])
 assert X==set(combined['family'])
 report=dict(status='PASS',seconds=time.monotonic()-start,classes=len(families),corner_sets=[dict(corners=list(C),count=sum(D==set(C) for D in corners)) for C in sorted({tuple(sorted(D)) for D in corners})],scope='Every discovered class is maximum VC3 and has exactly two corners; every class is reachable by saved legal exchanges. The component is not closed; no barrier or global exclusion follows.',verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
 (out/Path(__file__).name).write_bytes(Path(__file__).read_bytes());(out/'audit.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))

if __name__=='__main__':main()
