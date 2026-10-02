"""Inspect double-reduction leaves in the twelve known Hall controls.

This tests whether a leaf forces an exposed corner. All inputs are the
known299-word class, so it supplies no new obstruction or smaller bound.
"""
import collections,hashlib,itertools,json,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from explore_smallest_additional_peripheral_obstruction import HALL_CONCEPTS
from verify_full_base_defect import submasks,shadow

def main():
 out=Path(__file__).resolve().parents[1]/'evidence/final_mutation_tree';out.mkdir(parents=True,exist_ok=True)
 source=out.parent/'final_mutation_two_types/manifest.json';controls=json.loads(source.read_text())['Hall_controls'];H=set(HALL_CONCEPTS);rows=[];covers=collections.Counter();leafcounts=collections.Counter()
 for control in controls:
  r=control['removed'];order=control['coordinate_order'];X={sum((((v^r)>>j)&1)<<k for k,j in enumerate(order)) for v in H}
  groups={t:[int(j) for j,k in control['types'].items() if k==t] for t in [3,5]};small,=[g for g in groups.values() if len(g)==2];mask=sum(1<<j for j in small);a_words=list(submasks(mask))
  tree={v for v in X if not v&mask and all(v|a in X for a in a_words)}
  active=[j for j in range(12) if j not in small];edges=[(v,v^(1<<j),j) for v in tree for j in active if v<(v^(1<<j)) and v^(1<<j) in tree]
  assert len(tree)==11 and len(edges)==10 and sorted(j for _,_,j in edges)==active
  assert shadow(tree,12)=={0}|{1<<j for j in active}
  leaves=[]
  for v in sorted(tree):
   dirs=[j for j in active if v^(1<<j) in tree]
   if len(dirs)!=1:continue
   k=dirs[0];profiles={j:sorted(a for a in a_words if (v|a)^(1<<j) in X) for j in active if j!=k};profiles={j:E for j,E in profiles.items() if E}
   assert all(0<len(E)<4 for E in profiles.values())
   assert set().union(*(set(E) for E in profiles.values()))==set(a_words)
   minimum=next(size for size in range(1,len(profiles)+1) if any(set().union(*(set(profiles[j]) for j in js))==set(a_words) for js in itertools.combinations(profiles,size)))
   covers[minimum]+=1;leaves.append(dict(vertex=v,tree_direction=k,escape_profiles=profiles,minimum_cover=minimum))
  leafcounts[len(leaves)]+=1;rows.append(dict(removed=r,added=control['added'],reduction_directions=small,tree=sorted(tree),edges=edges,leaves=leaves))
 report=dict(status='PASS',controls=len(rows),leaf_count_histogram=dict(leafcounts),minimum_cover_histogram=dict(covers),rows=rows,sources={},scope='Known Hall controls only. Double-reduction trees and all proper escape profiles are exact. Does not exclude a dimension-sensitive lifting theorem.')
 for p in [Path(__file__),Path(__file__).with_name('verify_full_base_defect.py'),ROOT/'explore_smallest_additional_peripheral_obstruction.py',source]:
  raw=p.read_bytes();report['sources'][str(p.relative_to(ROOT))]=hashlib.sha256(raw).hexdigest();(out/('control_manifest.json' if p==source else p.name)).write_bytes(raw)
 (out/'manifest.json').write_text(json.dumps(report,indent=2)+'\n');print({k:v for k,v in report.items() if k not in ['rows','sources']})

if __name__=='__main__':main()
