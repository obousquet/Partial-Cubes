import sys,json,collections,hashlib,time
from pathlib import Path
sys.path.insert(0,'forbidden_pc_minor_campaign');sys.path.insert(0,'papers/ample-corners/scripts')
from explore_damp_maximum_vc3_corner_mutations import HALL_CONCEPTS,legal_swaps,build_state
from verify_maximum_mutation_criterion import check
from verify_full_base_defect import shadow
start=time.monotonic();X=frozenset(HALL_CONCEPTS);rows=[]
for r,a in legal_swaps(X,12):
 Y=(X-{r})|{a};state=build_state(Y,12);k,info=check(Y,X,12,3);assert k==0
 rows.append(dict(removed=r,added=a,corners=sorted(state.corners),reverse=info,union_shattered=len(shadow(X|{a},12)),corner_distances=[(u^v).bit_count() for u in state.corners for v in state.corners if u<v]))
p=Path('papers/ample-corners/evidence/hall_reverse_mutation_criterion');p.mkdir(parents=True,exist_ok=True)
report=dict(seconds=time.monotonic()-start,checked=len(rows),corner_histogram=dict(collections.Counter(len(r['corners']) for r in rows)),union_shadow_histogram=dict(collections.Counter(r['union_shattered'] for r in rows)),rows=rows,scope='All legal maximum VC3 one-word exchanges of the fixed Hall class. Reuses existing exchange enumeration; tests final-step mechanism, not smaller obstructions.')
report['sources']={}
for name in ['papers/ample-corners/scripts/verify_hall_reverse_mutations.py','papers/ample-corners/scripts/verify_maximum_mutation_criterion.py','papers/ample-corners/scripts/verify_full_base_defect.py','forbidden_pc_minor_campaign/explore_damp_maximum_vc3_corner_mutations.py','explore_smallest_additional_peripheral_obstruction.py']:
 raw=Path(name).read_bytes();report['sources'][name]=hashlib.sha256(raw).hexdigest();(p/Path(name).name).write_bytes(raw)
(p/'manifest.json').write_text(json.dumps(report,indent=2)+'\n');print({k:v for k,v in report.items() if k!='rows'})
