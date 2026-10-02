"""Check the mutation formula on all Q3 maximum exchanges and saved Q11 edges."""
import collections,hashlib,itertools,json,time
from pathlib import Path
from verify_full_base_defect import shadow,corners

def check(X,Y,n,d):
 r,=X-Y;a,=Y-X
 R={j for j in range(n) if r^(1<<j) in X};A={j for j in range(n) if a^(1<<j) in Y};assert R==A
 dx={v:sum(v^(1<<j) in X for j in range(n)) for v in X}
 dy={v:sum(v^(1<<j) in Y for j in range(n)) for v in Y}
 assert min(dx.values())>=d and min(dy.values())>=d
 U={v for v in X-{r} if (v^r).bit_count()==1 and (v^a).bit_count()!=1}
 V={v for v in X-{r} if (v^a).bit_count()==1 and (v^r).bit_count()!=1}
 old={v for v in X if dx[v]==d};new={v for v in Y if dy[v]==d}
 births={v for v in U if dx[v]==d+1};losses={v for v in V if dx[v]==d}
 assert len(new)-len(old)==len(births)-len(losses)
 assert (not new)==(len(R)>d and old<=V and not births)
 return len(new),dict(births=sorted(births),losses=sorted(losses),removed_degree=len(R))

def main():
 start=time.monotonic();toy=0
 for d in range(4):
  target={s for s in range(8) if s.bit_count()<=d}
  families=[set(C) for C in itertools.combinations(range(8),len(target)) if shadow(set(C),3)==target]
  for X in families:
   for Y in families:
    if len(X-Y)!=1:continue
    check(X,Y,3,d);toy+=1
    assert set(corners(X,3))=={v for v in X if sum(v^(1<<j) in X for j in range(3))==d}
 root=Path(__file__).resolve().parents[1];source=root/'evidence/maximum_two_corner_component/checks.json';data=json.loads(source.read_text());count=0;hist=collections.Counter()
 for row in data['processed']:
  X=set(data['families'][row['id']])
  for r,a,k,_ in row['exchanges']:
   new,info=check(X,(X-{r})|{a},11,3);assert new==k;count+=1;hist[(len(info['births']),len(info['losses']))]+=1
 out=root/'evidence/maximum_mutation_criterion';out.mkdir(parents=True,exist_ok=True)
 report=dict(status='PASS',seconds=time.monotonic()-start,all_Q3_maximum_directed_exchanges=toy,saved_Q11_directed_exchanges=count,birth_loss_histogram=[dict(births=a,losses=b,count=k) for (a,b),k in sorted(hist.items())],source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),scope='Finite verification of the independently proved formula. No enumeration of all Q11 classes or proof of a mutation barrier.')
 (out/Path(__file__).name).write_bytes(Path(__file__).read_bytes());(out/'manifest.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))

if __name__=='__main__':main()
