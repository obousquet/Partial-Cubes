"""Bounded test of prescribed sinks in elementary pc-minors of seed267.

Success stores a complete order for each endpoint. Greedy failures are unknown,
never negative. Stop at the first endpoint unresolved by three fixed priorities.
The decision is whether proper-minor rooted failures can supply new primitive
inputs, or this seed forces a compatibility failure beyond every such input.
"""
from pathlib import Path
import hashlib,json
from itertools import combinations
HERE=Path(__file__).resolve().parent
SOURCE=HERE.parent/'evidence/order267_obstruction/record.json'

def dirs(H,x,n):
    return sum(1<<i for i in range(n) if x^(1<<i) in H)

def iscorner(H,x,n):
    d=dirs(H,x,n);m=d
    while m:
        if x^m not in H:return False
        m=(m-1)&d
    return True

def order(H,root,n,priority):
    left=set(H);out=[];corners={x for x in left if iscorner(left,x,n)}
    while len(left)>1:
        allowed=corners-{root}
        if not allowed:return None
        if priority==0:x=max(allowed,key=lambda x:((x^root).bit_count(),x))
        elif priority==1:x=min(allowed)
        else:x=max(allowed)
        left.remove(x);out.append(x);corners.remove(x)
        affected={x^(1<<i) for i in range(n)}&left
        # Non-neighbours can only lose cornerhood after this deletion.
        # Recheck all previous corners, and neighbours that may gain it.
        candidates=corners|affected
        corners={y for y in candidates if iscorner(left,y,n)}
    return out+[root]

def verify(H,o,n):
    assert len(o)==len(H) and set(o)==H
    # Independent order check: outgoing sets are the later unit neighbours.
    # Every corresponding cube must lie in the suffix; also check nonclashing.
    times={x:i for i,x in enumerate(o)};maps={}
    for x in o:
        d=sum(1<<i for i in range(n) if times.get(x^(1<<i),-1)>times[x])
        maps[x]=d;m=d
        while m:
            assert times.get(x^m,-1)>=times[x],(x,m)
            m=(m-1)&d
    assert len(set(maps.values()))==len(H)
    for x,y in combinations(o,2):
        assert (x^y)&(maps[x]|maps[y])

def main():
    d=json.loads(SOURCE.read_text());X=set(d['family']);n=d['n'];rows=[];unknown=None
    assert len(X)==267 and not any(iscorner(X,x,n) for x in X)
    for e in range(n):
        bit=1<<e
        for op in ['contraction','zero','one']:
            H={x&~bit for x in X if op=='contraction' or bool(x&bit)==(op=='one')}
            orders={}
            for r in sorted(H):
                o=None
                for priority in range(3):
                    o=order(H,r,n,priority)
                    if o is not None:break
                if o is None:unknown=dict(coordinate=e,operation=op,root=r);break
                verify(H,o,n);orders[str(r)]=o
            rows.append(dict(coordinate=e,operation=op,family=sorted(H),orders=orders,complete=len(orders)==len(H)))
            print(e,op,len(H),len(orders),'complete' if len(orders)==len(H) else 'UNKNOWN',flush=True)
            if unknown:break
        if unknown:break
    out=dict(source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
             verifier_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             all_elementary_minors_complete=unknown is None and len(rows)==3*n,
             orders_verified=sum(len(r['orders']) for r in rows),unknown=unknown,rows=rows,
             scope='Every reported order is checked through suffix cubes, distinct out-sets and nonclashing. Unknown means three greedy priorities failed; no negative inference. Complete elementary coverage implies every sink of every nonempty proper pc-minor is attainable by minor transport.')
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
    print('RESULT',out['all_elementary_minors_complete'],out['orders_verified'],'unknown',unknown,flush=True)
if __name__=='__main__':main()
