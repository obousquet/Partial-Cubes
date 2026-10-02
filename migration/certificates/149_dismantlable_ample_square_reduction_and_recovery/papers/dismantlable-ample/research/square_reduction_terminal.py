"""Check the square reduction identities on positive and negative controls.

The universal equivalence and minimal-counterexample transfer are proved in
the main source. These controls do not establish existence of a terminal
forbidden square attachment or settle the endpoint-only converse.
"""
from pathlib import Path
import hashlib,json
from component_label_partitions import corner

HERE=Path(__file__).resolve().parent

def reduce(F,j):
    return {x for x in F if not x&(1<<j) and x^(1<<j) in F}

def replay(F,order,root,n):
    left=set(F)
    assert set(order)==left and len(order)==len(left) and order[-1]==root
    for x in order[:-1]:
        assert corner(left,x,n)
        left.remove(x)
    assert left=={root}

def crossing(F,n):
    return {(i,j) for i in range(n) for j in range(i+1,n)
            if len({((x>>i)&1)|(((x>>j)&1)<<1) for x in F})==4}

def recover(X,n,edges):
    universal=[i for i in range(n) if all((min(i,j),max(i,j)) in edges for j in range(n) if j!=i)]
    if len(universal)!=2:return dict(universal=universal,candidates=[])
    e,f=universal;em=(1<<e)|(1<<f)
    remaining=set(range(n))-set(universal);components=[]
    while remaining:
        start=min(remaining);part={start};todo=[start];remaining.remove(start)
        while todo:
            i=todo.pop()
            neighbours={j for j in remaining if (min(i,j),max(i,j)) in edges}
            remaining-=neighbours;part|=neighbours;todo.extend(neighbours)
        components.append(sorted(part))
    candidates=[]
    for part in components:
        lm=sum(1<<i for i in part);am=((1<<n)-1)^em^lm
        if not am:continue
        fibres={}
        code=lambda x:((x>>e)&1)|(((x>>f)&1)<<1)
        for x in X:fibres.setdefault(x&lm,set()).add(code(x))
        roots=[b for b,t in fibres.items() if t==set(range(4))]
        if len(roots)!=1:continue
        rb=roots[0];other=[t for b,t in fibres.items() if b!=rb]
        if not other or len(other[0])!=3 or any(t!=other[0] for t in other):continue
        ra_values={x&am for x in X if x&lm!=rb}
        if len(ra_values)!=1:continue
        ra=next(iter(ra_values));hole=next(iter(set(range(4))-other[0]))
        rows={t:{x&am for x in X if x&lm==rb and code(x)==t} for t in range(4)}
        if any(ra not in row for row in rows.values()):continue
        diagonals=[]
        for diagonal in [(0,3),(1,2)]:
            sides=sorted(set(range(4))-set(diagonal));union=rows[diagonal[0]]|rows[diagonal[1]]
            if hole not in diagonal and rows[sides[0]]==rows[sides[1]]==union:diagonals.append(list(diagonal))
        if diagonals:candidates.append(dict(second_coordinates=part,first_root=ra,second_root=rb,
            hole=hole,exceptional_diagonals=diagonals,second_factor_order=len(fibres)))
    return dict(universal=universal,private_components=components,candidates=candidates)

def main():
    names=['rooted_square_assembly.json','forced_prefix_rooted_inputs.json',
           'verify_edge_boundary_transfer.json','cube_boundary_retention.json',
           'square_attachment_sources.json']
    square,factors,edge,cube,punctures=[json.loads((HERE/n).read_text()) for n in names]
    bf=next(f for f in factors['factors'] if f['name']=='B289')
    B=set(bf['family']);replay(B,bf['ordinary_order'],bf['ordinary_order'][-1],12)
    assert [x for x in sorted(B) if corner(B,x,12)]==[0]
    A=set(square['first_factor']);assert A==B
    root2=factors['attainable_root_orders']['2'];replay(A,root2,2,12)
    for endpoint,order in zip(edge['edge'],edge['endpoint_orders']):
        replay(set(edge['family']),order,endpoint,13)
    rows=[]
    for label,A,n,z,endpoints in [
        ('negative_1601',A,12,1,[False,True]),
        ('positive_2029',set(edge['family']),13,12,[True,True]),
    ]:
        erase=lambda x:(x&((1<<z)-1))|((x>>(z+1))<<z)
        A0={erase(x) for x in A if not x&(1<<z)}
        A1={erase(x) for x in A if x&(1<<z)};C=A0|A1
        offset=n-1;e=offset+12;f=e+1;hole=2
        P=A0|{x|(3<<e) for x in A1}|{x|(t<<e) for x in C for t in [1,2]}
        K={t<<e for t in [0,1,3]}
        X=P|{k|(b<<offset) for k in K for b in B}
        checks=[]
        for old in range(n):
            if old==z:continue
            j=old if old<z else old-1
            assert reduce(X,j)==reduce(P,j)
            checks.append(dict(coordinate=j,type='first_factor',positive=True,
                               reason='Reduction of the positive block P'))
        for old in range(12):
            j=offset+old
            expected={k|(b<<offset) for k in K for b in reduce(B,old)}
            assert reduce(X,j)==expected
            checks.append(dict(coordinate=j,type='second_factor',positive=True,
                               reason='Product of K with a reduction of positive B'))
        for deleted,surviving,b in [(e,f,0),(f,e,1)]:
            expected={erase(x)|(((x>>z)&1)<<surviving) for x in A}
            expected|={x<<offset|(b<<surviving) for x in B}
            assert reduce(X,deleted)==expected
            checks.append(dict(coordinate=deleted,type='square',attachment_endpoint=b,
                               positive=endpoints[b],reason='Vertex-gluing criterion and certified endpoint status'))
        assert len(checks)==n+13
        assert len(X)==(1601 if label=='negative_1601' else 2029)
        cr=crossing(X,n+13)
        expected_cr=crossing(C,offset)|{(i+offset,j+offset) for i,j in crossing(B,12)}
        expected_cr|={(min(i,j),max(i,j)) for i in [e,f] for j in range(n+13) if i!=j}
        assert cr==expected_cr
        recovery=recover(X,n+13,cr)
        assert recovery['universal']==[e,f]
        assert len(recovery['candidates'])==1
        assert recovery['candidates'][0]['exceptional_diagonals']==[[0,3]]
        assert any(c['second_coordinates']==list(range(offset,offset+12)) and
                   c['hole']==2 and [0,3] in c['exceptional_diagonals'] for c in recovery['candidates'])
        ordinary=None
        if label=='negative_1601':
            assert not any(corner(X,x,n+13) for x in X)
        else:
            embed=lambda x:(x&4095)|((x>>12)<<e)
            assert {embed(x) for x in cube['family']}==P
            saved=next(r for r in punctures['punctured_boundary_certificates'] if r['hole']==hole)
            ordinary=[embed(x) for x in saved['relative_order']]
            ordinary += [(b<<offset)|(t<<e) for b in bf['ordinary_order'] for t in [0,3,1]]
            replay(X,ordinary,ordinary[-1],n+13)
        rows.append(dict(control=label,order=len(X),dimension=n+13,
                         reductions=checks,positive_reductions=sum(c['positive'] for c in checks),
                         ordinary_order=ordinary,unmarked_recovery=recovery))
    assert rows[0]['positive_reductions']==24 and rows[1]['positive_reductions']==26
    seed_path=HERE.parent/'evidence/order267_obstruction/record.json'
    growth_path=HERE/'protected_core_delay.json'
    seed=json.loads(seed_path.read_text());growth=json.loads(growth_path.read_text())
    exclusions=[]
    for label,F in [('seed267',set(seed['family'])),('terminal627',set(growth['base_family'])|{32}),
                    ('terminal630',set(growth['full_family']))]:
        n=max(F).bit_length();result=recover(F,n,crossing(F,n))
        exclusions.append(dict(control=label,order=len(F),dimension=n,recovery=result))
    assert all(not r['recovery']['candidates'] for r in exclusions)
    out=dict(all_checks_passed=True,
        source_sha256={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names},
        controls=rows,
        existing_seed_recovery=exclusions,
        recovery_source_sha256={str(p.relative_to(HERE.parent)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [seed_path,growth_path]},
        scope='All private and square reduction identities checked on two archived controls. The 1601 obstruction has exactly one negative coordinate reduction; the positive 2029 attachment has all reductions positive. Neither is a terminal forbidden square attachment.')
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print('PASS: all 51 reduction identities; 1601 control has 24 positive/1 negative, positive2029 has 26 positive.')
    print('No terminal forbidden square witness or universal endpoint-converse proof is claimed.')
    print('Existing seed universal coordinates:',[(r['control'],r['recovery']['universal']) for r in exclusions])
    print('Unmarked presentations recovered:',[(r['control'],len(r['unmarked_recovery']['candidates'])) for r in rows])

if __name__=='__main__':main()
