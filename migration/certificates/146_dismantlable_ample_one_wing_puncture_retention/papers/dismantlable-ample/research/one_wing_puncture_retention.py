"""Replay one-wing mixed-puncture retention on the archived 294-word class.

Both marked endpoints have old certified orders. Only one section retains
the overlap; the other has no initial exclusive corner. The new schedule
uses the successful wing, a whole-family endpoint order, and two rooted row
orders. It does not infer negativity from a greedy stall.
"""
from pathlib import Path
import hashlib,json
from component_label_partitions import corner
from four_vertex_path_obstruction import outmap
from four_vertex_path_minors import image,greedy_to

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]

def replay(F,order,target,n):
    left=set(F)
    assert set(order)==left-target and len(order)==len(left)-len(target)
    for v in order:
        assert corner(left,v,n),('not a corner',v)
        left.remove(v)
    assert left==target

def main():
    source=HERE/'endpoint_coordinate_stripping.json'
    base_source=ROOT/'forbidden_pc_minor_campaign/damp_one_corner_peelable_control.json'
    data=json.loads(source.read_text());base=set(json.loads(base_source.read_text())['family'])
    roots={r['original_root']:r for r in data['results']}
    u,v,z=1639,1655,4
    H={x^u for x in base};n=12
    order0=roots[u]['rooted_order']
    order1=[x^u^v for x in roots[v]['rooted_order']]
    assert order0[-1]==0 and order1[-1]==1<<z
    replay(H,order0[:-1],{0},n);replay(H,order1[:-1],{1<<z},n)
    erase=lambda x:(x&((1<<z)-1))|((x>>(z+1))<<z)
    A0={erase(x) for x in H if not x&(1<<z)}
    A1={erase(x) for x in H if x&(1<<z)}
    D=A0&A1;C=A0|A1
    saved=next(c for c in roots[v]['coordinates'] if c['coordinate']==z)['opposite_components']
    assert len(saved)==1 and saved[0]['status']=='positive'
    wing=saved[0]['relative_order'];replay(A0,wing,D,11)
    assert len(A0-D)==71 and len(A1-D)==93
    assert not any(corner(A1,x,11) for x in A1-D)
    witnesses=[]
    for x in sorted(A1-D):
        ds=[1<<j for j in range(11) if x^(1<<j) in A1]
        masks=[0]
        for d in ds:masks += [m|d for m in masks]
        missing=next(m for m in masks if x^m not in A1)
        witnesses.append(dict(word=x,directions=missing,missing=x^missing))
    # Restriction of an endpoint order and last-occurrence contraction.
    a1order=[erase(x) for x in order1 if x&(1<<z)]
    last={erase(x):i for i,x in enumerate(order0)}
    corder=sorted(last,key=last.get)
    replay(A1,a1order[:-1],{0},11);replay(C,corder[:-1],{0},11)
    row=lambda x,t:x|(t<<11)
    P=A0|{row(x,3) for x in A1}|{row(x,t) for x in C for t in [1,2]}
    target={row(0,t) for t in [0,1,3]}
    stages=[
        [row(x,2) for x in wing],
        [row(erase(x),2 if x&(1<<z) else 0) for x in order0[:-1]],
        [row(x,3) for x in a1order[:-1]],
        [row(x,1) for x in corder[:-1]],
    ]
    relative=sum(stages,[])
    assert len(P)==752 and list(map(len,stages))==[71,293,157,228]
    replay(P,relative,target,13)
    full=relative+[0,3<<11,1<<11]
    r=outmap(P,full,13)
    assert all((x^y)&(r[x]|r[y]) for x in P for y in P if x<y)
    # Scope check for the graph-only corollary on the existing fixed input.
    square_source=HERE/'rooted_square_assembly.json'
    first=set(json.loads(square_source.read_text())['first_factor'])
    graph_rows=[]
    for j in range(12):
        if j==1:continue
        for op in ['zero','contract']:
            M=image(first,j,op)
            sections=[{x for x in M if not x&2},{x&~2 for x in M if x&2}]
            overlap=sections[0]&sections[1];sides=[]
            for section in sections:
                left=section-overlap;components=[]
                while left:
                    seed=next(iter(left));vertices={seed};todo=[seed];left.remove(seed)
                    while todo:
                        x=todo.pop()
                        for k in range(12):
                            y=x^(1<<k)
                            if y in left:left.remove(y);vertices.add(y);todo.append(y)
                    edges=sum(x^(1<<k) in vertices for x in vertices for k in range(12))//2
                    components.append(dict(vertices=len(vertices),edges=edges))
                sides.append(components)
            qualifies=any(all(c['vertices']<=8 or c['edges']==c['vertices']-1 for c in side) for side in sides)
            graph_rows.append(dict(coordinate=j,operation=op,components=sides,qualifies=qualifies))
    assert sum(x['qualifies'] for x in graph_rows)==1
    # A different archived pair has two immobile wings, but retains the edge.
    # It therefore remains a positive control, not an eligible falsifier.
    H2={x^1655 for x in base};J2={0,256}
    edge_order=greedy_to(H2,J2,12);replay(H2,edge_order,J2,12)
    for endpoint in [1655,1911]:
        shore=next(c for c in roots[endpoint]['coordinates'] if c['coordinate']==8)
        assert all(t['status']=='negative_no_first_move' for t in shore['opposite_components'])
    report=dict(all_checks_passed=True,
        source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [source,base_source,square_source]},
        original_endpoints=[u,v],marked_coordinate=z,input_order=294,
        section_orders=[len(A0),len(A1)],overlap_order=len(D),contraction_order=len(C),
        successful_wing_order=wing,failed_wing_witnesses=witnesses,
        block_order=len(P),block_dimension=13,hole=2,target=sorted(target),
        stage_lengths=list(map(len,stages)),relative_order=relative,
        fixed_first_input_graph_audit=graph_rows,
        both_immobile_wings_control=dict(original_endpoints=[1655,1911],coordinate=8,
            normalization_xor=1655,full_edge_order=edge_order,
            scope='Both relative sections fail by archived initial-corner certificates, but full-edge retention makes mixed-puncture retention a consequence of the old block theorem.'),
        scope='Constructive check of the uniform one-wing puncture theorem on a case outside the two-relative-wing hypothesis. Both endpoint orders and one wing order are reused; the other wing has no initial corner. No universal endpoint-only converse or global generator is claimed.')
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS: both endpoints of H294; 71-word wing retains overlap, 93-word wing has no first move.')
    print('PASS: 749 deletions retain the mixed-punctured square in W2(H294), order752; independent nonclashing check.')
    print('SCOPE: the graph-only corollary qualifies only 1 of 22 minors of the fixed A289 input; it does not replace its existing certificates.')

if __name__=='__main__':main()
