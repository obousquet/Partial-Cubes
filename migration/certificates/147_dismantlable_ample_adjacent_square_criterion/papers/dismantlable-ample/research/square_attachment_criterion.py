"""Check the product-attachment interpretation of rooted square assemblies.

Reuses the attachment theorems: no new general attachment theorem is claimed.
Checks a noncorner-root second factor, the exact output corners, and a
positive 2029-word control that retains a punctured square but not the full
square. This last strictness is local, not a counterexample to necessity of
every strong input condition inside a forbidden adjacent-hole assembly.
"""
from pathlib import Path
import hashlib
import json
from component_label_partitions import corner
from four_vertex_path_minors import image

HERE = Path(__file__).resolve().parent

def replay(F, order, target, n):
    left = set(F)
    assert set(order) == left-target and len(order) == len(left)-len(target)
    for x in order:
        assert corner(left,x,n), ('not a corner',x)
        left.remove(x)
    assert left == target

def main():
    names = ['rooted_square_assembly.json','forced_prefix_rooted_inputs.json',
             'cube_boundary_retention.json','square_attachment_sources.json',
             'relative_wing_square.json']
    square, factors, cube, punctures, wings = [json.loads((HERE/n).read_text()) for n in names]
    A = set(square['first_factor'])
    factor = next(f for f in factors['factors'] if f['name']=='D289')
    B = set(factor['family'])
    replay(B,factor['ordinary_order'],set(),12)
    corners = [b for b in sorted(B) if corner(B,b,12)]
    assert corners == [256]
    assert [b for b in sorted(B-{256}) if corner(B-{256},b,12)] == [0]
    # The forced first two deletions prove failure to retain root 0.
    assert {(r['coordinate'],r['operation']) for r in factor['rooted_minors']} == {
        (j,op) for j in range(12) for op in ['zero','contract']}
    for r in factor['rooted_minors']:
        assert r['order'][-1]==0
        replay(image(B,r['coordinate'],r['operation']),r['order'][:-1],{0},12)

    erase = lambda x:(x&1)|((x>>2)<<1)
    A0 = {erase(x) for x in A if not x&2}
    A1 = {erase(x) for x in A if x&2}
    C = A0|A1
    P = A0|{x|(3<<23) for x in A1}|{x|(t<<23) for x in C for t in [1,2]}
    T = {t<<23 for t in range(4)}
    hole = 2
    K = T-{hole<<23}
    Q = {b<<11|(t<<23) for b in B for t in range(4) if t!=hole}|T
    X = P|Q
    attachment = P|{k|(b<<11) for k in K for b in B}
    assert X == attachment and len(X)==1601
    assert [v for v in P if corner(P,v,25)] == [0]
    expected_corners = {(b<<11)|(t<<23) for b in corners for t in [0,3]}
    assert {v for v in X if corner(X,v,25)} == expected_corners
    # Minimality uses the uniform theorem and the archived first-input orders;
    # the new seam is that B is root-minimal without being a unique-root-corner factor.
    for r in square['edge_relative_minor_orders']:
        replay(image(A,r['coordinate'],r['operation']),r['order'],{0,2},12)
    for r in square['projected_endpoint_orders']:
        replay(image(A,1,r['operation']),r['order'],{r['target']},12)

    Y = set(cube['family'])
    row = next(r for r in punctures['punctured_boundary_certificates'] if r['hole']==1)
    target = set(row['target'])
    replay(Y,row['relative_order'],target,14)
    control = wings['controls'][0]
    U = set(control['A'])|set(control['D'])
    # Y is the coordinate block of the established 581-word edge piece.
    block = set(control['A'])|{v|(3<<12) for v in control['D']}
    block |= {v|(t<<12) for v in U for t in [1,2]}
    assert Y == block and len(Y)==1165
    positive = Y|{k|(b<<14) for k in target for b in B}
    order = row['relative_order'] + [k|(b<<14) for b in factor['ordinary_order']
                                    for k in [0,3<<12,2<<12]]
    assert len(positive)==2029
    replay(positive,order,set(),26)
    report = dict(all_checks_passed=True,
        source_sha256={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names},
        noncorner_root_factor='D289',forced_prefix=[256,0],rooted_minor_orders_replayed=24,
        forbidden_assembly=dict(order=len(X),dimension=25,corners=sorted(expected_corners),
            proof='Uniform square theorem; first-factor relative orders and all second-factor rooted minor orders replayed. No output-minor search.'),
        weaker_relative_target=dict(block_order=1165,underlying_edge_piece_order=581,
            puncture=1,positive_attachment_order=2029,positive_attachment_dimension=26,
            ordinary_order=order,
            full_square_failure='Existing cube-boundary-retention theorem, via full-fibre reduction to its root-negative C289; not inferred from a greedy stall.'),
        scope='Exact product-attachment identification and generalized second input. Punctured-square retention is strictly weaker locally; necessity of the stronger first-input conditions on forbidden adjacent-hole assemblies remains unresolved.')
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS: noncorner-root B has 24 rooted minor orders; 1601-word assembly has exactly two corners.')
    print('PASS: punctured-square order gives a full peeling of the 2029-word attachment; full-square failure is the archived structural theorem.')

if __name__=='__main__':
    main()
