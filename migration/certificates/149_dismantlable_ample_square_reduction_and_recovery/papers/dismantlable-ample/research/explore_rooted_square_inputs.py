"""Bounded audit of marked edges on the two certified unique-corner factors.

Success is witnessed by replayed orders. A greedy stall is unknown, unless
the initial family has no corner outside its retained target. Stop at this
fixed input inventory; passing it is not a universal rooted-factor theorem.
"""
from pathlib import Path
import json
from component_label_partitions import corner
from four_vertex_path_minors import image

HERE = Path(__file__).resolve().parent

def test(F, target, n=12):
    left = set(F)
    initial = [v for v in sorted(left-target) if corner(left, v, n)]
    order = []
    while left != target:
        choices = [v for v in sorted(left-target) if corner(left, v, n)]
        if not choices:
            return dict(status='negative' if not initial else 'unknown',
                        reason='no initial outside corner' if not initial else 'greedy stall',
                        remaining=sorted(left), prefix=order)
        v = choices[0]
        left.remove(v)
        order.append(v)
    replay = set(F)
    for v in order:
        assert corner(replay, v, n)
        replay.remove(v)
    assert replay == target
    return dict(status='positive', order=order)

def main():
    source = json.loads((HERE/'forced_prefix_rooted_inputs.json').read_text())
    rows = []
    for factor in source['factors']:
        A = set(factor['family'])
        if [v for v in sorted(A) if corner(A,v,12)] != [0]:
            continue
        for z in range(12):
            if 1 << z not in A:
                continue
            edge = {0, 1 << z}
            checks = []
            for j in range(12):
                if j == z:
                    continue
                for op in ['zero', 'contract']:
                    checks.append(dict(coordinate=j, operation=op,
                                       retained=sorted(edge), **test(image(A,j,op),edge)))
            for op in ['zero','one','contract']:
                target = {1 << z} if op == 'one' else {0}
                checks.append(dict(coordinate=z,operation=op,
                                   retained=sorted(target),**test(image(A,z,op),target)))
            status = ('negative' if any(c['status']=='negative' for c in checks) else
                      'unknown' if any(c['status']=='unknown' for c in checks) else 'positive')
            rows.append(dict(factor=factor['name'],coordinate=z,status=status,checks=checks))
            print(factor['name'],z,status,
                  [(c['coordinate'],c['operation'],c['status']) for c in checks if c['status']!='positive'],flush=True)
    Path(__file__).with_suffix('.json').write_text(json.dumps(dict(rows=rows,
        scope='Only the certified unique-root-corner factors in the three-factor inventory. Negative means no possible first deletion, unknown means greedy stall.'),indent=2)+'\n')

if __name__ == '__main__':
    main()
