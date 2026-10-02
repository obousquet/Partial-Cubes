"""Check the finite endpoint logic used in peripheral-separator recovery.

Formal states only; no claim that every Boolean state has a realization.
The uniform geometric proof is cor:damp-peripheral-separator-recovery.
"""
from itertools import product
from pathlib import Path
import json


def main():
    states = list(product([0, 1], repeat=3))  # base layer, base sink, site sink
    report = []
    for m in [2, 3]:
        accepted = []
        rejected_peripheral = 0
        for pieces in product(states, repeat=m):
            # Every canonical piece of a forbidden amalgam fails edge retention.
            if any(p and q for b, p, q in pieces):
                continue
            # Contraction and both sections must be positive vertex gluings.
            if sum(not p for b, p, q in pieces) > 1:
                continue
            if any(sum(not (p if side == b else q) for b, p, q in pieces) > 1
                   for side in [0, 1]):
                continue
            if len({b for b, p, q in pieces}) == 1:
                rejected_peripheral += 1
                continue  # Both positive sections are globally nested.
            accepted.append(pieces)
        expected = [((0, 1, 0), (1, 1, 0)), ((1, 1, 0), (0, 1, 0))] if m == 2 else []
        assert sorted(accepted) == sorted(expected)
        report.append(dict(pieces=m, formal_assignments=len(states)**m,
                           surviving_assignments=accepted,
                           excluded_global_peripheral_assignments=rejected_peripheral))
    out = dict(status='PASS', cases=report,
               scope='Finite audit of necessary endpoint logic; not a geometric census. Only opposite bases with both base roots attainable and both site roots unattainable survive. Rooted minimality and unique recovery require the symbolic proof.')
    Path(__file__).with_suffix('.json').write_text(json.dumps(out, indent=2) + '\n')
    print('PASS: 576 formal assignments; two ordered two-piece states, no three-piece state.')


if __name__ == '__main__':
    main()
