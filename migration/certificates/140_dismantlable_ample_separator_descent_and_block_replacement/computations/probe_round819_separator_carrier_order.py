#!/usr/bin/env python3
"""Supply one positive order for the fixed terminal edge amalgam's carrier.

This bounded diagnostic uses a deterministic greedy order. Failure would
be unresolved, not a proof of nonpeelability. The separate verifier replays
the successful order and does not search.
"""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parent.parent


def corner(family, vertex, dimension):
    support = sum(1 << e for e in range(dimension)
                  if vertex ^ (1 << e) in family)
    sub = support
    while sub:
        if vertex ^ sub not in family:
            return False
        sub = (sub - 1) & support
    return True


def main():
    source = ROOT / 'damp_terminal_edge_amalgam_verification.json'
    data = json.loads(source.read_text())
    family, e = set(data['family']), 3
    carrier = {x for x in family if not x & (1 << e)
               and x ^ (1 << e) in family}
    left, order = set(carrier), []
    while left:
        v = next((x for x in sorted(left) if corner(left, x, 25)), None)
        if v is None:
            raise RuntimeError('Greedy attempt unresolved; no negative claim')
        order.append(v)
        left.remove(v)
    report = {
        'source': source.name,
        'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'generator_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'coordinate': e,
        'family': sorted(carrier),
        'order': len(carrier),
        'peeling': order,
    }
    out = ROOT / 'computations/round819_separator_carrier_order.json'
    out.write_text(json.dumps(report, indent=2) + '\n')
    print('Supplied positive order:', len(order), 'concepts; coordinate', e)


if __name__ == '__main__':
    main()
