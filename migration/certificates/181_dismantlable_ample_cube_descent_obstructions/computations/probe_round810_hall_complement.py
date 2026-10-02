#!/usr/bin/env python3
"""Find peelings of Hall's complement and its 36 fixed elementary minors.

These are fixed input certificates for the cube-descent counterexample;
no obstruction or large assembly is enumerated.
"""
import hashlib
import heapq
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from verify_damp_cut_vertex_obstruction import check_order
from verify_damp_layer_midpoint_cores import corner
from verify_damp_terminal_cut_vertex import elementary


def greedy(family):
    remaining = set(family)
    todo = [x for x in remaining if corner(remaining, x, 12)]
    heapq.heapify(todo)
    order = []
    initial_corners = sorted(todo)
    while todo:
        x = heapq.heappop(todo)
        if x not in remaining or not corner(remaining, x, 12):
            continue
        remaining.remove(x)
        order.append(x)
        # Only neighbors can acquire corner status when one word is removed.
        # Old candidates are rechecked on extraction, including nonneighbors
        # whose previously full corner cube has lost x.
        for e in range(12):
            y = x ^ (1 << e)
            if y in remaining and corner(remaining, y, 12):
                heapq.heappush(todo, y)
    assert not remaining, (len(order), len(remaining))
    check_order(family, order, 12)
    return initial_corners, order


def main():
    source_path = ROOT / 'computations/hall_corner_exchange/source_H.json'
    source = json.loads(source_path.read_text())
    bits = int(source['source_family_hex'], 16)
    H = {x for x in range(4096) if bits >> x & 1}
    K = set(range(4096)) - H
    assert len(H) == 299 and len(K) == 3797
    initial_corners, order = greedy(K)
    minor_orders = []
    for coordinate in range(12):
        for operation in ['root_section', 'opposite_section', 'contraction']:
            child = elementary(H, coordinate, operation)
            _, peeling = greedy(child)
            minor_orders.append({'coordinate': coordinate,
                                 'operation': operation, 'peeling': peeling})
    result = {'schema': 'damp-hall-complement-peeling-v2',
              'source_sha256': hashlib.sha256(source_path.read_bytes()).hexdigest(),
              'probe_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'source_order': len(H), 'complement_order': len(K),
              'dimension': 12, 'initial_corners': initial_corners,
              'peeling': order, 'seed_minor_peelings': minor_orders,
              'all_checks_passed': True}
    (ROOT / 'computations/round810_hall_complement_peeling.json').write_text(
        json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ['peeling', 'seed_minor_peelings']}, indent=2))
    print('Saved 36 seed-minor peelings and one complement peeling.')


if __name__ == '__main__':
    main()
