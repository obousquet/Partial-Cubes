#!/usr/bin/env python3
"""Check the fixed Hall three-repair example and its clause recovery.

No search for new families and no pc-minor census. The proof that the
constructed family is a terminal obstruction uses the facet-partition
and Hall-obstruction theorems. This verifier checks their repair inputs,
the explicit clause description, thickness, and the three candidate
layer sets. It records and checks actual peelings of the repair inputs.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def subsets(mask: int):
    sub = mask
    while True:
        yield sub
        if not sub:
            break
        sub = (sub - 1) & mask


def corner(row: set[int], x: int, dimension: int) -> bool:
    directions = sum(1 << e for e in range(dimension) if x ^ (1 << e) in row)
    return all(x ^ sub in row for sub in subsets(directions))


def peeling_with_first(row: set[int], first: int) -> list[int] | None:
    remaining = set(row)
    if not corner(remaining, first, 12):
        return None
    order = [first]
    remaining.remove(first)
    while remaining:
        found = next((x for x in sorted(remaining) if corner(remaining, x, 12)), None)
        if found is None:
            return None
        remaining.remove(found)
        order.append(found)
    return order


def main() -> None:
    path = ROOT / 'computations/hall_corner_exchange/source_H.json'
    source = json.loads(path.read_text())
    mask = int(source['source_family_hex'], 16)
    seed = {x for x in range(1 << 12) if mask >> x & 1}
    assert len(seed) == 299 and source['source_label'] == 'H'
    assert not any(corner(seed, x, 12) for x in seed)
    for triple in itertools.combinations(range(12), 3):
        support = sum(1 << e for e in triple)
        assert any(all(x ^ sub in seed for sub in subsets(support)) for x in seed)
    old_clauses = {}
    for coords in itertools.combinations(range(12), 4):
        support = sum(1 << e for e in coords)
        missing = set(subsets(support)) - {x & support for x in seed}
        assert len(missing) == 1
        old_clauses[support] = missing.pop()
    tips = [183, 315, 1836]
    expected_supports = [{2, 3, 6, 7}, {2, 3, 8, 9}, {4, 5, 10, 11}]
    acquired = []
    repair_proofs = []
    for tip, coords in zip(tips, expected_supports):
        support = sum(1 << e for e in range(12) if tip ^ (1 << e) in seed)
        assert support == sum(1 << e for e in coords)
        assert old_clauses[support] == tip & support
        assert all(tip ^ sub in seed for sub in subsets(support) if sub)
        row = seed | {tip}
        # Try only the other initial corners, then use the fixed greedy rule.
        order = next((order for first in sorted(seed) if corner(row, first, 12)
                      if (order := peeling_with_first(row, first)) is not None), None)
        assert order is not None, 'No fixed greedy continuation; do not infer nonpeelability.'
        remaining = set(row)
        for x in order:
            assert x in remaining and corner(remaining, x, 12)
            remaining.remove(x)
        assert not remaining and len(order) == 300
        repair_proofs.append({'tip': tip, 'support': sorted(coords), 'peeling': order})
        acquired.append(support)
    p, q = 1 << 12, 1 << 13
    layer_supports = [p | q, p, q]
    layer_missing = [p | q, 0, 0]
    rows = [{0, p, q}, {p, p | q}, {q, p | q}]
    family = {x | t for x in seed for t in [0, p, q, p | q]}
    for tip, row in zip(tips, rows):
        family.update(tip | t for t in row)
    assert len(family) == 1203
    clauses = [(s, trace) for s, trace in old_clauses.items() if s not in acquired]
    clauses += [(s | k, (tip & s) | z)
                for tip, s, k, z in zip(tips, acquired, layer_supports, layer_missing)]
    solutions = {x for x in range(1 << 14)
                 if all(x & support != trace for support, trace in clauses)}
    assert solutions == family and len(clauses) == 495
    incidence = [[trace >> e & 1 for support, trace in clauses if support >> e & 1]
                 for e in range(14)]
    balanced = [e for e, signs in enumerate(incidence) if sorted(signs) == [0, 1]]
    assert balanced == [12, 13]
    assert [len(signs) for signs in incidence] == [165] * 12 + [2, 2]
    thickness = []
    for e in range(14):
        bit = 1 << e
        sides = [{x & ~bit for x in family if (x >> e & 1) == i} for i in [0, 1]]
        sizes = [len(sides[i] - sides[1-i]) for i in [0, 1]]
        assert min(sizes) >= 1 and sum(sizes) >= 3
        thickness.append(sizes)
    candidates = []
    for marked in [p, q, p | q]:
        old = ((1 << 14) - 1) ^ marked
        old_family = {x for x in subsets(old)
                      if all(x & s != z for s, z in clauses if not s & marked)}
        mixed = [(s, z) for s, z in clauses if s & marked]
        fibres = [{x for x in old_family if x & (s & old) == z & old}
                  for s, z in mixed]
        sizes = [len(fibre) for fibre in fibres]
        success = all(size == 1 for size in sizes)
        assert success == (marked == p | q)
        if success:
            recovered_tips = set().union(*fibres)
            assert recovered_tips == set(tips) and old_family - recovered_tips == seed
        candidates.append({'marked': [e for e in range(14) if marked >> e & 1],
                           'trace_fibre_sizes': sizes, 'singleton_test': success})
    result = {'schema': 'damp-terminal-clause-recovery-v1',
              'source_file': str(path.relative_to(ROOT)),
              'source_file_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
              'source_mathematical_profile_sha256': source['mathematical_profile_sha256'],
              'seed_order': len(seed), 'order': len(family), 'dimension': 14,
              'clause_count': len(clauses), 'coordinate_clause_degrees': [len(s) for s in incidence],
              'balanced_coordinates': balanced, 'exclusive_section_orders': thickness,
              'layer_candidates': candidates, 'repair_proofs': repair_proofs,
              'scope': 'Fixed-witness verification; pc-minimality and terminality follow from the manuscript theorems.'}
    output = ROOT / 'damp_terminal_clause_recovery_verification.json'
    output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: result[k] for k in ['order', 'dimension', 'balanced_coordinates',
                                           'exclusive_section_orders', 'layer_candidates']}))


if __name__ == '__main__':
    main()
