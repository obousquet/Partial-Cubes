#!/usr/bin/env python3
"""Verify separation of metric constraints while an ample overlap is retained.

Exhaust all active splits and relative orders in Q3. Check the seven-vertex
face accounting for every new same-layer midpoint under contraction. Then
replay the existing elementary-minor certificates for two fixed forbidden
families and check where their contraction orders first enter the overlap.
This is a bounded theorem check, not an endpoint or obstruction search.
"""
import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path

from verify_damp_layer_midpoint_cores import corner, shadow, subsets
from verify_damp_relative_metric_peeling import isometric, relative_checks

ROOT = Path(__file__).resolve().parent


def intervals(family, dimension):
    neighbors = {v: {v ^ (1 << e) for e in range(dimension)
                     if v ^ (1 << e) in family} for v in family}
    return {(u, v): tuple(sorted(neighbors[u] & neighbors[v]))
            for u, v in itertools.combinations(sorted(family), 2)
            if (u ^ v).bit_count() == 2}


def retained_constraints(data, target):
    return {(u, v, common) for (u, v), common in data.items()
            if not set(common) & target}


def split_check(family, coordinate, dimension):
    bit = 1 << coordinate
    sides = [{v & ~bit for v in family if (v >> coordinate & 1) == b}
             for b in (0, 1)]
    A, C = sides[0] & sides[1], sides[0] | sides[1]
    wings = [side - A for side in sides]
    assert A and all(len(shadow(F, dimension)) == len(F)
                     for F in [A, C, *sides])
    assert not any((u ^ v).bit_count() == 1
                   for u in wings[0] for v in wings[1])
    data = [intervals(F, dimension) for F in [C, *sides]]
    assert all(common for table in data for common in table.values())
    active = [retained_constraints(table, A) for table in data]
    assert active[0] == active[1] | active[2]
    assert not active[1] & active[2]
    for b in (0, 1):
        for u, v, common in active[b + 1]:
            assert {u, v, *common} - A <= wings[b]
    cross = 0
    for (u, v), common in data[0].items():
        if (u in wings[0] and v in wings[1]) or (v in wings[0] and u in wings[1]):
            assert set(common) <= A
            cross += 1
    new_midpoints = []
    for b in (0, 1):
        for (u, v), old in data[b + 1].items():
            new = set(data[0][u, v]) - set(old)
            if not new:
                continue
            assert len(old) == len(new) == 1
            w, z = old[0], next(iter(new))
            assert {u, v, w} <= A and z in wings[1 - b]
            support = (u ^ v) | bit
            face = {u ^ t for t in subsets(support)}
            missing = z | (bit if b else 0)
            assert family & face == face - {missing}
            new_midpoints.append({'side': b, 'endpoints': [u, v],
                                  'old_midpoint': w, 'new_midpoint': z,
                                  'cube_support': support, 'missing_tip': missing})
    return sides, A, C, data, {
        'coordinate': coordinate, 'overlap_order': len(A),
        'wing_orders': list(map(len, wings)),
        'unprotected_constraints': list(map(len, active)),
        'cross_wing_pairs': cross, 'new_midpoints': new_midpoints,
    }


def check_order(family, order, dimension):
    remaining = set(family)
    assert len(order) == len(remaining) and set(order) == remaining
    for v in order:
        assert corner(remaining, v, dimension)
        remaining.remove(v)


def metric_prefix(family, order, dimension):
    remaining = set(family)
    for v in order:
        remaining.remove(v)
        if not isometric(remaining, dimension):
            return False
    return True


def small_checks():
    counts = Counter()
    for mask in range(1, 1 << 8):
        H = {v for v in range(8) if mask >> v & 1}
        if len(shadow(H, 3)) != len(H):
            continue
        counts['ample_q3_families'] += 1
        for e in range(3):
            if len({v >> e & 1 for v in H}) != 2:
                continue
            B, A, C, data, row = split_check(H, e, 3)
            counts['active_q3_splits'] += 1
            counts['q3_new_midpoints'] += len(row['new_midpoints'])
            counts['q3_cross_wing_pairs'] += row['cross_wing_pairs']
            constraints = [[(u, v, w) for (u, v), w in table.items()]
                           for table in data]
            for length in range(len(C - A) + 1):
                for prefix in itertools.permutations(sorted(C - A), length):
                    whole = metric_prefix(C, prefix, 3)
                    sides = [metric_prefix(B[b], [v for v in prefix if v in B[b]], 3)
                             for b in (0, 1)]
                    assert whole == all(sides)
                    counts['q3_partial_orders'] += 1
            for order in itertools.permutations(sorted(C - A)):
                whole = relative_checks(C, A, order, 3, constraints[0])
                wings = [relative_checks(B[b], A, [v for v in order if v in B[b]],
                                         3, constraints[b + 1]) for b in (0, 1)]
                assert len(set(whole)) == 1
                assert all(len(set(result)) == 1 for result in wings)
                assert whole[0] == all(result[0] for result in wings)
                counts['q3_relative_orders'] += 1
                counts['q3_valid_relative_orders'] += whole[0]
    return counts


def fixed_checks():
    source_name = 'damp_empty_core_counterexample.json'
    source = json.loads((ROOT / source_name).read_text())
    counts, rows = Counter(), []
    for record in source['elementary_minor_peelings']:
        H = set(record['family'])
        assert len(shadow(H, 12)) == len(H)
        if record['label'] == 'cornerless_predecessor':
            assert H == set(source['predecessor']) and len(H) == 291
            assert not any(corner(H, v, 12) for v in H)
        else:
            assert record['label'] == 'counterexample'
            assert H == set(source['counterexample']) and len(H) == 292
            assert [v for v in sorted(H) if corner(H, v, 12)] == [3959]
            assert H - {3959} == set(source['predecessor'])
        seen, contraction_orders = set(), {}
        for minor in record['minor_checks']:
            e, b = minor['coordinate'], minor['side']
            image = {v & ~(1 << e) for v in H if b is None or (v >> e & 1) == b}
            assert len(image) == minor['order']
            check_order(image, minor['peeling'], 12)
            seen.add((e, b))
            if b is None:
                contraction_orders[e] = minor['peeling']
            counts['replayed_elementary_orders'] += 1
        assert seen == {(e, b) for e in range(12) for b in (None, 0, 1)}
        family_rows = []
        for e in range(12):
            B, A, C, _, row = split_check(H, e, 12)
            # Produce and directly verify one fixed overlap peeling. There is
            # no branching or use of this order as construction admissibility.
            remaining, overlap_order = set(A), []
            while remaining:
                choices = [v for v in sorted(remaining) if corner(remaining, v, 12)]
                assert choices, (record['label'], e, 'fixed overlap peeling stalled')
                v = choices[0]
                overlap_order.append(v)
                remaining.remove(v)
            check_order(A, overlap_order, 12)
            row['overlap_peeling'] = overlap_order
            order = contraction_orders[e]
            first = next(i for i, v in enumerate(order) if v in A)
            exterior_prefix = set(order[:first])
            remaining_wings = [len(side - A - exterior_prefix) for side in B]
            assert all(remaining_wings)
            if record['label'] == 'cornerless_predecessor':
                assert first == 0
                created_corners = {v for v in C if corner(C, v, 12)}
                assert created_corners <= A
                assert all(not corner(A, v, 12) for v in created_corners)
                row['contraction_corners_noncorner_in_overlap'] = sorted(created_corners)
                counts['created_corners_noncorner_in_overlap'] += len(created_corners)
                assert len(shadow(A - {order[first]}, 12)) != len(A) - 1
                counts['first_deletions_make_overlap_nonample'] += 1
            row.update(first_overlap_index=first, first_overlap_vertex=order[first],
                       first_overlap_vertex_is_corner_in_overlap=corner(A, order[first], 12),
                       remaining_wing_orders=remaining_wings)
            family_rows.append(row)
            counts['fixed_coordinate_splits'] += 1
            counts['fixed_overlap_peelings'] += 1
            counts['fixed_new_midpoints'] += len(row['new_midpoints'])
        rows.append({'source_label': record['label'], 'order': len(H),
                     'coordinate_checks': family_rows})
    return source_name, counts, rows


def main():
    counts = small_checks()
    source_name, fixed_counts, rows = fixed_checks()
    counts.update(fixed_counts)
    report = {
        'schema': 'damp-overlap-metric-separation-v1', 'all_checks_passed': True,
        'input_sha256': {source_name: hashlib.sha256((ROOT / source_name).read_bytes()).hexdigest()},
        'verifier_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'counts': dict(counts), 'fixed_examples': rows,
        'scope': 'Relative constraints split while the entire overlap is retained. '
                 'In the gluing-critical case a contraction peeling must enter that '
                 'overlap before either wing is exhausted. For a cornerless source '
                 'its first deletion makes the overlap nonample. This rejects retaining '
                 'the canonical overlap as a general contraction-repair strategy; '
                 'it does not classify the required overlap deletions.',
    }
    (ROOT / 'damp_overlap_metric_separation_verification.json').write_text(
        json.dumps(report, indent=2) + '\n')
    print(json.dumps(dict(counts), indent=2))
    print('All overlap-separation, midpoint-face, and forced-entry checks passed.')


if __name__ == '__main__':
    main()
