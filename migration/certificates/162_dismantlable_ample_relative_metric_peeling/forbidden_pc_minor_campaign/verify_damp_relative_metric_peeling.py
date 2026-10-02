#!/usr/bin/env python3
"""Check metric one-concept extensions and relative peeling to an ample target.

Exhaustive small-cube extension and relative-order controls, an isometric
nonample target showing the hypothesis is needed, and a relative peeling
to the fixed cornerless 291-concept ample seed. No endpoint search.
"""
import hashlib
import itertools
import json
from pathlib import Path

from verify_damp_layer_midpoint_cores import corner, shadow, subsets
from verify_damp_rank_three_corner_neighbors import corners, directions

ROOT = Path(__file__).resolve().parent


def isometric(family, dimension):
    for source in family:
        distance = {source: 0}
        queue = [source]
        for v in queue:
            for e in range(dimension):
                w = v ^ (1 << e)
                if w in family and w not in distance:
                    distance[w] = distance[v] + 1
                    queue.append(w)
        if any(distance.get(v) != (source ^ v).bit_count() for v in family):
            return False
    return True


def metric_constraints(family):
    result = []
    for u, v in itertools.combinations(sorted(family), 2):
        if (u ^ v).bit_count() == 2:
            common = [w for w in family if (w ^ u).bit_count() == (w ^ v).bit_count() == 1]
            assert common
            result.append((u, v, common))
    return result


def relative_checks(family, target, order, dimension, constraints):
    assert set(order) == family - target and len(order) == len(family - target)
    current = set(family)
    actual, metric = True, True
    for v in order:
        actual &= corner(current, v, dimension)
        current.remove(v)
        metric &= isometric(current, dimension)
    position = {v: j for j, v in enumerate(order)}
    static = all(any(w in target for w in common) or
                 any(z not in target and position[z] < position[w]
                     for z in [u, v] for w in common)
                 for u, v, common in constraints)
    return bool(actual), bool(metric), static


def main():
    extension_count = ample_family_count = positive_extensions = 0
    for mask in range(1, 1 << 8):
        S = {v for v in range(8) if mask >> v & 1}
        sigma = shadow(S, 3)
        if len(sigma) != len(S):
            continue
        ample_family_count += 1
        for a in set(range(8)) - S:
            I = directions(S, a, 3)
            assert all(a ^ u in S for u in subsets(I) if u)
            H = S | {a}
            tests = [isometric(H, 3), len(shadow(H, 3)) == len(H),
                     I not in sigma, all((a ^ x) & I for x in S)]
            assert len(set(tests)) == 1
            if tests[0]:
                assert corner(H, a, 3)
                assert shadow(H, 3) == sigma | {I}
                positive_extensions += 1
            extension_count += 1

    H, A = set(range(8)), {0, 1}
    constraints = metric_constraints(H)
    relative_orders = valid_orders = 0
    for order in itertools.permutations(sorted(H - A)):
        actual, metric, static = relative_checks(H, A, order, 3, constraints)
        assert actual == metric == static
        valid_orders += actual
        relative_orders += 1
    assert relative_orders == 720 and 0 < valid_orders < relative_orders
    nonample = H - {0, 7}
    assert len(shadow(nonample, 3)) != len(nonample)
    assert relative_checks(H, nonample, [0, 7], 3, constraints) == (False, True, True)

    source_name = 'damp_fixed_291_repair_creation.json'
    source = json.loads((ROOT / source_name).read_text())
    seed = set(source['base'])
    extension = seed | {91}
    assert len(seed) == 291 and 91 not in seed and len(extension) == 292
    assert len(shadow(seed, 12)) == 291 and len(shadow(extension, 12)) == 292
    assert not corners(seed, 12)
    assert corner(extension, 91, 12)
    assert relative_checks(extension, seed, [91], 12,
                           metric_constraints(extension)) == (True, True, True)
    counts = {'nonempty_ample_q3_families': ample_family_count,
              'q3_one_concept_extensions': extension_count,
              'positive_q3_extensions': positive_extensions,
              'relative_cube_to_edge_orders': relative_orders,
              'valid_relative_cube_to_edge_orders': valid_orders,
              'nonample_target_counterexamples': 1,
              'cornerless_ample_target_orders': 1}
    report = {'schema': 'damp-relative-metric-peeling-v1', 'all_checks_passed': True,
              'input_sha256': {source_name: hashlib.sha256((ROOT / source_name).read_bytes()).hexdigest()},
              'verifier_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'counts': counts,
              'nonample_target_control': {'ambient': sorted(H), 'target': sorted(nonample),
                                         'order': [0, 7], 'isometric_prefixes': True,
                                         'corner_peeling': False},
              'nonpeelable_ample_target': {'ambient': sorted(extension), 'target': sorted(seed),
                                         'relative_order': [91], 'target_corners': []},
              'scope': 'Exact relative metric equivalence requires an ample target, but not a '
                       'peelable target. A relative peeling does not promise a complete peeling '
                       'with that target as its final block. No intrinsic classification is claimed.'}
    (ROOT / 'damp_relative_metric_peeling_verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(counts))
    print('All checks passed: one-concept and relative metric criteria, with both target-scope controls.')


if __name__ == '__main__':
    main()
