#!/usr/bin/env python3
"""Fixed checks for section deficits along ample contraction residuals.

Exhaust the ample residuals and corner transitions of all active Q3 splits.
Check every contraction corner of the saved cornerless 291-concept source
and one complete saved contraction order. A product example refutes the
unconditional repair that adjoins the entire contraction corner cube.
No obstruction or endpoint search is performed.
"""
import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path

from verify_damp_layer_midpoint_cores import corner, shadow, subsets

ROOT = Path(__file__).resolve().parent


def subfamilies(F):
    words = sorted(F)
    for mask in range(1 << len(words)):
        yield {v for i, v in enumerate(words) if mask >> i & 1}


def sh_within(F, allowed):
    return {mask for mask in allowed
            if len({v & mask for v in F}) == 1 << mask.bit_count()}


def split(H, e):
    bit = 1 << e
    B = [{v & ~bit for v in H if (v >> e & 1) == b} for b in (0, 1)]
    return B, B[0] & B[1], B[0] | B[1]


def state(H, B, A, K, e, sigma_C, sigma_H):
    sides = [side & K for side in B]
    overlap = A & K
    lift = {v for v in H if v & ~(1 << e) in K}
    sigma = sh_within(K, sigma_C)
    sigmas = [sh_within(F, sigma_C) for F in sides]
    sigma_A = sh_within(overlap, sigma_C)
    sigma_lift = sh_within(lift, sigma_H)
    assert len(sigma) == len(K)
    assert sigma == sigmas[0] | sigmas[1]
    deficits = [len(sigmas[b]) - len(sides[b]) for b in (0, 1)]
    overlap_deficit = len(sigma_A) - len(overlap)
    lift_deficit = len(sigma_lift) - len(lift)
    excess_intersection = len(sigmas[0] & sigmas[1]) - len(sigma_A)
    assert min(*deficits, overlap_deficit, lift_deficit, excess_intersection) >= 0
    assert sigmas[0] & sigmas[1] == sigma_A
    assert excess_intersection == 0
    assert lift_deficit == sum(deficits) == overlap_deficit
    return {'sigma': sigma, 'sigmas': sigmas, 'sigma_A': sigma_A,
            'sigma_lift': sigma_lift, 'deficits': deficits,
            'overlap_deficit': overlap_deficit, 'lift_deficit': lift_deficit,
            'excess_intersection': excess_intersection}


def check_transition(before, after, x, B, A):
    lost = [before['sigmas'][b] - after['sigmas'][b] for b in (0, 1)]
    multiplicity = sum(x in side for side in B)
    lost_contraction = before['sigma'] - after['sigma']
    assert len(lost_contraction) == 1
    assert lost_contraction <= lost[0] | lost[1]
    change = after['lift_deficit'] - before['lift_deficit']
    assert change == multiplicity - sum(map(len, lost))
    lost_overlap = before['sigma_A'] - after['sigma_A']
    if x in A:
        assert change == 1 - len(lost_overlap) <= 1
    else:
        assert change == 0 and before['deficits'] == after['deficits']
    return change, [sorted(S) for S in lost]


def local_holes(F, x, dimension):
    directions = sum(1 << e for e in range(dimension) if x ^ (1 << e) in F)
    return {T for T in subsets(directions) if T.bit_count() >= 2 and x ^ T not in F
            and all(x ^ U in F for U in subsets(T) if U != T)}


def first_step(H, B, A, C, x, e, dimension, sigma_C, sigma_H):
    assert x in A and corner(C, x, dimension) and not corner(A, x, dimension)
    support = sum(1 << f for f in range(dimension) if x ^ (1 << f) in C)
    Q = {x ^ T for T in subsets(support)}
    owners = [b for b in (0, 1) if Q <= B[b]]
    assert len(owners) == 1
    b = owners[0]
    assert corner(B[b], x, dimension) and not corner(B[1 - b], x, dimension)
    before = state(H, B, A, C, e, sigma_C, sigma_H)
    after = state(H, B, A, C - {x}, e, sigma_C, sigma_H)
    assert before['lift_deficit'] == 0 and after['lift_deficit'] == 1
    assert after['deficits'][b] == 0 and after['deficits'][1 - b] == 1
    assert after['overlap_deficit'] == 1 and after['excess_intersection'] == 0
    assert before['sigmas'][b] - after['sigmas'][b] == {support}
    assert before['sigmas'][1 - b] == after['sigmas'][1 - b]
    assert before['sigma_A'] == after['sigma_A']
    assert before['sigma_lift'] - after['sigma_lift'] == {support}
    holes = local_holes(A, x, dimension)
    assert holes and holes == local_holes(B[1 - b], x, dimension)
    records = []
    for T in sorted(holes):
        face = {x ^ U for U in subsets(T)}
        missing_tip = x ^ T
        assert missing_tip in B[b] - A
        for F in (A, B[1 - b]):
            assert (F - {x}) & face == face - {x, missing_tip}
        records.append({'support': T, 'other_missing_vertex': missing_tip})
    return {'coordinate': e, 'deleted_vertex': x, 'owner': b,
            'corner_support': support, 'minimal_nonample_faces': records}


def small_checks():
    ample = [F for F in subfamilies(set(range(8))) if len(shadow(F, 3)) == len(F)]
    counts = Counter()
    for H in ample:
        if not H:
            continue
        sigma_H = shadow(H, 3)
        for e in range(3):
            if len({v >> e & 1 for v in H}) < 2:
                continue
            B, A, C = split(H, e)
            sigma_C = shadow(C, 3)
            counts['active_q3_splits'] += 1
            for K in ample:
                if not K <= C:
                    continue
                before = state(H, B, A, K, e, sigma_C, sigma_H)
                for J in range(8):
                    traces = [{v & J for v in K & side} for side in B]
                    assert traces[0] & traces[1] == {v & J for v in A & K}
                counts['ample_q3_residuals'] += 1
                for x in sorted(K):
                    if not corner(K, x, 3):
                        continue
                    after = state(H, B, A, K - {x}, e, sigma_C, sigma_H)
                    check_transition(before, after, x, B, A)
                    counts['q3_corner_transitions'] += 1
                    if K == C and x in A and not corner(A, x, 3):
                        first_step(H, B, A, C, x, e, 3, sigma_C, sigma_H)
                        counts['q3_first_deficit_states'] += 1
    return counts


def fixed_checks():
    source_name = 'damp_empty_core_counterexample.json'
    source = json.loads((ROOT / source_name).read_text())
    record = next(row for row in source['elementary_minor_peelings']
                  if row['label'] == 'cornerless_predecessor')
    H = set(record['family'])
    sigma_H = shadow(H, 12)
    assert len(H) == len(sigma_H) == 291
    assert not any(corner(H, x, 12) for x in H)
    counts, first_states = Counter(), []
    for e in range(12):
        B, A, C = split(H, e)
        sigma_C = shadow(C, 12)
        for x in sorted(C):
            if corner(C, x, 12):
                first_states.append(first_step(H, B, A, C, x, e, 12, sigma_C, sigma_H))
                counts['fixed_first_deficit_states'] += 1
    e = 0
    B, A, C = split(H, e)
    sigma_C = shadow(C, 12)
    order = next(row['peeling'] for row in record['minor_checks']
                 if row['coordinate'] == e and row['side'] is None)
    assert len(order) == len(C) and set(order) == C
    K = set(C)
    before = state(H, B, A, K, e, sigma_C, sigma_H)
    profile = []
    for x in order:
        assert corner(K, x, 12)
        K.remove(x)
        after = state(H, B, A, K, e, sigma_C, sigma_H)
        change, lost = check_transition(before, after, x, B, A)
        profile.append({'deleted_vertex': x, 'in_original_overlap': x in A,
                        'remaining_order': len(K), 'side_deficits': after['deficits'],
                        'overlap_deficit': after['overlap_deficit'],
                        'lift_deficit': after['lift_deficit'],
                        'excess_intersection': after['excess_intersection'],
                        'deficit_change': change, 'lost_side_supports': lost,
                        'lost_overlap_supports': sorted(before['sigma_A'] - after['sigma_A'])})
        before = after
    assert not K and profile[0]['lift_deficit'] == 1 and profile[-1]['lift_deficit'] == 0
    counts['saved_contraction_transitions'] = len(profile)
    counts['saved_profile_maximum_lift_deficit'] = max(r['lift_deficit'] for r in profile)
    counts['saved_profile_maximum_overlap_deficit'] = max(r['overlap_deficit'] for r in profile)
    return source_name, counts, first_states, profile


def completion_control():
    base_C = {0, 1, 2, 3, 6, 7, 8, 9}
    base_A = {0, 2, 6, 7, 8}
    base_Q = {0, 1, 8, 9}
    assert len(shadow(base_C, 4)) == len(base_C) == 8
    assert len(shadow(base_A, 4)) == len(base_A) == 5
    assert len(shadow(base_A | base_Q, 4)) == 8 > len(base_A | base_Q) == 7
    C = {v | (t << 4) for v in base_C for t in range(4)}
    A = {v | (t << 4) for v in base_A for t in (0, 1, 2)}
    Q = {v | (t << 4) for v in base_Q for t in range(4)}
    x = 8
    support = sum(1 << e for e in range(6) if x ^ (1 << e) in C)
    assert Q == {x ^ T for T in subsets(support)}
    assert len(shadow(C, 6)) == len(C) == 32
    assert len(shadow(A, 6)) == len(A) == 15
    assert x in A <= C and corner(C, x, 6) and not corner(A, x, 6)
    U = A | Q
    assert len(U) == 25 and len(shadow(U, 6)) == 28
    # The original four-dimensional failure is retained in the 00 layer.
    assert {v for v in U if v < 16} == base_A | base_Q
    return {'base_upper': sorted(base_C), 'base_lower': sorted(base_A),
            'base_cube': sorted(base_Q), 'upper': sorted(C), 'lower': sorted(A),
            'vertex': x, 'corner_cube': sorted(Q), 'union': sorted(U),
            'union_order': 25, 'union_shadow_order': 28,
            'scope': 'The unconditional corner-cube completion rule fails even when '
                     'the marked vertex is a noncorner of the lower ample family. '
                     'No forbidden-minor assertion is made for this control.'}


def main():
    counts = small_checks()
    source_name, more_counts, first_states, profile = fixed_checks()
    counts.update(more_counts)
    control = completion_control()
    report = {'schema': 'damp-contraction-section-deficits-v1',
              'all_checks_passed': True,
              'input_sha256': {source_name: hashlib.sha256((ROOT / source_name).read_bytes()).hexdigest()},
              'verifier_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'counts': dict(counts), 'fixed_first_states': first_states,
              'saved_contraction_profile': profile, 'cube_completion_counterexample': control,
              'scope': 'Exact shadow-loss accounting and the first nonample state. '
                       'The scalar deficits do not replace attachment data or supply '
                       'intrinsic admissibility for the forbidden-minor classification.'}
    (ROOT / 'damp_contraction_section_deficits_verification.json').write_text(
        json.dumps(report, indent=2) + '\n')
    print(json.dumps(dict(counts), indent=2))
    print('All section-deficit, local co-pair, and failed cube-completion checks passed.')


if __name__ == '__main__':
    main()
