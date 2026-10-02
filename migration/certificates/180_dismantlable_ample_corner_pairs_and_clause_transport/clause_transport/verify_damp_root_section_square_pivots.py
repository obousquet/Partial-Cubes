#!/usr/bin/env python3
"""Replay a rooted-minimal adjacent clause change in an enlarged root section.

All sixteen occurrence rows are checked: fourteen are ample, eight give
rooted endpoint obstructions, and six have explicit endpoint orders.
Nonpeelability uses two fixed overlaps and a missing root neighbor,
not an exhaustive search. The 900-concept negative member has all 28
proper rooted elementary images checked by explicit orders.
"""
import hashlib
import json
from pathlib import Path

from verify_damp_cut_vertex_obstruction import check_order
from verify_damp_forbidden_clause_pivots import signing
from verify_damp_layer_midpoint_cores import corner, subsets
from verify_damp_paired_root_transfer import square
from verify_damp_rank_three_corner_neighbors import corners, directions
from verify_damp_single_clause_completions import small_shadow
from verify_damp_terminal_cut_vertex import elementary

ROOT = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def add_corner(family, sigma, word, dimension):
    assert word not in family
    support = directions(family, word, dimension)
    assert support not in sigma
    assert all(word ^ u in family for u in subsets(support) if u)
    family.add(word)
    sigma.add(support)
    assert corner(family, word, dimension)
    return support


def row_order(row):
    remaining, order = set(row), []
    while remaining:
        cs = corners(remaining, 2)
        assert cs
        v = min(cs)
        order.append(v)
        remaining.remove(v)
    return order


def root_lock(family, core, missing_neighbor):
    """Validate the hypotheses of the first-deletion obstruction proof."""
    sections = [{v >> 2 for v in family if v & 3 == z} for z in range(4)]
    assert sections[1] & sections[3] == core
    assert sections[2] & sections[3] == core
    assert corners(core, 12) == [0]
    assert all(0 in section for section in sections)
    assert missing_neighbor in core and missing_neighbor.bit_count() == 1
    assert missing_neighbor not in sections[0]
    return {'core': sorted(core), 'core_corners': [0],
            'opposite_edge_overlaps': [[1, 3], [2, 3]],
            'root_section': sorted(sections[0]),
            'missing_root_neighbor': missing_neighbor,
            'reason': 'Before the first deletion in a nonroot-layer core copy, '
                      'both displayed overlaps are the unchanged core. Such a deletion '
                      'would have to project to its only corner zero. At a zero copy, '
                      'the two layer directions and the missing-neighbor direction '
                      'would force that neighbor in the root section.'}


def lift_minor_orders(source, final, additions, source_orders, dimension):
    output = []
    for e in range(dimension):
        for kind in ['root_section', 'contraction']:
            base = elementary(source, e, kind)
            current, extra = set(base), []
            for v in additions:
                if kind == 'root_section' and v & (1 << e):
                    continue
                w = v & ~(1 << e)
                if w not in current:
                    assert w != 0
                    current.add(w)
                    extra.append(w)
            child = elementary(final, e, kind)
            assert current == child and len(child) < len(final)
            order = list(reversed(extra)) + source_orders[e, kind]
            check_order(child, order, dimension, last=0)
            output.append({'coordinate': e, 'operation': kind, 'peeling': order})
    return output


def main():
    names = ['damp_rooted_clause_pivots_verification.json',
             'damp_partial_tip_endpoints_verification.json']
    rooted, old_rows = [json.loads((ROOT / n).read_text()) for n in names]
    for data, script in [(rooted, 'verify_damp_rooted_clause_pivots.py'),
                         (old_rows, 'verify_damp_partial_tip_endpoints.py')]:
        assert data['all_checks_passed'] and data['verifier_sha256'] == digest(ROOT / script)
        for name, value in data['input_sha256'].items():
            assert digest(ROOT / name) == value
    first, second = rooted['square_lift_tower'][:2]
    B = set(first['negative_family'])
    alpha, delta, gamma, I = 128, 290, 604, 1632
    assert corners(B, 12) == [0]
    S = square(B, B | {alpha}, B | {delta})
    assert S == set(second['negative_family']) == set(old_rows['base_rooted_factor'])
    sigma_s = set(second['equal_shadow'])
    check_order(S, second['negative_ordinary_peeling'], 14)
    source_orders = {(r['coordinate'], r['operation']): r['peeling']
                     for r in second['rooted_minor_peelings']}
    assert len(source_orders) == 28
    face_mask = I | (gamma & ~I)
    D = {v for v in B if not v & ~face_mask}
    assert len(D) == 28 and face_mask == 1660
    assert {gamma ^ u for u in subsets(I) if u} <= D
    assert 2 in B and 2 not in D
    initial_order = next(r['peeling'] for r in first['rooted_minor_peelings']
                         if r['coordinate'] == 0 and r['operation'] == 'root_section')
    d_order = [v for v in initial_order if not v & ~face_mask]
    check_order(D, d_order, 12, last=0)
    D_sigma = small_shadow(D)
    assert len(D_sigma) == len(D)
    d_additions = [v << 2 for v in reversed(d_order[:-1])]
    T, sigma_t, d_gains = set(S), set(sigma_s), []
    for v in d_additions:
        gain = add_corner(T, sigma_t, v, 14)
        assert gain.bit_count() >= 3
        d_gains.append(gain)
    assert T == S | {v << 2 for v in D} and len(T) == 897
    assert sigma_t == sigma_s | {j << 2 | 3 for j in D_sigma}
    old_positive = {tuple(r['occurrence_row']): r['endpoint_peeling']
                    for r in next(c for c in old_rows['cases'] if c['tip'] == gamma)['rows']
                    if r.get('peels_to_root')}
    cases, families, additions_by_row = [], {}, {}
    for row_mask in range(16):
        row = {z for z in range(4) if row_mask >> z & 1}
        H = T | {gamma << 2 | z for z in row}
        row_sigma = small_shadow(row) if row else set()
        sigma = sigma_t | {I << 2 | j for j in row_sigma}
        ample = len(row_sigma) == len(row)
        assert ample == (row not in [{0, 3}, {1, 2}])
        record = {'occurrence_row': sorted(row), 'order': len(H), 'ample': ample}
        if not ample:
            assert len(sigma) == len(H) + 1
            assert all(len({v & j for v in H}) == 1 << j.bit_count() for j in sigma)
            record['shattered_support_count_lower_bound'] = len(sigma)
            cases.append(record)
            continue
        row_peeling = row_order(row)
        additions = d_additions + [gamma << 2 | z for z in reversed(row_peeling)]
        built, built_sigma = set(T), set(sigma_t)
        for z in reversed(row_peeling):
            assert add_corner(built, built_sigma, gamma << 2 | z, 14).bit_count() >= 4
        assert built == H and built_sigma == sigma
        sig = signing(H, sigma, 14)
        positive = {1, 3} <= row or {2, 3} <= row
        record.update(peels_to_root=positive, vc_dimension=max(map(int.bit_count, sigma)))
        if positive:
            old_row = tuple(sorted(row - {0}))
            order = ([gamma << 2] if 0 in row else [])
            order += [v << 2 for v in d_order[:-1]] + old_positive[old_row]
            check_order(H, order, 14, last=0)
            record['endpoint_peeling'] = order
        else:
            record['root_lock_certificate'] = root_lock(H, B, 2)
            order = [gamma << 2 | z for z in row_peeling]
            order += [v << 2 for v in d_order[:-1]] + second['negative_ordinary_peeling']
            check_order(H, order, 14)
            record['ordinary_peeling'] = order
            record['rooted_factor'] = True
        record['ample_addition_order_from_paired_factor'] = additions
        cases.append(record)
        families[tuple(sorted(row))] = H
        additions_by_row[tuple(sorted(row))] = additions
    negative_row = (0, 1, 2)
    negative = families[negative_row]
    minor_orders = lift_minor_orders(S, negative, additions_by_row[negative_row], source_orders, 14)
    common = families[(0, 1, 2, 3)]
    J, outside = I << 2 | 3, (gamma & ~I) << 2
    four = []
    common_sigma = None
    for omitted in range(4):
        row = tuple(z for z in range(4) if z != omitted)
        H = families[row]
        record = next(r for r in cases if tuple(r['occurrence_row']) == row)
        sigma = sigma_t | {I << 2 | j for j in [0, 1, 2]}
        sig = signing(H, sigma, 14)
        trace = (gamma & I) << 2 | omitted
        assert {v & ~J for v in common if v & J == trace} == {outside}
        assert H == common - {gamma << 2 | omitted}
        assert sig[J] == trace and len(H) == 900
        if common_sigma is None:
            common_sigma = sig
        else:
            assert {j for j in sig if sig[j] != common_sigma[j]} == {J}
        four.append({'omitted_row': omitted, 'omitted_trace': trace,
                     'deleted_vertex': gamma << 2 | omitted,
                     'peels_to_root': record['peels_to_root'],
                     'family': sorted(H),
                     'peeling': record.get('endpoint_peeling', record.get('ordinary_peeling'))})
    assert [r['peels_to_root'] for r in four] == [True, True, True, False]
    assert J == 6531 and outside == 112 and len(common) == 901
    assert (four[2]['deleted_vertex'] ^ four[3]['deleted_vertex']) == 1
    counts = {'occurrence_rows': 16, 'ample_rows': 14, 'nonample_rows': 2,
              'positive_endpoint_rows': 6, 'rooted_obstruction_rows': 8,
              'negative_member_rooted_elementary_orders': len(minor_orders),
              'displayed_same_shadow_signs': len(four)}
    assert sum(r.get('peels_to_root') is True for r in cases) == 6
    assert sum(r.get('rooted_factor') is True for r in cases) == 8
    report = {'schema': 'damp-root-section-square-pivots-v1',
              'all_checks_passed': True,
              'input_sha256': {n: digest(ROOT / n) for n in names},
              'verifier_sha256': digest(Path(__file__)),
              'core': sorted(B), 'root': 0, 'root_section': sorted(D),
              'root_section_coordinate_mask': face_mask, 'root_section_peeling': d_order,
              'missing_root_neighbor': 2, 'base_paired_factor': sorted(S),
              'root_section_addition_gains': d_gains,
              'cases': cases, 'common_completion': sorted(common),
              'changed_support': J, 'outside_family': [outside],
              'four_signs': four, 'common_shadow': sorted(sigma),
              'negative_member_rooted_elementary_peelings': minor_orders,
              'counts': counts,
              'scope': 'Refutes adjacent-sign endpoint invariance for rooted factors and the '
                       'single-face implication. The negative factors admit corner descents to '
                       'the known paired factor; no new terminal type or full classification is claimed.'}
    (ROOT / 'damp_root_section_square_pivots_verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(counts))
    print('All checks passed: adjacent 900-concept members have opposite endpoint status; the negative member is rooted-minimal.')


if __name__ == '__main__':
    main()
