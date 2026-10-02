#!/usr/bin/env python3
"""Replay ordinary forbidden/positive pairs from the rooted adjacent signs.

Glue the four 900-concept choices to the fixed 289-concept rooted factor,
then complete a common root square. Check all signed presentations, six
positive orders, and all 156 elementary orders of the two negative members.
Nonpeelability follows from the checked rooted-factor decomposition and
the proved gluing and square-completion theorems, without endpoint search.
"""
import hashlib
import json
from pathlib import Path

from verify_damp_cut_vertex_obstruction import check_order
from verify_damp_forbidden_clause_pivots import project_order, signing
from verify_damp_layer_midpoint_cores import corner
from verify_damp_terminal_cut_vertex import elementary

ROOT = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    names = ['damp_root_section_square_pivots_verification.json',
             'damp_rooted_clause_pivots_verification.json']
    lifted, old = [json.loads((ROOT / n).read_text()) for n in names]
    for data, script in [(lifted, 'verify_damp_root_section_square_pivots.py'),
                         (old, 'verify_damp_rooted_clause_pivots.py')]:
        assert data['all_checks_passed'] and data['verifier_sha256'] == digest(ROOT / script)
        for name, value in data['input_sha256'].items():
            assert digest(ROOT / name) == value
    base = old['square_lift_tower'][0]
    B = set(lifted['core'])
    assert B == set(base['negative_family'])
    b_order = base['negative_ordinary_peeling']
    b_minors = {(r['coordinate'], r['operation']): r['peeling']
                for r in base['rooted_minor_peelings']}
    negative = next(r for r in lifted['four_signs'] if not r['peels_to_root'])
    N, n_order = set(negative['family']), negative['peeling']
    n_minors = {(r['coordinate'], r['operation']): r['peeling']
                for r in lifted['negative_member_rooted_elementary_peelings']}
    shift, dimension = 14, 26
    B_shifted = {v << shift for v in B}
    shadow = set(lifted['common_shadow']) | {s << shift for s in base['equal_shadow']}
    square_support = (1 << 0) | (1 << 15)
    tip = square_support
    cases = []
    old_signatures = {}
    for completed in [False, True]:
        sigma = shadow | ({square_support} if completed else set())
        stage = []
        for row in lifted['four_signs']:
            H = set(row['family']) | B_shifted
            if completed:
                H.add(tip)
            sig = signing(H, sigma, dimension)
            if old_signatures.get(completed) is None:
                old_signatures[completed] = sig
            else:
                prior = old_signatures[completed]
                assert {s for s in sig if sig[s] != prior[s]} == {lifted['changed_support']}
            assert sig[lifted['changed_support']] == row['omitted_trace']
            record = {'omitted_row': row['omitted_row'], 'omitted_trace': row['omitted_trace'],
                      'order': len(H), 'peelable': row['peels_to_root']}
            assert len(H) == 1188 + completed
            if row['peels_to_root']:
                order = ([tip] if completed else [])
                order += row['peeling'][:-1] + [v << shift for v in b_order]
                check_order(H, order, dimension)
                record['peeling'] = order
            else:
                assert set(row['family']) == N
                record['forbidden_status_certificate'] = {
                    'rooted_factor_orders': [900, 289],
                    'root': 0, 'shared_vertex_only': True,
                    'square_completion_tip': tip if completed else None,
                    'basis': 'rooted-factor gluing, followed by support-two completion when present'}
            stage.append(record)
        cases.append({'square_completed': completed, 'common_shadow': sorted(sigma), 'signs': stage})
    X = N | B_shifted
    assert N & B_shifted == {0}
    assert 1 in N and (1 << 15) in B_shifted
    assert tip not in X and square_support not in shadow
    assert corner(X | {tip}, tip, dimension)
    elementary_orders = []
    for e in range(dimension):
        for kind in ['root_section', 'opposite_section', 'contraction']:
            child = elementary(X, e, kind)
            if e < shift:
                if kind == 'opposite_section':
                    order = project_order(n_order, e, kind)
                else:
                    order = n_minors[e, kind][:-1] + [v << shift for v in b_order]
            else:
                old_e = e - shift
                if kind == 'opposite_section':
                    order = [v << shift for v in project_order(b_order, old_e, kind)]
                else:
                    order = [v << shift for v in b_minors[old_e, kind][:-1]] + n_order
            check_order(child, order, dimension)
            complete_child = elementary(X | {tip}, e, kind)
            extra = complete_child - child
            assert child <= complete_child and len(extra) <= 1
            completed_order = sorted(extra) + order
            check_order(complete_child, completed_order, dimension)
            elementary_orders.append({'coordinate': e, 'operation': kind,
                                      'glued_peeling': order, 'completed_peeling': completed_order})
    assert len(elementary_orders) == 78
    report = {'schema': 'damp-forbidden-adjacent-sign-pivots-v1',
              'all_checks_passed': True,
              'input_sha256': {n: digest(ROOT / n) for n in names},
              'verifier_sha256': digest(Path(__file__)),
              'dimension': dimension, 'changed_support': lifted['changed_support'],
              'adjacent_opposite_status_signs': [2307, 2306],
              'common_square_tip': tip, 'stages': cases,
              'negative_glued_family': sorted(X),
              'negative_elementary_peelings': elementary_orders,
              'counts': {'signed_presentations': 8, 'positive_orders': 6,
                         'negative_elementary_orders': 156},
              'scope': 'One forbidden and three peelable choices on a square of adjacent signs. '
                       'The completed choices are two-connected by the rooted construction and '
                       'square-completion proof. Both negative classes have certified proper minors; '
                       'the full intrinsic classification remains open.'}
    (ROOT / 'damp_forbidden_adjacent_sign_pivots_verification.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report['counts']))
    print('All checks passed for orders 1188 and 1189; adjacent signs change forbidden-minor status.')


if __name__ == '__main__':
    main()
