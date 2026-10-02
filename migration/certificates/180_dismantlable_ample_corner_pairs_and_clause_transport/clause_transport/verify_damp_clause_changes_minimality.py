#!/usr/bin/env python3
"""Check clause transport and loss of minimality within nonpeelable outputs.

Replay one existing adjacent 1188-concept pair, not a new search. Check
nonpeelability by its two rooted pieces and a first-deletion obstruction,
all 78 proper-minor orders, and the positive partner's peeling. Construct
the compact rows and complete signed presentations of the two outputs of
the reduced-base universality theorem. The outputs themselves are never
materialized; their nonpeelability and minimality follow from that theorem.
"""
from pathlib import Path
import json

from verify_damp_reduced_base_universality import (
    check_order, clauses, corner, minor, sha, shadow,
)

ROOT = Path(__file__).resolve().parent


def check_record(data, verifier):
    assert data['all_checks_passed']
    assert data['verifier_sha256'] == sha(ROOT / verifier)
    for name, expected in data['input_sha256'].items():
        assert sha(ROOT / name) == expected, name


def main():
    names = [
        'damp_forbidden_adjacent_sign_pivots_verification.json',
        'damp_root_section_square_pivots_verification.json',
        'damp_rooted_clause_pivots_verification.json',
        'damp_reduced_base_universality_verification.json',
        'computations/hall_corner_exchange/source_H.json',
    ]
    pair, lifted, rooted, universal, source = [
        json.loads((ROOT / name).read_text()) for name in names
    ]
    for data, script in [
        (pair, 'verify_damp_forbidden_adjacent_sign_pivots.py'),
        (lifted, 'verify_damp_root_section_square_pivots.py'),
        (rooted, 'verify_damp_rooted_clause_pivots.py'),
        (universal, 'verify_damp_reduced_base_universality.py'),
    ]:
        check_record(data, script)

    # Check the obstruction to peeling each glued piece to their shared root.
    B = set(lifted['core'])
    base = rooted['square_lift_tower'][0]
    assert B == set(base['negative_family']) and len(B) == 289
    check_order(B, base['negative_ordinary_peeling'], 12)
    assert {x for x in B if corner(B, x, 12)} == {0}
    negative = next(row for row in lifted['four_signs']
                    if row['omitted_trace'] == 2307)
    N = set(negative['family'])
    assert len(N) == 900
    check_order(N, negative['peeling'], 14)
    sections = [{x >> 2 for x in N if x & 3 == z} for z in range(4)]
    assert all(0 in section for section in sections)
    assert all(B <= sections[z] for z in (1, 2, 3))
    assert sections[1] & sections[3] == B
    assert sections[2] & sections[3] == B
    missing_neighbor = lifted['missing_root_neighbor']
    assert missing_neighbor in B and missing_neighbor.bit_count() == 1
    assert missing_neighbor not in sections[0]
    # Before a first deletion in a nonroot B-copy, its two-layer carrier
    # forces the old word to be B's only corner, zero. The three retained
    # nonroot copies, retained root, and old neighbor then require the
    # missing root-section neighbor in that corner cube. This is impossible.
    G = set(pair['negative_glued_family'])
    B_shifted = {x << 14 for x in B}
    assert N & B_shifted == {0} and G == N | B_shifted
    n = pair['dimension']
    assert n == 26 and len(G) == 1188
    sg = shadow(G, n)
    assert all(1 << e in sg for e in range(n))
    seen = set()
    for row in pair['negative_elementary_peelings']:
        key = row['coordinate'], row['operation']
        assert key not in seen
        seen.add(key)
        H = minor(G, *key)
        assert 0 < len(H) < len(G)
        check_order(H, row['glued_peeling'], n)
    assert seen == {(e, op) for e in range(n)
                    for op in ('root_section', 'opposite_section', 'contraction')}
    stage = next(row for row in pair['stages'] if not row['square_completed'])
    positive = next(row for row in stage['signs'] if row['omitted_trace'] == 2306)
    K = (G - {2418}) | {2419}
    assert len(K) == 1188 and 0 in K and 0 in G
    assert G - K == {2418} and K - G == {2419}
    check_order(K, positive['peeling'], n)
    sk = shadow(K, n)
    assert sk == sg == set(stage['common_shadow'])
    ck, cg = clauses(K, sk, n), clauses(G, sg, n)
    I = pair['changed_support']
    assert I == 6531 and I.bit_count() == 6
    assert ck.keys() == cg.keys()
    assert {s for s in ck if ck[s] != cg[s]} == {I}
    assert (ck[I], cg[I]) == (2306, 2307)
    assert (ck[I] ^ cg[I]).bit_count() == 1
    print('Input pair checked: 78 proper-minor orders, positive peeling, '
          'two rooted obstructions, actual cube shadows and signed clauses.', flush=True)

    bits = int(source['source_family_hex'], 16)
    S = {x for x in range(4096) if bits >> x & 1}
    ss = shadow(S, 12)
    assert len(S) == 299
    assert ss == {s for s in range(4096) if s.bit_count() <= 3}
    tips, supports = [183, 315, 1836], [204, 780, 3120]
    cs = clauses(S, ss, 12)
    for tip, support in zip(tips, supports):
        assert tip not in S and corner(S | {tip}, tip, 12)
        assert sum(1 << e for e in range(12) if tip ^ (1 << e) in S) == support
        assert cs[support] == tip & support
    old = {s: t for s, t in cs.items() if s not in supports}
    assert len(old) == 492
    w, p, q, z = [1 << (n + e) for e in range(4)]
    m, full = n + 4, (1 << (n + 4)) - 1
    output_signings, rows, reports = [], [], []
    for name, H, ch, order in [('peelable', K, ck, positive['peeling']),
                               ('forbidden', G, cg, None)]:
        U = H | {x | w for x in H}
        L = U | {p, q, p | q}
        assert len(L) == 2379
        sl = shadow(L, m)
        c0 = clauses(L, sl, m)
        expected = dict(ch)
        expected.update({(1 << e) | t: (1 << e) | t
                         for e in range(n + 1) for t in (p, q)})
        expected[z] = z
        assert c0 == expected
        if order is not None:
            check_order(L, [p | q, p, q] +
                        [x | t for x in order for t in (0, w)], m)
        signed = dict(old)
        for J, trace in c0.items():
            signed[supports[0] | (J << 12)] = (tips[0] & supports[0]) | (trace << 12)
        signed[supports[1] | (full << 12)] = tips[1] & supports[1]
        signed[supports[2] | (full << 12)] = (tips[2] & supports[2]) | (full << 12)
        patterns = [tuple(bool(s >> e & 1) for s in signed) for e in range(m + 12)]
        assert all(any(v) for v in patterns) and len(set(patterns)) == m + 12
        assert {x for x in L if not x & (w | p | q | z)} == H
        assert max(3 + m, 4 + max(s.bit_count() for s in sl), 4 + m - 1) == m + 3
        output_signings.append(signed)
        rows.append(L)
        reports.append(dict(
            input_status=name, input_order=len(H), enlarged_row_order=len(L),
            enlarged_row_peeling_replayed=(order is not None),
            input_face_recovered=True,
            output_order=301 * (1 << m) + 2 * len(H) + 1,
            output_dimension=m + 12, output_vc_dimension=m + 3,
            output_dimension_difference=9, distinct_support_columns=len(patterns),
            output_clause_count=len(signed),
            minimum_degree_lower_bound_by_uniform_theorem=4,
            complete_crossing_graph_by_uniform_theorem=True,
            output_nonpeelable_by_uniform_theorem=True,
            output_forbidden_by_uniform_theorem=(order is not None),
            full_output_materialized=False, output_minor_orders_replayed=0,
        ))
    a, b = output_signings
    assert a.keys() == b.keys()
    changed = {s for s in a if a[s] != b[s]}
    J = supports[0] | (I << 12)
    assert changed == {J} and J.bit_count() == 10
    assert a[J] ^ b[J] == 1 << 12
    assert sum((a[s] ^ b[s]).bit_count() for s in a) == 1
    assert len(rows[0] - rows[1]) == len(rows[1] - rows[0]) == 2
    report = dict(
        schema=1, all_checks_passed=True, verifier_sha256=sha(Path(__file__)),
        helper_sha256={'verify_damp_reduced_base_universality.py':
                       sha(ROOT / 'verify_damp_reduced_base_universality.py')},
        input_sha256={name: sha(ROOT / name) for name in names},
        input_positive_orders_replayed=1, rooted_piece_ordinary_orders_replayed=2,
        input_negative_elementary_orders_replayed=78,
        negative_certificate='Two rooted pieces; unique corner in B and two fixed '
                             'overlaps with missing root neighbor in N.',
        fixed_seed_certificates='Hashes of the Round-823 verification and all its '
                                'inputs agree; its 36 minor and 3 repair orders '
                                'are reused, not replayed by this script.',
        input_changed_support=I, input_changed_traces=[ck[I], cg[I]],
        output_changed_support=J, output_changed_traces=[a[J], b[J]],
        output_changed_literal=12, output_changed_clause_count=1,
        concepts_changed_per_output=2, outputs=reports,
        scope='The symbolic clause-transport theorem covers all fixed-shadow inputs '
              'containing zero. This compact replay checks one existing adjacent pair '
              'and its loss of output minimality while nonpeelability persists. '
              'No output enumeration, output-minor replay, new search, improved '
              'cardinality bound, or exhaustive intrinsic classification is claimed.',
    )
    target = ROOT / 'damp_clause_changes_minimality_verification.json'
    target.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({'outputs': reports, 'changed_support': J,
                      'changed_traces': [a[J], b[J]]}), flush=True)
    print('All checks passed:', target, flush=True)


if __name__ == '__main__':
    main()
