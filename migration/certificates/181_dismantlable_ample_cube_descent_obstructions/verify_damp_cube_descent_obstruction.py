#!/usr/bin/env python3
"""Replay compact inputs for an obstruction with nonpeelable complement.

The 1,236,691-concept output is certified by the uniform common-seed
theorem, not by materializing it or enumerating its deletion states.
This file independently replays all fixed input orders and verifies the
row conditions, the complement section, and the output cardinality.
"""
import hashlib
import itertools
import json
from pathlib import Path

from verify_damp_cornered_maximal_seed import exact_shadow
from verify_damp_cut_vertex_obstruction import check_order
from verify_damp_layer_midpoint_cores import corner, subsets
from verify_damp_terminal_cut_vertex import elementary

ROOT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    inputs = ['computations/hall_corner_exchange/source_H.json',
              'damp_terminal_clause_recovery_verification.json',
              'computations/round810_hall_complement_peeling.json']
    source, repairs, certificate = [json.loads((ROOT / p).read_text())
                                    for p in inputs]
    hashes = {p: sha(ROOT / p) for p in inputs}
    assert certificate['source_sha256'] == hashes[inputs[0]]
    assert certificate['probe_sha256'] == sha(
        ROOT / 'computations/probe_round810_hall_complement.py')
    assert repairs['source_mathematical_profile_sha256'] == source[
        'mathematical_profile_sha256']
    bits = int(source['source_family_hex'], 16)
    cube, full = set(range(4096)), 4095
    H = {x for x in cube if bits >> x & 1}
    assert len(H) == 299
    expected_shadow = {s for s in cube if s.bit_count() <= 3}
    assert exact_shadow(H, expected_shadow, 12) == 495
    assert not any(corner(H, x, 12) for x in H)
    assert all(any(x ^ (1 << e) in H for x in H) for e in range(12))

    minor_keys = set()
    minor_sizes = []
    for row in certificate['seed_minor_peelings']:
        key = row['coordinate'], row['operation']
        assert key not in minor_keys
        minor_keys.add(key)
        child = elementary(H, *key)
        assert 0 < len(child) < len(H)
        check_order(child, row['peeling'], 12)
        minor_sizes.append(len(child))
    assert minor_keys == set(itertools.product(
        range(12), ['root_section', 'opposite_section', 'contraction']))

    tips = [183, 315, 1836]
    repair_orders = {r['tip']: r['peeling'] for r in repairs['repair_proofs']}
    supports = []
    for a in tips:
        assert a not in H
        P = sum(1 << e for e in range(12) if a ^ (1 << e) in H)
        assert P.bit_count() == 4
        assert all(a ^ s in H for s in subsets(P) if s)
        assert a & P not in {x & P for x in H}
        exact_shadow(H | {a}, expected_shadow | {P}, 12)
        check_order(H | {a}, repair_orders[a], 12)
        supports.append(P)
    assert len(set(supports)) == 3
    assert all((a ^ b).bit_count() > 1 for a, b in itertools.combinations(tips, 2))

    K = cube - H
    check_order(K, certificate['peeling'], 12)
    # Ampleness of K follows from complement duality for the verified ample H.
    rows = {tips[0]: K, tips[1]: cube - {0}, tips[2]: cube - {full}}
    check_order(rows[tips[1]], sorted(rows[tips[1]], key=lambda x: (x.bit_count(), x)), 12)
    check_order(rows[tips[2]], sorted(rows[tips[2]], key=lambda x: (-x.bit_count(), x)), 12)
    assert all(0 < len(row) < len(cube) for row in rows.values())
    facets = []
    for e in range(12):
        for sign in [0, 1]:
            facet = {t for t in cube if (t >> e & 1) == sign}
            owner = tips[1] if sign else tips[2]
            assert facet <= rows[owner]
            facets.append({'coordinate': e, 'sign': sign, 'owner': owner})

    def in_output(old, layer):
        return old in H or layer in rows.get(old, set())

    # A complete check of one 12-dimensional section, not of the full output.
    complement_section = {t for t in cube if not in_output(tips[0], t)}
    assert complement_section == H
    assert not in_output(tips[1], 0)  # the section is a proper pc-minor
    output_order = len(H) * len(cube) + sum(map(len, rows.values()))
    assert output_order == (len(H) + 3) * 4096 - len(H) - 2 == 1236691
    assert 2 ** 24 - output_order == 15540525

    result = {
        'schema': 'damp-cube-descent-obstruction-v1',
        'input_sha256': hashes,
        'verifier_sha256': sha(Path(__file__)),
        'helper_sha256': {p: sha(ROOT / p) for p in [
            'verify_damp_cornered_maximal_seed.py',
            'verify_damp_cut_vertex_obstruction.py',
            'verify_damp_layer_midpoint_cores.py',
            'verify_damp_terminal_cut_vertex.py']},
        'source': {'order': len(H), 'dimension': 12, 'vc_dimension': 3,
                   'corner_count': 0, 'elementary_minor_orders_replayed': 36,
                   'elementary_minor_sizes': minor_sizes},
        'repairs': [{'tip': a, 'support': [e for e in range(12) if P >> e & 1],
                     'peeling_length': len(repair_orders[a])}
                    for a, P in zip(tips, supports)],
        'rows': {'orders': [len(rows[a]) for a in tips],
                 'peelings_replayed': 3, 'signed_facet_owners': facets},
        'output': {'order': output_order, 'dimension': 24,
                   'complement_order': 2 ** 24 - output_order,
                   'complement_section_old_word': tips[0],
                   'complement_section_equals_H': True,
                   'common_seed_theorem_inputs_verified': True,
                   'full_output_materialized': False,
                   'full_output_minor_orders_replayed': 0},
        'conclusion_scope': 'Uniform theorem plus fixed input replay; no output census.',
        'all_checks_passed': True,
    }
    (ROOT / 'damp_cube_descent_obstruction_verification.json').write_text(
        json.dumps(result, indent=2) + '\n')
    print('PASS: Hall ample/cornerless, 36 seed-minor orders, 3 repairs, 3 row orders,')
    print('24 signed facets, exact Hall complement section; output order 1,236,691.')
    print(result['conclusion_scope'])


if __name__ == '__main__':
    main()
