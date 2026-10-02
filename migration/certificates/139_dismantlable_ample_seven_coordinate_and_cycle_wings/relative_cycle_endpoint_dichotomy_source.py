"""Exact small control: the relative leaf dichotomy fails in degree two.

This does not refute relative peeling for cycles. Both endpoints are ample,
the exclusive graph is C6, and an explicit relative order is supplied.
The designated vertex can occur neither first nor last in an addition order.
"""
from pathlib import Path
import json
from component_label_partitions import shadow, corner


def main():
    A = set(range(7)) | {8, 10, 11}
    H = {0, 8, 10, 11}
    D, d, n = A - H, 1, 4
    assert len(shadow(A, n)) == len(A)
    assert len(shadow(H, n)) == len(H)
    cycle = [1, 3, 2, 6, 4, 5]
    edges = {tuple(sorted((x, x ^ (1 << e)))) for x in D for e in range(n)
             if x ^ (1 << e) in D}
    assert edges == {tuple(sorted((cycle[i], cycle[(i+1) % 6]))) for i in range(6)}
    # In H+d, distance-two vertices 1,11 have neither intermediate 3 nor9.
    first = H | {d}
    assert {1, 11} <= first and not {3, 9} & first
    assert len(shadow(first, n)) != len(first)
    # In A-d, distance-two vertices 3,5 have neither intermediate 1 nor7.
    last = A - {d}
    assert {3, 5} <= last and not {1, 7} & last
    assert len(shadow(last, n)) != len(last)
    assert not corner(A, d, n)
    deletion = [5, 6, 4, 1, 3, 2]
    left = set(A)
    for x in deletion:
        assert x in left - H and corner(left, x, n)
        left.remove(x)
        assert len(shadow(left, n)) == len(left)
    assert left == H
    supports = set(shadow(A, n)) - set(shadow(H, n))
    isolated = [s for s in supports if all(
        s == t or s & t not in (s, t) for t in supports)]
    markers = {}
    complement = set(range(1 << n)) - A
    for s in isolated:
        other = ((1 << n) - 1) ^ s
        candidates = [x for x in D if x & s not in {h & s for h in H}
                      and x & other not in {y & other for y in complement}]
        assert len(candidates) == 1
        markers[s] = candidates[0]
    assert supports == {3, 4, 5, 6, 9, 10}
    assert markers == {3: 1, 9: 3, 10: 2}
    assert any(s != t and s & t == s for s in supports for t in supports)
    assert any((x ^ y).bit_count() == 1 for x in markers.values()
               for y in markers.values() if x != y)
    available = [x for x in sorted(D) if
                 len(shadow(H | {x}, n)) == len(H) + 1 or corner(A, x, n)]
    assert len(available) >= len(supports) - len(isolated)
    out = dict(status='PASS', A=sorted(A), H=sorted(H), designated_vertex=d,
               exclusive_cycle=cycle, relative_deletion_order=deletion,
               first_failure=dict(pair=[1,11], missing_intermediates=[3,9]),
               last_failure=dict(pair=[3,5], missing_intermediates=[1,7]),
               relative_supports=sorted(supports), isolated_markers=markers,
               available_vertices=available,
               nonisolated_support_count=len(supports)-len(isolated),
               scope='Refutes the degree-two vertexwise first-or-last dichotomy only. The pair is relatively peelable and supplies no Damp obstruction or counterexample to arbitrary-cycle relative peeling.')
    Path(__file__).with_suffix('.json').write_text(json.dumps(out, indent=2)+'\n')
    print('PASS: C6 exclusive graph; vertex1 cannot go first or last; relative peeling replayed.')


if __name__ == '__main__':
    main()
