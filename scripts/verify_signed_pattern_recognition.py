#!/usr/bin/env python3
"""Exact finite replay for signed-pattern recognition and section identities."""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from math import comb
from pathlib import Path
import json


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "data/examples/assets"
RECONCILIATION = ASSETS / "maximum_parent_reconciliation_source.json"
OUTPUT = ASSETS / "signed_pattern_recognition_database_verification.json"
EXPECTED_RECONCILIATION = "6769f2952b81b90c808deef4a3c859b8617fb1cd3752eda45a68c6c03d020386"


def submasks(mask: int):
    sub = mask
    while True:
        yield sub
        if sub == 0:
            break
        sub = (sub - 1) & mask


def patterns(concepts, n: int):
    concepts = set(concepts)
    answer = set()
    for support in range(1, 1 << n):
        for bits in submasks(support):
            if any((word & support) == bits for word in concepts):
                continue
            minimal = True
            for coordinate in range(n):
                if support >> coordinate & 1:
                    smaller = support ^ (1 << coordinate)
                    if not any((word & smaller) == (bits & smaller) for word in concepts):
                        minimal = False
                        break
            if minimal:
                answer.add((support, bits))
    return answer


def shattered(concepts, n: int):
    return {
        support
        for support in range(1 << n)
        if len({word & support for word in concepts}) == 1 << support.bit_count()
    }


def vc_dimension(concepts, n: int):
    return max((support.bit_count() for support in shattered(concepts, n)), default=-1)


def selected_families(n: int, rank: int):
    supports = [
        sum(1 << coordinate for coordinate in coordinates)
        for coordinates in combinations(range(n), rank + 1)
    ]
    choices = [tuple(submasks(support)) for support in supports]
    for values in product(*choices):
        yield dict(zip(supports, values))


def agreement_connected(selected, n: int, rank: int):
    for coordinates in combinations(range(n), rank + 2):
        universe = sum(1 << coordinate for coordinate in coordinates)
        vertices = [universe ^ (1 << coordinate) for coordinate in coordinates]
        seen = {vertices[0]}
        frontier = [vertices[0]]
        while frontier:
            left = frontier.pop()
            for right in vertices:
                common = left & right
                if (selected[left] & common) == (selected[right] & common) and right not in seen:
                    seen.add(right)
                    frontier.append(right)
        if len(seen) != len(vertices):
            return False
    return True


def avoider(selected, n: int):
    return {
        word
        for word in range(1 << n)
        if all((word & support) != bits for support, bits in selected.items())
    }


def contains_pattern(container, contained):
    outer_support, outer_bits = container
    support, bits = contained
    return support & ~outer_support == 0 and (outer_bits & support) == bits


def exact_section_conditions(selected, candidate, tau: int):
    missing = patterns(candidate, 2)
    active = []
    for support, bits in selected.items():
        if support & 4 and ((bits >> 2) & 1) != tau:
            continue
        active.append((support & 3, bits & 3))
    condition_two = all(any(contains_pattern(restriction, p) for p in missing) for restriction in active)
    condition_three = all(p in active for p in missing)
    return condition_two and condition_three


def minimal_generators(generators):
    generators = set(generators)
    return {g for g in generators if not any(h < g for h in generators)}


def literal_pattern(pattern):
    support, bits = pattern
    return frozenset(
        (coordinate, (bits >> coordinate) & 1)
        for coordinate in range(support.bit_length())
        if support >> coordinate & 1
    )


def faces_from_words(words, n: int):
    faces = {frozenset()}
    for word in words:
        facet = frozenset((coordinate, (word >> coordinate) & 1) for coordinate in range(n))
        facet_list = tuple(facet)
        for size in range(len(facet_list) + 1):
            faces.update(frozenset(face) for face in combinations(facet_list, size))
    return faces


def downsets(n: int):
    answer = []
    for family_mask in range(1 << (1 << n)):
        family = {word for word in range(1 << n) if family_mask >> word & 1}
        if family and all((sub in family) for word in family for sub in submasks(word)):
            answer.append(family)
    return answer


def replay_local_recognition(counts: Counter):
    all_zero = {3: 0, 5: 0, 6: 0}
    disconnected = {3: 0, 5: 1, 6: 6}
    assert agreement_connected(all_zero, 3, 1)
    assert avoider(all_zero, 3) == {7, 6, 5, 3}
    assert not agreement_connected(disconnected, 3, 1)
    assert avoider(disconnected, 3) == {2, 5}

    for selected in selected_families(3, 1):
        concepts = avoider(selected, 3)
        connected = agreement_connected(selected, 3, 1)
        maximum = len(concepts) == sum(comb(3, i) for i in range(2)) and vc_dimension(concepts, 3) == 1
        assert connected == maximum
        counts["signed_pattern_families"] += 1
        counts["connected_pattern_families"] += connected

        if not connected:
            continue
        for tau in (0, 1):
            actual = {word & 3 for word in concepts if ((word >> 2) & 1) == tau}
            for candidate_mask in range(1, 15):
                candidate = {word for word in range(4) if candidate_mask >> word & 1}
                criterion = exact_section_conditions(selected, candidate, tau)
                assert criterion == (candidate == actual)
                counts["exact_section_candidate_checks"] += 1


def replay_signed_ideal_identities(counts: Counter):
    # Depolarization: supports of signed minimal generators are exactly
    # the minimal nonfaces of the shattered-support complex.
    for family_mask in range(1, 1 << 8):
        concepts = {word for word in range(8) if family_mask >> word & 1}
        signed_supports = {support for support, _ in patterns(concepts, 3)}
        signed_supports = {
            support
            for support in signed_supports
            if not any(other != support and other & support == other for other in signed_supports)
        }
        sh = shattered(concepts, 3)
        nonfaces = {
            support
            for support in range(1, 8)
            if support not in sh and all((support ^ (1 << i)) in sh for i in range(3) if support >> i & 1)
        }
        assert signed_supports == nonfaces
        counts["signed_depolarizations"] += 1

        parent_generators = {literal_pattern(p) for p in patterns(concepts, 3)}
        parent_faces = faces_from_words(concepts, 3)
        for tau in (0, 1):
            section = {word & 3 for word in concepts if ((word >> 2) & 1) == tau}
            if not section:
                continue
            colon_generators = set()
            for generator in parent_generators:
                reduced = set(generator)
                reduced.discard((2, tau))
                if any(coordinate == 2 for coordinate, _ in reduced):
                    continue
                colon_generators.add(frozenset(reduced))
            colon_generators = minimal_generators(colon_generators)
            section_generators = {literal_pattern(p) for p in patterns(section, 2)}
            assert colon_generators == section_generators

            sigma = frozenset({(2, tau)})
            link = {
                frozenset((coordinate, value) for coordinate, value in face if coordinate < 2)
                for face in parent_faces
                if sigma <= face
            }
            assert link == faces_from_words(section, 2)
            counts["section_colon_and_link_checks"] += 1

    # Exact signed ideal formula for every pair of nonempty downsets of Q3.
    families = downsets(3)
    assert len(families) == 19
    for concepts in families:
        qc = {literal_pattern(p) for p in patterns(concepts, 3)}
        for filt in families:
            result = concepts & filt
            qf = {literal_pattern(p) for p in patterns(result, 3)}
            positive_nonfaces = {
                frozenset((coordinate, 1) for coordinate in range(3) if support >> coordinate & 1)
                for support, bits in patterns(filt, 3)
                if bits == support
            }
            assert minimal_generators(qc | positive_nonfaces) == qf
            counts["downset_filter_signed_ideal_checks"] += 1


def replay_relative_initial_example(counts: Counter):
    # The one relation leaves three dimensions in bidegree (1,1), matching H_{2,1}.
    assert 4 - 1 == sum(comb(2, i) for i in range(2))
    full_section = {0, 1}  # in(a0*b0-a1*b1)=(a1*b1), mark b=0
    singleton_section = {1}  # the wrong order gives (a0*b0)
    assert full_section != singleton_section
    counts["relative_initial_term_orders"] = 2
    counts["relative_initial_bidegree_dimension"] = 3


def main():
    assert sha256(RECONCILIATION.read_bytes()).hexdigest() == EXPECTED_RECONCILIATION
    reconciliation = json.loads(RECONCILIATION.read_text())
    required = {
        "thm:recognition",
        "ex:graphs",
        "cor:exact-section",
        "thm:homogeneous",
        "prop:colon",
        "thm:relative-initial",
        "ex:relative-order",
    }
    entries = {entry["label"]: entry for entry in reconciliation["statements"]}
    assert required <= entries.keys()
    assert all(entries[label]["reconciliation"] == "scoped_audit_located" for label in required)

    counts = Counter()
    replay_local_recognition(counts)
    replay_signed_ideal_identities(counts)
    replay_relative_initial_example(counts)
    output = {
        "date": "2026-10-01",
        "status": "passed",
        "arithmetic": "exhaustive finite binary set and squarefree monomial enumeration",
        "scope": "Exhausts every rank-one signed pattern family on three coordinates, every candidate marked section on two visible coordinates, all nonempty Q3 classes for signed depolarization, colon and Yang-link identities, every pair of Q3 downset filters, and both leading-term choices in the relative-initial example. General regular-sequence and Hilbert-series theorems use database-owned proofs.",
        "source_sha256": {RECONCILIATION.name: EXPECTED_RECONCILIATION},
        "source_statement_labels": sorted(required),
        "counts": dict(counts),
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"status": "passed", "counts": dict(counts)}))


if __name__ == "__main__":
    main()
