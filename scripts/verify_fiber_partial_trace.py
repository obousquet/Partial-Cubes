#!/usr/bin/env python3
"""Exact checks for fiber replacement, partial traces, and witness blocks."""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
from itertools import combinations, product
from math import comb
from pathlib import Path
import json

from verify_signed_pattern_recognition import patterns, shattered


ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "data/examples/assets"
RECONCILIATION = ASSETS / "maximum_parent_reconciliation_source.json"
OUTPUT = ASSETS / "fiber_partial_trace_database_verification.json"
EXPECTED_RECONCILIATION = "6769f2952b81b90c808deef4a3c859b8617fb1cd3752eda45a68c6c03d020386"


def vc_dimension(concepts, n: int):
    return max((support.bit_count() for support in shattered(concepts, n)), default=-1)


def is_maximum(concepts, n: int, rank: int):
    return (
        len(concepts) == sum(comb(n, i) for i in range(rank + 1))
        and vc_dimension(concepts, n) == rank
    )


def is_ample(concepts, n: int):
    return len(concepts) == len(shattered(concepts, n))


def weight_sequence(n: int):
    sequence = {-1: frozenset()}
    for rank in range(n + 1):
        sequence[rank] = frozenset(word for word in range(1 << n) if word.bit_count() <= rank)
    return sequence


def fiber_graph(sequence, a: int, k: int):
    layer = sequence[k] - sequence[k - 1]
    if k == 0:
        return layer, {x: {x} for x in layer}
    missing = [p for p in patterns(sequence[k - 1], a) if p[0].bit_count() == k]
    adjacency = {x: {x} for x in layer}
    for left, right in combinations(layer, 2):
        if any((left & support) == bits == (right & support) for support, bits in missing):
            adjacency[left].add(right)
            adjacency[right].add(left)
    components = {}
    unseen = set(layer)
    while unseen:
        root = unseen.pop()
        component = {root}
        frontier = [root]
        while frontier:
            word = frontier.pop()
            for neighbor in adjacency[word] - component:
                component.add(neighbor)
                unseen.discard(neighbor)
                frontier.append(neighbor)
        for word in component:
            components[word] = component
    return layer, components


def build_replacement(sequence_a, a: int, sequence_b, m: int, k: int, s: int, replacements):
    rank = k + s
    layer = sequence_a[k] - sequence_a[k - 1]
    answer = set()
    for x in range(1 << a):
        if x in layer:
            fiber = replacements[x]
        else:
            entry = min(level for level in range(a + 1) if x in sequence_a[level])
            index = rank - entry
            if index < 0:
                fiber = frozenset()
            elif index >= m:
                fiber = frozenset(range(1 << m))
            else:
                fiber = sequence_b[index]
        answer.update(x | (y << a) for y in fiber)
    return frozenset(answer)


def replacement_conditions(sequence_a, a: int, sequence_b, m: int, k: int, s: int, replacements):
    layer, components = fiber_graph(sequence_a, a, k)
    expected_size = sum(comb(m, i) for i in range(s + 1))
    if any(not is_maximum(replacements[x], m, s) or len(replacements[x]) != expected_size for x in layer):
        return False
    if any(replacements[x] != replacements[next(iter(components[x]))] for x in layer):
        return False
    if k > 0 and any(not replacements[x] <= sequence_b[s + 1] for x in layer):
        return False
    if k < a and any(not sequence_b[s - 1] <= replacements[x] for x in layer):
        return False
    return True


def replay_fiber_criterion(counts: Counter):
    # A non-weight rank-one sequence gives a connected changed layer.
    sequence_a = {
        -1: frozenset(),
        0: frozenset({0}),
        1: frozenset({0, 1, 3}),
        2: frozenset(range(4)),
    }
    sequence_b = weight_sequence(2)
    for k in range(3):
        layer, _ = fiber_graph(sequence_a, 2, k)
        for s in range(2):
            size = sum(comb(2, i) for i in range(s + 1))
            candidates = [frozenset(x) for x in combinations(range(4), size)]
            for chosen in product(candidates, repeat=len(layer)):
                replacements = dict(zip(sorted(layer), chosen))
                modified = build_replacement(sequence_a, 2, sequence_b, 2, k, s, replacements)
                conditions = replacement_conditions(sequence_a, 2, sequence_b, 2, k, s, replacements)
                assert is_maximum(modified, 4, k + s) == conditions
                counts["fiber_replacement_assignments"] += 1

    # The source's two-component example.
    old = weight_sequence(2)
    control = weight_sequence(3)
    layer = sorted(old[1] - old[0])
    replacements = {layer[0]: frozenset({0, 1, 3, 5}), layer[1]: frozenset({0, 2, 3, 6})}
    example = build_replacement(old, 2, control, 3, 1, 1, replacements)
    assert replacement_conditions(old, 2, control, 3, 1, 1, replacements)
    assert is_maximum(example, 5, 2) and len(example) == 16
    counts["distinct_fiber_example_words"] = len(example)

    # Complementing a middle fiber fails at both old levels for a=1.
    old_one = weight_sequence(1)
    middle_complement = frozenset(set(range(8)) - set(control[1]))
    for k in (0, 1):
        layer = old_one[k] - old_one[k - 1]
        replacements = {x: middle_complement for x in layer}
        modified = build_replacement(old_one, 1, control, 3, k, 1, replacements)
        assert not replacement_conditions(old_one, 1, control, 3, k, 1, replacements)
        assert not is_maximum(modified, 4, k + 1)
        counts["complementary_fiber_obstructions"] += 1


def partial_trace(concepts, n: int):
    answer = set()
    for mask in range(1 << n):
        specified = ((1 << n) - 1) ^ mask
        for visible in range(1 << n):
            if any((visible & specified) == (word & specified) for word in concepts):
                answer.add(visible | (mask << n))
    return frozenset(answer)


def replay_partial_traces(counts: Counter):
    n = 3
    for family_mask in range(1, 1 << (1 << n)):
        concepts = frozenset(word for word in range(1 << n) if family_mask >> word & 1)
        expanded = partial_trace(concepts, n)
        sh = shattered(concepts, n)
        predicted = {
            support
            for support in range(1 << (2 * n))
            if sum(
                (1 << i)
                for i in range(n)
                if (support >> i) & 1 and (support >> (n + i)) & 1
            ) in sh
        }
        assert shattered(expanded, 2 * n) == predicted
        assert vc_dimension(expanded, 2 * n) == n + vc_dimension(concepts, n)
        zero_section = {word & ((1 << n) - 1) for word in expanded if word >> n == 0}
        full_reduction = {
            visible
            for visible in range(1 << n)
            if all((visible | (mask << n)) in expanded for mask in range(1 << n))
        }
        assert zero_section == set(concepts) == full_reduction
        assert is_ample(expanded, 2 * n) == is_ample(concepts, n)
        counts["partial_trace_classes_q3"] += 1

    singleton = partial_trace({0}, 1)
    assert singleton == {0, 2, 3}
    assert is_maximum(singleton, 2, 1)
    star = partial_trace({0, 1, 2, 4}, 3)
    assert len(star) == 54
    assert vc_dimension(star, 6) == 4
    assert sum(comb(6, i) for i in range(5)) == 57
    assert is_ample(star, 6) and not is_maximum(star, 6, 4)
    counts["partial_trace_star_words"] = len(star)


def replay_disjoint_witness_example(counts: Counter):
    # Visible coordinates 0,1; auxiliary blocks {2},{3}. The path flips
    # b1,b2,x1,x2, so the zero-marked section is the visible singleton.
    parent = frozenset({0, 4, 12, 13, 15})
    assert is_maximum(parent, 4, 1)
    section = {word & 3 for word in parent if word & 12 == 0}
    assert section == {0}
    missing = patterns(parent, 4)
    assert (5, 1) in missing   # x1=1, b1=0
    assert (10, 2) in missing  # x2=1, b2=0
    visible_patterns = patterns({0}, 2)
    assert visible_patterns == {(1, 1), (2, 2)}
    block_total = sum(2 - support.bit_count() for support, _ in visible_patterns)
    assert block_total == 2
    counts["disjoint_witness_parent_words"] = len(parent)
    counts["disjoint_witness_blocks"] = block_total


def main():
    assert sha256(RECONCILIATION.read_bytes()).hexdigest() == EXPECTED_RECONCILIATION
    reconciliation = json.loads(RECONCILIATION.read_text())
    required = {
        "thm:fibers",
        "thm:partial-traces",
        "thm:disjoint-witnesses",
        "ex:fibers",
        "cor:complement-fiber",
    }
    entries = {entry["label"]: entry for entry in reconciliation["statements"]}
    assert required <= entries.keys()
    assert all(entries[label]["reconciliation"] == "scoped_audit_located" for label in required)

    counts = Counter()
    replay_fiber_criterion(counts)
    replay_partial_traces(counts)
    replay_disjoint_witness_example(counts)
    output = {
        "date": "2026-10-01",
        "status": "passed",
        "arithmetic": "exhaustive finite binary families, exact shattering sets, and explicit maximum parents",
        "scope": "Exhausts the fiber criterion for a connected two-word changed layer across all old and control levels, checks the distinct-fiber and complementary-fiber examples, verifies the full partial-trace theorem for every nonempty Q3 class, and replays an explicit disjoint-witness parent. Universal parent-rank and coordinate-copy constructions use database-owned proofs.",
        "source_sha256": {RECONCILIATION.name: EXPECTED_RECONCILIATION},
        "source_statement_labels": sorted(required),
        "counts": dict(counts),
    }
    OUTPUT.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"status": "passed", "counts": dict(counts)}))


if __name__ == "__main__":
    main()
