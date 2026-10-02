#!/usr/bin/env python3
"""Verify the complete order-nine two-unit opposite-wing comparison.

This is a normal-form audit, not a cardinality census.  Every feasible
two-unit wing over an overlap of order at most nine has one of five forms:

* the pure eta=2 form, obtained by adjoining two vertices to P7;
* the disjoint eta=1 (2,3)-shore form, obtained by adjoining one vertex
  to a one-vertex wedge of rooted ample gates of orders three and four;
* the vertex-free disjoint (2,4)- and (3,3)-shore wedges; or
* the vertex-free singleton-extra-shore wedge.

The shared-owner eta=1 (2,2)-shore form and the overlapping vertex-free
form are impossible by the structural results in the manuscript.

The verifier constructs the complete finite universe furnished by these
normal forms, rediscovers all exact two-unit markers, and compares each of
them with every sharp, uncounted one-unit, two-unit, and pair-rich
extra-owner marker in the same overlap.  Every comparison is resolved by
a repeated acquired support, passive shielding, a forced-owner conflict,
or a matching of three acquired pairs.  The pair-poor extra-owner branch
is handled separately by the rooted-claw theorem in the manuscript.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from collections import Counter
from itertools import combinations, permutations
from pathlib import Path


ROOT = Path(__file__).resolve().parent
MIXED_VERIFIER = ROOT / "verify_damp_pure_eta_two_mixed_profiles.py"


def load_mixed_verifier():
    spec = importlib.util.spec_from_file_location(
        "pure_eta_two_mixed_verifier", MIXED_VERIFIER
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {MIXED_VERIFIER}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


R = load_mixed_verifier()
M = R.M
V = R.V


def graph_from_edges(order, edges):
    graph = {vertex: set() for vertex in range(order)}
    for first, second in edges:
        M.add_edge(graph, first, second)
    return graph


def rooted_graph_key(graph, root):
    """Canonical adjacency key for a rooted graph of order at most five."""

    other_vertices = [vertex for vertex in graph if vertex != root]
    best = None
    for order in permutations(other_vertices):
        ordered_vertices = (root,) + order
        candidate = tuple(
            int(ordered_vertices[second] in graph[ordered_vertices[first]])
            for first in range(len(ordered_vertices))
            for second in range(first + 1, len(ordered_vertices))
        )
        if best is None or candidate < best:
            best = candidate
    assert best is not None
    return best


def rooted_ample_gates(order):
    """Return all rooted connected ample partial-cube graphs of given order."""

    possible_edges = tuple(combinations(range(order), 2))
    output = {}
    for edge_mask in range(1 << len(possible_edges)):
        graph = graph_from_edges(
            order,
            (
                edge
                for index, edge in enumerate(possible_edges)
                if (edge_mask >> index) & 1
            ),
        )
        if not M.connected(graph, set(graph)):
            continue
        embedding = M.partial_cube_embedding(graph)
        if embedding is None:
            continue
        words, edge_class = embedding
        dimension = 1 + max(edge_class.values(), default=-1)
        if not M.is_ample_embedding(words, dimension):
            continue
        for root in graph:
            output.setdefault(
                rooted_graph_key(graph, root), (graph, root)
            )
    return list(output.values())


def wedge_components(components):
    """Attach disjoint rooted active gates to one passive base vertex."""

    graph = {0: set()}
    next_vertex = 1
    for component, root in components:
        relabel = {
            vertex: next_vertex + index
            for index, vertex in enumerate(component)
        }
        next_vertex += len(component)
        for vertex in component:
            graph[relabel[vertex]] = set()
        for first, second in M.graph_edges(component):
            M.add_edge(graph, relabel[first], relabel[second])
        M.add_edge(graph, 0, relabel[root])
    return graph


def ample_embedding(graph):
    if not M.connected(graph, set(graph)):
        return None
    embedding = M.partial_cube_embedding(graph)
    if embedding is None:
        return None
    words, edge_class = embedding
    dimension = 1 + max(edge_class.values(), default=-1)
    if not M.is_ample_embedding(words, dimension):
        return None
    return words, edge_class


def exact_profile_universe():
    """Construct every overlap supplied by a feasible two-unit normal form."""

    path_two = (
        graph_from_edges(2, ((0, 1),)),
        0,
    )
    path_three = (
        graph_from_edges(3, ((0, 1), (1, 2))),
        0,
    )
    rooted_four = rooted_ample_gates(4)
    rooted_five = rooted_ample_gates(5)

    families = {}
    source_occurrences = Counter()

    def retain(graph, expected_profile):
        embedding = ample_embedding(graph)
        if embedding is None:
            return
        words, edge_class = embedding
        if expected_profile == "eta2-disjoint-2-2":
            if not V.pure_decompositions(graph, edge_class):
                return
        elif not any(
            marker["type"] == expected_profile
            for marker in R.exact_two_unit_markers(words, edge_class)
        ):
            return
        source_occurrences[expected_profile] += 1
        family = V.canonical_family(words)
        families.setdefault(family, (graph, words, edge_class))

    pure_spine = M.build_tree(1, (0, 0))
    for first_neighbours in V.subsets(pure_spine):
        first_extension = M.add_vertex(pure_spine, first_neighbours)
        for second_neighbours in V.subsets(first_extension):
            retain(
                M.add_vertex(first_extension, second_neighbours),
                "eta2-disjoint-2-2",
            )

    for large_gate in rooted_four:
        base = wedge_components([path_three, large_gate])
        for neighbours in V.subsets(base):
            retain(
                M.add_vertex(base, neighbours),
                "eta1-disjoint-2-3",
            )

    for large_gate in rooted_five:
        retain(
            wedge_components([path_three, large_gate]),
            "eta0-disjoint-2-4",
        )

    for first_index, first_gate in enumerate(rooted_four):
        for second_gate in rooted_four[first_index:]:
            retain(
                wedge_components([first_gate, second_gate]),
                "eta0-disjoint-3-3",
            )

    retain(
        wedge_components([path_three, path_three, path_two]),
        "eta0-extra-singleton",
    )

    return {
        "rooted_ample_gate_counts": {
            "order_4": len(rooted_four),
            "order_5": len(rooted_five),
        },
        "source_occurrences": dict(sorted(source_occurrences.items())),
        "families": families,
    }


def extra_owner_markers(graph, words, edge_class):
    """Recognize every pair-rich order-nine extra-owner decomposition."""

    all_edges = set(M.graph_edges(graph))
    markers = {}
    pair_rich_decompositions = set()
    pair_poor_decompositions = set()

    for core_tuple in combinations(graph, 2):
        core = set(core_tuple)
        if core_tuple[1] not in graph[core_tuple[0]]:
            continue
        exterior_components = M.components(graph, set(graph) - core)
        if sorted(map(len, exterior_components)) != [3, 4]:
            continue

        attachment_data = {}
        valid = True
        for component in exterior_components:
            attachment_edges = [
                (core_vertex, component_vertex)
                for core_vertex in core
                for component_vertex in component
                if component_vertex in graph[core_vertex]
            ]
            if len(attachment_edges) != 1:
                valid = False
                break
            attachment_data[len(component)] = (
                component,
                attachment_edges[0],
            )
        if not valid:
            continue

        small_component, (small_core, first) = attachment_data[3]
        small_degrees = {
            vertex: len(graph[vertex] & small_component)
            for vertex in small_component
        }
        if (
            sorted(small_degrees.values()) != [1, 1, 2]
            or small_degrees[first] != 1
        ):
            continue
        middle = next(iter(graph[first] & small_component))
        terminal = next(
            iter((graph[middle] & small_component) - {first})
        )
        small_edges = (
            M.canonical_edge(small_core, first),
            M.canonical_edge(first, middle),
            M.canonical_edge(middle, terminal),
        )
        small_directions = tuple(
            edge_class[edge] for edge in small_edges
        )
        small_passive = small_directions[0]
        ordinary_support = frozenset(
            (small_passive, small_directions[2])
        )

        large_component, (large_core, collar) = attachment_data[4]
        large_attachment = M.canonical_edge(large_core, collar)
        large_passive = edge_class[large_attachment]
        if small_passive == large_passive:
            continue
        if (
            sum(
                edge_class[edge] == small_passive for edge in all_edges
            )
            != 1
            or sum(
                edge_class[edge] == large_passive for edge in all_edges
            )
            != 1
        ):
            continue

        large_embedding = ample_embedding(
            {
                vertex: graph[vertex] & large_component
                for vertex in large_component
            }
        )
        if large_embedding is None:
            continue

        large_shore = set(large_component) - {collar}
        shore_edges = M.graph_edges(
            {
                vertex: graph[vertex] & large_shore
                for vertex in large_shore
            }
        )
        decomposition_key = (
            frozenset(core),
            frozenset(small_component),
            frozenset(large_component),
            small_passive,
            large_passive,
            ordinary_support,
        )
        if not shore_edges:
            pair_poor_decompositions.add(decomposition_key)
            continue

        pair_rich_decompositions.add(decomposition_key)
        collar_supports = {
            frozenset((small_passive, small_directions[1]))
        }
        collar_supports.update(
            frozenset(
                (
                    large_passive,
                    edge_class[M.canonical_edge(collar, neighbour)],
                )
            )
            for neighbour in graph[collar] & large_shore
        )
        collar_supports = frozenset(
            support
            for support in collar_supports
            if not V.shatters_pair(words, support)
        )

        for shore_edge in shore_edges:
            large_support = frozenset(
                (large_passive, edge_class[shore_edge])
            )
            if (
                V.shatters_pair(words, large_support)
                or not ordinary_support.isdisjoint(large_support)
            ):
                continue
            marker = {
                "type": "eta0-extra-owner",
                "core": min(core),
                "passive": frozenset(
                    (small_passive, large_passive)
                ),
                "supports": frozenset(
                    (ordinary_support, large_support)
                ),
                "collars": collar_supports,
            }
            passive_mask = (1 << small_passive) | (1 << large_passive)
            base_bits = words[min(core)] & passive_mask
            embedded_words = set(words.values())
            marker["shore_words"] = frozenset(
                normalized
                for word in embedded_words
                for normalized in ((word & ~passive_mask) | base_bits,)
                if normalized not in embedded_words
            )
            key = (
                decomposition_key,
                marker["supports"],
                marker["collars"],
            )
            markers[key] = marker

    return (
        list(markers.values()),
        len(pair_rich_decompositions),
        len(pair_poor_decompositions),
    )


def marker_signature(marker, words):
    return (
        marker["type"],
        words[marker["core"]],
        tuple(sorted(marker["passive"])),
        tuple(
            sorted(
                tuple(sorted(support))
                for support in marker["supports"]
            )
        ),
        tuple(
            sorted(
                tuple(sorted(support))
                for support in marker.get("collars", frozenset())
            )
        ),
    )


def run_audit():
    universe = exact_profile_universe()
    families = universe.pop("families")
    marker_counts = Counter()
    outcome_counts = Counter()
    comparison_keys = set()
    records = []
    extra_owner_decomposition_counts = Counter()

    def compare(
        comparison_class,
        family,
        first,
        second,
        graph,
        words,
        edge_class,
        symmetric=False,
    ):
        first_signature = marker_signature(first, words)
        second_signature = marker_signature(second, words)
        signatures = (first_signature, second_signature)
        if symmetric and second_signature < first_signature:
            signatures = (second_signature, first_signature)
        key = (family, comparison_class, signatures)
        if key in comparison_keys:
            return
        comparison_keys.add(key)
        outcome = R.classify_comparison(
            first, second, graph, words, edge_class
        )
        outcome_counts[(comparison_class, outcome)] += 1
        records.append(
            {
                "family": family,
                "comparison_class": comparison_class,
                "first": signatures[0],
                "second": signatures[1],
                "outcome": outcome,
            }
        )

    for family, (graph, words, edge_class) in families.items():
        two_unit_markers = R.exact_two_unit_markers(words, edge_class)
        sharp_markers = R.tree_arm_markers(
            graph, words, edge_class, 0, 3, "sharp"
        )
        eta_one_markers = R.tree_arm_markers(
            graph, words, edge_class, 1, 2, "eta1"
        )
        for marker in sharp_markers + eta_one_markers:
            marker["collars"] = R.acquired_collar_supports(
                marker, words
            )
        (
            extra_markers,
            pair_rich_count,
            pair_poor_count,
        ) = extra_owner_markers(graph, words, edge_class)
        extra_owner_decomposition_counts["pair_rich"] += pair_rich_count
        extra_owner_decomposition_counts["pair_poor"] += pair_poor_count

        for marker in (
            two_unit_markers
            + sharp_markers
            + eta_one_markers
            + extra_markers
        ):
            marker_counts[marker["type"]] += 1

        for first, second in combinations(two_unit_markers, 2):
            compare(
                "two-unit | two-unit",
                family,
                first,
                second,
                graph,
                words,
                edge_class,
                symmetric=True,
            )
        for first in two_unit_markers:
            for second in sharp_markers:
                compare(
                    "two-unit | sharp",
                    family,
                    first,
                    second,
                    graph,
                    words,
                    edge_class,
                )
            for second in eta_one_markers:
                compare(
                    "two-unit | eta1-one-unit",
                    family,
                    first,
                    second,
                    graph,
                    words,
                    edge_class,
                )
            for second in extra_markers:
                compare(
                    "two-unit | pair-rich-extra-owner",
                    family,
                    first,
                    second,
                    graph,
                    words,
                    edge_class,
                )

    unresolved = sum(
        count
        for (_comparison_class, outcome), count in outcome_counts.items()
        if outcome.endswith(":unresolved")
    )
    if unresolved:
        raise AssertionError(
            f"{unresolved} all-profile comparisons remain"
        )

    records.sort(
        key=lambda record: json.dumps(record, separators=(",", ":"))
    )
    digest = hashlib.sha256(
        json.dumps(records, separators=(",", ":")).encode()
    ).hexdigest()
    return {
        "scope": (
            "all order-nine comparisons with one feasible two-unit wing"
        ),
        "distinct_embedded_overlap_families": len(families),
        **universe,
        "marker_counts": dict(sorted(marker_counts.items())),
        "extra_owner_decomposition_counts": dict(
            sorted(extra_owner_decomposition_counts.items())
        ),
        "comparison_counts": {
            f"{comparison_class} | {outcome}": count
            for (comparison_class, outcome), count in sorted(
                outcome_counts.items()
            )
        },
        "coincident_marker_rule": (
            "Opposite wings using the same marker repeat both acquired "
            "terminal supports and are excluded by shadow intersection."
        ),
        "pair_poor_extra_owner_rule": (
            "The rooted-claw comparison theorem excludes every pair-poor "
            "extra-owner wing opposite a two-unit wing."
        ),
        "unresolved_comparisons": unresolved,
        "mathematical_digest": digest,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "damp_two_unit_all_opposite_profiles.json",
    )
    args = parser.parse_args()
    result = run_audit()
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
