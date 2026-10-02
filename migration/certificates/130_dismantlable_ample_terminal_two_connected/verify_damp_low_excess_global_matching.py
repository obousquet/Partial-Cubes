#!/usr/bin/env python3
"""Verify the global matching step when both wing excesses are at most one.

The local comparison theorems leave only three-pair matchings (the central
square is excluded separately by passive shielding).  This verifier applies
the forced-owner reconstruction to the complete eta-one normal-form universe
and to the two exceptional extra-owner paths.  Every eta-one matching and the
P9 extra-owner matching forces six active coordinates in a residual wing.
The sole five-coordinate residue is the positioned P8 pattern passed to the
five-coordinate completion verifier.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from itertools import combinations
from pathlib import Path

import explore_damp_two_arm_tree_decompositions as M
import verify_damp_pure_eta_two_mixed_profiles as R
import explore_damp_three_pair_global_compatibility as G


ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "damp_low_excess_global_matching.json"
P8_REQUIRED_TYPE = [0, 1, 2, 13, 29]


def marker_key(marker, words):
    return (
        words[marker["core"]],
        tuple(sorted(marker["passive"])),
        tuple(sorted(tuple(sorted(pair)) for pair in marker["supports"])),
    )


def required_exclusive(marker, graph, words, edge_class):
    owners = R.forced_owners(marker, graph, words, edge_class)
    if owners is None:
        return None
    if "shore_words" in marker:
        shores = set(marker["shore_words"])
    else:
        embedded = set(words.values())
        passive_mask = sum(1 << direction for direction in marker["passive"])
        base_bits = words[marker["core"]] & passive_mask
        shores = {
            (word & ~passive_mask) | base_bits for word in embedded
        } - embedded
    required = shores | owners
    active = G.varying_coordinates(required)
    return {
        "forced_owners": sorted(owners),
        "forced_shore_words": sorted(shores),
        "required_exclusive_words": sorted(required),
        "active_coordinates": sorted(active),
        "active_coordinate_count": len(active),
        "required_active_type": list(G.canonical_active_family(required)),
    }


def eta_normal_form_records():
    eta_eta = {}
    eta_extra = {}
    generated_graphs = M.eta_one_graphs()

    for graph, edge_class in generated_graphs:
        words, recovered_edge_class = M.partial_cube_embedding(graph)
        if recovered_edge_class != edge_class:
            raise AssertionError("noncanonical eta-one embedding")
        family = tuple(sorted(words.values()))
        eta_markers = R.tree_arm_markers(
            graph, words, edge_class, 1, len(graph) - 7, "eta1"
        )
        extra_decompositions = M.extra_owner_decompositions(graph)

        for first, second in combinations(eta_markers, 2):
            pairs = set(first["supports"]) | set(second["supports"])
            if len(pairs) != 4 or not R.matching_number_at_least_three(pairs):
                continue
            key = (family,) + tuple(sorted((
                marker_key(first, words), marker_key(second, words)
            )))
            eta_eta[key] = {
                "family": list(family),
                "first": required_exclusive(
                    first, graph, words, edge_class
                ),
                "second": required_exclusive(
                    second, graph, words, edge_class
                ),
                "support_pairs": sorted(
                    tuple(sorted(pair)) for pair in pairs
                ),
            }

        for eta in eta_markers:
            eta_record = required_exclusive(eta, graph, words, edge_class)
            for extra in extra_decompositions:
                ordinary = frozenset(
                    edge_class[edge] for edge in extra["small_support"]
                )
                if not all(
                    ordinary.isdisjoint(pair) for pair in eta["supports"]
                ):
                    continue
                extra_marker = extra_path_marker(
                    graph, words, edge_class, extra
                )
                key = (
                    family,
                    marker_key(eta, words),
                    tuple(sorted(ordinary)),
                )
                eta_extra[key] = {
                    "family": list(family),
                    "eta": eta_record,
                    "extra": required_exclusive(
                        extra_marker, graph, words, edge_class
                    ),
                    "eta_supports": sorted(
                        tuple(sorted(pair)) for pair in eta["supports"]
                    ),
                    "extra_ordinary_support": sorted(ordinary),
                }

    eta_eta_records = list(eta_eta.values())
    eta_extra_records = list(eta_extra.values())
    for record in eta_eta_records:
        for side in ("first", "second"):
            side_record = record[side]
            if side_record is None:
                raise AssertionError("a matching has a distant corner")
            if side_record["active_coordinate_count"] >= 6:
                continue
            if not (
                len(record["family"]) == 8
                and side_record["active_coordinate_count"] == 5
                and side_record["required_active_type"] == P8_REQUIRED_TYPE
            ):
                raise AssertionError("unexpected eta/eta equality residue")
    for record in eta_extra_records:
        for side in ("eta", "extra"):
            side_record = record[side]
            if side_record is None:
                raise AssertionError("a matching has a distant corner")
            if side_record["active_coordinate_count"] >= 6:
                continue
            if not (
                len(record["family"]) == 8
                and side_record["active_coordinate_count"] == 5
                and side_record["required_active_type"] == P8_REQUIRED_TYPE
            ):
                raise AssertionError("unexpected eta/extra equality residue")
    return len(generated_graphs), eta_eta_records, eta_extra_records


def path_graph(order):
    graph = {vertex: set() for vertex in range(order)}
    for vertex in range(order - 1):
        M.add_edge(graph, vertex, vertex + 1)
    return graph


def extra_path_marker(graph, words, edge_class, decomposition):
    core = set(decomposition["core"])
    components = M.components(graph, set(graph) - core)
    large = next(component for component in components if len(component) == 4)
    attachment = next(
        M.canonical_edge(core_vertex, component_vertex)
        for core_vertex in core
        for component_vertex in large
        if component_vertex in graph[core_vertex]
    )
    small_passive = edge_class[decomposition["small_edges"][0]]
    large_passive = edge_class[attachment]
    passive = frozenset((small_passive, large_passive))
    ordinary = frozenset(
        edge_class[edge] for edge in decomposition["small_support"]
    )
    second_edge = set(decomposition["small_edges"][1])
    third_edge = set(decomposition["small_edges"][2])
    middle = next(iter(second_edge & third_edge))
    terminal = next(iter(third_edge - {middle}))
    collar = next(vertex for vertex in attachment if vertex in large)
    shore_words = {
        words[middle] ^ (1 << small_passive),
        words[terminal] ^ (1 << small_passive),
    }
    shore_words.update(
        words[vertex] ^ (1 << large_passive)
        for vertex in large - {collar}
    )
    return {
        "type": "extra-owner",
        "core": min(core),
        "passive": passive,
        "supports": frozenset((ordinary,)),
        "shore_words": frozenset(shore_words),
    }


def extra_path_records():
    records = []
    for order in (8, 9):
        graph = path_graph(order)
        words, edge_class = M.partial_cube_embedding(graph)
        decompositions = M.extra_owner_decompositions(graph)
        if len(decompositions) != 2:
            raise AssertionError((order, len(decompositions)))
        side_records = []
        ordinary_supports = []
        for decomposition in decompositions:
            marker = extra_path_marker(
                graph, words, edge_class, decomposition
            )
            ordinary_supports.append(next(iter(marker["supports"])))
            side_records.append(
                required_exclusive(marker, graph, words, edge_class)
            )
        if ordinary_supports[0] == ordinary_supports[1]:
            raise AssertionError("the two path decompositions repeat an arm")
        records.append({
            "overlap_order": order,
            "family": sorted(words.values()),
            "ordinary_supports": sorted(
                tuple(sorted(pair)) for pair in ordinary_supports
            ),
            "sides": side_records,
        })

    p8, p9 = records
    if any(
        side is None or side["active_coordinate_count"] != 5
        for side in p8["sides"]
    ):
        raise AssertionError("P8 is not the five-coordinate equality residue")
    if any(
        side["required_active_type"] != P8_REQUIRED_TYPE
        for side in p8["sides"]
    ):
        raise AssertionError("P8 has an unexpected forced side type")
    if any(
        side is None or side["active_coordinate_count"] < 6
        for side in p9["sides"]
    ):
        raise AssertionError("P9 does not force six active coordinates")
    return records


def run():
    graph_count, eta_eta, eta_extra = eta_normal_form_records()
    paths = extra_path_records()
    records = {
        "eta_eta": eta_eta,
        "eta_extra": eta_extra,
        "extra_paths": paths,
    }
    five_states = {}

    def retain_five_state(family, side):
        if side is None or side["active_coordinate_count"] != 5:
            return
        state = {
            "family": family,
            "active": side["active_coordinates"],
            "required": side["required_exclusive_words"],
            "required_type": side["required_active_type"],
        }
        key = json.dumps(state, sort_keys=True, separators=(",", ":"))
        five_states[key] = state

    for record in eta_eta:
        for side in ("first", "second"):
            retain_five_state(record["family"], record[side])
    for record in eta_extra:
        for side in ("eta", "extra"):
            retain_five_state(record["family"], record[side])
    for record in paths:
        for side in record["sides"]:
            retain_five_state(record["family"], side)

    canonical = json.dumps(records, sort_keys=True, separators=(",", ":"))
    return {
        "scope": "all global three-pair residues with both excesses at most one",
        "eta_normal_form_graph_count": graph_count,
        "distinct_eta_eta_matching_count": len(eta_eta),
        "distinct_eta_extra_matching_count": len(eta_extra),
        "eta_matchings_force_six_active_coordinates_or_P8": True,
        "extra_owner_P9_forces_six_active_coordinates": True,
        "extra_owner_P8_is_unique_five_coordinate_residue": True,
        "five_coordinate_side_state_count": len(five_states),
        "five_coordinate_side_states": [
            five_states[key] for key in sorted(five_states)
        ],
        "records": records,
        "mathematical_digest": hashlib.sha256(canonical.encode()).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    record = run()
    arguments.output.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        key: value for key, value in record.items() if key != "records"
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
