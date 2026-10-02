#!/usr/bin/env python3
"""Verify every mixed comparison having one pure eta=2 two-apex wing.

The seven-profile excess identity gives a coordinate-level normal form for
the overlap B.  Once the two or three passive coordinates, their common
word, and the uncounted remainder are chosen, the remaining vertices are
exactly the collar vertices and the nonempty fibres of the shore-incidence
sets.  This verifier recognizes that identity directly inside every ample
two-vertex completion of the pure P7 spine.

For every opposite sharp, eta=1, or feasible eta=2 marker, the verifier
compares the two acquired terminal pairs.  A no-matching comparison is
resolved by one of four exact outcomes: a passive-distance-two corner, a
forced exclusive-owner collision, a repeated acquired collar, or a matching
of order three after acquired collars are included.  No mixed profile
survives.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from collections import Counter
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PURE_VERIFIER = ROOT / "verify_damp_pure_eta_two_connected_incidence.py"


def load_pure_verifier():
    spec = importlib.util.spec_from_file_location(
        "pure_eta_two_verifier", PURE_VERIFIER
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {PURE_VERIFIER}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


V = load_pure_verifier()
M = V.M


# name, passive-coordinate count, eta, sorted shore sizes, intersection size
TWO_UNIT_PROFILES = (
    ("eta2-disjoint-2-2", 2, 2, (2, 2), 0),
    ("eta1-overlap-2-2", 2, 1, (2, 2), 1),
    ("eta1-disjoint-2-3", 2, 1, (2, 3), 0),
    ("eta0-disjoint-2-4", 2, 0, (2, 4), 0),
    ("eta0-disjoint-3-3", 2, 0, (3, 3), 0),
    ("eta0-extra-singleton", 3, 0, (1, 2, 2), 0),
)


def is_ample_word_family(family, coordinates):
    coordinates = tuple(coordinates)
    shattered_count = 0
    for size in range(len(coordinates) + 1):
        for support in combinations(coordinates, size):
            mask = sum(1 << direction for direction in support)
            if len({word & mask for word in family}) == 1 << size:
                shattered_count += 1
    return shattered_count == len(family)


def active_edges(family):
    output = []
    for first, second in combinations(family, 2):
        difference = first ^ second
        if difference.bit_count() == 1:
            output.append((first, second, difference.bit_length() - 1))
    return output


def matching_number_at_least_three(pairs):
    return any(
        first.isdisjoint(second)
        and first.isdisjoint(third)
        and second.isdisjoint(third)
        for first, second, third in combinations(pairs, 3)
    )


def acquired_collar_supports(marker, words):
    output = set()
    for arm in marker["arms"]:
        support = frozenset(
            (arm["directions"][0], arm["directions"][1])
        )
        if not V.shatters_pair(words, support):
            output.add(support)
    return frozenset(output)


def tree_arm_markers(
    graph, words, edge_class, delete_count, core_order, name
):
    output = {}
    for deleted_tuple in combinations(graph, delete_count):
        deleted = frozenset(deleted_tuple)
        tree = {
            vertex: neighbours - deleted
            for vertex, neighbours in graph.items()
            if vertex not in deleted
        }
        if (
            len(M.graph_edges(tree)) != len(tree) - 1
            or not M.connected(tree, set(tree))
        ):
            continue
        for decomposition in M.decompositions(tree):
            if len(decomposition["core"]) != core_order:
                continue
            arms = []
            valid = True
            for arm in decomposition["arms"]:
                vertices = arm["vertices"]
                path_edges = tuple(
                    M.canonical_edge(vertices[index], vertices[index + 1])
                    for index in range(3)
                )
                directions = tuple(
                    edge_class[edge] for edge in path_edges
                )
                passive, _, terminal = directions
                support = frozenset((passive, terminal))
                other_edges = set(M.graph_edges(tree)) - {path_edges[0]}
                if any(
                    edge_class[edge] == passive for edge in other_edges
                ) or V.shatters_pair(words, support):
                    valid = False
                    break
                arms.append(
                    {
                        "vertices": tuple(vertices),
                        "directions": directions,
                        "support": support,
                        "shore_words": frozenset((
                            words[vertices[2]] ^ (1 << passive),
                            words[vertices[3]] ^ (1 << passive),
                        )),
                    }
                )
            if not valid:
                continue
            supports = frozenset(arm["support"] for arm in arms)
            if len(supports) != 2:
                continue
            first_support, second_support = tuple(supports)
            if not first_support.isdisjoint(second_support):
                continue
            marker = {
                "type": name,
                "core": min(decomposition["core"]),
                "passive": frozenset(
                    arm["directions"][0] for arm in arms
                ),
                "supports": supports,
                "arms": tuple(arms),
                "shore_words": frozenset().union(*(
                    arm["shore_words"] for arm in arms
                )),
            }
            key = (
                marker["type"],
                frozenset(decomposition["core"]),
                deleted,
                marker["passive"],
                supports,
            )
            output[key] = marker
    return list(output.values())


def exact_two_unit_markers(words, edge_class):
    dimension = 1 + max(edge_class.values(), default=-1)
    all_coordinate_mask = (1 << dimension) - 1
    embedded_words = set(words.values())
    markers = {}

    for (
        profile_name,
        passive_count,
        eta,
        shore_sizes,
        intersection_size,
    ) in TWO_UNIT_PROFILES:
        for passive_tuple in combinations(range(dimension), passive_count):
            passive = frozenset(passive_tuple)
            passive_mask = sum(1 << direction for direction in passive)
            active_mask = all_coordinate_mask ^ passive_mask
            active_coordinates = tuple(
                direction
                for direction in range(dimension)
                if direction not in passive
            )

            for core, core_word in words.items():
                base_slice = {
                    word
                    for word in embedded_words
                    if ((word ^ core_word) & passive_mask) == 0
                }
                if base_slice != {core_word}:
                    continue

                for uncounted_tuple in combinations(
                    embedded_words - {core_word}, eta
                ):
                    counted_words = embedded_words - set(uncounted_tuple)
                    by_active_word = {}
                    for word in counted_words:
                        active_word = word & active_mask
                        passive_flip = (word ^ core_word) & passive_mask
                        by_active_word.setdefault(active_word, set()).add(
                            passive_flip
                        )

                    core_active_word = core_word & active_mask
                    collar_flips = {
                        1 << direction for direction in passive
                    }
                    if by_active_word.get(core_active_word) != (
                        {0} | collar_flips
                    ):
                        continue

                    shores = {direction: set() for direction in passive}
                    valid = True
                    for active_word, observed_flips in by_active_word.items():
                        if active_word == core_active_word:
                            continue
                        if 0 in observed_flips:
                            valid = False
                            break
                        incidence_mask = 0
                        for flip in observed_flips:
                            incidence_mask |= flip
                        expected_flips = {
                            flip
                            for flip in range(1, 1 << dimension)
                            if (flip & ~incidence_mask) == 0
                        }
                        if (
                            observed_flips != expected_flips
                            or (incidence_mask & ~passive_mask)
                        ):
                            valid = False
                            break
                        for direction in passive:
                            if (incidence_mask >> direction) & 1:
                                shores[direction].add(active_word)
                    if not valid:
                        continue

                    if any(
                        not is_ample_word_family(
                            {core_active_word} | shores[direction],
                            active_coordinates,
                        )
                        for direction in passive
                    ):
                        continue

                    cycle_assignments = []
                    if passive_count == 2:
                        first, second = passive_tuple
                        if (
                            tuple(
                                sorted(
                                    (len(shores[first]), len(shores[second]))
                                )
                            )
                            != tuple(sorted(shore_sizes))
                            or len(shores[first] & shores[second])
                            != intersection_size
                        ):
                            continue
                        cycle_assignments.append(passive_tuple)
                    else:
                        if (
                            tuple(
                                sorted(
                                    len(shores[direction])
                                    for direction in passive_tuple
                                )
                            )
                            != shore_sizes
                            or any(
                                shores[first] & shores[second]
                                for first, second in combinations(
                                    passive_tuple, 2
                                )
                            )
                        ):
                            continue
                        extra_directions = [
                            direction
                            for direction in passive_tuple
                            if len(shores[direction]) == 1
                        ]
                        if len(extra_directions) != 1:
                            continue
                        cycle_assignments.append(
                            tuple(
                                direction
                                for direction in passive_tuple
                                if direction != extra_directions[0]
                            )
                        )

                    for cycle_directions in cycle_assignments:
                        lower_roots = {
                            direction: {
                                shore_word
                                for shore_word in shores[direction]
                                if (
                                    core_active_word ^ shore_word
                                ).bit_count()
                                == 1
                            }
                            for direction in passive
                        }
                        if any(
                            not lower_roots[direction]
                            or lower_roots[direction] == shores[direction]
                            for direction in cycle_directions
                        ):
                            continue
                        shared_cycle_words = (
                            shores[cycle_directions[0]]
                            & shores[cycle_directions[1]]
                        )
                        if any(
                            shared_cycle_words & lower_roots[direction]
                            for direction in cycle_directions
                        ):
                            continue
                        if passive_count == 3:
                            extra_direction = next(
                                direction
                                for direction in passive
                                if direction not in cycle_directions
                            )
                            if not lower_roots[extra_direction]:
                                continue

                        terminal_options = []
                        for direction in cycle_directions:
                            options = []
                            for _, _, active_direction in active_edges(
                                shores[direction]
                            ):
                                support = frozenset(
                                    (direction, active_direction)
                                )
                                if not V.shatters_pair(words, support):
                                    options.append(support)
                            terminal_options.append(options)
                        if (
                            len(terminal_options) != 2
                            or not terminal_options[0]
                            or not terminal_options[1]
                        ):
                            continue

                        collar_supports = set()
                        for direction in cycle_directions:
                            for shore_word in shores[direction]:
                                difference = core_active_word ^ shore_word
                                if difference.bit_count() != 1:
                                    continue
                                support = frozenset(
                                    (direction, difference.bit_length() - 1)
                                )
                                if not V.shatters_pair(words, support):
                                    collar_supports.add(support)

                        for first_support in terminal_options[0]:
                            for second_support in terminal_options[1]:
                                if not first_support.isdisjoint(second_support):
                                    continue
                                marker = {
                                    "type": profile_name,
                                    "core": core,
                                    "passive": passive,
                                    "supports": frozenset(
                                        (first_support, second_support)
                                    ),
                                    "collars": frozenset(collar_supports),
                                    "shore_words": frozenset(
                                        active_word
                                        | (core_word & passive_mask)
                                        for direction in passive
                                        for active_word in shores[direction]
                                    ),
                                }
                                key = (
                                    profile_name,
                                    core,
                                    passive,
                                    marker["supports"],
                                    marker["collars"],
                                )
                                markers[key] = marker
    return list(markers.values())


def forced_owners(marker, graph, words, edge_class):
    core_word = words[marker["core"]]
    owners = set()
    for vertex in graph:
        if not V.is_corner(vertex, graph, words, edge_class):
            continue
        differences = [
            direction
            for direction in marker["passive"]
            if ((words[vertex] ^ core_word) >> direction) & 1
        ]
        if len(differences) >= 2:
            return None
        if len(differences) == 1:
            owners.add(words[vertex] ^ (1 << differences[0]))
    return owners


def comparison_key(family, first, second):
    def marker_key(marker):
        return (
            marker["type"],
            marker["core"],
            tuple(sorted(marker["passive"])),
            tuple(
                sorted(tuple(sorted(pair)) for pair in marker["supports"])
            ),
            tuple(
                sorted(
                    tuple(sorted(pair))
                    for pair in marker.get("collars", frozenset())
                )
            ),
        )

    return (
        family,
        marker_key(first),
        marker_key(second),
    )


def classify_comparison(first, second, graph, words, edge_class):
    all_terminal_pairs = set(first["supports"]) | set(second["supports"])
    if len(all_terminal_pairs) < 4:
        return "repeated_terminal"
    if matching_number_at_least_three(all_terminal_pairs):
        return "terminal_matching"
    profile = V.incidence_type(first, second)
    if profile is None:
        raise AssertionError("unclassified four-pair comparison")

    first_owners = forced_owners(first, graph, words, edge_class)
    second_owners = forced_owners(second, graph, words, edge_class)
    if first_owners is None or second_owners is None:
        return f"{profile}:distant_corner"

    embedded_words = set(words.values())
    if (
        first_owners & second_owners
        or (first_owners | second_owners) & embedded_words
    ):
        return f"{profile}:forced_owner_conflict"

    common_collars = first["collars"] & second["collars"]
    if common_collars:
        return f"{profile}:repeated_collar"

    acquired_pairs = (
        all_terminal_pairs | set(first["collars"]) | set(second["collars"])
    )
    if matching_number_at_least_three(acquired_pairs):
        return f"{profile}:collar_augmented_matching"
    return f"{profile}:unresolved"


def run_audit():
    base = M.build_tree(1, (0, 0))
    seen_families = set()
    seen_comparisons = set()
    outcome_counts = Counter()
    marker_counts = Counter()
    records = []

    for first_neighbours in V.subsets(base):
        first_extension = M.add_vertex(base, first_neighbours)
        for second_neighbours in V.subsets(first_extension):
            graph = M.add_vertex(first_extension, second_neighbours)
            if not M.connected(graph, set(graph)):
                continue
            embedding = M.partial_cube_embedding(graph)
            if embedding is None:
                continue
            words, edge_class = embedding
            dimension = 1 + max(edge_class.values(), default=-1)
            if not M.is_ample_embedding(words, dimension):
                continue
            family = V.canonical_family(words)
            if family in seen_families:
                continue
            seen_families.add(family)

            pure_markers = V.pure_decompositions(graph, edge_class)
            for marker in pure_markers:
                marker["type"] = "pure-eta2"
                marker["collars"] = acquired_collar_supports(marker, words)

            opposite_markers = []
            opposite_markers.extend(
                tree_arm_markers(
                    graph, words, edge_class, 0, 3, "sharp"
                )
            )
            opposite_markers.extend(
                tree_arm_markers(
                    graph, words, edge_class, 1, 2, "eta1"
                )
            )
            opposite_markers.extend(
                exact_two_unit_markers(words, edge_class)
            )
            for marker in opposite_markers:
                marker_counts[marker["type"]] += 1
                if "collars" not in marker:
                    marker["collars"] = acquired_collar_supports(
                        marker, words
                    )

            for pure in pure_markers:
                for opposite in opposite_markers:
                    key = comparison_key(family, pure, opposite)
                    if key in seen_comparisons:
                        continue
                    seen_comparisons.add(key)
                    outcome = classify_comparison(
                        pure, opposite, graph, words, edge_class
                    )
                    outcome_counts[(opposite["type"], outcome)] += 1
                    records.append(
                        {
                            "family": family,
                            "opposite_type": opposite["type"],
                            "outcome": outcome,
                        }
                    )

    unresolved = sum(
        count
        for (profile, outcome), count in outcome_counts.items()
        if outcome.endswith(":unresolved")
    )
    if unresolved:
        raise AssertionError(f"{unresolved} mixed comparisons remain")

    records.sort(
        key=lambda record: (
            record["family"],
            record["opposite_type"],
            record["outcome"],
        )
    )
    digest = hashlib.sha256(
        json.dumps(records, separators=(",", ":")).encode()
    ).hexdigest()
    return {
        "scope": "all comparisons with one pure eta=2 two-apex wing",
        "distinct_ample_two_apex_overlap_embeddings": len(seen_families),
        "opposite_marker_counts": dict(sorted(marker_counts.items())),
        "outcome_counts": {
            f"{profile} | {outcome}": count
            for (profile, outcome), count in sorted(outcome_counts.items())
        },
        "unresolved_comparisons": unresolved,
        "mathematical_digest": digest,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "damp_pure_eta_two_mixed_profiles.json",
    )
    args = parser.parse_args()
    result = run_audit()
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
