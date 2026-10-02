#!/usr/bin/env python3
"""Verify the connected-incidence exclusion in the pure eta=2 normal form.

The mathematical input is the two-apex normal form from
Corollary damp-eta-two-two-apex-normal-form: deleting two distinguished
vertices from the nine-vertex ample overlap leaves a P7 whose centre is the
one-vertex active slice and whose two length-three sides are the counted
arms.  Every such labelled overlap is obtained by adjoining two vertices to
one fixed P7 and choosing their neighbour sets.

This verifier enumerates that constant-size normal form, not cube families
by cardinality.  It rediscovers every pure decomposition of every ample
partial-cube completion, compares their two terminal-support systems, and
checks the passive-radius-one shielding obstruction.  It proves that the
connected incidence graphs C4 and P5 always expose, for one decomposition,
an uncounted corner at passive distance two.  In the disconnected 2P3 case,
cornerlessness forces exclusive shield owners.  These collide between the
two wings except in one P9 role-reversal form, where the two wings acquire
the same collar-square support.  Thus the entire pure/pure alternating
branch is excluded.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parent
HELPER_PATH = ROOT / "explore_damp_two_arm_tree_decompositions.py"


def load_helper():
    spec = importlib.util.spec_from_file_location("two_arm_helper", HELPER_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {HELPER_PATH}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


M = load_helper()


def subsets(vertices):
    vertices = tuple(vertices)
    for size in range(len(vertices) + 1):
        yield from map(set, combinations(vertices, size))


def edge_on_path(vertices, index):
    return M.canonical_edge(vertices[index], vertices[index + 1])


def pure_decompositions(graph, edge_class):
    """Return every distinguished pure eta=2 decomposition of graph."""

    output = []
    vertex_set = set(graph)
    for extras_tuple in combinations(vertex_set, 2):
        extras = frozenset(extras_tuple)
        tree = {
            vertex: neighbours - extras
            for vertex, neighbours in graph.items()
            if vertex not in extras
        }
        if len(M.graph_edges(tree)) != 6 or not M.connected(tree, set(tree)):
            continue

        for decomposition in M.decompositions(tree):
            if len(decomposition["core"]) != 1:
                continue
            arms = []
            valid = True
            for arm in decomposition["arms"]:
                vertices = arm["vertices"]
                path_edges = tuple(
                    edge_on_path(vertices, index) for index in range(3)
                )
                directions = tuple(edge_class[edge] for edge in path_edges)
                passive, _, terminal = directions
                other_tree_edges = set(M.graph_edges(tree)) - {path_edges[0]}
                if any(
                    edge_class[edge] == passive for edge in other_tree_edges
                ):
                    valid = False
                    break
                arms.append(
                    {
                        "vertices": tuple(vertices),
                        "directions": directions,
                        "support": frozenset((passive, terminal)),
                    }
                )
            if not valid:
                continue
            passive_directions = frozenset(
                arm["directions"][0] for arm in arms
            )
            supports = frozenset(arm["support"] for arm in arms)
            if len(passive_directions) != 2 or len(supports) != 2:
                continue
            first_support, second_support = tuple(supports)
            if (
                len(first_support) != 2
                or len(second_support) != 2
                or not first_support.isdisjoint(second_support)
            ):
                continue
            output.append(
                {
                    "extras": extras,
                    "core": next(iter(decomposition["core"])),
                    "passive": passive_directions,
                    "supports": supports,
                    "arms": tuple(arms),
                }
            )
    return output


def matching_number_at_least_three(pairs):
    return any(
        first.isdisjoint(second)
        and first.isdisjoint(third)
        and second.isdisjoint(third)
        for first, second, third in combinations(pairs, 3)
    )


def incidence_type(first, second):
    pairs = set(first["supports"]) | set(second["supports"])
    if len(pairs) != 4 or matching_number_at_least_three(pairs):
        return None

    adjacency = {direction: set() for pair in pairs for direction in pair}
    for x, y in pairs:
        adjacency[x].add(y)
        adjacency[y].add(x)

    components = []
    unseen = set(adjacency)
    while unseen:
        root = min(unseen)
        reached = {root}
        frontier = [root]
        unseen.remove(root)
        edge_twice = 0
        while frontier:
            vertex = frontier.pop()
            edge_twice += len(adjacency[vertex])
            for neighbour in adjacency[vertex]:
                if neighbour in unseen:
                    unseen.remove(neighbour)
                    reached.add(neighbour)
                    frontier.append(neighbour)
        components.append((len(reached), edge_twice // 2))
    components.sort()
    if components == [(4, 4)]:
        return "C4"
    if components == [(5, 4)]:
        return "P5"
    if components == [(3, 2), (3, 2)]:
        return "2P3"
    raise AssertionError(f"unexpected four-edge support graph {components}")


def square_supports_at(vertex, graph, words, edge_class):
    inverse_words = {word: owner for owner, word in words.items()}
    directions = {
        edge_class[M.canonical_edge(vertex, neighbour)]
        for neighbour in graph[vertex]
    }
    squares = set()
    for first, second in combinations(directions, 2):
        opposite_word = words[vertex] ^ (1 << first) ^ (1 << second)
        if opposite_word in inverse_words:
            squares.add(frozenset((first, second)))
    return directions, squares


def is_corner(vertex, graph, words, edge_class):
    """Test cornerhood; all cubes in this normal form have rank at most two."""

    directions, squares = square_supports_at(
        vertex, graph, words, edge_class
    )
    maximal_supports = set(squares)
    maximal_supports.update(
        frozenset((direction,))
        for direction in directions
        if not any(direction in square for square in squares)
    )
    return len(maximal_supports) == 1


def passive_distance(vertex, decomposition, words):
    core_word = words[decomposition["core"]]
    vertex_word = words[vertex]
    return sum(
        ((core_word >> direction) & 1)
        != ((vertex_word >> direction) & 1)
        for direction in decomposition["passive"]
    )


def has_distant_corner(decomposition, graph, words, edge_class):
    return any(
        passive_distance(vertex, decomposition, words) == 2
        and is_corner(vertex, graph, words, edge_class)
        for vertex in decomposition["extras"]
    )


def forced_passive_owners(decomposition, graph, words, edge_class):
    """Return the exclusive concepts forced by distance-one overlap corners."""

    core_word = words[decomposition["core"]]
    owners = set()
    for vertex in graph:
        if not is_corner(vertex, graph, words, edge_class):
            continue
        differing_passive = [
            direction
            for direction in decomposition["passive"]
            if ((words[vertex] ^ core_word) >> direction) & 1
        ]
        if len(differing_passive) >= 2:
            return None
        if len(differing_passive) == 1:
            owners.add(words[vertex] ^ (1 << differing_passive[0]))
    return owners


def role_signature(first, second):
    roles = []
    for first_arm in first["arms"]:
        for second_arm in second["arms"]:
            common = first_arm["support"] & second_arm["support"]
            if not common:
                continue
            direction = next(iter(common))
            first_role = (
                "A" if first_arm["directions"][0] == direction else "T"
            )
            second_role = (
                "A" if second_arm["directions"][0] == direction else "T"
            )
            roles.append(first_role + second_role)
    return tuple(sorted(roles))


def collar_supports(decomposition):
    return {
        frozenset((arm["directions"][0], arm["directions"][1]))
        for arm in decomposition["arms"]
    }


def shatters_pair(words, support):
    mask = sum(1 << direction for direction in support)
    return len({word & mask for word in words.values()}) == 4


def is_nine_vertex_path(graph):
    return (
        len(graph) == 9
        and len(M.graph_edges(graph)) == 8
        and sorted(map(len, graph.values())) == [1, 1] + [2] * 7
    )


def canonical_family(words):
    """A deterministic embedded-family identifier for the audit digest."""

    return tuple(sorted(words.values()))


def run_audit():
    base = M.build_tree(1, (0, 0))
    candidate_count = 0
    ample_completion_count = 0
    comparison_counts = {"C4": 0, "P5": 0, "2P3": 0}
    survivor_counts = {"C4": 0, "P5": 0, "2P3": 0}
    disconnected_resolution_counts = {
        "forced_owner_conflict": 0,
        "repeated_collar": 0,
        "unresolved": 0,
    }
    unique_disconnected_resolutions = {
        "forced_owner_conflict": set(),
        "repeated_collar": set(),
        "unresolved": set(),
    }
    records = []

    for first_neighbours in subsets(base):
        first_extension = M.add_vertex(base, first_neighbours)
        for second_neighbours in subsets(first_extension):
            graph = M.add_vertex(first_extension, second_neighbours)
            candidate_count += 1
            if not M.connected(graph, set(graph)):
                continue
            embedding = M.partial_cube_embedding(graph)
            if embedding is None:
                continue
            words, edge_class = embedding
            dimension = 1 + max(edge_class.values(), default=-1)
            if not M.is_ample_embedding(words, dimension):
                continue
            ample_completion_count += 1

            decompositions = pure_decompositions(graph, edge_class)
            for first, second in combinations(decompositions, 2):
                profile = incidence_type(first, second)
                if profile is None:
                    continue
                comparison_counts[profile] += 1
                shielding_obstruction = (
                    has_distant_corner(first, graph, words, edge_class)
                    or has_distant_corner(second, graph, words, edge_class)
                )
                if not shielding_obstruction:
                    survivor_counts[profile] += 1
                disconnected_resolution = None
                if profile == "2P3":
                    first_owners = forced_passive_owners(
                        first, graph, words, edge_class
                    )
                    second_owners = forced_passive_owners(
                        second, graph, words, edge_class
                    )
                    if first_owners is None or second_owners is None:
                        raise AssertionError(
                            "a 2P3 profile unexpectedly has a distant corner"
                        )
                    embedded_words = set(words.values())
                    owner_conflict = bool(
                        first_owners & second_owners
                        or (first_owners | second_owners) & embedded_words
                    )
                    if owner_conflict:
                        disconnected_resolution = "forced_owner_conflict"
                    else:
                        common_collars = (
                            collar_supports(first) & collar_supports(second)
                        )
                        if (
                            is_nine_vertex_path(graph)
                            and role_signature(first, second) == ("AT", "TA")
                            and len(common_collars) == 1
                            and not shatters_pair(
                                words, next(iter(common_collars))
                            )
                        ):
                            disconnected_resolution = "repeated_collar"
                        else:
                            disconnected_resolution = "unresolved"
                    disconnected_resolution_counts[
                        disconnected_resolution
                    ] += 1
                    comparison_key = (
                        canonical_family(words),
                        tuple(
                            sorted(
                                tuple(sorted(pair))
                                for pair in first["supports"]
                            )
                        ),
                        tuple(
                            sorted(
                                tuple(sorted(pair))
                                for pair in second["supports"]
                            )
                        ),
                    )
                    unique_disconnected_resolutions[
                        disconnected_resolution
                    ].add(comparison_key)
                records.append(
                    {
                        "family": canonical_family(words),
                        "profile": profile,
                        "shielding_obstruction": shielding_obstruction,
                        "disconnected_resolution": disconnected_resolution,
                    }
                )

    records.sort(
        key=lambda record: (
            record["family"],
            record["profile"],
            record["shielding_obstruction"],
            str(record["disconnected_resolution"]),
        )
    )
    mathematical_digest = hashlib.sha256(
        json.dumps(records, separators=(",", ":")).encode()
    ).hexdigest()
    result = {
        "scope": "nine-vertex pure-eta-two two-apex normal form",
        "labelled_two_vertex_extensions": candidate_count,
        "ample_partial_cube_completions_with_multiplicity": (
            ample_completion_count
        ),
        "comparison_counts_with_multiplicity": comparison_counts,
        "passive_radius_survivor_counts_with_multiplicity": survivor_counts,
        "disconnected_resolution_counts_with_multiplicity": (
            disconnected_resolution_counts
        ),
        "unique_disconnected_resolution_counts": {
            key: len(value)
            for key, value in unique_disconnected_resolutions.items()
        },
        "mathematical_digest": mathematical_digest,
    }

    if survivor_counts["C4"] or survivor_counts["P5"]:
        raise AssertionError(
            "a connected terminal-incidence comparison escaped shielding"
        )
    if disconnected_resolution_counts["unresolved"]:
        raise AssertionError("a disconnected incidence comparison survived")
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "damp_pure_eta_two_connected_incidence.json",
    )
    args = parser.parse_args()
    result = run_audit()
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
