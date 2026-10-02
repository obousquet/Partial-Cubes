#!/usr/bin/env python3
"""Complete the sharp five-coordinate residue of the three-pair audit.

After passive shielding and the five-active-coordinate bound, every surviving
three-pair comparison has an overlap of order nine and forces one of two
six-vertex subsets in each exclusive wing: a square plus a disjoint edge, or
two disjoint three-vertex paths.  This diagnostic enumerates all induced
two-connected completions of those required vertices through relative order
eleven and tests the full upper family for ampleness and terminality.  A
forced shore need not display every coordinate of its eventual residual
graph, so the module also enumerates all active-coordinate super-fibres of
dimension at most five, including both orientations of a coordinate invisible
on the lower family.

The input is the complete global-compatibility diagnostic.  No ambient-family
cardinality is enumerated.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from itertools import combinations
from pathlib import Path

import verify_damp_pure_eta_two_mixed_profiles as R
import verify_damp_low_excess_global_matching as L


ROOT = Path(__file__).resolve().parent
DEFAULT_INPUT = ROOT / "damp_three_pair_global_compatibility.json"
DEFAULT_OUTPUT = ROOT / "damp_five_coordinate_matching_completions.json"


def cube_graph(vertices, dimension):
    vertices = set(vertices)
    return {
        vertex: {
            vertex ^ (1 << coordinate)
            for coordinate in range(dimension)
            if (vertex ^ (1 << coordinate)) in vertices
        }
        for vertex in vertices
    }


def connected_after_deleting(graph, deleted=None):
    remaining = set(graph)
    if deleted is not None:
        remaining.discard(deleted)
    if not remaining:
        return True
    seen = {next(iter(remaining))}
    stack = list(seen)
    while stack:
        vertex = stack.pop()
        for neighbour in graph[vertex] & remaining - seen:
            seen.add(neighbour)
            stack.append(neighbour)
    return seen == remaining


def is_two_connected(vertices, dimension):
    graph = cube_graph(vertices, dimension)
    return (
        len(graph) >= 3
        and connected_after_deleting(graph)
        and all(connected_after_deleting(graph, vertex) for vertex in graph)
    )


def is_ample(family, dimension):
    return R.is_ample_word_family(set(family), range(dimension))


def terminality_profile(lower, exclusive, dimension):
    upper = lower | exclusive
    addable = {
        word for word in exclusive if is_ample(lower | {word}, dimension)
    }
    deletable = {
        word for word in exclusive if is_ample(upper - {word}, dimension)
    }
    return addable, deletable


def state_key(state):
    return (
        tuple(state["family"]),
        tuple(state["active"]),
        tuple(state["required"]),
    )


def required_type(record, side):
    return tuple(record[f"{side}_required_active_type"])


def sharp_side_states(source):
    if "five_coordinate_side_states" in source:
        return [
            {"source": "two-unit matching residue", **state}
            for state in source["five_coordinate_side_states"]
        ]
    states = {}
    for record in source["shielding_compatible_records"]:
        if not (
            record["first_forced_active_coordinate_count"] == 5
            and record["second_forced_active_coordinate_count"] == 5
        ):
            continue
        for side in ("first", "second"):
            state = {
                "source": "two-unit matching residue",
                "family": record["family"],
                "active": record[f"{side}_forced_active_coordinates"],
                "required": record[f"{side}_required_exclusive_words"],
                "required_type": list(required_type(record, side)),
            }
            states[state_key(state)] = state
    return list(states.values())


def low_excess_path_equality_states():
    """All positioned P8 side states from the low-excess verifier."""

    source = L.run()
    return [
        {
            "source": "low-excess P8 equality residue",
            **state,
        }
        for state in source["five_coordinate_side_states"]
    ]


def varying_coordinates(words, dimension):
    words = set(words)
    return {
        coordinate
        for coordinate in range(dimension)
        if len({(word >> coordinate) & 1 for word in words}) == 2
    }


def completions(state, relative_order, *, require_all_active=False):
    lower = set(state["family"])
    required = set(state["required"])
    active = tuple(state["active"])
    dimension = max(
        max(lower | required).bit_length(),
        max(active, default=-1) + 1,
    )
    active_mask = sum(1 << coordinate for coordinate in active)
    fixed_mask = ((1 << dimension) - 1) ^ active_mask
    fixed_words = {word & fixed_mask for word in required}
    if len(fixed_words) != 1:
        raise AssertionError((state, fixed_words))
    fixed_word = next(iter(fixed_words))
    fibre = {
        fixed_word
        | sum(((choice >> index) & 1) << coordinate
              for index, coordinate in enumerate(active))
        for choice in range(1 << len(active))
    }
    if not required <= fibre or required & lower:
        raise AssertionError(state)
    candidates = sorted(fibre - lower - required)
    extra_count = relative_order - len(required)

    counts = Counter()
    examples = []
    for extra in combinations(candidates, extra_count):
        exclusive = required | set(extra)
        counts["candidate_subsets"] += 1
        if (
            require_all_active
            and varying_coordinates(exclusive, dimension) != set(active)
        ):
            counts["inactive_allowed_coordinate"] += 1
            continue
        counts["exact_active"] += 1
        if not is_two_connected(exclusive, dimension):
            continue
        counts["two_connected"] += 1
        upper = lower | exclusive
        if not is_ample(upper, dimension):
            continue
        counts["ample_upper"] += 1
        addable, deletable = terminality_profile(lower, exclusive, dimension)
        if addable:
            counts["has_addable_exclusive"] += 1
        else:
            counts["addition_terminal"] += 1
        if deletable:
            counts["has_deletable_exclusive"] += 1
        else:
            counts["deletion_terminal"] += 1
        if not addable and not deletable:
            counts["terminal"] += 1
            if len(examples) < 5:
                examples.append(sorted(exclusive))
    return dict(counts), examples


def active_superfibre_states(state, max_active=5):
    """Enumerate every positioned active-coordinate super-fibre.

    The words forced by a normal form need not vary on every coordinate of
    the eventual residual graph.  Coordinates already visible in the lower
    family retain their positions.  Coordinates invisible on both the lower
    family and the forced words are symmetric; after complementing them, the
    lower family is fixed at zero, while the forced shore may be fixed at an
    arbitrary binary pattern.  We therefore use canonical fresh coordinates
    and enumerate that pattern exactly.
    """

    lower = set(state["family"])
    required = set(state["required"])
    forced_active = tuple(sorted(state["active"]))
    base_dimension = max(
        max(lower | required).bit_length(),
        max(forced_active, default=-1) + 1,
    )
    if len(forced_active) > max_active:
        return []

    existing_inactive = tuple(
        coordinate
        for coordinate in range(base_dimension)
        if coordinate not in forced_active
    )
    records = []
    capacity = max_active - len(forced_active)
    for existing_count in range(min(capacity, len(existing_inactive)) + 1):
        for existing in combinations(existing_inactive, existing_count):
            fresh_capacity = capacity - existing_count
            for fresh_count in range(fresh_capacity + 1):
                fresh = tuple(
                    range(base_dimension, base_dimension + fresh_count)
                )
                for pattern in range(1 << fresh_count):
                    fresh_mask = sum(
                        ((pattern >> index) & 1) << coordinate
                        for index, coordinate in enumerate(fresh)
                    )
                    records.append({
                        **state,
                        "active": sorted(forced_active + existing + fresh),
                        "required": sorted(word | fresh_mask for word in required),
                        "added_existing_active_coordinates": list(existing),
                        "fresh_active_coordinate_count": fresh_count,
                        "fresh_required_pattern": pattern,
                    })
    return records


def completion_superfibre_record(state, relative_order, max_active=5):
    """Complete a forced state in every possible residual super-fibre."""

    records = []
    totals = Counter()
    terminal_examples = []
    for expanded_state in active_superfibre_states(state, max_active):
        counts, examples = completions(
            expanded_state,
            relative_order,
            require_all_active=True,
        )
        totals.update(counts)
        records.append({
            "state": expanded_state,
            "completion_counts": counts,
            "terminal_examples": examples,
        })
        for example in examples:
            if len(terminal_examples) < 5:
                terminal_examples.append({
                    "state": expanded_state,
                    "exclusive": example,
                })
    return {
        "forced_state": state,
        "active_superfibre_count": len(records),
        "completion_counts": dict(totals),
        "terminal_examples": terminal_examples,
        "active_superfibre_records": records,
    }


def run(input_path):
    source = json.loads(input_path.read_text())
    states = sharp_side_states(source) + low_excess_path_equality_states()
    records = []
    totals = Counter()
    for state in states:
        for relative_order in (10, 11):
            counts, examples = completions(state, relative_order)
            for name, count in counts.items():
                totals[(state["source"], relative_order, name)] += count
            records.append({
                **state,
                "relative_order": relative_order,
                "completion_counts": counts,
                "terminal_examples": examples,
            })
    records.sort(key=lambda record: json.dumps(record, sort_keys=True))
    for record in records:
        counts = record["completion_counts"]
        ample = counts.get("ample_upper", 0)
        if counts.get("terminal", 0):
            raise AssertionError("a terminal five-coordinate completion survives")
        if (
            counts.get("has_addable_exclusive", 0) != ample
            or counts.get("has_deletable_exclusive", 0) != ample
        ):
            raise AssertionError(
                "an ample completion lacks one of the two terminality failures"
            )
    canonical = json.dumps(records, sort_keys=True, separators=(",", ":"))
    return {
        "scope": (
            "all five-active-coordinate side states surviving the "
            "three-pair shielding audit, through relative order eleven"
        ),
        "side_state_count": len(states),
        "side_state_counts_by_source": dict(Counter(
            state["source"] for state in states
        )),
        "required_active_types": sorted({
            tuple(state["required_type"]) for state in states
        }),
        "totals": {
            f"{source} | order_{order}_{name}": count
            for (source, order, name), count in sorted(totals.items())
        },
        "unresolved_terminal_completions": sum(
            record["completion_counts"].get("terminal", 0)
            for record in records
        ),
        "records": records,
        "mathematical_digest": hashlib.sha256(canonical.encode()).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    record = run(arguments.input)
    arguments.output.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "side_state_count": record["side_state_count"],
        "required_active_types": record["required_active_types"],
        "totals": record["totals"],
        "unresolved_terminal_completions": record[
            "unresolved_terminal_completions"
        ],
        "mathematical_digest": record["mathematical_digest"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
