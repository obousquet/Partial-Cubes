#!/usr/bin/env python3
"""Test global shielding compatibility of the three-pair outcomes.

The complete two-unit comparison audit classifies a comparison as soon as
its terminal and acquired-collar supports contain a matching of order three.
This diagnostic deliberately continues past that stopping rule.  It applies
the passive-radius-one and forced-owner tests to every such comparison, in
order to determine whether the apparent matching cases can occur in a
cornerless two-layer family at all.

This is a normal-form diagnostic, not an order-by-order family census.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from itertools import combinations
from itertools import permutations
from pathlib import Path

import verify_damp_two_unit_all_opposite_profiles as A


ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = ROOT / "damp_three_pair_global_compatibility.json"
_CANONICAL_ACTIVE_CACHE = {}


def pair_key(pair):
    return tuple(sorted(pair))


def marker_record(marker, words):
    return {
        "type": marker["type"],
        "core_word": words[marker["core"]],
        "passive": sorted(marker["passive"]),
        "supports": sorted(pair_key(pair) for pair in marker["supports"]),
        "collars": sorted(
            pair_key(pair) for pair in marker.get("collars", frozenset())
        ),
    }


def forced_exclusive_words(marker, words):
    """Recover the shore words in the common passive layer of the wing."""

    if "shore_words" in marker:
        return set(marker["shore_words"])
    embedded = set(words.values())
    passive_mask = sum(1 << direction for direction in marker["passive"])
    base_bits = words[marker["core"]] & passive_mask
    output = set()
    for word in embedded:
        normalized = (word & ~passive_mask) | base_bits
        if normalized not in embedded:
            output.add(normalized)
    return output


def varying_coordinates(family):
    """Coordinates that are nonconstant on a nonempty word family."""

    family = set(family)
    if not family:
        return set()
    first = next(iter(family))
    variation_mask = 0
    for word in family:
        variation_mask |= first ^ word
    return {
        direction
        for direction in range(variation_mask.bit_length())
        if (variation_mask >> direction) & 1
    }


def canonical_active_family(family):
    """Canonicalize a word family on its nonconstant coordinates."""

    cache_key = tuple(sorted(family))
    if cache_key in _CANONICAL_ACTIVE_CACHE:
        return _CANONICAL_ACTIVE_CACHE[cache_key]
    family = set(cache_key)
    coordinates = sorted(varying_coordinates(family))
    compressed = {
        sum(((word >> coordinate) & 1) << index
            for index, coordinate in enumerate(coordinates))
        for word in family
    }
    best = None
    for origin in compressed:
        translated = [word ^ origin for word in compressed]
        for order in permutations(range(len(coordinates))):
            candidate = tuple(sorted(
                sum(((word >> old) & 1) << new
                    for new, old in enumerate(order))
                for word in translated
            ))
            if best is None or candidate < best:
                best = candidate
    _CANONICAL_ACTIVE_CACHE[cache_key] = best
    return best


def matching_witness(pairs):
    for selected in combinations(sorted(pairs, key=pair_key), 3):
        if all(
            first.isdisjoint(second)
            for first, second in combinations(selected, 2)
        ):
            return tuple(pair_key(pair) for pair in selected)
    return None


def continued_outcome(first, second, graph, words, edge_class):
    first_owners = A.R.forced_owners(first, graph, words, edge_class)
    second_owners = A.R.forced_owners(second, graph, words, edge_class)
    if first_owners is None or second_owners is None:
        return "distant_corner", first_owners, second_owners

    embedded_words = set(words.values())
    if (
        first_owners & second_owners
        or (first_owners | second_owners) & embedded_words
    ):
        return "forced_owner_conflict", first_owners, second_owners

    return "shielding_compatible", first_owners, second_owners


def run_diagnostic():
    universe = A.exact_profile_universe()
    families = universe["families"]
    seen = set()
    counts = Counter()
    viable_records = []

    def inspect(comparison_class, family, first, second, graph, words, edge_class,
                symmetric=False):
        first_signature = A.marker_signature(first, words)
        second_signature = A.marker_signature(second, words)
        signatures = (first_signature, second_signature)
        markers = (first, second)
        if symmetric and second_signature < first_signature:
            signatures = (second_signature, first_signature)
            markers = (second, first)
        key = (family, comparison_class, signatures)
        if key in seen:
            return
        seen.add(key)

        original = A.R.classify_comparison(first, second, graph, words, edge_class)
        if original not in {"terminal_matching", "C4:collar_augmented_matching",
                            "P5:collar_augmented_matching",
                            "2P3:collar_augmented_matching"}:
            return

        terminal_pairs = set(first["supports"]) | set(second["supports"])
        all_pairs = (
            terminal_pairs
            | set(first.get("collars", frozenset()))
            | set(second.get("collars", frozenset()))
        )
        witness = matching_witness(all_pairs)
        if witness is None:
            raise AssertionError((comparison_class, original))

        outcome, first_owners, second_owners = continued_outcome(
            first, second, graph, words, edge_class
        )
        counts[(comparison_class, original, outcome)] += 1
        if outcome == "shielding_compatible":
            first_shores = forced_exclusive_words(first, words)
            second_shores = forced_exclusive_words(second, words)
            first_required = first_shores | first_owners
            second_required = second_shores | second_owners
            first_active = varying_coordinates(first_required)
            second_active = varying_coordinates(second_required)
            first_required_type = canonical_active_family(first_required)
            second_required_type = canonical_active_family(second_required)
            viable_records.append({
                "family": list(family),
                "order": len(words),
                "dimension": 1 + max(edge_class.values(), default=-1),
                "comparison_class": comparison_class,
                "original_outcome": original,
                "matching": witness,
                "first": marker_record(markers[0], words),
                "second": marker_record(markers[1], words),
                "first_forced_owners": sorted(first_owners),
                "second_forced_owners": sorted(second_owners),
                "first_forced_shore_words": sorted(first_shores),
                "second_forced_shore_words": sorted(second_shores),
                "first_required_exclusive_words": sorted(first_required),
                "second_required_exclusive_words": sorted(second_required),
                "first_required_exclusive_count": len(first_required),
                "second_required_exclusive_count": len(second_required),
                "first_forced_active_coordinates": sorted(first_active),
                "second_forced_active_coordinates": sorted(second_active),
                "first_forced_active_coordinate_count": len(first_active),
                "second_forced_active_coordinate_count": len(second_active),
                "first_required_active_type": first_required_type,
                "second_required_active_type": second_required_type,
            })

    for family, (graph, words, edge_class) in families.items():
        two_unit = A.R.exact_two_unit_markers(words, edge_class)
        sharp = A.R.tree_arm_markers(graph, words, edge_class, 0, 3, "sharp")
        eta_one = A.R.tree_arm_markers(graph, words, edge_class, 1, 2, "eta1")
        for marker in sharp + eta_one:
            marker["collars"] = A.R.acquired_collar_supports(marker, words)
        extra, _, _ = A.extra_owner_markers(graph, words, edge_class)

        for first, second in combinations(two_unit, 2):
            inspect("two-unit | two-unit", family, first, second,
                    graph, words, edge_class, symmetric=True)
        for first in two_unit:
            for second in sharp:
                inspect("two-unit | sharp", family, first, second,
                        graph, words, edge_class)
            for second in eta_one:
                inspect("two-unit | eta1-one-unit", family, first, second,
                        graph, words, edge_class)
            for second in extra:
                inspect("two-unit | pair-rich-extra-owner", family, first, second,
                        graph, words, edge_class)

    viable_records.sort(key=lambda record: json.dumps(record, sort_keys=True))
    canonical = json.dumps(viable_records, sort_keys=True, separators=(",", ":"))
    active_coordinate_counts = Counter(
        (
            record["first_forced_active_coordinate_count"],
            record["second_forced_active_coordinate_count"],
        )
        for record in viable_records
    )
    six_active = sum(
        count
        for (first_count, second_count), count in active_coordinate_counts.items()
        if max(first_count, second_count) >= 6
    )
    five_coordinate_records = [
        record for record in viable_records
        if record["first_forced_active_coordinate_count"] == 5
        and record["second_forced_active_coordinate_count"] == 5
    ]
    if six_active + len(five_coordinate_records) != len(viable_records):
        raise AssertionError("a matching survives below five active coordinates")

    five_side_states = {}
    for record in five_coordinate_records:
        for side in ("first", "second"):
            state = {
                "family": record["family"],
                "active": record[f"{side}_forced_active_coordinates"],
                "required": record[f"{side}_required_exclusive_words"],
                "required_type": record[f"{side}_required_active_type"],
            }
            key = json.dumps(state, sort_keys=True, separators=(",", ":"))
            five_side_states[key] = state
    return {
        "scope": "continued shielding audit of every three-pair comparison",
        "comparison_counts": {
            " | ".join(key): value for key, value in sorted(counts.items())
        },
        "total_matching_comparisons": sum(counts.values()),
        "shielding_compatible_comparisons": len(viable_records),
        "active_coordinate_resolution_counts": {
            f"{first_count},{second_count}": count
            for (first_count, second_count), count in sorted(
                active_coordinate_counts.items()
            )
        },
        "six_active_coordinate_comparisons": six_active,
        "five_coordinate_comparisons": len(five_coordinate_records),
        "five_coordinate_side_state_count": len(five_side_states),
        "five_coordinate_side_states": [
            five_side_states[key] for key in sorted(five_side_states)
        ],
        "mathematical_digest": hashlib.sha256(canonical.encode()).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    arguments = parser.parse_args()
    record = run_diagnostic()
    arguments.output.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "total_matching_comparisons": record["total_matching_comparisons"],
        "shielding_compatible_comparisons": record[
            "shielding_compatible_comparisons"
        ],
        "active_coordinate_resolution_counts": record[
            "active_coordinate_resolution_counts"
        ],
        "six_active_coordinate_comparisons": record[
            "six_active_coordinate_comparisons"
        ],
        "five_coordinate_comparisons": record[
            "five_coordinate_comparisons"
        ],
        "five_coordinate_side_state_count": record[
            "five_coordinate_side_state_count"
        ],
        "comparison_counts": record["comparison_counts"],
        "mathematical_digest": record["mathematical_digest"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
