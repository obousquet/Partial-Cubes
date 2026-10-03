#!/usr/bin/env python3
"""Verify the rank-two carrier-tree profile of Hall's maximum class.

For every coordinate pair ``S``, the ``S``-cube carrier of a maximum
VC-dimension-three class is a maximum VC-one tree.  This verifier extracts
all 66 carrier trees of Hall's fixed 299-concept class and checks that each
is a path.  It also checks that every maximal 3-cube is an internal edge in
at least one of its three carrier paths, which supplies two opposite square
facets covering that cube.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import sys
from collections import Counter, defaultdict, deque
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from explore_smallest_additional_peripheral_obstruction import (  # noqa: E402
    HALL_CONCEPTS,
    shattered_coordinate_sets_mask,
)
from forbidden_pc_minor_campaign.verify_signed_clutter_parametrization import (  # noqa: E402
    signing,
)


DIMENSION = 12
OUTPUT = ROOT / "forbidden_pc_minor_campaign" / "damp_hall_carrier_trees.json"


def hall_maximal_cube_anchors() -> dict[int, int]:
    """Map each free triple support to its fixed ambient anchor word."""

    family = sum(1 << concept for concept in HALL_CONCEPTS)
    full_family = (1 << (1 << DIMENSION)) - 1
    complement = full_family ^ family
    complement_shadow = shattered_coordinate_sets_mask(complement, DIMENSION)
    fixed_signs = signing(complement, complement_shadow, DIMENSION)
    all_coordinates = (1 << DIMENSION) - 1
    anchors = {
        all_coordinates ^ fixed_support: fixed_word
        for fixed_support, fixed_word in fixed_signs.items()
    }
    expected_supports = {
        sum(1 << coordinate for coordinate in triple)
        for triple in itertools.combinations(range(DIMENSION), 3)
    }
    if set(anchors) != expected_supports:
        raise AssertionError("Hall's maximal-cube supports are not all triples")
    return anchors


def carrier_tree(pair: int, anchors: dict[int, int]) -> dict[str, object]:
    pair_coordinates = tuple(
        coordinate for coordinate in range(DIMENSION) if pair & (1 << coordinate)
    )
    outside = tuple(
        coordinate for coordinate in range(DIMENSION) if not pair & (1 << coordinate)
    )
    adjacency: dict[int, set[int]] = defaultdict(set)
    edge_rows = []
    internal_labels = []

    for coordinate in outside:
        triple = pair | (1 << coordinate)
        base_anchor = anchors[triple]
        left = base_anchor
        right = base_anchor | (1 << coordinate)
        adjacency[left].add(right)
        adjacency[right].add(left)
        edge_rows.append((coordinate, left, right))

    if len(edge_rows) != DIMENSION - 2:
        raise AssertionError("wrong carrier edge count")
    if len(adjacency) != DIMENSION - 1:
        raise AssertionError("wrong carrier vertex count")
    root = min(adjacency)
    seen = {root}
    queue = deque([root])
    while queue:
        vertex = queue.popleft()
        for neighbor in adjacency[vertex]:
            if neighbor not in seen:
                seen.add(neighbor)
                queue.append(neighbor)
    if len(seen) != len(adjacency):
        raise AssertionError("carrier is disconnected")
    degree_sequence = tuple(sorted(len(adjacency[vertex]) for vertex in adjacency))
    expected_degrees = (1, 1) + (2,) * (DIMENSION - 3)
    if degree_sequence != expected_degrees:
        raise AssertionError("Hall carrier is not a path")

    for coordinate, left, right in edge_rows:
        if len(adjacency[left]) > 1 and len(adjacency[right]) > 1:
            internal_labels.append(coordinate)

    return {
        "pair": list(pair_coordinates),
        "vertex_count": len(adjacency),
        "edge_count": len(edge_rows),
        "degree_sequence": list(degree_sequence),
        "internal_edge_labels": sorted(internal_labels),
        "internal_edge_count": len(internal_labels),
        "edges": [
            {"label": coordinate, "endpoints": [left, right]}
            for coordinate, left, right in sorted(edge_rows)
        ],
    }


def main() -> None:
    anchors = hall_maximal_cube_anchors()
    rows = []
    internal_ridge_count: Counter[int] = Counter()
    for pair_coordinates in itertools.combinations(range(DIMENSION), 2):
        pair = sum(1 << coordinate for coordinate in pair_coordinates)
        row = carrier_tree(pair, anchors)
        rows.append(row)
        for coordinate in row["internal_edge_labels"]:
            internal_ridge_count[pair | (1 << coordinate)] += 1

    expected_triples = set(anchors)
    if set(internal_ridge_count) != expected_triples:
        raise AssertionError("some Hall 3-cube lacks an internal carrier edge")
    histogram = Counter(internal_ridge_count.values())
    if histogram != Counter({1: 36, 2: 60, 3: 124}):
        raise AssertionError("unexpected Hall internal-ridge histogram")
    if any(row["internal_edge_count"] != DIMENSION - 4 for row in rows):
        raise AssertionError("a Hall carrier path has the wrong internal-edge count")

    mathematics = {
        "schema": "damp-hall-carrier-trees-v1",
        "source_class": {
            "dimension": DIMENSION,
            "order": len(HALL_CONCEPTS),
            "maximum_vc_dimension": 3,
            "maximal_cube_count": len(anchors),
        },
        "carrier_profile": {
            "pair_count": len(rows),
            "all_carriers_are_paths": True,
            "vertices_per_path": DIMENSION - 1,
            "edges_per_path": DIMENSION - 2,
            "internal_edges_per_path": DIMENSION - 4,
            "internal_ridge_count_per_triple_histogram": {
                str(key): value for key, value in sorted(histogram.items())
            },
            "every_maximal_cube_has_opposite_facet_cover": True,
        },
        "carrier_trees": rows,
    }
    digest = hashlib.sha256(
        json.dumps(mathematics, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    OUTPUT.write_text(
        json.dumps(
            {**mathematics, "mathematical_profile_sha256": digest},
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(f"carrier_paths={len(rows)}")
    print(f"internal_ridge_histogram={dict(sorted(histogram.items()))}")
    print(f"mathematical_profile_sha256={digest}")
    print(f"wrote {OUTPUT}")


if __name__ == "__main__":
    main()
