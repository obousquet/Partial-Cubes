#!/usr/bin/env python3
"""Extract balanced Yang order-sphere circuit carriers from Hall seeds.

The order-sphere theorem proves that every nonpeelable family has a
collection I of latent charts whose failure-relation digraph is a disjoint
union of strongly connected components and satisfies

    |I| + number_of_components = |H| + 1.

This worker diagnostic searches for the particularly transparent case in
which the selected relation digraph is strongly connected on all concepts.
It then has exactly |H| charts.  The computation is not used to prove the
theorem or nonpeelability of either seed.
"""

from __future__ import annotations

import hashlib
import itertools
import json
import random
from collections import Counter, deque
from pathlib import Path

from explore_smallest_additional_peripheral_obstruction import HALL_CONCEPTS


HERE = Path(__file__).resolve().parent
MODEL_291 = HERE / "damp_hall_carrier_defect_min8.json"
OUTPUT = HERE / "damp_yang_balanced_circuits.json"


def digest(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def suppress(vertex: int, coordinate: int) -> int:
    lower = (1 << coordinate) - 1
    return (vertex & lower) | ((vertex >> (coordinate + 1)) << coordinate)


def submasks(mask: int) -> tuple[int, ...]:
    answer = []
    current = mask
    while True:
        answer.append(current)
        if current == 0:
            return tuple(answer)
        current = (current - 1) & mask


def reconstruct_291() -> frozenset[int]:
    model = json.loads(MODEL_291.read_text(encoding="utf-8"))["mathematics"][
        "verified_model"
    ]
    projection = {suppress(concept, 0) for concept in HALL_CONCEPTS}
    overlap = set(map(int, model["overlap_members"]))
    red_supports = set(map(int, model["red_supports"]))
    lower = set(overlap)
    upper = set(overlap)
    for coordinates in itertools.combinations(range(11), 3):
        support = sum(1 << coordinate for coordinate in coordinates)
        anchor_counts = Counter(vertex & ~support for vertex in projection)
        anchors = [
            anchor
            for anchor, count in anchor_counts.items()
            if count == 1 << len(coordinates)
        ]
        if len(anchors) != 1:
            raise AssertionError((support, anchors))
        cube = {anchors[0] | word for word in submasks(support)}
        (lower if support in red_supports else upper).update(cube)
    family = lower | {word | (1 << 11) for word in upper}
    if len(family) != 291:
        raise AssertionError(len(family))
    return frozenset(family)


def latent_blocks(
    family: frozenset[int], dimension: int
) -> tuple[dict[str, object], ...]:
    blocks = []
    for order in range(2, min(4, dimension) + 1):
        for coordinates in itertools.combinations(range(dimension), order):
            support = sum(1 << coordinate for coordinate in coordinates)
            flips = tuple(flip for flip in submasks(support) if flip)
            for tip in range(1 << dimension):
                if not all(tip ^ flip in family for flip in flips):
                    continue
                antipode = tip ^ support
                protectors = tuple(
                    tip ^ (support ^ (1 << coordinate))
                    for coordinate in coordinates
                )
                arcs = []
                if tip in family:
                    arcs.append((tip, antipode))
                arcs.extend((antipode, protector) for protector in protectors)
                blocks.append(
                    {
                        "tip": tip,
                        "support": support,
                        "support_order": order,
                        "tip_present": tip in family,
                        "arcs": tuple(arcs),
                    }
                )
    return tuple(blocks)


def adjacency(
    vertices: tuple[int, ...],
    blocks: tuple[dict[str, object], ...],
    selected: set[int] | None = None,
) -> tuple[dict[int, list[tuple[int, int]]], dict[int, list[tuple[int, int]]]]:
    forward = {vertex: [] for vertex in vertices}
    reverse = {vertex: [] for vertex in vertices}
    indices = range(len(blocks)) if selected is None else sorted(selected)
    for index in indices:
        for source, target in blocks[index]["arcs"]:
            forward[source].append((target, index))
            reverse[target].append((source, index))
    return forward, reverse


def reached(
    start: int, graph: dict[int, list[tuple[int, int]]]
) -> set[int]:
    seen = {start}
    queue = [start]
    while queue:
        source = queue.pop()
        for target, _ in graph[source]:
            if target not in seen:
                seen.add(target)
                queue.append(target)
    return seen


def strongly_connected(
    vertices: tuple[int, ...],
    blocks: tuple[dict[str, object], ...],
    selected: set[int],
) -> bool:
    forward, reverse = adjacency(vertices, blocks, selected)
    root = vertices[0]
    return len(reached(root, forward)) == len(vertices) and len(
        reached(root, reverse)
    ) == len(vertices)


def tree_blocks(
    root: int,
    graph: dict[int, list[tuple[int, int]]],
    rng: random.Random,
) -> set[int]:
    seen = {root}
    queue = deque([root])
    chosen = set()
    while queue:
        source = queue.popleft()
        edges = list(graph[source])
        rng.shuffle(edges)
        for target, block in edges:
            if target in seen:
                continue
            seen.add(target)
            chosen.add(block)
            queue.append(target)
    return chosen


def extract_carrier(
    family: frozenset[int], dimension: int, label: str
) -> dict[str, object]:
    vertices = tuple(sorted(family))
    blocks = latent_blocks(family, dimension)
    full_forward, full_reverse = adjacency(vertices, blocks)
    root = vertices[0]
    if len(reached(root, full_forward)) != len(vertices) or len(
        reached(root, full_reverse)
    ) != len(vertices):
        raise AssertionError((label, "full failure-relation graph is not strong"))

    rng = random.Random(20260829)
    best = set(range(len(blocks)))
    trial_sizes = []
    for _ in range(96):
        selected = tree_blocks(root, full_forward, rng) | tree_blocks(
            root, full_reverse, rng
        )
        candidates = list(selected)
        rng.shuffle(candidates)
        for block in candidates:
            trial = selected - {block}
            if strongly_connected(vertices, blocks, trial):
                selected = trial
        trial_sizes.append(len(selected))
        if len(selected) < len(best):
            best = selected
        if len(best) <= len(vertices):
            break

    if len(best) > len(vertices):
        raise AssertionError((label, "carrier search too large", len(best)))
    selected = set(best)
    for index in range(len(blocks)):
        if len(selected) == len(vertices):
            break
        selected.add(index)
    if len(selected) != len(vertices):
        raise AssertionError((label, len(selected), len(vertices)))
    if not strongly_connected(vertices, blocks, selected):
        raise AssertionError((label, "padded carrier lost strong connectivity"))

    chosen = [blocks[index] for index in sorted(selected)]
    rank_counts = Counter(
        (block["support_order"], block["tip_present"]) for block in chosen
    )
    mathematics = {
        "label": label,
        "dimension": dimension,
        "order": len(vertices),
        "latent_block_count": len(blocks),
        "full_relation_graph_strongly_connected": True,
        "initial_strong_carrier_size": min(trial_sizes),
        "selected_block_count": len(chosen),
        "strong_component_count": 1,
        "balance_left": len(chosen) + 1,
        "balance_right": len(vertices) + 1,
        "balance_verified": len(chosen) + 1 == len(vertices) + 1,
        "selected_rank_tip_counts": [
            {
                "support_order": order,
                "tip_present": present,
                "multiplicity": multiplicity,
            }
            for (order, present), multiplicity in sorted(rank_counts.items())
        ],
        "selected_blocks": [
            {
                "tip": block["tip"],
                "support": block["support"],
                "support_order": block["support_order"],
                "tip_present": block["tip_present"],
                "arcs": [list(arc) for arc in block["arcs"]],
            }
            for block in chosen
        ],
    }
    mathematics["mathematical_sha256"] = digest(mathematics)
    return mathematics


def main() -> None:
    rows = [
        extract_carrier(frozenset(HALL_CONCEPTS), 12, "Hall-299"),
        extract_carrier(reconstruct_291(), 12, "Hall-291"),
    ]
    mathematics = {
        "theorem": "balanced Yang order-sphere circuit carrier",
        "rows": rows,
    }
    mathematics["mathematical_sha256"] = digest(mathematics)
    payload = {
        "schema": "damp-yang-balanced-circuits-v1",
        "mathematics": mathematics,
        "metadata": {"generator": Path(__file__).name},
    }
    OUTPUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    for row in rows:
        print(
            row["label"],
            f"latent={row['latent_block_count']}",
            f"initial={row['initial_strong_carrier_size']}",
            f"selected={row['selected_block_count']}",
            f"digest={row['mathematical_sha256']}",
        )
    print(f"mathematical_sha256={mathematics['mathematical_sha256']}")
    print(f"wrote {OUTPUT}")


if __name__ == "__main__":
    main()
