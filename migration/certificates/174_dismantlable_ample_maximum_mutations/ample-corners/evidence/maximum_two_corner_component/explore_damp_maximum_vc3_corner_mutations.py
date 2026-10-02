#!/usr/bin/env python3
"""Search for a smaller cornerless maximum VC-three class by exact mutations.

The search starts from coordinate projections of Hall's 299-concept class.
Each projection is a maximum VC-three class of order 232 in ``Q_11``.  A
proposal removes one concept and adds one missing word.  It is accepted as a
legal proposal only when every trace on at most three coordinates remains
realized; because the order is the Sauer bound, this is equivalent to
remaining a maximum VC-three class.

The objective is the exact number of corners, computed from the unique full
three-cube on every triple support.  The walk is a bounded construction
search.  Finding zero corners is a finite upper-bound certificate; failure to
find one is not evidence of nonexistence.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import random
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from explore_smallest_additional_peripheral_obstruction import (
    HALL_CONCEPTS,
    project_vertex,
)


CAMPAIGN = ROOT / "forbidden_pc_minor_campaign"


def subsets(mask: int) -> tuple[int, ...]:
    answer = []
    subset = mask
    while True:
        answer.append(subset)
        if subset == 0:
            return tuple(answer)
        subset = (subset - 1) & mask


def maximum_order(dimension: int) -> int:
    return sum(math.comb(dimension, size) for size in range(4))


def supports_at_most_three(dimension: int) -> tuple[int, ...]:
    return tuple(
        sum(1 << coordinate for coordinate in coordinates)
        for size in range(4)
        for coordinates in itertools.combinations(range(dimension), size)
    )


def triple_supports(dimension: int) -> tuple[int, ...]:
    return tuple(
        sum(1 << coordinate for coordinate in coordinates)
        for coordinates in itertools.combinations(range(dimension), 3)
    )


def quadruple_supports(dimension: int) -> tuple[int, ...]:
    return tuple(
        sum(1 << coordinate for coordinate in coordinates)
        for coordinates in itertools.combinations(range(dimension), 4)
    )


def hall_projection(coordinate: int) -> frozenset[int]:
    return frozenset(project_vertex(concept, coordinate) for concept in HALL_CONCEPTS)


@dataclass(frozen=True)
class MaximumState:
    family: frozenset[int]
    cube_anchor: dict[int, int]
    incidence: dict[int, int]
    corners: frozenset[int]


def build_state(family: frozenset[int], dimension: int) -> MaximumState:
    if len(family) != maximum_order(dimension):
        raise AssertionError((len(family), maximum_order(dimension)))
    cube_anchor: dict[int, int] = {}
    incidence = Counter()
    for support in triple_supports(dimension):
        outside_counts = Counter(concept & ~support for concept in family)
        full = [anchor for anchor, count in outside_counts.items() if count == 8]
        if len(full) != 1:
            raise AssertionError((support, full))
        anchor = full[0]
        cube_anchor[support] = anchor
        incidence.update(anchor | word for word in subsets(support))
    if set(incidence) != set(family):
        raise AssertionError("some concept is in no full three-cube")
    corners = frozenset(concept for concept in family if incidence[concept] == 1)
    return MaximumState(family, cube_anchor, dict(incidence), corners)


def critical_masks(family: frozenset[int], dimension: int) -> dict[int, int]:
    """Coordinates on which an added word must agree with the removed word."""

    answer = {concept: 0 for concept in family}
    for support in supports_at_most_three(dimension):
        counts = Counter(concept & support for concept in family)
        for concept in family:
            if counts[concept & support] == 1:
                answer[concept] |= support
    return answer


def legal_swaps(
    family: frozenset[int], dimension: int
) -> tuple[tuple[int, int], ...]:
    """Return all one-out/one-in swaps preserving maximum VC dimension three."""

    universe = frozenset(range(1 << dimension))
    outside = universe - family
    masks = critical_masks(family, dimension)
    missing_trace_bits = {word: 0 for word in universe}
    nonunique_trace_bits = {concept: 0 for concept in family}
    for index, support in enumerate(quadruple_supports(dimension)):
        counts = Counter(concept & support for concept in family)
        possible = tuple(subsets(support))
        missing = [trace for trace in possible if counts[trace] == 0]
        if len(missing) != 1:
            raise AssertionError((support, missing))
        missing_trace = missing[0]
        bit = 1 << index
        for word in universe:
            if word & support == missing_trace:
                missing_trace_bits[word] |= bit
        for concept in family:
            if counts[concept & support] > 1:
                nonunique_trace_bits[concept] |= bit
    return tuple(
        (removed, added)
        for removed in sorted(family)
        for added in sorted(outside)
        if (removed ^ added) & masks[removed] == 0
        and not (
            missing_trace_bits[added] & nonunique_trace_bits[removed]
        )
    )


def swap_corner_count(
    state: MaximumState,
    removed: int,
    added: int,
    dimension: int,
) -> int:
    """Compute the exact corner count after a legal swap without rebuilding."""

    family_after = (state.family - {removed}) | {added}
    incidence_delta = Counter()
    for support, old_anchor in state.cube_anchor.items():
        if removed & ~support != old_anchor:
            continue
        new_anchor = added & ~support
        new_cube = tuple(new_anchor | word for word in subsets(support))
        if any(vertex not in family_after for vertex in new_cube):
            raise AssertionError((removed, added, support, new_anchor))
        old_cube = tuple(old_anchor | word for word in subsets(support))
        incidence_delta.subtract(old_cube)
        incidence_delta.update(new_cube)

    affected = set(incidence_delta) | {removed, added}
    old_corner_count = len(state.corners)
    old_affected_corners = sum(vertex in state.corners for vertex in affected)
    new_affected_corners = 0
    for vertex in affected:
        if vertex not in family_after:
            continue
        new_incidence = state.incidence.get(vertex, 0) + incidence_delta[vertex]
        if new_incidence == 1:
            new_affected_corners += 1
        elif new_incidence < 1:
            raise AssertionError((vertex, new_incidence))
    return old_corner_count - old_affected_corners + new_affected_corners


def walk(
    initial: frozenset[int],
    dimension: int,
    steps: int,
    candidate_limit: int,
    rng: random.Random,
) -> tuple[MaximumState, list[dict[str, int]]]:
    state = build_state(initial, dimension)
    best = state
    trace = [{"step": 0, "corner_count": len(state.corners)}]
    seen = {state.family}

    for step in range(1, steps + 1):
        swaps = list(legal_swaps(state.family, dimension))
        if not swaps:
            break
        rng.shuffle(swaps)
        if candidate_limit and len(swaps) > candidate_limit:
            swaps = swaps[:candidate_limit]
        scored = []
        for removed, added in swaps:
            candidate_family = (state.family - {removed}) | {added}
            if candidate_family in seen:
                continue
            score = swap_corner_count(state, removed, added, dimension)
            scored.append((score, rng.random(), removed, added))
        if not scored:
            seen.clear()
            continue
        scored.sort()

        current_score = len(state.corners)
        improving = [row for row in scored if row[0] < current_score]
        if improving:
            _, _, removed, added = improving[0]
        else:
            # Diversify among the best few exact mutations at a local minimum.
            pool = scored[: min(16, len(scored))]
            _, _, removed, added = rng.choice(pool)
        state = build_state((state.family - {removed}) | {added}, dimension)
        seen.add(state.family)
        if len(state.corners) < len(best.corners):
            best = state
            trace.append({"step": step, "corner_count": len(best.corners)})
        if not best.corners:
            break
    return best, trace


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dimension", type=int, default=11)
    parser.add_argument("--steps", type=int, default=250)
    parser.add_argument("--starts", type=int, default=12)
    parser.add_argument("--candidate-limit", type=int, default=4000)
    parser.add_argument("--seed", type=int, default=20260820)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if args.dimension != 11:
        raise ValueError("the current bounded search uses Hall's Q_11 projections")

    rng = random.Random(args.seed)
    rows = []
    global_best: MaximumState | None = None
    for start in range(min(args.starts, 12)):
        initial = hall_projection(start)
        initial_state = build_state(initial, args.dimension)
        initial_swaps = legal_swaps(initial, args.dimension)
        initial_neighbor_histogram = Counter(
            swap_corner_count(initial_state, removed, added, args.dimension)
            for removed, added in initial_swaps
        )
        best, trace = walk(
            initial,
            args.dimension,
            args.steps,
            args.candidate_limit,
            rng,
        )
        rows.append(
            {
                "projected_coordinate": start,
                "initial_corner_count": trace[0]["corner_count"],
                "initial_legal_mutation_count": len(initial_swaps),
                "initial_neighbor_corner_histogram": {
                    str(corner_count): multiplicity
                    for corner_count, multiplicity in sorted(
                        initial_neighbor_histogram.items()
                    )
                },
                "best_corner_count": len(best.corners),
                "improvement_trace": trace,
                "best_family_sha256": hashlib.sha256(
                    json.dumps(sorted(best.family), separators=(",", ":")).encode()
                ).hexdigest(),
                "best_family": sorted(best.family) if not best.corners else None,
            }
        )
        if global_best is None or len(best.corners) < len(global_best.corners):
            global_best = best
        print(
            f"start={start} initial={trace[0]['corner_count']} "
            f"best={len(best.corners)}"
        )
        if global_best is not None and not global_best.corners:
            break

    assert global_best is not None
    mathematics = {
        "schema": "damp-maximum-vc3-corner-mutation-search-v1",
        "scope": (
            "Bounded exact same-shadow mutation search from the twelve Hall "
            "coordinate projections; a zero-corner row is a construction, "
            "whereas a positive minimum is only a search diagnostic."
        ),
        "dimension": args.dimension,
        "order": maximum_order(args.dimension),
        "steps_per_start": args.steps,
        "candidate_limit": args.candidate_limit,
        "seed": args.seed,
        "starts_completed": len(rows),
        "best_corner_count": len(global_best.corners),
        "found_cornerless_class": not global_best.corners,
        "rows": rows,
    }
    canonical = json.dumps(mathematics, sort_keys=True, separators=(",", ":"))
    payload = {
        **mathematics,
        "mathematical_profile_sha256": hashlib.sha256(canonical.encode()).hexdigest(),
    }
    output = args.output or CAMPAIGN / "damp_maximum_vc3_corner_mutation_d11.json"
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"best_corner_count={len(global_best.corners)}")
    print(f"found_cornerless_class={not global_best.corners}")
    print(f"mathematical_profile_sha256={payload['mathematical_profile_sha256']}")
    print(f"wrote {output}")


if __name__ == "__main__":
    main()
