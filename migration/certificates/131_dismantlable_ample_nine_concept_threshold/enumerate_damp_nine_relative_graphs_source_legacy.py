"""Enumerate the dimension-free nine-vertex induced-cube graph frontier.

The graph calculation is the nine-vertex analogue of YangDim's proved
eight-facet classifier.  A connected bipartite graph embeds as an induced
subgraph of a Boolean cube exactly when its edges can be partitioned into
matching cuts and every nonedge is separated by at least two cuts.

Only the abstract graph frontier is classified here.  Boundary realizability
and the terminal ample-interval conditions are separate mathematical steps.
"""

from __future__ import annotations

import json
import sys
from itertools import permutations
from pathlib import Path


YANGDIM_SCRIPTS = Path(__file__).resolve().parents[2] / "YangDim" / "scripts"
sys.path.insert(0, str(YANGDIM_SCRIPTS))

import enumerate_eight_deletion_graphs as cube_graphs  # noqa: E402


N = 9
cube_graphs.N = N
cube_graphs.PAIRS = tuple((i, j) for i in range(N) for j in range(i + 1, N))
cube_graphs.PAIR_INDEX = {
    pair: index for index, pair in enumerate(cube_graphs.PAIRS)
}


def canonical_bipartite_key(
    matrix: int, left: int, right: int
) -> tuple[int, ...]:
    """Canonicalize by permuting only the smaller bipartition.

    Once the ``left`` vertices are ordered, each right vertex is determined
    by its neighborhood bit mask on the left.  Sorting those masks removes
    every permutation of the right side.  Enumerating the at most ``4!``
    orders of the left side therefore gives an exact canonical form, instead
    of the ``right!`` loop used by the eight-vertex reference script.
    """

    rows = cube_graphs.matrix_rows(matrix, left, right)
    columns = tuple(
        sum(((rows[row] >> column) & 1) << row for row in range(left))
        for column in range(right)
    )
    return min(
        tuple(
            sorted(
                sum(((column >> old) & 1) << new for new, old in enumerate(order))
                for column in columns
            )
        )
        for order in permutations(range(left))
    )


def enumerate_types() -> dict[str, object]:
    abstract_types: dict[tuple[int, int, tuple[int, ...]], int] = {}
    labeled_counts: dict[str, int] = {}
    for left in range(1, N // 2 + 1):
        right = N - left
        labeled = 0
        for matrix in cube_graphs.relevant_matrices(left, right):
            labeled += 1
            canonical = canonical_bipartite_key(matrix, left, right)
            key = (left, right, canonical)
            abstract_types.setdefault(key, cube_graphs.matrix_graph_mask(matrix, left, right))
        labeled_counts[f"{left}+{right}"] = labeled

    cube_types = []
    for key, mask in abstract_types.items():
        embedding = cube_graphs.induced_cube_embedding(mask)
        if embedding is None:
            continue
        cube_types.append((key, mask, embedding))

    rows = []
    for number, (key, mask, embedding) in enumerate(
        sorted(
            cube_types,
            key=lambda item: (
                len(cube_graphs.edges_from_mask(item[1])),
                sorted(cube_graphs.degrees(item[1]), reverse=True),
                item[0],
            ),
        ),
        start=1,
    ):
        left, right, _matrix = key
        rows.append(
            {
                "type": number,
                "bipartition": [left, right],
                "edge_count": len(cube_graphs.edges_from_mask(mask)),
                "degree_sequence": sorted(
                    cube_graphs.degrees(mask), reverse=True
                ),
                "edges": [list(edge) for edge in cube_graphs.edges_from_mask(mask)],
                "embedding": list(embedding),
                "embedding_dimension": max(embedding).bit_length(),
            }
        )
    return {
        "status": "VERIFIED",
        "vertex_count": N,
        "labeled_candidates_by_bipartition": labeled_counts,
        "abstract_connected_bipartite_types_tested": len(abstract_types),
        "induced_cube_type_count": len(rows),
        "types": rows,
    }


if __name__ == "__main__":
    print(json.dumps(enumerate_types(), indent=2, sort_keys=True))
