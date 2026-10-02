#!/usr/bin/env python3
"""Enumerate the eight-vertex induced-cube graph frontier.

A minimal eight-facet counterexample to relative linear quotients must have
a connected relative facet graph of minimum degree at least two.  This
script enumerates, up to graph isomorphism, all connected bipartite graphs
with those properties that occur as induced subgraphs of a Boolean cube.

The calculation is dimension-free.  Connected bipartite graphs have a
unique bipartition up to exchange, so a bipartite adjacency matrix is
canonicalized by column permutations and sorting its rows (and by
transposition in the balanced case).  A cube embedding is then equivalent
to a family of coordinate cuts that partitions the edges and separates
every nonedge at least twice.

Only the Python standard library is used.
"""

from __future__ import annotations

from collections import Counter
from itertools import permutations


N = 8
PAIRS = tuple((i, j) for i in range(N) for j in range(i + 1, N))
PAIR_INDEX = {pair: index for index, pair in enumerate(PAIRS)}


def graph_mask(edges: list[tuple[int, int]]) -> int:
    result = 0
    for first, second in edges:
        if first > second:
            first, second = second, first
        result |= 1 << PAIR_INDEX[(first, second)]
    return result


def edges_from_mask(mask: int) -> tuple[tuple[int, int], ...]:
    return tuple(
        pair for index, pair in enumerate(PAIRS) if (mask >> index) & 1
    )


def degrees(mask: int) -> tuple[int, ...]:
    result = [0] * N
    for first, second in edges_from_mask(mask):
        result[first] += 1
        result[second] += 1
    return tuple(result)


def connected(mask: int) -> bool:
    adjacency = [[] for _ in range(N)]
    for first, second in edges_from_mask(mask):
        adjacency[first].append(second)
        adjacency[second].append(first)
    seen = {0}
    stack = [0]
    while stack:
        vertex = stack.pop()
        for neighbor in adjacency[vertex]:
            if neighbor not in seen:
                seen.add(neighbor)
                stack.append(neighbor)
    return len(seen) == N


def matrix_rows(matrix: int, left: int, right: int) -> tuple[int, ...]:
    row_mask = (1 << right) - 1
    return tuple(
        (matrix >> (row * right)) & row_mask for row in range(left)
    )


def rows_matrix(rows: tuple[int, ...], right: int) -> int:
    return sum(row << (index * right) for index, row in enumerate(rows))


def transpose(matrix: int, size: int) -> int:
    result = 0
    for row in range(size):
        for column in range(size):
            if (matrix >> (row * size + column)) & 1:
                result |= 1 << (column * size + row)
    return result


def canonical_matrix(matrix: int, left: int, right: int) -> int:
    """Canonicalize a connected bipartite adjacency matrix."""
    original_rows = matrix_rows(matrix, left, right)

    def one_orientation(rows: tuple[int, ...]) -> int:
        best = None
        for column_permutation in permutations(range(right)):
            permuted_rows = tuple(
                sorted(
                    sum(
                        ((row >> old_column) & 1) << new_column
                        for new_column, old_column in enumerate(
                            column_permutation
                        )
                    )
                    for row in rows
                )
            )
            candidate = rows_matrix(permuted_rows, right)
            if best is None or candidate < best:
                best = candidate
        assert best is not None
        return best

    best = one_orientation(original_rows)
    if left == right:
        best = min(
            best,
            one_orientation(matrix_rows(transpose(matrix, left), left, right)),
        )
    return best


def matrix_graph_mask(matrix: int, left: int, right: int) -> int:
    return graph_mask(
        [
            (row, left + column)
            for row in range(left)
            for column in range(right)
            if (matrix >> (row * right + column)) & 1
        ]
    )


def relevant_matrices(left: int, right: int):
    """Yield all labeled matrices with minimum degree two and connected."""
    for matrix in range(1 << (left * right)):
        rows = matrix_rows(matrix, left, right)
        if any(row.bit_count() < 2 for row in rows):
            continue
        if any(
            sum((row >> column) & 1 for row in rows) < 2
            for column in range(right)
        ):
            continue
        mask = matrix_graph_mask(matrix, left, right)
        if connected(mask):
            yield matrix


def cut_is_matching(cut_edges: int, edges: tuple[tuple[int, int], ...]) -> bool:
    used_vertices = 0
    for edge_number, (first, second) in enumerate(edges):
        if not (cut_edges >> edge_number) & 1:
            continue
        endpoints = (1 << first) | (1 << second)
        if used_vertices & endpoints:
            return False
        used_vertices |= endpoints
    return True


def induced_cube_embedding(mask: int) -> tuple[int, ...] | None:
    """Return a minimum-coordinate induced cube embedding, if one exists."""
    edges = edges_from_mask(mask)
    edge_index = {edge: index for index, edge in enumerate(edges)}
    all_edges = (1 << len(edges)) - 1
    nonedges = tuple(
        pair for index, pair in enumerate(PAIRS) if not (mask >> index) & 1
    )

    # Vertex zero is fixed on the zero side, removing complementary cuts.
    cuts = []
    for small_side in range(1, 1 << (N - 1)):
        side = small_side << 1
        cut_edges = 0
        for edge, edge_number in edge_index.items():
            first, second = edge
            if ((side >> first) ^ (side >> second)) & 1:
                cut_edges |= 1 << edge_number
        if cut_edges and cut_is_matching(cut_edges, edges):
            cuts.append((cut_edges, side))

    cuts_by_edge = [[] for _ in edges]
    for cut_number, (cut_edges, _side) in enumerate(cuts):
        for edge_number in range(len(edges)):
            if (cut_edges >> edge_number) & 1:
                cuts_by_edge[edge_number].append(cut_number)

    best: list[int] | None = None

    def separates_all_nonedges(chosen: list[int]) -> bool:
        for first, second in nonedges:
            separated = sum(
                ((cuts[cut_number][1] >> first)
                 ^ (cuts[cut_number][1] >> second))
                & 1
                for cut_number in chosen
            )
            if separated < 2:
                return False
        return True

    def exact_cover(covered: int, chosen: list[int]) -> None:
        nonlocal best
        if best is not None and len(chosen) >= len(best):
            return
        if covered == all_edges:
            if separates_all_nonedges(chosen):
                best = chosen.copy()
            return

        missing = all_edges ^ covered
        best_options = None
        for edge_number in range(len(edges)):
            if not (missing >> edge_number) & 1:
                continue
            options = [
                cut_number
                for cut_number in cuts_by_edge[edge_number]
                if not (cuts[cut_number][0] & covered)
            ]
            if not options:
                return
            if best_options is None or len(options) < len(best_options):
                best_options = options
        assert best_options is not None
        for cut_number in best_options:
            exact_cover(
                covered | cuts[cut_number][0], chosen + [cut_number]
            )

    exact_cover(0, [])
    if best is None:
        return None

    vertices = []
    for vertex in range(N):
        cube_vertex = sum(
            1 << coordinate
            for coordinate, cut_number in enumerate(best)
            if (cuts[cut_number][1] >> vertex) & 1
        )
        vertices.append(cube_vertex)
    return tuple(vertices)


def format_vertices(vertices: tuple[int, ...]) -> str:
    width = max(1, max(vertices).bit_length())
    return "{" + ",".join(f"{vertex:0{width}b}" for vertex in vertices) + "}"


def main() -> int:
    abstract_types: dict[tuple[int, int, int], int] = {}
    labeled_counts = {}
    for left in range(1, N // 2 + 1):
        right = N - left
        labeled = 0
        for matrix in relevant_matrices(left, right):
            labeled += 1
            canonical = canonical_matrix(matrix, left, right)
            key = (left, right, canonical)
            abstract_types.setdefault(
                key, matrix_graph_mask(canonical, left, right)
            )
        labeled_counts[(left, right)] = labeled

    cube_types = []
    for key, mask in abstract_types.items():
        embedding = induced_cube_embedding(mask)
        if embedding is not None:
            assert len(set(embedding)) == N
            realized = graph_mask([
                (i, j) for i, j in PAIRS
                if (embedding[i] ^ embedding[j]).bit_count() == 1
            ])
            assert realized == mask, "reported embedding changes graph adjacency"
            cube_types.append((key, mask, embedding))

    print("=== Eight-vertex induced cube graphs with minimum degree two ===")
    print(f"labeled candidates by bipartition: {labeled_counts}")
    abstract_counts = Counter(key[:2] for key in abstract_types)
    print(f"abstract candidates by bipartition: {dict(sorted(abstract_counts.items()))}")
    print(f"abstract connected bipartite types tested: {len(abstract_types)}")
    print(f"induced cube types: {len(cube_types)}")
    for number, (key, mask, embedding) in enumerate(
        sorted(
            cube_types,
            key=lambda item: (
                len(edges_from_mask(item[1])),
                sorted(degrees(item[1]), reverse=True),
                item[0],
            ),
        ),
        start=1,
    ):
        left, right, _matrix = key
        print(
            f"type {number}: bipartition={left}+{right}, "
            f"edges={len(edges_from_mask(mask))}, "
            f"degrees={sorted(degrees(mask), reverse=True)}, "
            f"edge_set={edges_from_mask(mask)}, "
            f"embedding={format_vertices(embedding)}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
