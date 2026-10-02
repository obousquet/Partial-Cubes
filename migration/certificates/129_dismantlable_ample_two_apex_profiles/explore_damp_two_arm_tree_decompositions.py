"""Test uniqueness of the small two-arm decompositions used in the paper.

A decomposition of a tree T consists of a connected core L such that T-L
has two components of order three.  Each component, together with its edge
to L, must be a pendant path of length three.  The square-support marker of
an arm is the pair formed by the attachment edge and the terminal edge.

The final audit treats the one-unit ``uncounted vertex'' case.  It adds one
vertex to a core-of-order-at-most-two two-arm tree, retains exactly the
graphs that are ample partial cubes, and compares every sharp and one-unit
decomposition of the same graph.  This is a bounded diagnostic for the
structural lemma in the manuscript, not an order-by-order obstruction
census.
"""

from itertools import combinations, combinations_with_replacement


def canonical_edge(u: int, v: int) -> tuple[int, int]:
    return (u, v) if u < v else (v, u)


Graph = dict[int, set[int]]


def add_edge(graph: Graph, u: int, v: int) -> None:
    graph[u].add(v)
    graph[v].add(u)


def connected(graph: Graph, vertices: set[int]) -> bool:
    if not vertices:
        return False
    reached = {next(iter(vertices))}
    frontier = list(reached)
    while frontier:
        u = frontier.pop()
        for v in graph[u] & vertices:
            if v not in reached:
                reached.add(v)
                frontier.append(v)
    return reached == vertices


def components(graph: Graph, vertices: set[int]) -> list[set[int]]:
    unseen = set(vertices)
    output = []
    while unseen:
        component = {next(iter(unseen))}
        frontier = list(component)
        unseen -= component
        while frontier:
            u = frontier.pop()
            new = (graph[u] & unseen)
            unseen -= new
            component |= new
            frontier.extend(new)
        output.append(component)
    return output


def graph_edges(graph: Graph) -> list[tuple[int, int]]:
    return [
        canonical_edge(u, v)
        for u in graph
        for v in graph[u]
        if u < v
    ]


def all_distances(graph: Graph) -> dict[int, dict[int, int]]:
    output = {}
    for source in graph:
        distance = {source: 0}
        frontier = [source]
        for u in frontier:
            for v in graph[u]:
                if v not in distance:
                    distance[v] = distance[u] + 1
                    frontier.append(v)
        if len(distance) != len(graph):
            raise ValueError("graph is disconnected")
        output[source] = distance
    return output


def bipartition(graph: Graph) -> dict[int, int] | None:
    colour: dict[int, int] = {}
    for source in graph:
        if source in colour:
            continue
        colour[source] = 0
        frontier = [source]
        for u in frontier:
            for v in graph[u]:
                wanted = 1 - colour[u]
                if v in colour and colour[v] != wanted:
                    return None
                if v not in colour:
                    colour[v] = wanted
                    frontier.append(v)
    return colour


def partial_cube_embedding(
    graph: Graph,
) -> tuple[dict[int, int], dict[tuple[int, int], int]] | None:
    """Return the canonical Theta-class embedding, or ``None``."""

    if bipartition(graph) is None:
        return None
    distance = all_distances(graph)
    edges = graph_edges(graph)

    def theta(first: tuple[int, int], second: tuple[int, int]) -> bool:
        a, b = first
        x, y = second
        return (
            distance[a][x] + distance[b][y]
            != distance[a][y] + distance[b][x]
        )

    classes: list[set[tuple[int, int]]] = []
    unseen = set(edges)
    while unseen:
        edge = min(unseen)
        block = {other for other in edges if theta(edge, other)}
        block.add(edge)
        # In a partial cube the Djokic relation is already transitive.
        if any(
            theta(first, second) != (second in block)
            for first in block
            for second in edges
        ):
            return None
        classes.append(block)
        unseen -= block

    edge_class = {
        edge: index
        for index, block in enumerate(classes)
        for edge in block
    }
    root = min(graph)
    words = {vertex: 0 for vertex in graph}
    for index, block in enumerate(classes):
        a, b = min(block)
        if distance[root][a] > distance[root][b]:
            a, b = b, a
        for vertex in graph:
            if distance[vertex][b] < distance[vertex][a]:
                words[vertex] |= 1 << index

    if len(set(words.values())) != len(graph):
        return None
    for u in graph:
        for v in graph:
            adjacent = v in graph[u]
            hamming_one = (words[u] ^ words[v]).bit_count() == 1
            if adjacent != hamming_one:
                return None
            if distance[u][v] != (words[u] ^ words[v]).bit_count():
                return None
    return words, edge_class


def is_ample_embedding(words: dict[int, int], dimension: int) -> bool:
    family = tuple(words.values())
    shattered = 0
    for support in range(1 << dimension):
        traces = {word & support for word in family}
        if len(traces) == 1 << support.bit_count():
            shattered += 1
    return shattered == len(family)


def support_classes(
    support: frozenset[tuple[int, int]],
    edge_class: dict[tuple[int, int], int],
) -> frozenset[int]:
    return frozenset(edge_class[edge] for edge in support)


def contains_three_matching(pairs: set[frozenset[int]]) -> bool:
    return any(
        first.isdisjoint(second)
        and first.isdisjoint(third)
        and second.isdisjoint(third)
        for first, second, third in combinations(pairs, 3)
    )


def build_tree(core_order: int, attachment_vertices: tuple[int, int]) -> Graph:
    graph = {vertex: set() for vertex in range(core_order + 6)}
    for vertex in range(core_order - 1):
        add_edge(graph, vertex, vertex + 1)
    next_vertex = core_order
    for root in attachment_vertices:
        first, middle, terminal = range(next_vertex, next_vertex + 3)
        next_vertex += 3
        add_edge(graph, root, first)
        add_edge(graph, first, middle)
        add_edge(graph, middle, terminal)
    return graph


def build_extra_owner_tree(
    core_order: int,
    small_root: int,
    large_root: int,
    large_shape: str,
    collar_orbit: str,
) -> Graph:
    graph = {vertex: set() for vertex in range(core_order + 7)}
    for vertex in range(core_order - 1):
        add_edge(graph, vertex, vertex + 1)

    first, middle, terminal = range(core_order, core_order + 3)
    add_edge(graph, small_root, first)
    add_edge(graph, first, middle)
    add_edge(graph, middle, terminal)

    large = list(range(core_order + 3, core_order + 7))
    if large_shape == "path":
        add_edge(graph, large[0], large[1])
        add_edge(graph, large[1], large[2])
        add_edge(graph, large[2], large[3])
        collar = large[0] if collar_orbit == "end" else large[1]
    elif large_shape == "star":
        add_edge(graph, large[0], large[1])
        add_edge(graph, large[0], large[2])
        add_edge(graph, large[0], large[3])
        collar = large[0] if collar_orbit == "center" else large[1]
    else:
        raise ValueError(large_shape)
    add_edge(graph, large_root, collar)
    return graph


def decompositions(tree: Graph) -> list[dict]:
    core_order = len(tree) - 6
    output = []
    for core_tuple in combinations(tree, core_order):
        core = set(core_tuple)
        if not connected(tree, core):
            continue
        exterior = set(tree) - core
        exterior_components = components(tree, exterior)
        if sorted(map(len, exterior_components)) != [3, 3]:
            continue

        arms = []
        valid = True
        for component in exterior_components:
            attachment_edges = [
                (u, v)
                for u in core
                for v in component
                if v in tree[u]
            ]
            if len(attachment_edges) != 1:
                valid = False
                break
            core_vertex, first = attachment_edges[0]
            component_degrees = {
                vertex: len(tree[vertex] & component)
                for vertex in component
            }
            if sorted(component_degrees.values()) != [1, 1, 2]:
                valid = False
                break
            if component_degrees[first] != 1:
                valid = False
                break
            middle = next(iter(tree[first] & component))
            terminal = next(iter((tree[middle] & component) - {first}))
            attachment_edge = canonical_edge(core_vertex, first)
            middle_edge = canonical_edge(first, middle)
            terminal_edge = canonical_edge(middle, terminal)
            arms.append(
                {
                    "vertices": (core_vertex, first, middle, terminal),
                    "support": frozenset({attachment_edge, terminal_edge}),
                    "edges": frozenset({attachment_edge, middle_edge, terminal_edge}),
                }
            )
        if valid:
            output.append(
                {
                    "core": frozenset(core),
                    "supports": frozenset(arm["support"] for arm in arms),
                    "arm_edges": frozenset().union(*(arm["edges"] for arm in arms)),
                    "arms": tuple(arms),
                }
            )
    return output


def extra_owner_decompositions(tree: Graph) -> list[dict]:
    core_order = len(tree) - 7
    output = []
    for core_tuple in combinations(tree, core_order):
        core = set(core_tuple)
        if not connected(tree, core):
            continue
        exterior_components = components(tree, set(tree) - core)
        if sorted(map(len, exterior_components)) != [3, 4]:
            continue

        attachment_data = {}
        valid = True
        for component in exterior_components:
            attachment_edges = [
                (u, v)
                for u in core
                for v in component
                if v in tree[u]
            ]
            if len(attachment_edges) != 1:
                valid = False
                break
            attachment_data[len(component)] = (
                component,
                attachment_edges[0],
            )
        if not valid:
            continue

        small_component, (core_vertex, first) = attachment_data[3]
        small_degrees = {
            vertex: len(tree[vertex] & small_component)
            for vertex in small_component
        }
        if sorted(small_degrees.values()) != [1, 1, 2]:
            continue
        if small_degrees[first] != 1:
            continue
        middle = next(iter(tree[first] & small_component))
        terminal = next(iter((tree[middle] & small_component) - {first}))
        small_support = frozenset(
            {
                canonical_edge(core_vertex, first),
                canonical_edge(middle, terminal),
            }
        )

        large_component, (_, collar) = attachment_data[4]
        large_degrees = sorted(
            len(tree[vertex] & large_component)
            for vertex in large_component
        )
        if large_degrees not in (
            [1, 1, 2, 2],
            [1, 1, 1, 3],
            [2, 2, 2, 2],
        ):
            continue
        output.append(
            {
                "core": frozenset(core),
                "small_support": small_support,
                "small_edges": (
                    canonical_edge(core_vertex, first),
                    canonical_edge(first, middle),
                    canonical_edge(middle, terminal),
                ),
                "large_collar_degree": len(tree[collar] & large_component),
                "large_degrees": tuple(large_degrees),
            }
        )
    return output


def induced_without(graph: Graph, deleted: int) -> Graph:
    return {
        vertex: neighbours - {deleted}
        for vertex, neighbours in graph.items()
        if vertex != deleted
    }


def one_vertex_two_arm_decompositions(
    graph: Graph,
    edge_class: dict[tuple[int, int], int],
) -> list[dict]:
    """Find all eta=1 decompositions of an ample graph."""

    core_order = len(graph) - 7
    if core_order not in (1, 2):
        return []
    output = []
    for extra in graph:
        tree = induced_without(graph, extra)
        if len(graph_edges(tree)) != len(tree) - 1 or not connected(
            tree, set(tree)
        ):
            continue
        for decomp in decompositions(tree):
            attachment_edges = [arm["vertices"][:2] for arm in decomp["arms"]]
            attachment_edges = [
                canonical_edge(*edge) for edge in attachment_edges
            ]
            attachment_classes = [edge_class[edge] for edge in attachment_edges]
            if len(set(attachment_classes)) != 2:
                continue
            # In the passive-gate normal form the two attachment edges
            # have the distinct passive directions z_1,z_2.  Every other
            # edge of the deleted tree is active, so neither passive class
            # may recur there.  A recurrence on an edge incident with the
            # added vertex is allowed and is precisely a completed square.
            other_tree_edges = set(graph_edges(tree)) - set(attachment_edges)
            if any(
                edge_class[edge] in attachment_classes
                for edge in other_tree_edges
            ):
                continue
            supports = frozenset(
                support_classes(support, edge_class)
                for support in decomp["supports"]
            )
            if len(supports) != 2 or any(len(pair) != 2 for pair in supports):
                continue
            output.append(
                {
                    "extra": extra,
                    "core": decomp["core"],
                    "supports": supports,
                }
            )
    return output


def sharp_decompositions_in_graph(
    graph: Graph,
    edge_class: dict[tuple[int, int], int],
) -> list[dict]:
    if len(graph_edges(graph)) != len(graph) - 1:
        return []
    output = []
    for decomp in decompositions(graph):
        output.append(
            {
                "core": decomp["core"],
                "supports": frozenset(
                    support_classes(support, edge_class)
                    for support in decomp["supports"]
                ),
            }
        )
    return output


def extra_owner_decompositions_in_graph(
    graph: Graph,
    edge_class: dict[tuple[int, int], int],
) -> list[dict]:
    output = []
    for decomp in extra_owner_decompositions(graph):
        support = support_classes(decomp["small_support"], edge_class)
        if len(support) != 2:
            continue
        output.append(
            {
                "core": decomp["core"],
                "ordinary_support": support,
                "large_degrees": decomp["large_degrees"],
            }
        )
    return output


def add_vertex(graph: Graph, neighbours: set[int]) -> Graph:
    output = {vertex: set(adjacent) for vertex, adjacent in graph.items()}
    new_vertex = max(output, default=-1) + 1
    output[new_vertex] = set()
    for vertex in neighbours:
        add_edge(output, new_vertex, vertex)
    return output


def eta_one_graphs() -> list[tuple[Graph, dict[tuple[int, int], int]]]:
    """Generate every ample partial-cube graph in the eta=1 normal form."""

    output = []
    seen = set()
    for core_order in (1, 2):
        for attachments in combinations_with_replacement(
            range(core_order), 2
        ):
            tree = build_tree(core_order, attachments)
            vertices = tuple(tree)
            for size in range(1, len(vertices) + 1):
                for neighbours_tuple in combinations(vertices, size):
                    graph = add_vertex(tree, set(neighbours_tuple))
                    embedding = partial_cube_embedding(graph)
                    if embedding is None:
                        continue
                    words, edge_class = embedding
                    dimension = 1 + max(edge_class.values(), default=-1)
                    if not is_ample_embedding(words, dimension):
                        continue
                    # Deduplicate labelled repetitions caused by symmetric
                    # attachment choices using the distance matrix and the
                    # canonical word set.  The decomposition comparison below
                    # still rediscovers all distinguished extra vertices.
                    normalized = tuple(sorted(words.values()))
                    if normalized in seen:
                        continue
                    seen.add(normalized)
                    output.append((graph, edge_class))
    return output


def main() -> None:
    print("=== Sharp two-arm tree decomposition audit ===")
    total_configurations = 0
    multiple = []
    support_ambiguous = []

    for core_order in range(1, 4):
        configuration_count = 0
        for index, attachment_vertices in enumerate(
            combinations_with_replacement(range(core_order), 2)
        ):
            configuration_count += 1
            total_configurations += 1
            tree = build_tree(core_order, attachment_vertices)
            decomps = decompositions(tree)
            if len(decomps) > 1:
                multiple.append(
                    (core_order, attachment_vertices, tree, decomps)
                )
            support_sets = {decomp["supports"] for decomp in decomps}
            if len(support_sets) > 1:
                support_ambiguous.append(
                    (core_order, attachment_vertices, tree, decomps)
                )

        print(
            f"core_order={core_order}: configurations={configuration_count}"
        )

    print(f"total configurations checked: {total_configurations}")
    print(f"trees with multiple cores: {len(multiple)}")
    print(f"trees with different support-pair systems: {len(support_ambiguous)}")

    for core_order, attachment_vertices, tree, decomps in support_ambiguous:
        edges = sorted(
            canonical_edge(u, v)
            for u in tree
            for v in tree[u]
            if u < v
        )
        print(
            "COUNTEREXAMPLE "
            f"core_order={core_order}, attachments={attachment_vertices}"
        )
        print(f"edges={edges}")
        for decomp in decomps:
            print(
                "  core=",
                sorted(decomp["core"]),
                "supports=",
                sorted(sorted(support) for support in decomp["supports"]),
            )

    print("=== One-unit extra-owner tree audit ===")
    extra_configurations = 0
    ambiguous_small_arm = []
    sharp_extra_disjoint = []
    for core_order in range(1, 3):
        for small_root in range(core_order):
            for large_root in range(core_order):
                for large_shape, collar_orbits in (
                    ("path", ("end", "internal")),
                    ("star", ("center", "leaf")),
                ):
                    for collar_orbit in collar_orbits:
                        extra_configurations += 1
                        tree = build_extra_owner_tree(
                            core_order,
                            small_root,
                            large_root,
                            large_shape,
                            collar_orbit,
                        )
                        extra_decomps = extra_owner_decompositions(tree)
                        small_supports = {
                            decomp["small_support"]
                            for decomp in extra_decomps
                        }
                        if len(small_supports) > 1:
                            ambiguous_small_arm.append(
                                (
                                    core_order,
                                    small_root,
                                    large_root,
                                    large_shape,
                                    collar_orbit,
                                    tree,
                                    extra_decomps,
                                )
                            )

                        sharp_decomps = decompositions(tree)
                        for extra_decomp in extra_decomps:
                            for sharp_decomp in sharp_decomps:
                                if all(
                                    extra_decomp["small_support"] != support
                                    for support in sharp_decomp["supports"]
                                ):
                                    sharp_extra_disjoint.append(
                                        (
                                            core_order,
                                            small_root,
                                            large_root,
                                            large_shape,
                                            collar_orbit,
                                            tree,
                                            extra_decomp,
                                            sharp_decomp,
                                        )
                                    )

    print(f"extra-owner configurations checked: {extra_configurations}")
    print(
        "trees with ambiguous ordinary-arm support: "
        f"{len(ambiguous_small_arm)}"
    )
    print(
        "sharp/extra decompositions with no repeated ordinary support: "
        f"{len(sharp_extra_disjoint)}"
    )
    for record in ambiguous_small_arm:
        (
            core_order,
            small_root,
            large_root,
            large_shape,
            collar_orbit,
            tree,
            extra_decomps,
        ) = record
        edges = sorted(
            canonical_edge(u, v)
            for u in tree
            for v in tree[u]
            if u < v
        )
        print(
            "AMBIGUOUS_EXTRA "
            f"core_order={core_order}, roots=({small_root},{large_root}), "
            f"large={large_shape}:{collar_orbit}"
        )
        print(f"edges={edges}")
        for decomp in extra_decomps:
            print(f"  {decomp}")
    for record in sharp_extra_disjoint:
        (
            core_order,
            small_root,
            large_root,
            large_shape,
            collar_orbit,
            tree,
            extra_decomp,
            sharp_decomp,
        ) = record
        edges = sorted(
            canonical_edge(u, v)
            for u in tree
            for v in tree[u]
            if u < v
        )
        print(
            "SHARP_EXTRA_COUNTEREXAMPLE "
            f"core_order={core_order}, roots=({small_root},{large_root}), "
            f"large={large_shape}:{collar_orbit}"
        )
        print(f"edges={edges}")
        print(f"extra={extra_decomp}")
        print(f"sharp={sharp_decomp}")

    print("=== One-unit uncounted-vertex audit ===")
    graphs = eta_one_graphs()
    eta_decomposition_total = 0
    eta_eta_unresolved = []
    eta_sharp_unresolved = []
    eta_extra_unresolved = []
    eta_eta_matching_records = []
    eta_eta_outcomes = {"repeated": 0, "matching3": 0, "exception": 0}
    eta_sharp_outcomes = {"repeated": 0, "matching3": 0}
    eta_extra_outcomes = {"repeated": 0, "matching3": 0, "unresolved": 0}
    profiles: dict[tuple[int, int, int], int] = {}
    for graph, edge_class in graphs:
        eta_decomps = one_vertex_two_arm_decompositions(graph, edge_class)
        sharp_decomps = sharp_decompositions_in_graph(graph, edge_class)
        extra_decomps = extra_owner_decompositions_in_graph(
            graph, edge_class
        )
        eta_decomposition_total += len(eta_decomps)
        squares = {
            frozenset(
                {
                    edge_class[canonical_edge(vertex, first)],
                    edge_class[canonical_edge(vertex, second)],
                }
            )
            for vertex in graph
            for first, second in combinations(graph[vertex], 2)
            if len(graph[first] & graph[second] - {vertex}) > 0
        }
        profile = (
            len(graph),
            len(graph[max(graph)]),
            len(squares),
        )
        profiles[profile] = profiles.get(profile, 0) + 1

        for first, second in combinations(eta_decomps, 2):
            pairs = set(first["supports"]) | set(second["supports"])
            if not set(first["supports"]).isdisjoint(second["supports"]):
                eta_eta_outcomes["repeated"] += 1
            elif contains_three_matching(pairs):
                eta_eta_outcomes["matching3"] += 1
                eta_eta_matching_records.append(
                    (graph, edge_class, first, second, pairs)
                )
            else:
                eta_eta_outcomes["exception"] += 1
                eta_eta_unresolved.append(
                    (graph, edge_class, first, second, pairs)
                )
        for eta_decomp in eta_decomps:
            for sharp_decomp in sharp_decomps:
                pairs = set(eta_decomp["supports"]) | set(
                    sharp_decomp["supports"]
                )
                if not set(eta_decomp["supports"]).isdisjoint(
                    sharp_decomp["supports"]
                ):
                    eta_sharp_outcomes["repeated"] += 1
                elif contains_three_matching(pairs):
                    eta_sharp_outcomes["matching3"] += 1
                else:
                    eta_sharp_unresolved.append(
                        (graph, edge_class, eta_decomp, sharp_decomp, pairs)
                    )
            for extra_decomp in extra_decomps:
                ordinary = extra_decomp["ordinary_support"]
                pairs = set(eta_decomp["supports"]) | {ordinary}
                if ordinary in eta_decomp["supports"]:
                    eta_extra_outcomes["repeated"] += 1
                elif all(
                    ordinary.isdisjoint(pair)
                    for pair in eta_decomp["supports"]
                ):
                    eta_extra_outcomes["matching3"] += 1
                else:
                    eta_extra_outcomes["unresolved"] += 1
                    eta_extra_unresolved.append(
                        (graph, edge_class, eta_decomp, extra_decomp, pairs)
                    )

    print(f"ample partial-cube graphs generated: {len(graphs)}")
    print(f"eta decompositions recovered: {eta_decomposition_total}")
    print(f"(order, added-degree, square-count) profiles: {profiles}")
    print(f"eta/eta comparison outcomes: {eta_eta_outcomes}")
    print(f"eta/sharp comparison outcomes: {eta_sharp_outcomes}")
    print(f"eta/extra comparison outcomes: {eta_extra_outcomes}")
    print(
        "eta/eta pairs with distinct supports and no three-matching: "
        f"{len(eta_eta_unresolved)}"
    )
    print(
        "eta/sharp pairs with distinct supports and no three-matching: "
        f"{len(eta_sharp_unresolved)}"
    )
    for kind, records in (
        ("ETA_ETA_MATCHING3", eta_eta_matching_records),
        ("ETA_ETA_UNRESOLVED", eta_eta_unresolved),
        ("ETA_SHARP_UNRESOLVED", eta_sharp_unresolved),
        ("ETA_EXTRA_UNRESOLVED", eta_extra_unresolved),
    ):
        for graph, edge_class, first, second, pairs in records:
            print(kind)
            print(f"  edges={sorted(graph_edges(graph))}")
            print(f"  theta={edge_class}")
            print(f"  first={first}")
            print(f"  second={second}")
            print(f"  support_pairs={sorted(map(sorted, pairs))}")


if __name__ == "__main__":
    main()
