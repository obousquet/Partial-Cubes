"""Interactive graph views for the Partial Cubes database."""

from collections import defaultdict, deque
from typing import Any


CATEGORY_STYLES = {
    "core": ("#355070", "#E7EEF5"),
    "metric_cubical": ("#2A6F97", "#E0F2F7"),
    "expansion_boundary": ("#9C6644", "#F4E9DF"),
    "set_family": ("#6A994E", "#EAF4E2"),
    "sign_vector_com": ("#7B2CBF", "#F1E5FA"),
    "symmetry": ("#BC6C25", "#FAEDCD"),
    "named_family": ("#6C757D", "#F1F3F5"),
}

INCLUSION_COLOR = "#475569"
STRICT_INCLUSION_COLOR = "#1D4ED8"


def _closure_style(entry: dict[str, Any]) -> dict[str, Any]:
    """Encode the three closure statuses independently on a class node."""
    p_status = entry["p_closure_status"]
    c_status = entry["c_closure_status"]
    pc_status = entry["pc_closure_status"]
    return {
        # Shape records projection closure.
        "shape": "box" if p_status == "closed" else "ellipse" if p_status == "not_closed" else "hexagon",
        # Line style records conditioning closure.
        "style": "filled" if c_status == "closed" else "filled,dashed" if c_status == "not_closed" else "filled,dotted",
        # A second outline records full pc-minor closure; three marks unresolved/N/A.
        "peripheries": 2 if pc_status == "closed" else 1 if pc_status == "not_closed" else 3,
    }


def _class_ref(entry: dict[str, Any]) -> str:
    return f'#classes/{entry["short_name"]}'


def _operation_name(cache, reference: str) -> str:
    _, operation = cache.lookup(reference)
    return operation.get("name", reference) if operation else reference


def _nodes(cache) -> tuple[list[dict[str, Any]], dict[str, str]]:
    nodes = []
    aliases = {}
    for entry in cache.get_table_entries("classes"):
        reference = _class_ref(entry)
        aliases[reference] = reference
        aliases[f'#classes/{entry["id"]}'] = reference
        color, fillcolor = CATEGORY_STYLES.get(
            entry.get("category"), ("#6C757D", "#F8F9FA")
        )
        node = {
            "id": reference,
            "ref": reference,
            "label": entry.get("graph_label") or entry["name"],
            "type": entry.get("category", "class"),
            "color": color,
            "fillcolor": fillcolor,
        }
        node.update(_closure_style(entry))
        nodes.append(node)
    return nodes, aliases


def _has_path(
    source: str,
    target: str,
    adjacency: dict[str, set[str]],
    omitted_edge: tuple[str, str],
) -> bool:
    pending = [source]
    seen = {source}
    while pending:
        current = pending.pop()
        for successor in adjacency.get(current, set()):
            if (current, successor) == omitted_edge:
                continue
            if successor == target:
                return True
            if successor not in seen:
                seen.add(successor)
                pending.append(successor)
    return False


def _reduced_inclusions(relations, aliases):
    """Return established inclusions after a conservative transitive reduction."""
    by_pair = {}
    for relation in relations:
        if (
            relation.get("status") == "established"
            and relation.get("relation_type") in {"inclusion", "strict_inclusion"}
        ):
            subject = aliases.get(relation["subject_id"])
            object_ = aliases.get(relation["object_id"])
            if subject and object_ and subject != object_:
                pair = (subject, object_)
                previous = by_pair.get(pair)
                if previous is None or relation["relation_type"] == "strict_inclusion":
                    by_pair[pair] = relation

    adjacency = defaultdict(set)
    for subject, object_ in by_pair:
        adjacency[subject].add(object_)
    return [
        relation
        for pair, relation in by_pair.items()
        if not _has_path(pair[0], pair[1], adjacency, pair)
    ]


def _assign_ranks(nodes, relations, aliases):
    """Put broad classes above their established subclasses."""
    children = defaultdict(set)
    indegree = {node["id"]: 0 for node in nodes}
    for relation in relations:
        subject = aliases.get(relation["subject_id"])
        object_ = aliases.get(relation["object_id"])
        if subject and object_ and subject not in children[object_]:
            children[object_].add(subject)
            indegree[subject] += 1
    queue = deque((node, 0) for node, degree in indegree.items() if degree == 0)
    ranks = {node: 0 for node in indegree}
    while queue:
        node, rank = queue.popleft()
        ranks[node] = max(ranks[node], rank)
        for child in children[node]:
            indegree[child] -= 1
            ranks[child] = max(ranks[child], rank + 1)
            if indegree[child] == 0:
                queue.append((child, ranks[child]))
    for node in nodes:
        node["rank"] = ranks[node["id"]]


def _relation_edges(relations, aliases, reduced):
    reduced_ids = {relation["id"] for relation in reduced}
    edges = []
    for relation in relations:
        source = aliases.get(relation.get("object_id"))
        target = aliases.get(relation.get("subject_id"))
        if not source or not target or source == target:
            continue
        relation_type = relation.get("relation_type")
        status = relation.get("status")
        reference = f'#relations/{relation["short_name"]}'
        if relation["id"] in reduced_ids:
            strict = relation_type == "strict_inclusion"
            edges.append({
                "source": source,
                "target": target,
                "ref": reference,
                "label": "strict inclusion" if strict else "inclusion",
                "color": STRICT_INCLUSION_COLOR if strict else INCLUSION_COLOR,
                "penwidth": 2.4 if strict else 1.4,
                "arrowhead": "normal" if strict else "vee",
            })
        elif status in {"open", "conjectured", "needs_verification"} and relation_type in {
            "inclusion", "strict_inclusion", "comparison"
        }:
            edges.append({
                "source": source,
                "target": target,
                "ref": reference,
                "label": status,
                "color": "#D97706",
                "style": "dashed",
                "arrowhead": "open",
                "constraint": False,
            })
        elif status == "established" and relation_type in {"noncontainment", "incomparable"}:
            edges.append({
                "source": target,
                "target": source,
                "ref": reference,
                "label": relation_type,
                "color": "#C1121F",
                "style": "dashed",
                "arrowhead": "tee",
                "dir": "both" if relation_type == "incomparable" else "forward",
                "constraint": False,
            })
        elif status == "established" and relation_type == "equality":
            edges.append({
                "source": source,
                "target": target,
                "ref": reference,
                "label": "equality",
                "color": "#5A189A",
                "style": "dotted",
                "dir": "both",
                "arrowhead": "none",
                "constraint": False,
            })
    return edges


def _legend(include_closure=False):
    legend = [
        {
            "type": "edge", "label": "", "text": "Established inclusion",
            "color": INCLUSION_COLOR, "arrowhead": "vee", "penwidth": 1.4,
        },
        {
            "type": "edge", "label": "", "text": "Established strict inclusion",
            "color": STRICT_INCLUSION_COLOR, "arrowhead": "normal", "penwidth": 2.4,
        },
        {
            "type": "edge", "label": "", "text": "Open or conjectured comparison",
            "color": "#D97706", "style": "dashed", "arrowhead": "open",
        },
        {
            "type": "edge", "label": "", "text": "Established negative comparison",
            "color": "#C1121F", "style": "dashed", "arrowhead": "tee",
        },
        {
            "type": "edge", "label": "", "text": "Established equality",
            "color": "#5A189A", "style": "dotted", "arrowhead": "none", "dir": "both",
        },
        {
            "type": "node", "label": "P", "text": "Rectangular node: P-closed",
            "shape": "box", "color": "#334155", "fillcolor": "#F8FAFC",
        },
        {
            "type": "node", "label": "not P", "text": "Elliptic node: not P-closed",
            "shape": "ellipse", "color": "#334155", "fillcolor": "#F8FAFC",
        },
        {
            "type": "node", "label": "C", "text": "Solid outline: C-closed",
            "shape": "box", "color": "#334155", "fillcolor": "#F8FAFC",
        },
        {
            "type": "node", "label": "not C", "text": "Dashed outline: not C-closed",
            "shape": "box", "color": "#334155", "fillcolor": "#F8FAFC", "style": "filled,dashed",
        },
        {
            "type": "node", "label": "PC", "text": "Double outline: PC-closed",
            "shape": "box", "color": "#334155", "fillcolor": "#F8FAFC", "peripheries": 2,
        },
        {
            "type": "node", "label": "P?", "text": "Hexagonal node: P closure open, unassessed, or N/A",
            "shape": "hexagon", "color": "#334155", "fillcolor": "#F8FAFC",
        },
        {
            "type": "node", "label": "C?", "text": "Dotted outline: C closure open, unassessed, or N/A",
            "shape": "box", "color": "#334155", "fillcolor": "#F8FAFC", "style": "filled,dotted",
        },
        {
            "type": "node", "label": "PC?", "text": "Triple outline: PC closure open, unassessed, or N/A",
            "shape": "box", "color": "#334155", "fillcolor": "#F8FAFC", "peripheries": 3,
        },
    ]
    if include_closure:
        legend.append({
            "type": "edge", "label": "", "text": "Exact operation closure",
            "color": "#0077B6", "style": "dotted", "arrowhead": "vee",
        })
    return legend


def generate_class_hierarchy(cache):
    nodes, aliases = _nodes(cache)
    relations = cache.get_table_entries("relations")
    reduced = _reduced_inclusions(relations, aliases)
    _assign_ranks(nodes, reduced, aliases)
    return {
        "nodes": nodes,
        "edges": _relation_edges(relations, aliases, reduced),
        "legend": _legend(),
        "layout": {"ranksep": 1.1, "nodesep": 0.35, "mclimit": 4},
    }


def generate_minor_closure_map(cache):
    graph = generate_class_hierarchy(cache)
    aliases = {
        reference: node["id"]
        for node in graph["nodes"]
        for reference in (node["id"],)
    }
    for entry in cache.get_table_entries("classes"):
        aliases[f'#classes/{entry["id"]}'] = _class_ref(entry)
    for result in cache.get_table_entries("operation_results"):
        if result.get("status") != "established" or result.get("result_type") != "closure_equals":
            continue
        source = aliases.get(result.get("class_id"))
        target = aliases.get(result.get("target_class_id"))
        if source and target:
            graph["edges"].append({
                "source": source,
                "target": target,
                "ref": f'#operation_results/{result["short_name"]}',
                "label": _operation_name(cache, result["operation_id"]),
                "color": "#0077B6",
                "style": "dotted",
                "arrowhead": "vee",
                "constraint": False,
            })
    graph["legend"] = _legend(include_closure=True)
    return graph
