#!/usr/bin/env python3
"""Explore the corner-exchange neighborhood of Hall's obstruction.

The computation has four bounded phases.

``source``
    Enumerate every one-concept ample extension of one of five fixed
    cornerless maximum classes: Hall's class ``H`` and its four first
    nonpeelable exchanges ``A``--``D``.  Classify all resulting maximum
    exchanges by corners and corner peelability.

``merge-discovery``
    Validate the five complete source profiles and construct the
    radius-two nonpeelable exchange graph.  Canonicalize every vertex and
    edge under cube isomorphism.

``audit``
    For a half-open interval of graph vertices, verify that every
    elementary section and contraction is corner-peelable.

``merge-audits``
    Validate an exact audit cover.  The audited cornerless vertices are
    forbidden pc-minors.  Every graph edge is an ample one-concept
    extension of either endpoint, so the one-concept extension lemma makes
    it a forbidden pc-minor as well.

The script deliberately stops at exchange radius two.  It proves a finite
family of obstruction types; it does not claim that the exchange component
or the full forbidden-minor basis is exhausted.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter, deque
from pathlib import Path
from typing import Iterable


REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))


from explore_smallest_additional_peripheral_obstruction import (
    HALL_CONCEPTS,
    family_from_vertices,
    irredundant_projection,
    is_ample_family,
    refined_canonical_cube_family,
    section,
    shattered_rank_histogram,
    vertices,
)
from explore_dismantlable_ample import (
    contraction,
    corners,
    is_dismantlable,
)


DIMENSION = 12
AMBIENT_ORDER = 1 << DIMENSION
CANONICAL_BYTES = AMBIENT_ORDER // 8
FIRST_EXCHANGES = {
    "A": (373, 1011),
    "B": (1105, 192),
    "C": (1360, 3120),
    "D": (1525, 252),
}
SOURCE_LABELS = ("H", "A", "B", "C", "D")


def encoded_digest(value: object) -> str:
    """Return the deterministic mathematical-profile digest."""

    return hashlib.sha256(
        json.dumps(
            value,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    ).hexdigest()


def write_payload(path: Path, mathematics: dict[str, object]) -> None:
    """Write one deterministic payload with its mathematical digest."""

    payload = {
        **mathematics,
        "mathematical_profile_sha256": encoded_digest(mathematics),
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(f"wrote {path}", flush=True)


def load_payload(path: Path) -> dict[str, object]:
    """Load and verify one deterministic payload."""

    payload = json.loads(path.read_text(encoding="utf-8"))
    recorded = payload.pop("mathematical_profile_sha256")
    actual = encoded_digest(payload)
    if recorded != actual:
        raise ValueError(
            f"{path}: mathematical digest mismatch: "
            f"{recorded} != {actual}"
        )
    payload["mathematical_profile_sha256"] = recorded
    return payload


def hall_family() -> frozenset[int]:
    return frozenset(HALL_CONCEPTS)


def source_family(label: str) -> frozenset[int]:
    """Return one of the five fixed cornerless source classes."""

    hall = hall_family()
    if label == "H":
        return hall
    if label not in FIRST_EXCHANGES:
        raise ValueError(f"unknown source label {label!r}")
    added, removed = FIRST_EXCHANGES[label]
    return frozenset((hall | {added}) - {removed})


def family_mask(members: Iterable[int]) -> int:
    return family_from_vertices(members)


def single_vertex_extension_is_isometric(
    base: frozenset[int],
    added: int,
) -> tuple[bool, tuple[int, ...]]:
    """Test the exact geodesic condition for ``base union {added}``.

    Since ``base`` is already isometric, the extension is isometric exactly
    when, for every base vertex ``x``, some base neighbor of ``added`` lies
    one Hamming step closer to ``x``.
    """

    neighbors = tuple(
        added ^ (1 << coordinate)
        for coordinate in range(DIMENSION)
        if added ^ (1 << coordinate) in base
    )
    if not neighbors:
        return False, neighbors
    return (
        all(
            any(
                (neighbor ^ target).bit_count()
                == (added ^ target).bit_count() - 1
                for neighbor in neighbors
            )
            for target in base
        ),
        neighbors,
    )


def source_profile(label: str) -> dict[str, object]:
    """Enumerate the complete one-concept neighborhood of one source."""

    base = source_family(label)
    base_mask = family_mask(base)
    if len(base) != 299:
        raise AssertionError(f"{label}: source order is not 299")
    if not is_ample_family(base_mask, DIMENSION):
        raise AssertionError(f"{label}: source is not ample")
    if corners(base_mask, DIMENSION):
        raise AssertionError(f"{label}: source is not cornerless")
    if is_dismantlable(base_mask, DIMENSION):
        raise AssertionError(f"{label}: source is unexpectedly peelable")

    extension_rows: list[dict[str, object]] = []
    for added in range(AMBIENT_ORDER):
        if added in base:
            continue
        isometric, neighbors = single_vertex_extension_is_isometric(
            base,
            added,
        )
        if not isometric:
            continue

        extension = frozenset(base | {added})
        extension_mask = family_mask(extension)
        if not is_ample_family(extension_mask, DIMENSION):
            raise AssertionError(
                f"{label}: isometric one-concept extension is not ample"
            )
        extension_corners = corners(extension_mask, DIMENSION)
        if added not in extension_corners:
            raise AssertionError(
                f"{label}: added concept is not a corner"
            )
        extension_dismantlable = is_dismantlable(
            extension_mask,
            DIMENSION,
        )

        swap_rows: list[dict[str, object]] = []
        for removed in extension_corners:
            if removed == added:
                continue
            exchanged = frozenset(extension - {removed})
            exchanged_mask = family_mask(exchanged)
            if not is_ample_family(exchanged_mask, DIMENSION):
                raise AssertionError(
                    f"{label}: deleting an extension corner is not ample"
                )
            exchanged_corners = corners(exchanged_mask, DIMENSION)
            swap_rows.append(
                {
                    "removed": removed,
                    "removed_bits": format(removed, f"0{DIMENSION}b"),
                    "family_hex": format(exchanged_mask, "x"),
                    "corners": list(exchanged_corners),
                    "corner_count": len(exchanged_corners),
                    "dismantlable": is_dismantlable(
                        exchanged_mask,
                        DIMENSION,
                    ),
                }
            )

        extension_rows.append(
            {
                "added": added,
                "added_bits": format(added, f"0{DIMENSION}b"),
                "neighbor_count": len(neighbors),
                "neighbors": list(neighbors),
                "family_hex": format(extension_mask, "x"),
                "corners": list(extension_corners),
                "corner_count": len(extension_corners),
                "dismantlable": extension_dismantlable,
                "swaps": swap_rows,
            }
        )

    profile_counts = Counter(
        (
            int(row["corner_count"]),
            bool(row["dismantlable"]),
        )
        for row in extension_rows
    )
    swap_counts = Counter(
        (
            int(swap["corner_count"]),
            bool(swap["dismantlable"]),
        )
        for row in extension_rows
        for swap in row["swaps"]
    )
    return {
        "schema": "hall-corner-exchange-source-v1",
        "source_label": label,
        "source_family_hex": format(base_mask, "x"),
        "source_order": len(base),
        "source_corner_count": 0,
        "source_dismantlable": False,
        "extensions": extension_rows,
        "counts": {
            "outside_concepts": AMBIENT_ORDER - len(base),
            "ample_one_concept_extensions": len(extension_rows),
            "extension_profiles": {
                f"corners={corner_count},dismantlable={str(damp).lower()}": count
                for (corner_count, damp), count in sorted(
                    profile_counts.items()
                )
            },
            "labeled_exchanges": sum(
                len(row["swaps"]) for row in extension_rows
            ),
            "exchange_profiles": {
                f"corners={corner_count},dismantlable={str(damp).lower()}": count
                for (corner_count, damp), count in sorted(
                    swap_counts.items()
                )
            },
        },
    }


def canonical_profile(mask: int) -> tuple[int, str]:
    """Return the exact cube canonical form and a compact digest."""

    canonical = refined_canonical_cube_family(mask, DIMENSION)
    digest = hashlib.sha256(
        canonical.to_bytes(CANONICAL_BYTES, "little")
    ).hexdigest()
    return canonical, digest


def merge_discovery(paths: list[Path]) -> dict[str, object]:
    """Merge the five complete source profiles into the radius-two graph."""

    payloads = [load_payload(path) for path in paths]
    by_label = {
        str(payload["source_label"]): payload for payload in payloads
    }
    if set(by_label) != set(SOURCE_LABELS):
        raise ValueError(
            "source profiles must cover exactly "
            f"{SOURCE_LABELS}, got {sorted(by_label)}"
        )

    for label in SOURCE_LABELS:
        payload = by_label[label]
        counts = payload["counts"]
        expected_extensions = 16 if label == "H" else 20
        expected_exchanges = 40 if label == "H" else 44
        if counts["ample_one_concept_extensions"] != expected_extensions:
            raise ValueError(f"{label}: unexpected extension count")
        if counts["labeled_exchanges"] != expected_exchanges:
            raise ValueError(f"{label}: unexpected exchange count")
        expected_extension_profiles = (
            {
                "corners=2,dismantlable=false": 4,
                "corners=4,dismantlable=true": 12,
            }
            if label == "H"
            else {
                "corners=2,dismantlable=false": 8,
                "corners=4,dismantlable=true": 12,
            }
        )
        if counts["extension_profiles"] != expected_extension_profiles:
            raise ValueError(f"{label}: unexpected extension profile")

    vertex_masks: set[int] = set()
    edge_occurrences: dict[int, list[dict[str, object]]] = {}
    for label in SOURCE_LABELS:
        payload = by_label[label]
        source_mask = int(str(payload["source_family_hex"]), 16)
        for row in payload["extensions"]:
            if bool(row["dismantlable"]):
                continue
            extension_mask = int(str(row["family_hex"]), 16)
            extension_corners = tuple(int(x) for x in row["corners"])
            if len(extension_corners) != 2:
                raise ValueError(
                    f"{label}: nonpeelable extension does not have two corners"
                )
            endpoints = tuple(
                extension_mask & ~(1 << corner)
                for corner in extension_corners
            )
            if source_mask not in endpoints:
                raise ValueError(
                    f"{label}: source is not an endpoint of its extension"
                )
            vertex_masks.update(endpoints)
            edge_occurrences.setdefault(extension_mask, []).append(
                {
                    "source_label": label,
                    "added": int(row["added"]),
                    "corners": list(extension_corners),
                    "endpoint_family_hex": [
                        format(endpoint, "x") for endpoint in endpoints
                    ],
                }
            )

    if len(vertex_masks) != 27 or len(edge_occurrences) != 32:
        raise ValueError(
            "unexpected radius-two inventory: "
            f"vertices={len(vertex_masks)}, edges={len(edge_occurrences)}"
        )

    raw_vertices: list[dict[str, object]] = []
    for index, mask in enumerate(sorted(vertex_masks)):
        if mask.bit_count() != 299:
            raise ValueError("exchange vertex does not have order 299")
        if corners(mask, DIMENSION):
            raise ValueError("exchange vertex is not cornerless")
        if is_dismantlable(mask, DIMENSION):
            raise ValueError("exchange vertex is unexpectedly peelable")
        if shattered_rank_histogram(mask, DIMENSION) != {
            0: 1,
            1: 12,
            2: 66,
            3: 220,
        }:
            raise ValueError("exchange vertex is not maximum VC-three")
        canonical, digest = canonical_profile(mask)
        raw_vertices.append(
            {
                "temporary_index": index,
                "family_hex": format(mask, "x"),
                "canonical_family_hex": format(canonical, "x"),
                "canonical_sha256": digest,
            }
        )

    canonical_masks = {
        int(str(row["canonical_family_hex"]), 16)
        for row in raw_vertices
    }
    if len(canonical_masks) != len(raw_vertices):
        raise ValueError("two radius-two vertices are cube-isomorphic")
    raw_vertices.sort(
        key=lambda row: int(str(row["canonical_family_hex"]), 16)
    )
    vertex_id_by_mask: dict[int, str] = {}
    vertices_payload: list[dict[str, object]] = []
    for index, row in enumerate(raw_vertices):
        vertex_id = f"V{index:02d}"
        mask = int(str(row["family_hex"]), 16)
        vertex_id_by_mask[mask] = vertex_id
        vertices_payload.append(
            {
                "id": vertex_id,
                "family_hex": row["family_hex"],
                "canonical_family_hex": row["canonical_family_hex"],
                "canonical_sha256": row["canonical_sha256"],
                "order": 299,
                "corner_count": 0,
                "dismantlable": False,
                "shattered_rank_histogram": {
                    "0": 1,
                    "1": 12,
                    "2": 66,
                    "3": 220,
                },
            }
        )

    edges_payload: list[dict[str, object]] = []
    edge_canonical_masks: set[int] = set()
    adjacency: dict[str, set[str]] = {
        vertex_id: set() for vertex_id in vertex_id_by_mask.values()
    }
    for extension_mask, occurrences in sorted(edge_occurrences.items()):
        extension_corners = corners(extension_mask, DIMENSION)
        endpoints = tuple(
            extension_mask & ~(1 << corner)
            for corner in extension_corners
        )
        if len(extension_corners) != 2:
            raise ValueError("radius-two edge does not have two corners")
        if any(endpoint not in vertex_id_by_mask for endpoint in endpoints):
            raise ValueError("radius-two edge has an endpoint outside inventory")
        if is_dismantlable(extension_mask, DIMENSION):
            raise ValueError("radius-two edge is unexpectedly peelable")
        if shattered_rank_histogram(extension_mask, DIMENSION) != {
            0: 1,
            1: 12,
            2: 66,
            3: 220,
            4: 1,
        }:
            raise ValueError("radius-two edge has the wrong shattering profile")
        canonical, digest = canonical_profile(extension_mask)
        edge_canonical_masks.add(canonical)
        endpoint_ids = sorted(vertex_id_by_mask[x] for x in endpoints)
        adjacency[endpoint_ids[0]].add(endpoint_ids[1])
        adjacency[endpoint_ids[1]].add(endpoint_ids[0])
        edges_payload.append(
            {
                "family_hex": format(extension_mask, "x"),
                "canonical_family_hex": format(canonical, "x"),
                "canonical_sha256": digest,
                "order": 300,
                "corner_count": 2,
                "corners": list(extension_corners),
                "dismantlable": False,
                "endpoints": endpoint_ids,
                "source_occurrences": occurrences,
                "shattered_rank_histogram": {
                    "0": 1,
                    "1": 12,
                    "2": 66,
                    "3": 220,
                    "4": 1,
                },
            }
        )
    if len(edge_canonical_masks) != len(edges_payload):
        raise ValueError("two radius-two edges are cube-isomorphic")
    edges_payload.sort(
        key=lambda row: int(str(row["canonical_family_hex"]), 16)
    )
    for index, row in enumerate(edges_payload):
        row["id"] = f"E{index:02d}"

    hall_mask = family_mask(hall_family())
    hall_id = vertex_id_by_mask[hall_mask]
    distances = {hall_id: 0}
    pending = deque([hall_id])
    while pending:
        current = pending.popleft()
        for neighbor in adjacency[current]:
            if neighbor not in distances:
                distances[neighbor] = distances[current] + 1
                pending.append(neighbor)
    distance_histogram = Counter(distances.values())
    if distance_histogram != {0: 1, 1: 4, 2: 22}:
        raise ValueError(
            f"unexpected Hall-distance profile {distance_histogram}"
        )

    return {
        "schema": "hall-corner-exchange-radius2-inventory-v1",
        "scope": {
            "dimension": DIMENSION,
            "exchange_radius": 2,
            "source_labels": list(SOURCE_LABELS),
            "source_profile_sha256": {
                label: by_label[label]["mathematical_profile_sha256"]
                for label in SOURCE_LABELS
            },
            "exhaustive_beyond_radius_two": False,
        },
        "hall_vertex_id": hall_id,
        "counts": {
            "cornerless_maximum_vertices": len(vertices_payload),
            "two_corner_extension_edges": len(edges_payload),
            "cube_isomorphism_types": (
                len(vertices_payload) + len(edges_payload)
            ),
            "hall_distance_histogram": {
                str(distance): count
                for distance, count in sorted(distance_histogram.items())
            },
        },
        "vertices": vertices_payload,
        "edges": edges_payload,
    }


def audit_interval(
    inventory_path: Path,
    start: int,
    stop: int,
) -> dict[str, object]:
    """Audit all elementary pc-minors of a vertex interval."""

    inventory = load_payload(inventory_path)
    if inventory["schema"] != (
        "hall-corner-exchange-radius2-inventory-v1"
    ):
        raise ValueError("unexpected inventory schema")
    vertex_rows = inventory["vertices"]
    if not 0 <= start <= stop <= len(vertex_rows):
        raise ValueError(
            f"invalid audit interval [{start},{stop}) "
            f"for {len(vertex_rows)} vertices"
        )

    audits: list[dict[str, object]] = []
    for index in range(start, stop):
        row = vertex_rows[index]
        mask = int(str(row["family_hex"]), 16)
        minor_rows: list[dict[str, object]] = []
        failures: list[dict[str, object]] = []
        for coordinate in range(DIMENSION):
            for bit in (0, 1):
                minor, minor_dimension = irredundant_projection(
                    section(mask, DIMENSION, coordinate, bit),
                    DIMENSION - 1,
                )
                dismantlable = is_dismantlable(
                    minor,
                    minor_dimension,
                )
                minor_row = {
                    "operation": "section",
                    "coordinate": coordinate,
                    "bit": bit,
                    "dimension": minor_dimension,
                    "order": minor.bit_count(),
                    "corner_count": len(corners(minor, minor_dimension)),
                    "dismantlable": dismantlable,
                }
                minor_rows.append(minor_row)
                if not dismantlable:
                    failures.append(minor_row)

            minor, minor_dimension = irredundant_projection(
                contraction(mask, DIMENSION, coordinate),
                DIMENSION - 1,
            )
            dismantlable = is_dismantlable(minor, minor_dimension)
            minor_row = {
                "operation": "contraction",
                "coordinate": coordinate,
                "bit": None,
                "dimension": minor_dimension,
                "order": minor.bit_count(),
                "corner_count": len(corners(minor, minor_dimension)),
                "dismantlable": dismantlable,
            }
            minor_rows.append(minor_row)
            if not dismantlable:
                failures.append(minor_row)

        audits.append(
            {
                "index": index,
                "vertex_id": row["id"],
                "canonical_sha256": row["canonical_sha256"],
                "elementary_minors": minor_rows,
                "failures": failures,
                "is_obstruction": not failures,
            }
        )
        print(
            f"audited {row['id']} ({index + 1}/{stop}); "
            f"failures={len(failures)}",
            flush=True,
        )
        is_dismantlable.cache_clear()
        corners.cache_clear()

    return {
        "schema": "hall-corner-exchange-radius2-audit-v1",
        "inventory_sha256": inventory["mathematical_profile_sha256"],
        "interval": {"start": start, "stop": stop},
        "audits": audits,
    }


def merge_audits(
    inventory_path: Path,
    audit_paths: list[Path],
) -> dict[str, object]:
    """Merge an exact audit cover and certify all radius-two types."""

    inventory = load_payload(inventory_path)
    audit_payloads = [load_payload(path) for path in audit_paths]
    expected_inventory_digest = inventory["mathematical_profile_sha256"]
    for payload in audit_payloads:
        if payload["inventory_sha256"] != expected_inventory_digest:
            raise ValueError("audit refers to a different inventory")

    audit_payloads.sort(
        key=lambda payload: int(payload["interval"]["start"])
    )
    cursor = 0
    audits_by_id: dict[str, dict[str, object]] = {}
    audit_digests: list[str] = []
    for payload in audit_payloads:
        interval = payload["interval"]
        start = int(interval["start"])
        stop = int(interval["stop"])
        if start != cursor:
            raise ValueError(
                f"audit cover gap/overlap at {cursor}: [{start},{stop})"
            )
        cursor = stop
        audit_digests.append(payload["mathematical_profile_sha256"])
        for row in payload["audits"]:
            vertex_id = str(row["vertex_id"])
            if vertex_id in audits_by_id:
                raise ValueError(f"duplicate audit for {vertex_id}")
            audits_by_id[vertex_id] = row
    if cursor != len(inventory["vertices"]):
        raise ValueError(
            f"audit cover stops at {cursor}, expected "
            f"{len(inventory['vertices'])}"
        )

    vertex_summaries: list[dict[str, object]] = []
    for vertex in inventory["vertices"]:
        vertex_id = str(vertex["id"])
        audit = audits_by_id[vertex_id]
        if audit["failures"] or not audit["is_obstruction"]:
            raise ValueError(f"{vertex_id}: elementary-minor audit failed")
        minor_rows = audit["elementary_minors"]
        vertex_summaries.append(
            {
                "id": vertex_id,
                "canonical_sha256": vertex["canonical_sha256"],
                "order": vertex["order"],
                "corner_count": vertex["corner_count"],
                "elementary_minors": len(minor_rows),
                "section_order_range": [
                    min(
                        int(row["order"])
                        for row in minor_rows
                        if row["operation"] == "section"
                    ),
                    max(
                        int(row["order"])
                        for row in minor_rows
                        if row["operation"] == "section"
                    ),
                ],
                "contraction_order_range": [
                    min(
                        int(row["order"])
                        for row in minor_rows
                        if row["operation"] == "contraction"
                    ),
                    max(
                        int(row["order"])
                        for row in minor_rows
                        if row["operation"] == "contraction"
                    ),
                ],
                "minor_corner_range": [
                    min(int(row["corner_count"]) for row in minor_rows),
                    max(int(row["corner_count"]) for row in minor_rows),
                ],
                "minor_failures": 0,
                "is_obstruction": True,
            }
        )

    edge_summaries = [
        {
            "id": edge["id"],
            "canonical_sha256": edge["canonical_sha256"],
            "order": edge["order"],
            "corner_count": edge["corner_count"],
            "endpoints": edge["endpoints"],
            "is_obstruction": True,
            "certificate": (
                "ample one-concept extension of an audited obstruction; "
                "apply the one-concept extension dichotomy"
            ),
        }
        for edge in inventory["edges"]
    ]

    all_digests = {
        str(row["canonical_sha256"])
        for row in vertex_summaries + edge_summaries
    }
    if len(all_digests) != len(vertex_summaries) + len(edge_summaries):
        raise ValueError("final obstruction types are not pairwise distinct")

    return {
        "schema": "hall-corner-exchange-radius2-obstructions-v1",
        "scope": {
            **inventory["scope"],
            "inventory_sha256": expected_inventory_digest,
            "audit_profile_sha256": audit_digests,
            "claim": (
                "all radius-two graph vertices and edges are distinct "
                "forbidden pc-minors for the corner-peelable ample class"
            ),
            "exhaustive_beyond_radius_two": False,
        },
        "counts": {
            **inventory["counts"],
            "audited_cornerless_obstructions": len(vertex_summaries),
            "extension_lemma_obstructions": len(edge_summaries),
            "total_distinct_obstruction_types": (
                len(vertex_summaries) + len(edge_summaries)
            ),
            "elementary_minor_audits": sum(
                int(row["elementary_minors"])
                for row in vertex_summaries
            ),
            "elementary_minor_failures": 0,
        },
        "vertices": vertex_summaries,
        "edges": edge_summaries,
    }


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="phase", required=True)

    source_parser = subparsers.add_parser("source")
    source_parser.add_argument(
        "--source",
        choices=SOURCE_LABELS,
        required=True,
    )
    source_parser.add_argument("--output", type=Path, required=True)

    discovery_parser = subparsers.add_parser("merge-discovery")
    discovery_parser.add_argument(
        "--inputs",
        type=Path,
        nargs="+",
        required=True,
    )
    discovery_parser.add_argument("--output", type=Path, required=True)

    audit_parser = subparsers.add_parser("audit")
    audit_parser.add_argument("--inventory", type=Path, required=True)
    audit_parser.add_argument("--start", type=int, required=True)
    audit_parser.add_argument("--stop", type=int, required=True)
    audit_parser.add_argument("--output", type=Path, required=True)

    merge_parser = subparsers.add_parser("merge-audits")
    merge_parser.add_argument("--inventory", type=Path, required=True)
    merge_parser.add_argument(
        "--audits",
        type=Path,
        nargs="+",
        required=True,
    )
    merge_parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    if arguments.phase == "source":
        mathematics = source_profile(arguments.source)
    elif arguments.phase == "merge-discovery":
        mathematics = merge_discovery(arguments.inputs)
    elif arguments.phase == "audit":
        mathematics = audit_interval(
            arguments.inventory,
            arguments.start,
            arguments.stop,
        )
    elif arguments.phase == "merge-audits":
        mathematics = merge_audits(
            arguments.inventory,
            arguments.audits,
        )
    else:  # pragma: no cover
        raise AssertionError(arguments.phase)
    write_payload(arguments.output, mathematics)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
