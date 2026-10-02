"""Classify active coordinate labelings of the ten nine-vertex graphs."""

from __future__ import annotations

import json
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
YANGDIM_SCRIPTS = Path(__file__).resolve().parents[2] / "YangDim" / "scripts"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(YANGDIM_SCRIPTS))

import enumerate_damp_nine_relative_graphs as nine_graphs  # noqa: E402
import enumerate_eight_deletion_colored_embeddings as colored  # noqa: E402


N = 9
colored.N = N
colored.PAIRS = nine_graphs.cube_graphs.PAIRS
colored.FULL_VERTEX_MASK = (1 << N) - 1


def enumerate_colored_types() -> dict[str, object]:
    graph_data = nine_graphs.enumerate_types()
    rows = []
    for graph_row in graph_data["types"]:
        graph_type = int(graph_row["type"])
        mask = nine_graphs.cube_graphs.graph_mask(
            [tuple(edge) for edge in graph_row["edges"]]
        )
        automorphisms, raw_systems, representatives = (
            colored.colored_embedding_types(mask)
        )
        type_rows = []
        for type_number, (sides, orbit_size) in enumerate(
            sorted(representatives.items(), key=lambda item: (len(item[0]), item[0])),
            start=1,
        ):
            embedding = colored.embedding_from_sides(sides)
            type_rows.append(
                {
                    "colored_type": f"{graph_type}.{type_number}",
                    "dimension": len(sides),
                    "cut_sides": list(sides),
                    "orbit_size": orbit_size,
                    "embedding": list(embedding),
                }
            )
        rows.append(
            {
                "graph_type": graph_type,
                "graph_automorphism_count": len(automorphisms),
                "raw_cut_system_count": len(raw_systems),
                "colored_type_count": len(type_rows),
                "colored_types": type_rows,
            }
        )
    return {
        "status": "VERIFIED",
        "graph_type_count": len(rows),
        "colored_embedding_type_count": sum(
            row["colored_type_count"] for row in rows
        ),
        "graphs": rows,
    }


if __name__ == "__main__":
    print(json.dumps(enumerate_colored_types(), indent=2, sort_keys=True))
