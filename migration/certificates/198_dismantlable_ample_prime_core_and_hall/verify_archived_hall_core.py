#!/usr/bin/env python3
"""Verify the archived Hall first-exchange and carrier-path claims."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def load_profile(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text())
    recorded = payload.pop("mathematical_profile_sha256")
    actual = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    assert recorded == actual, (path, recorded, actual)
    payload["mathematical_profile_sha256"] = recorded
    return payload


def main() -> None:
    campaign = ROOT / "forbidden_pc_minor_campaign"
    source = load_profile(campaign / "computations/hall_corner_exchange/source_H.json")
    radius = load_profile(campaign / "hall_corner_exchange_radius2_obstructions.json")
    carriers = load_profile(campaign / "damp_hall_carrier_trees.json")

    assert source["source_order"] == 299
    assert source["source_corner_count"] == 0
    assert source["source_dismantlable"] is False
    assert source["counts"]["outside_concepts"] == 3797
    assert source["counts"]["ample_one_concept_extensions"] == 16
    assert source["counts"]["extension_profiles"] == {
        "corners=2,dismantlable=false": 4,
        "corners=4,dismantlable=true": 12,
    }

    first = {}
    for extension in source["extensions"]:
        if extension["dismantlable"]:
            continue
        added = extension["added"]
        other = next(c for c in extension["corners"] if c != added)
        swap = next(s for s in extension["swaps"] if s["removed"] == other)
        assert extension["corner_count"] == 2
        assert swap["corner_count"] == 0
        assert swap["dismantlable"] is False
        first[added] = other
    assert first == {373: 1011, 1105: 192, 1360: 3120, 1525: 252}

    counts = radius["counts"]
    assert counts["elementary_minor_audits"] == 972
    assert counts["elementary_minor_failures"] == 0
    assert counts["cornerless_maximum_vertices"] == 27
    assert counts["cube_isomorphism_types"] == 59
    assert len(radius["scope"]["source_labels"]) == 5
    assert len(set(radius["scope"]["source_profile_sha256"].values())) == 5

    profile = carriers["carrier_profile"]
    assert carriers["mathematical_profile_sha256"] == (
        "65c031beb26e738618488674026fd5308367c99c513773edf2793e09c2ba7dcf"
    )
    assert profile == {
        "all_carriers_are_paths": True,
        "edges_per_path": 10,
        "every_maximal_cube_has_opposite_facet_cover": True,
        "internal_edges_per_path": 8,
        "internal_ridge_count_per_triple_histogram": {"1": 36, "2": 60, "3": 124},
        "pair_count": 66,
        "vertices_per_path": 11,
    }
    assert len(carriers["carrier_trees"]) == 66

    report = {
        "schema": "archived-hall-core-verification-v1",
        "hall_order": 299,
        "ample_extensions": 16,
        "peelable_four_corner_extensions": 12,
        "nonpeelable_two_corner_extensions": 4,
        "first_exchanges": {str(k): v for k, v in sorted(first.items())},
        "audited_elementary_minors": counts["elementary_minor_audits"],
        "elementary_minor_failures": counts["elementary_minor_failures"],
        "distinct_first_source_types": 5,
        "carrier_paths": profile["pair_count"],
        "carrier_internal_ridge_histogram": profile[
            "internal_ridge_count_per_triple_histogram"
        ],
        "every_maximal_cube_has_opposite_facet_cover": profile[
            "every_maximal_cube_has_opposite_facet_cover"
        ],
    }
    output = ROOT / "archived_hall_core_verification.json"
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, sort_keys=True))


if __name__ == "__main__":
    main()
