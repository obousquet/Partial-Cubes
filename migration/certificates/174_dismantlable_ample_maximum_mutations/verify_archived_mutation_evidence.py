#!/usr/bin/env python3
"""Replay the finite mutation controls archived with migration batch 174."""

from __future__ import annotations

import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parent
EVIDENCE = ROOT / "ample-corners" / "evidence"


def load(path: Path):
    return json.loads(path.read_text())


def masks(n: int, k: int):
    return [sum(1 << i for i in c) for c in itertools.combinations(range(n), k)]


def submasks(mask: int):
    sub = mask
    while True:
        yield sub
        if sub == 0:
            return
        sub = (sub - 1) & mask


def shadow(family: set[int], n: int):
    return {
        support
        for support in range(1 << n)
        if len({word & support for word in family}) == 1 << support.bit_count()
    }


def maximum_vc3_corners(family: set[int], n: int):
    assert len(family) == 1 + n + n * (n - 1) // 2 + n * (n - 1) * (n - 2) // 6
    incidence = Counter()
    for support in masks(n, 3):
        outside = Counter(word & ~support for word in family)
        anchors = [anchor for anchor, count in outside.items() if count == 8]
        assert len(anchors) == 1
        anchor = anchors[0]
        incidence.update(anchor | word for word in submasks(support))
    assert set(incidence) == family
    return {word for word in family if incidence[word] == 1}


def check_maximum_exchange_formula(old: set[int], new: set[int], n: int):
    removed, = old - new
    added, = new - old
    incident_old = {i for i in range(n) if removed ^ (1 << i) in old}
    incident_new = {i for i in range(n) if added ^ (1 << i) in new}
    assert incident_old == incident_new
    old_degree = {v: sum((v ^ (1 << i)) in old for i in range(n)) for v in old}
    new_degree = {v: sum((v ^ (1 << i)) in new for i in range(n)) for v in new}
    assert min(old_degree.values()) >= 3 and min(new_degree.values()) >= 3
    lose = {
        v for v in old - {removed}
        if (v ^ removed).bit_count() == 1 and (v ^ added).bit_count() != 1
    }
    gain = {
        v for v in old - {removed}
        if (v ^ added).bit_count() == 1 and (v ^ removed).bit_count() != 1
    }
    births = {v for v in lose if old_degree[v] == 4}
    losses = {v for v in gain if old_degree[v] == 3}
    assert len(maximum_vc3_corners(new, n)) - len(maximum_vc3_corners(old, n)) == len(births) - len(losses)
    return {
        "births": sorted(births),
        "losses": sorted(losses),
        "removed_degree": len(incident_old),
    }


def check_two_corner_component():
    folder = EVIDENCE / "maximum_two_corner_component"
    checks = load(folder / "checks.json")
    manifest = load(folder / "manifest.json")
    assert hashlib.sha256((folder / "checks.json").read_bytes()).hexdigest() == manifest["checks_sha256"]
    families = [set(row) for row in checks["families"]]
    assert len(families) == 505 == len({tuple(sorted(row)) for row in families})
    triple_masks, quadruple_masks = masks(11, 3), masks(11, 4)
    corners = []
    for family in families:
        assert len(family) == 232
        assert all(len({v & support for v in family}) == 8 for support in triple_masks)
        assert all(len({v & support for v in family}) < 16 for support in quadruple_masks)
        corner_set = maximum_vc3_corners(family, 11)
        assert len(corner_set) == 2
        corners.append(corner_set)
    reachable = {0}
    changed = True
    while changed:
        changed = False
        for row in checks["processed"]:
            assert set(row["corners"]) == corners[row["id"]]
            for removed, added, count, target in row["exchanges"]:
                if target is None:
                    continue
                assert count == 2
                assert families[target] == (families[row["id"]] - {removed}) | {added}
                if row["id"] in reachable and target not in reachable:
                    reachable.add(target)
                    changed = True
    assert reachable == set(range(505))
    combined = load(folder / "combined_corner_moves.json")
    family = set(families[0])
    for row in combined["steps"]:
        removed, added = row["exchange"]
        family = (family - {removed}) | {added}
        assert maximum_vc3_corners(family, 11) == set(row["corners"])
    assert family == set(combined["family"])
    assert (1537 ^ 2023).bit_count() == 6
    return {"classes": len(families), "processed": len(checks["processed"]), "combined_steps": len(combined["steps"])}


def check_hall_reverse_mutations():
    hall = set(load(ROOT / "hall_maximum_class.json")["members"])
    assert maximum_vc3_corners(hall, 12) == set()
    report = load(EVIDENCE / "hall_reverse_mutation_criterion" / "manifest.json")
    histogram = Counter()
    for row in report["rows"]:
        neighbor = (hall - {row["removed"]}) | {row["added"]}
        corners = maximum_vc3_corners(neighbor, 12)
        assert corners == set(row["corners"])
        histogram[len(corners)] += 1
        assert check_maximum_exchange_formula(neighbor, hall, 12) == row["reverse"]
        assert len(shadow(hall | {row["added"]}, 12)) == row["union_shattered"] == 300
        assert sorted((u ^ v).bit_count() for u in corners for v in corners if u < v) == row["corner_distances"]
    assert len(report["rows"]) == report["checked"] == 40
    assert {str(k): v for k, v in histogram.items()} == report["corner_histogram"]
    special = next(row for row in report["rows"] if row["removed"] == 2 and row["added"] == 2049)
    assert special["corners"] == [1026, 2050]
    return {"exchanges": len(report["rows"]), "corner_histogram": dict(sorted(histogram.items()))}


def check_final_mutation_trees():
    hall = set(load(ROOT / "hall_maximum_class.json")["members"])
    folder = EVIDENCE / "final_mutation_tree"
    controls = load(folder / "control_manifest.json")["Hall_controls"]
    report = load(folder / "manifest.json")
    by_key = {(row["removed"], row["added"]): row for row in report["rows"]}
    cover_histogram = Counter()
    for control in controls:
        removed = control["removed"]
        order = control["coordinate_order"]
        family = {
            sum((((word ^ removed) >> old) & 1) << new for new, old in enumerate(order))
            for word in hall
        }
        groups = {kind: [int(i) for i, value in control["types"].items() if value == kind] for kind in (3, 5)}
        small, = [group for group in groups.values() if len(group) == 2]
        mask = sum(1 << i for i in small)
        fibre_words = list(submasks(mask))
        tree = {v for v in family if not v & mask and all(v | a in family for a in fibre_words)}
        active = [i for i in range(12) if i not in small]
        edges = {(v, v ^ (1 << i), i) for v in tree for i in active if v < (v ^ (1 << i)) and v ^ (1 << i) in tree}
        assert len(tree) == 11 and len(edges) == 10
        assert sorted(i for _, _, i in edges) == active
        saved = by_key[(control["removed"], control["added"])]
        assert sorted(tree) == saved["tree"] and set(map(tuple, saved["edges"])) == edges
        assert sorted(small) == sorted(saved["reduction_directions"])
        saved_leaves = {row["vertex"]: row for row in saved["leaves"]}
        leaves = 0
        for v in sorted(tree):
            directions = [i for i in active if v ^ (1 << i) in tree]
            if len(directions) != 1:
                continue
            leaves += 1
            tree_direction = directions[0]
            profiles = {
                i: sorted(a for a in fibre_words if (v | a) ^ (1 << i) in family)
                for i in active if i != tree_direction
            }
            profiles = {i: values for i, values in profiles.items() if values}
            assert all(0 < len(values) < 4 for values in profiles.values())
            assert set().union(*(set(values) for values in profiles.values())) == set(fibre_words)
            minimum = next(
                size for size in range(1, len(profiles) + 1)
                if any(set().union(*(set(profiles[i]) for i in chosen)) == set(fibre_words)
                       for chosen in itertools.combinations(profiles, size))
            )
            saved_leaf = saved_leaves[v]
            assert saved_leaf["tree_direction"] == tree_direction
            assert {int(i): values for i, values in saved_leaf["escape_profiles"].items()} == profiles
            assert saved_leaf["minimum_cover"] == minimum == 2
            cover_histogram[minimum] += 1
        assert leaves == 2
    assert len(controls) == report["controls"] == 12
    assert {str(k): v for k, v in cover_histogram.items()} == report["minimum_cover_histogram"]
    return {"controls": len(controls), "leaves": sum(cover_histogram.values()), "minimum_cover_histogram": dict(cover_histogram)}


def main():
    result = {
        "status": "PASS",
        "two_corner_component": check_two_corner_component(),
        "hall_reverse_mutations": check_hall_reverse_mutations(),
        "final_mutation_trees": check_final_mutation_trees(),
        "scope": "Exact replay of the archived finite controls. The 505-state search is partial, and the Hall tree controls do not prove existence below the recorded orders.",
    }
    output = ROOT / "archived_mutation_evidence_verification.json"
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
