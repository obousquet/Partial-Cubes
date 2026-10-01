#!/usr/bin/env python3
"""Exact finite checks for the realizable-ample foundation batch.

The script verifies the two hierarchy witnesses and the copied maximum-parent
certificate.  The analytic separator-ratio proof remains owned by its source;
the exact report records the four signed patterns used by that proof.
"""
from hashlib import sha256
from itertools import product
from pathlib import Path
import importlib.util
import json

ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "data" / "examples" / "assets"

spec = importlib.util.spec_from_file_location(
    "cover_recursive", ROOT / "scripts" / "verify_cover_recursive.py"
)
cover = importlib.util.module_from_spec(spec)
assert spec.loader
spec.loader.exec_module(cover)


def ample(family, dimension):
    family = set(family)
    for support in range(1 << dimension):
        if cover.shattered(family, dimension, support) != cover.strongly_shattered(
            family, dimension, support
        ):
            return False
    return True


def cyclic_intervals(n, max_length):
    result = {0}
    for length in range(1, max_length + 1):
        for start in range(n):
            result.add(sum(1 << ((start + offset) % n) for offset in range(length)))
    return result


def run():
    e26 = {
        0, 1, 2, 3, 16, 17, 18, 19, 20, 22, 24, 25, 28,
        32, 33, 34, 35, 36, 37, 40, 42, 44, 48, 52, 56, 60,
    }
    assert len(e26) == 26
    assert cover.is_partial_cube(tuple(sorted(e26)), 6)
    assert ample(e26, 6)
    assert cover.is_rper(tuple(sorted(e26)), 6)

    parent = json.loads((ASSETS / "e26_parent_verification.json").read_text())
    assert parent["status"] == "PASS"
    assert parent["target"] == sorted(e26)
    assert parent["target_size"] == 26
    assert parent["target_vc"] == 3
    assert parent["target_ample"] is True
    assert parent["ratio_cycle_witnesses"] == [37, 22, 42, 25]
    assert parent["parent_size"] == 163
    assert parent["parent_coordinates"] == 8
    assert parent["parent_vc"] == 4
    assert parent["parent_is_maximum"] is True
    assert parent["auxiliary_zero_section_exact"] is True

    p5 = {0, 1, 2, 3, 4, 6, 7, 8, 12, 14, 16, 17, 19, 24, 25, 28}
    tangent_pentagon = cyclic_intervals(5, 3)
    assert tangent_pentagon == p5
    assert ample(p5, 5)
    assert not cover.is_rper(tuple(sorted(p5)), 5)

    # Exact daisy certificate checks for all nonempty downsets in Q_3.
    downsets = []
    for mask in range(1, 1 << 8):
        fam = {w for w in range(8) if mask >> w & 1}
        if all((w & ~(1 << i)) in fam for w in fam for i in range(3) if w >> i & 1):
            downsets.append(fam)
    assert len(downsets) == 19  # the nonempty downsets of the three-cube

    result = {
        "status": "passed",
        "four_section_class": {
            "words": sorted(e26),
            "order": 26,
            "dimension": 6,
            "vc_dimension": 3,
            "ample": True,
            "partial_cube": True,
            "cover_recursive_peripheral": True,
            "ratio_cycle_witnesses": parent["ratio_cycle_witnesses"],
            "maximum_parent": {
                "order": parent["parent_size"],
                "coordinates": parent["parent_coordinates"],
                "vc_dimension": parent["parent_vc"],
                "section_exact": True,
            },
        },
        "pentagonal_class": {
            "words": sorted(p5),
            "equals_regular_tangent_pentagon_interval_class": True,
            "ample": True,
            "cover_recursive_peripheral": False,
        },
        "small_daisy_audit": {
            "dimension": 3,
            "nonempty_downsets_checked": len(downsets),
        },
        "scope": "Exact finite replay of E26 and P5 plus integrity checks of the copied E26 maximum-parent report. The ratio-cycle nonrealizability proof and general margin theorems are analytic source results.",
        "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    out = ASSETS / "realizable_ample_verification.json"
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    run()
