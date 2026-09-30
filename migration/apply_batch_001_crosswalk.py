#!/usr/bin/env python3
"""Apply the reviewed Batch 001 canonical class identity decisions once."""

from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path

from common import canonical_json, load_json


MIGRATION_DIR = Path(__file__).resolve().parent


# source id: (canonical short name, canonical display name)
SURVEY_IDENTITIES = {
    "acc": ("acc_graphs", "ACC graphs"),
    "almost": ("almost_median_partial_cubes", "Almost-median partial cubes"),
    "ample": ("ample", "Ample classes"),
    "ampletwo": ("ample_rank_at_most_two", "Ample classes of rank at most two"),
    "antipodal": ("antipodal_partial_cubes", "Antipodal partial cubes"),
    "aom": ("aom_tope_graphs", "AOM tope graphs"),
    "avoidance": ("isometric_word_avoidance_cubes", "Isometric word-avoidance cubes"),
    "benzenoid": ("benzenoid_systems", "Benzenoid systems"),
    "bip": ("bipartite_graphs", "Bipartite graphs"),
    "bubblesort": ("bubble_sort_graphs", "Bubble-sort graphs"),
    "c4cactoid": ("c4_cactoids", "C4-cactoids"),
    "c4tree": ("c4_trees", "C4-trees"),
    "cayley": ("cayley_partial_cubes", "Cayley partial cubes"),
    "cellular": ("cellular_partial_cubes", "Cellular partial cubes"),
    "com": ("com_tope_graphs", "COM tope graphs"),
    "condanti": ("conditional_antimatroids", "Conditional antimatroids"),
    "cornercom": ("corner_peelable_coms", "Corner-peelable COMs"),
    "crosscomplete": ("complete_crossing_graph_partial_cubes", "Complete crossing-graph partial cubes"),
    "cubicvt": ("cubic_vertex_transitive_partial_cubes", "Cubic vertex-transitive partial cubes"),
    "daisy": ("daisy_cubes", "Daisy cubes"),
    "daisyres": ("daisy_resonance_graphs", "Daisy resonance graphs"),
    "damp": ("dismantlable_ample", "Dismantlable ample classes"),
    "deletedcubes": ("vertex_deleted_cubes", "Vertex-deleted cubes"),
    "diametrical": ("diametrical_partial_cubes", "Diametrical partial cubes"),
    "distbalanced": ("distance_balanced_partial_cubes", "Distance-balanced partial cubes"),
    "edgecritical": ("edge_critical_partial_cubes", "Edge-critical partial cubes"),
    "eiso": ("sequential_isometric_partial_cubes", "Sequential-isometric partial cubes"),
    "evencycleblocks": ("edge_even_cycle_block_partial_cubes", "Edge-or-even-cycle block partial cubes"),
    "evencycles": ("even_cycles", "Even cycles"),
    "fibonacci": ("fibonacci_cubes", "Fibonacci cubes"),
    "hypercellular": ("hypercellular_partial_cubes", "Hypercellular partial cubes"),
    "hypercubes": ("hypercubes", "Hypercubes"),
    "learning": ("learning_space_graphs", "Learning-space graphs"),
    "locreal": ("locally_realizable_coms", "Locally realizable COMs"),
    "lucas": ("lucas_cubes", "Lucas cubes"),
    "maximum": ("maximum", "Maximum classes"),
    "median": ("median", "Median classes"),
    "middlelevels": ("middle_levels_graphs", "Middle-levels graphs"),
    "mirror": ("mirror_graphs", "Mirror graphs"),
    "netlike": ("netlike_partial_cubes", "Netlike partial cubes"),
    "om": ("om_tope_graphs", "OM tope graphs"),
    "opphelly": ("opposite_semicube_helly_partial_cubes", "Opposite-semicube-Helly partial cubes"),
    "oppiso": ("opposite_semicube_isomorphic_partial_cubes", "Opposite-semicube-isomorphic partial cubes"),
    "pasch": ("pasch_partial_cubes", "Pasch partial cubes"),
    "paths": ("paths", "Paths"),
    "pc": ("partial_cubes", "Partial cubes"),
    "peano": ("peano_partial_cubes", "Peano partial cubes"),
    "planar": ("planar_partial_cubes", "Planar partial cubes"),
    "polat": ("polat_partial_cubes", "Polat partial cubes"),
    "prismfree": ("prism_free_almost_median_partial_cubes", "Prism-free almost-median partial cubes"),
    "propersemi": ("proper_semi_median_partial_cubes", "Proper semi-median partial cubes"),
    "purecom": ("pure_coms", "Pure COMs"),
    "qht": ("quasi_hypertori", "Quasi-hypertori"),
    "rankr": ("partial_cubes_rank_at_most_r", "Partial cubes of rank at most r"),
    "ranktwo": ("two_dimensional_partial_cubes", "Two-dimensional partial cubes"),
    "realample": ("realizable_ample", "Realizable ample classes"),
    "realcom": ("realizable_coms", "Realizable COMs"),
    "regular": ("regular_partial_cubes", "Regular partial cubes"),
    "resonance": ("connected_resonance_graphs", "Connected resonance graphs"),
    "rper": ("cover_recursive_peripheral", "Cover-recursively peripheral partial cubes"),
    "semimedian": ("semi_median_partial_cubes", "Semi-median partial cubes"),
    "shellcom": ("shellable_coms", "Shellable COMs"),
    "simplex": ("simplex_graphs", "Simplex graphs"),
    "squaregraph": ("squaregraphs", "Squaregraphs"),
    "stars": ("stars", "Stars"),
    "taucomplete": ("complete_tau_graph_partial_cubes", "Complete tau-graph partial cubes"),
    "thetacomplete": ("complete_theta_graph_partial_cubes", "Complete Theta-graph partial cubes"),
    "tiled": ("tiled_partial_cubes", "Tiled partial cubes"),
    "treelike": ("tree_like_partial_cubes", "Tree-like partial cubes"),
    "trees": ("trees", "Trees"),
    "treezone": ("tree_zone_partial_cubes", "Tree-zone partial cubes"),
    "twoarc": ("two_arc_transitive_partial_cubes", "2-arc-transitive partial cubes"),
    "varchenko": ("arrangement_factorizable_varchenko_partial_cubes", "Arrangement-factorizable Varchenko partial cubes"),
    "vt": ("vertex_transitive_partial_cubes", "Vertex-transitive partial cubes"),
    "wellzone": ("well_embedded_zone_partial_cubes", "Well-embedded zone partial cubes"),
    "zonoreal": ("zonotopally_realizable_coms", "Zonotopally realizable COMs"),
}


MINOR_CANONICAL = {
    "minimal-nonample": ("conditioning_minimal_non_ample", "Conditioning-minimal non-ample classes"),
    "maxmin": ("maximum_generated_classes", "Maximum-generated classes"),
    "peelable-parent": ("dismantlable_maximum_parent_sections", "Sections of dismantlable maximum parents"),
    "cube-extension": ("universal_cube_extension_ample", "Ample classes with universal cube extension"),
    "both-peelable": ("both_sides_dismantlable_ample", "Ample classes with both complementary sides dismantlable"),
    "represented": ("representation_map_ample", "Ample classes with representation maps"),
    "fractionally-represented": ("fractionally_representable_ample", "Fractionally representable ample classes"),
    "learning-closure": ("convex_geometry_ideals", "Convex-geometry ideals"),
    "parent-levels": ("maximum_parent_decomposition_level_k", "Maximum-parent decomposition level k"),
    "co-peelable": ("dismantlable_complement_ample", "Ample classes with dismantlable complement"),
    "auxiliary-levels": ("maximum_section_auxiliary_level_k", "Maximum-section auxiliary level k"),
}


SAME_AS = {
    "ample": "ample",
    "maximum": "maximum",
    "median": "median",
    "daisy": "daisy_cubes",
    "learning": "learning_space_graphs",
    "treelike": "tree_like_partial_cubes",
    "peelable": "dismantlable_ample",
    "cover": "cover_recursive_peripheral",
    "realizable": "realizable_ample",
}


PARAMETERIZED_MEMBERS = {
    "d1": ("yang_decomposition_level_k", "k=1"),
    "d2": ("yang_decomposition_level_k", "k=2"),
    "m1": ("maximum_section_auxiliary_level_k", "k=1"),
    "m2": ("maximum_section_auxiliary_level_k", "k=2"),
}


EXPLICIT_PARAMETERIZED = {
    "partial_cubes_rank_at_most_r": "integer r >= 3 in the atlas entry",
    "yang_decomposition_level_k": "positive integer decomposition level k",
    "maximum_parent_decomposition_level_k": "integer parent decomposition level k >= 0",
    "maximum_section_auxiliary_level_k": "integer auxiliary-coordinate bound k >= 0, with displayed coordinates retained",
}


def main() -> int:
    crosswalk_path = MIGRATION_DIR / "class_crosswalk.json"
    output_path = MIGRATION_DIR / "canonical_classes.json"
    if output_path.exists():
        raise SystemExit(f"Refusing to overwrite completed Batch 001 output: {output_path}")
    crosswalk = load_json(crosswalk_path)
    occurrences = crosswalk["occurrences"]
    expected_survey = {record["source_id"] for record in occurrences if record["source"] == "survey_atlas"}
    expected_minor = {record["source_id"] for record in occurrences if record["source"] == "minor_inventory"}
    handled_minor = set(MINOR_CANONICAL) | set(SAME_AS) | set(PARAMETERIZED_MEMBERS)
    if expected_survey != set(SURVEY_IDENTITIES):
        raise SystemExit(f"Survey mapping mismatch: {sorted(expected_survey ^ set(SURVEY_IDENTITIES))}")
    if expected_minor != handled_minor:
        raise SystemExit(f"Minor mapping mismatch: {sorted(expected_minor ^ handled_minor)}")

    identities: dict[str, dict] = {}
    by_identity: dict[str, list[str]] = defaultdict(list)
    for occurrence in occurrences:
        source_id = occurrence["source_id"]
        if occurrence["source"] == "survey_atlas":
            short_name, name = SURVEY_IDENTITIES[source_id]
            occurrence["disposition"] = "canonical"
            occurrence["canonical_short_name"] = short_name
            occurrence["decision_basis"] = (
                "The survey atlas supplies the canonical occurrence. A primary prose characterization "
                "is still required before database promotion."
            )
            identities.setdefault(short_name, {"short_name": short_name, "name": name})
        elif source_id in MINOR_CANONICAL:
            short_name, name = MINOR_CANONICAL[source_id]
            occurrence["disposition"] = "canonical"
            occurrence["canonical_short_name"] = short_name
            occurrence["decision_basis"] = (
                "The minor atlas supplies a distinct named family; no compatible survey occurrence was found."
            )
            identities.setdefault(short_name, {"short_name": short_name, "name": name})
        elif source_id in SAME_AS:
            short_name = SAME_AS[source_id]
            occurrence["disposition"] = "same_as"
            occurrence["canonical_short_name"] = short_name
            occurrence["decision_basis"] = (
                "The two structured atlases use the same name or symbol for the same graph/set-family "
                "dictionary class. Batch 001 records identity only; characterization equivalence remains "
                "a later proof-bearing record."
            )
        else:
            short_name, parameter_value = PARAMETERIZED_MEMBERS[source_id]
            occurrence["disposition"] = "parameterized_member"
            occurrence["canonical_short_name"] = short_name
            occurrence["parameter_value"] = parameter_value
            occurrence["decision_basis"] = (
                "Represent this fixed level as a member of one parameterized catalogue class rather than "
                "allocating an unrelated class identity."
            )
            if short_name == "yang_decomposition_level_k":
                identities.setdefault(
                    short_name,
                    {"short_name": short_name, "name": "Yang-complex decomposition level k"},
                )
        by_identity[short_name].append(occurrence["occurrence_id"])

    preferred = ["partial_cubes", "bipartite_graphs"]
    ordered_names = preferred + sorted(set(identities) - set(preferred))
    canonical = []
    for allocated_id, short_name in enumerate(ordered_names, start=1):
        identity = identities[short_name]
        canonical.append(
            {
                "allocated_id": allocated_id,
                "short_name": short_name,
                "name": identity["name"],
                "source_occurrence_ids": sorted(by_identity[short_name]),
                "parameterized": short_name in EXPLICIT_PARAMETERIZED,
                "parameters": EXPLICIT_PARAMETERIZED.get(short_name),
                "identity_state": "reviewed",
                "definition_state": "awaiting_primary_characterization",
                "promotion_state": "latex",
            }
        )

    disposition_counts = Counter(record["disposition"] for record in occurrences)
    crosswalk["counts"] = {
        "occurrences": len(occurrences),
        "canonical_identities": len(canonical),
        "by_disposition": dict(sorted(disposition_counts.items())),
        "unresolved": disposition_counts.get("unresolved", 0),
    }
    crosswalk["batch"] = "001_canonical_crosswalk"
    crosswalk["identity_rule"] = (
        "Identity review does not establish a definition, relation, or theorem. Each canonical class remains "
        "blocked from data/ until a primary characterization is reviewed."
    )
    crosswalk_path.write_text(canonical_json(crosswalk), encoding="utf-8")
    output = {
        "version": 1,
        "kind": "canonical_class_id_allocation",
        "batch": "001_canonical_crosswalk",
        "count": len(canonical),
        "allocation_rule": "partial_cubes and bipartite_graphs first, followed by canonical short name; IDs are reserved and must not be reused",
        "classes": canonical,
    }
    output_path.write_text(canonical_json(output), encoding="utf-8")
    print(f"Resolved {len(occurrences)} occurrences into {len(canonical)} canonical identities.")


if __name__ == "__main__":
    raise SystemExit(main())
