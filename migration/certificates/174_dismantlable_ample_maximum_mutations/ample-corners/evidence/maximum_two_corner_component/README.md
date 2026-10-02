# Bounded two-corner mutation component

Starting from Hall's projection at coordinate0, inspect every maximum VC3
single-word exchange of each processed state; enqueue states with at most
two corners. The128-processed-state cap was reached, with505 discovered
states and377 still pending. Thus status is PARTIAL, not a closed-component
certificate. All505 families are independently checked maximum and have
exactly two corners; five change the original corner pair. Two three-step
paths, saved in corner_change_paths.json, move the two corners separately.
Their concatenation, checked in combined_corner_moves.json, still has two corners. No zero/one-corner class was found.

checks.json stores every discovered family, processed transitions, and the
pending queue. manifest.json stores fingerprints and counts; audit.json
records the independent trace/cube/reachability check. Source snapshots
are retained. No larger search is justified by these finite observations
alone; the next question is structural corner-set change under mutations.
