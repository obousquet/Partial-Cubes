# Seven-coordinate unique-corner degree-three exclusion

All 185 proper-unit-fiber cases have fresh UNSAT solves and independent
`drat-trim` acceptance. Combined with the earlier exclusions of the other
corner degrees and the certified existence of a corner in every nonempty
ample class on at most seven coordinates, this proves that every such class
can peel to any prescribed vertex (Theorem E in `../../main.tex`).

The proper-fiber lemma selects a proper nontrivial ample Q4 fiber containing
zero at a unit neighbour of the root. Enumerating all odd family masks from
3 through 65533 by the exact shattered-support count gives 2,763 labelled
fibers. Their orbits under the 24 coordinate permutations give 185 cases.
The complete lists are in `../unique_corner_q7_degree3_fibers.json`.
No symmetry of the entire seven-coordinate family is assumed.

Each formula has 10,943 variables and 65,379 clauses. It combines Lawrence's
exact total-asymmetry criterion, exact unique-corner clauses, full support,
the proved order lower bound 33, nontrivial nonzero root-cube fibers, and 16
unit literals for the chosen fiber. Every certified case was solved afresh:
no solver assumptions or learned clauses from another case are used.

`manifest.json` records original and compressed-file SHA-256 hashes, source
fingerprints, checker identity, logs, and solve/check times. `.cnf.gz` and
`.drat.gz` contain the complete formula and proof. `sources/` contains the
exact generating scripts. Original CNF/proof hashes refer to decompressed
bytes. The exploratory report's one timeout is superseded by a checked proof
for that case; the separate unsplit 45-second test remains UNKNOWN.

To replay a stored case, decompress its files in a directory with sufficient
space, verify their SHA-256 values, and run, for example:

```bash
drat-trim /path/to/255.cnf /path/to/255.drat -t 60
```

Acceptance requires exit status zero and `s VERIFIED`. The checker used was
`/home/ec2-user/damp-proof-archive/checkers/drat-trim`; its hash is recorded.
To regenerate the complete cover from the retained source snapshots:

```bash
TMPDIR=/path/with/space uv run --with python-sat python sources/certify_unique_corner_degree3_fibers.py \
  --seconds 10 --checker /path/to/drat-trim --output /path/to/new-output
```

Apply the repository memory cap and low-priority execution rules. A timeout
in regeneration is unresolved; the retained proof can still be replayed.
The new result excludes rooted factors through seven coordinates and
cut-vertex forbidden obstructions through fifteen coordinates. It does not
improve the minimum-order interval 64--267 or finish primitive generation.
