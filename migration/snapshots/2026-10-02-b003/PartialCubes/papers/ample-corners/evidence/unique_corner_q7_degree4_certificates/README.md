# Seven-coordinate unique-corner degree-four exclusion

All 18 three-coordinate fiber orbits have fresh UNSAT solves and independent
`drat-trim` acceptance. The proper-fiber lemma in the paper shows that the 17
proper-fiber cases already cover the target; case 255 is a redundant check.
This proves the degree-four exclusion, not the degree-three case or a new
minimum-order bound.

The main proof is `thm:unique-corner-degree-four` in `../../main.tex`.
The case cover is every nontrivial ample Q3 family containing zero, modulo
coordinate permutations: 63 labelled families, 18 orbits. Selecting a root
neighbour and permuting the three fiber coordinates does not assume symmetry
of the full seven-coordinate family. A separate symbolic argument guarantees
a proper unit fiber and reduces the required cover to 62 families, 17 orbits.

Each formula combines Lawrence's exact total-asymmetry criterion, exact
unique-corner clauses, full coordinate support, the proved order lower bound
33, nontrivial nonzero root-cube fibers, and eight unit literals fixing the
chosen fiber. It has 10,943 variables and 65,378 clauses. There are no solver
assumptions or reused learned clauses in any certified case.

`manifest.json` records the original CNF/proof SHA-256 values, compressed-file
SHA-256 values, checker identity, log hashes and source fingerprints. `.cnf.gz`
and `.drat.gz` retain the complete formula/proof bytes. Source snapshots under
`sources/` preserve the exact scripts used. The exploratory incremental run is
recorded separately in `../unique_corner_q7_fibers.json`; its answer-only
statuses are superseded by these fresh certificate checks.

To replay a stored case, decompress its CNF and DRAT files into a directory
with sufficient space, then run:

```bash
drat-trim /path/to/3.cnf /path/to/3.drat -t 60
```

Check for exit status zero and `s VERIFIED`, and compare the decompressed
SHA-256 values against the manifest. The actual checker used here was
`/home/ec2-user/damp-proof-archive/checkers/drat-trim`; its hash is recorded.
To regenerate all cases, use the snapshot (with python-sat installed):

```bash
TMPDIR=/path/with/space uv run --with python-sat python sources/certify_unique_corner_fibers.py \
  --seconds 60 --checker /path/to/drat-trim --output /path/to/new-output
```

Apply the repository's memory cap and low-priority execution rules when
rerunning. Regeneration writes uncompressed files; the archived copies here
use deterministic gzip compression. A timeout in regeneration is unresolved,
not evidence against the stored checked certificate.
