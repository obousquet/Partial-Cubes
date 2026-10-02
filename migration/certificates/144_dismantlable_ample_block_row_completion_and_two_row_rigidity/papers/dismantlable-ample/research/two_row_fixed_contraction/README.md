# Six fixed-contraction rigidity certificates

For each pair (0,1), (2,3), (4,5), (6,7), (8,9), (10,11) of the labelled267
seed, retain its two large rows D0,D1 and allow every diagonal assignment
B0,B1 with union C=D0 intersection D1. There are exactly two cornerless
ample assignments: the original rows and their exchange. Both have267 words.
This has no output-size cutoff. It does not classify cornered negative
assemblies or assemblies with changed D0,D1, and is not a global lower bound.

The symbolic coverage proof is Proposition `prop:record-fixed-contraction-rigidity`
in `papers/ample-corners/main.tex`. The formula encodes ampleness by equality
of shattered and strongly shattered supports, and cornerlessness by missing
incident-cube vertices. No necessary-only degree test is substituted.

`report.json` records six prescribed267 baseline models and six checked
at-most266 refutations. `rigidity.json` records complete membership-model
enumeration with two models per pair. Each `*_rigidity.cnf` contains the
base formula plus exactly two blocking clauses. Its DRAT proof excludes
every other assignment. `verification.json` hashes the six rebuilt formulas,
proofs, and checker binary. Check logs are retained. The computational proof
uses python-sat's Glucose4 and the existing drat-trim checker at
`/home/ec2-user/damp-proof-archive/checkers/drat-trim`.

From the repository root, replay with:

```bash
mem_cap=$(awk '/MemAvailable:/ {printf "%.0f\n", $2 * 1024 * 0.75}' /proc/meminfo)
TMPDIR=/home/ec2-user/damp-validation nice -n 10 prlimit --as="$mem_cap" --cpu=120 -- \
  uv run --with python-sat python papers/dismantlable-ample/research/verify_two_row_fixed_contraction.py
```

Exploration scripts `two_row_fixed_contraction_sat.py` and
`two_row_fixed_rigidity.py` recreate all artifacts. Their SAT timeouts and
64-model cap return UNKNOWN or MODEL-CAP, not an exhaustiveness claim;
neither limit was reached in these certificates. Clause normalization only
removes duplicate literals and tautologies. Construction order is sorted
so replay is independent of set insertion order.
