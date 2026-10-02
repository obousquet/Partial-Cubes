# All ample subfamilies of one 96-word envelope peel

Source: the labelled order267 obstruction at coordinates (0,1), flipping bit1.
Its two-row common intersection C has48 words. The potential family is
C times Q1, with the fresh bit in position10. All96 memberships are free.
The exact query requires nonempty, ample and cornerless; it has no projection,
order, VC or support restriction. UNSAT excludes every cornerless ample
subfamily, hence every ample subfamily is dismantlable by corner induction.

The formula uses the previously audited shattered/strong-shattered and
missing-incident-cube encoder. `report.json` records source, CNF and rawproof
SHA256 hashes. `search.drat.gz` is the compressed raw proof. `check.log` and
`replay.log` both contain `s VERIFIED`. `verification.json` records the fresh
reconstruction and replay, recovering C directly from coordinate-zero edges.

Reproduce from repository root, under the prescribed memory/CPU limits:

    uv run --with python-sat python papers/ample-corners/scripts/explore_two_row_overlap_envelope.py
    uv run --with python-sat python papers/ample-corners/scripts/verify_two_row_overlap_envelope.py

Scope: this particular labelled envelope. No universal96-vertex lower bound,
no smaller obstruction and no resolution of the general two-row inverse.
Main source: ample-corners `prop:record-two-row-overlap-envelope`.
