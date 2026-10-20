# Agent handoff

## Goal and current checkpoint

This project is migrating the accumulated mathematics about partial-cube
classes from `~/latex/PartialCubes` into a structured database and generated
website in `~/code/Partial-Cubes`. The migration is incomplete and should
continue in small, reviewed batches.

At database commit `a5eb779` (`main`, synchronized with `origin/main`):

- the database contains 2,125 records across 12 tables;
- 816 of 1,143 audited manuscript claims are represented;
- the main remaining source is `dismantlable-ample` (380/683 covered, 303
  uncovered), followed by `ample-corners` (29/53 covered);
- all other audited workspaces are complete at this checkpoint;
- Batch 199 is the latest ledger; the next natural claims are
  `prop:damp-geodesic-disagreement-tax`,
  `lem:damp-least-cornerless-cardinality-minimality`, and
  `prop:damp-least-cornerless-every-shift-dismantlable`.

## What each checkout contains

`~/code/Partial-Cubes` is the normalized catalogue and website repository.

- `data/` is authoritative for accepted records. Its 12 tables cover classes,
  characterizations and their equivalences, relations, operations and closure
  results, invariants and values, obstruction families and bases, examples,
  and general results. One record is one numbered JSON file; references use
  `#table/short_name`.
- `data/latex/references.bib` owns bibliography entries shared by accepted
  database records.
- `docs/` is entirely generated, but committed. GitHub Pages deploys `main`
  from `/docs`; never edit these files by hand.
- `../math_database` supplies the generic validator-independent renderer and
  local server used by the Makefile.
- `migration/` contains frozen source manifests, byte snapshots, finite
  certificates, generated audits, review queues, and one immutable ledger per
  accepted batch. `migration/generated/` is extraction output, not database
  input.
- `sync/ownership.json` is the definitive ownership contract;
  `sync/generate_latex_catalog.py` is the one-way database-to-LaTeX interface.

`~/latex/PartialCubes` is the research and scholarly-writing workspace.

- `papers/` contains the survey, minor atlas, theorem papers, active research
  papers, bibliographies, figures, and paper-local scripts.
- `forbidden_pc_minor_campaign/` contains the dismantlable-ample research
  campaign, computations, exact reports, dashboards, and source material used
  by many later batches.
- `computations/`, `reviews/`, and `scripts/` contain supporting experiments,
  reviews, and tooling.
- Full proofs, paper-specific theorem wording and notation, exposition,
  diagrams, historical discussion, open questions, conjecture programmes,
  ledgers, and research state remain LaTeX-owned even after a normalized fact
  is entered in the database.

The live LaTeX checkout is active and currently has extensive unrelated
uncommitted work. Its HEAD was `d382927c4acbc4fa2e96e8204e0528c3f5b527c5` on
2026-10-20. Do not clean, reset, reformat, or casually edit that tree, and do
not treat its current bytes as the source for an existing batch.

## Ownership and synchronization

Ownership moves **record by record**, not paper by paper:

1. Before review, the fact is LaTeX-owned.
2. A reviewed database record may shadow the source while provenance and
   parity are checked.
3. Once accepted, the normalized catalogue fact and concise database proof are
   database-owned. The paper still owns its full proof and scholarly form.

Every accepted claim must have a stable `proof_source`, normally
`PartialCubes:path/to/file.tex#label`, plus citations and a concise proof when
useful. Do not use line numbers or a dashboard/research note as final proof
provenance. Broad questions stay as prose in LaTeX.

The current dismantlable-ample migration uses the frozen manuscript

```text
migration/snapshots/2026-10-02-b002/PartialCubes/
  papers/dismantlable-ample/sections/dismantlable_ample_classification.tex
```

at source commit `7ea344598fa468ab29cf94c7220ffeb989267805`, SHA-256
`5c575ab2b35cfae192bdd218167477b8da30a6d3f7f195000a3bf98c9bd1647e`.
Use this frozen file for the next claims. If newer live text is needed, make an
explicit supplemental snapshot and record its commit, byte hash, and state;
never silently replace the reviewed source with “latest file wins.”

Do not mirror accepted facts manually back into papers. Only a specifically
chosen catalogue consumer should import a generated fragment. The supported
one-way outputs are `class-table`, `class-definitions`, and `bibliography`, for
example:

```bash
python3 sync/generate_latex_catalog.py \
  --section class-definitions \
  --output ~/latex/PartialCubes/includes/generated/class_definitions.tex

python3 sync/generate_latex_catalog.py \
  --section class-definitions \
  --output ~/latex/PartialCubes/includes/generated/class_definitions.tex \
  --check
```

Only replace a LaTeX block classified as `generated_catalogue` in
`migration/latex_consumers.json`. Preserve mixed or paper-specific blocks, and
build and compare the affected paper before promoting the handoff.

## Continuing the migration

For each batch:

1. Start from a clean, synchronized `~/code/Partial-Cubes` checkout. Read the
   complete statement and proof in the frozen source, together with nearby
   definitions and dependencies. Check `migration/manuscript_result_audit.json`
   to confirm that each label is still uncovered.
2. Create atomic records under the appropriate `data/` tables. Preserve exact
   scope and hypotheses; avoid duplicating an existing stronger record. Use the
   next unused numeric ID and a stable descriptive `short_name`.
3. If a result is computer-assisted, copy the exact clean scripts, inputs, and
   reports into a new `migration/certificates/NNN_.../` packet. Add a manifest
   with byte hashes and, where practical, a compact independent verifier.
   Large or otherwise unrecreatable artifacts must be preserved rather than
   summarized away.
4. Run:

   ```bash
   make migration-extract
   make check
   make migration-check
   git diff --check
   ```

   `make check` validates the data, rebuilds `docs/`, normalizes the generated
   site, and checks both graph views. Commit the regenerated site with the data.
   `make migration-check` verifies frozen-source recoverability and that audits
   are current. If live-source drift causes a failure, investigate it; do not
   overwrite a snapshot or regenerate from unreviewed live changes.
5. Add `migration/batches/NNN_description.json` with source refs and packets,
   hashes of each complete statement-plus-proof block, record hashes,
   validation counts, ownership/preservation notes, and explicitly deferred
   work. Verify every recorded hash before committing.
6. Commit and push the complete batch to `main`, then confirm both a clean tree
   and `git rev-list --left-right --count origin/main...HEAD` equal to `0 0`.

The class graph renders the transitive reduction of established inclusions.
Stored redundant inclusions are allowed, but rendered edges must form the
Hasse diagram. Keep strict inclusion visually distinct, preserve closure
status through proof-bearing `operation_results`, ensure every class has an
audited relation to Partial Cubes, and retain every minor-atlas subclass of
Ample in the graph.

Read `README.md`, `MIGRATION_PLAN.md`, `migration/README.md`, and
`sync/ownership.json` before changing policy or creating a new source freeze.
