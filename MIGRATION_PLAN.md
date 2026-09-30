# Partial Cubes database migration plan

## Objective

Migrate accepted mathematical content from `~/latex/PartialCubes` into the
structured database without making any one manuscript the accidental source of
truth. The migration is organized around canonical classes and atomic claims,
not around copying documents in sequence.

During migration, the LaTeX workspace remains authoritative for material that
has not passed a database review. Authority moves record by record. The papers
retain exposition and long proofs throughout.

## Current source inventory

The migration starts from several overlapping sources with different roles.

| Source | Current content | Migration role |
|---|---:|---|
| Survey class atlas | 76 mathematical nodes | Candidate class identities, graph labels, aliases, and visual categories |
| Survey containment JSON | 134 asserted containments | Candidate binary relations; the diagram alone is not proof provenance |
| Minor-atlas inventory | 24 families | Candidate classes, closure facts, bases, and assessments |
| Minor-atlas relations | 48 relations: 43 proved, 5 open | Candidate inclusions, strictness facts, open comparisons, and non-containments |
| Minor-atlas source packets | 33 packets | Source paths, hashes, labels, and theorem locations to verify |
| Minor-atlas questions | 13 questions | Remain in LaTeX; checked only to ensure they are not accidentally imported |
| Survey TeX | 324 labels, including 19 definitions and 83 theorems | Definitions, characterizations, literature-backed relations, and terminology |
| Twelve paper workspaces | Over 2,000 labels and several maturity levels | Theorem packages, proofs, examples, obstructions, and invariants |
| Root bibliography | 168 entries | Main bibliography seed, subject to key and identity reconciliation |
| Five local paper bibliographies | 70 entries before deduplication | Additional entries and possible key conflicts |

The survey and minor atlas have nine likely shared class identities. Six use
the same IDs: `ample`, `daisy`, `learning`, `maximum`, `median`, and
`treelike`. Three appear to be aliases:

| Survey ID | Minor-atlas ID | Provisional canonical name |
|---|---|---|
| `damp` | `peelable` | `dismantlable_ample` |
| `rper` | `cover` | `cover_recursive_peripheral` |
| `realample` | `realizable` | `realizable_ample` |

If these nine matches survive review, the two structured inventories describe
about 91 canonical classes. This is a planning estimate, not an ID assignment.
Graph classes, set families, tope graphs, and parameterized hierarchies must not
be merged merely because their display names are similar.

## Migration principles

1. **Extract candidates before writing records.** Extractors write migration
   queues. They never write directly to `data/`.
2. **Create a class only with a primary characterization.** A TikZ node is
   enough to propose an identity, but not enough to create a database class.
3. **Promote atomic claims.** A theorem containing an inclusion, three closure
   facts, and a strictness witness becomes separate records with shared
   provenance.
4. **Do not infer proof strength from presentation.** A solid arrow, a field
   saying `yes`, or an inventory status saying `proved` is a lead. Promotion to
   `established` requires inspection of the cited theorem or proof.
5. **Preserve scope exactly.** Finite versus infinite, simple versus
   non-simple, graph versus set-family, fixed-rank versus parameterized, and
   ambient-universe qualifications remain explicit.
6. **Do not use “latest file wins.”** Conflicting formulations enter a review
   queue. Neither overwrites the other.
7. **Use stable labels, never line numbers.** Provenance points to a repository,
   relative path, and LaTeX label or exact published theorem.
8. **Keep broad questions out of the database.** Only a precise open relation,
   characterization, or basis claim may receive an `open` record.
9. **Do not copy long proofs.** The database stores a concise proof when useful
   and otherwise points to the canonical proof.
10. **Move ownership only after parity checks.** Existing LaTeX definitions are
    not replaced by generated fragments until the corresponding database batch
    has been reviewed.
11. **Preserve the papers as scholarly artifacts.** Migration must not flatten
    a paper into database fields or discard exposition, proof structure,
    notation, diagrams, history, examples, or local context that the database
    cannot reconstruct.
12. **Replace only explicitly shared catalogue text.** Finished and
    paper-specific manuscripts may remain entirely authored LaTeX. Generated
    inputs are intended for active catalogue or survey views, not as a reason
    to rewrite every source paper.

## Repository ownership contract

The final division follows the model used by Combinatorial Parameters, adapted
to a workspace containing many papers rather than one survey. The database
owns normalized catalogue facts. The LaTeX repository owns scholarly
documents. Generated files are the one-way interface from the database to a
LaTeX consumer.

| Information | Authority after promotion | Representation and rule |
|---|---|---|
| Canonical class identity, aliases, category, scope, and catalogue summary | Database | `data/classes/*.json`; the summary is concise catalogue prose, not a replacement for a paper's introduction |
| Canonical class criteria | Database | `data/characterizations/*.json`; an active survey may import generated definitions |
| Equivalence, relation, closure, invariant-value, obstruction-basis, and general result facts | Database | Typed records determine the website and graphs; do not maintain a second editable catalogue table in LaTeX |
| Operation, invariant, obstruction-family, and example catalogue entries | Database | Their normalized definitions and metadata live in the corresponding tables |
| Concise catalogue proofs | Database | Optional `proof` fields; generate a LaTeX appendix only when a manuscript deliberately consumes these proofs |
| Canonical full proofs | LaTeX | The source theorem and proof remain in their paper; database records use `proof_source` with a stable label |
| Paper-specific theorem formulation and notation | LaTeX | Retained even when the database has a normalized statement of the same mathematical fact |
| Narrative, motivation, history, diagrams, examples in context, and document organization | LaTeX | Never reconstructed from database records |
| Broad questions, conjecture programmes, ledgers, dashboards, and research notes | LaTeX | Precise catalogue claims may be promoted, but the surrounding research state stays in the source workspace |
| Shared bibliography used by generated catalogue material | Database | `data/latex/references.bib` generates `includes/generated/references.bib`; paper-only references may stay in local bibliographies |
| Website pages and graph layouts | Database | Derived artifacts in `docs/`; the survey's editorial TikZ layout remains LaTeX-owned |

Database ownership is **record-scoped during migration**. Content already
accepted into `data/` is database-owned; all other mathematical content remains
LaTeX-owned. A source theorem does not cease to belong to its paper when a
normalized database record cites it.

The database repository should eventually provide these shared generated
artifacts:

| Generated artifact | Database source | Generator section | Intended consumer |
|---|---|---|---|
| `includes/generated/class_table.tex` | `classes` plus closure summaries derived from `operation_results` | `class-table` | Survey/catalogue overview |
| `includes/generated/class_definitions.tex` | `classes` and `characterizations` | `class-definitions` | Survey/catalogue definition section |
| `includes/generated/references.bib` | `data/latex/references.bib` | `bibliography` | Any manuscript consuming generated catalogue text |

Additional operation, invariant, obstruction, or concise-proof fragments
should be added only when a named manuscript will actually consume them. A
database table does not need a LaTeX mirror merely because it exists.

As in Combinatorial Parameters, every generated file must name its generator,
say that it must not be edited, and be reproducible with `--check`. The database
owns the generator and source data. The LaTeX repository commits the generated
output so that each manuscript remains buildable on its own.

### Preserving non-reconstructible LaTeX

Before changing any LaTeX consumer, classify each affected source environment
or contiguous block as:

- `generated_catalogue`: all semantic content is intentionally database-owned
  and the block may be replaced by a generated input;
- `latex_only`: the block contains document-specific content and must remain;
- `mixed`: split the block so only its catalogue fragment is generated;
- `archival`: a stable or published artifact that will not be rewritten.

Record this classification in `migration/latex_consumers.json`, including the
paper, source path, labels, database records, intended action, and reviewed
source hash. A paper is not required to consume generated files. In particular,
published or stable theorem papers should normally remain unchanged and serve
as proof sources.

For every handoff batch:

1. Pin the exact pre-handoff LaTeX commit and file hashes and create a durable
   migration tag or equivalent archival reference in the LaTeX repository.
2. Build and retain the pre-handoff PDF or its checksum and normalized text
   extraction as a comparison artifact.
3. Replace only blocks marked `generated_catalogue`; split `mixed` blocks and
   retain all LaTeX-owned sentences, labels, diagrams, and proof environments.
4. Rebuild the document and compare labels, citations, page inventory,
   normalized PDF text, and selected rendered pages. Every removed passage must
   be either byte-preserved in the pinned source or deliberately represented by
   an accepted database record; any semantic difference requires review.
5. Record the database commit, LaTeX commit, generated-file hashes, builds, and
   reviewer decision in the batch ledger.

Git history is the recovery mechanism for exact old source, but it is not a
license to remove useful material from the live manuscripts. Narrative and
proof content stays in the working LaTeX documents. Generated inputs should be
small, local, and auditable; no paper is regenerated wholesale from the
database.

### Two-repository handoff

Each promoted batch moves through three states:

1. `latex`: the source workspace is authoritative; extracted database
   candidates have no public ownership.
2. `shadowed`: reviewed database records exist and generated output is compared
   with the LaTeX consumer, but the authored LaTeX block is still authoritative.
3. `database`: the relevant shared catalogue block imports checked generated
   output; the database record is authoritative for that normalized fact.

The batch ledger records the paired commits from both repositories. Promotion
to `database` occurs only after the database checks pass, every affected LaTeX
document builds, and preservation checks pass. Rollback restores the authored
block from the pinned source and marks the record `shadowed`; it never requires
reconstructing paper prose from database fields.

## Exact source snapshots

The LaTeX working tree is active and currently contains tracked and untracked
research changes. A Git commit alone is therefore insufficient to identify the
source used by a migration batch.

Before extraction, create `migration/source_manifest.json` containing:

- the Git repository name and `HEAD` commit;
- each source path and SHA-256 digest;
- whether the file is tracked, modified, or untracked;
- the labels used from that file;
- referenced external repositories and their commits;
- the extraction date and extractor version.

A modified or untracked manuscript may be inventoried, but its claims are
blocked from promotion until that exact text is reviewed and deliberately
accepted. Each migration batch records the manifest digest it used. Later
source edits then produce drift reports rather than silently changing database
meaning.

For recoverability, a clean tracked file is pinned by repository commit and
path. Every selected modified or untracked file is also copied byte-for-byte
into the migration snapshot before extraction. Extractors read the frozen
bytes, so an active research process may continue editing the live workspace
without changing an in-progress batch. Direct proof-source selection excludes
`research/` and `archive/` trees unless a later review deliberately promotes a
specific artifact.

Use portable locators of these forms:

```text
PartialCubes:papers/daisy-cubes/main.tex#thm:intrinsic
YangDim:minimal_nonample_classes.tex#thm:classification
CompressionSchemes:partitionable_yang_complexes.tex#thm:maxmin
published:CCMW#Theorem-6.8
```

The manifest is a migration aid, not a new public `sources` table. Final claim
records continue to use `references`, `proof_source`, and `proof`.

## Migration workspace

Create a `migration/` directory with reviewable, generated artifacts:

| File | Purpose |
|---|---|
| `source_manifest.json` | Exact input files, repository commits, hashes, and labels |
| `class_crosswalk.json` | Every source class occurrence and its proposed canonical class |
| `claim_queue.json` | Atomic candidate claims extracted from structured and TeX sources |
| `conflicts.json` | Incompatible names, scopes, statuses, endpoints, or statements |
| `bibliography_queue.json` | Citation keys, normalized identities, collisions, and dispositions |
| `coverage.json` | Candidate, promoted, deferred, rejected, and unresolved counts by source |
| `latex_consumers.json` | Block-level ownership, preservation action, and generated-input mapping for each manuscript |
| `batches/NNN_*.json` | Exact record list, sources, reviewer decisions, and validation for each batch |

These files must distinguish generated extraction from human decisions. A
regeneration may update source occurrences, but it must not overwrite a manual
canonicalization or disposition.

Initial read-only tools should be:

- `migration/extract_survey_atlas.py` for the 76 nodes and 134 containments;
- `migration/extract_minor_inventory.py` for families, operations, relations,
  source packets, and questions;
- `migration/index_latex_labels.py` for definitions and theorem environments;
- `migration/check_source_manifest.py` for source drift;
- `migration/audit_coverage.py` for reconciliation totals;
- `migration/audit_bibliography.py` for missing keys and identity collisions.

## Class crosswalk rules

Every class occurrence receives one of these dispositions:

- `canonical`: introduces a new canonical class;
- `same_as`: an exact synonym or an equivalent model already represented by a
  canonical class;
- `parameterized_member`: one member of a parameterized class record;
- `related_not_equal`: similar terminology but mathematically distinct;
- `presentation_only`: legend, heading, or diagram helper;
- `unresolved`: requires mathematical review;
- `excluded`: outside the intended catalogue, with a reason.

For every proposed merge, record:

- the two source definitions;
- ambient scopes and parameter conventions;
- the theorem or dictionary proving equivalence, when the models differ;
- the chosen canonical `short_name`;
- retained aliases and graph labels.

Numeric IDs are assigned only after the crosswalk decision. IDs are never
renumbered or reused. The stable public identifier is the `short_name`; numeric
IDs provide ordering and file naming.

## Status and admission rules

The database status describes the mathematical claim, while the migration
queue separately describes review progress.

| Database status | Admission rule |
|---|---|
| `established` | Exact statement and scope checked against a proof, a precise theorem source, or a published result |
| `open` | A precise, currently open mathematical statement with reviewed wording |
| `conjectured` | Explicitly proposed as a conjecture, not merely absent from the literature |
| `needs_verification` | Intentionally public claim whose proof or scope still needs checking |
| `refuted` | A useful historical or comparison claim with an explicit counterexample or correction |

Raw migration candidates normally stay outside `data/`. `needs_verification`
is not a bulk-import status. It is used only when displaying the unresolved
claim is itself useful and its uncertainty is accurately stated.

Manuscript review status is handled as follows:

| Source group | Default migration treatment |
|---|---|
| Survey atlas and minor inventory | Discovery and crosswalk sources; inspect their cited proofs before promotion |
| Daisy-cubes and gluing papers | First theorem packages to migrate because their current source bundles and focused internal reviews are stable |
| Learning-spaces paper | Migrate exact labeled classical equivalences and proved operation results; keep incomplete obstruction programmes out |
| VC–Littlestone note | Use as an index to primary sources and verified examples; avoid turning exposition into new provenance |
| Cover-recursive paper | Migrate only labels covered by a stable review, or wait for the documented delta defects and source refresh |
| Structural-invariants and realizability-margin papers | Migrate definitions first; promote laws and class characterizations in reviewed theorem batches |
| Coefficient-topology paper | Defer broad import; select only results needed by accepted class or invariant records |
| Ample-corners paper | Import bounds and examples only with their exact certificate or proof disposition |
| Dismantlable-ample compilation | Never migrate wholesale; use only focused, closed theorem packages extracted and reviewed elsewhere |
| Dashboards, ledgers, and `research/` files | Never direct proof sources; they may point to a theorem or a replayable certificate |

## Reconciliation rules

### Duplicate claims

Normalize a candidate fingerprint from its table, mathematical type, canonical
endpoints, ambient scope, hypotheses, and parameter range.

- Identical claims from several documents become one record with one canonical
  proof source and all useful published references.
- Equivalent statements in different mathematical languages become separate
  characterization records joined by an equivalence theorem.
- A stronger theorem and a weaker corollary are not automatically duplicate.
  Store the stronger direct fact unless the weaker form has independent value.
- Derived transitive inclusions are not stored merely to reproduce a diagram.

### Conflicts

Conflicts are never resolved by source date alone. The review packet shows both
statements, their source hashes, proof locations, and consequences. Possible
dispositions are:

- one statement has a narrower scope;
- one source is stale and should be corrected;
- both are valid and require separate records;
- one claim is refuted;
- the issue remains unresolved and neither claim is promoted.

### Bibliography

Import the citation closure of accepted records rather than copying every
bibliography immediately. Resolve entries by DOI, arXiv identifier, or
normalized title and authors. Preserve an existing key when it is unambiguous.
If one key names different works, block the batch until a deliberate rename and
consumer update are recorded.

The database bibliography becomes authoritative only after all generated
LaTeX consumers reproduce their previous citation keys and builds.

## Mapping source material to tables

| Source statement | Database representation |
|---|---|
| Named mathematical family | `classes` plus a primary `characterizations` record |
| Alternative membership criterion | `characterizations` |
| “The following conditions are equivalent” | One `characterization_equivalences` record joining all conditions |
| Inclusion, equality, or non-containment | `relations` |
| Strictness or failure witness | `examples`, referenced by the claim |
| Definition of projection, conditioning, product, expansion, or gluing | `operations` |
| Closure, failed closure, or exact closure | `operation_results`; class closure booleans are derived caches |
| Definition of a numerical or structural quantity | `invariants` |
| Formula, bound, range, or value | `invariant_values` with an explicit quantifier |
| Invariant threshold with a converse theorem | An invariant value, an invariant-condition characterization, and an equivalence record |
| Reusable forbidden family | `obstruction_families` |
| Complete or partial avoidance theorem | `obstruction_bases` |
| Avoidance as a class criterion | A forbidden-minor characterization linked to the basis |
| Recognition, construction, or implication not covered above | `results` |
| Broad programme or informal question | Remains in LaTeX |

## Migration phases

### Phase 0: freeze inputs and install audits

1. Create the migration workspace and source manifest.
2. Record the current LaTeX `HEAD`, dirty paths, and exact hashes.
3. Add extractors for the two structured atlases and a LaTeX label indexer.
4. Add coverage, drift, and bibliography audits.
5. Amend `sync/ownership.json` to express transitional authority: a database
   record owns a fact only after its batch is accepted.
6. Inventory every LaTeX document as `archival`, `latex_only`, or a potential
   generated consumer; do not alter any manuscript in this phase.

**Gate:** extractors are deterministic, do not change `data/`, and reproduce
76 survey nodes, 134 survey containments, 24 minor families, 48 minor
relations, 33 source packets, and 13 excluded questions.

### Phase 1: canonical class crosswalk

1. Reconcile all 100 class occurrences from the survey and minor atlas.
2. Review the nine likely overlaps listed above.
3. Resolve graph/set-family/tope-graph dictionary equivalences explicitly.
4. Decide parameterized records for rank classes and the `D_k`, `P_k`, and
   `M_k` hierarchies.
5. Allocate stable short names and IDs after review.
6. Build the bibliography queue for definitions used in the pilot.

**Gate:** every occurrence has a disposition; no unresolved merge creates a
database record; the provisional canonical count is explained by the
crosswalk rather than inferred from names.

### Phase 2: end-to-end pilot

Populate a small vertical slice that exercises every table before a broad
import. Recommended anchor classes are:

- partial cubes;
- ample classes;
- dismantlable ample classes;
- cover-recursive peripheral classes;
- daisy cubes;
- learning-space graphs;
- realizable ample classes;
- maximum classes and MaxMin.

The pilot should include:

- projection, conditioning, and pc-minor operation definitions;
- one primary characterization for each class;
- the standard partial-cube equivalence theorem;
- reviewed inclusions and at least one strictness witness;
- p-, c-, and pc-closure records and derived class summaries;
- one complete obstruction family and basis from the daisy paper;
- one structural invariant and value theorem;
- one named example and one general result.

**Gate:** all twelve tables contain a reviewed record, all validation passes,
both graphs render, provenance links resolve, and a reviewer can follow every
record back to the exact source statement.

### Phase 3: broad class dictionary

1. Migrate the remaining survey identities in small thematic batches:
   metric/cubical, expansion, set-family, sign-vector/COM, symmetry, and named
   families.
2. Take canonical definitions from survey prose or primary literature, not
   from the shortened TikZ labels.
3. Add aliases, graph labels, categories, ambient scopes, and parameter
   conventions.
4. Add one primary characterization with each class in the same commit.

**Gate:** every displayed class node resolves to a canonical class or an
explicit exclusion; no class lacks a primary characterization.

### Phase 4: survey hierarchy

1. Process the 134 containment candidates.
2. Locate a proof or exact citation for each inclusion.
3. Determine whether the inclusion is strict and attach a reviewed witness
   where one is claimed.
4. Store only asserted direct facts. Recompute the transitive reduction from
   the database and compare it with the survey atlas.
5. Put missing provenance, questionable scope, and redundant edges into the
   coverage report instead of guessing.

**Gate:** every original edge is promoted, derived, deferred, or rejected with
a reason; the generated hierarchy has no unexplained node or edge difference.

### Phase 5: operation and minor atlas

1. Reconcile all 24 family rows with canonical classes.
2. Convert each `p`, `c`, and `pc` cell into a proof-bearing operation result.
   This gives 72 candidate closure classifications.
3. Create `not_closed_under` claims only with a witness or exact theorem.
4. Normalize `c_closure` and `pc_closure` prose into `closure_equals` only when
   exact equality is really proved; otherwise use a weaker typed result.
5. Reconcile all 48 relations, including strictness, negative relations, and
   the five precise open comparisons.
6. Use the 33 source packets to locate proofs, but verify their recorded hashes
   and labels against the source manifest.
7. Confirm that all 13 broad questions remain outside the database.

**Gate:** every family field and relation has a disposition, class closure
summaries agree with operation results, and the minor/closure graph matches the
reviewed portion of the original map.

### Phase 6: characterizations and equivalence theorems

Process one canonical class at a time.

1. Split each multi-condition theorem into atomic characterizations.
2. Add one equivalence record for the theorem rather than pairwise duplicates.
3. Keep hypotheses common to the theorem on the equivalence record.
4. Verify that every established alternative connects to the primary
   characterization through established equivalence records.
5. Treat one-way implications as `results`, not characterizations.

Start with the survey's central dictionaries, then the daisy,
cover-recursive, and learning-space packages.

**Gate:** the characterization connectivity validator passes, theorem labels
resolve, and no criterion is promoted from a merely necessary condition.

### Phase 7: obstruction families and bases

Migrate complete, partial, and open basis statements from reviewed theorem
packages.

1. Separate the reusable obstruction family from the target class.
2. Record the ambient class and exact minor operation.
3. Distinguish complete, partial, necessary, and sufficient claims.
4. Record pc-, p-, c-, or convex-minimality only when proved.
5. Link a complete avoidance criterion to the class through a
   characterization and equivalence theorem.
6. Keep “basis not inventoried” and “no basis known” as absence of a record,
   not as mathematical claims.

Recommended order: daisy cubes, ample classes, cover-recursive classes, then
the proved learning-space subfamilies. Defer the open dismantlable-ample
generator programme.

**Gate:** each basis has compatible target, ambient, operation, families,
status, and proof provenance; completeness is never inferred from examples.

### Phase 8: invariants and values

Start with definitions whose normalization and domain are stable:

- ordering width;
- majority-completion defect;
- resolution slope and rounded slope defect;
- separation margin;
- selected decomposition or parent-depth quantities after their conventions
  are fixed.

For each invariant:

1. Migrate its definition, domain, codomain, parameters, and normalization.
2. Record exact values, formulas, bounds, or ranges with explicit quantifiers.
3. Attach named-object values to examples rather than class families.
4. Store product, deletion, gluing, and minor laws as typed results where they
   are not simply values.
5. Create a characterization only after checking the converse theorem.

**Gate:** values have exactly one subject, citation keys resolve, and no
universal class claim is confused with an individual example.

### Phase 9: remaining examples and results

Migrate only examples that are reproducibly defined and actually support a
catalogue claim. Preserve a certificate or script pointer in `verification`
when the fact is computational.

Use `results` last, for important recognition, construction, implication, and
invariant-law theorems that do not fit another table. If many records require a
new result type, review the schema before adding it.

**Gate:** every witness is referenced by a claim or documented as a standard
example; `results` has not become a prose dumping ground.

### Phase 10: LaTeX ownership handoff

This phase applies only to a manuscript deliberately selected as an active
catalogue consumer. Source papers marked `archival` or `latex_only` remain
unchanged and continue to supply full-proof provenance.

For each accepted class batch and each selected consumer:

1. Complete `migration/latex_consumers.json` and
   `sync/latex_mapping.json` with source blocks, existing labels, and deliberate
   overrides.
2. Pin the pre-handoff source and build its preservation comparison artifacts.
3. Generate class definitions and tables into `includes/generated/` in the
   LaTeX repository.
4. Compare generated statements and labels against authored source and split
   every `mixed` block so its LaTeX-only content stays in place.
5. Build every consumer before and after changing an input.
6. Replace only blocks classified and reviewed as `generated_catalogue`.
7. Generate the database bibliography and confirm citation-key parity without
   removing paper-local entries that are not part of the shared catalogue.
8. Compare normalized PDF text, labels, citations, diagrams, and selected page
   renders; account for every removed source passage.
9. Record paired database and LaTeX commits and promote the batch ownership
   state from `shadowed` to `database`.

Narrative, long proofs, and open-question sections remain authored LaTeX.

**Gate:** generated files are reproducible, `--check` detects no drift, all
consumers build, no shared catalogue definition remains independently editable,
and no LaTeX-only content or rendered scholarly artifact has been lost.

## Batch size and review packet

Use small commits: roughly 5–15 canonical classes or 10–30 claim records,
grouped by one mathematical theorem package. Do not combine broad extraction,
canonicalization, and promotion in one commit.

Each accepted batch records:

- source-manifest digest;
- candidate IDs and final record paths;
- canonicalization decisions;
- proof and bibliography checks;
- conflicts and deferrals;
- validator and generated-site results;
- graph or LaTeX parity results when applicable;
- pre/post LaTeX preservation comparison and any intentional rendered change;
- paired repository commits and the resulting ownership state;
- reviewer and date.

Required commands are:

```bash
python3 migration/check_source_manifest.py
python3 migration/audit_coverage.py
make validate
make check
```

When a batch changes database-owned LaTeX output, also run the relevant
generator with `--check` and build its consumers.

## Completion criteria

The initial migration is complete when:

- all 76 survey nodes and 24 minor-atlas families have reviewed crosswalk
  dispositions;
- all 134 survey containments and 48 minor-atlas relations are promoted,
  derived, deferred, or rejected with reasons;
- all 72 p/c/pc classifications have dispositions and every displayed closure
  summary is derived from an operation result;
- all 33 source packets resolve to pinned, checked locations;
- all 13 broad questions are confirmed to remain in LaTeX;
- every database class has a primary characterization;
- every established alternative characterization has an equivalence path;
- every established claim has proof or exact provenance;
- every citation key exists in the database bibliography;
- both graph views are generated entirely from database records;
- migration coverage contains no unexplained omissions;
- `make check` passes and the committed `docs/` tree is current;
- the accepted LaTeX consumers use checked generated definitions and
  bibliography files;
- every source paper has an ownership classification and pinned source
  snapshot, and preservation checks account for all modified LaTeX blocks.

The active research manuscripts can continue evolving after this milestone.
New facts then enter through the same candidate, review, and promotion process
instead of another bulk migration.
