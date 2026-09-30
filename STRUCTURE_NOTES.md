# Structural notes for a Partial Cubes database

These notes record the current repository architecture and the source material
that should inform the database schema. They deliberately stop short of fixing
the schema.

## Repositories and roles

### `~/code/Partial-Cubes`

- This is the intended database repository.
- It currently has no commits and no project files; only Git and environment
  metadata are present.
- Its remote is `git@github.com:obousquet/Partial-Cubes`.

### `~/code/Combinatorial-Parameters`

- This is the closest model for the new repository.
- `data/` is the authoritative structured database.
- `docs/` is generated static HTML committed for GitHub Pages.
- `sync/` contains cross-repository generators, mappings, validation, audits,
  and mathematical verification scripts.
- The current catalogue has 61 class records, 108 parameter records, 513
  relationship records, and 1,147 value records.
- Every logical row is stored as one numbered JSON file. Each table directory
  has a `schema.json`; records use integer `id` values and stable
  `short_name` values. Cross-references use strings such as
  `#parameters/vc_dimension`.
- `data/main.json` supplies site metadata, bibliography configuration, and a
  list of graph hooks. `data/make_graph.py` converts the records into the node,
  edge, legend, cluster, and layout dictionaries expected by the shared graph
  renderer.
- `data/latex/references.bib` is a database-owned bibliography source. The
  compatibility symlink `data/bibliography.bib` points to it.

The four current table shapes are:

- `classes`: identity, name, symbol, full definition, comments, graph label,
  and computed reverse references.
- `parameters`: identity and definitions; category; symmetry and several
  monotonicity classifications; evidence; computed reverse references.
- `relationships`: two parameter references, a typed mathematical relation,
  constants/variants, epistemic status, proof or proof pointer, references,
  and one- or two-direction separation witnesses.
- `values`: one class and one parameter, exact/asymptotic value, epistemic
  status, proof, and references.

The reference repository treats statuses and provenance as mathematical data,
not editorial decoration. Relationship records distinguish established,
conjectured, open, needs-verification, and refuted claims. Witness strength and
the precise kind of comparison are also explicit.

### `~/latex/CombinatorialParameters`

- This is a separate Git repository containing the survey narrative, long
  proofs, layout, research ledgers, and exploratory work.
- It consumes six committed generated assets under `includes/generated/`:
  parameter and value tables, parameter and class definitions, value proofs,
  and the bibliography.
- `includes/defs.tex` and `includes/cl_defs.tex` are thin inputs of generated
  files rather than independent definition sources.
- `sync/ownership.json` in the database repository defines the ownership
  boundary. Structured definitions, relationships, values, concise value
  proofs, and bibliography belong to the database. Survey narrative and
  canonical long relationship proofs belong to LaTeX.
- `sync/latex_mapping.json` maps stable database short names to stable LaTeX
  labels, including overrides and grouped definitions.
- `sync/generate_latex_catalog.py` writes a selected generated section to an
  explicit output path in the LaTeX checkout. `--check` detects drift. There is
  no Git submodule or live data symlink joining the two repositories; the link
  is an explicit, reviewable generation workflow.

### `~/code/math_database`

- This is the reusable editor and static-site generator. It is a separate Git
  repository rather than vendored code.
- `generate_website.py --data_dir ... --output_dir ...` discovers every
  subdirectory containing `schema.json`, emits table indexes and record pages,
  copies shared assets, renders a bibliography, and invokes configured graph
  hooks.
- `server.py --data-dir ...` provides a local Flask editor. Saving a record
  creates a filename of the form `NNN_short_name.json`.
- Supported editable field types are `string`, `integer`, `boolean`, `text`,
  `latex`, `enum`, and `array`; `reference` fields are computed reverse
  references and are omitted from forms.
- A `reference` column names a foreign table and foreign column; its displayed
  value is populated by finding records whose referenced field points back to
  the current record.
- Arbitrary inline links are resolved from `#table/short_name`, `#table/id`, or
  an unambiguous short reference. Bibliography links use `#bib/key`, and TeX
  prose can use `\cite{key}`.
- A table may provide `render_<table>.py`, but at present this only customizes
  card titles. Most rendering is schema-driven.
- The renderer contains some Combinatorial-Parameters-specific behavior for
  the exact table names `parameters`, `relationships`, and `values`. A Partial
  Cubes schema can use the generic behavior immediately, but rich class
  relation cards may require a small generalization in `math_database` or a
  new hook rather than copying the parameter-specific special cases.
- A graph hook is fully repository-specific. This is the natural place to
  compute a transitive reduction, ranks, overlays for open or negative
  relations, styling, and popup references.
- The shared tool loads JSON but is not itself a strict JSON-schema validator.
  Project-level validation and semantic audits belong in the data repository.

## Existing source material in `~/latex/PartialCubes`

The LaTeX repository is already a multi-paper research workspace. It has twelve
paper directories:

1. `ample-corners`
2. `coefficient-topology`
3. `cover-recursive`
4. `daisy-cubes`
5. `dismantlable-ample`
6. `gluing`
7. `learning-spaces`
8. `minor-atlas`
9. `realizability-margin`
10. `structural-invariants`
11. `survey`
12. `vc-littlestone`

The repository also has active ledgers, dashboards, verification scripts,
certificates, computational campaigns, and many compatibility symlinks. It is
not suitable as the database source of truth: authored papers and active
research state have different lifecycles from accepted catalogue facts.

The root `main.bib` contains 168 bibliography entries. The survey alone is
about 10,800 lines of TeX, and the minor-atlas paper is about 3,000 lines. The
survey covers several distinct dictionaries and kinds of class: graph metric,
set-family/VC, media, antimatroids, expansions, COMs and oriented matroids,
symmetry, named graph families, obstruction classes, and algorithmic results.

### Broad class inclusion atlas

`papers/survey/figures/class_inclusion_atlas.tex` contains 76 mathematical
class nodes. Their labels, visual categories, and some forbidden-minor
annotations are embedded in TikZ.

`papers/survey/figures/class_containments.json` contains 134 asserted
subclass-to-superclass pairs. The associated renderer validates endpoints,
rejects loops and duplicates, computes a transitive reduction and ranks, and
regenerates the solid arrows in the TikZ file. The current audit passes and the
134 displayed arrows are transitively irredundant.

This is useful seed data, but it is not yet a catalogue:

- class identities and display text live in TeX rather than JSON;
- an inclusion has no status, proof pointer, citation, scope, or strictness
  field;
- TikZ routing is mixed with mathematical data;
- open comparisons and non-containments are drawn separately rather than
  represented by the JSON relation;
- several nodes describe parameterized families or ambient classes whose
  semantics need more than a label.

### Projection/conditioning minor atlas

`papers/minor-atlas/inventory.json` is a narrower but richer structured source.
It currently contains:

- 24 family records;
- 48 relation records;
- 33 source records;
- 13 open questions;
- a separate boundary-study object.

Family records include closure under projection (`p`), conditioning (`c`), and
their combination (`pc`), named closures, obstruction/basis information,
assessment prose, and source IDs. Relation records include proved/open status,
source pointers, strictness, strict witnesses, negative relations, and diagram
visibility.

This inventory overlaps the broad survey atlas but serves a different purpose.
It should seed normalized database records, not be copied wholesale as one
opaque record. In particular, `p`, `c`, and `pc` closure claims are individual
mathematical facts with their own provenance, while an inclusion or
non-containment is a relation between two class records.

### Other manuscript content

The twelve manuscripts contain additional reusable kinds of facts:

- definitions and aliases;
- equivalent characterizations;
- closure under operations such as products, projections, conditionings,
  reductions, complement, expansion, and gluing;
- forbidden-minor and obstruction descriptions;
- recognition and algorithmic results;
- examples, counterexamples, and strictness witnesses;
- numerical or structural invariants;
- open questions and research status.

These should not all be forced into inclusion edges. The schema needs a small
set of normalized fact types and explicit source/proof links.

### Characterizations, invariants, and obstruction families

The survey is explicitly organized around theorem-level characterizations.
Many classes have several genuinely different membership criteria. For
example, partial cubes can be defined by isometric hypercube embeddings,
recognized by the Djokovi\'c--Winkler relation, and represented as state graphs
of well-graded families or media. The manuscripts commonly state several
conditions in one “the following are equivalent” theorem. A class therefore
needs atomic characterization records plus separate equivalence-theorem
records; a list of alternate-definition strings would lose the proof joining
them.

The structural-invariants workspace supplies a second normalized pattern. An
invariant is defined on an individual concept class or partial cube, has its own
codomain and normalization, and has laws under pc-minors, products, deletion,
or gluing. Values attached to a named class can mean different things: an exact
formula on a parameterized family, a universal bound, an attained range, or a
zero/positivity criterion. Examples in the current papers include ordering
width, majority-completion defect, resolution slope and its rounded defect,
and separation margin. In particular, an invariant condition can itself be a
class characterization; this converse assertion is stronger than merely
recording that every member of the class has a given value.

Forbidden-minor results introduce a third pattern. The target class, the
containment/minor notion, the ambient universe, the obstruction family, and the
completeness of the claimed basis are separate data. Existing results include:

- an explicit parameterized family, such as the antipodally punctured cubes
  for ample partial cubes;
- a union of several families, as in the daisy-cube basis;
- a recursive or predicate-defined complete family, as in the current
  cover-recursive obstruction theorem;
- proved infinite subfamilies whose completeness remains open, as in parts of
  the learning-space programme.

Consequently an obstruction family should not be stored only as prose on its
target class. One obstruction family may occur in several bases, one basis may
be the union of several families, and the same target can have equivalent
coordinate, intrinsic, recursive, or generative descriptions of its basis.

## Architectural lessons for the schema design

1. **Database authority should be explicit from the beginning.** Accepted
   class definitions and structured facts should move into
   `~/code/Partial-Cubes`; papers should retain exposition and long proofs.
   A checked `sync/ownership.json` should state this boundary.

2. **Class identity must be separate from presentation.** Stable short names
   should drive cross-references. TeX labels, graph labels, aliases, visual
   category, and layout hints are attributes or mappings, not identifiers.

3. **Relations need epistemic and mathematical types.** At minimum the design
   must distinguish inclusion, equality/equivalence, non-containment,
   incomparability, and an open comparison. It must also distinguish direct
   stored facts from transitive consequences and record strictness and
   witnesses without encoding them in diagram style.

4. **Unary facts deserve their own representation.** Closure under an
   operation, recognition complexity, existence of a basis, and possession of
   a representation are properties of one class. Modeling them as prose on a
   class record would make provenance, status, and querying difficult.

5. **Operations and obstructions are likely first-class entities.** Projection,
   conditioning, pc-minors, products, complements, expansions, and reductions
   recur across many papers. Likewise, named forbidden graphs/classes recur in
   basis statements. Normalizing these avoids an ever-growing set of Boolean
   columns tied to today's survey.

6. **Sources and proof locations must be structured.** A citation key alone is
   insufficient for local results. Records should be able to point to a paper,
   theorem/definition label, repository-relative path, and optionally a concise
   database proof. Long canonical proofs can remain in LaTeX.

7. **Open research state should remain distinct from accepted catalogue
   facts.** The database may display a precise open relation or candidate
   characterization, but broader questions, route ledgers, dashboards, and
   computational campaign state should stay in the LaTeX research workspace.

8. **Layout data should be downstream.** The database graph should derive its
   Hasse backbone, ranks, and overlays from records. Curator hints may still be
   useful, but TikZ paths and coordinates should not be the only store of a
   mathematical assertion.

9. **Migration needs reconciliation rather than extraction from one file.** A
   first seed pass should reconcile the 76-node survey atlas with the 24-family
   minor atlas, assign canonical IDs and aliases, and then attach definitions
   and proof locations from the manuscripts. Conflicts or weaker/stronger
   formulations should become a review queue.

10. **Validation is part of the database design.** Useful initial checks
    include unique IDs and aliases, valid references and citation keys,
    relation endpoint/type rules, acyclicity of established proper inclusions,
    consistency of strictness witnesses, proof/source requirements for
    established facts, and drift checks for generated LaTeX.

## Remaining design decisions

- How much proof text belongs in the database versus a stable pointer into a
  LaTeX manuscript.
- Which existing diagram is the first published graph: the 76-node survey
  hierarchy, the operation-focused minor atlas, or both as separate graph
  views over one catalogue.
- Whether the shared renderer should first be generalized so relation and fact
  cards are table-agnostic, or whether the initial site can use generic cards
  plus a Partial-Cubes-specific graph hook.

## Decisions from the schema discussion

- Keep `p_closed`, `c_closed`, and `pc_closed` on each class as compact,
  filterable summaries. Derive and validate them from proof-bearing operation
  result records so the booleans do not become a second mathematical source.
- Do not create a `sources` table initially. Keep the BibTeX bibliography and
  place `references`, an exact `proof_source`, and any concise `proof` directly
  on each mathematical claim.
- Do not create a `questions` table initially. Broad questions and research
  status remain prose in the LaTeX papers, ledgers, and dashboards. Precisely
  formulated open relations or candidate characterizations may still appear
  in their mathematical tables with an open status.
