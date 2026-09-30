# Partial Cubes Database

This repository is the structured-data and website counterpart of the LaTeX
workspace at `~/latex/PartialCubes`. It catalogues classes of partial cubes,
their characterizations, inclusions, closure properties, invariants, forbidden
pc-minor bases, examples, and supporting results.

The database is intentionally empty at first. Its schemas and validation are in
place so the existing mathematical material can be migrated record by record.

## Repository structure

- `data/` contains the authoritative JSON database and BibTeX bibliography.
- `data/*/schema.json` defines the twelve database tables.
- `data/make_graph.py` builds the class-hierarchy and minor/closure graph views.
- `scripts/validate_data.py` performs structural and mathematical consistency
  checks beyond the editor's field schemas.
- `sync/` records the ownership boundary with `~/latex/PartialCubes` and
  generates database-owned LaTeX catalogue fragments.
- `docs/` is the generated static site committed on `main`. GitHub Pages
  publishes it from the `/docs` folder, matching Combinatorial-Parameters.

## Shared renderer

The website uses [`math_database`](https://github.com/obousquet/math_database).
With the sibling checkout at `~/code/math_database`, run:

```bash
make validate
make build
make serve
```

The local server is then available at <http://localhost:8080>. Override the
renderer location with `MATH_DATABASE_DIR=/path/to/math_database` when needed.
Run `make check` before committing database changes, then commit the regenerated
`docs/` tree together with the source records.

## Data conventions

- Each record is one `NNN_short_name.json` file.
- Cross-references use `#table/short_name`.
- Mathematical prose is Markdown with dollar-delimited TeX.
- Bare symbols use the `latex` field type; renderers add math delimiters.
- Claim records distinguish mathematical type from epistemic `status`.
- `references` contains BibTeX citations in `\cite{key}` or `[@key]` form;
  `proof_source` uses a stable theorem
  locator such as
  `PartialCubes:papers/survey/main.tex#thm:partial-cube-wellgraded-medium`.
- The class fields `p_closed`, `c_closed`, and `pc_closed` are cached summaries
  checked against proof-bearing `operation_results` records.
- Broad open questions and research campaign state remain in the LaTeX
  workspace.

See [SCHEMA_PROPOSAL.md](SCHEMA_PROPOSAL.md) for the design rationale and
[STRUCTURE_NOTES.md](STRUCTURE_NOTES.md) for the repository survey.

## LaTeX synchronization

The database owns accepted, normalized catalogue entries and their shared
bibliography. Unmigrated content remains LaTeX-owned. The LaTeX repository
continues to own every paper as a scholarly artifact, including exposition,
paper-specific statements, notation, full proofs, diagrams, and research
ledgers. Only a manuscript deliberately selected as an active catalogue
consumer imports generated fragments; stable theorem papers may remain
unchanged. Generate a catalogue fragment explicitly, for example:

```bash
python3 sync/generate_latex_catalog.py \
  --section class-definitions \
  --output ~/latex/PartialCubes/includes/generated/class_definitions.tex
```

Use `--check` to detect drift without writing the output.

See [MIGRATION_PLAN.md](MIGRATION_PLAN.md) for the record-by-record migration,
preservation checks, and two-repository handoff.
