# Migration workspace

This directory contains the review boundary between the LaTeX research
workspace and the public database. Extractors only write candidate artifacts
under `migration/generated/`; they never write records under `data/`.

The generated files are reproducible from the pinned inputs in
`source_manifest.json`. Human decisions belong in reviewed batch ledgers and,
from Batch 001 onward, in the class crosswalk and claim dispositions. A later
extraction may refresh candidates but must not overwrite those decisions.

Freeze a new explicitly named snapshot, then run the extraction pipeline with:

```bash
make migration-freeze
make migration-extract
make migration-check
```

`LATEX_ROOT` defaults to `~/latex/PartialCubes` and may be overridden on either
command. Change `snapshot_id` in `config.json` before freezing a later source
state. `migration-freeze` refuses to overwrite an existing byte archive.
The check target is read-only.

Clean tracked inputs are recoverable from the repository commit recorded in
`source_manifest.json`. Every selected modified or untracked input is copied
byte-for-byte under `snapshots/<snapshot_id>/`. Extractors read the frozen
version even if the live research workspace changes later. The direct source
corpus excludes `research/` and `archive/` directories; those artifacts remain
in the LaTeX repository but cannot silently serve as proof provenance.

The ownership states are:

- `latex`: extracted candidate; the source document remains authoritative;
- `shadowed`: reviewed database representation exists but has not replaced a
  shared LaTeX catalogue block;
- `database`: the normalized fact is database-owned and any selected LaTeX
  consumer uses checked generated output.

Source papers remain scholarly artifacts in every state. Full proofs,
paper-specific statements, narrative, notation, diagrams, and organization
stay in the LaTeX repository.
