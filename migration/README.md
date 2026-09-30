# Migration workspace

This directory contains the review boundary between the LaTeX research
workspace and the public database. Extractors only write candidate artifacts
under `migration/generated/`; they never write records under `data/`.

The generated files are reproducible from the pinned inputs in
`source_manifest.json`. Human decisions belong in reviewed batch ledgers and,
from Batch 001 onward, in the class crosswalk and claim dispositions. A later
extraction may refresh candidates but must not overwrite those decisions.

Run the first-batch pipeline with:

```bash
make migration-extract
make migration-check
```

`LATEX_ROOT` defaults to `~/latex/PartialCubes` and may be overridden on either
command. The check target is read-only.

The ownership states are:

- `latex`: extracted candidate; the source document remains authoritative;
- `shadowed`: reviewed database representation exists but has not replaced a
  shared LaTeX catalogue block;
- `database`: the normalized fact is database-owned and any selected LaTeX
  consumer uses checked generated output.

Source papers remain scholarly artifacts in every state. Full proofs,
paper-specific statements, narrative, notation, diagrams, and organization
stay in the LaTeX repository.
