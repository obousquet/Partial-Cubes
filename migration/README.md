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

When a reviewed claim changes after the main freeze, its exact current source
file may be preserved as a supplemental snapshot under a new snapshot ID.
`supplemental_source_archive.json` records these targeted refreshes, and
`check_source_manifest.py` verifies their archived bytes on every migration
check. Supplemental entries overlay the same logical path in extraction and
coverage indexing, with the latest listed entry taking precedence. Batch
ledgers must state which claims use the main snapshot and which use a
supplemental refresh.

The ownership states are:

- `latex`: extracted candidate; the source document remains authoritative;
- `shadowed`: reviewed database representation exists but has not replaced a
  shared LaTeX catalogue block;
- `database`: the normalized fact is database-owned and any selected LaTeX
  consumer uses checked generated output.

Source papers remain scholarly artifacts in every state. Full proofs,
paper-specific statements, narrative, notation, diagrams, and organization
stay in the LaTeX repository.

## Current minor-atlas updates

Batch 200 reviews the atlas state of 20 October 2026 in the separate byte
snapshot `snapshots/2026-10-20-b200/`. Its manifest records source commits,
working-tree state and file hashes. `minor_atlas_graph_classes.json` points
to that inventory: graph checks cover all 35 Ample-family rows, all 80
relationships at their stated strength, and every affirmative or negative
minor-closure cell. Signed bouquets reuse the existing conditional-antimatroid
identity; ACC likewise reuses its existing class.

The original 24-family seed inventory, numbered JSON locators, review queues
and manuscript-label audit deliberately keep their historical frozen inputs.
Replacing that inventory in place would change the meaning of earlier
`#/relations/N` provenance. Consequently the seed audit's 37/37 minor-atlas
labels do not assert coverage of every theorem in the current expanded paper.
The current atlas graph coverage is checked independently by `make check`.
