# 4. Judge the adaptation; export with restraint

An upstream `ok` result checks validity, not the four independent
engineering questions. Keep the exact targets used by the low-level
solve; check validity **and** a chosen quality floor, a containing
face-connected cell corridor, lattice-chain continuity, and physical
rim motion. Run `python -m fieldbook.cli run 10_four_gates`.

{{ example:10_four_gates.py }}

The `fieldbook.diagnostics.evaluate` helper **calls upstream checks**;
it does not implement a second morph or infer targets from `adapt_grid`
(which does not expose them). A chain can be continuous without its
fault face geometry satisfying your field-specific requirements.
Its rim tolerance is a **fixed-rim policy**: even tangential sliding
fails. If sliding is intentionally allowed, measure departure from
the actual physical outline instead (lesson `09_sheared_rim`), not
just the distance every rim pillar moved.

Only after review, export a **separate** grid INCLUDE. Run:
`python -m fieldbook.cli run 11_export_includes --out output/one`.

{{ example:11_export_includes.py }}

With `regenerate_faults=False` there is no FAULTS INCLUDE. The
master deck is never rewritten. You must review INCLUDE order,
ACTNUM, MAPAXES, other deck keywords and simulator parsing yourself.
Neither a generated grid INCLUDE nor four passing gates establishes
simulator-readiness.

## Start with a read-only field-deck intake

Use `python -m fieldbook.cli run 13_field_deck_intake --deck FIELD.GRDECL`
or provide an `.EGRID` path. This prints dimensions, ACTNUM, coordinate
frame, quality on both planes and recognized/unrecognized property names;
it does **not** build full topology, rewrite a master deck or declare
the input valid. The code is shown below:

{{ example:13_field_deck_intake.py }}

## Visual inspection uses the same scale for before and after

Run `python -m fieldbook.cli run 14_visualize --out output/figures`.
The drawings show top-plane pillar lines and the authored fault;
they do not show bottom-plane quality, simulator transmissibility
or guarantee fault-chain connectivity. The example prints all
four gates beside the SVG paths:

{{ example:14_visualize.py }}
