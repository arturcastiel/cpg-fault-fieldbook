# Workflows: which tool, in which order

The [Learn](../learn/grids.md) pages teach one tool at a time. This
section starts from the job instead. Each workflow is a flowchart of
the calls in order. It links to the lessons and pages that show each
step, and shows a computed before-and-after picture.

## Which workflow?

Answer two questions: should the grid's geometry change, and if so,
what should move?

![Which workflow? A decision chart from a grid and a fault to the eight workflows](../assets/workflows/chooser.png)

| # | Workflow | Changes | Use it when | Lessons |
| --- | --- | --- | --- | --- |
| 1 | [Stair-step projection](stair-step.md) | FAULTS only | you have a map trace and the grid must not move | 03, 04, 28 |
| 2 | [Full 3D projection](full-3d.md) | FAULTS only (X, Y and Z faces) | you have a 3D surface: a polygon, facets or a quad mesh | 05, 13, 25 |
| 3 | [2D / 2.5D deformation](deformation.md) | COORD (+ regenerated FAULTS) | the grid should follow the fault: in plan, or tilted into it | 06–12, 16, 30 C |
| 4 | [Fault slide](slide.md) | ZCORN only | the layers should be offset across a fault (the throw) | 19–22, 23–29 |
| 5 | [Slide + pillar shift](slide-pillars.md) | COORD, then ZCORN | the pillars should move and the layers slide along them | 21 B, 30 A, 30 C, 31 I |
| 6 | [Set a dip angle](dip-angle.md) | COORD only | a fault should dip at a given angle | 31, 32 |
| 7 | [Fault zones](fault-zones.md) | nothing, or a local edit | you need the cells within ±D of a fault, or want 4–6 kept inside them | 33 |
| 8 | [A field deck, end to end](field-deck.md) | separate INCLUDE files | any of 3–7 on a deck on disk | 11, 12, 14, 15, 17, 18, 29, 30, 31 J, 33 K |

What each workflow writes into the grid:

| Keyword | 1 | 2 | 3 | 4 | 5 | 6 | 7 (select) | 7 (confine) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| COORD (pillars) | – | – | changes | – | changes | changes | – | changes, in the zone |
| ZCORN (corner depths) | – | – | – | changes | changes | – | – | changes, in the zone (a slide) |
| FAULTS | written | written (3D) | regenerated | kept | regenerated (`adapt_deck`) or kept (`dip_grid`) | kept | kept | kept |
| a region keyword | – | – | – | – | – | – | `FLTZONE` | – |

## The workflows at a glance

The same synthetic grid (16 × 8 × 6 cells, 50 m by 15 m) run through
every workflow. Each panel is the grid a call returned. The sections
cut along the south faces of the cells in row J = 4 (y = 150 m), at
true scale, with the same depth range in every panel.

![The workflows at a glance: one grid through every workflow](../assets/workflows/at-a-glance.png)

## How to read the flowcharts

| Shape | Meaning |
| --- | --- |
| blue rounded box | what you start from |
| orange slanted box | data or a file |
| white box | a library call; blue text is the call itself |
| purple cut-corner box | a check or gate |
| yellow diamond | a decision |
| green rounded box | the result |
| dashed grey box | a note, or where to read more |

Grey italic text names the lessons and pages: "lesson 31 E" is case E
of lesson 31. Each lesson runs with
`python -m fieldbook.cli run NN_name`; the lesson table on each
workflow page gives the exact command.

## How these figures are made

`python tools/figures/workflows.py` draws every figure in this
section into `docs/assets/workflows/`. The flowcharts are drawn with
matplotlib. The before-and-after pictures are not sketches: the
script calls `fault_faces_for`, `project_fault_surface`, `adapt_grid`,
`slide_grid`, `dip_grid` and `fault_zone` on the synthetic grid,
draws the grids they return, and prints the measured numbers into
the titles. It asserts what it draws: the adapted grid is valid, the
slide opens only FAULTS faces, and nothing outside the zone moves.
