# 8. Stress-test the fault dip: 22 hard cases, judged independently

Lesson 31 shows what `dip_grid` does on clean grids. This lesson tries to
break it. It mixes the dip with curved and S-shaped traces, slides, growth
faults, blind tips, holes, junctions, skewed columns, shallow angles, random
angles, seeded random grids and a field deck.

**The judge is not the library.** Every grid is rebuilt from raw COORD and
ZCORN by a cell checker in the lesson itself (`audit`), which shares no code
with `fault_projections`. For every cell it measures:

- the scaled Jacobian of the trilinear hexahedron at 45 points;
- the volume;
- whether the plan quadrilateral of the column's four pillars folds, at
  depths where all four carry corners.

A cell is **broken** when any of these holds:

- it inverts;
- its worst Jacobian falls below 0.05 *and* below what it was (a 1 % wobble
  allowed);
- it loses 90 % of its volume;
- its pillars fold in plan.

Only active cells are judged. A fold under a column with no active cell is
reported separately.

Dips with the default `boundary="fixed"` are also checked, from raw COORD, to
hold the grid's outline: no rim pillar moves across its side, and the corners
do not move.

**Every accepted grid below has zero broken cells.** Every refusal is shown,
together with what the checker says about the grid the gate refused.

Run `python -m fieldbook.cli run 32_fault_dip_stress` (about 80 s). Add
`--out DIR` for the figures, or `--case B,C` to run only those cases.

> **Library calls.** `fault_projections.Dip`, `dip_grid`, `Throw`,
> `slide_grid` and `adapt_grid`. The method reference is
> [Fault dip](../functions/fault-dip.md).

## How to read the figures

Each figure shows the plan BEFORE and AFTER, coloured by each column's worst
3-D scaled Jacobian (1 is a perfect box; 0 or less is inverted). It also
shows sections through the grid at each step: the original mesh, each dip
and slide, and the final mesh. The fault is drawn in black: solid at the
top, dashed at the base. **A red outline marks a broken cell**; no figure of
an accepted grid has one.

## The cases

| Case | What it tries | Result |
|---|---|---|
| A | One rim-to-rim fault at 89, 75, 60, 45, 35 and 30° (`boundary="free"`) | Exact at every angle; hanging-wall volumes change by at most 4e-15; min 3-D Jacobian 0.998 → 0.457 as the angle falls; 0 broken |
| B | An arc adapted into the pillars, dipped 60°, slid 20 m | 60.00° all along; throw 20.000 m along 341 m of arc; 0 broken |
| C | S-shaped trace, angle table 50 → 80°, throw table 8 → 20 m | 50.0° at one end, 65.3° mid-way, 80.0° at the other; along row 10 the largest step between neighbours is 0.48 m (no tear inside the bend); 0 broken |
| D | Stair-stepped diagonal FAULTS (23 faces, not adapted), 50°, slid 15 m | Dips along the smoothed strike; 0 broken |
| E | A fault with two tips; an angle per pillar 55 → 80 → 55°; throw split 50/50 | Decay mode, fade 110 m; throw 0 at the tips and 15 m mid-fault; 0 broken |
| F | Graben and horst: four faults | At 60° the graben faults cross 5.2 m below the cells, and the slide is **refused** (36 cells would fold; the checker agrees). Widening ×16 cannot help. At 75° they cross 99 m below and everything passes |
| G | Domino: four parallel faults at 50, 55, 60, 65°, 10 m each | The blocks step down 0, 10, 20, 30, 40 m, measured exactly; 0 broken |
| H | T-junction: a fault ending on another, both dipped | The shared pillar holds both angles: 65.0000° across MAIN and 70.0000° across BRANCH; 0 broken |
| I | Bedding dipping 22° along strike (the pivot depth varies by 120 m) | Hanging-wall volumes change by at most 4.6e-14 (`boundary="free"`) |
| J | Every pillar pre-tilted 35 m east; the fault turned to dip west | 74.9° east → 60.000000° west; 0 broken |
| K | Dip 60, slide 15, re-dip 75, slide 10, re-dip 50 | Total throw 25.000000 m: a dip moves pillars, never depths; 0 broken |
| L | Thin and pinched-out layers, a zero-thickness layer, 15 % inactive cells | 408 zero-volume and 270 inactive cells, none judged, none broken |
| M | Skewed, uneven columns (3.7 to 39.5 m wide, sheared) | The fault strikes 21.8° off north and dips along its own normal; 0 broken |
| N | A dome with a ring fault, dipped 70° inward, collapsed 25 m (a caldera) | Rigid inside the ring; throw 25.000 m all round; 0 broken |
| O | Shallow angles on bedding dipping 14° toward the fault | See below: 15° with the rim fixed is **refused**; 30° passes |
| P | Random angles 40–85° along a curved trace | Noise at each pillar is **refused**: the fault's base runs backwards along strike on 8 of 16 segments, so it folds onto itself. Five random control angles, interpolated, are exact |
| Q | 18 m on a short 40° fault in 5 m layers | Without the slide's cell gate, 2 cells invert, and the checker finds them. The gate refuses them; bisect manages 87.4 %. `dip_grid(slide=...)` widens the fades ×2 and fits the full 18 m |
| R | Bedding rising 35° toward a fault dipping 30° | The gate refuses (264 plan folds), but the checker finds the 3-D cells sound. **A conservative refusal**: the plan view of a cell's top flips where bedding is steeper than the fault |
| S | Growth fault: throw 10 m at the top growing to 16 m at the base | Hanging-wall layers thicken by 10.0, 4.0, 10.0, 5.0, 5.0 %, exactly as computed by hand; 0 broken |
| T | Blind fault cutting layers 3–5 only | Offset 0, 0, 0, 3, 9, 12 m at horizons 0–5; 0 broken |
| U | Seeded fuzz: 2 seeds × 20 random grids; dip → slide → re-dip | 90 steps accepted, every one passed by the checker. 11 refused, and the checker agrees they break cells. 4 refused where the checker would have passed (conservative). **0 disagreements the other way** |
| V | Norne: E_01 and BC set to 60° | E_01 is exact on 28 pillars and 13,015 cells move, with 0 broken. BC is **refused** (1 degenerate cell), and the checker agrees (2 degenerate) |

## A. Angle sweep

![Angle sweep](../assets/fault-dip-stress/32_A_angle_sweep.png)

## B. A curved trace, dipped and slid

![Curved arc](../assets/fault-dip-stress/32_B_curved_arc.png)

## C. An S-shaped trace with angle and throw tables

![S trace](../assets/fault-dip-stress/32_C_s_trace_angle_table.png)

Arc length runs from the chain end with the smaller (I, J), here the north
end. A rigid slide on a curved fault used to read each pillar's throw at its
nearest point on the trace. Inside a bend that point jumps across the
medial axis, which tore the hanging wall by 9 m over one column. This case
found the bug; the arc length is now spread smoothly.

## D. A stair-stepped diagonal fault

![Staircase](../assets/fault-dip-stress/32_D_staircase.png)

## E. Two tips, an angle per pillar

![Tips](../assets/fault-dip-stress/32_E_tips_angle_per_pillar.png)

## F. Graben and horst

![Graben and horst](../assets/fault-dip-stress/32_F_graben_horst.png)

With the middle pivot at 60°, graben pillars 6 and 7 cross at 2114 m, 5.2 m
below the cells. A slide moves the cells down onto that crossing:

- `slide_grid` refuses it (36 folded);
- without the cell gate, the checker finds 36 inverted;
- `dip_grid(..., slide=...)` widens ×16 and still refuses, because the
  faults themselves meet within the throw.

At 75° they cross 99 m below, and the slide passes. Throws hold exactly:
15 m on the graben faults and 10 m on the horst faults.

## G. Domino

![Domino](../assets/fault-dip-stress/32_G_domino.png)

## H. T-junction

![T-junction](../assets/fault-dip-stress/32_H_t_junction.png)

This case found a bug. The dip damped the other fault's pillars to zero,
along its strike too, which pinched MAIN's 25 m edges to about 1 m at the
junction, and the dip was refused at any widening. Only the component across
the other fault is removed now.

## I. Bedding dipping along strike

![Bedding along strike](../assets/fault-dip-stress/32_I_bedding_along_strike.png)

An earlier version spread the pivot off the fault with a weighted mean,
pulled toward the fault's middle near its ends. That strained the "rigid"
hanging wall. The spread is now a local linear fit along strike.

## J. Pre-tilted pillars, dip flipped

![Pre-tilted](../assets/fault-dip-stress/32_J_pretilted_flip.png)

## K. Re-dip sequence

![Re-dip](../assets/fault-dip-stress/32_K_redip_sequence.png)

## L. Pinch-outs and holes

![Pinch-outs](../assets/fault-dip-stress/32_L_pinchout_holes.png)

## M. Skewed columns

![Skewed](../assets/fault-dip-stress/32_M_skewed_columns.png)

## N. A caldera

![Caldera](../assets/fault-dip-stress/32_N_dome_ring_caldera.png)

## O. Shallow angles, and the fixed outline

![Shallow](../assets/fault-dip-stress/32_O_shallow_widen_bisect.png)

- **15° with the default fixed outline:** refused even widened ×16. The
  pillars must shift ±187 m at the top and base, and the east side, parallel
  to the fault and 300 m from it, holds still.
- **30°, default:** fades widened ×2, every cell valid, the outline held.
- **15° with `boundary="free"`:** fades widened ×2 (one warning), footwall
  fade 1,172 m, every cell valid.
- **15°, free, `widen=0`:** `DipError`.
- **15°, free, `widen=0, on_invalid="bisect"`:** 76.4 % of the change. The
  fault dips **19.3°, not 15°**, and `scale=0.764` says so.
- **8°, free:** refused even widened ×16.

The fixed outline costs room: a shallow fault close to a parallel side may be
refused where `boundary="free"` passes.

## P. Random angles

![Random angles](../assets/fault-dip-stress/32_P_random_angles.png)

## Q. The slide's cell gate

![Slide gate](../assets/fault-dip-stress/32_Q_slide_cell_gate.png)

Corners slide along their pillars, so on leaning pillars a throw also moves
them sideways. A large throw near a fault tip can turn cells inside out while
ZCORN stays in order. The fuzzer used to build this lesson found such grids,
and `slide_grid` gained `min_jacobian` (3-D cell shape, on by default) to
refuse them. The full 18 m fits only after `dip_grid(slide=...)` widens the
fades. The cells then change volume by up to 43 %.

## R. A conservative refusal

![Conservative](../assets/fault-dip-stress/32_R_conservative_refusal.png)

The gate is allowed to be cautious, but not wrong in the other direction.
Here it refuses a grid whose 3-D cells are sound (min Jacobian 0.63), because
the plan view of a cell's top flips.

## S. Growth fault

![Growth](../assets/fault-dip-stress/32_S_growth_fault.png)

## T. Blind fault

![Blind](../assets/fault-dip-stress/32_T_blind_fault.png)

## U. Seeded fuzz

![Fuzz](../assets/fault-dip-stress/32_U_seeded_fuzz.png)

| Seed | Accepted | Refused | Conservative | Disagree |
|---|---|---|---|---|
| 2024 | 45 | 5 | 2 | 0 |
| 7 | 45 | 6 | 2 | 0 |

"Disagree" would be a grid the gate accepted and the checker breaks: there
are none. In the larger fuzzer used while building this lesson, the
pivot-spread fix (case I) raised the share of dips accepted from about 91 %
to 99 %.

## V. Norne

![Norne](../assets/fault-dip-stress/32_V_norne.png)

- **E_01 → 60°:** exact on 28 pillars; decay mode; no widening. The checker
  finds 0 broken cells (and 1 fold in all-inactive columns, not judged).
  13,015 cells moved; their min 3-D Jacobian goes from 0.279 to 0.442.
- **BC → 60°:** refused. One cell that was already poor (corner Jacobian
  0.021) would drop to 0.0187. The checker agrees and finds 2 degenerate
  cells.

{{ example:32_fault_dip_stress.py }}
