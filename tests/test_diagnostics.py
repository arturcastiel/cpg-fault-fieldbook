"""Do not confuse a valid grid with a quality-approved or rim-safe grid."""

from dataclasses import replace

import numpy as np

from cpg2unstructured import morph_pillars, read_grdecl
from fault_projections import adapt_grid, targets_from_polyline
from fieldbook.diagnostics import evaluate
from fieldbook.fixtures import quickstart


def test_polish_passes_all_gates():
    grid = read_grdecl(quickstart())
    trace = np.array([[1., 1.], [3., 3.]])
    targets = targets_from_polyline(grid, trace)
    result = morph_pillars(grid, targets, boundary="pinned", w0=0.05,
                           polish="barrier", on_invalid="raise")
    report = evaluate(grid, result.grid, trace, targets)
    assert report.ok and report.quality_top >= .3 and report.quality_bottom >= .3


def test_valid_unpolished_grid_fails_explicit_quality_gate():
    grid = read_grdecl(quickstart())
    trace = np.array([[1., 1.], [3., 3.]])
    targets = targets_from_polyline(grid, trace)
    result = morph_pillars(grid, targets, boundary="pinned", w0=.05,
                           on_invalid="raise")
    report = evaluate(grid, result.grid, trace, targets)
    assert report.valid and not report.quality_floor_ok and not report.ok


def test_free_boundary_can_pass_solver_but_fail_rim_gate():
    grid = read_grdecl(quickstart())
    trace = np.array([[1., 1.], [3., 3.]])
    targets = targets_from_polyline(grid, trace)
    result = morph_pillars(grid, targets, boundary="free", w0=.05,
                           polish="barrier", on_invalid="raise")
    report = evaluate(grid, result.grid, trace, targets)
    assert result.ok and report.rim_max_displacement > 0
    assert not report.rim_within_tolerance and not report.ok


def test_sheared_slide_is_not_outline_safe():
    original = read_grdecl(quickstart())
    coords = original.coord.copy()
    coords[:, 0] += .5 * coords[:, 1]
    coords[:, 3] += .5 * coords[:, 4]
    grid = replace(original, coord=coords)
    trace = np.array([[1.1, 2.2], [6.7, 3.4]])
    result = adapt_grid(grid, trace, boundary="slide", w0=.05,
                        regenerate_faults=False, on_invalid="raise")
    assert result.ok
    # Physical outline check in lesson 09 is more specific than this
    # displacement metric: even tangential rim motion can be acceptable.
    assert np.linalg.norm(result.grid.coord - grid.coord) > 0
