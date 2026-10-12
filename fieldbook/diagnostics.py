"""Report distinct adaptation gates without replacing the upstream checks."""

from dataclasses import dataclass

import numpy as np
from cpg2unstructured import (
    active_columns,
    boundary_mask,
    cell_quality,
    check_validity,
    pillar_lattice,
)
from fault_projections.adapt import fault_cell_map
from fault_projections.corridor import cell_corridor


@dataclass(frozen=True)
class GateReport:
    valid: bool
    quality_top: float
    quality_bottom: float
    quality_floor_ok: bool
    corridor_ok: bool
    chain_continuous: bool
    rim_max_displacement: float
    rim_within_tolerance: bool

    @property
    def ok(self) -> bool:
        return (
            self.valid
            and self.quality_floor_ok
            and self.corridor_ok
            and self.chain_continuous
            and self.rim_within_tolerance
        )


def evaluate(original, adapted, trace, targets, *, min_sj: float = 0.3,
             rim_tolerance: float = 1e-6) -> GateReport:
    """Check grid, quality, corridor, target chain and rim, each separately.

    ``targets`` must be the actual targets used by the low-level solver.
    This cannot infer them from ``adapt_grid``'s result, which does not
    expose its target object. The chosen floors are caller policy.
    The rim gate enforces a *fixed* rim, so it deliberately rejects even
    tangential sliding. For permitted slide, separately measure deviation
    from the original physical outline as in lesson 09.
    """
    if min_sj <= 0 or rim_tolerance < 0:
        raise ValueError("min_sj must be positive and rim_tolerance nonnegative")
    before = pillar_lattice(original)
    after = pillar_lattice(adapted)
    validity = check_validity(after, before, active=active_columns(original))
    strict = check_validity(
        after, before, active=active_columns(original),
        min_scaled_jacobian=min_sj,
    )
    rim = boundary_mask(original, scope="active")
    displacement = np.linalg.norm(after.top[rim] - before.top[rim], axis=-1)
    max_move = float(displacement.max()) if displacement.size else 0.0
    return GateReport(
        valid=validity.ok,
        quality_top=cell_quality(after, plane="top").sj_min,
        quality_bottom=cell_quality(after, plane="bot").sj_min,
        quality_floor_ok=strict.ok,
        corridor_ok=cell_corridor(adapted, trace).ok,
        chain_continuous=fault_cell_map(targets, original.cart_dims).is_continuous,
        rim_max_displacement=max_move,
        rim_within_tolerance=max_move <= rim_tolerance,
    )
