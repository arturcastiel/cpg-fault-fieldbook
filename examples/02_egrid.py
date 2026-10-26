"""Read binary EGRID, compare MAPAXES-applied and raw-local coordinates."""

import warnings

import numpy as np
from cpg2unstructured import from_grdecl, read_egrid
from fieldbook.fixtures import egrid

plain = read_egrid(egrid())
world = read_egrid(egrid(rotated=True))
with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    local = read_egrid(egrid(rotated=True), apply_mapaxes=False)
assert plain.cart_dims == (2, 2, 2) and int(plain.actnum.sum()) == 8
assert world.mapaxes_applied and not local.mapaxes_applied
assert caught and any("MAPAXES" in str(item.message) for item in caught)
assert not np.allclose(world.coord, local.coord)
assert from_grdecl(world).cart_dims == world.cart_dims
print("raw local first pillar:", local.coord[0, :2])
print("MAPAXES world first pillar:", world.coord[0, :2])
print("grid source:", world.source, "geometry dtype:", world.coord.dtype)
