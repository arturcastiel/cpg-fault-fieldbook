"""Read a text deck; inspect logical and active dimensions."""

from cpg2unstructured import read_grdecl
from fieldbook.fixtures import quickstart

grid = read_grdecl(quickstart())
assert grid.cart_dims == (5, 5, 3)
assert grid.n_pillars == 36
assert int(grid.actnum.sum()) == 70
assert grid.coord.shape == (36, 6)
print("cells:", grid.cart_dims, "active:", grid.actnum.sum(), "pillars:", grid.n_pillars)
print("properties:", sorted(grid.props))
