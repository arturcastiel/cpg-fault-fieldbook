"""Convert a corner-point grid to explicit faces, properties and a cell path."""

from cpg2unstructured import (
    attach_props,
    build_topology,
    cell_adjacency,
    from_grdecl,
    read_grdecl,
    shortest_path,
)
from fieldbook.fixtures import quickstart

source = read_grdecl(quickstart())
mesh = attach_props(build_topology(from_grdecl(source)), source)
assert mesh.n_cells == 70
assert mesh.faces.neighbours.shape[1] == 2  # British spelling in upstream API
assert "PORO" in mesh
distance, path = shortest_path(cell_adjacency(mesh, weights="unit"), 0, mesh.n_cells - 1)
assert path[0] == 0 and path[-1] == mesh.n_cells - 1
print("nodes:", mesh.nodes.shape[0], "faces:", mesh.faces.neighbours.shape[0])
print("active cells:", mesh.n_cells, "path length:", distance)
