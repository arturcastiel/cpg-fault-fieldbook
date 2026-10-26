"""Project a 2D authored fault onto an unchanged mesh and emit FAULTS text."""

import numpy as np
from cpg2unstructured import build_topology, from_grdecl, projected_topology, read_grdecl
from fault_projections import (
    build_projected_graph,
    faults_text,
    path_to_faults,
    project_polylines,
    stitch_paths,
)
from fieldbook.fixtures import quickstart

mesh = build_topology(from_grdecl(read_grdecl(quickstart())))
top = projected_topology(mesh, "xy")
graph = build_projected_graph(top)
line = np.array([[1., 1.], [3., 3.]])
paths = stitch_paths(project_polylines([line], graph), graph)
faces = [path_to_faults(path, top) for path in paths]
text = faults_text(faces, nz=mesh.cart_dims[2])
assert "FAULTS" in text and "FAULTDIM" in text and "'Y'" in text
print(text)
