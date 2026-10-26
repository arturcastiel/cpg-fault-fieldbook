"""Project a planar polygon onto 3D cell faces (no grid deformation)."""

import numpy as np
from cpg2unstructured import build_topology, from_grdecl, read_grdecl
from fault_projections import faults_text_3d, project_fault_surface
from fieldbook.fixtures import quickstart

mesh = build_topology(from_grdecl(read_grdecl(quickstart())))
polygon = np.array([[2., 0., 0.], [2., 5., 0.], [2., 5., 3.], [2., 0., 3.]])
surface = project_fault_surface(polygon, mesh)
assert surface.face_ids.size == 14
text = faults_text_3d({"F1": surface}, mesh)
assert "FAULTS" in text
print("projected 3D cell faces:", surface.face_ids.size)
