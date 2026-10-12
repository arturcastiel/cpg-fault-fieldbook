"""Reproducible case metadata; no interpretation of simulator readiness."""

import hashlib
import importlib.metadata
import platform
from dataclasses import asdict
from pathlib import Path

import numpy as np
import scipy

from fieldbook.cli import upstream_revisions


def case_record(input_path: Path, trace: np.ndarray, settings: dict,
                morph, gates) -> dict:
    return {
        "schema": "cpg-fault-fieldbook-case-v1",
        "input": {
            "name": input_path.name,
            "sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
            "coordinate_units": "toy fixture: arbitrary consistent units",
        },
        "versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "cpg2unstructured": importlib.metadata.version("cpg2unstructured"),
            "fault-projections": importlib.metadata.version("fault-projections"),
            "revisions": upstream_revisions(),
        },
        "fault_trace_xy": trace.tolist(),
        "settings": settings,
        "solver": {
            "ok": morph.ok,
            "weight_final": morph.w,
            "n_moved": morph.n_moved,
            "residual_mean": float(np.mean(morph.residuals)),
        },
        "gates": asdict(gates),
        "verdict": "synthetic_geometry_gates_pass" if gates.ok else "blocked_by_geometry_gate",
        "simulator_validated": False,
    }
