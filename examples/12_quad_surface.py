"""Run the upstream coherent quad example without writing into the checkout."""

import argparse
import subprocess
import sys
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
upstream = Path(__file__).resolve().parents[2] / "fault-projections"
script = upstream / "examples" / "quad_mesh_fault.py"
if not script.is_file():
    parser.error("this advanced lesson needs the sibling fault-projections checkout")
subprocess.run([sys.executable, str(script), "--out", str(args.out)], check=True)
assert (args.out / "quad_mesh_fault.FAULTS").is_file()
assert (args.out / "quad_mesh_fault.vtk").is_file()
print("Coherent quad example produced FAULTS and VTK.")
