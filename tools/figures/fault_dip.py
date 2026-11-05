"""Regenerate the fault-dip figures of the wiki from lesson 31.

usage: python tools/figures/fault_dip.py [OUTDIR]

Runs examples/31_fault_dip.py --out into a temporary directory (the lesson draws the figures
itself, from the grids it checks) and copies its PNGs to OUTDIR, by default
docs/assets/fault-dip/. Needs matplotlib.
"""

import shutil
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[2]
NAMES = {
    "31_A_rigid_60.png": "rigid-vs-decay.png",
    "31_B_why_neighbours_move.png": "why-neighbours-move.png",
    "31_D_tip_55.png": "tip-fade.png",
    "31_G_graben.png": "graben.png",
}


def main(out: Path) -> int:
    out.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory() as temporary:
        work = Path(temporary)
        subprocess.run([sys.executable, str(ROOT / "examples" / "31_fault_dip.py"),
                        "--out", str(work)], check=True)
        for source, target in NAMES.items():
            shutil.copyfile(work / source, out / target)
            print(out / target)
    return 0


if __name__ == "__main__":
    target = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "docs" / "assets" / "fault-dip"
    raise SystemExit(main(target))
