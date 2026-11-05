"""Regenerate the fault-dip and fault-zone figures of the wiki from lessons 31, 32 and 33.

usage: python tools/figures/fault_dip.py [ASSETS_DIR]

Runs each lesson with --out into a temporary directory (the lessons draw the figures
themselves, from the grids they check) and copies the PNGs under ASSETS_DIR, by default
docs/assets/: lesson 31 to fault-dip/ (renamed), lesson 32 to fault-dip-stress/ and lesson 33
to fault-zone/. Needs matplotlib; about 5 minutes (Norne cases run when it is present).
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
LESSONS = [("31_fault_dip.py", "fault-dip"), ("32_fault_dip_stress.py", "fault-dip-stress"),
           ("33_fault_zone.py", "fault-zone")]


def main(assets: Path) -> int:
    for script, folder in LESSONS:
        out = assets / folder
        out.mkdir(parents=True, exist_ok=True)
        with TemporaryDirectory() as temporary:
            work = Path(temporary)
            subprocess.run([sys.executable, str(ROOT / "examples" / script), "--out", str(work)],
                           check=True)
            pngs = sorted(work.rglob("*.png"))
            for source in pngs:
                if script.startswith("31") and source.name not in NAMES:
                    continue
                target = out / NAMES.get(source.name, source.name)
                shutil.copyfile(source, target)
                print(target)
    return 0


if __name__ == "__main__":
    target = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "docs" / "assets"
    raise SystemExit(main(target))
