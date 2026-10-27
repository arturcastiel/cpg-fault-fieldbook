"""Regenerate the Z-slide and gallery figures (06-16) from the lessons, in one command.

usage: python tools/figures/make_figures.py [OUTDIR] [--only 07 10 ...]

Runs each lesson with --out OUTDIR/<figure> to write its decks, then draws OUTDIR/<figure>.png
and OUTDIR/<figure>_3d.html (interactive). OUTDIR defaults to output/figures. Needs the
figures extra: pip install -e ".[figures]".
"""

import argparse
import ast
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent

# figure stem, lesson that writes its decks (None: the script builds them), drawing script
FIGURES = [
    ("06_two_sequential_faults", None, "sequential_faults.py"),
    ("07_partial_fault", "20_z_slide_partial_fault", "partial_fault.py"),
    ("08_oblique_fault", "21_z_slide_oblique_fault", "oblique_fault.py"),
    ("09_throw_control", "22_z_slide_throw_control", "throw_control.py"),
    ("10_structural_styles", "23_gallery_structural_styles", "gallery.py"),
    ("11_fault_interactions", "24_gallery_fault_interactions", "gallery.py"),
    ("12_geometry_inputs", "25_gallery_geometry_inputs", "gallery.py"),
    ("13_throw_definitions", "26_gallery_throw_definitions", "gallery.py"),
    ("14_reservoir_applications", "27_gallery_reservoir_applications", "gallery.py"),
    ("15_regridding", "28_gallery_regridding", "gallery.py"),
    ("16_export_workflows", "29_gallery_export_workflows", "gallery.py"),
]


def gallery_title(lesson):
    """'Gallery 3: from fault data ...' -> 'Gallery 3 (lesson 25): from fault data ...'."""
    doc = ast.get_docstring(ast.parse(lesson.read_text(encoding="utf-8")))
    head, rest = doc.splitlines()[0].split(":", 1)
    return f"{head} (lesson {lesson.name[:2]}):{rest}"


def run(*cmd):
    print("$", " ".join(str(c) for c in cmd), flush=True)
    subprocess.run([sys.executable, *map(str, cmd)], check=True, cwd=ROOT)


def main(out, only):
    out.mkdir(parents=True, exist_ok=True)
    for stem, lesson, script in FIGURES:
        if only and stem[:2] not in only:
            continue
        decks = out / stem
        if lesson is not None:
            lesson_file = ROOT / "examples" / f"{lesson}.py"
            run(lesson_file, "--out", decks)
        if script == "gallery.py":
            run(HERE / script, decks, lesson_file, out / stem, gallery_title(lesson_file))
        else:
            run(HERE / script, out)
    print(f"figures in {out}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("out", nargs="?", type=Path, default=ROOT / "output" / "figures")
    parser.add_argument("--only", nargs="*", default=[], help="figure numbers, e.g. 07 10")
    args = parser.parse_args()
    unknown = set(args.only) - {stem[:2] for stem, _, _ in FIGURES}
    if unknown or any(not re.fullmatch(r"\d\d", n) for n in args.only):
        parser.error(f"--only takes figure numbers 06-16, got {sorted(unknown)}")
    main(args.out.resolve(), set(args.only))
