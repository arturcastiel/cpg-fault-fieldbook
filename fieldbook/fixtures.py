"""Locate the bundled, source-attributed offline teaching grid."""

from importlib.resources import files
from pathlib import Path


def quickstart() -> Path:
    return Path(files("fieldbook").joinpath("data/quickstart.grdecl"))


def egrid(*, rotated: bool = False) -> Path:
    name = "cart_2x2x2_rotated.EGRID" if rotated else "cart_2x2x2_no_actnum.EGRID"
    return Path(files("fieldbook").joinpath("data", name))
