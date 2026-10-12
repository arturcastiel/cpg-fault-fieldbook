"""Minimal SVG lattice illustration; numerical gates remain authoritative."""

from html import escape
from pathlib import Path

import numpy as np
from cpg2unstructured import pillar_lattice


def render_pair(before_grid, after_grid, trace: np.ndarray, out_dir: Path) -> tuple[Path, Path]:
    """Draw same-scale top-plane grids and authored trace to separate SVG files."""
    before, after = pillar_lattice(before_grid).top, pillar_lattice(after_grid).top
    points = np.vstack((before.reshape(-1, 2), after.reshape(-1, 2), trace))
    low, high = points.min(axis=0), points.max(axis=0)
    span = max(float(np.max(high - low)), 1e-9)

    def point(p):
        x = 32 + 336 * (float(p[0]) - low[0]) / span
        y = 368 - 336 * (float(p[1]) - low[1]) / span
        return f"{x:.2f},{y:.2f}"

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for label, lattice, name in (
        ("Before", before, "before.svg"), ("After", after, "after.svg")
    ):
        lines = [
            '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 400 400" '
            'role="img" aria-label="Top-plane grid and authored fault trace">',
            '<rect width="400" height="400" fill="#f9fbff"/>',
            f'<text x="18" y="22" font-size="16" fill="#20324a">{escape(label)}</text>',
        ]
        for row in lattice:
            lines.append(f'<polyline fill="none" stroke="#52708e" stroke-width="1.4" '
                         f'points="{" ".join(point(p) for p in row)}"/>')
        for row in lattice.transpose(1, 0, 2):
            lines.append(f'<polyline fill="none" stroke="#52708e" stroke-width="1.4" '
                         f'points="{" ".join(point(p) for p in row)}"/>')
        lines.append(f'<polyline fill="none" stroke="#c12736" stroke-width="3" '
                     f'points="{" ".join(point(p) for p in trace)}"/>')
        lines.append("</svg>")
        destination = out_dir / name
        destination.write_text("\n".join(lines) + "\n", encoding="utf-8")
        paths.append(destination)
    return paths[0], paths[1]
