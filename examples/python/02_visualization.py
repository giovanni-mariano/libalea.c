#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0
"""Slice and render an inline model; save PPM images using only NumPy/pyalea.

Usage: python3 examples/python/02_visualization.py [output_directory]
PPM is a simple RGB image format; no plotting or image packages are needed.
"""
import argparse
from pathlib import Path

import numpy as np
import pyalea


def write_ppm(path, rgb):
    height, width, channels = rgb.shape
    assert channels == 3
    with path.open("wb") as stream:
        stream.write(f"P6\n{width} {height}\n255\n".encode("ascii"))
        stream.write(np.asarray(rgb, dtype=np.uint8).tobytes())
    print("Saved", path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_directory", type=Path, nargs="?", default=Path("."))
    args = parser.parse_args()
    args.output_directory.mkdir(parents=True, exist_ok=True)
    system = pyalea.load_mcnp_string("""Concentric spheres
1 1 -10 -1
2 2 -1 1 -2
3 0 2

1 so 5
2 so 10

""")
    system.build_universe_index()
    system.prepare_query_acceleration()
    width, height = 160, 120
    grid = system.find_cells_grid_z(0, -15, 15, -11.25, 11.25, width, height)
    ids = np.asarray(grid["material_ids"]).reshape(height, width)
    palette = np.array([[245, 245, 245], [220, 100, 70], [80, 150, 210]], dtype=np.uint8)
    # Grid rows increase in Y; image rows increase downwards.
    image = palette[np.clip(ids, 0, 2)][::-1]
    assert set(np.unique(ids)) == {0, 1, 2}
    write_ppm(args.output_directory / "slice.ppm", image)
    for name, clips in [("render", []), ("cutaway", [(0, 0, -1, 0)])]:
        result = system.render_3d(
            width=width, height=height, eye=(25, 25, 15), target=(0, 0, 0),
            color_by="material", shadows=True, edges=True, clips=clips,
        )
        assert result["rgb"].shape == (height, width, 3)
        write_ppm(args.output_directory / f"{name}.ppm", result["rgb"])


if __name__ == "__main__":
    main()
