#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0
"""Build, query, export and reload a shielded sphere without external data.

Usage: python3 examples/python/01_geometry.py [output_directory]
Requires pyalea and NumPy. The output directory is created if necessary.
"""
import argparse
from pathlib import Path

import pyalea


def build_model():
    system = pyalea.System()
    fuel = system.add_material(1)
    moderator = system.add_material(2)
    _, outside_inner, inside_inner = system.sphere_surface(1, 0, 0, 0, 5)
    _, outside_outer, inside_outer = system.sphere_surface(2, 0, 0, 0, 10)
    system.add_cell(1, inside_inner, fuel, 10.0)
    shell = system.create_intersection(outside_inner, inside_outer)
    system.add_cell(2, shell, moderator, 1.0)
    system.add_cell(3, outside_outer)  # Omitted material means void.
    system.build_universe_index()
    system.prepare_query_acceleration()
    return system


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_directory", type=Path, nargs="?", default=Path("."))
    args = parser.parse_args()
    args.output_directory.mkdir(parents=True, exist_ok=True)
    system = build_model()
    for point, expected in [((0, 0, 0), (1, 1)), ((7, 0, 0), (2, 2)), ((12, 0, 0), (3, 0))]:
        result = system.find_cell_at(*point)
        assert result == expected, (point, result)
        print(f"{point}: cell {result[0]}, material {result[1]}")
    ray = system.raycast(-15, 0, 0, 1, 0, 0, t_max=30)
    print("Ray segments:", ray["segments"])
    path = args.output_directory / "shielded_sphere.i"
    system.export_mcnp(str(path))
    restored = pyalea.load_mcnp(str(path))
    restored.build_universe_index()
    assert restored.find_cell_at(0, 0, 0) == (1, 1)
    print("Saved and reloaded", path)


if __name__ == "__main__":
    main()
