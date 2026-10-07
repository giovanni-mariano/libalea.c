#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0
"""Check plot sampling methods and explicit diagnostic overlays.

Run after make tools: python3 tools/check_mc_plotter.py
"""
import argparse
import os
from pathlib import Path
import subprocess
import tempfile


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plotter", type=Path, default=root / "bin/mc_plotter")
    args = parser.parse_args()
    executable = str(args.plotter.resolve())
    env = {key: value for key, value in os.environ.items()
           if not key.startswith(("ALEA_PLOT_", "ALEA_GRID_"))}
    env["ALEA_NUM_THREADS"] = "2"
    with tempfile.TemporaryDirectory(prefix="alea-plot-check-") as directory:
        tmp = Path(directory)
        model = tmp / "spheres.i"
        model.write_text("Spheres\n1 1 -1 -1\n2 2 -1 1 -2\n3 0 2\n\n1 so 2\n2 so 4\n\n")

        def plot(name, options=(), source=model, environment=env):
            output = tmp / f"{name}.bmp"
            result = subprocess.run(
                [executable, str(source), "Z", "0", "-5", "5", "-5", "5", "96",
                 str(output), *options], env=environment, capture_output=True,
                text=True, timeout=30, check=True)
            return output.read_bytes(), result.stdout

        grid, log = plot("grid")
        assert "Grid query:" in log and "Grid coverage query:" not in log
        ray, log = plot("ray", ["--method=ray"])
        assert "Ray raster query:" in log
        assert grid == ray, "Valid sphere geometry changed between sampling methods"
        stats, log = plot("stats", environment=dict(env, ALEA_PLOT_ERROR_STATS="1"))
        assert "Grid coverage query:" in log and stats == grid

        overlap = root / "tests/data/mcnp_overlap.mcnp"
        plain, _ = plot("plain-overlap", source=overlap)
        diagnostics, log = plot("diagnostics", ["--errors"], source=overlap)
        assert "Grid coverage query:" in log and diagnostics != plain
        ray_errors, log = plot("ray-errors", ["--method=ray", "--errors"], source=overlap)
        assert "Diagnostics use grid coverage" in log and diagnostics == ray_errors

        batch = tmp / "batch.txt"
        batch_output = tmp / "batch.bmp"
        batch.write_text(f"Z 0 -5 5 -5 5 96 {batch_output} method=ray\n")
        subprocess.run([executable, str(model), f"--batch={batch}"],
                       env=env, capture_output=True, check=True, timeout=30)
        assert batch_output.read_bytes() == ray
        bad = subprocess.run([executable, str(model), "Z", "0", "-5", "5", "-5", "5", "96",
                              "--method=unknown"], env=env, capture_output=True, timeout=30)
        assert bad.returncode != 0
    print("Plotter grid/ray equivalence, diagnostics, batch, and method validation: PASS")


if __name__ == "__main__":
    main()
