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

        # Read uncompressed 24-bit BMP pixels in file row order. Label checks
        # compare with the unlabeled image, so they also cover the grid Y flip.
        def bitmap(data):
            import struct
            offset = struct.unpack_from("<I", data, 10)[0]
            w, h = struct.unpack_from("<ii", data, 18)
            assert struct.unpack_from("<H", data, 28)[0] == 24
            stride = (w * 3 + 3) & ~3
            return [[tuple(data[offset + y * stride + x * 3:
                                offset + y * stride + x * 3 + 3])
                     for x in range(w)] for y in range(abs(h))]

        baseline, _ = plot("label-baseline", ["--color=materials"])
        base = bitmap(baseline)
        for mode in ("cells", "materials", "cells,materials"):
            rendered, _ = plot("label-" + mode,
                               ["--color=materials", "--labels=" + mode])
            rows = bitmap(rendered)
            changed = 0
            for y in range(96):
                for x in range(96):
                    if rows[y][x] == base[y][x]:
                        continue
                    changed += 1
                    assert 0 < x < 95 and 0 < y < 95
                    # Every ink pixel retains a padding pixel before a contour.
                    assert all(base[ny][nx] == base[y][x]
                               for ny in range(y - 1, y + 2)
                               for nx in range(x - 1, x + 2)), mode
            assert changed > 0, mode

        # A thin annulus is large enough to label by area, but cannot contain
        # the text footprint. It must remain untouched in every label mode.
        thin = tmp / "thin.i"
        thin.write_text(model.read_text().replace("1 so 2", "1 so 3.7"))
        thin_base, _ = plot("thin-base", ["--color=materials"], source=thin)
        thin_rows = bitmap(thin_base)
        ring_color = base[48][76]
        assert sum(color == ring_color for row in thin_rows for color in row) > 100
        for mode in ("cells", "materials", "cells,materials"):
            rendered, _ = plot("thin-" + mode,
                               ["--color=materials", "--labels=" + mode], source=thin)
            rows = bitmap(rendered)
            assert any(rows[y][x] != thin_rows[y][x] for y in range(96) for x in range(96))
            assert all(rows[y][x] == thin_rows[y][x]
                       for y in range(96) for x in range(96)
                       if thin_rows[y][x] == ring_color), mode

        high_ids = tmp / "high-ids.i"
        high_ids.write_text(model.read_text().replace("1 1 -1 -1", "10001 1 -1 -1")
                            .replace("2 2 -1 1 -2", "20002 2 -1 1 -2"))
        regular, _ = plot("regular-materials", ["--color=materials", "--labels=materials"])
        large, _ = plot("large-materials", ["--color=materials", "--labels=materials"], source=high_ids)
        assert regular == large, "Material labels depend on cell ID magnitude"

        # Adjacent cells with the same material must produce one material label.
        same_material = tmp / "same-material.i"
        same_material.write_text(model.read_text().replace("2 2 -1", "2 1 -1"))
        merged = tmp / "merged.i"
        merged.write_text("Merged sphere\n1 1 -1 -2\n3 0 2\n\n2 so 4\n\n")
        split, _ = plot("split-material", ["--color=materials", "--labels=materials"], source=same_material)
        whole, _ = plot("whole-material", ["--color=materials", "--labels=materials"], source=merged)
        def material_ink(data):
            return {(x, y) for y, row in enumerate(bitmap(data))
                    for x, color in enumerate(row) if color == (100, 200, 255)}
        assert material_ink(split) == material_ink(whole)

        batch = tmp / "batch.txt"
        batch_output = tmp / "batch.bmp"
        batch.write_text(f"Z 0 -5 5 -5 5 96 {batch_output} method=ray\n")
        subprocess.run([executable, str(model), f"--batch={batch}"],
                       env=env, capture_output=True, check=True, timeout=30)
        assert batch_output.read_bytes() == ray
        bad = subprocess.run([executable, str(model), "Z", "0", "-5", "5", "-5", "5", "96",
                              "--method=unknown"], env=env, capture_output=True, timeout=30)
        assert bad.returncode != 0
    print("Plotter grid/ray equivalence, diagnostics, labels, batch, and method validation: PASS")


if __name__ == "__main__":
    main()
