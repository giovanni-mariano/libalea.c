#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0
"""Check selected segments and points against analytic fills and the full CSG tracer."""
import argparse
import json
import math
import os
import random
from pathlib import Path
import subprocess
import sys
import tempfile


def check_spatial_index(adapter, root):
    """Exercise indexed CSG, aliases, overlaps and the out-of-domain fallback."""
    model = root / "regional.i"
    cards = ["1 0 -101 102 -103", "2 0 (102 -104):(105 -103)", "3 0 -106"]
    surfaces = ["101 so 4", "102 px -3", "103 px 3", "104 px -2", "105 px 2",
                "106 gq 1 1 1 0 0 0 0 0 0 -25", "107 px 3"]
    # A translated, rotated ellipsoid exercises every mixed term in the
    # quadric enclosure. Its complement tests sign folding in another cell.
    cards.extend(["4 0 -108 -109", "5 0 #4 -110"])
    surfaces.extend(["108 gq 2 1.75 3 0.8 -0.4 1.2 -42 2.5 -25.2 234.15",
                     "109 p 0.3 0.4 -0.6 8", "110 so 14"])
    for i in range(80):
        cards.append(f"{10+i} 0 {200+i} -{201+i} -107")
        surfaces.append(f"{200+i} py {i-40}")
    surfaces.append("280 py 40")
    model.write_text("Mixed bounded and unbounded cells\n" + "\n".join(cards) +
                     "\n\n" + "\n".join(surfaces) + "\n\n")
    workload = root / "regional.rays"
    rays = ["p 0 0 0", "r -9 0 0 1 0 0 18", "r 9 0 0 -1 0 0 18",
            "r 3 -15 0 0 1 0 30", "r -15 4 0 1 0 0 30",
            "r -10 -9 -8 1 0.9 0.75 29", "r 4.9999999999 -8 0 0 1 0 16",
            "r -40 2 0 1 0 0 90", "r 2 0 0 0 0 -1 7"]
    rng = random.Random(135231)
    for _ in range(150):
        origin = [rng.uniform(-40, 40) for _ in range(3)]
        direction = [rng.uniform(-1, 1) for _ in range(3)]
        rays.append("p " + " ".join(map(str, origin)))
        rays.append("r " + " ".join(map(str, [*origin, *direction, rng.uniform(1, 100)])))
    for x in (-120, -96, -32, -5, -3, 0, 3, 5, 32, 96, 120):
        rays.append(f"p {x} 0 0")
        rays.append(f"r {x} 0 0 1 0 0 240")
    workload.write_text("\n".join(rays) + "\n")
    command = [adapter.resolve(), "fresh", model, workload, root / "unused"]
    original = json.loads(subprocess.check_output(command))
    indexed = json.loads(subprocess.check_output(command + ["--point-grid", *(["-32", "32"] * 3)]))
    assert indexed["ray_algorithm"] == "parent-interval-regional"
    assert indexed["ray_index_bytes"] > original["ray_index_bytes"]
    for key in ("points", "rays"):
        assert indexed[key] == original[key], key


def check_lattices_and_coincident_cards(adapter, root):
    fixture = Path(__file__).resolve().parents[1] / "tests" / "data"
    alias = root / "coincident.i"
    alias.write_text("""Coincident planes on an inactive union branch
1 0 (5 -6):(1 -6 7 -8)

1 px 0
5 px 0
6 px 2
7 py 5
8 py 6

""")
    workload = root / "boundary.workload"
    rng = random.Random(18472)
    queries = ["r -2 0 0 1 0 0 14", "r 12 0 0 -1 0 0 14",
               "r 3 1 0 1 0 0 7", "r 3 1 0 1 0 0.0000001 7"]
    for _ in range(160):
        x, y, z = rng.uniform(-1, 11), rng.uniform(-2, 2), rng.uniform(-1.2, 1.2)
        queries.extend((f"p {x} {y} {z}", f"r {x} {y} {z} 1 0.04 -0.08 8"))
    workload.write_text("\n".join(queries) + "\n")
    for model in (alias, fixture / "mcnp_nested_lattice_transformed.mcnp",
                  fixture / "mcnp_lattice_repeating.mcnp"):
        cmd = [adapter.resolve(), "fresh", model, workload, root / "unused"]
        reference = json.loads(subprocess.check_output(cmd))
        indexed = json.loads(subprocess.check_output(cmd + ["--point-grid", *(["-8", "8"] * 3)]))
        for key in ("rays", "points"):
            assert indexed[key] == reference[key], (model.name, key)
        if model == alias:
            assert reference["rays"][0][0]["surface_out"] == 1


def check_radial_surfaces(adapter, root):
    model = root / "radial.i"
    surfaces = ["1 tz 0 0 0 3 0.5 1", "2 tx 12 0 0 4 1.5 0.7",
                "3 ty -12 4 3 2 0.7 0.5", "4 c/x 8 3 0.8", "5 c/y -7 4 2",
                "6 c/z 5 -3 1", "7 s 4 5 -5 2"]
    cells = [f"{i} 0 -{i}" for i in range(1, 8)]
    for i in range(25):
        cells.append(f"{100+i} 0 {200+i} -{201+i} 7 1 2 3 4 5 6")
        surfaces.append(f"{200+i} py {i-12}")
    surfaces.append("225 py 13")
    model.write_text("Radial predicates and complements\n" + "\n".join(cells) +
                     "\n\n" + "\n".join(surfaces) + "\n\n")
    queries = []
    for y in (0, 2, 3, 4, math.nextafter(4, math.inf), math.nextafter(4, -math.inf),
              4+1e-10, 4-1e-10, 8, 12):
        for z in (0, 0.5, math.nextafter(0.5, 0), math.nextafter(0.5, 1), 0.5-1e-12, 0.5+1e-12, 3):
            queries.extend((f"r -20 {y} {z} 1 0 0 40", f"r 20 {y} {z} -1 0 0 40",
                            f"p 0 {y} {z}", f"p 12 {y} {z}"))
    rng = random.Random(4006)
    for _ in range(200):
        point = [rng.uniform(-18, 18) for _ in range(3)]
        direction = [rng.uniform(-1, 1) for _ in range(3)]
        queries.append("r " + " ".join(map(str, [*point, *direction, 40])))
        queries.append("p " + " ".join(map(str, point)))
    path = root / "radial.workload"
    path.write_text("\n".join(queries) + "\n")
    cmd = [adapter.resolve(), "write", model, path]
    full_image, regional_image = root / "radial-full.bin", root / "radial-prepared.bin"
    full = json.loads(subprocess.check_output([*cmd, full_image]))
    regional = json.loads(subprocess.check_output(
        [*cmd, regional_image, "--point-grid", *(["-16", "16"] * 3)]))
    assert full_image.read_bytes() == regional_image.read_bytes()
    if sys.platform.startswith("linux"):
        launcher = Path(__file__).with_name("run_interval_cache.py")
        cache = root / "index-cache"
        for options in ([], ["--selective"]):
            restore = [sys.executable, launcher, "--cache", cache, *options, adapter.resolve(),
                       "snapshot", model, path, full_image, "--point-grid", *(["-16", "16"] * 3)]
            for hit in (False, True):
                run = subprocess.run(restore, capture_output=True, check=True)
                assert (b"verified regional image" in run.stderr) == hit
                record = json.loads(run.stdout)
                assert all(record[k] == full[k] for k in ("points", "rays"))
        shifted_domain = restore[:-6] + ["-15", "17", "-16", "16", "-16", "16"]
        shifted = subprocess.run(shifted_domain, capture_output=True, check=True)
        assert b"verified regional image" not in shifted.stderr
        assert all(json.loads(shifted.stdout)[k] == full[k] for k in ("points", "rays"))
        images = list(cache.glob("*/index"))
        assert len(images) == 3
        assert full_image.read_bytes() == regional_image.read_bytes()
        for image in images:
            with image.open("r+b") as f:
                f.truncate(512)
        assert subprocess.run(restore, capture_output=True).returncode != 0
    for key in ("rays", "points"):
        assert full[key] == regional[key], ("radial", key)



def check_point_tile_precedence(adapter, root):
    model = root / "tile-order.i"
    model.write_text("""Earlier partial owners and later solid owners
10 0 -1
20 0 -2
30 0 -3
40 0 -4
50 0 -5

1 s -2 0 0 1
2 s 2 0 0 1
3 px 0
4 py 0
5 so 20

""")
    xs = [-8, -2, -1, math.nextafter(-1, 0), math.nextafter(0, -1), 0,
          math.nextafter(0, 1), math.nextafter(1, 0), 1, 2, 8]
    queries = [f"p {x} {y} {z}" for x in xs for y in [-1, 0, 1] for z in [-1, 0, 1]]
    center_row = [i * 9 + 4 for i in range(len(xs))]
    workload = root / "tile-order.workload"
    workload.write_text("\n".join(queries + ["r -8 0 0 1 0 0 16"]) + "\n")
    cmd = [adapter.resolve(), "fresh", model, workload, root / "unused"]
    full = json.loads(subprocess.check_output(cmd))
    indexed = json.loads(subprocess.check_output(cmd + ["--point-grid", *(["-8", "8"] * 3)]))
    assert [full["points"][i]["cell"] for i in center_row] == [
        30, 10, 10, 10, 30, 30, 40, 20, 20, 20, 40]
    assert indexed["points"] == full["points"] and indexed["rays"] == full["rays"]


def check_adjacent_roots(adapter, root):
    # A locally redundant slanted plane cuts a ray within one representable
    # step of the active x plane. Match the native FP64 reference at the cut.
    model = root / "adjacent.i"
    model.write_text("""Adjacent roots from distinct planes
10 0 1 : (2 3)
20 0 -4
30 0 3

1 px 1000
2 p 1 1 0 1000
3 py 100
4 so 1600

""")
    work = root / "adjacent.workload"
    work.write_text("p 0 0 0\nr -572.5672551759482 -2.842170943040401e-14 0 1 0 0 2000\n")
    cmd = [adapter.resolve(), "fresh", model, work, root / "unused"]
    ref = json.loads(subprocess.check_output(cmd))
    grid = ["--point-grid", "-2400", "2400", "-16", "16", "-16", "16"]
    regional = json.loads(subprocess.check_output(cmd + grid))
    assert len(ref["rays"][0]) == 2 and ref["rays"][0][0]["surface_out"] in (1, 2)
    assert ref["rays"] == regional["rays"]


def check_fill_world_distance(adapter, root):
    # The rounded xy rotation leaves z and dz unchanged. The child's
    # z=10 plane must coincide with the parent's in either direction.
    model = root / "coincident-fill.i"
    deck = """Coincident planes across a rounded fill rotation
1 0 -1 fill=1 (58)
2 0 1
101 0 -2 u=1
102 0 2 -3 u=1
103 1 -1 3 u=1

1 pz 10
2 pz 0
3 pz 10

tr58 0 0 {shift} 0.86603 0.50000 0 0.50000 -0.86603 0 0 0 1
m1 1001 1
"""
    work = root / "coincident-fill.workload"
    work.write_text("p 0 0 1\nr 0 0 -114 4 0 3 500\nr 400 0 186 -4 0 -3 500\n")
    cmd = [adapter.resolve(), "fresh", model, work, root / "unused"]
    # Moving the child down by 2^-43 cm opens a genuine material layer.
    for shift in (0, -2.0**-43):
        model.write_text(deck.format(shift=shift))
        # An exact tie follows the incident cell's card, including on the
        # entry of the following segment. Check this against the public walk.
        env = dict(os.environ, ALEA_INTERVAL_COMPARE_PUBLIC="surfaces") if shift == 0 else None
        full = json.loads(subprocess.check_output(cmd, env=env))["rays"]
        indexed = json.loads(subprocess.check_output(cmd + ["--point-grid", *(["-600", "600"] * 3)]))["rays"]
        assert full == indexed and len(full) == 2
        for row, origin, dz in zip(full, [-114, 186], [0.6, -0.6]):
            cells = [101, 102, 103, 2] if shift else [101, 102, 2]
            surfaces = [2, 3, 1] if shift else [2, 3 if dz > 0 else 1]
            heights = [shift, 10 + shift, 10] if shift else [0, 10]
            if dz < 0:
                cells.reverse()
                surfaces.reverse()
                heights.reverse()
            ends = [0, *((z - origin)/dz for z in heights), 500]
            assert [s["label"]["cell"] for s in row] == cells
            assert [s["label"]["material"] for s in row] == [int(c == 103) for c in cells]
            assert [s["enter"] for s in row] == ends[:-1] and [s["exit"] for s in row] == ends[1:]
            assert [s["surface_in"] for s in row] == [-1, *surfaces]
            assert [s["surface_out"] for s in row] == [*surfaces, -1]


def check_close_quadric_gap(adapter, root):
    # The two roots are ~1.19e-12 cm apart. In the first deck the half-space
    # -3 makes the cut by surface 2 redundant; in the second it bounds cell 10.
    deck = """Two nearly coincident quadrics
10 1 -1 {first_cell}
20 2 -1 -1

1 gq 0.030153689607081 0.969846310392919 1 0.342020143325863 0 0
     -206.814337848850040 -1172.902394277684600 -102.000808796413540
     216968.840459326430000
2 gq 0.030153689607081 0.969846310392919 1 0.342020143325863 0 0
     -206.814337848838900 -1172.902394277686900 -102.000808796413590
     216968.840459326510000
3 px 5000

m1 1001 1
m2 1001 1
"""
    model = root / "close-quadrics.i"
    work = root / "close-quadrics.workload"
    work.write_text("p 0 0 0\nr -262.4563849952908 1157.0080930303711 -403.9923819206924 "
                    "0.8565508818580871 -0.45145463514230855 0.7019354298985057 1000\n"
                    "r 453.77634785815815 779.5096934187713 182.9536383230443 "
                    "-0.8565508818580871 0.45145463514230855 -0.7019354298985057 1000\n")
    cmd = [adapter.resolve(), "fresh", model, work, root / "unused"]
    for expression, expected in [("1 (2 : -3)", [10, 20]), ("2", [10, -1, 20])]:
        model.write_text(deck.format(first_cell=expression))
        for options in [[], ["--point-grid", *(["-2000", "2000"] * 3)]]:
            forward, reverse = json.loads(subprocess.check_output(cmd + options))["rays"]
            assert [s["label"]["cell"] for s in forward] == expected
            assert [s["label"]["cell"] for s in reverse] == expected[::-1]
            assert forward[-1]["surface_in"] == reverse[0]["surface_out"] == 1
            assert abs(forward[-1]["enter"] - 576.1743444620345) < 3e-12
            for row, surfaces in ((forward, (2, 1)), (reverse, (1, 2))):
                assert row[0]["enter"] == 0 and row[-1]["exit"] == 1000
                assert all(a["exit"] == b["enter"] for a, b in zip(row, row[1:]))
                if len(expected) == 3:
                    gap = row[1]
                    assert (gap["surface_in"], gap["surface_out"]) == surfaces
                    assert abs((gap["exit"] - gap["enter"]) - 1.18624e-12) <= math.ulp(gap["enter"])


def check_substep_fill(adapter, root):
    # Rounding the plane numerator before dividing by dx puts these cuts
    # one step apart. Both exact roots round to the same world distance.
    deck = root / "substep.i"
    deck.write_text("""Parallel planes in a fill
1 0 -1 fill=7
11 0 -2 u=7
12 0 3 u=7

1 so 500000000
2 px 5.960464477539063e-9
3 px 8.940696716308593e-9

""")
    work = root / "substep.workload"
    work.write_text("p 0 0 0\nr -102400000.00000001 0 0 3 4 0 200000000\n")
    cmd = [adapter.resolve(), "fresh", deck, work, root / "unused"]
    for options in ([], ["--point-grid", *(["-500000000", "500000000"] * 3)]):
        ray, = json.loads(subprocess.check_output(cmd + options))["rays"]
        assert [s["label"]["cell"] for s in ray] == [11, 12]
        assert ray[0]["exit"] == ray[1]["enter"] == float.fromhex("0x1.4585555555557p+27")


def check_overlap_scoring(adapter, root):
    model = root / "overlap-score.i"
    deck = """Earlier fill placement inside a later material cell
900 0 {region} fill=7
100 2 -2 -3
15 1 -1 -2 -3 u=7

1 px 1
2 px 2
3 so 20
4 px 3
5 px 4

m1 1001 1
m2 1001 1
"""
    work = root / "overlap-score.workload"
    work.write_text("p 1.5 0 0\nr 0 0 0 1 0 0 5\nr 0 0 0 -1 0 0 1\n")
    cmd = [adapter.resolve(), "fresh", model, work, root / "unused"]
    env = dict(os.environ, ALEA_INTERVAL_COMPARE_PUBLIC="overlaps")
    # A second occurrence of the fill is undefined throughout [3,4].
    # It keeps the entire ray in the failure count despite the overlap at [1,2].
    for region, failures in [("1 -2", 0), ("(1 -2):(4 -5)", 1)]:
        model.write_text(deck.format(region=region))
        for options in ([], ["--point-grid", *(["-8", "8"] * 3)]):
            run = subprocess.run(cmd + options, env=env, capture_output=True)
            assert run.returncode == failures, run.stderr
            score = json.loads(run.stdout)["geometry_comparison"]
            assert score == dict(rays=2, matched=1, overlap_only=1-failures, remaining=failures), score


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("adapter", type=Path)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp)
        check_adjacent_roots(args.adapter, root)
        check_fill_world_distance(args.adapter, root)
        check_close_quadric_gap(args.adapter, root)
        check_substep_fill(args.adapter, root)
        check_overlap_scoring(args.adapter, root)
        check_spatial_index(args.adapter, root)
        check_lattices_and_coincident_cards(args.adapter, root)
        check_radial_surfaces(args.adapter, root)
        check_point_tile_precedence(args.adapter, root)
        model = root / "fills.i"
        # The second occurrence places surface 5 at x=1, inside the first
        # parent. Only the first occurrence's copy at x=0 can bound it.
        model.write_text("""Two disjoint parents sharing one universe
1 0 1 -2 fill=1
2 0 3 -4 fill=1 (1 0 0)
101 0 -5 u=1
102 0 5 u=1

1 px -2
2 px 2
3 px 9
4 px 11
5 px 0

""")
        workload = root / "rays"
        workload.write_text("p 0.5 0 0\nr -3 0 0 1 0 0 15\n"
                            "r 12 0 0 -1 0 0 15\nr 0.5 0 0 1 0 0 1\n")
        result = json.loads(subprocess.check_output(
            [args.adapter.resolve(), "fresh", model, workload, root / "unused"]))
        indexed_fill = json.loads(subprocess.check_output(
            [args.adapter.resolve(), "fresh", model, workload, root / "unused",
             "--point-grid", *(["-4", "4"] * 3)]))
        assert all(result[k] == indexed_fill[k] for k in ("points", "rays"))
        # A rounded rotation changes the length of the inverse-transformed
        # unit ray. The local plane x'=1 intersects the world x axis at
        # (a*a + b*b)/a for the matrix [[a,b],[b,-a]] on this card.
        model.write_text("""Rounded transform and a local plane
1 0 -1 fill=1 (58)
101 0 -5 u=1
102 0 5 u=1

1 so 10
5 px 1

tr58 0 0 0 0.86603 0.50000 0 0.50000 -0.86603 0 0 0 1
""")
        a, b = 0.86603, 0.5
        crossing = (a * a + b * b) / a
        workload.write_text(f"p {crossing-1e-7} 0 0\np {crossing+1e-7} 0 0\n"
                            "r 0 0 0 1 0 0 2\nr 2 0 0 -1 0 0 2\n")
        rounded = json.loads(subprocess.check_output(
            [args.adapter.resolve(), "fresh", model, workload, root / "unused"]))
        indexed_rounded = json.loads(subprocess.check_output(
            [args.adapter.resolve(), "fresh", model, workload, root / "unused",
             "--point-grid", *(["-4", "4"] * 3)]))
        assert all(rounded[k] == indexed_rounded[k] for k in ("points", "rays"))
    cells = [-1, 101, 102, -1, 102, -1]
    boundaries = [-1, 1, 5, 2, 3, 4, -1]
    cases = [
        ([0, 1, 3, 5, 12, 14, 15], cells, boundaries),
        ([0, 1, 3, 10, 12, 14, 15], cells[::-1], boundaries[::-1]),
        ([0, 1], [102], [-1, -1]),
    ]
    assert result["ray_algorithm"] == "parent-interval-full-expression"
    for segments, (ts, owners, surfaces) in zip(result["rays"], cases, strict=True):
        expected = list(zip(ts[:-1], ts[1:], owners, surfaces[:-1], surfaces[1:], strict=True))
        actual = [(s["enter"], s["exit"], s["label"]["cell"],
                   s["surface_in"], s["surface_out"]) for s in segments]
        assert actual == expected, (actual, expected)
    assert [p["cell"] for p in rounded["points"]] == [101, 102]
    for segments, owners, t in zip(rounded["rays"], [[101, 102], [102, 101]],
                                   [crossing, 2 - crossing], strict=True):
        assert [s["label"]["cell"] for s in segments] == owners
        assert segments[0]["enter"] == 0 and segments[1]["exit"] == 2
        assert segments[0]["surface_out"] == segments[1]["surface_in"] == 5
        assert segments[0]["exit"] == segments[1]["enter"]
        assert math.isclose(segments[0]["exit"], t, abs_tol=1e-12, rel_tol=0)
    print("Analytic fills, randomized CSG, curved surfaces, coincident cards and nested lattices passed.")


if __name__ == "__main__":
    main()
