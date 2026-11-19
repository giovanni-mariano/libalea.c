#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0
"""Reproduce the E-lite prepared-query replays without reading model geometry."""
import argparse
from pathlib import Path
import random
from tempfile import TemporaryDirectory

from model_lifecycle_benchmark import default_workload


def row(kind, values):
    return " ".join([kind, *map(str, values)])


def additional_rays(seed=1568290, count=160):
    rng = random.Random(seed)
    for i in range(count):
        origin = [rng.uniform(-1300, 1300) for _ in range(3)]
        if i % 4 == 0:
            origin[i % 3] = rng.choice([-512, -64, 0, 64, 512])
        direction = ([rng.uniform(-1, 1) for _ in range(3)] if i % 3
                     else [1, 0, 0])
        yield row("r", [*origin, *direction, rng.choice([100, 500, 1000, 2000])])


def historical_suite():
    with TemporaryDirectory() as tmp:
        path = Path(tmp) / "axis.workload"
        default_workload(path, (500, 1500, -500, 500, -500, 500), 256)
        axis = path.read_text().splitlines()
    points = [line for line in axis if line.startswith("p ")] * 2
    rays = [line for line in axis if line.startswith("r ")]
    for i, line in enumerate(rays.copy()):
        direction = [1.0] * 3
        direction[(i + 1) % 3] = 0.17 if i % 2 == 0 else -0.17
        direction[(i + 2) % 3] = 0.11 if i % 4 < 2 else -0.11
        rays.append(row("r", [*line.split()[1:4], *direction, (50, 350, 700, 1000)[i % 4]]))
    rng = random.Random(8742026)
    for i in range(1024):
        coords = ([rng.uniform(350, 1700), rng.uniform(-650, 650), rng.uniform(-650, 650)]
                  if i % 2 else [rng.uniform(-2600, 2600) for _ in range(3)])
        points.append(row("p", coords))
    for _ in range(16):
        origin = [rng.uniform(350, 1700), rng.uniform(-650, 650), rng.uniform(-650, 650)]
        rays.append(row("r", [*origin, *(rng.uniform(-1, 1) for _ in range(3)),
                              rng.choice([50, 350, 700, 1000])]))
    return points + rays


def fresh_points(seed=97812037, count=2048):
    rng = random.Random(seed)
    for i in range(count):
        bounds = [(-1800, 1800)] * 3 if i % 2 else [(400, 1600), (-600, 600), (-600, 600)]
        yield row("p", [rng.uniform(a, b) for a, b in bounds])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workload", choices=("suite", "combined", "fresh-points", "first-use"))
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    if args.workload == "first-use":
        lines = [*fresh_points(415307, 8192), *additional_rays(630743, 640)]
    elif args.workload == "fresh-points":
        lines = [*fresh_points(), *additional_rays()]
    else:
        lines = historical_suite()
        if args.workload == "combined":
            lines += ["p 0 0 0", *additional_rays()]
    args.output.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
