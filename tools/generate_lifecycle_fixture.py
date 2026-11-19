#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0
"""Generate disjoint cylinders with oblique caps and a fixed point/ray workload."""
import argparse
from pathlib import Path


def generate(model, workload, count):
    with model.open("w") as stream:
        stream.write("Oblique capped cylinders for certified mixed-bound replay\n")
        for i in range(count):
            s = 3 * i + 1
            stream.write(f"{i+1} 0 -{s} {s+1} -{s+2} imp:n=1\n")
        stream.write("\n")
        for i in range(count):
            s, x = 3 * i + 1, 4 * i
            stream.write(f"{s} c/z {x} 0 1\n{s+1} p -1 0 1 {-x}\n{s+2} p -1 0 1 {3-x}\n")
        stream.write("\nmode n\n")
    with workload.open("w") as stream:
        for i in range(count):
            x = 4 * i
            for z in (-15, -3, 0.5, 1.5, 2.5, 5, 15):
                for y in (0, 0.5, 2):
                    stream.write(f"p {x} {y} {z}\n")
            stream.write(f"r {x} 0 -20 0 0 1 40\n")
        # Cross the entire assembly above/below all caps: no cell candidates
        # should survive certified z bounds.
        for z in (-15, 15):
            stream.write(f"r -2 0 {z} 1 0 0 {4*count+4}\n")


def generate_quadric(model, workload, count):
    """Rotated cylinders closed by oblique slabs; independent fixed probes."""
    with model.open("w") as stream:
        stream.write("Rotated quadric cylinders with certified oblique slabs\n")
        for i in range(count):
            s = 3 * i + 1
            stream.write(f"{i+1} 0 -{s} {s+1} -{s+2} imp:n=1\n")
        stream.write("\n")
        for i in range(count):
            s, x = 3 * i + 1, 20 * i
            # (X-x-Y)^2 + 2*Z^2 <= 8; x <= X+Y <= x+10.
            stream.write(f"{s} gq 1 1 2 -2 0 0 {-2*x} {2*x} 0 {x*x-8}\n")
            stream.write(f"{s+1} p 1 1 0 {x}\n{s+2} p 1 1 0 {x+10}\n")
        stream.write("\nmode n\n")
    with workload.open("w") as stream:
        for i in range(count):
            x = 20 * i
            for t in (-100, -5, 0, 3, 5, 10, 15, 100):
                for d in (-4, 0, 1, 4):
                    for z in (0, 3):
                        stream.write(f"p {x+(t+d)/2} {(t-d)/2} {z}\n")
            stream.write(f"r {x-20} 2 0 1 0 0 60\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", type=Path)
    parser.add_argument("workload", type=Path)
    parser.add_argument("--cells", type=int, default=128)
    parser.add_argument("--rotated-quadrics", action="store_true")
    args = parser.parse_args()
    if not 1 <= args.cells <= 100000:
        parser.error("cells must be in [1, 100000]")
    (generate_quadric if args.rotated_quadrics else generate)(args.model, args.workload, args.cells)


if __name__ == "__main__":
    main()
