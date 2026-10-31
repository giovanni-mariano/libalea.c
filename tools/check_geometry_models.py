#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0
"""Compare representative geometry workloads. Private inputs are never copied.

Raw logs and summaries default to the Git-ignored build directory. Output hashes
are regression evidence; independent analytic tests remain the accuracy oracle.
"""
import argparse
import json
import math
import os
from pathlib import Path
import statistics
import subprocess


def parse(text):
    rows = []
    for line in text.splitlines():
        kind, _, tail = line.partition(" ")
        if kind not in ("load", "prepare", "tracks", "navigation", "raster"):
            continue
        row = {"kind": kind}
        for token in tail.split():
            key, sep, value = token.partition("=")
            if sep:
                row[key] = value if key.endswith("hash") else float(value)
        rows.append(row)
    return rows


def signature(rows):
    # Timing and allocator variation do not affect correctness comparison.
    return [{k: v for k, v in row.items() if k not in ("ms", "peak_mib")}
            for row in rows if row["kind"] in ("tracks", "navigation", "raster")]


def complete_output(rows, frames):
    """A successful process must still supply every requested measurement."""
    expected = [("load",), ("prepare",), ("navigation",)]
    expected += [("tracks", paths, surfaces, frame)
                 for paths in range(2) for surfaces in range(2)
                 for frame in range(-1, frames)]
    expected += [("raster", frame) for frame in range(-1, frames)]
    actual = []
    for row in rows:
        kind = row["kind"]
        if not all(math.isfinite(row.get(k, math.nan)) for k in ("ms", "peak_mib")):
            return False
        if kind in ("tracks", "navigation", "raster"):
            hashes = ("hash", "surface_hash") if kind == "tracks" else ("hash",)
            if any(len(row.get(k, "")) != 16 or
                   any(c not in "0123456789abcdef" for c in row[k]) for k in hashes):
                return False
        if kind == "tracks":
            actual.append((kind, row.get("paths"), row.get("surfaces"), row.get("frame")))
        elif kind == "raster":
            actual.append((kind, row.get("frame")))
        else:
            actual.append((kind,))
    return len(actual) == len(expected) and set(actual) == set(expected)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", required=True, type=Path)
    parser.add_argument("--candidate", required=True, type=Path)
    parser.add_argument("--private-model", type=Path)
    parser.add_argument("--private-rows", type=int, default=32)
    parser.add_argument("--private-view", type=float, nargs=5,
                        metavar=("XMIN", "XMAX", "YMIN", "YMAX", "Z"),
                        help="locally supplied slice coordinates; required with --private-model")
    parser.add_argument("--output", type=Path, default=Path("build/geometry-model-validation"))
    parser.add_argument("--pairs", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=180)
    parser.add_argument("--expect-navigation-repair", action="store_true",
                        help="accept navigation-only differences when failed baseline histories complete cleanly")
    args = parser.parse_args()
    if args.pairs < 1 or not math.isfinite(args.timeout) or args.timeout <= 0 or not 1 <= args.private_rows <= 4096:
        parser.error("pairs and timeout must be positive; private rows must be 1 to 4096")
    if bool(args.private_model) != bool(args.private_view):
        parser.error("--private-model and --private-view must be supplied together")
    if args.private_view:
        xmin, xmax, ymin, ymax, z = args.private_view
        if not all(map(math.isfinite, args.private_view)) or not (xmin < xmax and ymin < ymax):
            parser.error("private view requires finite coordinates with increasing x/y bounds")
    if args.private_model:
        check = subprocess.run(["git", "check-ignore", "-q",
                                str((args.output / "private-baseline-0.log").resolve())])
        if check.returncode:
            parser.error("private diagnostics require a Git-ignored output directory")
    args.output.mkdir(parents=True, exist_ok=True)
    cases = [
        ("repeating", Path("tests/data/mcnp_lattice_repeating.mcnp"), [1024, 5, -8, 8, -.8, .8, 0], args.pairs),
        ("deep_fills", Path("tests/data/mcnp_deep_translated_fills.mcnp"), [1024, 5, -4, 6, -1, 1, 0], args.pairs),
        ("transformed_lattice", Path("tests/data/mcnp_nested_lattice_transformed.mcnp"), [1024, 5, 3, 9, -.8, .8, 0], args.pairs),
        ("transformed_navigation", Path("tests/data/mcnp_nested_lattice_transformed.mcnp"), [16, 2, 3, 9, -.8, .8, 0], args.pairs),
    ]
    if args.private_model:
        cases.append(("private", args.private_model, [args.private_rows, 3, *args.private_view, 4, 16], 1))
    results = {}
    regression = False
    unresolved = False
    for name, model, config, pairs in cases:
        measurements = {"baseline": [], "candidate": []}
        complete = True
        for pair in range(pairs):
            order = ("baseline", "candidate") if pair % 2 == 0 else ("candidate", "baseline")
            for label in order:
                command = [str(getattr(args, label).resolve()), str(model), *map(str, config)]
                print(f"Running {name} {label} pair {pair + 1}", flush=True)
                try:
                    run = subprocess.run(command, text=True, capture_output=True,
                                         env={**os.environ, "ALEA_NUM_THREADS": "1"},
                                         timeout=args.timeout)
                    stdout, stderr, code = run.stdout, run.stderr, run.returncode
                except subprocess.TimeoutExpired as error:
                    stdout = error.stdout or b""
                    stderr = error.stderr or b""
                    stdout = stdout.decode() if isinstance(stdout, bytes) else stdout
                    stderr = stderr.decode() if isinstance(stderr, bytes) else stderr
                    code = "timeout"
                (args.output / f"{name}-{label}-{pair}.log").write_text(stdout + "\n" + stderr)
                rows = parse(stdout)
                measurements[label].append({"returncode": code, "rows": rows})
                if code != 0 or not complete_output(rows, config[1]):
                    complete = False
        a, b = measurements["baseline"], measurements["candidate"]
        equal = complete and all(signature(x["rows"]) == signature(y["rows"])
                                 for x, y in zip(a, b))
        summary = {"complete": complete, "outputs_equal": equal, "timings": {}}
        for label in ("baseline", "candidate"):
            grouped = {}
            for run in measurements[label]:
                for row in run["rows"]:
                    if row.get("frame", 0) < 0:
                        continue
                    key = row["kind"]
                    if key == "tracks":
                        key += f"_paths{int(row['paths'])}_surfaces{int(row['surfaces'])}"
                    grouped.setdefault(key, []).append(row["ms"])
            summary["timings"][label] = {k: {"median_ms": statistics.median(v),
                                                    "min_ms": min(v), "max_ms": max(v),
                                                    "samples": len(v)} for k, v in grouped.items()}
        summary["navigation"] = {label: [r for run in measurements[label]
                                                for r in run["rows"] if r["kind"] == "navigation"]
                                 for label in ("baseline", "candidate")}
        def clean(rows):
            return bool(rows) and all(row["sources"] > 0 and row["events"] > 0
                                      and row["failures"] == 0 and row["ambiguous"] == 0
                                      for row in rows)

        summary["baseline_navigation_clean"] = complete and clean(summary["navigation"]["baseline"])
        summary["navigation_clean"] = complete and clean(summary["navigation"]["candidate"])
        # A correctness fix intentionally changes the event digest. Require
        # explicit opt-in, unchanged non-navigation results, and completion of
        # the same source set in every pair. Never accept a timeout as repair.
        repaired = complete and args.expect_navigation_repair and all(
            [r for r in signature(x["rows"]) if r["kind"] != "navigation"] ==
            [r for r in signature(y["rows"]) if r["kind"] != "navigation"]
            and len(nx := [r for r in x["rows"] if r["kind"] == "navigation"]) == 1
            and len(ny := [r for r in y["rows"] if r["kind"] == "navigation"]) == 1
            and nx[0]["failures"] + nx[0]["ambiguous"] > 0 and clean(ny)
            and nx[0]["sources"] == ny[0]["sources"]
            and ny[0]["events"] >= nx[0]["events"]
            for x, y in zip(a, b))
        summary["navigation_repaired"] = repaired
        regression |= not (equal or repaired)
        unresolved |= not summary["navigation_clean"]
        summary["peak_mib"] = {label: max((row["peak_mib"] for run in measurements[label]
                                             for row in run["rows"]), default=0)
                               for label in ("baseline", "candidate")}
        results[name] = summary
        (args.output / "summary.json").write_text(json.dumps(results, indent=2) + "\n")
        print(f"{name}: complete={complete}, outputs_equal={equal}, "
              f"navigation_clean={summary['navigation_clean']}, "
              f"navigation_repaired={repaired}", flush=True)
    return int(regression or unresolved)


if __name__ == "__main__":
    raise SystemExit(main())
