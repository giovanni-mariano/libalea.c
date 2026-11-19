#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0
"""Alternate native CPU adapter builds and require identical exported results."""
import argparse
import json
import os
from pathlib import Path
import platform
import statistics
import subprocess

from model_lifecycle_benchmark import sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", type=Path)
    parser.add_argument("workload", type=Path)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=5)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("positive repeats required")
    args.output_directory.mkdir(parents=True, exist_ok=False)
    report = dict(platform=platform.platform(), threads=1, mode="fresh", complete=False,
                  model_sha256=sha256(args.model), workload_sha256=sha256(args.workload),
                  binaries={k: sha256(getattr(args, k)) for k in ("baseline", "candidate")},
                  records=[], summary={}, results_identical=True,
                  note="OS cache/frequency uncontrolled; sequential alternating processes; JSON excluded from query timers")
    path = args.output_directory / "report.json"
    reference = None
    for iteration in range(args.repeats):
        for engine in (("baseline", "candidate") if iteration % 2 == 0 else ("candidate", "baseline")):
            print(engine, iteration+1, flush=True)
            stem = args.output_directory / f"{engine}-{iteration+1}"
            raw, log = stem.with_suffix(".json"), stem.with_suffix(".log")
            with raw.open("w") as out, log.open("w") as err:
                subprocess.run([str(getattr(args, engine).resolve()), "fresh", str(args.model.resolve()),
                                str(args.workload.resolve()), "-"], stdout=out, stderr=err,
                               env=dict(os.environ, ALEA_CACHE="off", OMP_NUM_THREADS="1"),
                               timeout=1800, check=True)
            result = json.loads(raw.read_text())
            answers = (result["points"], result["rays"])
            if reference is None:
                reference = answers
            matches = answers == reference
            report["results_identical"] &= matches
            row = {k: v for k, v in result.items() if k not in ("points", "rays")}
            row.update(build=engine, iteration=iteration+1, raw=str(raw), log=str(log), matches_reference=matches)
            report["records"].append(row)
            path.write_text(json.dumps(report, indent=2)+"\n")
            if not matches:
                raise RuntimeError(f"Result changed: {raw}; comparison stopped")
    for engine in ("baseline", "candidate"):
        rows = [r for r in report["records"] if r["build"] == engine]
        report["summary"][engine] = {
            key: dict(median=statistics.median(r[key] for r in rows),
                      minimum=min(r[key] for r in rows), maximum=max(r[key] for r in rows))
            for key in ("ready_s", "point_s", "ray_s", "peak_rss_bytes")}
    report["complete"] = True
    path.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report["summary"], indent=2))


if __name__ == "__main__":
    main()
