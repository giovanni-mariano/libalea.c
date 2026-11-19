#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0
"""Run isolated lifecycle scenarios and gate exact point/ray replay results."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import statistics
import subprocess
import tempfile
import time

STAGE = re.compile(r"\[alea-load-profile\]\s+(\S+)\s+([0-9.eE+-]+) s")
METRICS = ("ready_s", "load_s", "prepare_s", "bounds_s", "previous_load_s",
           "first_query_s", "point_s", "ray_s", "peak_rss_bytes",
           "points_per_s", "rays_per_s", "point_candidates", "ray_point_candidates")


def sha256(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def default_workload(path, bounds, count):
    """Deterministic samples from a caller-fixed box, independent of model bounds."""
    import random
    rng = random.Random(20261117)
    lo, hi = bounds[::2], bounds[1::2]
    with path.open("w") as stream:
        for _ in range(count):
            p = [rng.uniform(a, b) for a, b in zip(lo, hi)]
            stream.write("p " + " ".join(format(v, ".17g") for v in p) + "\n")
        for i in range(max(1, count // 16)):
            axis = i % 3
            p = [rng.uniform(a, b) for a, b in zip(lo, hi)]
            p[axis] = lo[axis]
            direction = [float(k == axis) for k in range(3)]
            values = p + direction + [hi[axis] - lo[axis]]
            stream.write("r " + " ".join(format(v, ".17g") for v in values) + "\n")


def digest(record):
    return record["point_digest"], record["ray_digest"]


def replay_matches(first, second):
    if digest(first) != digest(second):
        return False
    if "acceleration_digest" in first and "acceleration_digest" in second:
        return first["acceleration_digest"] == second["acceleration_digest"]
    return True  # Older reports contain only query replay; record this below.


def summarize(records):
    result = {}
    for scenario in sorted({r["scenario"] for r in records}):
        rows = [r for r in records if r["scenario"] == scenario]
        result[scenario] = {}
        for key in METRICS:
            values = [r[key] for r in rows if r.get(key) is not None]
            if values:
                result[scenario][key] = dict(median=statistics.median(values),
                                             minimum=min(values), maximum=max(values))
        for key in sorted({k for r in rows for k in r["load_stages_s"]}):
            values = [r["load_stages_s"][key] for r in rows if key in r["load_stages_s"]]
            result[scenario]["stage:" + key] = dict(median=statistics.median(values),
                                                    minimum=min(values), maximum=max(values))
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", type=Path)
    parser.add_argument("--binary", type=Path, default=Path("bin/model_lifecycle"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workload", type=Path)
    parser.add_argument("--bounds", type=float, nargs=6, default=(-32, 32, -32, 32, -32, 32))
    parser.add_argument("--points", type=int, default=4096)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--edited-model", type=Path)
    parser.add_argument("--mixed-bounds", action="store_true")
    parser.add_argument("--candidate-profile", action="store_true",
                        help="attribute BLAS point candidates to cells; adds profiling time/memory")
    parser.add_argument("--ray-profile", action="store_true",
                        help="attribute ray costs; requires an ALEA_RAY_PROFILE build and adds timing overhead")
    parser.add_argument("--ray-api", choices=("global", "hierarchical"), default="hierarchical")
    parser.add_argument("--snapshot-io", choices=("buffered", "mapped", "both"), default="buffered")
    parser.add_argument("--baseline", type=Path)
    parser.add_argument("--max-ready-ratio", type=float)
    parser.add_argument("--max-memory-ratio", type=float)
    parser.add_argument("--build-label", default="unspecified; use an isolated release build")
    parser.add_argument("--timeout", type=float, default=1800)
    args = parser.parse_args(argv)
    import math
    if (args.repeats < 1 or args.points < 1 or args.points > 900000 or
            any(not math.isfinite(v) for v in args.bounds) or
            any(a >= b for a, b in zip(args.bounds[::2], args.bounds[1::2]))):
        parser.error("positive repeats/points and finite ordered bounds required")
    if any(v is not None and (not math.isfinite(v) or v <= 0)
           for v in (args.max_ready_ratio, args.max_memory_ratio)):
        parser.error("performance ratios must be finite and positive")
    if not math.isfinite(args.timeout) or args.timeout <= 0:
        parser.error("timeout must be finite and positive")
    if (args.max_ready_ratio is not None or args.max_memory_ratio is not None) and not args.baseline:
        parser.error("performance gates require --baseline")
    inputs = {p.resolve() for p in (args.model, args.binary, args.workload,
                                    args.edited_model, args.baseline) if p is not None}
    outputs = {args.output.resolve()}
    generated_workload = args.output.with_name(args.output.name + ".workload")
    if args.workload is None:
        outputs.add(generated_workload.resolve())
    if inputs & outputs:
        parser.error("report/generated workload must not overwrite an input or baseline")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    report = dict(schema=2, model=str(args.model.resolve()), source_sha256=sha256(args.model),
                  binary_sha256=sha256(args.binary), build_label=args.build_label,
                  platform=platform.platform(), machine=platform.machine(),
                  cpu_count=os.cpu_count(), query_threads=1, repeats=args.repeats,
                  mixed_bounds=args.mixed_bounds, ray_api=args.ray_api,
                  candidate_profile=args.candidate_profile,
                  ray_profile=args.ray_profile,
                  snapshot_io=args.snapshot_io,
                  environment={k: v for k, v in os.environ.items()
                               if k.startswith(("ALEA_", "OMP_", "TINYPAR_"))},
                  page_cache="uncontrolled; warm cache is primed, not an OS-cold measurement",
                  records=[], summary={}, complete=False)
    if args.edited_model:
        report["edited_source_sha256"] = sha256(args.edited_model)
        report["edited_model"] = str(args.edited_model.resolve())
    def save():
        args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    save()
    # Place potentially large snapshots on the output filesystem, not /tmp tmpfs.
    with tempfile.TemporaryDirectory(prefix="lifecycle-", dir=args.output.parent) as temp:
        temp = Path(temp)
        workload = args.workload
        if workload is None:
            workload = generated_workload
            default_workload(workload, args.bounds, args.points)
        report["workload_sha256"] = sha256(workload)
        report["workload"] = str(workload.resolve())
        snapshot, cache = temp / "model.abin", temp / "cache.abin"
        env = dict(os.environ, ALEA_PROFILE_LOAD="1", ALEA_CACHE="on",
                   ALEA_SNAPSHOT_IO="buffered" if args.snapshot_io == "both" else args.snapshot_io)
        env.pop("ALEA_CACHE_DIR", None)
        env.pop("ALEA_CACHE_FILE", None)
        report["effective_environment"] = {k: v for k, v in env.items()
                                            if k.startswith(("ALEA_", "OMP_", "TINYPAR_"))}
        def run(scenario, mode, model, aux, mixed=False, snapshot_io=None):
            command = [str(args.binary.resolve()), mode, str(model.resolve()),
                       str(workload.resolve()), str(aux)]
            if mixed:
                command.append("--mixed-bounds")
            if args.ray_api == "hierarchical":
                command.append("--hierarchical-rays")
            if args.candidate_profile:
                command.append("--candidate-profile")
            if args.ray_profile:
                command.append("--ray-profile")
            print(f"{scenario}: {model}", flush=True)
            start = time.monotonic()
            log = args.output.with_suffix(f".{len(report['records']):03d}-{scenario}.log")
            if log.resolve() in inputs | outputs:
                raise RuntimeError(f"log would overwrite an input/report: {log}")
            run_env = dict(env)
            if snapshot_io:
                run_env["ALEA_SNAPSHOT_IO"] = snapshot_io
            with log.open("w") as errors:
                process = subprocess.run(command, text=True, stdout=subprocess.PIPE,
                                         stderr=errors, env=run_env, timeout=args.timeout)
            if process.returncode:
                raise RuntimeError(f"{scenario} exited {process.returncode}; see {log}")
            record = json.loads(process.stdout)
            record["points_per_s"] = record["points"] / record["point_s"] if record["point_s"] > 0 else None
            record["rays_per_s"] = record["rays"] / record["ray_s"] if record["ray_s"] > 0 else None
            record.update(scenario=scenario, wall_s=time.monotonic() - start,
                          log=str(log), command=command, load_stages_s={},
                          snapshot_io_requested=run_env["ALEA_SNAPSHOT_IO"],
                          snapshot_io_used=re.findall(r"\[alea-snapshot-io\] (\w+)",log.read_text()))
            for name, elapsed in STAGE.findall(log.read_text()):
                stages = record["load_stages_s"]
                stages[name] = stages.get(name, 0) + float(elapsed)
            report["records"].append(record)
            report["summary"] = summarize(report["records"])
            save()
            return record
        try:
            reference = run("snapshot-write", "write", args.model, snapshot)
            prime = run("cache-prime", "cache", args.model, cache)
            report["snapshot_bytes"] = snapshot.stat().st_size
            report["cache_bytes"] = cache.stat().st_size if cache.exists() else None
            if not replay_matches(reference, prime) or prime["cache_hit"]:
                raise RuntimeError("cache priming replay mismatch or unexpected hit")
            edited_reference = None
            if args.edited_model:
                edited_reference = run("edited-reference", "fresh", args.edited_model, "-")
            scenarios = [("fresh", "fresh", args.model, "-", False),
                         ("warm-cache", "cache", args.model, cache, False),
                         ("snapshot", "snapshot", args.model, snapshot, False)]
            if args.snapshot_io == "both":
                scenarios += [("warm-cache-mapped", "cache", args.model, cache, False, "mapped"),
                              ("snapshot-mapped", "snapshot", args.model, snapshot, False, "mapped")]
            if args.mixed_bounds:
                scenarios.append(("mixed-bounds", "fresh", args.model, "-", True))
            if args.edited_model:
                scenarios += [("edited-fresh", "fresh", args.edited_model, "-", False),
                              ("edit-with-previous", "edit", args.edited_model, snapshot, False)]
            for repeat in range(args.repeats):
                for scenario in scenarios if repeat % 2 == 0 else reversed(scenarios):
                    record = run(*scenario)
                    wanted = edited_reference if scenario[0].startswith(("edit-", "edited-")) else reference
                    if not replay_matches(record, wanted):
                        raise RuntimeError(f"replay mismatch: {scenario[0]}")
                    if scenario[0].startswith("warm-cache") and not record["cache_hit"]:
                        raise RuntimeError("warm-cache scenario missed the cache")
            if args.baseline:
                old = json.loads(args.baseline.read_text())
                if not old.get("complete"):
                    raise RuntimeError("baseline is incomplete")
                for key in ("source_sha256", "edited_source_sha256", "workload_sha256", "mixed_bounds", "ray_api"):
                    if old.get(key) != report.get(key):
                        raise RuntimeError(f"baseline input differs: {key}")
                old_records = {r["scenario"]: r for r in old["records"]}
                report["baseline_acceleration_checked"] = all(
                    "acceleration_digest" in r for r in old["records"] + report["records"])
                for r in report["records"]:
                    if r["scenario"] not in old_records or not replay_matches(r, old_records[r["scenario"]]):
                        raise RuntimeError(f"baseline replay mismatch: {r['scenario']}")
                for scenario, metrics in report["summary"].items():
                    for metric, ratio in (("ready_s", args.max_ready_ratio),
                                          ("peak_rss_bytes", args.max_memory_ratio)):
                        if ratio is not None:
                            if metric not in metrics or metric not in old["summary"].get(scenario, {}):
                                raise RuntimeError(f"performance metric unavailable: {scenario} {metric}")
                            base = old["summary"][scenario][metric]["median"]
                            if metrics[metric]["median"] > base * ratio:
                                raise RuntimeError(f"performance gate failed: {scenario} {metric}")
            report["acceleration_checked"] = all("acceleration_digest" in r for r in report["records"])
            report["complete"] = True
        except (RuntimeError, subprocess.TimeoutExpired, ValueError, OSError) as error:
            report["error"] = str(error)
            save()
            raise SystemExit(str(error)) from error
        save()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
