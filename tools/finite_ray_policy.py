#!/usr/bin/env python3
# SPDX-License-Identifier: MPL-2.0
"""Replay finite geometry calls on a fresh model and check their ordered segments."""
import argparse
import json
import math
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
import pyalea

try:
    import resource
except ImportError:  # Windows does not expose getrusage.
    resource = None


def read_workload(args):
    inputs = []
    for line in args.workload.read_text().splitlines():
        fields = line.split()
        if not fields or fields[0].startswith("#"):
            continue
        if len(fields) != 8 or fields[0] != "r":
            raise ValueError("expected finite rays with a shared distance limit")
        inputs.append([float(value) for value in fields[1:]])
    rows = np.array(inputs)
    reference = json.loads(args.reference.read_text())["rays"]
    if (not inputs or len(rows) != len(reference) or not np.isfinite(rows).all()
            or rows[0, 6] <= 0 or not (rows[:, 6] == rows[0, 6]).all()):
        raise ValueError("workload and reference must cover the same finite batch")
    indices = list(range(len(rows)))
    if args.batch_size is None:
        indices *= args.repeat
    if args.batch_size is None and args.limit is not None:
        indices = indices[:args.limit]
    # Sorting uses known reference lengths only to create adversarial input orders.
    if args.order != "input":
        indices.sort(key=lambda i: len(reference[i]), reverse=args.order == "dense-first")
    if args.batch_size is not None:
        # Consecutive calls consume the ordered stream one pass at a time.
        indices *= args.repeat
        if args.limit is not None:
            indices = indices[:args.limit]
    rays = rows[indices]
    origins = np.ascontiguousarray(rays[:, :3])
    directions = np.ascontiguousarray(rays[:, 3:6])
    return origins, directions, rows[0, 6], indices, reference


def read_calibration(args, distance):
    calibration = None
    if args.calibration:
        calibration = json.loads(args.calibration.read_text())
        if (calibration["mode"] != "regional" or calibration["model"] != str(args.model.resolve())
                or calibration["distance"] != distance or calibration["cold"] != args.cold
                or calibration["rays"] <= 0
                or not all(math.isfinite(calibration[k]) and calibration[k] >= 0
                           for k in ("build_s", "tail_s"))):
            raise ValueError("expected matching explicit regional calibration")

    return calibration


def verify_rows(result, indices, reference):
    max_delta = 0.
    segments = 0
    for row, (ray, i) in enumerate(zip(result, indices, strict=True)):
        if len(ray) != len(reference[i]):
            raise ValueError(f"row {row}, input {i}: segment count differs")
        for got, expected in zip(ray, reference[i], strict=True):
            lo, hi, *owner = expected
            delta = max(abs(got["t_enter"] - lo), abs(got["t_exit"] - hi))
            if (not all(math.isfinite(t) for t in (lo, hi, got["t_enter"], got["t_exit"]))
                    or delta > 1e-6 or owner != [got[k] for k in
                    ("cell_id", "material_id", "density", "resolution_flags")]):
                raise ValueError(f"row {row}, input {i}: {got} disagrees with {expected}")
            max_delta = max(max_delta, delta)
            segments += 1
    return segments, max_delta


def peak_bytes():
    if resource is None:
        return None
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return peak if sys.platform == "darwin" else peak * 1024


def replay_calls(args, model, trace, indices, reference, distance, load_s, setup_s):
    calls = []
    cumulative = load_s + setup_s
    build_s = 0.
    for lo in range(0, len(indices), args.batch_size):
        hi = min(lo + args.batch_size, len(indices))
        begin = perf_counter()
        prepare = args.mode == "regional" and len(calls) == args.prepare_after
        if prepare:
            model.prepare_ray_acceleration()
        ready = perf_counter()
        result = trace(slice(lo, hi))
        end = perf_counter()
        if prepare:
            build_s = ready - begin
        cumulative += end - begin
        peak = peak_bytes()
        segments, delta = verify_rows(result, indices[lo:hi], reference)
        del result
        calls.append(dict(call=len(calls) + 1, rays=hi-lo, segments=segments,
                          build_s=build_s if prepare else 0., trace_s=end-ready,
                          source_s=cumulative, peak_bytes=peak, max_endpoint_delta_cm=delta))
    return dict(model=str(args.model.resolve()), mode=args.mode, order=args.order,
                distance=distance, cold=args.cold, rays=len(indices), batch_size=args.batch_size,
                load_s=load_s, setup_s=setup_s, build_s=build_s, prepare_after=args.prepare_after,
                first_result_s=calls[0]["source_s"], source_s=cumulative,
                peak_bytes=calls[-1]["peak_bytes"], calls=calls)


def replay(args):
    origins, directions, distance, indices, reference = read_workload(args)
    count = len(indices)
    if args.batch_size and args.mode == "regional" and args.prepare_after * args.batch_size >= count:
        raise ValueError("regional preparation must occur before the last call returns")
    calibration = read_calibration(args, distance)
    start = perf_counter()
    model = pyalea.load_mcnp(str(args.model))
    loaded = perf_counter()
    if not args.cold:
        model.prepare_query_acceleration()
    ready = perf_counter()

    def trace(selection):
        return model.raycast_hier_fast_batch(origins[selection], directions[selection],
                                             distance, include_surfaces=False)

    if args.batch_size is not None:
        return replay_calls(args, model, trace, indices, reference, distance,
                            loaded-start, ready-loaded)

    prefix_count = 0
    prefix_s = build_s = 0.
    estimated_saving = None
    prepare = args.mode == "regional"
    result = []
    remainder = slice(None)
    scattered = args.mode == "stratified" and count > 16
    if args.mode in ("prefix", "stratified") and count > 16:
        # The cold first ray can build caches. It belongs in caller time but
        # its one-time setup must not be multiplied by the remaining ray count.
        result = trace(slice(0, 1))
        if scattered:
            # One reproducible sample per stratum, even for a sorted batch.
            edges = np.linspace(1, count, 16, dtype=np.int64)
            sample = np.random.default_rng(42).integers(edges[:-1], edges[1:])
            mask = np.ones(count, dtype=bool)
            mask[0] = False
            mask[sample] = False
            remainder = np.flatnonzero(mask)
        else:
            sample = slice(1, 16)
            remainder = slice(16, None)
        sample_start = perf_counter()
        result.extend(trace(sample))
        sample_s = perf_counter() - sample_start
        prefix_count = 16
        prefix_s = perf_counter() - ready
        regional_per_ray = calibration["tail_s"] / calibration["rays"]
        estimated_saving = (count - prefix_count) * (sample_s / 15 - regional_per_ray)
        prepare = estimated_saving > calibration["build_s"]
    if prepare:
        begin = perf_counter()
        model.prepare_ray_acceleration()
        build_s = perf_counter() - begin
    begin = perf_counter()
    if args.mode == "lazy-chunks":
        # Single-worker hierarchy control. Seven rows never trigger automatic
        # recovery preparation; include all calls and list assembly in the timer.
        tail = []
        for lo in range(0, count, 7):
            tail.extend(trace(slice(lo, lo + 7)))
    else:
        tail = trace(remainder)
    if scattered:
        ordered = [None] * count
        ordered[0] = result[0]
        for index, row in zip(sample, result[1:], strict=True):
            ordered[index] = row
        for index, row in zip(remainder, tail, strict=True):
            ordered[index] = row
        result = ordered
    elif prefix_count:
        result.extend(tail)
    else:
        result = tail
    end = perf_counter()
    peak = peak_bytes()
    segments, max_delta = verify_rows(result, indices, reference)
    return dict(model=str(args.model.resolve()), distance=distance, mode=args.mode,
                cold=args.cold, order=args.order, rays=count, segments=segments,
                load_s=loaded-start, setup_s=ready-loaded, call_s=end-ready, source_s=end-start,
                prefix_count=prefix_count, prefix_s=prefix_s, estimated_saving_s=estimated_saving,
                requested_regional=prepare, build_s=build_s, tail_s=end-begin,
                after_load_s=end-loaded, peak_bytes=peak,
                max_endpoint_delta_cm=max_delta)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", type=Path)
    parser.add_argument("workload", type=Path)
    parser.add_argument("reference", type=Path, help="ordered geometry JSON from model_lifecycle public")
    parser.add_argument("mode", choices=("hierarchy", "regional", "prefix", "stratified", "lazy-chunks"))
    parser.add_argument("--calibration", type=Path, help="previous regional replay on the same model")
    parser.add_argument("--repeat", type=int, default=1)
    parser.add_argument("--batch-size", type=int, help="split repeated passes into consecutive public calls")
    parser.add_argument("--prepare-after", type=int, default=0,
                        help="complete this many calls before explicit regional preparation")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--order", choices=("input", "sparse-first", "dense-first"), default="input")
    parser.add_argument("--cold", action="store_true", help="leave initial hierarchy preparation to the query")
    args = parser.parse_args()
    if args.batch_size is not None and (args.batch_size < 1 or args.mode not in ("hierarchy", "regional")):
        parser.error("--batch-size requires a positive size and hierarchy or regional mode")
    if args.prepare_after < 0 or (args.prepare_after and (not args.batch_size or args.mode != "regional")):
        parser.error("--prepare-after requires consecutive regional calls and a nonnegative count")
    if args.mode == "lazy-chunks" and not args.cold:
        parser.error("use --cold for the lazy hierarchy control")
    if (args.repeat < 1 or (args.limit is not None and args.limit < 1)
            or (args.mode in ("prefix", "stratified")) != bool(args.calibration)):
        parser.error("use positive batch sizes and supply a calibration only for sampled modes")
    print(json.dumps(replay(args), indent=2))
