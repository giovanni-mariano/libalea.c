#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0

"""Run a reproducible bounded slice-error workload on an external MCNP model."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import subprocess
import time
from pathlib import Path
from typing import Any

import pyalea

try:
    import resource
except ImportError:  # Windows does not provide the POSIX resource module.
    resource = None


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _current_rss_bytes() -> int | None:
    try:
        for line in Path("/proc/self/status").read_text().splitlines():
            if line.startswith("VmRSS:"):
                return int(line.split()[1]) * 1024
    except (OSError, ValueError, IndexError):
        pass
    return None


def _peak_rss_bytes() -> int | None:
    if resource is None:
        return None
    value = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(value if platform.system() == "Darwin" else value * 1024)


def _git_revision(root: Path) -> str | None:
    try:
        return subprocess.run(
            ["git", "rev-parse", "HEAD"], cwd=root, check=True,
            text=True, stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return None


def _parse_pages(value: str, page_count: int) -> list[int]:
    if value == "all":
        return list(range(page_count))
    pages: list[int] = []
    for part in value.split(","):
        part = part.strip()
        if not part:
            continue
        if ":" in part:
            start_text, stop_text = part.split(":", 1)
            start, stop = int(start_text), int(stop_text)
            pages.extend(range(start, stop))
        else:
            pages.append(int(part))
    if len(set(pages)) != len(pages):
        raise ValueError("page selection contains duplicates")
    if any(page < 0 or page >= page_count for page in pages):
        raise ValueError(f"page selection must lie in [0, {page_count})")
    return pages


def run_benchmark(model_path: Path, *, value: float = 80.0,
                  bounds=(500.0, 1500.0, -500.0, 500.0),
                  tile_columns: int = 22, tile_rows: int = 20,
                  pages: str = "0", interior_probes_per_axis: int = 3,
                  max_index_bytes: int = 32 * 1024 * 1024) -> dict[str, Any]:
    root = Path(__file__).resolve().parents[1]
    model_path = model_path.resolve()
    started = time.perf_counter()
    rss_before_load = _current_rss_bytes()
    load_started = time.perf_counter()
    system = pyalea.load_mcnp(str(model_path))
    load_seconds = time.perf_counter() - load_started
    rss_after_load = _current_rss_bytes()
    query_started = time.perf_counter()
    query = system.slice_error_query(
        (0.0, 0.0, value), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0),
        bounds, bounds, tile_columns=tile_columns, tile_rows=tile_rows,
        options={"interior_probes_per_axis": interior_probes_per_axis},
        max_index_bytes=max_index_bytes)
    query_seconds = time.perf_counter() - query_started
    rss_after_query = _current_rss_bytes()
    selected_pages = _parse_pages(pages, query.page_count)
    page_records = []
    first_witness_seconds = None
    scan_started = time.perf_counter()
    try:
        for page_index in selected_pages:
            page_started = time.perf_counter()
            page = query.run_page(page_index)
            page_wall_seconds = time.perf_counter() - page_started
            receipt = dict(page["receipt"])
            witnesses = list(page.get("witnesses", ()))
            if witnesses and first_witness_seconds is None:
                first_witness_seconds = time.perf_counter() - scan_started
            page_records.append({
                "page_index": page_index,
                "wall_seconds": page_wall_seconds,
                "receipt": receipt,
                "witness_count": len(witnesses),
                "gap_witness_count": sum(
                    item["kind"] == "gap" for item in witnesses),
                "overlap_witness_count": sum(
                    item["kind"] == "overlap" for item in witnesses),
                "witnesses": witnesses,
                "unresolved": list(page.get("unresolved", ())),
            })
    finally:
        query.close()
    scan_seconds = time.perf_counter() - scan_started
    receipts = [item["receipt"] for item in page_records]
    return {
        "schema": "libalea.geometry-error-model-benchmark.v1",
        "metadata": {
            "model_path": str(model_path),
            "model_size_bytes": model_path.stat().st_size,
            "model_sha256": _sha256(model_path),
            "libalea_version": pyalea.version(),
            "git_revision": _git_revision(root),
            "python_version": platform.python_version(),
            "platform": platform.platform(),
            "processor": platform.processor(),
            "cpu_count": os.cpu_count(),
        },
        "workload": {
            "preset": "elite_z80",
            "axis": "z", "value": value,
            "bounds": list(bounds),
            "tile_columns": tile_columns, "tile_rows": tile_rows,
            "total_page_count": tile_columns * tile_rows,
            "selected_pages": selected_pages,
            "interior_probes_per_axis": interior_probes_per_axis,
            "max_index_bytes": max_index_bytes,
        },
        "timing": {
            "load_seconds": load_seconds,
            "query_preparation_seconds": query_seconds,
            "scan_seconds": scan_seconds,
            "first_confirmed_witness_seconds": first_witness_seconds,
            "total_seconds": time.perf_counter() - started,
        },
        "memory": {
            "rss_before_load_bytes": rss_before_load,
            "rss_after_load_bytes": rss_after_load,
            "rss_after_query_bytes": rss_after_query,
            "peak_process_rss_bytes": _peak_rss_bytes(),
            "query_index_bytes": max(
                (item["query_index_bytes"] for item in receipts), default=0),
            "peak_page_scratch_bytes": max(
                (item["peak_scratch_bytes"] for item in receipts), default=0),
        },
        "summary": {
            "completed_page_count": len(page_records),
            "classified_page_count": sum(
                item["scope_classified"] for item in receipts),
            "unresolved_page_count": sum(
                not item["scope_classified"] for item in receipts),
            "confirmed_gap_witness_count": sum(
                item["confirmed_gap_witness_count"] for item in receipts),
            "confirmed_overlap_witness_count": sum(
                item["confirmed_overlap_witness_count"] for item in receipts),
            "confirmation_attempt_count": sum(
                item["confirmation_attempt_count"] for item in receipts),
            "confirmation_failure_count": sum(
                item["confirmation_failure_count"] for item in receipts),
            "interior_probe_count": sum(
                item["interior_probe_count"] for item in receipts),
            "candidate_curve_count": sum(
                item["candidate_curves"] for item in receipts),
            "candidate_pair_count": sum(
                item["candidate_pairs_tested"] for item in receipts),
            "symbolic_one_sided_interval_count": sum(
                item["symbolic_one_sided_intervals"] for item in receipts),
        },
        "pages": page_records,
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--pages", default="0",
                        help="comma-separated ordinals, half-open ranges, or all")
    parser.add_argument("--value", type=float, default=80.0)
    parser.add_argument("--bounds", type=float, nargs=4,
                        default=(500.0, 1500.0, -500.0, 500.0))
    parser.add_argument("--tile-columns", type=int, default=22)
    parser.add_argument("--tile-rows", type=int, default=20)
    parser.add_argument("--interior-probes-per-axis", type=int, default=3)
    parser.add_argument("--max-index-bytes", type=int,
                        default=32 * 1024 * 1024)
    args = parser.parse_args(argv)
    if (args.tile_columns < 1 or args.tile_rows < 1 or
            not 0 <= args.interior_probes_per_axis <= 64 or
            args.max_index_bytes < 0):
        parser.error("tile counts, interior probes, or index budget are invalid")
    report = run_benchmark(
        args.model, value=args.value, bounds=tuple(args.bounds),
        tile_columns=args.tile_columns, tile_rows=args.tile_rows,
        pages=args.pages,
        interior_probes_per_axis=args.interior_probes_per_axis,
        max_index_bytes=args.max_index_bytes)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
