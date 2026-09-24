#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0

"""Inject and detect a known point defect in an external geometry model.

The input file is never modified.  A deterministic interior probe with one
complete terminal owner is selected from a requested slice page.  Removing
that owner seeds a gap; duplicating its CSG root in the same universe seeds an
overlap.  The tool verifies the point classification before and after the
mutation, then runs the normal slice-error page on both states.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import time
from pathlib import Path
from typing import Any

import pyalea


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _page_bounds(bounds: tuple[float, float, float, float], columns: int,
                 rows: int, page_index: int) -> tuple[float, float, float, float]:
    page_count = columns * rows
    if columns < 1 or rows < 1 or not 0 <= page_index < page_count:
        raise ValueError("invalid page layout or page index")
    u0, u1, v0, v1 = bounds
    col, row = page_index % columns, page_index // columns
    return (
        u0 + (u1 - u0) * col / columns,
        u0 + (u1 - u0) * (col + 1) / columns,
        v0 + (v1 - v0) * row / rows,
        v0 + (v1 - v0) * (row + 1) / rows,
    )


def _probe_points(page_bounds: tuple[float, float, float, float],
                  probes_per_axis: int, value: float):
    u0, u1, v0, v1 = page_bounds
    for j in range(probes_per_axis):
        v = v0 + (j + 0.5) * (v1 - v0) / probes_per_axis
        for i in range(probes_per_axis):
            u = u0 + (i + 0.5) * (u1 - u0) / probes_per_axis
            yield (u, v, value)


def _select_seed(system: Any,
                 page_bounds: tuple[float, float, float, float],
                 probes_per_axis: int, value: float) -> dict[str, Any]:
    rejected: list[dict[str, Any]] = []
    for point in _probe_points(page_bounds, probes_per_axis, value):
        coverage = system.find_all_cells_coverage(*point, max_hits=4096)
        owners = coverage["owners"]
        reason = None
        if coverage["truncated"]:
            reason = "truncated"
        elif coverage["kind"] != "unique" or len(owners) != 1:
            reason = str(coverage["kind"])
        else:
            owner = dict(owners[0])
            cell = dict(system.get_cell_by_index(owner["cell_index"]))
            if cell.get("fill_universe", 0) > 0 or cell.get("lat_type", 0):
                reason = "container_owner"
            elif not system.point_inside(
                    cell["root_node"], owner["local_x"],
                    owner["local_y"], owner["local_z"]):
                reason = "root_predicate_disagrees"
            else:
                return {
                    "world_point": list(point),
                    "coverage": coverage,
                    "owner": owner,
                    "cell": cell,
                    "rejected_probe_count": len(rejected),
                    "rejected_probes": rejected,
                }
        rejected.append({"world_point": list(point), "reason": reason})
    raise RuntimeError(
        "no complete unique terminal owner exists at the configured probes")


def _new_cell_id(system: Any) -> int:
    candidate = 2_000_000_000
    while system.cell_find(candidate) is not None:
        candidate -= 1
    return candidate


def _new_surface_id(system: Any) -> int:
    used = {system.surface_id_at(index)
            for index in range(system.surface_count)}
    for candidate in range(1, len(used) + 2):
        if candidate not in used:
            return candidate
    raise RuntimeError("failed to allocate a compact synthetic surface ID")


def _apply_mutation(system: Any, seed: dict[str, Any], kind: str,
                    radius: float) -> dict[str, Any]:
    owner = seed["owner"]
    cell = seed["cell"]
    surface_id = _new_surface_id(system)
    _, _, inside = system.sphere_surface(
        surface_id, owner["local_x"], owner["local_y"], owner["local_z"],
        radius)
    if kind == "gap":
        carved_region = system.create_difference(cell["root_node"], inside)
        system.cell_set_region(owner["cell_index"], carved_region)
        return {
            "kind": kind,
            "operation": "carve_sphere_from_terminal_cell",
            "source_cell_id": owner["cell_id"],
            "source_cell_index": owner["cell_index"],
            "inserted_surface_id": surface_id,
            "local_center": [owner["local_x"], owner["local_y"],
                             owner["local_z"]],
            "radius": radius,
            "affected_universe_id": owner["universe_id"],
            "shared_universe_may_affect_other_occurrences":
                owner["universe_id"] != 0,
        }
    if kind == "overlap":
        new_cell_id = _new_cell_id(system)
        new_index = system.add_cell(
            new_cell_id, inside, -1, 0.0,
            owner["universe_id"])
        return {
            "kind": kind,
            "operation": "add_competing_sphere_cell",
            "source_cell_id": owner["cell_id"],
            "inserted_cell_id": new_cell_id,
            "inserted_cell_index": new_index,
            "inserted_surface_id": surface_id,
            "local_center": [owner["local_x"], owner["local_y"],
                             owner["local_z"]],
            "radius": radius,
            "affected_universe_id": owner["universe_id"],
            "shared_universe_may_affect_other_occurrences":
                owner["universe_id"] != 0,
        }
    raise ValueError("mutation kind must be gap or overlap")


def _scan_page(system: Any, *, value: float,
               bounds: tuple[float, float, float, float], columns: int,
               rows: int, page_index: int, probes_per_axis: int,
               max_index_bytes: int) -> dict[str, Any]:
    started = time.perf_counter()
    with system.slice_error_query(
            (0.0, 0.0, value), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0),
            bounds, bounds, tile_columns=columns, tile_rows=rows,
            options={"interior_probes_per_axis": probes_per_axis},
            max_index_bytes=max_index_bytes) as query:
        prepared_seconds = time.perf_counter() - started
        scan_started = time.perf_counter()
        page = query.run_page(page_index)
        scan_seconds = time.perf_counter() - scan_started
    return {
        "query_preparation_seconds": prepared_seconds,
        "scan_seconds": scan_seconds,
        "receipt": dict(page["receipt"]),
        "witnesses": list(page.get("witnesses", ())),
        "unresolved": list(page.get("unresolved", ())),
    }


def _matching_seed_witnesses(scan: dict[str, Any], seed_point: list[float],
                             kind: str) -> list[dict[str, Any]]:
    u, v, _ = seed_point
    scale = max(1.0, abs(u), abs(v))
    tolerance = 128.0 * 2.220446049250313e-16 * scale
    return [witness for witness in scan["witnesses"]
            if witness["kind"] == kind and
            abs(witness["uv"][0] - u) <= tolerance and
            abs(witness["uv"][1] - v) <= tolerance]


def run_seeded_system(system: Any, *, mutation: str, value: float = 80.0,
                      bounds=(500.0, 1500.0, -500.0, 500.0),
                      tile_columns: int = 22, tile_rows: int = 20,
                      page_index: int = 0, probes_per_axis: int = 3,
                      max_index_bytes: int = 32 * 1024 * 1024,
                      defect_radius: float = 0.0) -> dict[str, Any]:
    bounds = tuple(float(item) for item in bounds)
    core_bounds = _page_bounds(bounds, tile_columns, tile_rows, page_index)
    if defect_radius <= 0.0:
        defect_radius = min(
            core_bounds[1] - core_bounds[0],
            core_bounds[3] - core_bounds[2]) / (8.0 * probes_per_axis)
    if defect_radius <= 0.0:
        raise ValueError("defect radius must be positive")
    seed_started = time.perf_counter()
    seed = _select_seed(system, core_bounds, probes_per_axis, value)
    seed_seconds = time.perf_counter() - seed_started
    baseline_scan = _scan_page(
        system, value=value, bounds=bounds, columns=tile_columns,
        rows=tile_rows, page_index=page_index,
        probes_per_axis=probes_per_axis, max_index_bytes=max_index_bytes)
    mutation_started = time.perf_counter()
    mutation_record = _apply_mutation(system, seed, mutation, defect_radius)
    system.prepare_query_acceleration()
    mutation_seconds = time.perf_counter() - mutation_started
    point = seed["world_point"]
    verification_started = time.perf_counter()
    mutated_coverage = system.find_all_cells_coverage(*point, max_hits=4096)
    verification_seconds = time.perf_counter() - verification_started
    mutated_scan = _scan_page(
        system, value=value, bounds=bounds, columns=tile_columns,
        rows=tile_rows, page_index=page_index,
        probes_per_axis=probes_per_axis, max_index_bytes=max_index_bytes)
    matching = _matching_seed_witnesses(mutated_scan, point, mutation)
    expected_point_kind = mutation
    return {
        "schema": "libalea.geometry-error-seeded-benchmark.v1",
        "workload": {
            "axis": "z", "value": value, "bounds": list(bounds),
            "tile_columns": tile_columns, "tile_rows": tile_rows,
            "page_index": page_index, "core_bounds": list(core_bounds),
            "interior_probes_per_axis": probes_per_axis,
            "defect_radius": defect_radius,
            "max_index_bytes": max_index_bytes,
        },
        "seed": seed,
        "mutation": mutation_record,
        "oracle": {
            "scope": "verified_point",
            "world_point": point,
            "baseline_expected_kind": "unique",
            "mutated_expected_kind": expected_point_kind,
            "baseline_observed_kind": seed["coverage"]["kind"],
            "mutated_observed_kind": mutated_coverage["kind"],
            "root_predicate_checked_in_owner_local_coordinates": True,
        },
        "baseline_scan": baseline_scan,
        "mutated_scan": mutated_scan,
        "phase_timing": {
            "seed_selection_seconds": seed_seconds,
            "mutation_seconds": mutation_seconds,
            "mutated_point_verification_seconds": verification_seconds,
        },
        "result": {
            "point_mutation_verified":
                mutated_coverage["kind"] == expected_point_kind and
                not mutated_coverage["truncated"],
            "seed_witness_detected": bool(matching),
            "matching_seed_witnesses": matching,
        },
    }


def run_benchmark(model_path: Path, **kwargs: Any) -> dict[str, Any]:
    model_path = model_path.resolve()
    started = time.perf_counter()
    load_started = time.perf_counter()
    system = pyalea.load_mcnp(str(model_path))
    load_seconds = time.perf_counter() - load_started
    report = run_seeded_system(system, **kwargs)
    report["metadata"] = {
        "model_path": str(model_path),
        "model_size_bytes": model_path.stat().st_size,
        "model_sha256": _sha256(model_path),
        "libalea_version": pyalea.version(),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
    }
    report["timing"] = {
        "load_seconds": load_seconds,
        "total_seconds": time.perf_counter() - started,
    }
    return report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("model", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mutation", choices=("gap", "overlap"), required=True)
    parser.add_argument("--page", type=int, default=0)
    parser.add_argument("--value", type=float, default=80.0)
    parser.add_argument("--bounds", type=float, nargs=4,
                        default=(500.0, 1500.0, -500.0, 500.0))
    parser.add_argument("--tile-columns", type=int, default=22)
    parser.add_argument("--tile-rows", type=int, default=20)
    parser.add_argument("--interior-probes-per-axis", type=int, default=3)
    parser.add_argument("--max-index-bytes", type=int,
                        default=32 * 1024 * 1024)
    parser.add_argument("--defect-radius", type=float, default=0.0,
                        help="local sphere radius; zero selects a page-relative default")
    args = parser.parse_args(argv)
    if args.interior_probes_per_axis < 1 or args.interior_probes_per_axis > 64:
        parser.error("interior probe count must be between 1 and 64")
    report = run_benchmark(
        args.model, mutation=args.mutation, value=args.value,
        bounds=tuple(args.bounds), tile_columns=args.tile_columns,
        tile_rows=args.tile_rows, page_index=args.page,
        probes_per_axis=args.interior_probes_per_axis,
        max_index_bytes=args.max_index_bytes,
        defect_radius=args.defect_radius)
    args.output.write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0 if (report["result"]["point_mutation_verified"] and
                 report["result"]["seed_witness_detected"]) else 1


if __name__ == "__main__":
    raise SystemExit(main())
