#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0

"""Known-answer acceptance harness for slice geometry defect witnesses.

The expected regions below use direct analytic arithmetic.  They deliberately
do not call libalea containment or coverage predicates.  The output is a
machine-readable baseline, with one result per expected geometric defect rather
than one result per probe emitted by the validator.
"""

from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable

import pyalea


PointPredicate = Callable[[float, float], bool]


@dataclass(frozen=True)
class ExpectedDefect:
    defect_id: str
    kind: str
    oracle: dict[str, Any]
    contains: PointPredicate
    owner_cell_ids: tuple[int, ...]
    target_depth: int
    require_distinct_occurrences: bool = False


@dataclass(frozen=True)
class Fixture:
    name: str
    family: str
    system: Any
    view_bounds: tuple[float, float, float, float]
    required_bounds: tuple[float, float, float, float]
    tile_columns: int
    tile_rows: int
    expected: tuple[ExpectedDefect, ...]


def _rectangle(u_min: float, u_max: float, v_min: float, v_max: float):
    return lambda u, v: u_min <= u <= u_max and v_min <= v <= v_max


def _outside_circle(cx: float, cy: float, radius: float, margin: float = 0.0):
    radius2 = (radius + margin) ** 2
    return lambda u, v: (u - cx) ** 2 + (v - cy) ** 2 > radius2


def _adjacent_fixture(name: str, mode: str) -> Fixture:
    system = pyalea.System()
    if mode == "clean":
        _, positive, negative = system.plane_surface(100, 1, 0, 0, 0)
        system.add_cell(101, negative)
        system.add_cell(102, positive)
        expected = ()
    elif mode == "gap":
        _, _, left = system.plane_surface(100, 1, 0, 0, 0.2)
        _, right, _ = system.plane_surface(101, 1, 0, 0, -0.2)
        system.add_cell(101, left)
        system.add_cell(102, right)
        expected = (ExpectedDefect(
            "central-gap", "gap",
            {"type": "rectangle", "uv_min": [-0.2, -0.5],
             "uv_max": [0.2, 0.5]},
            _rectangle(-0.2, 0.2, -0.5, 0.5), (), 0),)
    elif mode == "overlap":
        _, _, left = system.plane_surface(100, 1, 0, 0, -0.2)
        _, right, _ = system.plane_surface(101, 1, 0, 0, 0.2)
        system.add_cell(101, left)
        system.add_cell(102, right)
        expected = (ExpectedDefect(
            "central-overlap", "overlap",
            {"type": "rectangle", "uv_min": [-0.2, -0.5],
             "uv_max": [0.2, 0.5]},
            _rectangle(-0.2, 0.2, -0.5, 0.5), (101, 102), 0),)
    else:  # pragma: no cover - fixture construction is internal
        raise ValueError(mode)
    return Fixture(name, "adjacent_halfspaces", system,
                   (-1.0, 1.0, -0.5, 0.5),
                   (-1.0, 1.0, -0.5, 0.5), 2, 1, expected)


def _narrow_offset_gap_fixture() -> Fixture:
    system = pyalea.System()
    _, _, left = system.plane_surface(200, 1, 0, 0, -0.123)
    _, right, _ = system.plane_surface(201, 1, 0, 0, -0.124)
    system.add_cell(201, left)
    system.add_cell(202, right)
    expected = ExpectedDefect(
        "off-grid-narrow-gap", "gap",
        {"type": "rectangle", "uv_min": [0.123, -0.4],
         "uv_max": [0.124, 0.6]},
        _rectangle(0.123, 0.124, -0.4, 0.6), (), 0)
    return Fixture("narrow_offset_gap", "narrow_defect", system,
                   (-1.0, 1.0, -0.4, 0.6),
                   (-1.0, 1.0, -0.4, 0.6), 3, 2, (expected,))


def _viewport_inside_gap_fixture() -> Fixture:
    fixture = _adjacent_fixture("viewport_inside_gap", "gap")
    return Fixture(fixture.name, "enclosed_or_hidden_defect", fixture.system,
                   (-0.1, 0.1, -0.25, 0.25),
                   (-0.1, 0.1, -0.25, 0.25), 1, 1,
                   (ExpectedDefect(
                       "viewport-gap", "gap",
                       {"type": "rectangle", "uv_min": [-0.1, -0.25],
                        "uv_max": [0.1, 0.25]},
                       _rectangle(-0.1, 0.1, -0.25, 0.25), (), 0),))


def _transformed_fill_gap_fixture() -> Fixture:
    system = pyalea.System()
    _, _, parent_region = system.plane_surface(300, 1, 0, 0, -10)
    _, _, child_region = system.plane_surface(301, 1, 0, 0, -2)
    parent = system.add_cell(300, parent_region)
    system.add_cell(301, child_region, universe_id=8)
    system.add_transform(300, (3.0, 0.0, 0.0))
    system.set_fill(parent, 8, 300)
    expected = ExpectedDefect(
        "translated-child-gap", "gap",
        {"type": "rectangle", "uv_min": [5.5, -0.5],
         "uv_max": [6.0, 0.5]},
        _rectangle(5.5, 6.0, -0.5, 0.5), (), 1)
    return Fixture("transformed_fill_gap", "transformed_hierarchy", system,
                   (5.5, 6.0, -0.5, 0.5), (5.5, 6.0, -0.5, 0.5),
                   1, 1, (expected,))


def _repeated_fill_overlap_fixture() -> Fixture:
    system = pyalea.System()
    _, _, parent_region = system.plane_surface(400, 1, 0, 0, -10)
    _, _, child_region = system.plane_surface(401, 1, 0, 0, -2)
    first = system.add_cell(400, parent_region)
    second = system.add_cell(402, parent_region)
    system.add_cell(401, child_region, universe_id=8)
    system.add_transform(400, (3.0, 0.0, 0.0))
    system.add_transform(402, (4.0, 0.0, 0.0))
    system.set_fill(first, 8, 400)
    system.set_fill(second, 8, 402)
    expected = ExpectedDefect(
        "repeated-child-overlap", "overlap",
        {"type": "rectangle", "uv_min": [4.0, -0.5],
         "uv_max": [4.5, 0.5]},
        _rectangle(4.0, 4.5, -0.5, 0.5), (401, 401), 1, True)
    return Fixture("repeated_fill_overlap", "repeated_occurrence", system,
                   (4.0, 4.5, -0.5, 0.5), (4.0, 4.5, -0.5, 0.5),
                   1, 1, (expected,))


def _curved_gap_fixture(name: str, cylinder: bool) -> Fixture:
    system = pyalea.System()
    if cylinder:
        _, _, inside = system.cylinder_z_surface(500, 0, 0, 0.5)
        family = "curved_boundary_incomplete_page_proof"
    else:
        _, _, inside = system.sphere_surface(500, 0, 0, 0, 0.5)
        family = "curved_boundary"
    system.add_cell(501, inside)
    expected = ExpectedDefect(
        "curved-exterior-gap", "gap",
        {"type": "outside_circle", "center_uv": [0.0, 0.0],
         "radius": 0.5, "strict_margin": 1e-9},
        _outside_circle(0.0, 0.0, 0.5, 1e-9), (), 0)
    return Fixture(name, family, system, (-1.0, 1.0, -1.0, 1.0),
                   (-1.0, 1.0, -1.0, 1.0), 1, 1, (expected,))


def _off_center_enclosed_overlap_fixture() -> Fixture:
    system = pyalea.System()
    _, _, broad = system.sphere_surface(600, 0, 0, 0, 10.0)
    _, _, enclosed = system.sphere_surface(601, 2.0 / 3.0, 0, 0, 0.12)
    system.add_cell(601, broad)
    system.add_cell(602, enclosed)
    expected = ExpectedDefect(
        "off-center-enclosed-overlap", "overlap",
        {"type": "inside_circle", "center_uv": [2.0 / 3.0, 0.0],
         "radius": 0.12},
        lambda u, v: (u - 2.0 / 3.0) ** 2 + v ** 2 < 0.12 ** 2,
        (601, 602), 0)
    return Fixture("off_center_enclosed_overlap",
                   "enclosed_or_hidden_defect", system,
                   (-1.0, 1.0, -1.0, 1.0),
                   (-1.0, 1.0, -1.0, 1.0), 1, 1, (expected,))


def _rotated_fill_fixture(name: str, gap: bool) -> Fixture:
    system = pyalea.System()
    _, _, parent_region = system.sphere_surface(700, 0, 0, 0, 10.0)
    _, _, child_region = system.plane_surface(701, 1, 0, 0, -2.0)
    parent = system.add_cell(700, parent_region)
    system.add_cell(701, child_region, universe_id=8)
    system.add_transform(700, (
        0, 0, 0,
        0, -1, 0,
        1, 0, 0,
        0, 0, 1))
    system.set_fill(parent, 8, 700)
    if gap:
        bounds = (-0.5, 0.5, -2.5, -1.5)
        expected = (ExpectedDefect(
            "rotated-child-gap", "gap",
            {"type": "rectangle", "uv_min": [-0.5, -2.5],
             "uv_max": [0.5, -2.0]},
            _rectangle(-0.5, 0.5, -2.5, -2.0), (), 1),)
    else:
        bounds = (-0.5, 0.5, -0.5, 0.5)
        expected = ()
    return Fixture(name, "rotated_hierarchy", system, bounds, bounds,
                   1, 1, expected)


def _macrobody_gap_fixture(name: str, macrobody: str) -> Fixture:
    system = pyalea.System()
    if macrobody == "box":
        _, _, inside = system.box_surface(
            800, -0.5, 0.5, -0.5, 0.5, -1.0, 1.0)
        oracle = {"type": "outside_rectangle", "uv_min": [-0.5, -0.5],
                  "uv_max": [0.5, 0.5]}
        contains = lambda u, v: abs(u) > 0.5 or abs(v) > 0.5
    elif macrobody == "rcc":
        _, _, inside = system.rcc_surface(
            800, 0, 0, -1.0, 0, 0, 2.0, 0.5)
        oracle = {"type": "outside_circle", "center_uv": [0.0, 0.0],
                  "radius": 0.5, "strict_margin": 1e-9}
        contains = _outside_circle(0.0, 0.0, 0.5, 1e-9)
    else:  # pragma: no cover - fixture construction is internal
        raise ValueError(macrobody)
    system.add_cell(801, inside)
    expected = ExpectedDefect(
        f"{macrobody}-exterior-gap", "gap", oracle, contains, (), 0)
    return Fixture(name, "native_macrobody", system,
                   (-1.0, 1.0, -1.0, 1.0),
                   (-1.0, 1.0, -1.0, 1.0), 1, 1, (expected,))


def _empty_required_domain_fixture() -> Fixture:
    system = pyalea.System()
    expected = ExpectedDefect(
        "required-domain-gap", "gap",
        {"type": "rectangle", "uv_min": [-1.0, -1.0],
         "uv_max": [1.0, 1.0]},
        _rectangle(-1.0, 1.0, -1.0, 1.0), (), 0)
    return Fixture("empty_required_domain", "required_coverage", system,
                   (-1.0, 1.0, -1.0, 1.0),
                   (-1.0, 1.0, -1.0, 1.0), 1, 1, (expected,))


def _valid_lattice_seam_fixture() -> Fixture:
    model_path = (Path(__file__).resolve().parents[1] / "tests" / "data" /
                  "mcnp_lattice_eval.mcnp")
    system = pyalea.load_mcnp(str(model_path))
    bounds = (0.9, 1.1, -0.5, 0.5)
    return Fixture("valid_rectangular_lattice_seam", "repeated_lattice",
                   system, bounds, bounds, 1, 1, ())


def _defect_beside_unsupported_fixture() -> Fixture:
    # The plane cells form an ordinary overlap at x=[-1.4, -0.6].  A remote
    # BOX macrobody is referenced by an analytically empty cell (inside AND
    # outside the same primitive), so it cannot change point ownership.  Its
    # contour arrangement remains unsupported and must not discard the
    # independently verified overlap evidence.
    system = pyalea.System()
    _, _, left = system.plane_surface(900, 1, 0, 0, 0.6)
    _, right, _ = system.plane_surface(901, 1, 0, 0, 1.4)
    system.add_cell(901, left)
    system.add_cell(902, right)
    _, outside_box, inside_box = system.box_surface(
        902, 2.0, 3.0, -0.25, 0.25, -1.0, 1.0)
    empty_macrobody = system.create_intersection(outside_box, inside_box)
    system.add_cell(903, empty_macrobody)
    expected = ExpectedDefect(
        "overlap-beside-unsupported-macrobody", "overlap",
        {"type": "rectangle", "uv_min": [-1.4, -0.5],
         "uv_max": [-0.6, 0.5]},
        _rectangle(-1.4, -0.6, -0.5, 0.5), (901, 902), 0)
    return Fixture("defect_beside_unsupported_macrobody",
                   "mixed_supported_and_unsupported", system,
                   (-2.0, 4.0, -0.5, 0.5),
                   (-2.0, 4.0, -0.5, 0.5), 6, 1, (expected,))


def fixtures() -> list[Fixture]:
    return [
        _adjacent_fixture("adjacent_clean", "clean"),
        _adjacent_fixture("adjacent_gap", "gap"),
        _adjacent_fixture("adjacent_overlap", "overlap"),
        _narrow_offset_gap_fixture(),
        _viewport_inside_gap_fixture(),
        _transformed_fill_gap_fixture(),
        _repeated_fill_overlap_fixture(),
        _rotated_fill_fixture("rotated_fill_clean", False),
        _rotated_fill_fixture("rotated_fill_gap", True),
        _valid_lattice_seam_fixture(),
        _defect_beside_unsupported_fixture(),
        _off_center_enclosed_overlap_fixture(),
        _empty_required_domain_fixture(),
        _curved_gap_fixture("sphere_exterior_gap", False),
        _curved_gap_fixture("cylinder_exterior_gap", True),
        _macrobody_gap_fixture("box_macrobody_exterior_gap", "box"),
        _macrobody_gap_fixture("rcc_macrobody_exterior_gap", "rcc"),
    ]


def _valid_evidence(witness: dict[str, Any], expected: ExpectedDefect) -> bool:
    if witness.get("evidence_scope") != "verified_point":
        return False
    if not witness.get("owners_complete", False):
        return False
    owners = witness.get("owners", ())
    cell_ids = tuple(sorted(int(owner["cell_id"]) for owner in owners))
    if cell_ids != tuple(sorted(expected.owner_cell_ids)):
        return False
    if int(witness.get("target_depth", -1)) != expected.target_depth:
        return False
    if expected.require_distinct_occurrences:
        keys = [int(owner["occurrence_key"]) for owner in owners]
        if len(keys) != len(set(keys)):
            return False
    return True


def run_fixture(fixture: Fixture) -> dict[str, Any]:
    started = time.perf_counter()
    with fixture.system.slice_error_query(
            (0.0, 0.0, 0.0), (0.0, 0.0, 1.0), (0.0, 1.0, 0.0),
            fixture.view_bounds, fixture.required_bounds,
            tile_columns=fixture.tile_columns,
            tile_rows=fixture.tile_rows) as query:
        pages = [query.run_page(index) for index in range(query.page_count)]
    wall_seconds = time.perf_counter() - started
    witnesses = [witness for page in pages
                 for witness in page.get("witnesses", ())]
    matched: dict[str, list[int]] = {
        expected.defect_id: [] for expected in fixture.expected}
    invalid_evidence: list[int] = []
    false_positive_witnesses: list[int] = []
    for index, witness in enumerate(witnesses):
        u, v = map(float, witness["uv"])
        candidates = [expected for expected in fixture.expected
                      if witness.get("kind") == expected.kind and
                      expected.contains(u, v)]
        valid = [expected for expected in candidates
                 if _valid_evidence(witness, expected)]
        if valid:
            matched[valid[0].defect_id].append(index)
        elif candidates:
            invalid_evidence.append(index)
        else:
            false_positive_witnesses.append(index)
    missing = [expected.defect_id for expected in fixture.expected
               if not matched[expected.defect_id]]
    duplicate_count = sum(max(0, len(indices) - 1)
                          for indices in matched.values())
    clean_confirmed_records = 0
    if not fixture.expected:
        clean_confirmed_records = len(witnesses)
        for page in pages:
            clean_confirmed_records += sum(
                region.get("kind") in ("gap", "overlap")
                for region in page.get("regions", ()))
            clean_confirmed_records += sum(
                any(interval.get(side) in ("gap", "overlap") for side in
                    ("negative_side_kind", "positive_side_kind"))
                for interval in page.get("intervals", ()))
            clean_confirmed_records += sum(
                any(circle.get(side) in ("gap", "overlap") for side in
                    ("inside_kind", "outside_kind"))
                for circle in page.get("circles", ()))
    receipts = [page["receipt"] for page in pages]
    passed = (not missing and not invalid_evidence and
              not false_positive_witnesses and clean_confirmed_records == 0)
    return {
        "name": fixture.name,
        "family": fixture.family,
        "passed": passed,
        "expected_defects": [{
            "id": expected.defect_id,
            "kind": expected.kind,
            "oracle": expected.oracle,
            "owner_cell_ids": list(expected.owner_cell_ids),
            "target_depth": expected.target_depth,
            "require_distinct_occurrences":
                expected.require_distinct_occurrences,
        } for expected in fixture.expected],
        "matched_defect_ids": [key for key, value in matched.items() if value],
        "missing_defect_ids": missing,
        "false_positive_witness_indices": false_positive_witnesses,
        "invalid_evidence_witness_indices": invalid_evidence,
        "duplicate_witness_count": duplicate_count,
        "clean_confirmed_record_count": clean_confirmed_records,
        "witness_count": len(witnesses),
        "witnesses": witnesses,
        "page_count": len(pages),
        "unresolved_page_count": sum(
            not receipt["scope_classified"] for receipt in receipts),
        "confirmation_attempt_count": sum(
            receipt["confirmation_attempt_count"] for receipt in receipts),
        "confirmation_failure_count": sum(
            receipt["confirmation_failure_count"] for receipt in receipts),
        "candidate_curve_count": sum(
            receipt["candidate_curves"] for receipt in receipts),
        "candidate_pair_count": sum(
            receipt["candidate_pairs_tested"] for receipt in receipts),
        "reported_page_seconds": sum(
            receipt["elapsed_seconds"] for receipt in receipts),
        "wall_seconds": wall_seconds,
    }


def run_acceptance() -> dict[str, Any]:
    started = time.perf_counter()
    cases = [run_fixture(fixture) for fixture in fixtures()]
    expected_count = sum(len(case["expected_defects"]) for case in cases)
    missed_count = sum(len(case["missing_defect_ids"]) for case in cases)
    false_positive_count = sum(
        (len(case["false_positive_witness_indices"])
         if case["expected_defects"] else
         case["clean_confirmed_record_count"])
        for case in cases)
    invalid_evidence_count = sum(
        len(case["invalid_evidence_witness_indices"]) for case in cases)
    return {
        "schema": "libalea.geometry-error-acceptance.v1",
        "metadata": {
            "libalea_version": pyalea.version(),
            "python_version": platform.python_version(),
            "platform": platform.platform(),
        },
        "summary": {
            "passed": all(case["passed"] for case in cases),
            "fixture_count": len(cases),
            "expected_defect_count": expected_count,
            "detected_defect_count": expected_count - missed_count,
            "missed_defect_count": missed_count,
            "false_positive_count": false_positive_count,
            "invalid_evidence_count": invalid_evidence_count,
            "duplicate_witness_count": sum(
                case["duplicate_witness_count"] for case in cases),
            "unresolved_page_count": sum(
                case["unresolved_page_count"] for case in cases),
            "wall_seconds": time.perf_counter() - started,
        },
        "cases": cases,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path,
                        help="write JSON to this path instead of stdout")
    parser.add_argument("--pretty", action="store_true",
                        help="indent the JSON output")
    parser.add_argument("--allow-fail", action="store_true",
                        help="return success even when acceptance cases fail")
    args = parser.parse_args(argv)
    report = run_acceptance()
    payload = json.dumps(report, indent=2 if args.pretty else None,
                         sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(payload, encoding="utf-8")
    else:
        sys.stdout.write(payload)
    return 0 if report["summary"]["passed"] or args.allow_fail else 1


if __name__ == "__main__":
    raise SystemExit(main())
