# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

import gc
import json
import weakref
from pathlib import Path

import numpy as np
import pytest

import pyalea


@pytest.fixture
def populated_system():
    system = pyalea.System()
    _, _, inside = system.sphere_surface(1, 0.0, 0.0, 0.0, 2.0)
    material = system.add_material(7)
    system.add_cell(101, inside, material, 1.0)
    system.build_universe_index()
    return system


def test_native_module_is_private_package_member():
    assert pyalea.System.__module__ == "pyalea._alea"
    assert pyalea.version().count(".") == 2
    assert pyalea.__version__ == pyalea.version()


def test_public_names_are_forwarded():
    assert "System" in pyalea.__all__
    assert pyalea.System is pyalea._alea.System
    assert "_alea" not in pyalea.__all__
    assert set(pyalea.__all__) == {
        name for name in dir(pyalea._alea) if not name.startswith("_")
    }


def test_rhp_constructor_accepts_documented_arity():
    system = pyalea.System()
    index, positive, negative = system.rhp_surface(
        17, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0,
        2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    assert index == 0
    assert positive != negative


def test_surface_boundary_metadata_round_trip():
    system = pyalea.System()
    system.sphere_surface(1, 0.0, 0.0, 0.0, 2.0)
    system.surface_set_boundary(1, "reflective")
    assert system.surface_get_boundary(1) == "reflective"


def test_slice_error_page_reports_verified_boundary_and_unresolved_geometry():
    system = pyalea.System()
    _, _, negative = system.plane_surface(10, 1.0, 0.0, 0.0, 0.0)
    system.add_cell(1, negative)
    args = ((0, 0, 0), (0, 0, 1), (0, 1, 0),
            (-0.25, 0.25, -0.25, 0.25), (-1, 1, -1, 1))

    page = system.slice_error_page(*args)
    assert page["receipt"]["scope_classified"] is True
    assert page["receipt"]["output_complete"] is True
    assert page["receipt"]["query_index_bytes"] > 0
    assert page["receipt"]["close_crossing_observations"] == 0
    assert page["receipt"]["symbolic_one_sided_intervals"] == 0
    assert page["receipt"]["boundary_analysis_policy_version"] == 4
    assert page["receipt"]["close_crossing_relative_tolerance"] == 1e-6
    assert page["receipt"]["numerical_unsafe_probe_intervals"] == 0
    assert page["receipt"]["numerical_unrepresentable_probe_intervals"] == 0
    assert page["receipt"]["numerical_inconsistent_probe_intervals"] == 0
    assert page["receipt"]["core_uv_min"] == (-1, -1)
    assert page["receipt"]["core_uv_max"] == (1, 1)
    assert len(page["context_findings"]) == page["receipt"]["contextual_finding_count"]
    assert all(finding["evidence_scope"] == "context"
               for finding in page["context_findings"])
    assert page["unresolved"] == []
    assert len(page["intervals"]) == 1
    assert page["intervals"][0]["surface_id"] == 10
    assert page["intervals"][0]["negative_owner_cell_ids"] == [1]
    assert page["intervals"][0]["positive_side_kind"] == "gap"
    assert page["regions"][0]["kind"] == "gap"
    assert page["receipt"]["confirmation_attempt_count"] > 0
    assert page["receipt"]["interior_probe_count"] == 9
    assert page["receipt"]["confirmation_failure_count"] == 0
    assert page["receipt"]["confirmed_gap_witness_count"] > 0
    assert page["receipt"]["confirmed_overlap_witness_count"] == 0
    assert page["receipt"]["confirmation_seconds"] >= 0.0
    assert page["receipt"]["elapsed_seconds"] >= (
        page["receipt"]["confirmation_seconds"])
    gap_witnesses = [item for item in page["witnesses"]
                     if item["kind"] == "gap"]
    assert gap_witnesses
    witness = next(item for item in gap_witnesses
                   if item["source"] == "boundary_probe")
    assert witness["evidence_scope"] == "verified_point"
    assert witness["source"] == "boundary_probe"
    assert witness["target_depth"] == 0
    assert witness["owner_count"] == 0
    assert witness["owner_count_lower_bound"] == 0
    assert witness["owners_complete"] is True
    assert witness["owners"] == []

    without_index = system.slice_error_page(*args, max_index_bytes=1)
    assert without_index["receipt"]["scope_classified"] is True
    assert without_index["receipt"]["query_index_bytes"] == 0
    assert without_index["intervals"] == page["intervals"]

    one_probe = system.slice_error_page(
        *args, options={"interior_probes_per_axis": 1})
    assert one_probe["receipt"]["interior_probe_count"] == 1
    with pytest.raises(ValueError):
        system.slice_error_page(
            *args, options={"interior_probes_per_axis": 65})

    sphere = pyalea.System()
    _, _, inside = sphere.sphere_surface(2, 0, 0, 0, 1)
    sphere.add_cell(2, inside)
    unresolved = sphere.slice_error_page(*args)
    assert unresolved["receipt"]["scope_classified"] is False
    assert unresolved["receipt"]["unresolved_reason"] == "unsupported_primitive"
    assert unresolved["unresolved"][0]["reason"] == "unsupported_primitive"
    assert unresolved["intervals"] == []
    assert unresolved["unresolved"]


def test_slice_error_page_reports_symbolic_probe_resolution_separately():
    system = pyalea.System()
    _, xlo_positive, _ = system.plane_surface(1310, 1, 0, 0, -999)
    _, _, xhi_negative = system.plane_surface(1311, 1, 0, 0, -1001)
    _, ylo_positive, _ = system.plane_surface(1312, 0, 1, 0, -1000)
    _, _, yhi_negative = system.plane_surface(
        1313, 0, 1, 0, -(1000 + 1e-12))
    _, zlo_positive, _ = system.plane_surface(1314, 0, 0, 1, 1)
    _, _, zhi_negative = system.plane_surface(1315, 0, 0, 1, -1)
    narrow = system.create_intersection(ylo_positive, yhi_negative)
    broad = system.create_union(narrow, xhi_negative)
    system.add_cell(1310, system.create_intersection_many([
        broad, xlo_positive, zlo_positive, zhi_negative]))

    page = system.slice_error_page(
        (0, 0, 0), (0, 0, 1), (0, 1, 0),
        (999, 1001, 999, 1001), (999, 1001, 999, 1001),
        options={"critical_relative_distance_tolerance": 1e-12})
    assert page["receipt"]["scope_classified"] is False
    assert page["receipt"]["scan_stop_reason"] == "numerical_unresolved"
    assert page["receipt"]["close_crossing_observations"] > 0
    assert page["receipt"]["symbolic_one_sided_intervals"] > 0
    assert page["receipt"]["numerical_unrepresentable_probe_intervals"] > 0
    assert page["regions"]
    assert all(region["kind"] == "unresolved" and
               region["numerical_cause"] == "unrepresentable_probes"
               for region in page["regions"])


def test_slice_error_proximity_tolerance_does_not_reject_a_thin_gap():
    system = pyalea.System()
    _, _, left = system.plane_surface(1316, 1, 0, 0, 0)
    _, right, _ = system.plane_surface(1317, 1, 0, 0, -1e-7)
    system.add_cell(1316, left)
    system.add_cell(1317, right)
    args = ((0, 0, 0), (0, 0, 1), (0, 1, 0),
            (-0.1, 0.1, -0.1, 0.1), (-0.1, 0.1, -0.1, 0.1))

    # The threshold records proximity but does not erase the distinct,
    # representable interval or turn it into a numerical failure.
    default = system.slice_error_page(*args)
    assert default["receipt"]["scope_classified"] is True
    assert default["receipt"]["unresolved_reason"] == "resolved"
    assert default["regions"][0]["kind"] == "gap"

    resolved = system.slice_error_page(
        *args, options={"critical_relative_distance_tolerance": 1e-9})
    assert resolved["receipt"]["scope_classified"] is True
    assert resolved["receipt"]["unresolved_reason"] == "resolved"
    assert resolved["regions"][0]["kind"] == "gap"
    assert resolved["regions"][0]["uv_max"][0] - (
        resolved["regions"][0]["uv_min"][0]) == pytest.approx(1e-7)

    with pytest.raises(ValueError):
        system.slice_error_page(
            *args,
            options={"critical_relative_distance_tolerance": -1e-6})


def test_slice_error_thin_valid_layer_and_overlap_keep_their_ownership():
    width = 1e-7
    args = ((0, 0, 0), (0, 0, 1), (0, 1, 0),
            (-1, 1, -1, 1), (-1, 1, -1, 1))

    valid = pyalea.System()
    _, low_positive, low_negative = valid.plane_surface(
        1320, 1, 0, 0, -0.2)
    _, high_positive, high_negative = valid.plane_surface(
        1321, 1, 0, 0, -(0.2 + width))
    valid.add_cell(1320, low_negative)
    valid.add_cell(1321, valid.create_intersection(
        low_positive, high_negative))
    valid.add_cell(1322, high_positive)
    for columns in (1, 2):
        with valid.slice_error_query(
                *args, tile_columns=columns) as query:
            pages = [query.run_page(i) for i in range(query.page_count)]
        assert all(page["receipt"]["scope_classified"] for page in pages)
        assert not any(page["regions"] for page in pages)

    overlap = pyalea.System()
    _, low_positive, _ = overlap.plane_surface(1330, 1, 0, 0, -0.2)
    _, _, high_negative = overlap.plane_surface(
        1331, 1, 0, 0, -(0.2 + width))
    overlap.add_cell(1330, high_negative)
    overlap.add_cell(1331, low_positive)
    page = overlap.slice_error_page(*args)
    assert page["receipt"]["scope_classified"] is True
    assert page["regions"][0]["kind"] == "overlap"
    assert page["regions"][0]["uv_max"][0] - (
        page["regions"][0]["uv_min"][0]) == pytest.approx(width)


def test_slice_error_page_keeps_plane_boundary_with_separated_sphere():
    system = pyalea.System()
    _, _, sphere_inside = system.sphere_surface(1060, 0, 0, 0, 10)
    _, _, plane_negative = system.plane_surface(1061, 1, 0, 0, 0)
    system.add_cell(1060, system.create_intersection(
        sphere_inside, plane_negative))
    page = system.slice_error_page(
        (0, 0, 0), (0, 0, 1), (0, 1, 0),
        (-1, 1, -1, 1), (-1, 1, -1, 1))
    assert page["receipt"]["scope_classified"] is True
    assert page["unresolved"] == []
    assert len(page["intervals"]) == 1
    assert page["intervals"][0]["surface_id"] == 1061
    assert page["regions"][0]["kind"] == "gap"


def test_slice_error_page_reports_verified_isolated_circle():
    system = pyalea.System()
    _, _, inside = system.sphere_surface(1063, 0, 0, 0, 0.5)
    system.add_cell(1063, inside)
    page = system.slice_error_page(
        (0, 0, 0), (0, 0, 1), (0, 1, 0),
        (-1, 1, -1, 1), (-1, 1, -1, 1))
    assert page["receipt"]["scope_classified"] is True
    assert page["receipt"]["verified_circle_count"] == 1
    assert page["intervals"] == []
    assert page["regions"] == []
    circle, = page["circles"]
    assert circle["surface_id"] == 1063
    assert circle["center_uv"] == (0.0, 0.0)
    assert circle["radius"] == pytest.approx(0.5)
    assert circle["inside_kind"] == "unique"
    assert circle["outside_kind"] == "gap"


def test_slice_error_page_reports_transverse_circle_arcs():
    system = pyalea.System()
    _, _, left = system.sphere_surface(1070, -0.3, 0, 0, 0.6)
    _, _, right = system.sphere_surface(1071, 0.3, 0, 0, 0.6)
    system.add_cell(1070, left)
    system.add_cell(1071, right)
    page = system.slice_error_page(
        (0, 0, 0), (0, 0, 1), (0, 1, 0),
        (-1, 1, -1, 1), (-1, 1, -1, 1))
    assert page["receipt"]["scope_classified"] is True
    assert page["receipt"]["verified_circle_count"] == 4
    assert sum(arc["inside_kind"] == "overlap"
               for arc in page["circles"]) == 2
    assert sum(arc["outside_kind"] == "gap"
               for arc in page["circles"]) == 2
    assert all(arc["start_angle"] < arc["end_angle"]
               for arc in page["circles"])


def test_slice_error_page_reports_mixed_plane_circle_boundaries():
    system = pyalea.System()
    _, _, plane_negative = system.plane_surface(1074, 1, 0, 0, 0)
    _, _, sphere_inside = system.sphere_surface(1075, 0, 0, 0, 0.6)
    system.add_cell(1074, plane_negative)
    system.add_cell(1075, sphere_inside)
    page = system.slice_error_page(
        (0, 0, 0), (0, 0, 1), (0, 1, 0),
        (-1, 1, -1, 1), (-1, 1, -1, 1))
    assert page["receipt"]["scope_classified"] is True
    assert page["receipt"]["verified_interval_count"] == 3
    assert page["receipt"]["verified_circle_count"] == 2
    assert {interval["surface_id"] for interval in page["intervals"]} == {1074}
    assert {arc["surface_id"] for arc in page["circles"]} == {1075}
    assert {arc["inside_kind"] for arc in page["circles"]} == {
        "unique", "overlap"}


def test_slice_error_page_reports_nested_circle_faces():
    system = pyalea.System()
    _, _, outer = system.sphere_surface(1072, 0, 0, 0, 0.8)
    _, _, inner = system.sphere_surface(1073, 0, 0, 0, 0.3)
    system.add_cell(1072, outer)
    system.add_cell(1073, inner)
    page = system.slice_error_page(
        (0, 0, 0), (0, 0, 1), (0, 1, 0),
        (-2, 2, -2, 2), (-2, 2, -2, 2))
    assert page["receipt"]["scope_classified"] is True
    assert page["receipt"]["verified_circle_count"] == 2
    assert [arc["inside_kind"] for arc in page["circles"]] == [
        "unique", "overlap"]
    assert [arc["outside_kind"] for arc in page["circles"]] == [
        "gap", "unique"]


def test_slice_error_page_exposes_three_line_gap():
    system = pyalea.System()
    _, x_positive, x_negative = system.plane_surface(1170, 1, 0, 0, 0)
    _, y_positive, y_negative = system.plane_surface(1171, 0, 1, 0, 0)
    _, diagonal_positive, _ = system.plane_surface(1172, 1, 1, 0, -0.8)
    system.add_cell(1170, x_negative)
    system.add_cell(1171, system.create_intersection(x_positive, y_negative))
    system.add_cell(1172, system.create_intersection(
        system.create_intersection(x_positive, y_positive),
        diagonal_positive))
    page = system.slice_error_page(
        (0, 0, 0), (0, 0, 1), (0, 1, 0),
        (-1, 1, -1, 1), (-1, 1, -1, 1))
    assert page["receipt"]["scope_classified"] is True
    assert len(page["regions"]) == 1
    assert page["regions"][0]["kind"] == "gap"
    assert len(page["intervals"]) == 3
    assert {interval["surface_id"] for interval in page["intervals"]} == {
        1170, 1171, 1172}


def test_slice_error_page_exposes_verified_oblique_polygon():
    system = pyalea.System()
    _, _, negative = system.plane_surface(1050, 1.0, 1.0, 0.0, -0.2)
    system.add_cell(1050, negative)
    page = system.slice_error_page(
        (0, 0, 0), (0, 0, 1), (0, 1, 0),
        (-1, 1, -1, 1), (-1, 1, -1, 1),
    )
    assert page["receipt"]["scope_classified"] is True
    assert page["receipt"]["region_count"] == 1
    assert page["receipt"]["verified_interval_count"] == 1
    region = page["regions"][0]
    assert region["kind"] == "gap"
    assert len(region["polygon_uv"]) == 3
    assert len(region["polygon_uv_uncertainty"]) == 3
    assert region["uv_min"] == pytest.approx((-0.8, -0.8))
    assert region["uv_max"] == pytest.approx((1.0, 1.0))
    interval = page["intervals"][0]
    assert interval["axis"] == -1
    assert interval["negative_side_kind"] == "unique"
    assert interval["positive_side_kind"] == "gap"

    slab = pyalea.System()
    _, _, left = slab.plane_surface(1060, 1.0, 1.0, 0.0, 0.2)
    _, right, _ = slab.plane_surface(1061, 1.0, 1.0, 0.0, -0.2)
    slab.add_cell(1060, left)
    slab.add_cell(1061, right)
    band = slab.slice_error_page(
        (0, 0, 0), (0, 0, 1), (0, 1, 0),
        (-1, 1, -1, 1), (-1, 1, -1, 1),
    )
    assert band["receipt"]["scope_classified"] is True
    assert band["receipt"]["region_count"] == 1
    assert band["receipt"]["verified_interval_count"] == 2
    assert band["regions"][0]["kind"] == "gap"
    assert len(band["regions"][0]["polygon_uv"]) == 6

    crossing = pyalea.System()
    _, first_pos, first_neg = crossing.plane_surface(
        1070, 1.0, 1.0, 0.0, -0.2)
    _, second_pos, second_neg = crossing.plane_surface(
        1071, 1.0, -1.0, 0.0, 0.1)
    crossing.add_cell(1070, crossing.create_intersection(
        first_neg, second_neg))
    crossing.add_cell(1071, crossing.create_intersection(
        first_pos, second_neg))
    crossing.add_cell(1072, crossing.create_intersection(
        first_neg, second_pos))
    quadrant = crossing.slice_error_page(
        (0, 0, 0), (0, 0, 1), (0, 1, 0),
        (-1, 1, -1, 1), (-1, 1, -1, 1),
    )
    assert quadrant["receipt"]["scope_classified"] is True
    assert quadrant["receipt"]["region_count"] == 1
    assert quadrant["receipt"]["verified_interval_count"] == 2
    assert len(quadrant["regions"][0]["polygon_uv"]) == 4


def test_slice_error_page_classifies_axis_and_oblique_crossing():
    system = pyalea.System()
    _, axis_pos, axis_neg = system.plane_surface(
        1090, 1.0, 0.0, 0.0, -0.2)
    _, oblique_pos, oblique_neg = system.plane_surface(
        1091, 1.0, 1.0, 0.0, -0.3)
    system.add_cell(1090, system.create_intersection(axis_neg, oblique_neg))
    system.add_cell(1091, system.create_intersection(axis_pos, oblique_neg))
    system.add_cell(1092, system.create_intersection(axis_neg, oblique_pos))
    with system.slice_error_query(
        (0, 0, 0), (0, 0, 1), (0, 1, 0),
        (-1, 1, -1, 1), (-1, 1, -1, 1),
        tile_columns=2, tile_rows=2,
    ) as query:
        assert query.page_count == 4
        assert all(query.run_page(i)["receipt"]["scope_classified"]
                   for i in range(query.page_count))
    page = system.slice_error_page(
        (0, 0, 0), (0, 0, 1), (0, 1, 0),
        (-1, 1, -1, 1), (-1, 1, -1, 1),
    )
    assert page["receipt"]["scope_classified"] is True
    assert page["receipt"]["region_count"] == 1
    assert page["receipt"]["verified_interval_count"] == 2
    assert page["regions"][0]["kind"] == "gap"
    assert len(page["regions"][0]["polygon_uv"]) >= 3


def test_slice_error_query_reuses_identity_and_survives_partial_pages():
    system = pyalea.System()
    _, _, negative = system.plane_surface(20, 1.0, 0.0, 0.0, 0.0)
    system.add_cell(20, negative)
    args = ((0, 0, 0), (0, 0, 1), (0, 1, 0),
            (-1, 1, -1, 1), (-1, 1, -1, 1))
    query = system.slice_error_query(*args, tile_columns=2, tile_rows=2)
    assert query.page_count == 4
    pages = [query.run_page(i) for i in (2, 0, 3, 1)]
    assert {page["receipt"]["page_index"] for page in pages} == {0, 1, 2, 3}
    assert len({page["receipt"]["query_id"] for page in pages}) == 1
    assert all(page["receipt"]["scope_classified"] for page in pages)
    repeated = query.run_page(0)
    for receipt in (repeated["receipt"], pages[1]["receipt"]):
        receipt.pop("confirmation_seconds")
        receipt.pop("elapsed_seconds")
    assert repeated == pages[1]
    with pytest.raises(IndexError):
        query.run_page(4)
    with pytest.raises(IndexError):
        query.run_page(-1)
    query.close()
    query.close()
    with pytest.raises(RuntimeError):
        query.run_page(0)


def test_slice_error_query_retains_system_and_rejects_replacement():
    system = pyalea.System()
    args = ((0, 0, 0), (0, 0, 1), (0, 1, 0),
            (-1, 1, -1, 1), (-1, 1, -1, 1))
    query = system.slice_error_query(*args)
    del system
    gc.collect()
    assert query.run_page(0)["receipt"]["scope_classified"] is True
    query.close()

    system = pyalea.System()
    query = system.slice_error_query(*args)
    system.__init__()
    with pytest.raises(RuntimeError, match="replaced"):
        query.run_page(0)
    query.close()

    query = system.slice_error_query(*args)
    system.plane_surface(21, 1.0, 0.0, 0.0, 0.0)
    with pytest.raises(RuntimeError):
        query.run_page(0)
    query.close()

    with system.slice_error_query(*args) as query:
        assert query.run_page(0)["receipt"]["scope_classified"] is True
    with pytest.raises(RuntimeError):
        query.run_page(0)


@pytest.mark.parametrize(
    ("primitive_type", "parameters", "inside", "outside"),
    [
        (pyalea.PRIMITIVE_SPHERE, (0, 0, 0, 2), (0, 0, 0), (3, 0, 0)),
        (pyalea.PRIMITIVE_BOX,
         (0, 0, 0, 2, 0, 0, 0, 3, 0, 0, 0, 4),
         (1, 1, 1), (3, 1, 1)),
        (pyalea.PRIMITIVE_REC,
         (0, 0, 0, 0, 0, 4, 2, 0, 0, 0, 1, 0),
         (0, 0, 2), (3, 0, 2)),
        (pyalea.PRIMITIVE_WED,
         (0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 4),
         (0.25, 0.25, 2), (2, 2, 2)),
        (pyalea.PRIMITIVE_RHP,
         (0, 0, 0, 0, 0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0),
         (0, 0, 2), (3, 0, 2)),
    ],
)
def test_standalone_primitive_evaluation(
        primitive_type, parameters, inside, outside):
    assert pyalea.primitive_evaluate(primitive_type, parameters, inside) < 0
    assert pyalea.primitive_evaluate(primitive_type, parameters, outside) > 0


def test_standalone_primitive_evaluation_validates_inputs():
    with pytest.raises(ValueError, match="exactly 4"):
        pyalea.primitive_evaluate(pyalea.PRIMITIVE_SPHERE, (0, 0, 0), (0, 0, 0))
    with pytest.raises(ValueError, match="invalid or unsupported"):
        pyalea.primitive_evaluate(pyalea.PRIMITIVE_SPHERE, (0, 0, 0, 0), (0, 0, 0))
    with pytest.raises(ValueError, match="finite"):
        pyalea.primitive_evaluate(
            pyalea.PRIMITIVE_SPHERE, (0, 0, 0, 2), (float("nan"), 0, 0))


def test_parallel_runtime_value_stays_synchronized():
    original = pyalea.parallel_max_threads()
    try:
        assert pyalea.set_parallel_threads(2) == 2
        assert pyalea.PARALLEL_MAX_THREADS == 2
        assert pyalea._alea.PARALLEL_MAX_THREADS == 2
    finally:
        pyalea.set_parallel_threads(original)


def test_slice_error_query_runs_scratch_bounded_parallel_pages():
    system = pyalea.System()
    _, _, negative = system.plane_surface(22, 1.0, 0.0, 0.0, 0.0)
    system.add_cell(22, negative)
    args = ((0, 0, 0), (0, 0, 1), (0, 1, 0),
            (-2, 2, -1, 1), (-2, 2, -1, 1))
    with system.slice_error_query(
            *args, tile_columns=4, tile_rows=1) as query:
        batch = query.run_pages(
            [3, 1, 0, 2], workers=4,
            max_parallel_scratch_bytes=8 * 1024 * 1024)
        assert [page["receipt"]["page_index"]
                for page in batch["pages"]] == [3, 1, 0, 2]
        assert batch["stats"]["page_count"] == 4
        assert batch["stats"]["completed_page_count"] == 4
        assert 1 <= batch["stats"]["actual_workers"] <= 2
        assert batch["stats"]["reserved_parallel_scratch_bytes"] <= (
            8 * 1024 * 1024)
        serial = query.run_pages(
            [0, 1], workers=4, max_parallel_scratch_bytes=0)
        assert serial["stats"]["actual_workers"] == 1
        with pytest.raises(IndexError):
            query.run_pages([4])
        with pytest.raises(TypeError):
            query.run_pages([True])
        with pytest.raises(ValueError):
            query.run_pages([0], workers=-1)


def test_packaged_build_metadata_when_present():
    metadata_path = Path(pyalea.__file__).with_name("_build_info.json")
    if not metadata_path.exists():
        pytest.skip("build metadata is generated by the package target")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    assert metadata["libalea_version"] == pyalea.version()
    assert len(metadata["libalea_commit"]) == 40
    assert set(metadata["libalea_commit"]) <= set("0123456789abcdef")
    assert Path(pyalea.__file__).with_name("__init__.pyi").is_file()
    assert Path(pyalea.__file__).with_name("py.typed").is_file()


def test_numpy_grid_arrays_have_stable_layout_and_dtypes(populated_system):
    grid = populated_system.find_cells_grid_z(
        0.0, -3.0, 3.0, -3.0, 3.0, 9, 7,
        error_mode="fast", _as_buffers=True,
    )

    expected = {
        "cell_ids": np.dtype(np.intc),
        "material_ids": np.dtype(np.intc),
        "secondary_cell_ids": np.dtype(np.intc),
        "coverage": np.dtype(np.uint8),
        "errors": np.dtype(np.uint8),
    }
    for name, dtype in expected.items():
        array = grid[name]
        assert isinstance(array, np.ndarray)
        assert array.dtype == dtype
        assert array.shape == (63,)
        assert array.flags.c_contiguous
        assert array.flags.writeable

    assert 101 in grid["cell_ids"]
    assert 7 in grid["material_ids"]


def test_numpy_arrays_own_native_storage_after_system_destruction(populated_system):
    grid = populated_system.find_cells_grid_z(
        0.0, -3.0, 3.0, -3.0, 3.0, 8, 6, _as_buffers=True,
    )
    cells = grid["cell_ids"]
    expected = cells.copy()
    assert cells.base is not None

    del grid
    del populated_system
    gc.collect()

    np.testing.assert_array_equal(cells, expected)
    cells[0] = 12345
    assert cells[0] == 12345


def test_numpy_grid_releases_all_arrays_with_result(populated_system):
    grid = populated_system.find_cells_grid_z(
        0.0, -3.0, 3.0, -3.0, 3.0, 8, 6,
        error_mode="fast", _as_buffers=True,
    )
    references = [weakref.ref(value) for value in grid.values()
                  if isinstance(value, np.ndarray)]
    assert references

    del grid
    gc.collect()

    assert all(reference() is None for reference in references)


def test_reinitializing_owned_types_replaces_native_state():
    system = pyalea.System()
    system.__init__()
    assert system.cell_count == 0

    material = pyalea.NucMaterial()
    material.__init__()

    multigroup = pyalea.Multigroup([10.0, 1.0, 0.1])
    multigroup.__init__([20.0, 2.0])
    assert multigroup.n_groups == 1


def transformed_fill_system(transform_factory):
    system = pyalea.System()
    _, _, sphere = system.sphere_surface(1, 0.0, 0.0, 0.0, 1.5)
    _, _, container_region = system.box_surface(
        2, -50.0, 50.0, -50.0, 50.0, -50.0, 50.0,
    )
    material = system.add_material(7)
    system.add_cell(10, sphere, material, -1.0, universe_id=1)
    container = system.add_cell(1, container_region)
    transform_id = transform_factory(system)
    system.set_fill(container, 1, transform_id)
    system.build_universe_index()
    return system, transform_id


def test_named_transform_positions_programmatic_fill():
    system, transform_id = transformed_fill_system(
        lambda system: system.add_transform(41, (10.0, 0.0, 0.0)),
    )

    assert transform_id == 41
    assert system.material_at(10.0, 0.0, 0.0) == 7

    assert system.add_transform(41, (-10.0, 0.0, 0.0)) == 41
    assert system.material_at(-10.0, 0.0, 0.0) == 7


def test_inline_transform_deduplicates_and_positions_programmatic_fill():
    def add_inline(system):
        first = system.add_inline_transform(
            (10.0, 0.0, 0.0), cell_id=1, role="fill",
        )
        duplicate = system.add_inline_transform(
            (10.0, 0.0, 0.0), cell_id=2, role="fill",
        )
        assert duplicate == first
        return first

    system, transform_id = transformed_fill_system(add_inline)

    assert transform_id > 0
    assert system.material_at(10.0, 0.0, 0.0) == 7


def test_transform_values_are_validated():
    system = pyalea.System()
    with pytest.raises(ValueError, match="between 3 and 13"):
        system.add_transform(1, (1.0, 2.0))
    with pytest.raises(ValueError, match="rotation"):
        system.add_inline_transform((0.0, 0.0, 0.0, 0.5))


def test_in_memory_exports(populated_system):
    mcnp = populated_system.export_mcnp_string()
    openmc = populated_system.export_openmc_string()
    serpent = populated_system.export_serpent_string()

    assert "Exported by CSG Library" in mcnp
    assert "<geometry>" in openmc
    assert "exported by Alea" in serpent


def test_compact_validation_reports_trace_reuse(populated_system):
    result = populated_system.validate_ray_slice_compact(
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 1.0),
        (0.0, 1.0, 0.0),
        -3.0, 3.0, -3.0, 3.0, 8,
    )

    assert isinstance(result["executed_trace_mask"], int)
    assert isinstance(result["reused_trace_mask"], int)
    assert result["executed_trace_mask"] != 0


def test_numpy_buffers_survive_repeated_allocation_and_collection():
    retained = []
    for size in range(2, 10):
        system = pyalea.System()
        _, _, inside = system.sphere_surface(1, 0.0, 0.0, 0.0, 2.0)
        material = system.add_material(7)
        system.add_cell(101, inside, material, 1.0)
        system.build_universe_index()
        grid = system.find_cells_grid_z(
            0.0, -2.5, 2.5, -2.5, 2.5, size, size,
            _as_buffers=True,
        )
        materials = grid["material_ids"]
        retained.append((materials, materials.copy()))
        del grid
        del system
    gc.collect()

    for size, (array, expected) in zip(range(2, 10), retained):
        assert array.shape == (size * size,)
        np.testing.assert_array_equal(array, expected)


def test_numpy_grid_rejects_invalid_dimensions(populated_system):
    with pytest.raises(ValueError):
        populated_system.find_cells_grid_z(
            0.0, -1.0, 1.0, -1.0, 1.0, 0, 4,
            _as_buffers=True,
        )
