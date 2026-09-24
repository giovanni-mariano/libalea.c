# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0

import importlib.util
import sys
from pathlib import Path


def _load_harness():
    path = Path(__file__).parents[3] / "tools" / "geometry_error_acceptance.py"
    spec = importlib.util.spec_from_file_location(
        "geometry_error_acceptance", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_known_answer_geometry_error_acceptance():
    report = _load_harness().run_acceptance()
    assert report["schema"] == "libalea.geometry-error-acceptance.v1"
    assert report["summary"]["passed"] is True
    assert report["summary"]["fixture_count"] == 17
    assert report["summary"]["expected_defect_count"] == 14
    assert report["summary"]["detected_defect_count"] == 14
    assert report["summary"]["missed_defect_count"] == 0
    assert report["summary"]["false_positive_count"] == 0
    assert report["summary"]["invalid_evidence_count"] == 0
    cases = {case["name"]: case for case in report["cases"]}
    assert cases["adjacent_clean"]["clean_confirmed_record_count"] == 0
    assert cases["viewport_inside_gap"]["witness_count"] > 0
    assert cases["repeated_fill_overlap"]["witnesses"][0][
        "owners_complete"] is True
    assert any(witness["source"] == "interior_probe" and
               witness["uv"][0] != 0.0
               for witness in cases["off_center_enclosed_overlap"]["witnesses"])
    assert cases["rotated_fill_clean"]["clean_confirmed_record_count"] == 0
    assert cases["valid_rectangular_lattice_seam"][
        "clean_confirmed_record_count"] == 0
    mixed = cases["defect_beside_unsupported_macrobody"]
    assert mixed["matched_defect_ids"] == [
        "overlap-beside-unsupported-macrobody"]
    assert mixed["unresolved_page_count"] > 0
    assert cases["box_macrobody_exterior_gap"]["unresolved_page_count"] > 0
    assert cases["rcc_macrobody_exterior_gap"]["unresolved_page_count"] > 0
    assert cases["cylinder_exterior_gap"]["unresolved_page_count"] > 0
