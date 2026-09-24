# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0

import importlib.util
import sys
from pathlib import Path

import pytest
import pyalea


def _load_harness():
    path = (Path(__file__).parents[3] / "tools" /
            "geometry_error_seeded_benchmark.py")
    spec = importlib.util.spec_from_file_location(
        "geometry_error_seeded_benchmark", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _partitioned_system():
    system = pyalea.System()
    _, positive, negative = system.plane_surface(1700, 1, 0, 0, 0)
    system.add_cell(1701, negative)
    system.add_cell(1702, positive)
    return system


@pytest.mark.parametrize("mutation", ["gap", "overlap"])
def test_seeded_benchmark_detects_programmatic_point_defect(mutation):
    report = _load_harness().run_seeded_system(
        _partitioned_system(), mutation=mutation, value=0.0,
        bounds=(-1.0, 1.0, -0.5, 0.5), tile_columns=1, tile_rows=1,
        page_index=0, probes_per_axis=3, max_index_bytes=1024 * 1024)
    assert report["schema"] == "libalea.geometry-error-seeded-benchmark.v1"
    assert report["oracle"]["baseline_observed_kind"] == "unique"
    assert report["oracle"]["mutated_observed_kind"] == mutation
    assert report["result"]["point_mutation_verified"] is True
    assert report["result"]["seed_witness_detected"] is True
    assert report["result"]["matching_seed_witnesses"]
    assert report["mutation"]["inserted_surface_id"] <= 2
