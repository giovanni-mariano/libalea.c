# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0

import importlib.util
import sys
from pathlib import Path


def _load_benchmark():
    path = (Path(__file__).parents[3] / "tools" /
            "geometry_error_model_benchmark.py")
    spec = importlib.util.spec_from_file_location(
        "geometry_error_model_benchmark", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_model_benchmark_emits_reproducible_receipt_fields():
    module = _load_benchmark()
    model = Path(__file__).parents[3] / "tests" / "data" / "simple_box.mcnp"
    report = module.run_benchmark(
        model, value=0.0, bounds=(-2.0, 2.0, -2.0, 2.0),
        tile_columns=2, tile_rows=1, pages="0:2",
        interior_probes_per_axis=1, max_index_bytes=1024 * 1024)
    assert report["schema"] == "libalea.geometry-error-model-benchmark.v1"
    assert report["metadata"]["model_sha256"]
    assert report["workload"]["selected_pages"] == [0, 1]
    assert report["summary"]["completed_page_count"] == 2
    assert report["summary"]["interior_probe_count"] == 2
    assert len(report["pages"]) == 2
    assert all("receipt" in page for page in report["pages"])
