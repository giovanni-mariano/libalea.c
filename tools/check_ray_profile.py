#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0
"""Check ray profiling against ordinary replay on a public generated fixture."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile

from generate_lifecycle_fixture import generate


def main():
    normal = Path(sys.argv[1] if len(sys.argv) > 1 else "bin/model_lifecycle").resolve()
    profiled = Path(sys.argv[2] if len(sys.argv) > 2 else
                    "bin/ray-profile/model_lifecycle").resolve()
    with tempfile.TemporaryDirectory() as directory:
        temp = Path(directory)
        model, workload = temp / "model.i", temp / "queries.txt"
        generate(model, workload, 16)
        arguments = ["fresh", str(model), str(workload), "-", "--hierarchical-rays"]

        def run(binary, *options):
            result = subprocess.run([str(binary), *arguments, *options],
                                    check=True, capture_output=True, text=True)
            return json.loads(result.stdout)

        ordinary = run(normal)
        disabled = run(profiled)
        enabled = run(profiled, "--ray-profile")
        assert "ray_profile" not in ordinary and "ray_profile" not in disabled
        for key in ("point_digest", "ray_digest", "acceleration_digest", "segments"):
            assert ordinary[key] == disabled[key] == enabled[key], key
        costs = enabled["ray_profile"]
        stages = {s["name"]: s for s in costs["stages"]}
        assert set(stages) == {"other", "hierarchy", "membership", "intersection",
                               "recovery", "recovery_scan", "recovery_cache"}
        assert all(s["seconds"] >= 0 and s["calls"] >= 0 for s in stages.values())
        # Start/stop enclose the replay clock by a few clock calls.
        assert abs(sum(s["seconds"] for s in stages.values()) - enabled["ray_s"]) < .01
        assert stages["membership"]["calls"] > 0
        assert costs["csg_nodes"] >= costs["csg_primitives"] > 0
        assert costs["dropped_root_samples"] == 0
        assert 0 < len(costs["top_roots"]) <= 32
        assert all(r["calls"] > 0 and r["seconds"] >= 0 for r in costs["top_roots"])
        assert costs["dropped_list_samples"] == 0
        assert len(costs["scan_types"]) == 32
        assert all(t["tests"] >= 0 and t["window_roots"] >= 0 for t in costs["scan_types"])
        events = costs["owner_probe_recoveries"]
        assert len(events) + costs["dropped_recoveries"] == stages["recovery"]["calls"]
        for e in events:
            assert 0 <= e["ray"] < enabled["rays"]
            assert 1 <= e["attempt"] <= 2
            assert e["source"] in {"cached_path", "neighbor", "path_descent", "root", "full_lookup"}
            assert e["current_t"] < e["proposed_t"] and e["current_t"] < e["recovered_t"]
            assert len(e["occurrence"]) == 16
        fills = costs["fill_traversal"]
        assert fills["nearer_hits"] <= fills["nearest_passes"]
        assert sum(fills[k] for k in ("duplicate_skips", "bbox_rejections", "descents")) <= fills["child_visits"]
        assert sum(t["tests"] for t in costs["root_build_types"]) == costs["recovery_root_tests"]
        assert costs["cached_non_torus_candidates"] <= costs["recovery_root_queries"]
        assert costs["indexed_root_lists"] <= costs["recovery_root_builds"]
        assert costs["roots_sorted"] <= sum(t["window_roots"] for t in costs["root_build_types"])
        lists = costs["top_scan_lists"]
        assert len(lists) <= 32
        assert sum(r["tests"] for r in lists) <= sum(t["tests"] for t in costs["scan_types"])
        assert all(r["calls"] > 0 and r["tests"] <= r["calls"] * r["surface_count"] for r in lists)
        rejected = subprocess.run([str(normal), *arguments, "--ray-profile"],
                                  capture_output=True, text=True)
        assert rejected.returncode == 2 and "ALEA_RAY_PROFILE" in rejected.stderr

        # Exercise the Python driver's forwarding and report metadata too.
        report = temp / "report.json"
        subprocess.run([sys.executable, str(Path(__file__).with_name("model_lifecycle_benchmark.py")),
                        str(model), "--binary", str(profiled), "--workload", str(workload),
                        "--output", str(report), "--repeats", "1", "--ray-profile"],
                       check=True, stdout=subprocess.PIPE)
        data = json.loads(report.read_text())
        assert data["complete"] and data["ray_profile"]
        assert all("ray_profile" in record for record in data["records"])
    print("Ray profiling replay, accounting, and driver checks passed")


if __name__ == "__main__":
    main()
