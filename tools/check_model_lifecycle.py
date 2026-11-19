#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0
"""Exercise replay gates, cache hits, edited-model isolation, and input rejection."""
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from generate_lifecycle_fixture import generate, generate_quadric


def main():
    binary = Path(sys.argv[1] if len(sys.argv) > 1 else "bin/model_lifecycle").resolve()
    runner = Path(__file__).with_name("model_lifecycle_benchmark.py")
    with tempfile.TemporaryDirectory() as temp:
        temp = Path(temp)
        model, workload = temp / "model.i", temp / "queries.txt"
        edited, report = temp / "edited.i", temp / "report.json"
        generate(model, workload, 16)
        edited.write_text(model.read_text().replace("3 p -1 0 1 3\n", "3 p -1 0 1 3.25\n"))
        command = [sys.executable, str(runner), str(model), "--binary", str(binary),
                   "--workload", str(workload), "--output", str(report), "--repeats", "1",
                   "--mixed-bounds", "--edited-model", str(edited)]
        subprocess.run(command, check=True, stdout=subprocess.PIPE)
        data = json.loads(report.read_text())
        assert data["complete"]
        rows = {r["scenario"]: r for r in data["records"]}
        assert rows["warm-cache"]["cache_hit"] == 1
        assert rows["mixed-bounds"]["bounds_changed"] == 16
        assert rows["mixed-bounds"]["point_candidates"] < rows["fresh"]["point_candidates"]
        assert rows["edited-fresh"]["ray_digest"] != rows["fresh"]["ray_digest"]
        assert rows["edit-with-previous"]["ray_digest"] == rows["edited-fresh"]["ray_digest"]
        good_baseline = temp / "baseline.json"
        good_baseline.write_text(json.dumps(data))
        subprocess.run(command + ["--candidate-profile", "--baseline", str(good_baseline)],
                       check=True, stdout=subprocess.PIPE)
        profiled = json.loads(report.read_text())
        assert profiled["candidate_profile"]
        for record in profiled["records"]:
            cells = record["candidate_profile"]
            assert sum(c["point_candidates"] for c in cells) == record["point_candidates"]
            assert sum(c["ray_candidates"] for c in cells) == record["ray_point_candidates"]
            assert len({c["cell_index"] for c in cells}) == len(cells)
            assert all(c["primitive_mask"] == (1 << 1) | (1 << 5) for c in cells)
        subprocess.run(command + ["--snapshot-io", "mapped", "--baseline", str(good_baseline)],
                       check=True, stdout=subprocess.PIPE)
        mapped = json.loads(report.read_text())
        assert all(r["snapshot_io_used"] == ["mapped"] for r in mapped["records"]
                   if r["scenario"] in ("warm-cache", "snapshot", "edit-with-previous"))
        subprocess.run(command + ["--snapshot-io", "both"], check=True, stdout=subprocess.PIPE)
        both = json.loads(report.read_text())
        rows = {r["scenario"]: r for r in both["records"]}
        assert rows["warm-cache-mapped"]["cache_hit"] == 1
        assert rows["snapshot-mapped"]["snapshot_io_used"] == ["mapped"]
        assert rows["snapshot"]["snapshot_io_used"] == ["buffered"]
        original = model.read_bytes()
        failed = subprocess.run(command + ["--output", str(model)], capture_output=True)
        assert failed.returncode == 2 and model.read_bytes() == original
        gate = json.loads(good_baseline.read_text())
        for metrics in gate["summary"].values():
            metrics["ready_s"]["median"] = 1e-100
        impossible = temp / "impossible-gate.json"
        impossible.write_text(json.dumps(gate))
        failed = subprocess.run(command + ["--baseline", str(impossible), "--max-ready-ratio", "1"],
                                capture_output=True)
        assert failed.returncode != 0
        assert "performance gate failed" in json.loads(report.read_text())["error"]
        # Query samples alone cannot certify exact adjacency/surface ordering.
        structural = json.loads(good_baseline.read_text())
        structural["records"][0]["acceleration_digest"] = "bad"
        structural_baseline = temp / "bad-acceleration.json"
        structural_baseline.write_text(json.dumps(structural))
        failed = subprocess.run(command + ["--baseline", str(structural_baseline)], capture_output=True)
        assert failed.returncode != 0
        assert "baseline replay mismatch" in json.loads(report.read_text())["error"]
        # A baseline with altered results must never pass a speed-only gate.
        data["records"][0]["point_digest"] = "bad"
        baseline = temp / "bad-baseline.json"
        baseline.write_text(json.dumps(data))
        failed = subprocess.run(command + ["--baseline", str(baseline)], capture_output=True)
        assert failed.returncode != 0
        assert "baseline replay mismatch" in json.loads(report.read_text())["error"]
        for text in ("p 1 2 3 garbage\n", "r 0 0 0 0 0 0 5\n", "p nan 1 2\n"):
            workload.write_text(text)
            failed = subprocess.run([str(binary), "fresh", str(model), str(workload), "-"],
                                    capture_output=True)
            assert failed.returncode == 2
        generate_quadric(model, workload, 16)
        subprocess.run([sys.executable, str(runner), str(model), "--binary", str(binary),
                        "--workload", str(workload), "--output", str(report),
                        "--repeats", "1", "--mixed-bounds", "--candidate-profile"],
                       check=True, stdout=subprocess.PIPE)
        quadrics = json.loads(report.read_text())
        assert quadrics["complete"]
        rows = {r["scenario"]: r for r in quadrics["records"]}
        assert rows["mixed-bounds"]["bounds_changed"] == 16
        assert rows["mixed-bounds"]["point_candidates"] < rows["fresh"]["point_candidates"]
        assert all(c["primitive_mask"] == (1 << 1) | (1 << 10)
                   for c in rows["mixed-bounds"]["candidate_profile"])
    print("Lifecycle regression checks passed")


if __name__ == "__main__":
    main()
