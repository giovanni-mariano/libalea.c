#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0
"""Exercise replay gates, cache hits, edited-model isolation, and input rejection."""
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from compare_navigator_replays import compare_crossings, compare_runs
from generate_lifecycle_fixture import generate, generate_quadric


def check_adaptive_lattice(binary, model, workload):
    fixture = Path(__file__).resolve().parents[1] / "tests/data/mcnp_nested_lattice_transformed.mcnp"
    # Make the lattice the last candidate in a frequently scanned child.
    # Widen the world so the transformed samples also fall inside the local index.
    cards = "".join(f"{10000 + i} 0 -999 u=20 imp:n=1\n" for i in range(32))
    text = fixture.read_text().replace("100 0  1 -2", cards + "100 0  1 -2")
    model.write_text(text.replace("1  px -1", "999 s 200 0 0 1\n1  px -1").replace("21 px 3", "21 px -4"))
    samples = [(5, 0), (4, .25), (6.1, .25), (8.5, 0), (5.5, .5), (20, 0)]
    workload.write_text("".join(f"p {x} {y} 0\n" for x, y in samples) * 2000)
    command = [binary, "fresh", model, workload, "-"]
    reference = json.loads(subprocess.check_output(command + ["--hierarchy-points"]))
    for regions in ([], ["--child-regions", "10"]):
        result = json.loads(subprocess.check_output(command + ["--root-points", *regions,
                                                               "--adaptive-child-bounds", "--repeat-points"]))
        assert result["point_digest"] == reference["point_digest"]
        assert result["repeat_point_s"] > 0
        child = next(row for row in result["root_scans_by_universe"] if row["universe"] == 20)
        assert child["bounded_cells"] == 33
        first_bytes = result["selected_points"]["bytes"]
        assert result["repeat_point_index_bytes"] >= first_bytes > result["point_index_bytes_at_ready"]
        if regions:
            assert result["child_universes"] == [10]
            assert result["child_index"] == "regions-and-adaptive-bounds"
            for row in result["root_scans_by_universe"]:
                if row["universe"] == 10:
                    assert "bounded_cells" not in row

    workload.write_text(workload.read_text() * 6)
    control = json.loads(subprocess.check_output(command + ["--root-points", "--study-promotion-ms", "0"]))
    env = dict(os.environ, ALEA_PROFILE_CHILD_POINTS="1")
    for strategy in ([], ["--study-child-only"]):
        run = json.loads(subprocess.check_output(command + ["--root-points", "--study-promotion-ms", "1",
                                                            *strategy], env=env))
        assert run["point_digest"] == control["point_digest"]
        assert run["promoted_universe"] == 20 and 4096 <= run["after_points"] < run["points"]
        assert run["rebuild_s"] > 0
    args = command + ["--root-points", "--study-promotion-ms", "1", "--study-background-child"]
    single = subprocess.run(args, capture_output=True, env=dict(env, ALEA_NUM_THREADS="1"))
    assert single.returncode == 2 and b"needs two workers" in single.stderr
    background = subprocess.run(args, capture_output=True, env=dict(env, ALEA_NUM_THREADS="2"))
    if background.returncode == 2 and b"needs two workers" in background.stderr:
        return  # Serial-only build.
    background.check_returncode()
    result = json.loads(background.stdout)
    assert result["point_digest"] == control["point_digest"]
    assert result["promoted_universe"] == 20 and result["next_result_s"] > 0
    assert len(result["windows_s"]) == len(control["windows_s"])


def check_navigator_replay(binary, workload):
    fixture = Path(__file__).resolve().parents[1] / "tests/data/simple_sphere.mcnp"
    workload.write_text("r -60 0 0 2 0 0 120\nr 0 0 0 2 0 0 10\n")
    for mode in ("fast", "strict", "interval"):
        output = subprocess.check_output([binary, "navigator", fixture, workload, mode])
        run = json.loads(output)
        assert run["failures"] == 0 and run["mode"] == mode
        assert run["first_event_s"] >= 0 and "peak_rss_bytes" in run
        outside, inside = run["flights"]
        assert outside["stop"] == "vacuum" and inside["stop"] == "distance-limit"
        exit_event, = outside["events"]
        assert exit_event[:2] == [3, 10] and exit_event[6][0] == 6
        limit, = inside["events"]
        assert limit[:2] == [1, 10] and limit[5] == limit[6] == inside["source"]

    fixtures = Path(__file__).resolve().parents[1] / "benchmarks/workloads/transport"
    # Analytic expected event kind, surface, distance and world position.
    # Includes both cloned histories, periodic destinations and one-ULP ties.
    expected = [
        (0, -1, .5, [.5, 0, 0]), (1, -1, 1, [-.5, 0, 0]),
        (1, -1, .5, [1, 0, 0]), (4, 2, 1, [2, 0, 0]),
        (0, -1, .25, [-1.75, 0, 0]), (4, 1, 1.5, [-2, 0, 0]),
        (4, 3, 2, [0, -2, 0]), (4, 4, 4, [0, 2, 0]),
        (3, 6, 2.5, [0, .5, 2]), (0, -1, 2-2**-52, [2-2**-52, 0, 0]),
        (4, 2, 2, [2, 0, 0]),
    ]
    reference = None
    for mode in ("strict", "fast", "interval"):
        run = json.loads(subprocess.check_output([
            binary, "navigator-script", fixtures / "box.i", fixtures / "box.ops", mode]))
        events = [op["event"] for op in run["operations"] if "event" in op]
        assert [(e[0], e[2], e[1], e[7]) for e in events] == expected
        boundaries = {-1: (0, [0, 0, 0]), 1: (3, [1, 0, 0]), 2: (3, [1, 0, 0]),
                      3: (1, [0, 1, 0]), 4: (2, [0, 1, 0]), 6: (4, [0, 0, 1])}
        source = run["operations"][0]["location"]
        assert source[:4] == [0, 1, 1, 1] and source[4] != 0
        for event in events:
            assert (event[3], event[4]) == boundaries[event[2]]
            assert event[5] == source
            if event[0] == 3:
                assert event[6][0] == 6  # vacuum leaks the particle
            else:
                assert event[6] == source
        assert [op["position"] for op in run["operations"] if op["op"] == "periodic"] == [
            [-2, 0, 0], [2, 0, 0]]
        if reference is None:
            reference = run
        else:
            compare_runs(reference, run)

    for mode in ("strict", "fast", "interval"):
        regional = json.loads(subprocess.check_output([
            binary, "navigator-script", fixtures / "box.i", fixtures / "box.ops", mode,
            "regional", "-12", "12", "-12", "12", "-12", "12"]))
        assert regional["producer"] == "regional"
        compare_runs(reference, regional)

    # The comparison must catch an incorrect physical surface even when the
    # distance, endpoint, material and subsequent trajectory still agree.
    run["operations"][6]["event"][2] = 99
    try:
        compare_runs(reference, run)
    except ValueError as error:
        assert "operations[6] advance parent" in str(error) and "surface:" in str(error)
    else:
        raise AssertionError("surface mismatch was accepted")

    lattice = fixtures.parents[2] / "tests/data/mcnp_lattice_repeating.mcnp"
    workload.write_text("r 0 .5 0 1 0 0 8\nr 4 .5 0 -1 0 0 8\n")
    for mode in ("fast", "strict", "interval"):
        command = [binary, "navigator", lattice, workload, mode]
        baseline = json.loads(subprocess.check_output(command))
        spans = json.loads(subprocess.check_output(command + [
            "regional", "-12", "12", "-12", "12", "-12", "12"]))
        assert spans["regional_spans"] > 0 and spans["failures"] == 0
        assert compare_crossings(baseline, spans)["matched_events"] == 10
        for flight in spans["flights"]:
            first = flight["events"][0]
            assert first[2] == 0 and first[5][1] == first[6][1] and first[5][4] != first[6][4]

    # An interior action card still reflects a particle even though the
    # selected cell owns both sides of the plane.
    internal = workload.with_name("internal-action.i")
    internal.write_text("Internal reflecting plane\n1 1 -1 -2 (-1 : 1) imp:n=1\n"
                        "2 0 2 imp:n=0\n\n*1 px 0\n2 so 10\n\nm1 1001 1\n")
    workload.write_text("restart parent -1 0 0 1 0 0\nadvance parent inf 5\n"
                        "reflect parent\nadvance parent inf 2\n")
    command = [binary, "navigator-script", internal, workload, "strict"]
    baseline = json.loads(subprocess.check_output(command))
    reflected = json.loads(subprocess.check_output(command + ["regional"]))
    compare_runs(baseline, reflected)
    assert reflected["operations"][1]["event"][:4] == [4, 1, 1, 1]
    assert reflected["operations"][-1]["event"][7] == [-2, 0, 0]

    check_entering_surface(binary, workload)

    # Two physical actions at the same corner cannot be disambiguated.
    workload.write_text("r 0 0 0 1 1 0 5\n")
    failed = subprocess.run([binary, "navigator", fixtures / "box.i", workload, "strict"],
                            capture_output=True)
    assert failed.returncode == 1
    trace = json.loads(failed.stdout)
    flight, = trace["flights"]
    assert trace["first_event_s"] is None and trace["failures"] == 1
    assert flight["failure_distance"] == 0 and flight["stop"] == "advance-error"
    assert "ambiguous" in flight["diagnostic"]

    # An invalid physical operation records its input and returns nonzero.
    workload.write_text("restart parent 0 0 0 1 0 0\nreflect parent\n")
    failed = subprocess.run([binary, "navigator-script", fixtures / "box.i", workload,
                             "strict"], capture_output=True)
    assert failed.returncode == 1
    last = json.loads(failed.stdout)["operations"][-1]
    assert last["op"] == "reflect" and last["rc"] != 0


def check_entering_surface(binary, workload):
    # The earlier cell starts at x=0. The departing cell's distinct, tilted
    # plane intersects just past x=0, inside the surface-group tolerance.
    model = workload.with_name("entering-surface.i")
    deck = ("Entering face attribution\n1 1 -1 1 -3\n2 2 -1 -2 -3\n3 0 3\n\n"
            "1 px 0\n2 p 1 .25 0 5e-11\n3 so 10\n\nm1 1001 1\nm2 1001 1\n")
    start = "restart parent -1 0 0 1 0 0\nadvance parent inf 3\n"
    command = [binary, "navigator-script", model, workload, "strict", "regional",
               "-12", "12", "-12", "12", "-12", "12"]
    model.write_text(deck)
    workload.write_text(start + "advance parent inf 2\n")
    ops = json.loads(subprocess.check_output(command))["operations"]
    event = ops[1]["event"]
    assert event[:5] == [2, 1, 1, 0, [1, 0, 0]]
    assert (event[5][1], event[6][1], event[7]) == (2, 1, [0, 0, 0])
    assert ops[-1]["event"][7] == [2, 0, 0]

    # A mirror on either face makes this two-card group ambiguous.
    for surface in ("1 px", "2 p"):
        model.write_text(deck.replace(surface, "*" + surface))
        failed = subprocess.run(command, capture_output=True)
        assert failed.returncode == 1
        op = json.loads(failed.stdout)["operations"][-1]
        assert op["op"] == "advance" and "ambiguous" in op["diagnostic"]

    # Only the entering face lies within tolerance. Its local X normal
    # becomes world Y in this rotated fill and reflection keeps the owner.
    model.write_text("Entering mirror in a rotated fill\n"
                     "4 0 -3 fill=1 (4 0 0 0 1 0 -1 0 0 0 0 1)\n"
                     "1 1 -1 1 u=1\n2 2 -1 -2 u=1\n3 0 3\n\n"
                     "*1 px 0\n2 p 1 .25 0 1e-9\n3 so 10\n\nm1 1001 1\nm2 1001 1\n")
    workload.write_text("restart parent 4 -1 0 0 1 0\nadvance parent inf 3\n"
                        "reflect parent\nadvance parent inf 2\n")
    ops = json.loads(subprocess.check_output(command))["operations"]
    event = ops[1]["event"]
    assert event[:5] == [4, 1, 1, 1, [0, 1, 0]] and event[5] == event[6]
    assert ops[-1]["event"][7] == [4, -2, 0]


def check_ray_root_candidates(binary, model, workload):
    fixtures = Path(__file__).resolve().parents[1] / "tests/data"
    # Leave part of the parent fill uncovered, so resolution flags and the
    # child's first crossing within a gap participate in the segment digest.
    model.write_text((fixtures / "simple_fill.mcnp").read_text().replace("3 0 20 U=1\n", ""))
    workload.write_text("r 0 0.25 0 1 0 0 12\nr 10 0.25 0 -1 0 0 12\n"
                        "r 5 0.25 0 1 0 0 12\nr 0 0 0 1 0 0 12\n")
    for fixture in (model, fixtures / "mcnp_nested_lattice_transformed.mcnp",
                    fixtures / "mcnp_lattice_repeating.mcnp", fixtures / "mcnp_deep_translated_fills.mcnp"):
        command = [binary, "fresh", fixture, workload, "-"]
        reference = json.loads(subprocess.check_output(command + ["--hierarchical-rays"]))
        regional = json.loads(subprocess.check_output(command + ["--study-ray-root"]))
        assert regional["regional_root_hits"] > 0, fixture
        assert regional["segments"] == reference["segments"] > 0, fixture
        assert regional["ray_digest"] == reference["ray_digest"], fixture
        paired = json.loads(subprocess.check_output(command + ["--compare-ray-root"]))
        assert paired["rays_per_method"] == reference["rays"], fixture
        assert paired["segments_per_method"] == reference["segments"], fixture


def check_public_query_study(binary, workload):
    fixture = Path(__file__).resolve().parents[1] / "tests/data/mcnp_lattice_repeating.mcnp"
    command = [binary, "public", fixture, workload]
    workload.write_text("p 0 0 0\np 0.5 0 0\n" * 16)
    reference = None
    for mode in ("points", "lookup", "adaptive", "auto-points"):
        record = json.loads(subprocess.check_output(command + [mode]))
        assert [row["n"] for row in record["batches"]] == [1, 16, 32]
        digests = [row["digest"] for row in record["batches"]]
        if reference is None:
            reference = digests
        assert digests == reference
        if mode == "auto-points":
            assert record["promoted_after"] == 0 and record["policy_build_s"] == 0
    workload.write_text("p 0 0 0\np 0.5 0 0\n" * 128)
    reference = json.loads(subprocess.check_output(command + ["lookup"]))
    promoted = json.loads(subprocess.check_output(command + ["auto-points", "0"]))
    assert promoted["promoted_after"] == 64
    assert [r["digest"] for r in promoted["batches"]] == [r["digest"] for r in reference["batches"]]
    workload.write_text("r -2 0 0 1 0 0 4\nr 2 0 0 -1 0 0 4\n")
    hierarchy = json.loads(subprocess.check_output(command + ["rays"]))["rays"]
    for options in (["ready-rays"], ["prepared", "-4", "4", "-4", "4", "-4", "4"]):
        candidate = json.loads(subprocess.check_output(command + options))["rays"]
        assert len(hierarchy) == len(candidate) == 2
        for expected, actual in zip(hierarchy, candidate):
            assert len(expected) == len(actual) > 0
            for a, b in zip(expected, actual):
                assert a[2:] == b[2:] and max(abs(a[i] - b[i]) for i in (0, 1)) < 1e-12


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
        check_adaptive_lattice(binary, model, workload)
        check_navigator_replay(binary, workload)
        check_ray_root_candidates(binary, model, workload)
        check_public_query_study(binary, workload)
    print("Lifecycle regression checks passed")


if __name__ == "__main__":
    main()
