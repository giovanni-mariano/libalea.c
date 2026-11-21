#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0
"""Report the first difference between two navigator traces, ignoring timings."""
import argparse
import json
import math
from itertools import zip_longest
from pathlib import Path


FIELDS = ("kind", "distance", "surface", "boundary type", "normal", "before", "after", "position")


def same(context, field, a, b):
    if a != b:
        raise ValueError(f"{context}: {field}: expected {a!r}, got {b!r}")


def compare_event(context, a, b):
    for i, (left, right) in enumerate(zip_longest(a, b)):
        same(context, FIELDS[i], left, right)


def compare_flight(context, a, b):
    if a == b:
        return
    same(context, "input ray", a.get("ray"), b.get("ray"))
    same(context, "source", a["source"], b["source"])
    for i, (left, right) in enumerate(zip_longest(a["events"], b["events"])):
        where = f"{context}, event {i}, ray {a.get('ray')}"
        if left is None or right is None:
            same(where, "event", left, right)
        compare_event(where, left, right)
    same(context, "stop", a["stop"], b["stop"])
    for field in ("failure_distance", "diagnostic"):
        same(context, field, a.get(field), b.get(field))


def compare_runs(reference, candidate):
    key = "flights" if "flights" in reference else "operations"
    if not reference.get(key) or key not in candidate:
        raise ValueError("missing or incompatible navigator traces")
    for i, (left, right) in enumerate(zip_longest(reference[key], candidate[key])):
        context = f"{key}[{i}]"
        if left is None or right is None:
            raise ValueError(f"{context}: trace length differs")
        if key == "flights":
            compare_flight(context, left, right)
        else:
            context += f" {left['op']} {left['particle']} {left['args']}"
            for field in ("op", "particle", "args", "rc", "location", "position", "diagnostic"):
                same(context, field, left.get(field), right.get(field))
            if "event" in left and "event" in right:
                compare_event(context, left["event"], right["event"])
            else:
                same(context, "event", left.get("event"), right.get("event"))
    same(key, "failures", reference.get("failures"), candidate.get("failures"))


def ownership_crossings(events):
    distance = 0.0
    selected = []
    for event in events:
        distance += event[1]
        # Preserve collisions, limits and all boundary actions. Only a
        # transmissive boundary with identical occurrence and state is skipped.
        if event[0] == 2 and event[3] == 0 and event[5] == event[6]:
            continue
        selected.append((distance, event))
    return selected


def compare_crossings(reference, candidate):
    left, right = reference["flights"], candidate["flights"]
    same("navigator", "validation mode", reference["mode"], candidate["mode"])
    same("navigator", "flight count", len(left), len(right))
    counts = [0, 0]
    total = 0
    maximum = [0.0, 0.0, 0.0]
    for i, (a, b) in enumerate(zip(left, right)):
        where = f"flight {i}"
        for field in ("ray", "source", "stop", "failure_distance", "diagnostic"):
            same(where, field, a.get(field), b.get(field))
        x, y = ownership_crossings(a["events"]), ownership_crossings(b["events"])
        same(where, "remaining event count", len(x), len(y))
        total += len(x)
        counts[0] += len(a["events"]) - len(x)
        counts[1] += len(b["events"]) - len(y)
        for j, ((ta, ea), (tb, eb)) in enumerate(zip(x, y)):
            where = f"flight {i} crossing {j}"
            for k in (0, 2, 3, 5, 6):
                same(where, FIELDS[k], ea[k], eb[k])
            for k, (label, u, v, tol) in enumerate((
                    ("world distance", [ta], [tb], 1e-9),
                    ("position", ea[7], eb[7], 1e-9),
                    ("normal", ea[4], eb[4], 1e-12))):
                if not all(math.isclose(x0, y0, rel_tol=0, abs_tol=tol)
                           for x0, y0 in zip(u, v, strict=True)):
                    raise ValueError(f"{where}: {label}: expected {u}, got {v}")
                maximum[k] = max(maximum[k], *(abs(x0-y0) for x0, y0 in zip(u, v)))
    return {"matched_events": total, "same_owner_transmissions": counts,
            "max_distance_delta": maximum[0], "max_position_delta": maximum[1],
            "max_normal_delta": maximum[2]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ownership-crossings", action="store_true",
                        help="compare cumulative flight endpoints and report omitted same-owner transmissions")
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    args = parser.parse_args()
    try:
        a, b = json.loads(args.reference.read_text()), json.loads(args.candidate.read_text())
        if args.ownership_crossings:
            print(json.dumps(compare_crossings(a, b)))
        else:
            compare_runs(a, b)
    except (ValueError, KeyError, OSError) as error:
        parser.exit(1, f"{error}\n")
