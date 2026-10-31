<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO
SPDX-License-Identifier: MPL-2.0
-->

# Remaining recovery work after prepared universe lists

For the consolidated decision and current limitations, see the
[geometry query status](GEOMETRY_QUERY_STATUS.md). This report retains the
baseline and conclusions of its own experiment.

Profile of `5a43f40`, following [prepared universe recovery](GEOMETRY_UNIVERSE_RECOVERY.md).
The production implementation is unchanged in this investigation.

## Method

An ignored instrumented copy of `raycast.c` was linked into a separate static
library and the existing row profiler. Instrumentation activates only during
actual selected-owner recovery. Timers surround root enumeration, filled
universe enumeration, lattice enumeration, and complete indexed-list calls;
there is no timer per primitive. Thread-local counters record reference visits,
primitive types, normals, occurrence scans, and fill-cell visits. Indexed-list
times are nested inside stage times and must not be added to them.

The established 32 systematic row samples at each of 512 and 1,024 resolution
run sequentially with one worker. All 64 completed rows match the existing
uninstrumented results in geometry digest, segment/interval counts, snapshot
counts, and ownership query counters. The profile observes 21 actual recovery
calls at 512 and 27 at 1,024. Timings include instrumentation and are for locating
cost, not a new speedup comparison or a full-frame prediction.

A second pass captures every 1,024th indexed solve per type into a bounded
64-entry reservoir per type per recovery. It replays these ray/payload pairs
through the existing intersection routine in batches of 128 repetitions, three
passes per batch. The reservoir selection is deterministic. These are warm
payload microbenchmarks; they exclude payload lookup/copy, memory locality of
the real model, traversal, normals, and publication. They establish relative
arithmetic cost, not each type's exact share of production time. Rare types
have insufficient samples for useful timing claims. Replay output is ignored
and cannot change the geometry results. All 32 replay-pass rows match.

A third pass computes existing primitive boxes for positive-width tori with
nonnegative major radius and positive minor radius, and tests them against the
same local ray and enumeration window. It still performs every intersection,
counts roots inside that window, and records any rejection/root conflict.
All 32 rows in this pass also match. This is a diagnostic screen, not an enabled
optimization or a proof of conservative floating-point rejection.

Private source data, paths, coordinates, hashes, captured payloads, instrumented
sources, binaries, and raw logs remain in ignored `build/recovery-remaining/`
or outside the repository. This report contains aggregate observations only.

## Findings

| Measurement | Resolution 512, 32 rows | Resolution 1,024, 32 rows |
| --- | ---: | ---: |
| Total sampled row time | 8.602 s | 10.565 s |
| Actual recovery time | 7.821 s | 9.791 s |
| Indexed-list time, included above | 6.984 s | 8.752 s |
| Indexed primitive solves | 233.0 million | 299.9 million |
| Indexed-list calls | 83,602 | 107,564 |
| Prepared-list calls | 83,602 | 107,564 |
| Fill-cell visits | 201.1 million | 258.9 million |
| Normals calculated in indexed scans | 426 | 479 |

Recovery accounts for approximately 91–93% of sampled row time. Approximately
89% of recovery time is inside indexed-list scanning. Filled-universe
enumeration accounts for about 99% of recovery stage time, including its indexed
scans; root enumeration is about 0.1% and lattice enumeration about 1.2%.
All observed universe surface scans use prepared lists, with no fallback-cell
scans and no bitset duplicate skips. The remaining solve volume comes from
many occurrences, even after repeated card references within a universe have
been removed. This does not establish whether some occurrences share exactly
the same local ray.

In the 1,024 sample set, planes account for 129.3 million indexed solves,
general quadrics for 85.0 million, and tori for 46.6 million. Median replay costs
are approximately 4 ns for planes, 9 ns for quadrics, and 81 ns for the dominant
torus orientation. Tori are approximately 15.5% of indexed solves but have much
higher arithmetic cost. Normal calculation is infrequent and is not the first
target suggested by this profile.

The torus-box diagnostic would reject 46,557,971 of 46,591,802 eligible solves:
99.93%. There are zero observed conflicts with roots in the enumeration window.
In fact, none of these eligible torus solves returns a root in that window in
this sampled view. This is a useful opportunity measurement with a narrow
scope: it provides no independent validation of torus intersections or bound
safety when crossings, tangencies, or difficult coordinates are present.

## Next experiment

Test a conservative finite-torus rejection bound before the quartic solver in
nearest recovery. Reuse the current local ray and existing occurrence transform;
keep the existing solver for every uncertain, degenerate, invalid, unbounded,
or non-finite case. Bounds must enclose the complete primitive surface, including
roots outside the owning cell's region. Cell-region boxes are insufficient for
this enumeration contract. Preserve card order, aliases, distances, and normals.

First validate independent crossing and miss cases on all torus axes, elliptic
tori, tangencies and nearby rays, parallel and tiny nonzero directions, short
windows near boundaries, large coordinates, repeated transformed occurrences,
and degenerate parameters. The rejection must be conservative under rounding;
a raw existing box/slab test should not be enabled merely because this private
sample has no observed conflicts. Ambiguous tests should execute the solver.
Then compare exact completed outputs, cancellation, preparation/memory cost,
and repeated heavy/control timings. Retain only a repeatable benefit without
control regression above the existing 5% screening threshold.

This bounded experiment needs no new persistent scene hierarchy. The remaining
fill-cell traversal cost and possible exact local-ray reuse can be investigated
after measuring how much torus rejection actually saves.

The subsequent [torus rejection implementation](GEOMETRY_TORUS_RECOVERY.md)
records its conservative envelope, compatibility limits, tests, and measurements.
