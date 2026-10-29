<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO
SPDX-License-Identifier: MPL-2.0
-->

# Nearest-boundary recovery across shared cells

For the consolidated decision and current limitations, see the
[geometry query status](GEOMETRY_QUERY_STATUS.md). This report retains the
baseline and conclusions of its own experiment.

Follow-up on `experiment/geometry-recovery-search`, based on `bbf414d`.
The [preceding experiment](GEOMETRY_QUERY_EXPERIMENT.md) establishes the
navigation repair and per-cell indexed enumeration baseline.

## Profile and hypothesis

The current baseline's expensive ray was profiled with stage timers in a local
diagnostic build. Each of three recovery calls spends approximately 2.37 seconds
in filled-universe enumeration, performing about 140 million primitive tests.
Root enumeration takes approximately 0.4–1 ms and lattice processing about 4 ms
per call. These are attribution measurements from one selected ray, not a model
latency distribution. Diagnostic instrumentation and private parameters remain
in ignored local storage.

Prepared cell surface lists already remove repeated leaves within a cell.
Different cells still refer to many of the same cards. The experiment tests
whether solving a card once per universe occurrence reduces recovery time while
preserving first-hit behavior. Acceptance was set before comparison: at least
a twofold median reduction on the bounded heavy workload, matching result
digests, no persistent geometry storage, and investigation of repeatable control
regressions above 5%.

## Implementation and contracts

Nearest-only universe enumeration uses an optional bitset indexed by surface
index. All indexed cells in one invocation use the same local ray, while its
nearest-hit bound only shrinks. A card tested earlier therefore cannot supply
a newly eligible closer hit when another cell references it. First occurrence
order remains authoritative for ties. Distinct cards sharing a primitive keep
separate bits, preserving their identities.

The bitset belongs to one invocation and is freed before returning. Different
transforms and lattice elements get independent storage. Scratch is capped at
1 MiB per active universe scan; larger models and allocation failure use the
existing scan. Cells without prepared indexes retain their tree fallback.
There is no new geometry cache, public API, timer in production traversal, or
change to general full-breakpoint enumeration.

The independent regression creates many diagnostic cells sharing a plane, then
a second translated occurrence of the same universe whose plane is closer.
It checks analytic distance, surface identity and normal, tree-fallback parity,
and reduced solve counts. The deliberately overlapping diagnostic cells test
mathematical crossing enumeration, not physical ownership validity. Existing
regressions retain coincident aliases, nearby distinct planes, repeated and
transformed lattices, and the short-flight navigation fix.

## Validation

Native unit, integration and Lua suites and serial unit/integration suites pass.
Three alternating public-model process pairs retain all compared track, path,
surface, navigation and raster digests. The targeted changing-direction case
continues to complete 512 events without ambiguity or failure.

## Bounded comparison and decision

Baseline and candidate use the same release compiler settings
(`-O3 -march=native -DNDEBUG`), one worker, full paths and requested surface IDs.
Three process pairs alternate execution order. Each process performs one warm-up;
heavy processes have one measured query and controls have 30. Loading and
preparation are excluded from query timing. Processes run sequentially under a
90-second limit. All complete, and every compared output digest matches,
including endpoints, ownership, full hierarchy paths and requested boundary IDs.

| Targeted private workload | Baseline median (range) | Candidate median (range) |
| --- | --- | --- |
| Recovery-heavy query, 3 measured samples each | 7.265 s (7.222–7.284) | 1.830 s (1.817–1.854) |
| Recovery-free control, 90 measured samples each | 23.390 ms (22.789–24.042) | 23.190 ms (22.270–24.333) |

The heavy median improves by 3.97 times, meeting the predeclared twofold gate.
Control variation is approximately -0.9%; this is not a separate speedup claim.
The new regression fails against the immediate baseline's work count and passes
with the optimization. Retain the scoped change for review. These results apply
to the selected workloads; complete useful-resolution large-model slices and
interactive rendering remain open validation milestones. A 1.83-second query
still warrants profiling before choosing another optimization.

Private models, workload coordinates, comparison clients and raw logs remain
outside publishable files. Only aggregate observations are recorded here.

The [complete-slice validation](GEOMETRY_SLICE_VALIDATION.md) measures the next
end-to-end milestone: public frames and material-containing private detail views
complete at 512 and 1024 square, while a broader focused view remains bounded
by a tracing timeout.

The subsequent [prepared-universe recovery experiment](GEOMETRY_UNIVERSE_RECOVERY.md)
records the implementation and validation of this next step.
