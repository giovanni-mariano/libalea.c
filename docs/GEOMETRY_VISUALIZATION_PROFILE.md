<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO
SPDX-License-Identifier: MPL-2.0
-->

# Broad-view baseline after ordered point queries

Measurement-only follow-up at `a8491e2` on 2026-10-29. This resumes the
visualization investigation after [recovery context reuse](GEOMETRY_RECOVERY_CONTEXT_REUSE.md)
and the separate [public point-query fix](GEOMETRY_POINT_QUERIES.md).
No production code changes in this step.

## Release measurements

The existing slice benchmark is rebuilt against the current release libraries
with Apple Clang 21, -O3 -march=native -DNDEBUG, on macOS ARM64 using native
TinyPar. One process runs at a time. Output budgets remain one million segments
and 256 MiB of trace storage. Hashing and display/upload are outside frame time.

All 32 public control frames match prior cell/material/density/flag digests:
four synthetic models at 512 and 1,024, split and fused, first and warm frame.
Both private 32-row sets also complete and retain their earlier output digests.
Their single-process total times are 2.142 s at 512 and 2.697 s at 1,024;
these are compatibility/baseline checks, not paired speedup evidence.

The same material-containing broad private 512 view has these results:

| Configured workers and mode | First frame | Warm frame | Load | Preparation |
| --- | ---: | ---: | ---: | ---: |
| 1, split | 36.046 s | 35.576 s | 8.468 s | 7.904 s |
| 1, fused | 35.144 s | Not measured | 6.724 s | 5.355 s |
| 4, split | 12.529 s | 11.761 s | 7.070 s | 5.375 s |

Each row represents one process, not a repeated timing distribution. The split
warm-series processes have a 120-second external timeout; the fused first-frame
screen retains the 60-second limit. Warm-series completion does not redefine
the earlier first-frame screen. Loading and preparation stay outside frame time.
One-worker split first-frame setup plus frame work totals about 52.4 s; the
four-worker equivalent is about 25.0 s. These exclude application display.

All five broad frames match each other and the preceding recovery build's
completed frame in all digest fields, with 88,597 material pixels and zero
resolution flags. Four configured workers give roughly a threefold warm-frame
ratio in this sample. This is useful existing CPU parallelism, not a repeatably
established scaling curve. Process peak RSS ranges roughly from 2.0 to 2.8 GiB;
these few runs establish no memory improvement.

The earlier recovery build's broad frame took about 45 s. These new runs were
not paired with that build, so their lower times cannot be attributed to the
public point-query change. The one-second broad warm-frame objective remains
unmet with either worker setting. No 1,024 broad frame, 3D renderer, device
backend, or application display latency was measured here.

## Full-frame diagnostic profile

An ignored copy of the current raycast source adds monotonic timers and
thread-local counters. One-worker diagnostics cover a complete broad split
frame, rather than extrapolating sampled rows. Four synthetic diagnostic frames
and the complete private diagnostic frame match release output. Instrumented
tracing takes 40.223 s; its times describe distribution, not release throughput.

| Diagnostic work | Time or count |
| --- | ---: |
| Actual nearest recovery calls | 349 |
| Recovery time | 36.946 s |
| Indexed scan time, inside recovery | 21.498 s |
| Selected-walker ownership-resolution block | 2.987 s |
| Filled-universe stage, including indexed scans | 35.042 s |
| Lattice stage | 1.726 s |
| Root stage | 0.173 s |
| Fill traversal cell-entry visits | 3,343,768,470 |
| Ordinary positive-fill candidates | 2,635,648 |
| Actual indexed scans | 256,516 |
| Surface references in actual indexed scans | 815,376,379 |
| Prepared universe-list calls | 1,389,891 |
| Reused cached contexts | 1,133,375 of 1,384,656 cache lookups |

Recovery accounts for about 92% of diagnostic tracing, and actual indexed scans
account for about 58% of recovery. All observed universe scans use prepared
lists; there are no fallback-cell scans or bitset duplicate skips. Actual solves
include 494.7 million planes and 226.4 million general quadrics; normal evaluation
is infrequent (5,457 indexed normal calculations). Torus quartics are no longer
the dominant arithmetic volume in this workload.

Ownership timing covers the selected walker's neighbor/path/root/full lookup
block; it excludes later containment verification and crossing work. Cached-helper
time includes its underlying scans. Stage and scan timers are nested and must
not be added together. The recovery time outside indexed scans includes fill
iteration, transforms, cache handling, and other work; it is not all removable
loop overhead. Positive-fill candidate counts exclude lattice-only cells, so
they must not be used to claim an exact percentage of terminal-cell visits.

## Next bounded experiment

Prepare an optional per-universe list of cells that can lead to children:
`fill_universe > 0` or an active lattice (`lat_type != 0 && lat_fill`). Preserve
original cell order. Nearest recovery can then iterate that subset in its
filled-universe recursion instead of walking every terminal cell. The existing
loop already skips ordinary non-fill/non-lattice cells without contributing
hits. This is a structural filtering hypothesis justified by the visit volume;
its speedup has not been measured or implemented here.

Retain all existing parent-box tests, transforms, child/lattice traversal,
primitive enumeration, first-hit ties, aliases, interruption, and error handling.
Do not substitute region bounds for complete primitive bounds or introduce new
approximate rejection. Scope initial consumption to nearest recovery. Missing,
incomplete, invalidated, oversized, or failed optional storage keeps the original
loop. Use the existing preparation/invalidation ownership conventions and cap
additional retained indices at 8 MiB per model, with no per-query list allocation.

Validate independently constructed mixed terminal/fill/lattice universes,
repeated and transformed children, childless universes, undefined fills,
coincident aliases, edits, clones, invalidation, and fallback. Check native and
serial integration and exact completed public/private results. Compare three
alternating uninstrumented process pairs; account for preparation, retained
memory, cancellation, and controls (no repeatable regression above 5%). Recheck
complete first and warm broad frames at one and four workers. The current
threefold worker ratio also needs repeated measurements before a scaling claim.

This experiment can remove avoidable CPU traversal before considering another
representation or hardware backend. Even eliminating all measured work outside
indexed scans would leave substantial tracing cost, so it must not be presented
as sufficient to reach the one-second target.

## Privacy and validation scope

Private models, paths, view coordinates, raw digests, diagnostic sources, logs,
and binaries remain outside Git or under ignored build/visualization-profile/.
Only aggregate private observations appear here. No new full test-suite run is
needed for this measurement-only change; the production suites at a8491e2
already passed, and fresh output compatibility checks are described above.
