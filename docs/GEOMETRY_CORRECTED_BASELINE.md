<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO
SPDX-License-Identifier: MPL-2.0
-->

# Corrected geometry baseline and remaining recovery work

Measured on 2026-10-29 at `2a3f8b2`, whose production code is unchanged from
the elliptical-torus correction at `0e9451e`. See the
[consolidated status](GEOMETRY_QUERY_STATUS.md) for the preceding decisions.
This investigation changes no production code and implements no optimization.

## Method

A separate release build uses Apple Clang 21.0.0 on macOS ARM64, with
`-O3 -march=native -DNDEBUG`, the native TinyPar backend, and
`ALEA_NUM_THREADS=1`. This establishes a one-worker CPU baseline; it does not
measure existing parallel scaling or GPU execution. The exact CPU brand was
unavailable through the sandbox; local metadata records the compiler and OS.

Complete slices use the existing benchmark, capped at one million segments
and 256 MiB of trace output. Each completing configuration has three sequential
processes, each producing a first frame and two warm frames. Times exclude
hashing and display/upload. Split mode measures tracing and rasterization
separately; fused mode reports their combined time. Preparation and loading are
separate from frame timing. Process peak RSS includes the whole process.

Private processes run sequentially. The broad 512 view has a 60-second external
whole-process timeout. Bounded row sets have 32 systematic samples at each
resolution, a 3-second per-row cooperative budget, a 90-second row-set budget,
and a 125-second external process timeout. These row sets are not random
samples and are not extrapolated into full-frame timings.

## Completed slices

All 144 public frames complete across four fixtures, two resolutions, two modes,
and three processes per configuration. Repeating, translated-fill, and
transformed-lattice fields match the preceding corrected build. Elliptical
fields agree across modes and with four separate pixel-by-pixel analytical
checks at both resolutions. Compared fields are cell ID, material ID, density,
and resolution flags. Public warm median frame times range from approximately
1.8 to 8.0 ms at 512 and 3.9 to 15.8 ms at 1,024. This is a baseline measurement,
not a speedup comparison.

All 36 private detail frames complete and match the prior completed-view
digests. Approximately 91.7% of pixels contain material, with no flagged pixels.

| Private detail configuration | First-frame median | Warm median (range, six frames) | Warm raster median |
| --- | ---: | ---: | ---: |
| 512 split | 2.910 s | 3.030 s (2.896–4.264 s) | 0.311 ms |
| 512 fused | 2.962 s | 2.901 s (2.886–3.653 s) | Included in frame |
| 1,024 split | 6.065 s | 5.849 s (5.788–5.925 s) | 1.172 ms |
| 1,024 fused | 5.941 s | 5.763 s (5.716–5.846 s) | Included in frame |

Across these configurations, median loading is 6.7–7.0 s and median preparation
is 5.3–5.5 s. One preparation run takes 9.3 s. Observed peak RSS varies roughly
from 1.8 to 2.9 GiB. Timing outliers and RSS variation preclude claims of a
small fused-mode advantage or a memory improvement. Tracing dominates the
measured split frames; rasterization accounts for less than 0.1% of their time.

The corrected broad 512 view still times out at 60 seconds without completing
its first frame. Loading takes 7.31 s and preparation 5.39 s in that screen.
There is no completed-field digest or full-frame timing for the interrupted
process. A broad 1,024 frame is not attempted after the 512 screen fails.
The one-second warm broad-slice target remains unmet; warming alone does not
make even the measured detail view meet it.

## Rows and responsiveness

Both uninstrumented row sets complete and match all 64 earlier row digests and
their segment/interval, snapshot, and ownership-query counters.

| Resolution | Total time for 32 rows | Median row | Slowest row |
| --- | ---: | ---: | ---: |
| 512 | 6.343 s | 255.9 ms | 704.8 ms |
| 1,024 | 7.443 s | 264.5 ms | 747.4 ms |

These are single baseline processes per row set, not repeated speedup evidence.
The cancellation check requests 100 ms per row within a one-second tracing
budget. It completes seven rows, cancels eight, and reports no query errors.
Cancelled calls return in 102.9–105.0 ms; the total tracing budget returns in
1.003 s. Completed queries after cancellations remain usable and cancelled
rows publish no digest. This meets the sampled 150 ms responsiveness target,
without establishing a universal hard deadline.

## Corrected recovery profile

Diagnostic copies of the current raycast source add stage timers and counters
only under ignored `build/`. A coarse build measures actual recovery calls,
indexed scans, fill traversal, and the selected walker's ownership-resolution
block. Ownership timing covers neighbor/path/root/full lookup in that block;
it does not include every later containment or interval-verification operation.
No timer is added per primitive. A second diagnostic counts repeated contexts
without skipping any solve or traversal. All eight diagnostic public control
frames and all 96 diagnostic private rows match release outputs.

| Coarse diagnostic, 32 rows | 512 | 1,024 |
| --- | ---: | ---: |
| Actual recovery calls | 21 | 27 |
| Row time | 6.284 s | 7.441 s |
| Recovery time | 5.480 s | 6.669 s |
| Indexed scan time, included in recovery | 4.553 s | 5.564 s |
| Ownership-resolution block | 0.785 s | 0.754 s |
| Indexed surface references | 233.0 million | 299.9 million |
| Fill-cell visits | 201.1 million | 258.9 million |

Recovery accounts for about 87–90% of diagnostic row time. Indexed scanning
accounts for about 83% of recovery time. Every observed universe scan uses its
prepared list; there are no fallback-cell scans or per-cell duplicate skips.
At 1,024, torus rejection removes 46.544 million of 46.592 million checked
torus calls (99.90%). Remaining indexed solves include 129.3 million planes and
85.0 million quadrics. Torus quartics are no longer the obvious first target.
These diagnostic times describe cost distribution, not release throughput.

## One next hypothesis: repeated identical recovery contexts

The shadow diagnostic resets a bounded 16,384-entry table for each actual
recovery call. Its exact key contains ordered-list identity and length, the
local ray including inverse directions, and the complete transform/inverse
values and flags. Hash collisions are resolved with full key comparison.
A repeat is counted only when its search window is contained in the previously
recorded window for that key. No work is skipped and no geometry is approximated.

In the 1,024 row set, 87,719 of 107,564 indexed calls repeat such a context:
81.6% of calls, accounting for 236.8 million of 299.9 million surface references
(79.0%). There are at most 735 distinct recorded contexts in an individual
recovery, and no table-capacity failures. This resolves the earlier uncertainty
about exact local-ray reuse for this workload. It does not establish how much
release time an implementation would save or prove a general cache contract.

The next bounded experiment should test skipping a repeated indexed scan
within one nearest-recovery call after an identical context has been scanned
successfully over a containing window. Preserve first-hit tie behavior, aliases,
normals, transforms, cancellation, and error propagation. Do not share a record
across recovery calls, rays, model edits, different transforms, or wider windows.
Do not prune recursive child traversal merely because one surface list repeats.

Use at most 1 MiB of optional scratch per recovery, with the existing scan on
allocation failure, table saturation, or uncertain inputs. Establish the rule
on independently constructed repeated-occurrence fixtures before private timing.
Retain only with three alternating release baseline/candidate process pairs,
matching completed results, a benefit above variability, no repeatable public
control regression above 5%, and explicit memory/cancellation accounting.
Recheck the broad-frame screen; a sampled gain alone does not meet its target.

No production cache has been implemented in this measurement step. A generic
GPU abstraction is still unsupported by these results; removing duplicated
queries is the narrower experiment now justified by evidence.

## Confidentiality and reproducibility

Private inputs, model paths, view coordinates, raw digests, extracted diagnostic
data, scripts, logs, and instrumented sources remain outside Git or under
ignored `build/corrected-baseline/`. Only aggregate results appear here. Public
fixtures remain independently specified. No models or diagnostics are published.
Existing production tests were not rerun for this measurement-only step;
the previous full-suite evidence remains scoped to the unchanged production
checkpoint, while the fresh output/oracle checks are described above.
