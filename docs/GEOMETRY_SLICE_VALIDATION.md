<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO
SPDX-License-Identifier: MPL-2.0
-->

# Complete CPU slice validation

This pass tests complete square rasters after the occurrence-scoped recovery
optimization at `8c8789f`, compared with its `bbf414d` baseline. See the
[recovery findings](GEOMETRY_RECOVERY_SEARCH.md) for the isolated query results.

## Workload and measurement

The [slice benchmark](../benchmarks/geometry_slice_queries.c) loads one model,
prepares acceleration once, and produces cell/material IDs, density and resolution
flags at 512 by 512 and 1024 by 1024. Split mode times compact row tracing and
rasterization separately; fused mode exercises the production combined API.
Timing excludes digest calculation, field statistics, image encoding and GUI
display. Split total is the sum of trace and raster phases, excluding its progress
notification. The first frame includes cold query-buffer/page costs; the second
uses the prepared system and reused buffers.

Both builds use identical release settings (`-O3 -march=native -DNDEBUG`) and
one worker. Public fixtures use three process pairs in alternating build order,
with two frames in each mode. Every expected full-frame result is required;
cell/material/density/flag digests must match across both builds, modes and
repeats. Agreement is regression evidence. A separate analytic pixel check on
the public deep-fill fixture independently verifies the x=2 plane, positive mass
density, material and zero flags at every pixel in both sizes and modes.

Caller pixel storage is approximately 4.25 MiB at 512 square and 17 MiB at 1024.
The benchmark caps dimensions at 2048, trace output at 256 MiB and segments at
one million. These limits bound caller buffers and published track output;
source geometry, preparation and internal traversal storage are additional.
Process peak RSS is recorded, not constrained by a whole-process memory cap.

The transformed-lattice viewport includes undefined-fill pixels. Their flags
are compared explicitly; matching images do not establish unique physical
ownership throughout that viewport.

## Reproduction

Build and run from the repository root:

```sh
make RELEASE=1 BUILD_DIR=build/slice-release BIN_DIR=build/slice-release-bin \
    build/slice-release-bin/geometry_slice_queries
ALEA_NUM_THREADS=1 build/slice-release-bin/geometry_slice_queries \
    tests/data/mcnp_deep_translated_fills.mcnp 1024 2 -4 6 -1 1 0 split
```

Use `fused` instead of `split` to measure the combined API. Link the same benchmark
source against a separately built baseline archive to compare the two revisions.
Private models and view coordinates must be supplied locally and all derived
output kept in ignored storage. This benchmark writes field digests and timing,
not images or geometry exports.

## Completed public frames

Candidate medians of the three second (warm) frames, in milliseconds:

| Public fixture | 512 square split total | 1024 square split total | 1024 square raster phase | 1024 square fused total |
| --- | --- | --- | --- | --- |
| Repeating lattice | 6.456 | 13.600 | 1.096 | 13.654 |
| Eight translated fills | 1.900 | 4.015 | 0.643 | 3.898 |
| Transformed nested lattice | 7.456 | 15.326 | 1.034 | 14.975 |

All fields match the baseline and fused/split modes at both resolutions. Public
peak process RSS is approximately 7–21 MiB. These workloads establish completed
raster behavior and regression compatibility; they do not predict large-model
latency or establish a new optimization speedup.

## Private completion screens

Private processes ran sequentially, one worker, with one cold frame per run.
Candidate processes were limited to 60 seconds and baseline processes to 30,
including loading and preparation. These asymmetric limits screen completion;
they are not a controlled speedup comparison. No image files were written.

A focused large-model view timed out at both 512 square and 1024 square in both
builds. Logs reach the start of tracing but never its completion, so no raster
phase or full-field comparison exists for those four runs. Their last RSS record
is from preparation, not a final query memory measurement.

A narrow initial detail view completes with matching fields in approximately
1.24 seconds at 512 square and 2.48 seconds at 1024 square, but contains zero
material pixels. It is retained as a void-region diagnostic and does not satisfy
useful geometry-view validation.

An additional material-containing detail view was selected from a material
interval in the previously used center track, with padding to include its
interfaces. Its coordinates and selection diagnostics remain local. It completes
in both builds and both sizes, with matching all-field digests:

| Material-containing detail | Baseline split total | Candidate trace | Candidate raster | Candidate split total |
| --- | --- | --- | --- | --- |
| 512 square | 2.720 s | 3.184 s | 0.908 ms | 3.184 s |
| 1024 square | 5.618 s | 5.582 s | 3.568 ms | 5.586 s |

Approximately 91.7% of pixels contain material and no pixel has an undefined-fill
flag. Completed private process peaks are approximately 2.3–2.8 GiB across these
runs. Loading and preparation add roughly 11–12 seconds per fresh process,
separately from query timings. Single cold runs and process/allocator variation
prevent attributing timing or RSS differences to the recovery optimization.
Baseline agreement is not an independent truth oracle for private geometry.

## Decision

Useful-resolution material-containing local views now have completed field
comparisons. The broader focused view still fails its declared deadline.
Rasterization is negligible relative to tracing: increasing pixel width alone
is cheap compared with adding traced rows. Image encoding, UI latency and 3D
rendering are outside this measurement.

The next performance work should measure per-row tracing costs and their
distribution in the failing view, with bounded work and progress/cancellation.
This pass establishes the stage containing the bottleneck; it does not determine
how much is concentrated in a few slow rows versus distributed across the view.
It provides no evidence for changing raster storage or adding another renderer
backend. Private inputs, coordinates, raw digests and diagnostic output remain
local and ignored; public fixtures are independently authored.

## Bounded row profiling

The benchmark also accepts `rows SAMPLES ROW_MS TOTAL_MS` in place of its final
mode, with frame count 1. For example, the public analytic view can be profiled
with:

```sh
ALEA_NUM_THREADS=1 build/slice-release-bin/geometry_slice_queries \
    tests/data/mcnp_deep_translated_fills.mcnp 1024 1 -4 6 -1 1 0 \
    rows 32 3000 90000
```

It selects one actual row center from each equal index stratum, preserving the
production row-origin arithmetic, direction, span and requested fields. A
single prepared model is retained. Completed rows report their segment-output
digest and existing internal interval/snapshot/ownership work counters. It
prints row-start and row-end progress. The profiler uses the public signal-safe
interrupt operation to request cancellation via its standalone SIGALRM handler.
An alarm is disarmed before clearing the flag, so it cannot poison the following
row. This is cooperative cancellation; a parent process must additionally impose
a hard timeout. Exit status 3 means cancellation or an incomplete sample set,
1 means error, and 2 means invalid arguments. Partial rows have no output digest.

Public analytic sampled rows retain the independently specified plane crossing
at distance 6 and final endpoint 10 in both resolutions. A mixed diagnostic run
cancels its first row and completes the next seven, confirming successful reuse
after cancellation. Existing split/fused complete-frame modes remain available.

Private sampling uses 32 rows per resolution, a 3-second per-row cancellation
request, a 90-second tracing budget and a 125-second hard parent-process timeout.
Preparation is outside the tracing budget. Both current-checkpoint sample sets
complete without cancellations or errors:

| Current checkpoint | Median row | 90th percentile row | Maximum row | Rows over 0.5 s | Total sampled query time |
| --- | --- | --- | --- | --- | --- |
| 512-row grid, 32 sampled | 0.622 s | 1.245 s | 1.277 s | 17 / 32 | 13.620 s |
| 1024-row grid, 32 sampled | 0.673 s | 1.279 s | 1.909 s | 23 / 32 | 20.451 s |

Percentiles use the nearest-rank convention. These systematic samples use
slightly different physical row centers at different resolutions and are not
random samples with a statistical confidence bound. One-row API calls also
include their own setup/publication overhead, so sums are not an exact timing
decomposition or prediction of the full-frame batch.

Rows without a recovery snapshot have median latency around 23 ms. Completed
rows with one snapshot have medians around 0.63–0.71 seconds, and those with two
around 1.25–1.30 seconds. Snapshot counts are an existing walker diagnostic,
not a direct measurement of time inside recovery. They identify the association
for subsequent profiling without introducing production timers. The four
slowest sampled rows account for about 37% and 29% of sampled query time,
respectively. Many sampled rows are expensive; a few extreme outliers alone do
not explain the sampled cost.

The same row sets were checked against `bbf414d` with identical limits. The
baseline completes 28/32 and 29/32 rows; its remaining seven rows are cancelled.
All 57 mutually completed rows match in output digest, segment/interval counts,
snapshot counts, and ownership query counters. Cancelled baseline rows are
censored and have no output-comparison claim. Both baseline runs finish their
sample attempts within the total tracing budget and hard process limit.

Cancellation remains coarse inside exceptional recovery. Baseline rows asked
to stop at three seconds return after approximately 4.96–6.08 seconds. A separate
current-checkpoint test requests cancellation at 100 ms and observes returns
at approximately 662–665 ms, with two cancelled rows, no query errors, and an
explicitly incomplete sample set. No partial row digest is published. This
explains why the parent process limit is necessary; the cooperative deadline
is a request, not a hard per-row bound. Alarm delivery is confined to the profiling
thread by blocking SIGALRM before preparation creates workers and unmasking it
only after the profiler installs its handler.

The next experiment should separate repeated reference scanning from remaining
primitive solves inside current recovery. A prepared first-occurrence-ordered
universe surface list is a candidate for avoiding the repeated per-cell scans,
subject to a measured memory/preparation budget and the same alias/occurrence
checks. Recovery also needs finer interruption points before interactive
progress can rely on short cancellation latency. Load balancing alone does not
remove the substantial recovery work observed across the sampled rows.

The subsequent [prepared-universe recovery experiment](GEOMETRY_UNIVERSE_RECOVERY.md)
records the implementation and validation of this next step.
