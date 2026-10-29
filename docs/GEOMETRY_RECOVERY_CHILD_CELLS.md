<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO
SPDX-License-Identifier: MPL-2.0
-->

# Ordered child cells for nearest recovery

Follow-up to the [complete-frame profile](GEOMETRY_VISUALIZATION_PROFILE.md),
compared with `a4d07b8` (production code at `a8491e2`). That profile counted
3.34 billion cell-entry visits inside fill recursion. Most loop entries did
not contribute a fill or lattice traversal.

## Implementation and contract

Preparation builds an optional ordered list per universe containing positive
fills and active lattices. Nearest recovery iterates this list; ordinary
all-hit enumeration continues to use the original universe entries. Parent
bounding-box checks, transforms, recursion, surface aliases, nearest-hit ties,
and failure handling keep their existing order and behavior. No new geometric
rejection or intersection formula is introduced.

Additional index storage is capped at 8 MiB per model, allocated during cache
preparation. There is no per-query list allocation. Two passes count and fill
an exact allocation. Empty lists have an explicit readiness flag. Invalid
indices, allocation failure, and insufficient budget leave that universe on
the original scan. The global cache bit is published under the existing build
lock only after preparation completes. Interruption leaves completed lists
available for a retry without publishing an incomplete cache.

Universe invalidation and direct index rebuilds discard the lists. Clones
prepare their own storage, and destruction frees it. Fill and region setters
continue to invalidate through the existing lifecycle. Concurrent geometry
mutation remains unsupported. Preparation and fill-loop cancellation checks
are cooperative, not hard real-time deadlines.

## Correctness and integration

Four new regression tests cover ordered mixed universes, empty child lists,
nested and repeated fills, translation and rotation, aliases, rectangular and
hexagonal lattices, missed parent bounds, undefined fills, geometry edits,
clone independence, direct universe rebuilds, allocation failure, budget
exhaustion, malformed indices, interruption, and subsequent reuse. Prepared
and fallback enumeration compare every returned hit field in nearest and
all-hit modes on synthetic fixtures. Analytic plane cases independently check
distance, alias precedence, normal, clipping, and failure behavior.

The clipping oracle accounts for the diagnostic API's global surface pass:
an unplaced plane remains visible even when a translated occurrence is clipped.
This is breakpoint enumeration, not proof of physical ownership or unique
material coverage.

Final native unit/integration/Lua and serial unit/integration suites pass.
The comparison includes 288 public frames at 512 and 1,024 in split and fused
modes, three alternating process pairs and three frames per process. All
cell/material/density/flag digests match. Median warm control changes range
from -4.3% to +4.4%, within the 5% regression screen; these small differences
are not speedup claims.

## Release measurements

Baseline and candidate use Apple Clang 21 on macOS ARM64 with
`-O3 -march=native -DNDEBUG` and native TinyPar. Processes run sequentially.
The baseline executable was preserved before edits. Both use the existing
one-million-segment and 256 MiB trace-storage budgets. Frame timings exclude
loading, preparation, hashing, and application display/upload.

Three alternating process pairs measure the complete broad 512 split view,
with a first and warm frame per process, separately at one and four configured
workers. The parent timeout is 120 seconds per two-frame process, as in the
preceding baseline warm series. The separate fused first-frame screen keeps
its 60-second timeout.

| Workers | Build | First-frame median | Warm-frame median | Warm range |
| --- | --- | ---: | ---: | ---: |
| 1 | Baseline | 35.028 s | 35.275 s | 34.159–35.405 s |
| 1 | Candidate | 24.814 s | 24.961 s | 24.450–27.841 s |
| 4 | Baseline | 12.489 s | 12.779 s | 12.294–13.615 s |
| 4 | Candidate | 7.932 s | 9.608 s | 8.196–13.517 s |

All 24 paired broad frames match the preceding completed frame in every digest
field, with 88,597 material pixels and no resolution flags. The one-worker
warm median improves by 29.2%; all three candidate warm frames are faster than
all three baseline warm frames. Four-worker first frames also consistently
improve. The four-worker warm median improves by 24.8%, but one candidate run
is slower than its paired baseline, and the ranges overlap. These observations
do not establish a uniform four-worker warm-frame gain or a scaling curve.

The separate one-worker fused first frame completes in 27.048 s with matching
fields. It is a single completion check, not a paired fused timing distribution.
Loading and preparation are additional. Their medians are:

| Workers | Build | Load | Preparation |
| --- | --- | ---: | ---: |
| 1 | Baseline | 6.512 s | 5.157 s |
| 1 | Candidate | 6.650 s | 5.256 s |
| 4 | Baseline | 6.698 s | 5.193 s |
| 4 | Candidate | 6.549 s | 5.341 s |

One-worker candidate preparation ranges from 5.171 to 6.415 s, versus
5.078–5.170 s at baseline; the four-worker-configured processes range from
5.180 to 5.511 s, versus 5.192–5.314 s. These are whole preparation times,
not an isolated measurement of child-list construction. Paired broad process
peak RSS spans about 2.13–2.86 GiB; allocator and process variation prevent
attributing an RSS improvement to this change.

The private storage audit retains 44,344 bytes of indices for 5,543 child entries.
All 198 universe lists are available, including 188 empty lists. The new fields
add 6,144 bytes across the universe vector's allocated capacity on this ABI.
The 8 MiB cap applies to index allocations; struct metadata is additional.


The same 32 sampled private rows at resolution 1,024 take a median 2.739 s
(range 2.710–2.773 s) at baseline and 1.885 s (1.880–2.790 s) with child lists.
The first candidate process did not improve, so retain that variability rather
than claiming every paired row run wins. The median reduction is 31.2%.
A separate 512 pair takes 2.261 s and 1.589 s. All 256 comparison rows match
in digest, segment count, and interval count. Eight private detail frames
also match at 512 and 1,024 in both modes.

Fresh independent ellipse pixel oracles pass at both resolutions and in both
modes. Sixteen existing geometry regressions, twelve ellipse regressions, the
four new child-list tests, and analytic characteristic tracks pass.

The one-second cancellation workload completes eleven rows, cancels seven,
and reports zero query errors, returning after 1.005 s. Full 100 ms requests
return in approximately 102–111 ms; the final shortened-budget cancellation
takes 31.9 ms. Cancelled rows have no digest, and completed rows after
cancellation match their reference outputs. This preserves observed reuse and
cooperative responsiveness, not a hard 100 ms deadline.

## Remaining work in the complete-frame profile

An ignored instrumented copy of the final implementation preserves release
output on four synthetic frames and the complete broad private frame. The
broad diagnostic takes 30.190 s, including 30.188 s of tracing; diagnostic
times are not release throughput measurements.

Fill-loop cell visits fall from 3,343,768,470 to 2,635,648, a 99.92% reduction.
The other recorded work counters are unchanged: 349 actual recoveries,
1,388,145 fill-recursion calls, 2,635,648 positive-fill candidates, 1,247,852
parent-box rejections, 256,516 indexed scans, 815,376,379 indexed references,
and 1,133,375 reused scan contexts. This confirms removal of terminal loop
entries without pruning recursive children or changing the indexed solves.
The equality between child visits and positive-fill candidates is specific to
this measured fill-recursion workload; synthetic tests also cover lattice-only
entries.

Recovery still takes 26.264 s, about 87% of diagnostic tracing time. Indexed
scans take 23.054 s, about 88% of recovery and 76% of total tracing. The measured
ownership-resolution block takes 3.569 s. Nested recovery, scan, and stage
timers overlap and must not be added together. The remaining indexed surface
work is now the clearest performance target.

## Decision and limits

Retain this change. The one-worker complete-frame benefit exceeds the observed
variation, all completed fields agree, independent geometry checks and CPU
suites pass, and public controls stay within the agreed regression screen.
Retained storage and fallback are bounded, and lifecycle/cancellation coverage
passes. The four-worker warm measurements remain more variable and should not
be presented as a guaranteed gain for every frame.

The one-second broad warm-frame target remains unmet. Further acceleration
should examine the remaining indexed surface scans and recovery strategy with
the same ordering, transform, and geometric correctness contracts. No new GPU
backend, 3D-renderer benchmark, display-latency claim, or complete transport
coverage validation is established here.

Only aggregate private observations belong in this report. Inputs, paths,
coordinates, extracted diagnostics, raw hashes, and logs remain outside Git
or under ignored `build/`. The branch remains local and unmerged; this report
is not added to the public documentation site's allowlist.
