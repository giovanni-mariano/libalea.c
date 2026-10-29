<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO
SPDX-License-Identifier: MPL-2.0
-->

# Direct plane calls in nearest indexed scans

Status: retained locally, measured on 2026-10-29 against production `9375ed9`
(documentation checkpoint `05cef84`). The small implementation follows the
[captured-context replay](GEOMETRY_SCAN_COST_REPLAY.md). Three paired runs show
a consistent one-worker complete-frame improvement: median warm time falls
from 23.815 s to 21.491 s, or 9.8%. Four-worker warm results are mixed and do
not establish a repeatable gain. No cache or prepared table is added.

## Implementation and contract

In `raycast_indexed_surface_scan`, nearest-only plane references use the existing
validated payload accessor and call `ray_intersect_plane` directly. The payload
is copied into the local union only if a root lies inside the scan window and
the existing normal calculation needs it. This bypasses an otherwise repeated
union copy and generic intersection dispatch. Other primitive types and ordinary
all-hit queries retain their existing path.

The scalar plane formula, coefficients, parallel threshold, original card
order, aliases, inclusive window endpoints, nearest-hit ties, occurrence
transforms, and entry-bound clipping are unchanged. The scan does not tighten
its upper bound during iteration. Index/payload validation, failure flags,
seen-card handling, tested-surface counts, and cancellation checks retain their
semantics. No public API, persistent allocation, preparation, ownership, or
invalidation mechanism changes. The normal path always receives an initialized
payload before reading it.

## Validation

Native unit/integration/Lua suites and serial unit/integration suites pass,
including characteristic tracks and the existing geometry and elliptical-torus
regressions. Four new focused tests pass on both backends and on the retained
baseline library. They exercise inclusive endpoints and adjacent ULPs in both
directions, near-parallel thresholds, signed zeros, extreme and nonfinite
payloads, invalid indexed cards and missing payloads, transformed normals,
terminal-cell-box misses, cancellation, and subsequent reuse.

The fixture rotates a child plane so that the global ray is parallel to the
unplaced surface. The expected crossing must therefore come from the child
occurrence; a global scan cannot mask a faulty indexed result. Analytic cases
check expected distances and normals. Extreme-payload cases compare root and
normal fields bit-for-bit with ordinary generic/tree enumeration, preserving
its existing behavior rather than defining new rules for malformed geometry.
Existing tests cover aliases, coincident ordering, later nearer occurrences,
repeated contexts, and transformed lattice paths.

Fresh independent ellipse pixel oracles pass at 512 and 1,024 in split/fused
modes. All 288 public comparison frames, 256 private comparison rows, eight
detail frames, twenty-four paired broad frames, and the additional fused broad
frame match their reference outputs. Frame comparisons cover cell, material,
density, and flags. These are bounded compatibility checks, not a proof of
physical ownership throughout arbitrary imported geometry.

## Complete-frame measurements

Apple Clang 21 on macOS ARM64, release `-O3 -march=native -DNDEBUG`, native
TinyPar. Three alternating process pairs at each worker count use the same
broad private 512x512 split view and output budgets. Each process measures a
first and warm frame; processes run sequentially. Frame time includes tracing
and rasterization and excludes loading, preparation, hashing, and display.

| Workers | Build | First-frame median | Warm-frame median | Warm range |
| --- | --- | ---: | ---: | ---: |
| 1 | Baseline | 24.904 s | 23.815 s | 23.105–23.982 s |
| 1 | Candidate | 21.913 s | 21.491 s | 20.570–21.727 s |
| 4 | Baseline | 7.277 s | 7.874 s | 7.004–8.094 s |
| 4 | Candidate | 6.563 s | 7.334 s | 7.203–8.006 s |

One-worker warm pairs improve by 14.2%, 7.0%, and 8.8%. Their ranges do not
overlap. First frames also improve in every pair, although setup and execution
variability make these noisier observations.

Four-worker warm medians differ by 6.9%, but the paired results are two
regressions of 1.7% and 2.8%, followed by a 9.4% improvement. Ranges overlap.
Those results are retained without attributing their variation to a specific
cause. The lower median alone does not justify a reliable four-worker speedup
claim. The separate one-worker fused first frame takes 21.620 s; it verifies
completion and matching fields, not a paired fused speedup.

The 32-row workload at resolution 1,024 has baseline/candidate medians of
1.799/1.789 s, only 0.5% apart. Ranges are 1.784–2.055 s and 1.575–1.853 s.
All three paired totals improve, but their gains vary substantially. The
single 512 row pair takes 1.461/1.318 s. These measurements do not replace the
complete-frame acceptance evidence.

Public warm-frame median changes span -4.5% to +4.6%, within the 5% regression
screen. Loading/preparation medians are 6.567/5.084 s at baseline and
6.597/5.053 s for the candidate with one worker; corresponding four-worker
values are 6.334/5.088 s and 6.891/5.204 s. No preparation work was added.
Paired-process peak RSS spans about 2.47–2.97 GiB; this establishes no memory
reduction. The change adds no persistent storage.

The one-second cancellation workload completes nineteen rows, cancels four,
and reports zero query errors, returning in 1.005 s. Full 100 ms requests
return in about 103–105 ms; the shortened final request returns in 47.1 ms.
Cancelled queries publish no digest, and subsequent successful queries match
the reference. Cancellation remains cooperative rather than a hard deadline.

## Decision and next boundary

Retain this local change for its consistent one-worker complete-frame benefit,
small implementation, and passing correctness/public-regression checks. The
measured benefit does not require a new index or lifecycle complexity. The
four-worker result remains explicitly inconclusive; the replay's 16% reduction
must not be presented as a complete-frame speedup.

The one-second warm broad-slice target remains unmet. Before another performance
implementation, profile the retained complete-frame path again and investigate
why nearest recovery needs so many surface references. Do not infer that wider
typed dispatch, packed metadata, SIMD, or GPU acceleration will pay off from
this plane-only result. Transport geometry validation remains a separate track.

The local plan is `plans/geometry-direct-plane-scans.md`. Candidate snapshots,
executables, raw measurements, and private diagnostics remain ignored under
`build/direct-plane/`. Private inputs, paths, coordinates, geometry, and hashes
are not tracked. Only aggregate observations are recorded here. This report is
not added to the public documentation-site allowlist; nothing is published.
