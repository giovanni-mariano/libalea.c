<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO
SPDX-License-Identifier: MPL-2.0
-->

# Conservative torus rejection in nearest recovery

Implementation following the [remaining recovery profile](GEOMETRY_REMAINING_RECOVERY.md),
compared with production code at `5a43f40` (documentation checkpoint `0319988`).

## Change and compatibility

Nearest-only indexed enumeration can reject a finite torus before invoking the
quartic solver. The internal `ray_torus_window_may_hit` helper checks the axial
slab first, then the two radial slabs. A linear coordinate stays between its
values at the two window endpoints. Rejection requires both endpoints to lie
strictly beyond the same expanded slab. The slab encloses the complete primitive
surface, including roots outside any containing Boolean region.

Endpoint evaluation uses FMA and a deliberately loose guard of 64 machine
epsilons times the largest relevant origin, center, endpoint product, extent,
or one. The guard covers rounding in endpoint and extent calculations. The
check uses no reciprocal direction or near-parallel cutoff. Boundary and
ambiguous cases execute the existing solver. Non-finite inputs, invalid windows,
non-positive widths, horn/spindle/degenerate tori, extent overflow, and very thin
tubes retain that solver. The thin-tube gate is 64 times the square root of
machine epsilon times major radius: the existing quartic can lose the squared
tube radius or underflow it and publish roots that a geometric envelope would
otherwise reject. A regression case preserves that behavior through fallback.

An existing limitation was identified during implementation: the ray quartic
uses the circular minor radius in all directions, while the geometry evaluator
and boxes also represent an elliptical axial semiwidth. The axial envelope uses
the larger width, enclosing both surfaces. This preserves existing ray behavior;
it does not repair or validate elliptical-torus intersection physics. That
solver/evaluator discrepancy remains a separate correctness issue for future
physics readiness.

Only nearest recovery uses rejection. Unpruned tree enumeration and ordinary
all-hit intersection paths keep their existing behavior. Surface aliases,
ordering, selected distances, normals, transforms, and public APIs are preserved.
The change adds no persistent geometry storage or preparation pass. Existing
interrupt checkpoints remain active even when an intersection is skipped.

## Validation

Thirteen focused geometry tests pass, including four added tests covering:

- Independently specified four-crossing rays on all torus axes; misses,
  tangencies, ULP perturbations, and short windows around known roots.
- Tiny nonzero directions, large coordinates, both elliptical-width directions,
  non-finite/overflow/degenerate inputs, and the thin-tube solver fallback.
- Rotated and translated occurrences, forward/reverse queries, and equality of
  distance, surface identity, and all normal components with unpruned trees.
- A deterministic 9,216-ray sweep over all axes and six scales from 1e-9 to 1e6,
  checking that rejection never removes an existing solver root in the window.

The occurrence fixture exercises mathematical breakpoint enumeration with
intentionally overlapping diagnostic parents; it does not establish unique
physical ownership. All independent analytic crossings in these tests are
circular tori. Non-circular cases test envelope compatibility, not ellipse roots.

Public complete slices match every pixel field at both 512 and 1,024 resolution,
in split and fused modes, over three alternating process pairs. Warm total-time
changes across the tested controls range from approximately -2% to +2%, below
the 5% regression threshold. These small timing differences are noise-scale.

All 64 initial private row samples match their previous geometry digests,
segment/interval counts, snapshot counts, and ownership query counters. After
adding the thin-tube fallback, three alternating baseline/candidate pairs at
resolution 1,024 complete all 192 row queries with exactly matching results.

Median total time for the same 32 sampled rows falls from 10.126 s to 6.893 s:
approximately 32% less tracing time, or 1.47 times the throughput for this
workload. Baseline totals range from 10.104 to 10.205 s; candidate totals from
6.779 to 7.195 s. Pooled row medians are approximately 369 ms and 247 ms.
This is a sampled recovery improvement, not a full-frame speedup claim.

Preparation remains approximately five seconds in these runs and has no new
work in this change. Process RSS variation does not establish a memory benefit.
The private model, coordinates, raw output digests, logs, scripts, and binaries
remain ignored under `build/torus-recovery/` or outside the repository.

The full native unit/integration/Lua suite and serial unit/integration suite
pass. Complete material-containing private slices match every field at both
resolutions, with approximately 91.7% material pixels and no flagged pixels.
These are single cold-frame completion checks, not repeated frame benchmarks.
The broader focused private view at resolution 512 still exceeds its 60-second
whole-process screen without completing the trace; no full-field comparison
or final peak-memory claim is made for that cancelled process.

The final-build cancellation check requests 100 ms on expensive rows and observes
approximately 103–105 ms returns. Its one-second tracing budget finishes in
1.002 s, with seven completed rows, nine cancelled rows, and no query errors.
All completed row hashes match the previous checkpoint, including successful
queries after cancellations; no partial digest is published for cancelled rows.
Cancellation remains cooperative.

Retain this bounded optimization: it has a repeated sampled recovery benefit,
matching completed outputs, passing independent/control checks, and no new
persistent storage. The broader view is still not interactive. The identified
elliptical-torus solver/evaluator discrepancy should receive a separate
correctness investigation before broader physics-readiness claims.

The subsequent [elliptical torus correction](GEOMETRY_ELLIPTICAL_TORUS.md) resolves
the identified finite-width solver/evaluator discrepancy and updates the axial
envelope, normals, and analytical slice curves accordingly.
