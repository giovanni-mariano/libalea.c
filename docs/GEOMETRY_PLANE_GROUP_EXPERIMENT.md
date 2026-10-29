<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO
SPDX-License-Identifier: MPL-2.0
-->

# Plane-group rejection experiment: not retained

The candidate was implemented and validated against `c11a95b` (production code
at `9375ed9`), following the [indexed scan attribution](GEOMETRY_INDEXED_SCAN_PROFILE.md).
It is not retained. Repeated complete-frame measurements do not demonstrate a
large enough improvement relative to their variability to justify the added
cache and lookup machinery. Production source is restored byte-for-byte to the
preceding checkpoint. The candidate source, executables, and raw evidence remain
under ignored `build/plane-families/` for local review.

## Candidate and numerical contract

The implementation prepared per-universe groups of at least sixteen planes
with identical finite normal coefficients. Opposite signed zeros were kept in
separate groups. Each group stored minimum/maximum offset representatives and
an auxiliary map from original card position to group. Only nearest recovery
used the filter, after a recovery-context cache miss; the original surface order
and ordinary all-hit enumeration remained intact.

Both endpoint planes were tested with the existing scalar plane kernel. A group
was rejected only if both calls returned parallelism under that kernel's rule,
or both finite roots lay strictly on the same outside side of the scan-entry
window. Endpoint equality remained eligible. For an identical normal and local
ray, the arithmetic involving the offset is a composition of monotone rounded
additions and division by a fixed denominator, possibly reversing order. The
implementation guarded finite inputs, projection, denominator, and endpoint
roots, falling back when uncertain. Fast-math, finite-only assumptions, and
excess-intermediate-precision builds disabled the certificate. No new root
formula, normal approximation, cell-box rejection, or intra-scan bound tightening
was introduced.

Preparation allocations were capped at 8 MiB per model, with separately capped
8 MiB temporary preparation storage and at most 4,096 groups per universe.
A non-inlined nearest-only helper used 4 KiB of rejection flags, without query
heap allocation. On the measured compiler its complete stack frame was 4,272
bytes; the ordinary scan and cached helper frames were 560 and 448 bytes.
Missing, incomplete, oversized, or allocation-failed metadata retained the
original scan. Publication, dependent invalidation, direct index rebuilds,
clones, interruption, and destruction were implemented explicitly.

## Validation

The candidate passes native unit/integration/Lua and serial unit/integration
suites, existing geometry and ellipse regressions, characteristic tracks, and
fresh independent ellipse pixel oracles at 512/1,024 in split/fused modes.
Seven new focused tests pass in both native and serial builds, including final
extensions for repeated contexts, later nearer occurrences, and undefined fills.
They cover endpoint equality and adjacent ULPs, near-parallel thresholds,
signed zeros, extreme scales, uncertain arithmetic, aliases, transformed fills,
rectangular/hexagonal lattices, the terminal-cell-box counterexample, cache
lifecycle, injected allocation failure, malformed/oversized metadata, and the
aggregate model storage cap. Separate build-mode probes confirm strict-mode
rejection and fallback under fast-math and finite-only compilation.

All 288 public comparison frames match their baseline cell/material/density/flag
fields. Median warm timing changes range from about -2.0% to +3.2%, within the
5% public regression screen. All 256 private comparison rows, eight detail
frames, twenty-four paired broad frames, and the separate fused broad frame
match. Cancelled rows publish no digest; subsequent completed queries match.
These checks support the tested geometric and lifecycle behavior, not universal
physical coverage validation of imported models.

## Paired release measurements

Both builds use Apple Clang 21 on macOS ARM64 with `-O3 -march=native -DNDEBUG`
and native TinyPar. Processes run sequentially with identical output budgets.
Three alternating process pairs measure the broad private 512 split view at
one and four configured workers, with a first and warm frame per process.
The two-frame process timeout is 120 seconds. Loading, preparation, hashing,
and application display are excluded from frame time.

| Workers | Build | First-frame median | Warm-frame median | Warm range |
| --- | --- | ---: | ---: | ---: |
| 1 | Baseline | 23.244 s | 23.804 s | 23.332–23.937 s |
| 1 | Candidate | 22.292 s | 23.282 s | 22.198–30.170 s |
| 4 | Baseline | 7.235 s | 7.236 s | 6.698–7.642 s |
| 4 | Candidate | 6.995 s | 7.041 s | 6.961–9.789 s |

Warm medians improve by only 2.2% and 2.7%. The first two candidate pairs improve
at each worker count, but the third candidate's warm frame is slower than its
paired baseline by 26.0% with one worker and 46.1% with four. Those slower runs
are retained in the evidence. Their cause was not isolated; they cannot be
attributed to system noise or discarded to claim a reliable gain. The observed
variation is large relative to the small median improvements.

The same 32-row workload at resolution 1,024 has a median total of 1.778 s at
baseline (range 1.775–1.901 s) and 1.675 s for the candidate (1.654–1.711 s),
a 5.8% reduction. A separate 512 pair takes 1.407 s and 1.319 s. These modest
sampled-row gains do not override the complete-frame acceptance result.
The candidate's separate one-worker fused first frame completes in 23.218 s;
this is an output/completion check, not a paired fused speedup measurement.

| Workers | Build | Median load | Median preparation |
| --- | --- | ---: | ---: |
| 1 | Baseline | 6.358 s | 4.922 s |
| 1 | Candidate | 6.582 s | 5.067 s |
| 4 | Baseline | 6.372 s | 5.008 s |
| 4 | Candidate | 5.960 s | 5.036 s |

These are whole-process preparation measurements, not isolated group-build
costs. The private metadata audit retains 2,575,404 bytes of maps/descriptors
for 2,484 groups in 169 universes, plus 6,144 bytes of universe-struct fields.
The largest universe needs 137 rejection flags. Signed-zero separation changes
the group count slightly from the earlier observational grouping. Paired-process
peak RSS spans approximately 2.46–2.98 GiB; these observations establish no
whole-process memory improvement.

The one-second cancellation workload completes nineteen rows, cancels four,
and reports no query errors, returning in 1.003 s. Full 100 ms requests return
in about 100–102 ms; the last shortened-budget cancellation takes 60.7 ms.
Cancellation remains cooperative. Successful queries follow cancelled ones.

## Work reduction and decision

A final diagnostic copy matches four synthetic frames and the complete private
frame. Its 24.083 s frame time is instrumented, not release throughput. Across 349 actual recoveries it visits 815,382,047 indexed references,
performs 3,010,477 group tests, and skips 145,527,303 per-card plane tests—about
17.85% of visits. Group tests themselves may call the plane solver twice and
are additional work, not included in the skipped per-card count. The number of
in-window normal calculations remains 5,457. Thus the intended work reduction
occurs, but does not by itself establish a worthwhile frame-time improvement.

The agreed acceptance criterion was net release benefit beyond variability,
after paying for preparation, storage, group probes, and map reads. This run
does not meet that criterion convincingly. The implementation and its dependent
tests were removed from the working source after their complete candidate
snapshot was saved locally. No previous production optimization was removed.

This result does not prove that plane grouping can never help. It establishes
that this particular filter is not justified for retention by the current
measurements. Further work should first distinguish primitive-kernel cost from
indexed dispatch/data-access cost, or investigate why recovery requires such
broad scans, before adding another prepared index. The exploratory one-second
target remains a prepared 512x512 broad 2D slice, excluding setup and display.
No new GPU, 3D-renderer, or transport-solver claim follows from this experiment.

Private inputs, paths, coordinates, raw hashes, and diagnostics stay ignored.
Only aggregate results are recorded here; this report is not added to the public
documentation-site allowlist, and nothing was pushed or published.
