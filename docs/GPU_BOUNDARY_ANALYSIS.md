<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO
SPDX-License-Identifier: MPL-2.0
-->

# GPU rendering experiments and accuracy findings

Status, 17 October 2026: the uniform-region compiler and surface-patch renderer
were removed after failing their performance checks. Their dedicated trace/replay
code, tests, benchmark targets, public patch API, and viewer options were also
removed. This report preserves the findings; it does not describe available tools.
The working CPU/Metal ray renderer, profiling, and its regression tests remain.

## Retired surface-patch renderer

The prototype rasterized conservative screen rectangles and evaluated analytic
surfaces and Boolean trimming masks in fragments. It required manual display
bounds and expanded bounded lattice placements.

On the development Apple M1, ITER blanket cell 427004 at 256x192 took
1,006/1,011 ms per warm patch frame versus 102/100 ms for ray rendering.
Preparation took 0.583 seconds and the cache occupied 9.426 MiB for 99,667
surface patches, 10,049 cap regions, and 639 terminal instances. All 49,152
pixels agreed in IDs, depth, normals, and RGB, but only 350 were solid hits.
The tested bounds were `180,420,180,420,-260,-140`, with a central Y cut and
four subdivisions per axis. Broad projected coverage and trimming work made
this implementation roughly ten times slower on that view.

## Retired regional analysis

The host analyzer split automatic universe domains at the midpoint of their
longest axis. It used conservative interval proofs, deduplicated residual CSG
programs, and bounded compilation work and storage. Regions were processed in
deterministic tree order; subdivision was not directed by measured regional cost.

Bounded captures recorded selected progressive rays. Replay checked their
original expression evaluations, then applied regional programs only where
proof coverage and recorded predicates were sufficient. Missing predicates,
uncovered intervals, and exhausted budgets retained the original work count.
Replay preserved the original expression checkpoints; it did not predict a
changed traversal schedule or frame time. The continuation threshold was a
fourfold reduction in replayed instructions without extra primitive solves.

## Recorded ITER experiment

The model was `E-lite_R250630.mcnp`; the packed-scene fingerprint was
`668f641b19166421`. Four complete AA1 traces captured pixels 0 and 85 of the
16 by 12 cutaway, pixel 85 of the opposite cutaway, and pixel 85 of the exterior
view. All four progressive images matched ordinary Metal rendering in cell and
material IDs, depth, normals, and RGB. Diagnostic timings are excluded from
performance conclusions.

Selection used the archived full-frame reference cell profile. Cell IDs,
materials, and original program sizes were checked; that legacy profile does
not carry the new geometry fingerprint. Universes 0 and 20040 account for
133,463,230 of its 151,820,554 instructions, or 87.91 percent. The traces cover
four rays and cannot establish full-image performance or spatial coverage.

Each run used a 32 MiB cache budget and a 67,108,864 source-instruction proof
budget. The domain/compiler scratch estimate was 45,158,361 bytes, separately
from trace input, replay storage, original geometry, and the persistent cache.
The final depth-six host run used 765,606,968 bytes of packed geometry, completed
analysis in 2.404 seconds after loading and preparation, and measured a peak
process RSS of 3,157,884,928 bytes (about 2.94 GiB), including the original MCNP
model, packed scene, trace input, and analysis allocations.

| Binary subdivision depth | Regions | Usable local regions | Budget-limited regions | Replay instructions | Baseline divided by replay |
| --- | --- | --- | --- | --- | --- |
| 0 | 2 | 2 | 0 | 2,499,746 | 0.996x |
| 3 | 16 | 7 | 9 | 2,369,852 | 1.051x |
| 6 | 128 | 49 | 79 | 2,282,330 | 1.091x |

The replay baseline was 2,490,900 instructions and 453,145 primitive solves.
Depth six used 415,651 solves, with zero proof disagreements. It still retained
84.42 percent of baseline instruction work as fallback, including 4,186 windows
whose alternative programs required uncaptured predicates. No selected ITER
region was proved empty, uniform, or a direct analytic boundary. Successful
depth-six local programs had median/p95/maximum sizes of 67/429/8,668 instructions.
The cache filled its 32 MiB allowance. These results identify substantial residual
complexity, program storage pressure, and limited replay coverage; they do not
rule out a better region-selection or proof strategy.

The sparse synthetic workload did produce eight direct-boundary
candidates in its child universe, showing that the simple case is reachable.
The removed analysis fixtures checked empty space, uniform ownership, thin shells,
overlapping owners, fills, repeated lattices, budgets, interruption, malformed
traces, and ray-identity requirements.

## Primitive root reuse

The four captures made 453,145 primitive requests for 27,677 distinct inputs,
where an input is a primitive plus its bitwise local ray origin and direction.
Distinct counts are summed per capture. Each simulated LRU starts empty for
each ray, with no reuse across captures or AA samples.

| Simulated entries | Hits | Hit fraction |
| --- | --- | --- |
| 8 | 210,989 | 46.56 percent |
| 16 | 333,618 | 73.62 percent |
| 32 | 385,588 | 85.09 percent |

Planes account for 334,924 requests and quadrics for 98,936. Avoiding cheap plane
solves may not repay lookup and continuation-state traffic. A future bounded root-cache experiment would need to compare all-primitive and
expensive-primitive policies using actual frame time. No root cache was implemented.

The boundary compiler needs a better locality and storage result before GPU
integration. A next compiler experiment should preserve existing adaptive
specializations, prioritize regions using measured work, and measure the extra
predicate data needed for replay. Increasing the uniform depth or memory cap
alone is not supported by this experiment.

## Retained evidence

Local diagnostic outputs remain under `build/`: `boundary-iter-final-d0-*`,
`boundary-iter-final-d3-*`, `boundary-iter-final-d6-*`, and
`boundary-iter-final-memory-*`, along with the four `boundary-iter-*.trace`
captures and existing patch benchmark images. These ignored artifacts are local
records, not a maintained replay interface.

Before cleanup, the source and documentation were saved locally in
`build/retired-gpu-experiments.tar.gz`, since the experimental source had not
been committed. Current builds no longer include those experiments. Their old
implementation plan was removed along with the runnable reproduction commands.

## Consolidation measurements, 17 October 2026

The manual 4x4x4 root-universe specialization was removed. Automatic adaptive
specialization remains, including its conservative primitive proofs. Removing
the obsolete per-cell field saves about 1.72 MiB on 450,648 ITER cells; this is
not a measured frame-time claim.

Candidate masks were tested with automatic spatial programs, selective lazy
clipping and slot refilling enabled, and rebuilt BVHs disabled. Only masks were
toggled. Each view used one warm-up and five alternating warm frames at 16x12,
AA1, on Apple M1. All output channels matched exactly.

| View | Masks off, median (range), s | Masks on, median (range), s | Change |
| --- | --- | --- | --- |
| Cutaway | 19.783 (19.642–19.981) | 20.544 (20.437–20.815) | 3.85% slower |
| Opposite | 21.808 (21.522–22.072) | 22.764 (22.404–22.967) | 4.38% slower |

The cache used 142 domains and 2.26 MiB of metadata, generated 140,392 reduced
masks and exhausted its 32-million-proof budget. Candidate masks and their
controls, buffers and profile column were removed. These measurements isolate
the old mask implementation; they are not a measured before/after speedup from
the subsequent structural cleanup. Weighted BVH rebuilding remains separately
opt-in.

Local records: `build/candidate-gate-cutaway.log`,
`build/candidate-gate-opposite.log`; diagnostic executable FNV1a64
`2441b2fb16c68ef0`. The source model was E-lite_R250630.mcnp, 410,298,522 bytes.
The pre-consolidation source is archived in
`build/gpu-baseline-before-consolidation.tar.gz`.

## ITER accuracy triage

The consolidated renderer was compared at 32x24, AA1, eye (1800,-2200,1400),
target (0,0,0), retaining Y>0. There were 14 cell-ID disagreements and 11
matching-ID depth differences exceeding 0.1 model units. Diagnostic timings
were collected alongside other GPU activity and are not performance evidence.

Seven rays were checked with double-precision point ownership queries and full
CPU interval queries restricted to a window around the disputed depths. These
checks establish local interval behavior, not absence of an earlier hit along
the entire camera ray. Pixel indices are zero-based, row-major.

| Pixel | GPU cell / depth | Bounded full CPU interval entry | Finding |
| --- | --- | --- | --- |
| 25 | 524241 / 3438.385498 | 524241 / 3438.382886 | CPU viewer depth is late at 3438.496826 |
| 63 | 960740 / 2522.746094 | 960740 / 2522.745980 | CPU viewer skips earlier material, reports 960471 |
| 94 | 962127 / 2796.740723 | 962127 / 2796.738045 | CPU viewer skips earlier material, reports 962041 |
| 491 | 151000 / 4261.583008 | 151000 / 4261.582500 | CPU viewer skips earlier material, reports 150906 |
| 634 | 962398 / 3252.429688 | 962398 / 3252.430351 | CPU viewer skips earlier material, reports 962377 |

Pixel 172 differs between cells 307927 and 307933, both material 110. The full
CPU query finds only about 0.000386 model units of the former after the clipping
plane. This is consistent with sensitivity to FP32 boundary rounding, but the
precise arithmetic source has not been established.

Pixel 266 is a material discrepancy: Metal shows cell 123045 (material 400,
water), while CPU interval and point queries find 124939 (material 100, steel).
Both cells belong to universe 161 under parent 459162. The disputed depth is
3664.2568359375. At the corresponding local point
(863.57041810763246, 152.27170774538098, 207.94093005024814), a flat two-cell
Metal classifier reports **two owners**, while the double CPU evaluator reports
steel alone. The same holds at offsets ±0.001, ±0.01 and ±0.1 along the ray.
This reproduces the discrepancy without hierarchy traversal or parent clipping;
the discrepancy was traced to quadric packing. Rounded coefficients left tiny
nonzero axial terms, causing the stationary-point solver to choose origins near
(138518784, 37116004, -52353060) and (147006208, 39390200, -55560872).
Subtracting these origins destroyed local FP32 detail.

The packer now stops solving translation directions whose pivot is at most
32 times FLT_EPSILON relative to the quadratic coefficient scale. It retains
all quadratic coefficients and residual linear terms: no surface term is dropped.
The resulting origins stay near the geometry. All seven isolated samples now
classify steel alone. A maintained regression uses both original quadrics and
inside/outside points; it fails with the old packing threshold and passes with
the fix. The full 32x24 progressive ITER render also changes pixel 266 from
123045 to 124939 at the same reported depth, matching the bounded CPU interval.
This addresses a concrete conditioning defect, not every FP32 error.

Records are local ignored artifacts: `build/iter-accuracy-probes.log`,
`build/iter-accuracy-gpu.bin`, and `build/iter-pixel266-isolated.log`.
The fixed isolated result is in `build/iter-pixel266-fixed.log`.
Do not treat CPU viewer agreement, or
agreement between two Metal configurations, as an independent accuracy oracle.

After consolidation and the quadric fix, 46 Metal-build tests, seven CPU compute
tests, one interrupted
BVH builder test and one resource test covering six upload-allocation failure
points passed. SDL worker/publication, controls, save and cancellation checks
also passed. These tests do not resolve the full-model accuracy cases above.

The accuracy-fix capture is `build/iter-accuracy-fixed-gpu.bin`, with its CPU
comparison in `build/iter-accuracy-fixed.log`. The full-model verification was
not a controlled timing run. Historical timings above predate the packing fix;
repeat the timing baseline before claiming performance for the final renderer.

The final comparison has **10 cell-ID disagreements**, down from 14 before the
quadric fix, and still 11 matching-ID depth differences above 0.1. Four pixel IDs
changed (266, 330, 440, 503), all to the CPU viewer's IDs; only pixel 266 received
the independent two-cell investigation described above. The remaining cases
must not be assumed to be GPU errors or CPU errors without further checks.

## Completed disputed-ray checks, 17 October 2026

All 21 pixels flagged after the quadric fix (ten ID differences and eleven
matching-ID depth differences above 0.1) now have bounded double-precision
interval checks. Six of the ten ID disagreements and all eleven large depth
comparisons involve the CPU viewer returning a later material entry. Metal
matches the earlier material found by the full CPU interval query. This includes
pixel 667, where the CPU viewer's returned depth is over 105 units later.

The table reports original-camera double-precision reference entries and the
saved FP32 GPU frame. Positive error means Metal reports a later entry. The
largest absolute local entry error in these 17 cases is **0.021598 model units**
at pixel 26; this is a measured sample maximum, not a global accuracy bound.

| Pixel | Earlier material cell | CPU viewer depth | Metal depth | Full interval entry | Metal − interval |
| --- | --- | --- | --- | --- | --- |
| 25 | 524241 | 3438.496826 | 3438.385498 | 3438.382886 | +0.002612 |
| 26 | 524220 | 3547.346924 | 3538.398926 | 3538.420524 | -0.021598 |
| 28 | 961901 | 3269.797852 | 3266.667725 | 3266.666519 | +0.001206 |
| 48 | 751019 | 3500.536133 | 3478.289062 | 3478.293987 | -0.004925 |
| 52 | 751012 | 3292.519531 | 3277.352783 | 3277.354658 | -0.001874 |
| 62 | 961901 | 2443.663330 | 2442.408203 | 2442.408938 | -0.000735 |
| 63 | 960740 | 2533.741943 | 2522.746094 | 2522.745980 | +0.000113 |
| 93 | 962127 | 2777.095215 | 2776.390625 | 2776.389898 | +0.000727 |
| 94 | 962127 | 2808.625977 | 2796.740723 | 2796.738045 | +0.002678 |
| 299 | 586245 | 3759.628662 | 3759.399902 | 3759.399813 | +0.000090 |
| 381 | 962109 | 2869.158447 | 2867.942627 | 2867.942967 | -0.000340 |
| 491 | 151000 | 4302.787598 | 4261.583008 | 4261.582500 | +0.000508 |
| 572 | 962763 | 3215.534668 | 3214.682129 | 3214.681386 | +0.000743 |
| 603 | 159254 | 2831.628906 | 2827.143555 | 2827.136927 | +0.006628 |
| 634 | 962398 | 3256.990479 | 3252.429688 | 3252.430351 | -0.000663 |
| 636 | 962328 | 3303.626709 | 3297.603516 | 3297.603307 | +0.000209 |
| 667 | 962769 | 3362.766357 | 3257.134277 | 3257.134190 | +0.000088 |

The other four cases occur at the Y clipping plane. The full CPU query finds a
very short interval of one cell followed by the cell Metal selects; both have
material 110. Metal and CPU depths differ by at most one FP32 ULP in these
cases, but their cell IDs differ.

| Pixel | CPU → Metal cell | Material on both sides | CPU interval width after clipping |
| --- | --- | --- | --- |
| 172 | 307927 → 307933 | 110 | 0.000386310 |
| 297 | 310509 → 310511 | 110 | 0.000459670 |
| 394 | 309265 → 309273 | 110 | 0.000465840 |
| 492 | 309283 → 309282 | 110 | 0.000457540 |

These four cases are consistent with FP32 loss of a thin clipped interval.
The individual contributions of camera rounding, transform rounding and boundary
evaluation have not been isolated. They remain ownership limitations, even when
material coloring happens to agree; do not silently relabel them as exact hits.

The checks used point ownership on either side of disputed depths, then full
interval queries around those depths. Later queries restricted the window to
Metal depth ±0.05 (clamped at the clipping plane), avoiding traversal of the
entire CPU/GPU depth gap. These windows establish that an earlier material
interval exists and its local boundary. They do not prove that no still-earlier
material exists elsewhere along the complete ray, and the unchanged pixels were
not exhaustively revalidated. Normals were not independently checked in this
full-model diagnostic.

Evidence: `build/iter-accuracy-all.log` (intentionally stopped after useful early
windows), `build/iter-accuracy-remaining.log` (completed), the earlier
`build/iter-accuracy-fixed.log`, and cached CPU/GPU ID/depth frames. These are
local diagnostic artifacts; no tracing subsystem or runtime workaround was added.
The CPU viewer's fast entry path needs its own regression and repair; using it
as the accuracy oracle would incorrectly condemn the earlier GPU hits.

## Corrected-renderer timing baseline

Measured after the quadric fix at source commit `a3e66ec`, on Apple M1 with
16 GiB RAM. Each view used a fresh process/context, one completed warm-up, then
five timed progressive frames at 16x12, AA1. Defaults were explicit by clearing
all overrides: automatic spatial programs, selective lazy parents and refilling
on; rebuilt BVHs and profiling off. Diagnostic ray queries had finished before
these serial GPU runs. This is a baseline for the corrected renderer, not a
controlled attribution of speedup to cleanup or numerical changes.

Scene loading and preparation are separate from frame time. Timed frames include
GPU submission, synchronization, readback and progress callbacks; they exclude
SDL display and PNG encoding. Every timed frame is compared bitwise against its
warm-up in cell IDs, material IDs, depth, normals and RGB. These repeated-frame
checks establish determinism, not independent geometric correctness.

| View | Median frame, s | Frame range, s | Median first pixels, s | Largest warm callback gap, ms |
| --- | --- | --- | --- | --- |
| Cutaway | 19.511 | 19.460–19.788 | 0.581 | 27.8 |
| Opposite cutaway | 21.571 | 21.536–22.698 | 0.328 | 56.2 |

Cutaway load/preparation: 8.224/10.574 s. Opposite-view load/preparation:
7.758/9.811 s. Both views produced identical output on all five timed frames.

The tiny 192-ray workload still takes tens of seconds and does not demonstrate
GPU utilization at ordinary image sizes. It remains unsuitable for smooth full
ITER interaction. These timings do not establish benefits from primitive-root
caching; any prototype must beat this baseline without changing the independent
geometry references or silently erasing the documented ownership limitations.

Local reproduction used `build/render_default_baseline.c`, a diagnostic adaptation
of the maintained progressive benchmark: the ordinary-render comparison was
replaced with five timed repeats and bitwise comparison to the warm-up. It is
not a new maintained benchmark interface. Executable SHA-256:
`7b751280e5ed5a6fdb65528a92d709ef046cb85ae0aa907dddd6712261aab2b0`.
Commands were `build/render_default_baseline MODEL 16 12 0 1` and then
`build/render_default_baseline MODEL 16 12 0 2`, with the same E-lite model as
above. Logs: `build/renderer-baseline-view1.log` and
`build/renderer-baseline-view2.log`.
