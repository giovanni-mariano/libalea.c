<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO
SPDX-License-Identifier: MPL-2.0
-->

# Native GPU execution audit: first iteration

Date: 2026-10-30. Apple M1, four CPU threads. Baseline: `e9231ef` plus the
working-tree CPU undefined-fill fix and benchmark changes. This is an initial
E0 comparison and source audit, not production acceptance or a demonstrated new
speedup. The initial audit did not change GPU execution. The follow-up below adds
an opt-in diagnostic; ordinary intersection calculations remain unchanged.

## CPU/device responsibility map

| Stage | Current implementation | Assessment |
| --- | --- | --- |
| Import, construction, validation | Core and format modules on CPU | Keep on CPU |
| Snapshot analysis/packing | `compute_scene.c`, `compute_render_scene.c`, `compute_spatial.c`, `compute_point.c` | Already CPU preparation; measure its lifetime and eliminated work before expansion |
| Ordered/all-owner points | `select.metal`, `classify.metal`, interpreter in `compute.metal` | Native batch execution; root candidate option, child scans remain ordered enumeration |
| First-visible rays | `render.metal`, `intersect.metal` | Native BVH/hierarchy traversal and analytic CSG intervals |
| Progressive rays | `resume.metal`, command management in `compute_metal.m` | Existing bounded continuations; retain complete-query state on device |
| Output/presentation | Device hit/shading results, host tools/render support | Geometry remains analytic; host handles UI, labels and image encoding |

The library is already split around a computational geometry core. The next
question is how much repeated work that core performs per query. Porting more
CPU subsystems does not address the measured problem.

## Reproducible comparison tool

`benchmarks/native_gpu_compare.py` accepts a JSON array of named variants with
explicit argv and environment overrides. It serializes runs, alternates variant
order by block, removes inherited GPU experiment settings, sets four CPU threads,
limits each process to at most 90 seconds, and terminates its process group on
timeout. Any launch failure, timeout or nonzero benchmark exit stops the run.
Existing evidence directories cannot be overwritten.

It records revision, executable digests, working-tree patch/status, configuration, process times and
raw logs under a Git-ignored `build/` directory. Keep configurations referring to
private models there too. The runner does not turn process times into GPU times,
validate arbitrary executable outputs, or declare a statistical speedup.
Each executable owns its warmups, semantic oracle and within-process timings.

Public example:

```sh
make RELEASE=1 USE_METAL=1 bin/compute/metal/benchmark_gpu_geometry
python3 benchmarks/native_gpu_compare.py \
  --config benchmarks/native_gpu_compare.example.json \
  --output build/native-gpu/comparison-new --blocks 2
```

At this checkpoint the expanded analytic benchmark reports CPU renderer oracle
failures, so this command intentionally stops with failed comparison status.
Do not ignore that exit status to claim acceptance. Fixes must be independently
validated; the report below retains diagnostic observations with that limitation.

Runner checks exercised alternating A/B then B/A order, successful completion,
process timeout/termination, and failure-stop behavior from the real benchmark.

## Fresh public matched operation measurements

The existing 64 translated sphere/void fixture now measures CPU rendering as well
as Metal, alternating backend order within each pass. It retains analytic point
ownership and sphere-intersection expectations, tests one initial and five warm
calls, and checks every rendered pixel outside the timed region, including explicit
nonfinite depth/normal rejection. Both CPU and GPU
ray outputs are checked for IDs/materials, depth and source-oriented normals.

One wide-format audit process produced these warm medians:

| Operation | CPU ms | Metal ms | CPU comparison |
| --- | ---: | ---: | --- |
| All-owner points, 512x512 | 47.297 | 33.021 | FP64 compute enumeration reference |
| Selected points, 512x512 | 37.216 | 28.323 | Prepared core lookup, serial |
| Exterior, 512x384 | 110.826 | 7.984 | Production renderer, four threads |
| Cutaway, 512x384 | 48.678 | 4.700 | Production renderer, four threads |
| Close-up, 512x384 | 163.927 | 12.375 | Production renderer, four threads |

Both point operations pass their oracle. Metal rendering has zero mismatching
pixels in all three views. CPU rendering reports one exterior and two close-up
pixels outside the existing oracle contract; cutaway has none. Across six calls
per view that yields 18 pixel mismatches. Their precise mechanism has not been
classified. The benchmark reports the counts, finishes the diagnostic views and
returns nonzero; the comparison runner stops before another variant/block.

These numbers are a single diagnostic process, not a completed wide/compact
comparison or release timing distribution. The CPU point reference is not the
optimized grid renderer. The few CPU ray discrepancies prevent an unconditional
claim of equivalent outputs, although the independent oracle supports Metal on
this fixture. No performance target for the private workload follows from this
small, regular public scene.

## Private CPU work versus existing GPU attribution

A local CPU-only probe calls the production `render_pixel_solid` helper with the
same stored camera configurations, 4x3 dimensions and pixel centers used by the
saved GPU attribution. It includes the current CPU correction. Each view uses
its own process and completes within 90 seconds. No private model or camera data
is included in this report.

| View, 12 pixels | CPU indexed surface tests | CPU walk iterations | Saved GPU primitive solves | Saved GPU interval steps |
| --- | ---: | ---: | ---: | ---: |
| Exterior | 31,694 | 172 | 2,672,665 | 338,332 |
| Cutaway | 119,467 | 482 | 1,981,525 | 200,871 |
| Close-up | 1,143 | 42 | 1,113,849 | 148,827 |

CPU terminal/ancestor indexed tests are 28,695/2,999 exterior,
37,293/82,174 cutaway and 669/474 close-up. These counters omit work such as
point-containment evaluation and some global searches. GPU solves count repeated
primitive intersection calls. CPU/GPU camera precision and clip-entry treatment
differ, and correctness differences remain open. Thus these are corresponding
view diagnostics, not identical-bit-ray measurements, equal-coverage counters,
CPU total-work counts or direct speed ratios. The GPU observations are the earlier
profiled runs, not fresh timing runs. Their source arithmetic is unchanged by the
subsequent auxiliary-normal fix.

Saved parent-check/cache-hit counts are 44,616/2 exterior, 41,623/2 cutaway and
48,226/0 close-up. Parent checks include eager work as well as lazy validation;
these totals do not define a cache hit rate. Existing attribution places 87-88%
of exterior/cutaway instructions in containers and 48% in close-up. Exact reuse
frequency across placed ray frames has not yet been measured.

## Source evidence and ranked hypotheses

In `render.metal::ray_expression` and `resume.metal::resume_expression_step`, each
executed leaf calls `primitive_roots`. Roots depend on primitive data and the
local ray origin/direction; cursor/end determine which returned root is used.
The current code recomputes roots when a leaf is revisited. Whether a small cache
would capture enough repetition to offset its lookup/state cost needs measurement.
The existing parent cache stores one truth segment for a parent frame; it is not
a primitive-root cache. Merely increasing that cache or adding a second BVH is not
yet justified.

| Rank | Hypothesis | Evidence and next decision |
| --- | --- | --- |
| 1 | Reuse intersection results within a placed ray | Repeated root-solving call sites and large solve counts; measure actual reuse distance and primitive mix before selecting a cache |
| 2 | Prepare smaller residual CSG programs/candidate domains | Root point cost and container ray cost remain large; measure work-weighted proof coverage, including budget-limited and ineligible cases |
| 3 | Change scheduling/queues | Existing continuations work; no hardware attribution yet establishes scheduling as the main limiter; defer extra queue machinery |

The parked bounds prototype's 529 certified cells cannot be judged by cell count
alone. Its work-weighted effect and controlled timing remain unknown. The public
hierarchy fixture's zero certificates also requires an eligibility audit before
concluding its proofs failed. Do not reuse its sequential timings as a speedup.

## One discriminating follow-up

Instrument a bounded diagnostic of exact primitive-root requests in the existing
native ray path, without changing roots or evaluation order. For selected public
rays and the same small private view samples, count:

- eligible requests by primitive family;
- repeat requests keyed by primitive ID and exact local-ray coordinate bits,
  scoped to the snapshot and ray/AA sample;
- hits/evictions for small candidate capacities, with memory accounted explicitly;
- error/nonfinite requests, which cannot be hidden by successful cached values.

Use per-ray bounded scratch or sampled diagnostics; avoid a full private trace or
unbounded hash table. Timing with diagnostic instrumentation is not acceptance
performance. Check instrumented output against the unchanged native renderer.
Retain all raw keys/coordinates privately. Resume and frame-reuse lifetime must
be explicit; array depth alone does not identify a placed ray.

If repetition is substantial in the expensive cases, specify one bounded root
reuse prototype for ordinary/resumable Metal. Preserve original arithmetic,
surface provenance, per-interval sign probes, errors and fallback on eviction.
Otherwise move to the residual-program hypothesis. No cache is promoted from
counts alone: complete frames must pass the plan's 20% improvement and regression,
resource, correctness and cancellation gates.

This first E0 iteration ends with that named measurement. Full E0 closure still
needs the reuse evidence, all-owner work attribution, and a completed controlled
baseline once the new CPU oracle failures are classified. CUDA and navigation
status are unchanged. Continue architecture work on independently validated
fixtures while keeping accuracy failures visible.


## Bounded root-reuse diagnostic

The ordinary Metal ray path now supports `ALEA_METAL_ROOT_REUSE=/local/output.csv`,
read when the context is created. It simulates four independent direct-mapped
caches with 8, 32, 128 and 512 entries. It still executes every original root solve,
interval selection and sign probe. With the option unset, a Metal function constant
specializes the diagnostic away and no diagnostic buffer is allocated.

The key is the packed primitive ID plus the exact six FP32 local-ray coordinate
bits. The immutable snapshot is implicit in context ownership. Different placed
rays cannot falsely match; identical local rays may reuse the same immutable
primitive regardless of placement. Hash collisions require full key comparison.
Keys reset for every pixel/AA sample; counters reset for every frame. Zero-root
successes are eligible. Solver errors and successful calls returning a nonfinite
root are counted separately and excluded from simulated hits/insertion. The
ordinary evaluator's numerical/error behavior is preserved.

Scratch is 22,380 bytes per pixel, capped at 1,024 pixels (21.86 MiB), independent
of cell count. The report contains aggregate counts by primitive family and cache
capacity, allocated scratch bytes, image/sample dimensions and saturation status.
It contains no primitive IDs, keys, coordinates or geometry. Family numbers follow
`compute_layout.h`. Misses count eligible requests minus hits; evictions count a
successful incoming request replacing a different valid key. Saturated counters
invalidate exact accounting. This is a fixed-policy simulation, not a full reuse
trace or an optimal-cache upper bound.

Only successful complete ordinary frames emit a report, overwriting the chosen
file. Check the render return code; an earlier report can remain after a failed
frame. Larger images and progressive rendering explicitly reject this diagnostic
mode. Progressive reuse, 2D points and production cache performance are not measured
by these counters. Unset the variable for normal rendering and timing.

Validation exercises bitwise keys (including signed zero and changed placement
coordinates), collisions/eviction, failed/nonfinite requests, saturation, frame and
AA reset, translated/rotated fills, clipping, report errors, allocation failure and
retry, and explicit size/progressive rejection. Public Metal images are compared
bitwise for IDs, materials, depth, normal and RGB, with independent analytic checks.
The public three-view 16x12 fixture passes in both code formats with seven-pixel
submission batches. Its 162 / 64 / 738 requests have 45.68% / 46.88% / 43.77% reuse
respectively, already captured by eight entries. This differs substantially from
the complex private workload and must remain in the regression set.

### Private result and decision

A fresh ordinary-render pass used the three established 4x3 views, AA1, compact
code, one pixel per submission and one context per process. Separate baseline and
diagnostic processes ran serially with a 90-second limit. These are work counts,
not timing comparisons or resolution-scaling evidence.

| View, 12 pixels | Root requests | 8-entry hits | 32-entry hits | 128-entry hits | 512-entry hits |
| --- | ---: | ---: | ---: | ---: | ---: |
| Exterior | 2,672,665 | 42.24% | 70.75% | 85.64% | 90.74% |
| Cutaway | 1,981,525 | 41.11% | 68.85% | 84.44% | 89.88% |
| Close-up | 1,113,849 | 44.66% | 68.20% | 79.10% | 84.11% |

Every request total agrees with the independent existing ray-profile counter in
both baseline and instrumented runs. All output fields are bit-identical between
them: IDs, materials, depth, normals and RGB. There are no solver errors, returned
nonfinite roots or saturated counters in these samples. Diagnostic scratch is
268,560 bytes per view. Across the final private processes, maximum observed
snapshot preparation was 10.134 s, process time 38.8 s, process peak 3.10 GiB and
Metal buffers 426.61 MiB, within the declared guards. Raw reports and images stay
in ignored local build storage.

Planes account for 71.1–77.1% of requests, quadrics 15.7–24.5%, and tori 0.48–0.78%.
Counts are not cost weights: even a high plane hit rate may lose to hashing and
memory traffic. A 512-entry cache captures only about five more percentage points
than 128 entries while requiring four times the entries. The evidence therefore
supports **one 128-entry root-cache prototype**, beginning with the measured exact
key and admission policy. It does not establish that root solving dominates frame
time or that any cache will pass the performance gate.

The experiment specification is saved in the local architecture plan. It covers
ordinary and resumable rays, retains original interval/sign work and errors,
limits scratch by resident rays, and requires controlled complete-frame timing.
The gate remains a repeatable >=20% improvement with no unexplained >5% regression,
plus resource, correctness and cancellation checks. Stop after at most two measured
revisions if that case cannot be made. No production cache has been implemented
or enabled by this diagnostic, and no CUDA or 2D acceleration claim follows.

Validation of the final implementation passes 74 compute tests, five builder tests
and ten resource/diagnostic tests in each code format. The diagnostic confirms
output preservation, not absolute correctness of the private geometry: existing
CPU/GPU accuracy discrepancies remain open. E0's all-owner attribution and the
controlled baseline blocked by the public CPU oracle failures are still pending.

## Root-cache prototype: implemented, tested, parked

The proposed 128-entry cache was implemented in both ordinary and resumable Metal
rays. Exact keys, finite-root admission, zero-root results, AA/ray resets, persistent
resume slots and original interval/sign work followed the specification. Cache
storage was separate from the copied continuation state. It used 48-byte entries
plus 24 bytes of optional counters per resident ray, with at most 4,096 slots;
resume allocation included the cache in the existing 64 MiB scratch allowance.

The prototype was then **parked and removed from the active source tree** because
it failed the performance gate. Its source, patch, binaries and raw measurements
are preserved in ignored local build storage. The prior renderer and bounded reuse
diagnostic are restored byte-for-byte. No root-cache option is retained in the
active implementation. Existing CPU changes are unaffected.

### Actual work eliminated, but no useful private timing gain

Two serialized baseline/cache then cache/baseline blocks used the established
private exterior 4x3 view, compact code, ordinary rendering and one context per
process. Each process completed one first-use frame and five warm frames. Both
profilers were disabled for timing; setup and output checks were outside frame
timing. Every process completed within 90 seconds.

| Process order | Baseline warm median | Cache warm median | Improvement |
| --- | ---: | ---: | ---: |
| Baseline then cache | 10.995764 s | 10.931805 s | 0.58% |
| Cache then baseline | 11.304467 s | 11.296184 s | 0.07% |

A separate profiled run confirmed 2,672,665 requests, 2,288,805 hits and only
383,860 actual solves: the predicted 85.64% reuse was realized. There were no
solver errors, returned nonfinite roots or saturated counters. Cache scratch for
these twelve rays was 74,016 bytes. All output fields were bit-identical to the
baseline, including RGB and auxiliary buffers. Thus failure to achieve useful
speedup was not caused by the cache failing to find repeated requests.

The timing variation is larger than the observed gain, and neither block comes
close to the required 20%. These observations do not establish a speedup. They
also do not isolate root-solver cost: hash/key comparison, global-memory accesses
and cache initialization can offset the work removed. CSG interpretation and
interval iteration remain unchanged. A high hit percentage is insufficient to
choose an execution representation.

### Public regressions in both execution modes

The independent 64-placement analytic fixture used 512x384 images, AA1, three
views, one first-use and five warm frames per view. Ordinary and progressive
baseline/cache variants ran in forward then reversed process order. Every frame
passed the analytic oracle. Cache/baseline output buffers matched bitwise in both
modes. The ranges below are the two per-process warm medians, not percentile
estimates or confidence intervals.

| Operation | Baseline warm medians, ms | Cache warm medians, ms |
| --- | ---: | ---: |
| Ordinary exterior | 8.43–8.87 | 55.21–66.27 |
| Ordinary cutaway | 4.61–4.63 | 53.93–62.58 |
| Ordinary close-up | 11.39–15.07 | 58.64–71.77 |
| Progressive exterior | 52.16–75.27 | 87.21–99.48 |
| Progressive cutaway | 48.95–57.51 | 84.48–118.03 |
| Progressive close-up | 59.79–69.85 | 97.79–114.90 |

Ordinary cache runs required bounded 4,096-pixel submissions instead of the
baseline's whole-image submission on this small scene. The comparison therefore
includes the cost of that required storage policy; it cannot attribute the whole
regression to root lookup. Progressive runs retained 4,096 resident slots in both
variants and also regressed. Cache scratch added 24.09 MiB at that capacity.

The native cache passed 74 compute tests in each code format, including existing
hierarchy/lattice/error and cancellation checks. Thirteen resource/diagnostic
tests passed, including new bitwise ordinary/resumable output, AA/frame reset,
allocation retry, invalid configuration and error/report handling tests. The
restored renderer passes the original 74 compute, five builder and ten resource
tests. Maximum private process time was 84.4 s, snapshot preparation 10.87 s,
process peak 3.29 GiB and Metal buffers 410.95 MiB, within the declared guards.

### Decision and continuation

One measured cache revision was sufficient to reject promotion: the declared
exterior target failed and every public view regressed. The optional second
revision was not used. Private cutaway/close-up timing and private progressive
cache timing were not expanded after this failure; earlier reuse counts for
those views remain diagnostics, not cache performance results.

The next investigation returns to executed CSG work and preparation coverage:
use existing 2D all-owner and 3D per-cell profiles to rank expensive original and
spatial programs, quantify work-weighted specialization coverage, and identify a
bounded way to reduce interpreted instructions or interval visits. If those
profiles cannot distinguish memory traffic from interpreter/solver cost, obtain
that attribution before another implementation. Do not add more caches or port
another CPU subsystem based on repetition counts alone.

The independent private accuracy blocker, public CPU oracle discrepancies,
full-resolution performance targets and hardware-gated CUDA work remain open.

## Parallel-query scaling follow-up

The [query-scaling and camera-motion evaluation](NATIVE_GPU_SCALING.md) measures
fixed point/ray workloads across submission sizes, records native command times,
and separately tests preview/detail and cancellation. Larger batches matter
substantially, including on the private slice. A larger native point-cap diagnostic
also triggered Metal's Impacting Interactivity error, so no production cap was
changed. The private exact-ray preview remains far outside the motion target.
The new all-owner work profile quantifies 3.67 billion executed instructions and
99.82% rejecting work on the fixed slice. Aggregate attribution is available;
work-weighted preparation opportunities still need classification.

## Motion-preview representation follow-up

The [cached-depth reprojection evaluation](NATIVE_GPU_REPROJECTION.md) adds only a
benchmark-local native display kernel. Projection is fast on the public fixture,
but sparse source frames and newly exposed surfaces prevent adoption as the
standalone motion strategy. No production renderer or viewer default changed.

## Prepared-bounds follow-up

The [native bounds-raster experiment](NATIVE_GPU_BOUNDS_PREVIEW.md) completed with
sub-millisecond median draws, but private bounds clutter and full-scene preparation
prevent a useful standalone first preview. It remains benchmark-only. Return to
work-weighted specialization coverage before adding another display representation;
useful geometry generation and analytic traversal remain the unresolved costs.

## Work-weighted specialization result

The [per-cell work audit](NATIVE_GPU_WORK_AUDIT.md) now distinguishes call counts
from executed cost. Outside-domain fallback consumes 50.30% of the fixed private
slice's instructions, terminal programs 14.21%, and unproven finite certification
only 2.30%. Public accepted outside-domain calls rule out treating those domains
as rejection bounds. Work is spread across many cells. A widespread static chain
of identical short-circuit skips makes preparation-time jump threading the next
bounded candidate; no speedup or changed production behavior is claimed yet.

## Jump-threading implementation result

The [opt-in bytecode pass](NATIVE_GPU_JUMP_THREADING.md) is implemented and measured.
It preserves retained private outputs exactly, reduces fixed-slice instructions by
40.4% without changing primitive counts, improves slice time by 32.6–35.6%, and
improves the bounded exterior ordinary-ray time by about 22%. Cutaway, close-up and
resumable views also improve. The switch remains disabled by default: public
microbenchmark order sensitivity and a paired preparation-overhead result leave
the promotion gates unresolved. See the report for scope, validation and commands.
