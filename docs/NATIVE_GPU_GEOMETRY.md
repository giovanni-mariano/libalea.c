<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO
SPDX-License-Identifier: MPL-2.0
-->

# Native GPU geometry: Metal classification

Implemented on `experiment/native-gpu-geometry`, starting from `6cc16d8`, on
2026-10-29. The local roadmap is `plans/native-gpu-geometry.md`.

Metal, CUDA, and future backends own their kernels, device representations,
memory management, scheduling, synchronization, and numerical capabilities.
The common boundary is libalea's C query contract and semantic tests. There is
no GPU portability framework or common device execution layer. The current
`compute_*` packed structures are private to Metal, not a future CUDA ABI.
CUDA implementation and validation require an actual NVIDIA device and remain
future work; no unavailable backend is advertised as implemented.

## First hierarchy increment (`97ed101`)

`alea_compute_classify()` now supports nested transformed fills on CPU FP64 and
native Metal FP32. Each point counts terminal owning occurrences from universe
zero, including terminal void, under all containing ancestors. Two overlapping
placements of the same cell definition count twice. Container cells do not
count as terminal owners. A gap in a container's contents remains a gap.
Unreferenced universe definitions do not participate. IDs are returned only
when exactly one occurrence owns the point.

Flat classification keeps its existing implementation. The new Metal kernel
uses one thread per point and a bounded stack of 32 universe levels. A compact
eight-byte link record per source cell provides original-order universe
enumeration and fill descent. It shares primitive/CSG buffers and the transform
upload with the existing 3D renderer, without duplicating the model or flattening
placements. Query metadata includes void; the visibility-only render BVH does
not, and remains unchanged. Render-specific spatial programs are not used as
point-ownership predicates.

Query links upload lazily, with ownership published only after all allocations
succeed. Their host copy is then released. Classification and ordinary or
progressive rendering may alternate on the same owned snapshot, including after
the source model has been destroyed. Concurrent use of one context remains
unsupported. New query buffers retain at most 4,096 points/results (96 KiB);
geometry, hierarchy links, shared transforms, and caller arrays are additional.
Classification checks cancellation between submissions. Work per submission
depends on the scene, so this is not a latency guarantee.

This first traversal deliberately enumerates each visited universe. Spatial
pruning needs a separate FP32 ownership argument and measurement. The CPU FP64
reference independently enumerates the hierarchy; it is not the accelerated
core point-query algorithm. CPU cycles/depth exhaustion fail at query time;
Metal packing rejects reachable fill cycles and excessive depth before dispatch.
Nonfinite device evaluation, transformed-coordinate overflow, and count/stack
exhaustion produce errors, not apparently valid ownership results.

Lattice classification is explicitly rejected in this increment. Existing 3D
rectangular-lattice rendering remains supported. There is no plotter/viewer
integration change, ray-interval slice implementation, CUDA implementation, or
transport-grade GPU query claim. FP32 can still alter boundary ownership and
lose small features. Point sampling does not prove complete coverage or preserve
features between samples; these limitations are distinct from traversal defects.

## Validation

The CPU-only compute build passes ten tests. The native Metal compute build
passes fifty-two tests, the builder test, and two resource-lifecycle tests on
an Apple M1. Existing primitive, Boolean, lattice-rendering, clipping, normal,
ordinary/progressive, cancellation, and numerical regressions remain passing.

New analytic fixtures cover nested translation and rotation, ancestor clipping,
terminal void, gaps, repeated overlapping placements, unreferenced definitions,
snapshot lifetime, buffer growth/reuse across 4,096-point submissions, interruption
and reuse, and the 32/33-level boundary. Separate tests alternate classification
and rendering in both orders, comparing 3D output with an independently expanded
world-space reference, including clipping. Injected upload failures verify retry,
retained shared-transform ownership, and destruction without leaked buffers.

No core geometry source changed; the earlier CPU recovery fixes are retained.
The broader core suites were not repeated for this optional compute-module
increment. GPU tests ran on a real device, not a CPU emulator. Apple's missing
Metal compiler component was installed before shader compilation and testing.

## Public synthetic measurements

`benchmarks/gpu_geometry.c` constructs 64 translated occurrences of one sphere
and terminal-void universe. A 512x512 planar grid checks ownership against
independent circle equations. Exterior and cutaway frames at 512x384 check
first-hit material/cell, depth, and normals against analytic sphere-ray roots.
All measured outputs pass those checks.

Apple M1, Apple Clang 21, native release build, `ALEA_NUM_THREADS=1`. One first
call and three warm calls per workload, with alternating CPU/Metal point-query
order, in one process. Times include public API validation, conversion, command
submission, synchronization and output copies. Oracle checks and display are
excluded. Scene creation is outside the timers. Context creation takes 0.061 ms
for CPU and 57.912 ms for Metal; lazy resources are included in first-call times.

| Workload | Backend | First | Warm median | Warm range |
| --- | --- | ---: | ---: | ---: |
| 512x512 point-sampled slice | CPU enumeration | 165.955 ms | 154.607 ms | 152.807–155.318 ms |
| 512x512 point-sampled slice | Metal | 34.720 ms | 30.774 ms | 29.994–32.854 ms |
| 512x384 exterior | Existing Metal renderer | 14.956 ms | 8.030 ms | 7.627–9.958 ms |
| 512x384 cutaway | Existing Metal renderer | 4.970 ms | 4.680 ms | 4.451–6.763 ms |

These observations establish a working two-consumer path on a simple public
hierarchy. They do not establish improvement over the core CPU spatial lookup,
a new 3D rendering speedup, complex-model performance, or GPU utilization.
The cutaway follows exterior rendering on the same context, so its first frame
does not include cold renderer setup. This sample is not the private broad
ray-based slice whose one-second target remains unmet.

Reproduce with:

```sh
make test-compute USE_METAL=0 RELEASE=1
make test-compute USE_METAL=1 RELEASE=1
ALEA_NUM_THREADS=1 make benchmark-gpu-geometry USE_METAL=1 RELEASE=1
```

## 2D plotter integration (`3dda308`)

The follow-up connects this classifier to `mc_plotter --backend=metal`, for
static PNG/BMP output and the SDL viewer. CPU remains the default. The tool owns
one snapshot across frames and classifies pixel centers in batches of 4,096,
checking worker cancellation between batches and before complete-frame publication.
Cancellation does not set a process-wide interrupt, and the same snapshot can be
reused immediately. No partial slice is published or saved.

Sampled terminal void, gaps and overlaps remain distinct: white, orange and
magenta respectively. Cell/material colors, contours, labels and ticks use the
returned grid. Unsupported CPU diagnostics, ray sampling, surface labels, batch
mode and lattices fail explicitly. The Metal viewer hides CPU-only controls.
This is approximate FP32 visualization, with the classifier's existing sampling
and performance limits; the query kernels are unchanged.

Validation on the Apple M1 passes the existing CPU plotter regression suite,
native Metal image tests, and the SDL dummy-driver viewer tests. Public fixtures
exercise all three axes and image orientation, translated fills, styles,
void/gap/overlap colors, option/environment/lattice rejection without output,
zoom/reset, axis/slice movement, save and quit. A direct plotter test cancels
before the first submission, between submissions and before publication, checks
that no image escapes, and compares the next frame with the original baseline.
The existing CPU and Metal 3D viewer checks remain passing. The CPU-only plotter
also builds and rejects a Metal request explicitly.

The new device check is `python3 tools/check_metal_plotter.py`; the shared
viewer/cancellation check is `USE_METAL=1 python3 tools/check_sdl_viewers.py`
(with `SDL_CONFIG` pointing to the local SDL installation if needed). No private
geometry was used. This increment introduces no new performance claim.

## Rectangular-lattice classification (`6570a4d`)

Native Metal FP32 and the CPU FP64 enumeration reference now classify finite
and repeating rectangular lattices procedurally. The point traversal selects an
element before descending, rather than testing the unit-element CSG as if it
were the full lattice envelope. Ancestor CSG retains its clipping role. Nested
lattices, transformed ancestor fills, outer universes, zero-element coordinates,
terminal void, empty fills and overlapping occurrences preserve the terminal
ownership contract. Unsupported lattice types, lattice-local fill transforms,
malformed dimensions, pitches and fill tables fail Metal validation before CSG
packing or spatial specialization.

Point-specific lattice metadata uses the same dimensions and pitch as rendering,
but its fill and outer links address complete universe cell lists, including
void. Renderer BVHs intentionally omit void and cannot provide these links.
The query uploads are lazy, independent of renderer lattice uploads, and retain
shared transform ownership in either initialization order. No new public API,
core point-lookup change, CUDA implementation, or render-kernel change is added.

Tiled axes select with floor; shared element faces belong to the upper index.
Finite singleton Y/Z windows include their upper face as in core point queries.
Untiled Z remains unshifted in center coordinates. Nonpositive fill entries are
gaps. Metal reports an error for relative/logical indices beyond the exactly
adjacent FP32 integer range (magnitude 16,777,215); CPU uses checked 32-bit indices.
Ordinary boundary rounding remains approximate and is not cured by these guards.

Validation passes 15 CPU compute tests, 62 native Metal compute tests, the builder
test and three upload-resource tests on Apple M1. Analytic public fixtures cover
nested/transformed/finite/repeating lattices, outer void universes, zero-element
coordinates, gaps/overlaps, faces on both sides of a boundary, finite singleton
layers, unbounded Z, axial repetition, owned snapshots, batches crossing 4,096
points, interruption/reuse and distant-index rejection. Existing rectangular
3D render comparisons also classify on the same snapshot in both initialization
orders. Injected allocation failures check retries and release of query buffers
before and after render preparation. Plotter image and SDL navigation tests
include public repeated and transformed nested lattices.

A bounded private E-lite check completes with 256/256 ownership results agreeing
between Metal FP32 and the CPU FP64 enumeration reference. It samples both void
and material and includes overlaps. Metal snapshot preparation takes 9.828 s;
first/warm queries take 1.465/1.356 s. CPU snapshot preparation takes 0.312 s and
its query 0.678 s. These are one process, one first and one warm Metal query on a
16x16 sample; they are not a full-resolution benchmark or independent proof of
private-model correctness. The CPU result is faster for this sample. The lattice
blocker is removed, but large interactive slices remain too expensive without
point-specific spatial pruning. Raw inputs, coordinates, images, logs and the
private probe remain ignored under `build/native-gpu/` or outside the repository.
The rebuilt `mc_plotter` also completes a 16x16 E-lite PNG slice.
Only aggregate observations are recorded here.

## Certified point-program selection

The Metal classifier now chooses an existing spatial CSG program when the
FP32 point lies inside its domain and a point-specific certificate has proved
all original primitive evaluations finite throughout that domain. It retains
the original expression outside the domain or when certification is unavailable.
All candidate cells remain eligible, including void and overlapping occurrences;
no visibility BVH or child-content bound is used to reject a point owner.
Lattice element selection and ancestor transforms retain their existing rules.

The finite certificate matters because a Boolean simplification can otherwise
hide a nonfinite operand that the original point evaluator would report as an
error. It conservatively checks even operands that might short-circuit. Interval
arithmetic encloses rounded FP32 intermediate values, and proof work is capped
at 32 million source instructions per snapshot; remaining cells simply retain
their original programs. No new CSG bytecode or spatial tree is built specifically
for classification. The per-cell point link grows from 8 to 12 bytes. Transform
and spatial-program GPU uploads are shared with the renderer, preserving retry
and ownership in either initialization order. CPU queries are unchanged.

An initial point-only candidate BVH was tested and discarded: broad conservative
bounds increased the private sample’s warm time from 1.405 to 2.225 seconds.
The final change retains the cell scan and reduces expression evaluation instead.
This avoids enabling a traversal strategy that regressed the representative case.
The prototype and raw diagnostics are retained only under ignored `build/`.

Validation passes 15 CPU tests, 64 Metal compute tests, two builder tests and
three resource tests. New comparisons exercise spatial selection on/off, voids,
overlaps, transformed lattices, disconnected parents, bad source-bound hints,
region faces and adjacent FP32 values, out-of-domain fallback, arithmetic errors,
interruption and retry. Public image and SDL tests cover the actual plotter and
both 2D/3D resource consumers.

Apple M1 release build, one process, identical private 16x16 sample, two snapshots,
one first query each followed by five timed warm queries per mode. Execution
order alternates on each repetition; GPU calls run serially. API timing includes
conversion, submissions, synchronization and copies, but excludes preparation.
`ALEA_NUM_THREADS=1`. Both modes agree with the CPU FP64 enumeration reference
on all 256 ownership results; this is agreement on a sample, not independent
proof of complete private-model correctness.

| Point mode | Snapshot preparation | First query | Warm median | Warm range |
| --- | ---: | ---: | ---: | ---: |
| Original scan | 9.570 s | 1.570 s | 1.456 s | 1.419–1.463 s |
| Certified spatial programs | 10.723 s | 1.357 s | 1.175 s | 1.170–1.216 s |

Warm query time decreases about 19.3% (1.24x throughput); preparation increases
in this run. CPU enumeration takes 0.594 s for the same sample and remains faster.
This is not the 512x512 one-second target. The public 64-occurrence analytic
benchmark at 512x512 passes its independent equations in both modes: Metal warm
medians are 32.898 ms off and 33.911 ms on, with overlapping ranges. No public
simple-scene or 3D rendering speedup is claimed.

The optimization is enabled by default. Reproduce the original point evaluator
with `ALEA_METAL_POINT_SPATIAL=0`; setting this does not suppress the renderer’s
spatial preprocessing. Private geometry, coordinates, images, probe sources and
raw logs remain ignored or outside Git. Only aggregate observations are recorded.

## Point-query cost attribution

Opt-in `ALEA_METAL_POINT_PROFILE` now reports aggregate versioned JSONL for
complete successful queries. The Metal shader function constant removes counters
from normal pipelines. The profiler covers flat and hierarchical classification,
counts actual bytecode/primitive execution, separates original and spatial
programs, records original-program selection reasons and CSG rejection work,
and reports bridge wall time separately from available GPU timestamps. Counters
reset per query and aggregate across the existing bounded submissions. Report
failures surface as errors without invalidating the context. See
[the profiler contract](GPU_COMPUTE.md#point-query-profiling).

Validation on Apple M1: 15 CPU, 67 Metal, 2 builder and 3 resource tests pass.
New checks compare profiled and normal ownership, cover spatial and lattice
accounting, multi-submission aggregation, interruption/nonfinite errors, empty
calls, repeated use and report-write failure/retry. Public 512x512 point and
512x384 exterior/cutaway analytic oracles pass with profiling off and on. The
public warm slice medians are 30.658 ms off and 37.853 ms on (about 23.5% overhead);
these are three warm calls per mode, measured sequentially. Profiling is a
diagnostic option, not a performance improvement. No new 3D speedup is claimed.

A private comparison keeps spatial selection enabled in both contexts,
alternates profiling off/on across one initial and five warm queries per mode,
and uses the same bounded 256-point sample as before. Only aggregate observations
are retained here:

| Measurement | Normal | Profiled |
| --- | ---: | ---: |
| Preparation | 9.541 s | 9.585 s |
| First query | 1.386 s | 1.365 s |
| Warm median | 1.226 s | 1.234 s |
| Warm range | 1.207–1.369 s | 1.229–1.300 s |
| Differences from CPU FP64 enumeration | 0 / 256 | 0 / 256 |

The roughly 0.7% median overhead is smaller than the observed timing spread and
is not an overhead guarantee. CPU enumeration takes 0.646 s. Profiled GPU time
has a 1.230 s warm median; bridge work outside the GPU is small in this sample.
All six reports parse as JSON, have identical work counts, and report no counter
saturation:

| Work over 256 points | Count |
| --- | ---: |
| Cell visits | 1,705,697 |
| CSG rejections | 1,705,144 |
| Original / spatial program calls | 942,296 / 763,401 |
| No program / uncertified | 404,961 / 35,328 |
| Outside certified domain / unchanged cache region | 474,849 / 27,158 |
| Original / spatial instructions | 40,534,967 / 16,521,851 |
| Rejected / accepted instructions | 56,962,030 / 94,788 |
| Original / spatial primitive evaluations | 7,186,991 / 2,926,989 |
| Universe descents / maximum depth | 269 / 2 |
| Lattice visits | 0 |

**99.83% of executed CSG instructions reject cells.** This is an instruction
fraction, not a measured time attribution. It supports targeting rejection work
before arithmetic tuning. The sample visits no lattice cells, so it does not
establish lattice traversal cost; public lattice fixtures validate the counters
and semantics separately. It also cannot predict full-resolution slice latency.

The user-facing Metal plotter is rebuilt. CPU plotter regressions and Metal image
checks pass, including 52 valid aggregate reports with profiling enabled. SDL
dummy-driver checks pass for 2D/3D navigation, publication, cancellation and reuse.

## Opt-in tile candidate filtering

The candidate prototype is retained behind `ALEA_METAL_POINT_CANDIDATES=1`
(default off). It proves root-cell rejection over the bounding box of each
converted FP32 input batch, reuses primitive interval proofs within that batch,
and supplies an ordered candidate list to the native Metal classifier. Original
CSG operands must all remain finite; unknown proofs, boundary contacts and
lattice containers stay in the list. All-owner semantics, void and ancestor
clipping are preserved. Child-universe scans are unchanged.

The host reads the existing shared Metal geometry buffers. Scratch is one byte
per packed primitive and four bytes per packed cell; lists are rebuilt per
batch rather than cached by view. A 32-million-source-instruction budget bounds
proof work, with unexamined cells retained. Cancellation checkpoints occur
between cells and every 1,024 instructions. `mc_plotter` stages 64x64 tiles when
the option is enabled and scatters results to the original pixels. The generic
point API retains input order and caps candidate batches at 4,096 points.

Validation: 15 CPU, 69 Metal, 4 builder and 4 resource tests pass. Candidate tests
cover primitive/CSG semantics, transformed occurrences, void/overlap ownership,
rectangular/repeating/nested lattices, exact boundaries, empty scenes, all-rejected
lists, nonfinite operands, proof-budget fallback, interruption/retry, allocation
failure and rebuilding lists between batches. The public Metal image checker
compares scanning and candidates, including partial edge tiles. CPU plotter
regressions and SDL 2D/3D navigation, cancellation, publication and reuse checks
pass with candidates enabled. The user-facing Metal plotter is rebuilt. These checks
retain the original FP32 accuracy limits; this is not transport-grade geometry.

### Complete public plotter comparison

`python3 benchmarks/gpu_plotter.py` generates the public 64-occurrence sphere/void
fixture, checks complete image equality and measures the actual CPU grid/ray
plotters against both Metal modes. Four fresh-process trials per mode, alternating
order, Apple M1, four CPU threads, 512x512, profiling disabled:

| Plotter path | Query median | Query range | Process wall median |
| --- | ---: | ---: | ---: |
| CPU grid | 6.450 ms | 6.300–6.500 ms | 16.2 ms |
| CPU ray raster | 2.900 ms | 2.300–3.700 ms | 12.2 ms |
| Metal scanner | 35.950 ms | 35.600–37.900 ms | 70.8 ms |
| Metal tile candidates | 24.050 ms | 21.600–24.700 ms | 56.0 ms |

All images agree. Metal query time improves about 33%, including candidate proof,
staging and all submissions. The optimized CPU plotters remain faster. Process
wall includes loading, acceleration preparation, drawing, output and shutdown;
query timing does not. This comparison replaces neither the analytic compute
oracle nor the separate 3D regression workload.

### Bounded private frame comparison

The private check was scaled from 256 sample points to a 32x32 view and then
128x128; no 512x512 private run was attempted. At 32x32 the two Metal images match,
and both classifiers match CPU FP64 enumeration on all 1,024 points. The CPU grid
image differs because 85 samples have multiple terminal owners: Metal marks an
overlap, whereas CPU grid plotting selects an owner. All 939 unique samples have
the same CPU grid IDs/materials as the enumeration reference. There are no gaps.
This is a semantic difference between plotting paths, not a candidate regression.

At 128x128, an actual CLI run measures CPU grid at 0.385 s, Metal scanning at
2.777 s and Metal candidates at 2.283 s; fresh-process wall times are respectively
15.823, 22.453 and 21.776 s. The two Metal images are identical. Loading/preparation
dominate fresh-process wall time.

A separate in-process comparison retains each snapshot, alternates modes across
one initial and five warm frames, and includes coordinate generation, staging,
candidate proofs, all four submissions and scattering. Scanner input uses the
existing row order; candidate input uses four 64x64 tiles. Profiling is disabled.

| Measurement | Metal scanner | Metal tile candidates |
| --- | ---: | ---: |
| Snapshot preparation | 9.631 s | 9.547 s |
| First frame | 2.870 s | 2.233 s |
| Warm median | 2.722 s | 2.198 s |
| Warm range | 2.658–2.928 s | 2.152–2.339 s |
| Differences from CPU FP64 enumeration | 0 / 16,384 | 0 / 16,384 |

The warm frame improves about 19.3%. Optimized CPU grid sampling has a 0.523 s
warm median (0.443–0.677 s), but implements selected-owner semantics. Independent
CPU all-owner compute enumeration takes 12.391 s. These CPU paths must not be
conflated when making speed claims. The private geometry, coordinates, images,
probe sources and raw logs remain ignored/outside Git.

A separate instrumented frame comparison attributes the reduction:

| Work per 128x128 frame | Scanner | Tile candidates |
| --- | ---: | ---: |
| GPU cell visits | 107,584,636 | 40,078,460 |
| GPU executed CSG instructions | 3,671,369,563 | 1,678,030,131 |
| Host root candidates examined / rejected | 0 / 0 | 22,648 / 16,481 |
| Host source instructions / primitive proofs | 0 / 0 | 7,848,296 / 131,072 |
| Candidate proof time | 0 | 0.092 s |
| Largest GPU submission | 1.079 s | 0.723 s |

All counters are unsaturated. Cell visits fall about 63% and CSG instruction
execution about 54%, but elapsed-time improvement is smaller; instruction counts
are not time fractions. The remaining GPU submissions still prevent a tight
cancellation-latency guarantee.

A smaller diagnostic split also demonstrates why coherent tiles alone are not
enough: splitting the old 256-point sample into four 64-point calls takes 3.786 s
with scanning and 2.220 s with candidates, versus roughly 1.2 s for the original
single batch. This is not a reason to subdivide already small queries. The
plotter's 64x64 tiles retain up to 4,096 points per submission.

## Next increment

Keep the option experimental and the optimized CPU plotter as the practical
performance baseline. Investigate expensive surviving candidates and child
universe scans before extending the proof machinery; use profiling only for
attribution and disabled profiling for timing. Any extension must include proof
construction, memory, complete frame latency and cancellation behavior in its
acceptance criteria. The private 512x512 one-second target is still unmet.

Maintain 2D and 3D as regression consumers. Establish separate 3D attribution on
exterior/cutaway/close-up views before changing ray traversal; the point candidate
list is not a ray visibility structure.

Private models, coordinates, extracted geometry, raw hashes, and diagnostics
remain ignored or outside Git. Tracked fixtures use public synthetic geometry;
private validation contributes only aggregate observations. Local logs and images
are under ignored `build/native-gpu/`. This experimental report is not added to the public
documentation-site allowlist.
