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

## Next increment

Reduce point-query candidate visits using tighter conservative proofs that can
handle broad and uncertain cells without adding a traversal penalty. Establish fresh 3D
attribution on representative exterior/cutaway/close-up views before changing
its traversal. Both consumers remain acceptance workloads throughout.

Private models, coordinates, extracted geometry, raw hashes, and diagnostics
remain ignored or outside Git. Tracked fixtures use public synthetic geometry;
private validation contributes only aggregate observations. Local logs and images
are under ignored `build/native-gpu/`. This experimental report is not added to the public
documentation-site allowlist.
