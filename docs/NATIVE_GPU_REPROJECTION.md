# Cached-depth motion-preview experiment

Status: benchmark-only native Metal prototype, 2026-10-30. **Useful for small
camera movements after a dense frame exists; insufficient as the standalone
motion-preview strategy.** No viewer, public API, geometry-query or rendering
default changed. CUDA implementation and validation remain hardware-gated.

## Scope and implementation

`benchmarks/metal_reprojection.m` caches world positions and linear RGB from a
completed analytic AA1 frame. Four kernels in `benchmarks/reprojection.metal`
clear scratch, project points, choose the nearest source sample, and resolve color.
Equal projected depths choose the lowest source index deterministically. Each
sample contributes to at most one target pixel. There is no surface triangulation,
splatting, interpolation or hole filling.

Only display color and an unknown mask are candidate preview outputs. Internal
projected depths/source indices are inspected by the benchmark to measure quality;
they are not geometry answers. Uncovered pixels, including background, have alpha
zero and a visible checkerboard. The cache excludes misses and the fixture's void
material. That material convention is specific to this experiment, not a general
library visibility contract.

Cache construction uses CPU ray reconstruction and writes shared Metal memory.
There is no separate discrete-GPU upload on this Apple M1. Projection includes
command encoding, four dispatches, completion wait and RGBA readback. Source
rendering, source-cache construction, comparison oracles and PNG writing are
outside projection timing. These are isolated benchmark timings, not
request-to-visible SDL latency or concurrent-render measurements.

The prototype accepts bounded dimensions/counts before dispatch, allocates scratch
once, and skips empty input dispatches. Its conservative preview-storage bound is
4.875 MiB at 256x192, including point staging/readback allowances; the allocator
rejects capacities above 64 MiB. Analytic model, source/reference framebuffers and
driver overhead are separate. The prototype has no source-lifetime or scene-revision
integration: those remain requirements before a viewer can consume it.

## Reproduction

Public synthetic geometry only: 64 translated spheres with an independent
intersection oracle. Run from the repository root on a Metal-capable Mac:

```sh
make -j4 RELEASE=1 USE_METAL=1 \
  BUILD_DIR=build/native-gpu/release BIN_DIR=build/native-gpu/bin \
  build/native-gpu/bin/compute/metal/benchmark_metal_reprojection

python3 benchmarks/native_gpu_compare.py \
  --config benchmarks/metal_reprojection.example.json \
  --output build/native-gpu/reprojection/example-runs --blocks 2 --timeout 60
```

The output directory must be new. The public comparison image is written under
ignored `build/native-gpu/reprojection/`: source, 2-degree orbit reprojection,
fresh analytic target, left to right. It is overwritten on each run. Logs and
machine-specific artifacts remain ignored. The runner serializes GPU processes
and applies timeouts; it does not establish statistical performance guarantees.

## Measurements

Apple M1, release build, 256x192 source and target, AA1, shadows disabled.
Two ordinary processes, each with one initial projection plus ten warm projections
per camera; two additional API-validation processes passed. The table gives the
range of the two ordinary warm medians. Validation timings are excluded.

| Camera change | Warm projection ms | Target-hit coverage | Missing hit pixels | Ghost pixels |
| --- | ---: | ---: | ---: | ---: |
| Perspective identity | 0.316–0.328 | 100% | 0 | 0 |
| 2-degree orbit | 0.305–0.424 | 94.93% | 331 | 227 |
| Pan by 0.25 scene units | 0.335–0.444 | 95.59% | 285 | 263 |
| Move 5% closer | 0.294–0.372 | 86.48% | 982 | 220 |
| 20-degree orbit | 0.323–0.327 | 82.70% | 1,113 | 279 |
| Orthographic identity | 0.295–0.350 | 100% | 0 | 0 |
| Orthographic pan | 0.286–0.381 | 94.38% | 408 | 401 |
| 16x12 source, 256x192 target, same camera | 0.262–0.271 | 0.40% | 6,484 | 4 |

Coverage means the fraction of **fresh analytic hit pixels** receiving a projected
sample. It does not imply correct depth or color. A ghost is a projected sample
where the fresh ray misses. Pixel-center sampling alone can produce boundary
ghosts, even without camera motion when the source and target resolutions differ.
Unknown background is excluded from hit coverage.

For the 2-degree orbit, 236 covered pixels differ in forward depth by more than
one target pixel's world-space height at that depth. Mean absolute RGB error on
covered hits is 0.0241 on the raw 0–1 scale. For the 20-degree orbit it is 0.0619.
These are diagnostic tolerances, not geometry-accuracy acceptance criteria.
Cached shading is view-dependent and does not become correct merely because a
surface sample projects to the right pixel.

Normal-run median GPU execution was 0.021–0.089 ms; maximum observed warm wall
time was 2.559 ms. Pipeline/buffer setup took 27.9–50.2 ms. CPU cache construction
and shared-buffer writes took 0.001–0.241 ms, including the sparse case. These
small samples establish feasibility, not a p95 responsiveness guarantee.

## Validation and limitations

All four final processes passed nearest-depth and deterministic tie selection,
invalid/behind/offscreen rejection, reset, empty input, capacity rejection and
translated/rotated camera tests in both projection modes. Synthetic projection
tests use non-multiple dispatch dimensions and independently constructed pixel
centers. Fresh analytic source and target renders pass the sphere oracle.
Perspective and orthographic same-resolution identity views retain every hit,
with no ghosts, no depth-tolerance failures and identical cached RGB. Metal API
validation reported no errors.

Visual inspection of the public comparison shows recognizable geometry during a
small orbit, with holes and shifted silhouettes visible against the unknown mask.
The prototype deliberately exposes those limitations.

The existing ignored private completed 16x12 frame has 192 valid samples. Even if
every sample lands on a different pixel, this implementation can populate at most
192 / 49,152 = **0.390625% of a 256x192 display**. This is an analytical upper bound
from the retained source data, not a private GPU timing or target-quality result.
No new private model load, target render or private image was needed. The public
sparse-source experiment independently demonstrates the same practical problem.

Point projection cannot reveal newly visible surfaces, create a first frame,
refresh view-dependent shading or recover detail missing from the source. It also
has not been tested while analytic work occupies the GPU. A sub-millisecond
isolated command does not guarantee that it can preempt a long analytic command.

## Decision and next boundary

Retain this small benchmark; do not enable it as the viewer's motion-preview path.
The latency gate is promising, but the representation fails sparse-source coverage
and cannot meet the private first-image requirement. The [scaling evaluation](NATIVE_GPU_SCALING.md)
remains the evidence for slow analytic preview generation.

Before viewer integration, evaluate a bounded preview representation that can
produce a useful first image independently of a completed dense analytic frame.
Prefer reusing prepared geometry to duplicating CPU services. Its experiment must
measure preparation time, memory, visible structure and disocclusion quality as
well as drawing time. Simply enlarging projected points hides unknown space and
is not evidence that the representation is sufficient. Any bounds-only proxy must
be visibly identified as bounds, rather than presented as a geometry surface.

If a suitable source/proxy becomes available, this projection kernel can be a
small-motion display bridge. Integration then needs immutable completed-source
snapshots, scene/clip/appearance invalidation, AA sample-position metadata,
camera-revision rejection, bounded scheduling alongside analytic work and explicit
unknown display. Never feed a reprojected frame back into the cache. Settled detail
and exports must still use fresh analytic rays; physics queries must remain exact
according to their existing backend contracts. Measure complete UI latency and
cancellation before declaring the interaction target met.

This experiment does not reduce analytic CSG work or settle existing accuracy
issues. Native point continuations and work-weighted preparation opportunities
remain separate implementation tasks; no additional CPU subsystem was ported.

### Prepared-bounds follow-up

The [prepared-bounds raster experiment](NATIVE_GPU_BOUNDS_PREVIEW.md) now evaluates
an independent first-image source. Drawing is fast, but overlapping root bounds
obscure the private model and still require full scene preparation. It was rejected
as a standalone preview. Both experiments point back to the cost of obtaining
useful geometry, rather than display-kernel throughput.
