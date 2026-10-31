# Prepared-bounds navigation experiment

Status: benchmark-only native Metal rasterization, 2026-10-30. **Rejected as a
standalone model preview.** Bounds draw quickly, but the private model produces an
unreadable mass of overlapping edges. No viewer or public API was changed.

## Representation

The prepared render scene already contains root-universe BVH leaf bounds in world
coordinates. `benchmarks/metal_bounds_preview.m` extracts finite, nonempty bounds
from that root tree. It does not descend into child universes, transform child
bounds, expand lattices or generate geometry from CSG. Filled containers therefore
appear as containing boxes; their internal surfaces are absent.

`benchmarks/bounds_preview.metal` uses a native Metal vertex/fragment pipeline:
24 line vertices per box for the navigation view, and 36 triangle vertices for a
separate filled-coverage diagnostic. Corners are generated in the vertex shader;
there is no CPU vertex mesh. Hidden-line removal and depth-based occlusion are
intentionally absent: this is an x-ray **bounds** view. The filled image is not a
geometry surface, a usable detailed image or a query result.

The benchmark limits extraction to 65,536 boxes and fails instead of truncating.
Unbounded, nonfinite and degenerate bounds are counted and omitted. These omissions,
root-only visibility and use of tightened render bounds mean this display must
not be advertised as a universal conservative envelope. The public fixture's
coverage is independently measured; no private conservativeness proof is claimed.
User clip planes, material filtering and source revisions are not implemented.
The near clip is 0.001 scene units; orthographic far clip is 1e8. Those experiment
constants would need scene-aware treatment in a product.

## Reproduction and validation

```sh
make -j4 RELEASE=1 USE_METAL=1 \
  BUILD_DIR=build/native-gpu/release BIN_DIR=build/native-gpu/bin \
  build/native-gpu/bin/compute/metal/benchmark_metal_bounds_preview

python3 benchmarks/native_gpu_compare.py \
  --config benchmarks/metal_bounds_preview.example.json \
  --output build/native-gpu/bounds-preview/example-runs --blocks 2 --timeout 60
```

Use a new output directory. The example includes ordinary and Metal API-validation
runs. The public images under ignored `build/native-gpu/bounds-preview/` show
wire bounds, filled bounds and fresh analytic geometry, left to right. The repeated
ordinary run overwrites these public comparison images.

Two ordinary and two API-validation processes passed. Independent CPU slab
intersections check every pixel of the synthetic filled-box masks in perspective
and orthographic modes, including behind-camera boxes and near-plane crossing.
Other checks cover empty/reset behavior, capacity rejection, root-only extraction
and invalid-bounds omission. The capacity test intentionally emits a rejection
message. Metal API validation reported no errors. Fresh public source renders pass
the independent sphere intersection oracle, and no analytic hit is lost by the
public filled-bounds masks in the three evaluated views.

The allocation guard limits dimensions to 1024x1024, widths to multiples of 64 for
the readback row alignment, and preview scratch to 64 MiB. At the measured 512x384
size, the conservative preview-storage bound is 4.5 MiB: the maximum box buffer,
one render texture, a readback buffer and two host RGBA images. Optional oracle
framebuffers and image-encoding scratch are separate.

## Public measurements

Apple M1, release build, 512x384, 64 translated spheres. Each process draws one
initial wire frame plus ten warm frames for each view. The following ranges are
the two ordinary process medians; API-validation timings are excluded.

| View | Warm wire draw/readback ms | Analytic hits retained by filled bounds | Filled pixels that actually hit geometry |
| --- | ---: | ---: | ---: |
| Exterior | 0.411–0.435 | 100% | 67.78% |
| 20-degree orbit | 0.335–0.461 | 100% | 67.34% |
| Orthographic | 0.277–0.382 | 100% | 69.05% |

Drawing time includes encoding, GPU completion, texture-to-buffer copy and host
RGBA copy. Median GPU time was 0.068–0.110 ms. Scene construction, extraction,
pipeline creation, analytic reference renders and image encoding are outside warm
drawing time. These are isolated measurements, not interactive viewer latency.

Full native scene preparation took 26.2–54.1 ms, raster pipeline/buffer setup
0.218–0.400 ms, and extraction/shared-memory upload 0.008–0.014 ms. The first
readback after an already prepared scene took 2.173–4.158 ms, including raster
setup and extraction. A reused OS/driver shader cache may affect pipeline timing;
these are first calls within each process, not guaranteed cold-machine timings.

All 64 root leaves produced boxes, using 1,536 bytes of actual box data. Visual
inspection shows recognizable object layout, but filled boxes merge gaps and
replace curved surfaces with block silhouettes. About 31–33% of the filled pixels
are empty space according to the analytic renderer. A fast rasterizer does not
make that representation faithful.

## Private bounded evaluation

One serialized process reused the local private model and camera specification.
No analytic target rays were traced; this run measures bounds startup, drawing and
clutter only. Raw logs, exact cameras, source probe, images and model paths remain
ignored. Only aggregate findings are recorded here.

- Model load: 7.866 s; full native scene preparation: 10.295 s.
- 5,590 root leaves; 5,583 finite boxes; seven omitted leaves.
- Actual box data: 133,992 bytes; extraction/upload: 0.047 ms.
- Raster setup: 1.112 ms; prepared-scene to first readback: 11.414 ms.
- Warm wire medians: exterior 0.933 ms, orbit 0.821 ms, orthographic 0.551 ms.
- Wire pixels occupied 72.9–75.6% of the display; filled bounds occupied 100%
  in every view. Those percentages measure clutter, not geometry accuracy.
- Peak process memory: approximately 2.45 GiB; native model buffers: approximately
  365.2 MiB. Total process time: 18.261 s, within the experiment's 90 s guard.

The existing 15 s preparation, 4 GiB process and 1 GiB model-buffer checks passed.
No extra private runs were needed after the usefulness failure. Visual inspection
of the exterior image confirms dense overlapping lines that obscure structure.
The blank third panel in private artifacts means **no analytic reference was
rendered**, not an empty model.

The first bounds image still depends on full native preparation: approximately
18.17 s from starting model load in this observation. Its sub-millisecond warm
draw time therefore does not meet a sub-200 ms initial-preview target. Drawing
while analytic commands occupy the GPU, SDL presentation, cancellation and stale
camera handling have not been tested by this experiment.

## Decision and continuation

Retain the bounded benchmark and reject this representation as the standalone
motion-preview path. It may eventually be useful as an explicit selection/bounds
overlay, but that is a different task. Do not enable it automatically or mistake
filled-bounds coverage for visible-surface quality.

Together with the [reprojection experiment](NATIVE_GPU_REPROJECTION.md), this
isolates the remaining issue: drawing cached/display geometry is cheap; obtaining
useful visible geometry before a dense analytic frame is available is not solved.
Neither adding more display kernels nor flattening the CPU geometry machinery is
justified by these results.

The next priority is the work-weighted native CSG preparation audit already
identified in the [execution audit](NATIVE_GPU_EXECUTION_AUDIT.md) and
[scaling evaluation](NATIVE_GPU_SCALING.md): rank expensive unspecialized programs,
determine why they retain rejecting instructions, then select one bounded,
conservative reduction. This benefits both parallel 2D queries and parallel 3D
rays, and any later sampled preview representation. Measure removed executed work
and end-to-end latency; preserve ownership, hierarchy, overlaps and boundary rules.

A surface proxy remains a separate possible representation, with an explicit
construction budget and missing/thin-feature policy. The existing `alea_mesh`
module samples CSG onto a structured grid; it is not an already available cheap
CSG-to-triangle converter. Reusing it would still require measuring sampling cost,
resolution and visibility quality. Do not port that subsystem to the GPU merely
to obtain a preview. CUDA and the existing analytic accuracy issues remain open.
