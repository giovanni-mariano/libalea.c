<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO
SPDX-License-Identifier: MPL-2.0
-->

# Experimental GPU rendering and batch compute

`alea_compute.h` provides an optional CPU/Metal geometry snapshot and rendering
API. CPU evaluation uses double precision; Metal uses approximate FP32 analytic
CSG. Existing core APIs retain their CPU implementation. CUDA and HIP are not
implemented. GPU transport is outside this module.

The Metal backend supports useful interactive rendering on smaller scenes.
**Full ITER E-lite is still too expensive for smooth navigation.** The recorded
16x12 default cutaway takes about 20 seconds. Responsive cancellation does not
mean that a useful new image is available at interactive speed.

## Build and use

On macOS, install Xcode with the Metal toolchain. SDL2 is required only for the
interactive viewer; see [Installation](INSTALLATION.md) for the local fallback.

```sh
make compute USE_METAL=0 RELEASE=1
make compute USE_METAL=1 RELEASE=1
make -C tools ../bin/render3d USE_METAL=1 USE_SDL=1 RELEASE=1
bin/render3d examples/data/doc_shield.i --backend metal --interactive \
    --width 800 --height 600 --eye 20,-30,20 --target 0,0,0 --clip-y=0
```

Archives are `bin/compute/cpu/libalea_compute.a` and
`bin/compute/metal/libalea_compute.a`. The module is separate from the default
and full archives. Link the selected compute archive before `libalea.a`; Metal
also needs Foundation and Metal frameworks. `make install-compute` installs the
selected archive and header. `--backend cpu` is the viewer default. A requested
backend never silently falls back to another backend or precision.

Create an owned snapshot once and reuse it across frames:

```c
alea_compute_t *gpu = alea_compute_create(sys, ALEA_COMPUTE_METAL,
                                        ALEA_COMPUTE_FP32);
if (!gpu) { /* inspect alea_error() */ }
/* Prepare camera/framebuffer with the existing render APIs. */
int rc = alea_compute_render(gpu, &config, &camera, framebuffer);
alea_compute_destroy(gpu);
```

The source may change or be freed after snapshot creation. Recreate the context
to incorporate geometry changes. Serialize calls on a context; independent
contexts can be used independently. Failed calls return an error rather than a
partial success. CPU contexts use `ALEA_COMPUTE_FP64`.

## Geometry and precision

Metal rendering supports planes, spheres, axis-aligned cylinders and cones,
general quadrics, RPP boxes, RCC cylinders, and circular/elliptical tori. It
supports transformed nested universes and finite/repeating rectangular lattices.
Universe geometry is shared across placements and lattice traversal is procedural.
Parent CSG intervals clip children, retaining parent-created surfaces and normals.

Supported image features are solid/depth/cell-ID/material-ID output, material or
cell colors, perspective/orthographic cameras, one clipping plane, and AA1–AA4.
CPU postprocessing supplies edges and image encoding. Shadows, X-ray, selection
filters, density/universe colors, and adaptive AA are rejected. Hexagonal
lattices and lattice-local fill transforms are unsupported; ordinary transformed
fills may contain a rectangular lattice. Overlapping solid cells have no defined
ownership resolution policy.

FP32 can change ownership near boundaries and erase small features. This is an
approximate visualization backend, not a transport or verified geometry engine.
Conservative specialization preserves the packed FP32 problem; it does not
restore the original double-precision coefficients. ITER accuracy triage found
that the CPU viewer misses earlier material on several inspected rays, and
identified a quadric-packing defect that created false water/steel overlaps. The packer now avoids distant origins in
nearly singular directions; both an isolated regression and the full ITER render
confirm that the disputed pixel now shows steel. The 32x24 comparison improves
from 14 to 10 cell-ID differences. Bounded interval checks attribute six to later
CPU viewer hits and four to same-material clipped intervals roughly 0.0004 units
wide, consistent with FP32 boundary sensitivity. In 17 disputed earlier-material
cases, the largest measured local Metal entry error is 0.0216 model units;
this is not a global error bound or an exhaustive image validation. Agreement
between two Metal configurations is not an independent accuracy check.

Source expression depth is limited to 4096. The device uses a 64-value expression
stack, 32 hierarchy levels, and 256 BVH deferred entries. Traversal/work limits,
nonfinite arithmetic, unsupported geometry and device allocation failures are
reported explicitly. Controlled tests include analytic geometry and independent
flattened references; finite test coverage does not prove all-model accuracy.

## Interactive and progressive rendering

The viewer starts with a coarse AA1 preview, then refines after the preview
completes and 150 ms passes without input. The initial longest edge is 96 pixels,
or 16 for models with more than 10,000 cells; subsequent previews adapt between
16 and 128. A coarse preview can miss thin features. Full refinement evaluates
every requested pixel. CPU viewing retains its half-resolution preview.

- Left drag/arrows orbit; right or middle drag pans; wheel or +/- zooms.
- C changes material/cell colors; E toggles edges.
- X toggles clipping; [ and ] move the first clip plane.
- R resets; S saves a full-resolution snapshot; H shows help; Q/Escape quits.

One worker owns rendering; the UI consumes immutable copies. Camera/color/clip
changes cancel obsolete work. Stale jobs cannot replace a newer image. Saves
retain their own camera and configuration and finish independently of navigation.

`alea_compute_render_progressive()` requires a Metal context and a callback:

```c
int rc = alea_compute_render_progressive(gpu, &config, &camera,
                                        framebuffer, on_progress, user);
```

The callback receives `(user, framebuffer, completed_pixels, total_pixels)` on
the calling thread, initially at zero and after submissions. The count can stay
unchanged while rays are suspended. Returning nonzero cancels this call. Return
values are 0 for completion, 1 for callback cancellation, and -1 for error.
Initialize unfinished pixels before displaying them and copy data before
publishing it to another thread. A cancelled context can immediately render again.

Up to 4096 rays occupy a bounded continuation pool. Completed slots receive new
pixels between submissions. State includes the hierarchy, expression evaluator,
intervals, closest hit and AA accumulation. Geometry remains uploaded across
camera changes. Continuation storage is capped at 64 MiB; geometry and image
buffers are additional. Work allowance adapts toward 20 ms submissions. This is
a scheduling target, not a latency guarantee.

## Acceleration and diagnostics

Automatic content bounds and adaptive container expressions are enabled by
default. Domains are inferred from finite geometry and used as optimization
hints: geometry outside them remains queryable through the original evaluator.
Only conservative proofs remove predicates or tighten bounds. Poorly conditioned
rays retain original bounds/programs. Adaptive trees are bounded to 10 MiB and
additional bytecode to 64 MiB; source geometry and compiler scratch are extra.

Selective deferred parent clipping searches child bounds first where the cost
estimate justifies it, then intersects material candidates with ancestor truth
intervals. It preserves parent-created boundaries. Expensive child universes and
lattice descendants keep eager clipping. Slot refilling is enabled by default.

Controls are read at context creation:

| Variable | Values | Default |
| --- | --- | --- |
| `ALEA_METAL_SPATIAL` | `0` original; `bounds` content bounds; `1` bounds and adaptive programs | `1` |
| `ALEA_METAL_LAZY` | `0` eager; `1` selective deferred clipping; `all` bypass cost selector for eligible parents | `1` |
| `ALEA_METAL_REFILL` | `0` finish each tile; `1` refill completed slots | `1` |
| `ALEA_METAL_BVH` | `0` median; `sah` surface-area splits; `weighted` instruction-weighted splits | `0` |
| `ALEA_METAL_RESUME_WORK` | Fixed 1–65536 operations per ray/submission | Adaptive |
| `ALEA_METAL_BATCH_PIXELS` | Ordinary-render batch size, 1–4294967295 | Automatic |
| `ALEA_METAL_PROFILE` | Output prefix for instrumented CSVs/heatmap | Disabled |
| `ALEA_BENCH_TRACE` | Set to print scheduling and buffer diagnostics | Disabled |

Weighted/SAH BVHs remain opt-in experiments. Unrestricted lazy clipping increased
ITER work substantially; keep the selective default for normal use. The manual
root-universe grid, candidate-region masks, surface-patch rasterizer and uniform
boundary compiler have been retired. Their environment variables and viewer
options are no longer supported. No primitive-root cache is implemented.

Profiling writes `PREFIX-rays.csv`, `PREFIX-cells.csv`, and
`PREFIX-heatmap.png`, replacing them per completed frame. Counters include BVH
visits, cell queries, instructions, primitive solves, interval/hierarchy steps,
parent checks, cache hits and rejected parent ranges. They saturate at UINT32_MAX.
Per-cell accumulation adds overhead; disable profiling for timing comparisons.
The retired candidate-rejection column is no longer written.

## Measurements and validation

Historical M1 measurements establish the following scope:

| Workload | Recorded result |
| --- | --- |
| 64 spheres, 800x600 | Metal 9.92 ms versus CPU 148.93 ms/frame; no cell-ID differences, some CPU normal discrepancies |
| Sparse expensive parent, 128x96 | Selective deferred clipping: about 24 to 3 ms |
| Mixed ray costs, 256x192 | Slot refilling: about 19 to 14 ms |
| ITER default cutaway, 16x12 | 19.619 s; first pixel 0.603 s; maximum callback gap 29 ms |
| ITER cancellation requested at 2 s | Returned at 2.003 s with 11 of 192 pixels complete |

These are specific workloads, not predictions for other cameras/resolutions.
The recorded ITER preparation is about nine seconds after roughly eight seconds
of model loading. The source contains 450,648 cells, 633,916 surfaces and 198
universes. High GPU utilization and normal image resolution remain separate
from the small 192-ray diagnostic workload.

The 17 October candidate-mask decision used one warm-up and five alternating
warm frames per configuration, keeping normal spatial/lazy/refill defaults:

| ITER 16x12 view | No masks, median | Masks, median |
| --- | --- | --- |
| Cutaway | 19.783 s | 20.544 s |
| Opposite cutaway | 21.808 s | 22.764 s |

All image channels matched. Masks were removed after slowing both tested views;
combined historical experiments do not establish their individual benefit.
These measurements precede the structural cleanup and quadric-conditioning fix;
they do not measure the final renderer’s speed.

The corrected renderer at commit `a3e66ec` was then measured with default
spatial/lazy/refill controls, one warm-up and five serial progressive frames per
view on the M1. The 16x12 cutaway median is 19.511 s (19.460–19.788 s); the
opposite cutaway median is 21.571 s (21.536–22.698 s). Every frame matched its warm-up
bitwise in all output channels. First-pixel medians were 0.581 s and 0.328 s
after scene preparation. Loading took about eight seconds and preparation about
ten seconds per fresh context. These tiny-image timings are a descriptive
baseline, not a speedup attribution or a normal-resolution utilization test.

```sh
make test-compute USE_METAL=0 RELEASE=1
make test-compute USE_METAL=1 RELEASE=1
make bin/compute/metal/benchmark_render_model USE_METAL=1 RELEASE=1
make bin/compute/metal/benchmark_render_progressive USE_METAL=1 RELEASE=1
bin/compute/metal/benchmark_render_model synthetic:mixed p-refill 256 192 5
bin/compute/metal/benchmark_render_model /path/to/E-lite_R250630.mcnp p-balanced 16 12 5 1
bin/compute/metal/benchmark_render_progressive /path/to/E-lite_R250630.mcnp 16 12 2000
```

Model benchmark modes are `auto`, `bounds`, `current`, `lazy`, `sah`, `weighted`,
`refill`, `balanced` and `all`; prefix with `p-` for progressive rendering.
`auto`/`bounds` use an unspecialized reference; other modes use automatic spatial
programs in both configurations. `all` combines weighted BVHs, selective lazy
clipping and refilling. Cameras 1–4 are cutaway, opposite cutaway, exterior and
component close-up. The benchmark reports its binary identity and model metadata.

Regression tests cover primitive geometry, CSG, nested transforms, lattices,
clipping, conservative proofs, original-program fallback, ordinary/progressive
agreement, cancellation/reuse, and upload failure/retry. GPU tests need a real
Metal device. Viewer checks use SDL's dummy driver:

```sh
SDL_CONFIG="$PWD/build/deps/SDL2-2.32.10/install/bin/sdl2-config" \
    USE_METAL=1 python3 tools/check_sdl_viewers.py
```

## Point classification

`alea_compute_classify()` classifies interleaved XYZ batches on a prepared
snapshot. It currently requires flat universe-0 geometry and scans all cells;
this is separate from hierarchical rendering and the core accelerated lookup.
An owner count of zero denotes a gap, one a unique owner, and larger counts an
overlap. IDs are -1 unless ownership is unique. Empty batches are valid. Inputs
must be finite and fit the selected backend's precision and allocation limits.

`make benchmark-compute USE_METAL=1` compares complete CPU/Metal batch costs.
`make benchmark-render-compute USE_METAL=1` compares the 64-sphere renderer.
Batch agreement is sampled evidence, not a geometry proof.
