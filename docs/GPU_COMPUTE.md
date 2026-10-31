<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO
SPDX-License-Identifier: MPL-2.0
-->

# GPU rendering and batch compute

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
make -C tools ../bin/render3d ../bin/mc_plotter USE_METAL=1 USE_SDL=1 RELEASE=1
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

Auxiliary render normals retain the source-surface orientation used by the CPU
renderer; shading applies camera-facing orientation to a separate copy. Clipping
faces use the negative clipping-plane normal. These visualization outputs do not
establish outward CSG boundary orientation for physics navigation.

## Intersection-group execution

The compute package includes native Metal execution of pure intersection subtrees
as ordered signed-leaf lists. Enable it before creating a context or launching a
tool:

```sh
ALEA_METAL_CSG_GROUPS=1 bin/mc_plotter model.mcnp Z 0 -10 10 -10 10 512 slice.bmp --backend=metal
```

This applies to `alea_compute_classify()` and `alea_compute_select()`, including
`mc_plotter`'s all-owner and selected ownership modes. It is a supported package
option; no benchmark harness or extra metallib is needed. Group buffers are owned
by each immutable context, so different contexts remain independent after the
source model is changed or freed.

The default is `0`: preparation and up to 64 MiB of additional maps/lists can be
worthwhile for large CSG models, while small expressions can be slower. Values
other than `0` or `1` fail context creation. If no eligible groups exist or optional
storage cannot fit its budget/device limits, the same Metal context uses generic
bytecode. Cancellation and malformed geometry still report errors. Grouping is
compatible with compact code, jump threading, spatial programs and point
candidates, and preserves evaluation order, finite errors, void and overlaps.

Point and selection profile records include `csg_groups` and `instruction_unit`.
When groups are active, `bytecode_and_group_ops` counts remaining bytecode plus
one group dispatch and its visited list entries. These counts are not directly
comparable to bytecode-only instruction totals. Primitive-evaluation counts keep
the same meaning. The ray evaluator remains independent of this point option.
Ray-group contract and coverage are documented in
[the ray report](NATIVE_GPU_RAY_GROUPS.md).

See [integration and validation](NATIVE_GPU_CSG_GROUPS.md) and the earlier
[controlled measurements](NATIVE_GPU_GROUP_EXECUTION.md).

## Native intersection groups for ordinary rays

`ALEA_METAL_RAY_GROUPS=1` at context creation enables supported native group
execution in `alea_compute_render()` and static `render3d --backend metal` output.
Its default is `0`. It is independent of the point-group option; both share one
bounded context-owned data set. Generic Metal evaluation handles missing groups
or unavailable optional storage. Progressive rendering keeps its generic
resumable evaluator. See [the package guide](NATIVE_GPU_RAY_GROUPS_PACKAGE.md) for
ownership, profiling units, validation and measured tradeoffs.

## Experimental compact Metal instructions

Set `ALEA_METAL_COMPACT_CODE=1` before snapshot creation to store each CSG
instruction in one 32-bit word instead of two. The default remains `0` for
comparison. Classification, selection, host candidate proofs and both render
paths use the context's chosen format. Independent contexts may use different
formats; changing the environment does not modify existing contexts.

The encoding preserves instruction order, logical addresses and arguments. It
uses 4 opcode bits and 28 argument bits. If any instruction does not fit, the
whole context retains the wide format without truncation or a new geometry
rejection. Both formats use the same FP32 arithmetic and retain its accuracy
limits. Invalid settings fail context creation.

On the measured private workload, compact storage halves the code buffer and
reduces total context buffers by about 45%. Host preprocessing still builds wide
instructions, so process peak memory does not halve. Query speed is workload
specific; use complete queries with profiling disabled for comparisons. Aggregate
measurements and remaining work are in [the M0 report](NATIVE_GPU_M0.md).

## Ordered batch selection

`alea_compute_select()` is a separate CPU FP64 / Metal FP32 operation for selecting
one cell in original model order. It returns `alea_compute_selection_t` with cell
and material IDs, root-based depth and a per-point status:

- `ALEA_COMPUTE_FOUND`: terminal material or terminal void.
- `ALEA_COMPUTE_NOT_FOUND`: no containing root cell; IDs/depth are -1.
- `ALEA_COMPUTE_UNDEFINED_FILL`: the winning container has no child at the point;
  return that container and stop, including empty rectangular lattice slots.

Overlapping siblings never override the first winning path. This differs from
`alea_compute_classify()`, which counts every terminal owner and reports a hole
in a fill as a gap. Each operation has its own result type. Selection currently
uses ordered CPU enumeration and a native Metal single-path traversal; it is
connected to `mc_plotter --ownership=selected`. Matched point and complete-frame
measurements are recorded in the native GPU implementation report.

Snapshots, explicit precision and numerical limits remain as described above.
Rectangular lattices are supported; hexagonal lattices and lattice-local fill
transforms remain rejected. Metal selection reuses the geometry, transforms,
spatial programs and optional root candidate lists; its pipeline is prepared
lazily. Selection and classification may alternate on a serialized context.
Result buffers grow for either layout and remain reusable after errors. Selection
submits at most 4,096 points at a time, with no per-submission time guarantee.
`ALEA_METAL_POINT_PROFILE` continues to report classification calls only.

```c
alea_compute_selection_t result;
double xyz[3] = {0, 0, 0};
if (alea_compute_select(gpu, xyz, 1, &result) == 0) {
    /* Inspect result.status before consuming cell_id/material_id/depth. */
}
```

## 2D Metal slices

The plotter now uses the native point classifier with explicit backend selection:

```sh
bin/mc_plotter examples/data/doc_shield.i Z 0 -12 12 -12 12 800x600 \
    slice.png --backend=metal --interactive
```

Omit `--interactive` to write a PNG or BMP directly. CPU remains the default;
`--backend metal` and `--backend=metal` are equivalent. Metal requires a build
with `USE_METAL=1`; SDL is needed only for interactive use.

Choose the ownership operation explicitly:

```sh
bin/mc_plotter model.mcnp Z 0 -10 10 -10 10 512 slice.bmp --backend=metal --ownership=selected
bin/mc_plotter model.mcnp Z 0 -10 10 -10 10 512 slice.bmp --backend=metal --ownership=all
```

Metal defaults to `all`, preserving terminal-owner counting and overlap display.
`selected` follows the first containing path in model order, including terminal
void. Undefined fills return their container and appear yellow; uncontained
points remain orange. Overlapping siblings do not override the winning path and
are not marked magenta in selected mode. The viewer title, help and console
identify the selected operation. The mode is fixed for that viewer session;
navigation, saving and snapshot reuse preserve it. Both `--ownership=value` and
`--ownership value` are accepted. CPU plotting supports explicit `selected`;
requesting `all` on CPU is rejected. Existing CPU diagnostics remain separately
selected through their normal options.

Each pixel samples its center using approximate FP32 ownership. Unique material
cells use the selected cell/material palette, terminal void is white, a sampled
gap is orange, and a sampled overlap is magenta. These are point samples, not
verified geometry diagnostics: boundaries may differ from CPU evaluation and
thin features can be missed. Contours, cell/material labels and ticks are drawn
on the CPU from the classified grid. Pan, zoom, axis and slice controls reuse
one geometry snapshot; `H` shows controls and the color legend.

Metal slices support finite/repeating rectangular lattices and outer fills.
They reject hexagonal lattices, lattice-local fill transforms, batch files, ray sampling, surface
labels, `--errors`, `--debug`, `--trace`, and the `ALEA_PLOT_*` coverage/error
diagnostic settings. These requests fail explicitly without CPU fallback.
The Ray and Errors controls are absent from the Metal viewer.

The worker checks cancellation between batches of 4,096 samples and before
publishing a complete image. One GPU submission has no latency bound; no partial
slice is published. Saving uses a complete frame for the requested view.
The classifier now selects certified spatial CSG programs to reduce expression
work. It still enumerates the cells of each visited universe. Complex-model
interactive performance remains unresolved.

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
python3 tools/check_metal_plotter.py
```

## Point classification

`alea_compute_classify()` classifies interleaved XYZ batches on an owned
snapshot. CPU FP64 and native Metal FP32 support nested transformed fills through
32 universe levels. Classification counts terminal owning occurrences, including
void, under every containing ancestor. Repeated placements count separately;
containers do not count as owners, and gaps in their contents remain gaps.
An owner count of zero denotes a gap, one a unique owner, and larger counts an
overlap. IDs are -1 unless ownership is unique. Unreferenced definitions are
ignored. Finite/repeating rectangular lattices are traversed procedurally, including
nested lattices, outer fills and zero-element coordinates. Element selection
uses floor on tiled axes, assigning a shared face to the upper element. Finite
singleton Y/Z windows include their upper face, matching core point lookup;
untiled Z remains unshifted in center coordinates. Nonpositive element fills
are gaps, while void-only fill universes remain unique owners. Ancestor fill
CSG still clips the lattice. Hexagonal lattices and lattice-local fill transforms
are rejected. Metal validates lattice type, dimensions, pitch and fill tables
before CSG packing and rendering specialization. Spatial CSG programs are prepared by snapshot creation. Classification now
selects them only after a separate finite-arithmetic certificate establishes
that simplification cannot hide an error in the original expression.

Metal rejects logical or relative element indices with magnitude greater than
16,777,215 rather than silently collapsing adjacent indices in FP32. CPU lookup
uses checked 32-bit indices. These guards do not remove ordinary FP32 boundary
sensitivity.

Metal shares primitive/CSG, transform and spatial-program buffers with its renderer, but uses
separate point-traversal and lattice-fill links that include void. These links
upload lazily and survive either classification/rendering initialization order.
Inside a certified region, the kernel selects a shorter spatial CSG program;
outside it, or when certification is unavailable, it evaluates the original
expression. Certification inspects every original primitive, including potentially
short-circuited operands, using outward-rounded FP32 interval bounds. A bounded
certification budget leaves remaining cells on their original programs. The
renderer’s visibility BVH and content bounds do not reject point-query owners.

`ALEA_METAL_POINT_SPATIAL=0` disables point-program selection at snapshot creation
for comparison; the default is enabled. This does not disable the renderer’s
spatial preprocessing. `ALEA_METAL_SPATIAL=0` disables that preprocessing itself,
which also leaves no spatial programs for classification. Neither switch changes
backend or precision. The simple public 512x512 fixture remains around 33–34 ms;
the measured private 256-point sample improves about 19% in warm query time, but
is still too slow to establish interactive performance on large slices.
Hierarchical queries submit up to 4,096 points at a time and check cancellation
between submissions; per-submission execution time is not bounded. Both point
implementations currently enumerate visited universes. They are separate from
the core accelerated point lookup. Empty batches are valid. Inputs must be
finite and fit the backend's precision/allocation limits. FP32 boundary and
thin-feature limitations remain; sampled ownership is not interval coverage.

### Experimental point candidates

`ALEA_METAL_POINT_CANDIDATES=1` at snapshot creation enables conservative
root-universe candidate filtering. It is **off by default**; `0` retains the
scanner. For each batch of up to 4,096 points, the host bounds the actual converted
FP32 inputs, proves which root cells cannot contain any of them, and uploads a
linked candidate list through shared Metal storage. The kernel reuses that list
across the batch. Child universes still use the existing traversal and certified
spatial CSG programs. Both options can be used independently.

The proof propagates strict signs through the original CSG bytecode and requires
every original operand to evaluate finitely over the region, including operands
that might be short-circuited. Unknown proofs, boundary-touching cells and lattice
containers remain candidates. Void and overlapping cells stay eligible. Primitive
proofs are reused within each batch; a 32-million-source-instruction budget leaves
remaining programs unfiltered. Interrupt checks occur between cells and every
1,024 source instructions. An individual GPU submission still has no latency
guarantee.

Candidate lists are rebuilt per batch and query; there is no view cache to
invalidate. Scratch storage is four bytes per packed cell plus one byte per
packed primitive, allocated lazily per context. Existing shared code/primitive
buffers supply the host proofs without another geometry copy. Allocation failure
returns an error and allows retry. Different contexts own independent scratch.

For `mc_plotter`, the same option stages Metal pixels in 64x64 tiles, including
partial edge tiles, and scatters results back into the image. This makes each
proof cover a compact region. Generic `alea_compute_classify()` callers keep
their supplied point order; broad or incoherent batches may gain little, and
splitting a small batch can lose GPU parallelism. Profile and benchmark the
complete query, including proof construction, before selecting this option.

```sh
ALEA_METAL_POINT_CANDIDATES=1 bin/mc_plotter model.mcnp Z 0 -10 10 -10 10 512 --backend=metal slice.bmp
python3 benchmarks/gpu_plotter.py
```

The benchmark generates public geometry and compares complete images and timings
for CPU grid/ray plotting and Metal scanning/candidates. It uses four CPU threads
by default and reports query time separately from fresh-process wall time. Unlike
the compute benchmark, its CPU measurements exercise the optimized plotter paths.
CPU selected-owner plotting and Metal all-owner sampling have different overlap
semantics; image equality on the public nonoverlapping fixture does not imply
that arbitrary models will produce identical images.

### Selected-query attribution

`ALEA_METAL_SELECT_PROFILE=build/selected.jsonl` at context creation enables
selection-only cell/depth attribution. CPU contexts ignore it. The normal Metal
selection pipeline specializes instrumentation out and allocates no attribution
buffer. This setting is independent of `ALEA_METAL_POINT_PROFILE`, which continues
to report classification only.

Each successful nonempty selection appends one version-1 JSON record. Calls that
fail or are interrupted before reporting emit no partial work. Report I/O errors
return failure and may leave an incomplete line; fixing the path permits retry.
Use separate files for concurrently used contexts. Counters reset per public call
and accumulate across its internal submissions.

`cells` contains visited packed cell indices, original program sizes and counts
of calls, executed instructions, primitive evaluations, CSG rejections, changed
spatial-program selections and lattice visits. `depths` reports the same counters
by universe depth (root = 0). CSG instructions include executed skip instructions
but exclude skipped operands. A lattice visit has zero CSG work; an outside
lattice window is not a CSG rejection. Per-record counters saturate at UINT32_MAX;
`counter_limit_records` counts saturated cell/depth records, whose values then
provide lower bounds. Host analysis must sum into 64-bit values.

Reports include bridge wall time, available command-buffer GPU time and valid
timestamp-pair count, candidate preparation time, context buffer lengths, and
selected geometry/counter buffer sizes. `device_process_allocated_bytes` is the
Metal device's process-wide allocation figure, not a per-context or physical RAM
measurement. Context buffer totals count each owned shared buffer once and
exclude pipeline/driver memory and host geometry. Apple unified memory figures
must not be added together as independent CPU and GPU memory pools.

Attribution uses saturating GPU atomics and can be very expensive. On the bounded
private 128x128 sample, normal candidate frames take about 1.3 s in the attribution
run while instrumented frames take about 29–30 s. These times diagnose profiling
overhead; use disabled profiling for acceptance measurements. Attribution storage
is 24 bytes per packed cell plus 32 depth records, allocated lazily. The files can
contain model-sensitive program/cell metadata: private-model reports belong in
ignored `build/` or outside Git, with only aggregates retained in tracked reports.

`make benchmark-metal-counters USE_METAL=1 RELEASE=1` enumerates native counter
sets and samples supported sets on public synthetic work. It verifies execution
and checks timestamp validity. A device exposing only timestamp counters does
not supply occupancy or stall attribution through this API. See Apple's
[counter sampling documentation](https://developer.apple.com/documentation/metal/sampling-gpu-data-into-counter-sample-buffers)
for the pass-boundary sampling mechanism.

### Point-query profiling

Set `ALEA_METAL_POINT_PROFILE=build/point-profile.jsonl` before creating a Metal
context (or launching `mc_plotter --backend=metal`) to append aggregate JSONL
reports. Create the output directory first. Unset or empty disables profiling;
the normal pipeline specializes out the shader counters and allocates no counter
buffer. CPU contexts ignore this option. It is separate from the 3D renderer's
`ALEA_METAL_PROFILE` setting and does not instrument rays.

Each successful nonempty classification appends one version-1 record covering
all its internal submissions. Empty, interrupted and numerically failed queries
do not report partial work. Report I/O failure returns an error; results are
then unspecified, and an I/O failure can leave a partial report line. The context
can be reused after fixing the path. Use separate paths for concurrently used
contexts. Reports contain no coordinates, cell IDs or model identifiers.

| Fields | Meaning |
| --- | --- |
| `points`, `submissions` | Input points and completed GPU submissions |
| `cells`, `lattices` | Cell visits, and the subset selecting lattice elements rather than evaluating CSG |
| `universes`, `descents`, `max_depth` | Root entries plus successful descents, descent count, and maximum universe depth (root = 1) |
| `original_calls`, `spatial_calls` | Original versus changed spatial CSG program evaluations |
| `no_program`, `uncertified`, `outside_domain`, `unchanged_region` | Reasons for original evaluations: no cache; no point certificate (also includes deliberately disabled selection); outside certified domain; selected cache leaf has the original program |
| `original_instructions`, `spatial_instructions` | Actually executed bytecode, including skip instructions but excluding skipped operands |
| `original_primitives`, `spatial_primitives` | Executed primitive evaluations |
| `rejected`, `owners` | Positive CSG evaluations and terminal owning occurrences; lattice gaps are not CSG rejections |
| `rejected_instructions`, `accepted_instructions` | Executed bytecode in rejecting versus containing CSG evaluations; containing ancestors are included |
| `candidate_seconds` | Host time bounding inputs and building candidates, included in `wall_seconds` |
| `candidate_cells`, `candidate_rejected` | Root cells examined by the host and cells removed from batch lists (not multiplied by point count) |
| `candidate_instructions`, `candidate_primitives` | Host source-instruction visits and primitive interval proofs; these differ from executed GPU instructions |
| `candidate_budget_fallbacks` | Root programs retained because the per-batch proof budget was insufficient |
| `counter_limit_points` | Points with at least one saturated 32-bit counter; nonzero means counts can be lower bounds |
| `wall_seconds` | Metal bridge time, including lazy uploads, conversion, waits, result copies and counter aggregation; excludes public input validation and final report writing |
| `gpu_seconds`, `gpu_timed_submissions` | Sum of available command-buffer GPU durations and number of valid timestamp pairs; zero valid pairs means GPU time is unavailable |

Host sums use 64-bit counters; `max_depth` is a maximum. Counter volume and
instrumentation overhead depend on the workload: the buffer uses 76 bytes per
submitted point, capped at 4,096 points for hierarchical queries. Flat queries
retain their whole-batch submission unless candidate filtering is enabled. Measure performance with profiling disabled;
in one public 512x512 run instrumentation increased the warm median from 30.658
to 37.853 ms. Instruction fractions describe work, not elapsed-time fractions.

`make benchmark-compute USE_METAL=1` compares complete flat CPU/Metal batch costs.
`make benchmark-render-compute USE_METAL=1` compares the 64-sphere renderer.
`make benchmark-gpu-geometry USE_METAL=1 RELEASE=1` checks and times hierarchical
512x512 point-sampled slices and 512x384 exterior/cutaway frames against analytic
references. Metal owns its kernels, device layouts and scheduling. Planned CUDA
and other backends will have their own native implementations, sharing libalea's
query contract and semantic fixtures rather than a GPU portability framework.
