<!-- SPDX-FileCopyrightText: 2026 Giovanni MARIANO -->
<!-- SPDX-License-Identifier: MPL-2.0 -->

# Native GPU implementation: M0 contracts and baseline

Checkpoint: `d3b7bff`, plus the M0 benchmark extension. Updated 2026-10-29.
This records the first implementation work package from the local native GPU
roadmap. Selected-owner semantics below are now implemented by the follow-up recorded at
the end of this report; the capability API remains a design decision. Metal is active; CUDA device work awaits NVIDIA hardware.

## Audited geometry and feature gaps

The implementation inventory must describe geometry after import as well as the
public primitive enumeration. The Metal packer consumes reachable CSG nodes in
`src/compute/compute_scene.c`. MCNP macrobody conversion already calls
`alea_expand_macrobody_immediate()` in `src/mcnp/conversion/surface_conv.c`.
`src/core/alea_macrobody.c` implements analytic expansion for all macrobody
families. Consequently, an unsupported raw macrobody tag in the Metal packer does
not prove that an imported model using that macrobody is unsupported.

| Area | Current Metal behavior | Required work / evidence |
| --- | --- | --- |
| Planes, spheres, axis cylinders/cones/tori, quadrics, RPP, RCC | Direct packed evaluation and analytic ray support | Retain independent primitive/boundary/root tests and numerical limits |
| BOX, SPH, TRC, ELL, REC, WED, RHP, ARB | Raw tags rejected by packer's default case; MCNP import can analytically lower them to supported CSG | Add imported-versus-direct-construction fixture pairs for every family; verify facet provenance before declaring parity |
| CSG | Union/intersection/difference/complement, short-circuit bytecode | Preserve original error behavior in reductions; control duplication/code size |
| Transformed fills | Packed inverse matrices, shared definitions | Translation/rotation fixtures exist; audit remaining accepted matrix conventions and FP32 range limits |
| Rectangular lattices | Finite/repeating/nested, outer fills, empty slots, zero coordinates | Preserve face/window/index rules in all query modes |
| Hexagonal lattices | Rejected | Implement procedural element lookup, windows, ray crossings and numerical tests |
| Lattice-local fill transforms | Rejected | Define transform application order against core; implement both point and ray paths |
| Selected-owner compute query | Implemented by follow-up on CPU/Metal | Benchmark against optimized CPU selection and integrate the chosen frontend mode |
| All-owner compute query | Implemented; includes terminal void and repeated occurrences | Accelerate child-universe scans; preserve counts and ancestor clipping |
| 2D interval/ray-row GPU slicing | Absent | Dedicated query and rasterization operation with explicit overlap policy |
| 2D static/interactive sampling | Implemented, arbitrary planes; optional 64x64 candidate tiles | Bound long submissions and complete remaining frontend workflows |
| GPU plotter batch files / surface labels | Rejected | Batch scheduling and original surface/facet provenance |
| Exact coverage diagnostics | CPU operations; GPU samples approximate | Keep explicit operation/backend provenance; sampled coverage is not exact diagnostics |
| 3D solid/depth/cell/material output | Implemented; normals available | Retain analytic first-hit/depth/normal checks, ordinary/progressive consistency |
| Perspective / orthographic | Implemented | Retain camera and clipping tests |
| Clipping | Zero or one plane | Core permits up to 16 and AND/OR combinations; add multi-plane combinations |
| Filters / density range | Rejected | Add cell/material include/exclude and density-range semantics |
| Color modes | Material/cell plus custom host palette | Add universe/density modes with occurrence-aware metadata |
| Shadows / X-ray / adaptive AA | Rejected | Implement separately with tests; AA 1..4 is currently supported |
| Progressive rays | Implemented, including CSG continuation | Submission latency is adaptive, not guaranteed; validate expensive primitive tails |
| Progressive point/slice work | No continuation inside long query | Add resumable state and generation-safe publication |
| Capability/device options | Backend availability only; experimental environment settings | Typed feature/precision/limit query and validated options |
| Physics boundary/transition operations | Not exposed by compute API | Separate navigation contract and numerical validation before use in solvers |

Source limits: 64 CSG value-stack entries, 32 universe levels, source CSG depth
4096; Metal rejects rectangular element indices beyond magnitude 16,777,215.
Actual buffer/device allocation limits also apply. Continuation, spilling or
explicit failure must replace any future capacity exhaustion; never truncate.
The importer/direct-construction fixture matrix remains an M1 acceptance item,
not a claim of tested macrobody coverage from this source audit alone.

## Operation contracts frozen for the next implementation

### Selected-owner points

Add a separate batch operation; do not change `alea_compute_classify()` to stop
at the first owner. C name: `alea_compute_select()`. Its per-point
result has cell ID, material ID, depth and resolution status. The public layout is in `alea_compute.h`; the following semantics are fixed:

1. Start at universe zero. Visit candidates in the original universe/model order.
2. Choose the first containing cell, then resolve its fill/lattice path. Do not
   choose a later overlapping sibling because its descendant is deeper.
3. If the winning containing fill has no child at the point, return that container
   with `UNDEFINED_FILL`. Do not continue to a sibling or invent terminal void.
   This follows `find_cell_recursive()` in `src/core/alea_universe.c`.
4. Terminal material and terminal void both return `FOUND`. A genuinely uncontained
   point returns `NOT_FOUND`, IDs -1 and depth -1. Not-found is a per-point result,
   not a failed batch. Undefined fill returns the selected container IDs/depth.
5. Invalid geometry, arithmetic, stack/capacity or execution failures return an
   error; outputs are unspecified on failed batches. Empty batches remain valid.
6. Coordinates are transformed at each occurrence; lattice selection follows the
   same documented core convention. Acceleration preserves ordered precedence.

The CPU oracle is prepared core selected lookup, with undefined-fill metadata
checked through the internal/public resolution interfaces. Existing core APIs
remain compatible. FP32 versus FP64 differences near boundaries are separately
classified; strict selected-owner equivalence is required on well-separated cases.

### All-owner points

The existing `alea_compute_classify()` contract stays authoritative: count every
terminal occurrence under all containing ancestors, including terminal void.
Repeated placements count separately. Containers do not count; a hole in their
children remains a gap. IDs are supplied only for exactly one owner. This contract
can produce different output from selected-owner lookup for both overlaps and
undefined fills; neither result is a substitute for the other.

CPU compute enumeration and independent analytic fixtures remain the correctness
oracles. Public comparisons label exhaustive enumeration versus optimized core
queries. All-owner acceleration must retain uncertain/unbounded candidates unless
a conservative proof justifies removing them, including finite-error behavior.

### Rays, intervals and future navigation

First visible material, complete terminal-region intervals and next physical
boundary are separate operations. Visibility can omit void; navigation cannot.
Interval outputs must specify ordering, overlap representation, clipping, finite
ray range and original surface/occurrence provenance. The existing image API is
not promoted into an implicit boundary-navigation API.

For N1, define normal orientation, parameter units, tangency/coincident-face and
boundary-start behavior before exposing an ABI. An occurrence identifier must be
unambiguous within its documented lifetime; an unchecked hash is insufficient.
Numerical uncertainty/error statuses must have a documented detection guarantee.
FP32 visualization and FP64 arithmetic alone do not establish transport validity.

### Capability and lifetime boundary

The future capability query reports backend/device identity, numerical modes,
operation/geometry features and published capacity limits. Features become
advertised only after their implementation and real-device tests pass. Requesting
unsupported operations or precision fails explicitly. Snapshot ownership and
per-context serialization remain explicit; job cancellation cannot free in-flight
native resources. CUDA's packer, queues and layouts will be native CUDA designs.

## Workload manifest and fresh public baseline

`benchmarks/native_gpu_workloads.json` is the machine-readable workload inventory.
Commands are argument arrays relative to the repository root. Setup/build happens
before measurement; entry points retain their own validation. Profiles are disabled
for timing and captured separately. CUDA entries are pending hardware, not runnable
claims. Private view coordinates and commands belong only in ignored local files.

The analytic benchmark now covers a close-up in addition to exterior and cutaway,
and uses one initial plus five warm calls. It checks every sampled owner against
independent equations and every ray against analytic first-hit/depth/normals.
Fresh Apple M1 run, four CPU threads, release build, native Metal, FP32, candidate
filtering off and default spatial programs enabled:

| Task | Backend | First | Warm median | Warm range |
| --- | --- | ---: | ---: | ---: |
| 512x512 all-owner points | CPU FP64 enumeration | 53.433 ms | 44.113 ms | 42.598–45.804 ms |
| 512x512 all-owner points | Metal | 37.073 ms | 33.894 ms | 32.420–35.080 ms |
| 512x384 exterior | Metal | 11.837 ms | 7.892 ms | 7.188–9.558 ms |
| 512x384 cutaway | Metal | 4.575 ms | 5.155 ms | 4.618–6.866 ms |
| 512x384 close-up | Metal | 13.894 ms | 11.901 ms | 9.881–12.888 ms |

All analytic checks pass. Setup is 0.066 ms CPU / 59.173 ms Metal. The three ray
views share a context, so cutaway/close-up first times are not cold-context times.
These simple public scenes do not establish complex-model interactivity or an
optimized CPU plotting speed advantage. Existing complete-plotter comparisons and
bounded private evidence are recorded in `NATIVE_GPU_GEOMETRY.md`.

## Native profiling availability and evidence limits

The installed Xcode tools include `xctrace` and the Metal System Trace template.
A short capture of `synthetic:mixed`, progressive rendering at 128x96, completed
successfully. Trace and exports are kept in ignored `build/native-gpu/`.

The trace lists GPU timing/counter/shader-profiler schemas, but exports contain
**zero rows** for `gpu-counter-info`, `gpu-counter-value` and
`gpu-shader-profiler-sample`. This capture establishes functioning trace recording,
not hardware counter availability or a diagnosis of occupancy, register spills,
cache misses or divergence. It does not prove the device cannot provide those
metrics. A configured native counter capture remains required.

Current point counters establish excessive surviving CSG work; they do not isolate
hardware stalls. Current private 128x128 evidence still reports approximately
1.68 billion executed instructions and 40 million cell visits after candidate
filtering, with 2.198 s warm frame time. Preserve that measured baseline.

## M0 status and first implementation decisions

Completed: source-level feature inventory; selected/all-owner resolution contracts;
public workload manifest; five-warm-call analytic 2D/three-view 3D baseline; native
trace availability check; CUDA hardware dependency recorded.

Still open before closing M0: configure usable hardware counters or record a
specific tooling limitation; attribute expensive surviving programs/child scans;
record packed/resident/scratch memory and preparation cost by workload; obtain
fresh bounded private 3D exterior/cutaway/close-up evidence; freeze per-device
resource budgets and validate proposed latency/performance targets. Existing
aggregate private 2D evidence remains useful but cannot replace these 3D checks.

The first functional implementation item is the selected-owner batch operation,
with an explicit result status and tests for ordered overlap and undefined-fill
behavior. It provides the missing matched comparison with CPU production lookup.
It does not change the all-owner target or explain away its remaining cost.
The first acceleration change is chosen after survivor/hardware attribution;
no additional BVH or expression-generation strategy is selected from instruction
counts alone. Geometry parity starts with imported/direct macrobody fixture pairs,
using existing analytic lowering before adding duplicate primitive implementations.

CPU-only/Metal builds remain independent of CUDA tooling. No unavailable CUDA
backend is advertised. Private models, paths, coordinates, hashes, captures,
generated geometry-specific code, images and raw logs remain ignored/outside Git.
Only aggregate private observations enter tracked reports.

## Selected-owner API implementation follow-up

`alea_compute_select()` now implements the selected-owner contract on CPU FP64
and native Metal FP32. Its public result is `alea_compute_selection_t`, with
cell/material IDs, zero-based depth and FOUND/NOT_FOUND/UNDEFINED_FILL status.
The CPU implementation is an ordered recursive reference. The Metal path follows
one winning occurrence without sibling backtracking, shares existing geometry
and optional candidates, and prepares its selection pipeline lazily. Current
snapshot geometry/precision limits remain applicable. Plotter integration and
matched performance measurements are still pending.

Selection testing exposed a core ordered-lookup defect: a zero lattice fill
could recurse into root universe zero. Ordered CPU lookup now returns the lattice
container as undefined fill immediately, matching its existing undefined-child
rule. A core regression checks both the returned cell and root depth.

Validation: 16 CPU and 72 Metal compute tests, 4 builder tests, 4 resource tests,
6 core ordered-point tests and 3 universe-BVH tests pass. Selection comparisons span
4,105 points
per scene, transformed overlaps, finite/repeating lattices, empty slots, source
lifetime, buffer reuse across both operations, candidate modes, cancellation,
invalid storage/counts/coordinates, and short-circuiting before later nonfinite
programs. Selection does not emit classification profile records.

The user-facing Metal plotter is rebuilt. CPU and Metal public image checks and
SDL 2D/3D navigation, cancellation and reuse checks pass after this follow-up.
