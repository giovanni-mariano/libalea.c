<!-- SPDX-FileCopyrightText: 2026 Giovanni MARIANO -->
<!-- SPDX-License-Identifier: MPL-2.0 -->

# Native GPU implementation: M0 contracts and baseline

Checkpoint: `d3b7bff`, plus the M0 benchmark extension. Updated 2026-10-30.
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

## Explicit plotter operation and matched benchmarks

`mc_plotter --backend=metal --ownership=selected|all` now dispatches the requested
operation. Default Metal behavior remains all-owner classification. Selected
undefined fills are yellow, uncontained samples orange, and terminal void white.
Selected overlaps follow the first model-order path rather than displaying
all-owner magenta. The initial geometry probe also uses the requested operation.
Viewer publication, saving and navigation retain the chosen mode. CPU accepts
explicit selected mode and rejects an all-owner request.

The analytic benchmark now additionally compares the same 512x512 point array
through prepared CPU core lookup (serial) and Metal selection, alternating one
initial and five warm calls. IDs/materials match the analytic ownership fixture;
Metal status/depth are checked. In this run, warm medians are 36.911 ms for serial
CPU core and 26.547 ms for Metal. This measures generic point lookup, not optimized
CPU grid traversal or complete plotting.

The six-mode public complete-plotter benchmark uses five fresh-process trials,
four CPU threads, 512x512 and profiling disabled. All complete images match on
the nonoverlapping, fully populated public fixture:

| Mode | Query median | Query range | Process wall median |
| --- | ---: | ---: | ---: |
| CPU grid | 6.200 ms | 5.800–7.100 ms | 15.1 ms |
| CPU ray | 2.700 ms | 2.400–3.500 ms | 12.4 ms |
| Metal all, scanner | 35.600 ms | 32.400–37.200 ms | 68.6 ms |
| Metal all, candidates | 21.600 ms | 18.800–24.300 ms | 54.6 ms |
| Metal selected, scanner | 28.300 ms | 26.000–30.400 ms | 59.9 ms |
| Metal selected, candidates | 21.700 ms | 20.100–24.600 ms | 51.2 ms |

Selected/all image equality on this fixture does not generalize to models with
overlaps or undefined fills. The image checker tests those differing policies
explicitly, including undefined fills winning before a later containing sibling.

On the existing bounded private 128x128 view, one initial and five alternating
warm selected frames include coordinate generation, staging, query preparation,
all submissions and scattering. Scanner and candidate warm medians are 1.285 s
(range 1.230–1.364) and 1.144 s (1.114–1.269). First frames take 1.341/1.207 s;
snapshot preparation takes 11.791/9.828 s. All 16,384 result structs, including
status and depth, match prepared CPU ordered lookup. Optimized CPU grid has a
0.393 s warm median (0.376–0.418), while one serial CPU core-query pass takes
19.691 s. The private results reinforce the need to distinguish point-call
throughput from production grid traversal. Earlier all-owner timings are from
a separate run and are not an alternating selected-versus-all speed ratio.

No full-resolution private speed target is claimed. Private coordinates, images,
probe sources and raw logs remain ignored/outside Git.

The actual bounded private CLI check produces identical selected Metal images
with candidates off/on. CPU grid, Metal selected scanning and Metal selected
candidates query times are 377.5, 1302.4 and 1148.8 ms. Fresh-process wall times
are 15.761, 21.999 and 22.826 s, respectively, in one run each. The candidate
query improvement does not imply a cold-process improvement; loading/preparation
dominate and vary between runs. No private CPU/Metal image-equality claim is
made for undefined-fill presentation, which uses different colors.

Validation: public analytic point/ray oracles, complete-image checks for both
operations and candidate modes, overlap precedence, undefined-fill colors, option
validation, CPU plotter regressions, both-mode cancellation/reuse and selected
2D navigation/save checks pass. Existing 3D SDL regressions also pass. The Metal
plotter is rebuilt.

## Selected-query attribution and preparation memory

`ALEA_METAL_SELECT_PROFILE` now enables selection-only per-cell and per-depth
counters. Normal pipelines compile them out. The public counter tests check exact
counts, matching output, resets, interruption/numerical failures, report I/O retry,
hierarchy/lattice behavior and counter allocation failure/retry. Raw per-cell
reports and all private measurements below remain ignored; only aggregates appear
here. The report schema and profiling limitations are in `GPU_COMPUTE.md`.

### Native hardware counter result

The new public `benchmark-metal-counters` probe enumerates the device's exposed
counter sets and actually samples a compute pass. On this Apple M1 it reports
stage-boundary sampling supported, dispatch-boundary sampling unavailable, and
one set containing only `GPUTimestamp`. Two timestamp samples resolve successfully;
the sample delta agrees with the short command duration. The synthetic kernel's
output is independently verified.

This closes the question of whether native timestamp collection works. The public
API on this device exposes no occupancy, cache-miss or stall set. The earlier
empty System Trace exports do not prove the GPU did no work. Detailed hardware
stall attribution remains unavailable through the tested interface; no claim of
memory-bound or divergence-bound execution is made. Apple's
[counter-sampling documentation](https://developer.apple.com/documentation/metal/sampling-gpu-data-into-counter-sample-buffers)
describes the pass attachment used by this probe.

### Where selected work remains

The bounded private 128x128 view was run twice with the same candidate tiles,
comparing a normal context against an instrumented context. All 16,384 selected
result structs agree between the two modes on both runs. Both profiles have
identical counts, valid JSON, matching cell/depth totals and no saturated records.

| Executed work per frame | Root depth 0 | Child depth 1 | Total |
| --- | ---: | ---: | ---: |
| Cell visits | 16,318,416 | 7,394,707 | 23,713,123 |
| CSG instructions | 789,703,494 | 163,160,969 | 952,864,463 |
| Primitive evaluations | 140,748,688 | 30,801,341 | 171,550,029 |
| CSG rejections | 16,302,032 | 7,379,631 | 23,681,663 |
| Changed spatial-program selections | 11,904,354 | 0 | 11,904,354 |
| Lattice visits | 0 | 0 | 0 |

Root traversal accounts for **82.88% of executed instructions**. Overall 99.87% of
cell visits end in a positive CSG rejection on this sample. These are work
fractions, not elapsed-time fractions; the sample still visits no lattice cells.
There are 111,444 distinct visited packed cells across the frame. The most
expensive cell contributes 0.393% of instructions; the top 5, 10, 50 and 100
contribute respectively 1.730%, 2.794%, 7.978% and 13.017%. The work is broadly
distributed, rather than concentrated in a few pathological cells.

Normal frames take 1.246/1.303 s in this diagnostic run. Instrumented frames take
30.065/28.795 s; summed GPU durations are 29.636/28.538 s and candidate preparation
0.053/0.053 s. Global saturating atomics substantially perturb execution. These
timings quantify profiling overhead, not the performance of the default kernels.
No speed regression conclusion is drawn by comparing this two-frame diagnostic
run with earlier five-frame medians.

### Preparation and memory

A separate single-context run reads Mach task resident and maximum-resident byte
counts around load, snapshot creation, one instrumented point and destruction.
Load takes 9.645 s and snapshot preparation 9.732 s. The resident size after load
is 2,079,703,040 bytes; the maximum observed by Mach during the process reaches
3,033,710,592 bytes (**2.825 GiB**). Resident memory after the query is
2,693,103,616 bytes; allocator/framework retention after destruction is not
classified as a leak from this observation alone.

| Owned Metal buffers after one instrumented point | Bytes |
| --- | ---: |
| CSG bytecode (original and spatial programs) | 672,009,512 |
| Primitives | 30,718,968 |
| Cells | 16,223,328 |
| Spatial tree | 6,832,840 |
| Point hierarchy links | 5,407,776 |
| Attribution counters | 10,816,320 |
| Total distinct context buffers, including other query buffers | 744,077,548 |
| Metal process-wide allocated size | 744,570,880 |

The bytecode buffer alone is **640.88 MiB**, about 90% of the measured context
buffer total. Counter storage is approximately 10.32 MiB. These are allocated
buffer lengths; unified-memory CPU/GPU figures must not be added as separate
physical pools. Host geometry, compiler/pipeline allocations and peak preparation
scratch are not included in the context buffer sum; the Mach high-water mark
captures the process-level peak. The 128x128 capture's larger staging arrays bring
the context-buffer total to 744,192,208 bytes.

Source inspection shows that `compute_scene.c` sums each cell's full compiled
program count and emits a separate contiguous program per cell. This establishes
where storage is produced, not how much could be eliminated safely. Exact shared
roots and common subexpressions need an explicit reuse-opportunity audit before a
new representation is selected. Small spatial-code measurements from earlier
reports were not total bytecode-buffer measurements.

### Decision and remaining M0 work

Prioritize reducing broad root candidate/CSG work. Child-only optimization cannot
remove the observed root instruction volume. Specializing a handful of expensive
cells also leaves most work intact. The first representation experiment should
measure reusable program ranges/subexpressions and its effect on preparation,
code size and complete queries; instruction counts do not by themselves justify
runtime code generation or another per-point BVH.

For the named M1 private profile, use provisional regression budgets of 1 GiB in
distinct context buffer lengths, 4 GiB process resident high-water mark for one
loaded model/one context, and 15 s snapshot preparation (load reported separately).
These are benchmark guardrails for subsequent experiments, not enforced allocation
limits or universal backend guarantees. Record revisions explicitly. The measured
baseline fits these guards but still misses the one-second full-slice objective.

Fresh bounded private 3D exterior/cutaway/close-up profiling and final latency
acceptance remain open before M0 closes. CUDA hardware remains unavailable.

Validation for this checkpoint: 16 CPU and 73 Metal compute tests, 4 builder
tests, 5 resource tests and public image checks with selected attribution enabled
pass. All 99 emitted public JSON records parse and have matching cell/depth
totals. The counter probe target collects valid timestamp samples. The plotter
is rebuilt; normal kernels retain disabled instrumentation.

## Representation audit and compact Metal storage

Follow-up measured on Apple M1, 2026-10-30. The source-root audit finds 450,648
packed cells with 450,648 distinct source roots: identity-based whole-program
sharing saves nothing on this private workload. Grouping original emitted
programs by their hashes suggests an upper bound of about **0.72%** original
bytecode savings. Hash matches alone are not proof of equality; this estimate
would require exact comparisons before implementing sharing. Common-subexpression
reuse remains unaudited. These results do not justify a global interning pass.
Only aggregates are recorded; geometry, hashes, coordinates, probes and raw logs
remain ignored or outside Git.

The retained experiment is lossless instruction compaction, enabled by
`ALEA_METAL_COMPACT_CODE=1` at context creation (default `0`). Each instruction
uses one 32-bit word: 4 opcode bits and 28 argument bits. A complete fit check
precedes in-place host encoding. A wider opcode/argument retains the existing
8-byte format for the entire context. Logical instruction addresses, skip counts,
primitive references and instruction order remain unchanged. Host candidate proofs
and all four Metal pipeline variants decode the chosen representation. The CPU
backend and FP32 numerical contract are unchanged.

### Single-context memory comparison

Separate fresh processes load the same private model, prepare one context and
query one instrumented selected point with candidates enabled. This is a buffer
and memory observation, not a query-speed benchmark. Load times are 9.659 s wide
and 9.584 s compact. Profile counters are included in both buffer totals.

| Measurement | Wide | Compact |
| --- | ---: | ---: |
| Bytecode buffer bytes | 672,009,512 | 336,004,756 |
| Bytecode buffer MiB | 640.88 | 320.44 |
| Total distinct context buffer bytes | 744,077,548 | 408,072,792 |
| Total distinct context buffers MiB | 709.61 | 389.17 |
| Metal process-wide allocated bytes | 744,570,880 | 408,567,808 |
| Process resident high-water mark bytes | 3,834,413,056 | 3,500,474,368 |
| Snapshot preparation | 9.212 s | 9.131 s |

Code storage halves; total context buffers decrease **45.16%**. Process peak
falls only about 8.7%: host packing still uses the wide representation and other
model/preparation allocations remain. These process peaks differ from the earlier
checkpoint, so use this paired observation for the comparison rather than
subtracting results across increments. Both fit the provisional 4 GiB process,
1 GiB context-buffer and 15 s preparation guards. Buffer lengths and resident
memory are different measures of shared system memory and must not be added.

### Complete-query comparison

A separate retained-context benchmark prepares both formats, alternates them
through one initial and five warm 128x128 selected frames, and disables all
profiling. Both use the same four 64x64 candidate tiles. Timing includes coordinate
generation, proof construction, staging, all submissions and result scattering.

| Private selected frame | Wide | Compact |
| --- | ---: | ---: |
| Context preparation | 9.429 s | 9.243 s |
| First frame | 1.324 s | 1.204 s |
| Warm median | 1.234 s | 1.178 s |
| Warm range | 1.202-1.251 s | 1.164-1.351 s |
| Differences from prepared CPU core lookup | 0 / 16,384 | 0 / 16,384 |

The two Metal result arrays also match exactly. The modest median improvement
has overlapping ranges; this is not a robust speedup claim. In the same run the
optimized CPU grid median is 0.323 s and serial core lookup takes 20.250 s.
Compaction reduces storage rather than the number of executed CSG instructions.
The full-slice latency objective remains unmet.

The public analytic benchmark independently checks 512x512 all-owner/selected
points and 512x384 exterior/cutaway/close-up depths, normals and IDs in both
formats. One initial and five warm calls produce these Metal medians:

| Public operation | Wide ms | Compact ms |
| --- | ---: | ---: |
| All-owner points | 36.868 | 35.203 |
| Selected points | 26.368 | 27.017 |
| Exterior render | 8.480 | 7.953 |
| Cutaway render | 4.581 | 4.630 |
| Close-up render | 13.112 | 11.673 |

All independent checks pass. These small workloads do not establish a universal
performance improvement. To reproduce either format, run
`ALEA_METAL_COMPACT_CODE=0` or `=1` with `make benchmark-gpu-geometry USE_METAL=1
RELEASE=1`; the benchmark reports its complete initial/warm ranges.

Keep compaction opt-in. The next performance work remains broad root rejection
and ordered traversal, with exact boundary, overlap, unbounded/nonfinite and
undefined-fill behavior preserved. Another point hierarchy or subexpression
interpreter needs evidence beyond these storage measurements. Fresh bounded
private 3D attribution remains an M0 gate; CUDA execution awaits NVIDIA hardware.

Validation: CPU-only build passes 16 tests. Each Metal format passes 73 compute,
5 builder and 6 resource tests. Added checks cover exact decode at encoding
limits, skip arguments, in-place encoding, untouched wide fallback, interruption
and retry, real-buffer size, independent context formats and invalid-option
recovery. Compact-mode public plotter image checks and real SDL 2D/3D navigation,
progressive publication, cancellation and reuse checks pass. Both user-facing
viewers are rebuilt. Private artifacts remain ignored.
