# Plan: native hierarchical lattice batch DDA

**Date:** 2026-08-07

**Status:** implementation in progress — phases 0–3 complete; phase 4
substantially complete; phase 5 initial baseline complete

**Target branch:** `feature/unified-ray-query`

**Primary implementation areas:** `src/raycast/raycast.c`,
`src/raycast/raycast.h`, `src/raycast/raycast_api.c`, and
`src/core/alea_spatial_hier.*`

**Related plans:**

- `plans/UNIFIED_RAY_QUERY_ARCHITECTURE.md`
- `plans/NATIVE_COMPACT_RAYCAST_BATCH.md`
- `plans/PLAN_LATTICE_AWARE_GRID_COHERENCE.md`

## Executive decision

Implement native batch `FIRST_VISIBLE`/`ANY_HIT` on top of the existing,
nested-lattice-capable hierarchical stepper. First refactor that scalar stepper
into a resumable, interval-at-a-time state machine. Then execute a bounded
packet of those same lane states per worker and write completed answers
directly into the existing structure-of-arrays batch results.

The batch implementation must not introduce a second ownership resolver or a
second DDA algorithm. Scalar full traces, scalar early-stop queries, and batch
early-stop queries must consume the same verified interval step primitive.

No public API change is required.

## Current baseline

The current branch already provides:

- canonical rectangular and hexagonal finite-lattice location and DDA helpers;
- forward and reverse scalar lattice parity for ordinary one-level lattices;
- scalar hierarchical lattice `FIRST_VISIBLE` and `ANY_HIT` early stopping;
- transformed lattice-placement entry validation against enclosing fills;
- synthetic lattice boundary semantics (`surface_id == 0`);
- physical surface/primitive/normal preservation;
- `t_min` cross-section semantics (`surface_id == -1` when starting inside an
  already-active interval);
- transactional SoA batch result allocation and publication;
- OpenMP orchestration with deterministic output indexing;
- rect/hex scalar and batch policy parity tests.

The former batch traversal was scalar dispatch:

```text
for each ray assigned to a worker:
    initialize ray
    call alea_raycast_hier_first_visible_nocache()
    copy one scalar result into the SoA output
```

`ANY_HIT` likewise delegates to a scalar query. This is parallel scalar
execution, not a native batch stepper.

Both early-stop batch APIs now use the native packet executor. The packet
width is currently a private compile-time constant of four lanes. Packets
cover contiguous input ranges and retire lanes independently; this initial
implementation does not refill an emptied lane from a later packet range.

The nested-lattice point fixture and hierarchical full/early-stop ray paths
resolve correctly. The legacy `GLOBAL` lattice pipeline remains unable to
compose nested lattice-element transforms because it scans lattice cells in
world coordinates. It is not the implementation reference for this program:
the hierarchical fast-forward trace is the nested-lattice oracle.

## Goals

1. Make the hierarchical ray traversal resumable at verified-interval
   boundaries.
2. Preserve the existing exact hierarchy occurrence identity when restarting
   below a changed lattice entry.
3. Keep nested-lattice full, early-stop, and batch behavior on the same
   hierarchical interval semantics.
4. Make scalar full trace, scalar `FIRST_VISIBLE`, scalar `ANY_HIT`, batch
   `FIRST_VISIBLE`, and batch `ANY_HIT` share the same interval stepper.
5. Remove per-ray scalar entry-point dispatch from the early-stop batch APIs.
6. Avoid complete hit/segment materialization for early-stop batch queries.
7. Bound traversal memory by worker count and packet width, not total ray
   count.
8. Preserve deterministic output under different OpenMP schedules and thread
   counts.
9. Preserve interruption, byte-budget, error, and transactional-publication
   behavior.
10. Establish performance evidence before increasing packet width or adding
    explicit SIMD intrinsics.

## Non-goals

- Do not redesign the public ray query or batch result APIs.
- Do not replace the compact segment CSR executor in this program.
- Do not flatten or materialize lattice elements globally.
- Do not add a separate batch-only cell ownership algorithm.
- Do not require AVX/SIMD in the first native batch implementation.
- Do not reorder input rays globally; contiguous packetization may be used,
  but outputs must remain indexed by original ray order.
- Do not weaken overlap/deck-order ownership rules.
- Do not change lattice parsing, coordinate conventions, or finite/repeating
  fill semantics.
- Do not make timing thresholds ordinary correctness-test assertions.

## Required semantic invariants

### Verified intervals

- A query answer is published only from a completed, ownership-verified open
  interval `[t_current, t_next)`.
- The next boundary is the nearest valid participant among the terminal
  physical surface, enclosing ordinary-fill surfaces, and every active
  lattice DDA level.
- Equal-distance physical and synthetic events retain existing tolerance and
  ordering rules.
- A retry may replace ownership only with an equally deep or deeper valid
  path; it must not silently climb to a container.

### Nested lattice paths

- Each lattice state is keyed by path-entry index, not cell index.
- Repeated occurrences of the same cell definition remain distinguishable.
- Crossing lattice level `n` invalidates/rebuilds only path entries below
  level `n`; entries above it remain validated.
- A child-lattice transition must retain the parent element transform.
- A parent-lattice transition invalidates all descendant lattice states.
- Element changes remain ownership transitions even if they resolve to the
  same terminal cell and material.

### First-visible and any-hit

- `FIRST_VISIBLE` skips void and material-zero intervals.
- Material filtering is applied only after interval ownership verification.
- `ANY_HIT` is the same acceptance policy with boolean-only publication.
- Physical entries preserve surface ID, primitive ID, and optional world-space
  normal.
- Synthetic lattice entries report surface ID `0`, invalid primitive, and zero
  normal.
- If `t_min` clips into an active interval, the reported entry has surface ID
  `-1`, invalid primitive, and zero normal.

### Batch behavior

- Result order equals input order.
- A lane may finish independently of other lanes.
- Failure or interruption never publishes a partially updated result object.
- `ray_count == 0` succeeds without allocating lane state.
- Output allocation remains proportional to requested fields.
- Per-lane traversal state is bounded by fixed packet width.

## Internal design

### 1. Resumable lane state

Introduce a private state type in `raycast.c`; expose it in `raycast.h` only if
`raycast_api.c` must allocate it directly. Prefer an opaque internal batch
entry point so detailed state remains owned by the traversal module.

Conceptual structure:

```c
typedef struct {
    alea_ray_t ray;
    double t_current;
    double t_max;
    double visible_t_min;

    int previous_cell_index;
    int previous_surface_id;
    int material_filter;

    alea_hier_ray_path_t path;
    alea_ray_boundary_event_t enter_event;
    double pending_interior_sample;

    lattice_level_state_t lattice_levels[ALEA_HIER_RAY_PATH_MAX];
    uint8_t lattice_level_count;

    uint8_t wants_normal;
    uint8_t done;
    uint8_t failed;
    uint8_t found;

    alea_ray_first_visible_result_t visible;
    alea_raycast_result_t* optional_trace;
} ray_query_lane_t;
```

Do not retain every transient variable from the current monolithic function.
Only state that survives from one completed interval to the next belongs in
the lane. Per-step cell candidates, transforms, probe values, and retry
snapshots remain stack-local inside the step function.

### 2. Interval-step result

The state machine should return an explicit outcome instead of relying on
multiple early returns:

```c
typedef enum {
    RAY_STEP_INTERVAL,
    RAY_STEP_ACCEPTED,
    RAY_STEP_FINISHED,
    RAY_STEP_ERROR
} ray_step_status_t;
```

One call performs at most one ownership interval:

```c
ray_step_status_t ray_query_lane_step(
    alea_system_t* sys,
    ray_query_lane_t* lane,
    ray_interval_t* out_interval);
```

`ray_interval_t` contains the verified owner, range, resolution flags, path,
and entry/exit boundary metadata needed by full-trace materialization or
early-stop acceptance. It must not own dynamic memory.

### 3. Per-lane hierarchy context

The existing hierarchical stepper already performs nested lattice descent and
suffix restart through `alea_hier_ray_path_t`. A packet lane must retain that
exact path between interval steps; it must not reconstruct hierarchy state
from a terminal cell ID or from the deepest lattice entry alone.

The first native batch implementation reuses the current path-aware boundary
selection unchanged. Any later lattice-level boundary cache is an optimization
only and must be keyed by path-entry index plus the concrete transform.

### 4. Scalar adapters

Keep current internal entry points, but implement them as adapters:

```text
initialize one lane
repeat lane_step
    full trace: materialize every interval/event
    first visible: apply acceptance and stop
    any hit: apply acceptance and stop with bool
finish lane
```

This maintains compatibility for current callers and makes scalar behavior the
oracle for the batch kernel without maintaining a second traversal body.

### 5. Native packet executor

Use a fixed-size packet per worker. The current compile-time default is four
lanes and remains private so it can be benchmarked and changed without ABI
impact.

```text
for each worker-owned contiguous ray range:
    fill packet with input rays
    while packet contains active lanes:
        step every active lane once
        publish accepted/finished lanes directly by original ray index
```

This is a native batch state machine even if the inner lane operation is
initially scalar C. It removes repeated scalar API dispatch and makes lane
state, scheduling, field computation, and output publication batch-aware.

Do not allocate one `alea_hier_ray_path_t` per input ray. Packet state is
worker-local and reused until the worker's range is exhausted.

### 6. Direct SoA publication

The existing batch output allocation and validation in `raycast_api.c` should
remain. After successful allocation, pass non-owning output spans to the
raycast batch kernel:

```c
typedef struct {
    uint8_t* found;
    double* t;
    int32_t* cell_ids;
    int32_t* material_ids;
    double* densities;
    int32_t* surface_ids;
    uint32_t* primitive_ids;
    uint8_t* resolution_flags;
    double* normals_xyz;
} first_visible_output_span_t;
```

Each lane writes only its original ray index, so no locks or atomics are
needed. The outer API publishes `next` into the caller's result only after all
workers succeed.

### 7. Error and interruption handling

Use one operation-local atomic failure/interruption flag. Workers check it:

- before starting and while advancing a packet;
- at a bounded step cadence;
- after any lane returns `RAY_STEP_ERROR`.

Per-ray status storage is unnecessary if the first error cancels the entire
operation and the failing lane records an operation-local error category.
The public wrapper remains responsible for setting a stable error detail after
the parallel region, because thread-local error strings written by a worker
are not reliably observable by the calling thread.

## Implementation phases

### Phase 0: lock the hierarchical nested contract

1. Add focused nested-lattice ray fixtures for forward and reverse rays.
2. Verify separately:
   - full hierarchical trace without hits;
   - full hierarchical trace with hits;
   - scalar `FIRST_VISIBLE`;
   - scalar `ANY_HIT`;
   - batch `FIRST_VISIBLE`;
   - batch `ANY_HIT`.
3. Use `FAST_FORWARD` hierarchical tracing as the nested-lattice full-trace
   oracle; document that the legacy `GLOBAL` lattice pipeline is not nested-
   transform-aware.
4. Establish expected owners using canonical point samples on each open
   interval.

Gate: hierarchical scalar paths are covered in both directions, and their
expected owner sequence is established from canonical point samples.

### Phase 1: extract the scalar lane state machine

1. Introduce `ray_query_lane_t`, `ray_interval_t`, and `ray_step_status_t`.
2. Move persistent locals from `raycast_cell_aware_impl()` into lane state.
3. Extract interval resolution and nearest-boundary selection into one step.
4. Preserve current midpoint verification/retry behavior.
5. Adapt existing scalar full-trace and early-stop functions.
6. Run all existing raycast, hierarchy, grid, and integration tests.

Gate: no observable output change for existing non-lattice and one-level
lattice cases; no performance claim yet.

### Phase 2: introduce the native packet kernel

1. Add a private batch-kernel input/output descriptor.
2. Allocate one fixed packet per OpenMP worker-owned contiguous range.
3. Initialize, step, and retire lanes by original ray index.
4. Write `FIRST_VISIBLE` results directly to SoA output.
5. Skip normal/primitive calculations unless requested.
6. Preserve transactional result publication in `raycast_api.c`.
7. Remove the per-ray call to
   `alea_raycast_hier_first_visible_nocache()` from the batch executor.

Gate: code inspection and tests prove that the batch loop no longer invokes
the scalar public/internal query entry point per ray.

### Phase 3: route `ANY_HIT` through the packet kernel

1. Reuse the same lane initialization and interval acceptance.
2. Publish only a boolean output.
3. Do not calculate normals, surface provenance, or density.
4. Remove lattice/non-lattice dispatch differences from the existing batch
   wrapper where semantics are identical.

Gate: any-hit output equals first-visible `found` for matching filters/ranges,
including nested lattices.

### Phase 4: determinism, resource, and failure tests

1. Compare scalar and batch outputs for packet-boundary ray counts:
   `0, 1, W-1, W, W+1, 2W+1`. **Complete for `W = 4`:** both early-stop APIs
   cover `0, 1, 3, 4, 5, 9` and reuse one result object across calls.
2. Repeat with one thread and multiple OpenMP thread counts. **Complete for
   focused raycast and MCNP integration suites** at one and four threads.
   Multi-ray packet parity now covers rectangular, hexagonal, and nested
   lattices, including a translated nested-lattice placement, clipped starts,
   filters, reverse directions, and void misses.
3. Interrupt a large batch and verify no partial result publication.
   **Partially complete:** pre-start interruption preserves the previously
   published result; deterministic mid-flight interruption coverage remains.
4. Force output-byte-budget rejection before traversal. **Complete.**
5. Reuse a result after success, failure, empty input, and interruption.
   **Complete for early-stop results:** success, byte-budget failure, empty
   input, and pre-start interruption are covered.
6. Verify optional fields remain unallocated when not requested. **Complete
   for existing API allocation coverage.**
7. Run sanitizer builds where available. **Complete for focused ASan/UBSan
   runs:** raycast unit and MCNP integration suites pass. LeakSanitizer cannot
   run in the hosted ptrace environment, so leak detection remains a CI/local
   follow-up.

Gate: deterministic byte-for-byte integer/ID fields and tolerance-equivalent
floating fields under all tested packet/thread configurations.

### Phase 5: benchmark and tune

**Initial baseline complete:** `perf_first_visible_native_packet_vs_scalar_20_shells`
compares the prepared scalar early-stop loop and native packet batch call with
identical rays and fields. On the recorded four-thread development run, the
scalar loop took 3.26 us/ray and the native four-lane packet path took 1.02
us/ray, with 3,114 accepted rays in each case. This is an informational
measurement, not a timing assertion.

The same benchmark shape now covers the rectangular MCNP lattice fixture. On
the recorded four-thread development run, its scalar loop took 5.06 us/ray
and the native packet path took 1.60 us/ray, with all 10,000 answers
matching. Nested-lattice, hexagonal-lattice, filtered-late-hit, and clipped
range workloads remain useful benchmark extensions.

`ANY_HIT` is also measured against its scalar loop on the shell workload. The
recorded four-thread run took 3.20 us/ray scalar and 0.96 us/ray native, with
3,114 matching hits. Its boolean-only result publishes 10,000 bytes for the
10,000-ray request.

The packet-width macro is private and build-overridable for measurement. A
separate optimized four-thread run measured widths `4`, `8`, and `16` at
approximately `0.74`, `0.77`, and `0.73` us/ray respectively on the same
shell workload. This difference is within ordinary run-to-run noise and does
not cover lattice divergence, so the committed default remains `4` pending
lattice-specific measurements.

Measure separately:

- cache preparation;
- batch output allocation;
- traversal;
- optional normal computation;
- total batch call;
- scalar loop baseline using identical rays and requested fields.

Use workloads covering:

- non-lattice shallow geometry;
- transformed ordinary fills;
- rectangular lattice;
- hexagonal lattice;
- nested rectangular lattice;
- mostly-miss/void rays;
- early visible hits;
- material-filtered late hits;
- clipped `t_min` starts.

Benchmark packet widths `4`, `8`, and `16`. Select one default based on wall
time and peak resident memory. Do not add explicit SIMD until profiles show
that lane stepping, rather than hierarchy/CSG resolution, dominates.

Gate: native batch is not slower than the current OpenMP scalar-dispatch
baseline on representative lattice workloads, and memory remains bounded by
worker packet state plus final output arrays. The initial non-lattice shell
baseline meets the performance direction; lattice-specific workload coverage
and width comparison remain.

## Test matrix

| Dimension | Required cases |
|---|---|
| Geometry | non-lattice, rect, hex, transformed lattice, nested rect |
| Direction | forward, reverse, oblique where stable fixture exists |
| Start | outside, inside terminal, inside outer lattice, inside inner lattice |
| Range | `t_min=0`, clipped `t_min`, finite `t_max`, no accepted interval |
| Policy | full trace, first-visible, any-hit |
| Filter | none, first material, later material, absent material |
| Entry | physical, synthetic DDA, clipped/no entry |
| Fields | IDs only, density, surfaces, primitive, normal, flags |
| Batch size | empty, singleton, packet edge, multiple packets |
| Parallel | one thread, several threads, repeated deterministic runs |

For nested full-trace tests, do not use forward/reverse agreement alone as the
oracle. Sample an interior point of every returned open interval with the
canonical deepest-path point resolver and compare cell, material, placement
path, and lattice element origins.

## Expected codebase impact

### `src/raycast/raycast.c`

Largest change. The current monolithic stepper becomes state initialization,
one-interval stepping, policy adapters, and a packet driver. Existing shared
lattice boundary helpers remain canonical.

### `src/raycast/raycast.h`

Add only private declarations needed by `raycast_api.c`. Avoid exposing lane
layout if the API module can call one opaque native batch entry point.

### `src/raycast/raycast_api.c`

Keep validation, requested-field allocation, byte budgets, and transactional
publication. Replace OpenMP scalar dispatch with the native batch kernel.

### `src/core/alea_spatial_hier.*`

Expected small changes: exact entry-indexed suffix restart or path-to-lattice-
stack helpers. Ownership resolution remains in the spatial core. Do not move
packet scheduling into this module.

### Renderer and slice consumers

No required source changes. They continue through existing query APIs and
benefit once those adapters use the native kernel.

## Risks and mitigations

### State-machine drift

Risk: scalar and batch paths acquire different semantics.

Mitigation: both consume the same `ray_query_lane_step()` implementation;
policy adapters only decide whether to materialize or stop.

### Excessive lane-state size

Risk: a full 64-entry hierarchy path per lane pressures cache.

Mitigation: fixed worker-local packets, benchmarked width, no per-ray state
array. A later compact path representation requires separate evidence.

### Packet divergence

Risk: rays take different numbers of steps and reduce SIMD efficiency.

Mitigation: retire lanes independently. The current contiguous-packet
implementation deliberately does not refill a completed lane; evaluate
refilling only after benchmark evidence. Treat the first implementation as
packet scheduling and direct materialization, not guaranteed vector SIMD.

### Nested path corruption

Risk: restarting by cell definition confuses repeated recursive occurrences.

Mitigation: retain the existing path-aware hierarchical restart logic inside
each lane; never identify a restart target by cell definition alone.

### OpenMP error propagation

Risk: worker TLS error details are lost or results are partially published.

Mitigation: operation-local atomic status plus caller-thread error reporting;
publish the prepared `next` result only after successful completion.

### Performance regression for small batches

Risk: packet setup costs more than a direct scalar call for one or two rays.

Mitigation: benchmark a small-batch threshold. A scalar adapter may remain for
tiny batches if it uses the same lane state machine and does not create a
separate semantic path.

## Suggested commit sequence

1. `Cover nested hierarchical ray policy parity`
2. `Refactor hierarchical ray traversal into lane steps`
3. `Add native first-visible packet traversal`
4. `Route batch any-hit through packet traversal`
5. `Cover packet boundaries and threaded determinism`
6. `Document and benchmark native lattice batch traversal`

Every commit after the initial failing-test characterization must leave the
normal test suite green. If the repository policy does not allow a committed
expected-failure test, keep the reproducer in the same commit as its fix.

## Completion criteria

The program is complete when all of the following are true:

- nested lattice full traces succeed in forward and reverse directions;
- full-trace intervals match canonical point/path ownership samples;
- scalar lattice `FIRST_VISIBLE`/`ANY_HIT` match full-trace projection;
- batch lattice `FIRST_VISIBLE`/`ANY_HIT` match scalar outputs;
- the first-visible and any-hit batch loops no longer call scalar query entry
  points per ray;
- batch state memory is bounded by workers and packet width;
- synthetic, physical, clipped, filtered, transformed, and nested semantics
  are covered;
- results are deterministic across packet edges and thread counts;
- interruption and failure preserve transactional output publication;
- benchmarks show no representative lattice regression against the current
  OpenMP scalar-dispatch baseline;
- related architecture-plan status text is updated to describe the completed
  implementation accurately.
