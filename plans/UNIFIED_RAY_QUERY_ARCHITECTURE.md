# Unified Ray Query Architecture Plan

**Status:** In progress — Phase 1 landed; later phases remain design work  
**Owning module:** `raycast`  
**Primary consumers:** 3D renderer, slice boundary provenance, compact slice
validation, transport and analysis APIs  
**Related plans:**

- `plans/NATIVE_COMPACT_RAYCAST_BATCH.md`
- `plans/NATIVE_BIDIRECTIONAL_RAY_SLICE_VALIDATION.md`
- `plans/SURFACE_BOUNDARY_PROVENANCE_LABELS.md`
- `plans/PARALLEL_SURFACE_BOUNDARY_MAP.md`

## Summary

Libalea currently has several ray paths optimized for different historical
uses: global surface collection, cell-aware traversal, hierarchical segment
stepping, hit-producing hierarchical tracing, compact batch tracing, first-cell
queries, and occlusion tests. Consumers also differ substantially:

- solid 3D rendering needs the first visible non-void hit and its normal;
- shadow rays need only an any-hit answer;
- X-ray rendering and transport need ordered material segments;
- surface boundary provenance needs ordered ownership-changing boundary events;
- directional validation needs comparable forward/reverse ownership intervals;
- diagnostic and analysis users may need hierarchy paths and occurrence keys.

The goal is not to replace these algorithms with one monolithic function. The
goal is to define one semantic query layer over shared traversal primitives,
with explicit stop policy, requested fields, reusable storage, and compact
batch execution. Existing specialized traversal code should be reused and
converged incrementally behind that contract.

Every migration is gated by output-equivalence tests and benchmarks. Existing
public APIs remain available until their replacements are proven and adopted.

## Current observations

### Non-lattice 3D rendering

The renderer already uses a sound execution pattern:

- OpenMP tile parallelism;
- one thread-local `alea_raycast_result_t` per worker;
- prebuilt query caches;
- internal no-cache hierarchical tracing;
- a specialized early-out occlusion query for shadows.

Solid rendering nevertheless builds a full hit/segment trace and then selects
the first visible segment. A first-visible query may avoid traversal and
materialization behind the winning hit.

### Lattice 3D rendering and boundary provenance

These paths use public `alea_raycast()`. That function frees result vectors on
every call. Calling it with a nominally reusable result therefore still causes
repeated allocation. A private canonical reusable entry point can benefit both
consumers immediately without changing public behavior.

### Compact batch and directional validation

`alea_raycast_hier_batch()` already provides:

- parallel ray orchestration;
- compact CSR segment storage;
- requested optional fields;
- hierarchy/occurrence projection;
- memory and segment budgets.

`alea_validate_ray_slice_compact()` builds on it for U-directed forward/reverse
ownership comparison. Its result is a diagnostic interval product, not a full
surface-provenance source. The underlying compact traces are reusable
infrastructure and should become shareable without conflating validation with
provenance.

### Boundary identity

Surface provenance is currently inferred from segment enter/exit fields. The
architecture needs a documented distinction between:

1. every mathematical surface intersection;
2. a physical surface event retained by traversal;
3. an ownership-changing boundary event;
4. a synthetic lattice/DDA transition;
5. the surface chosen for reporting when equivalent/coincident primitives are
   present.

Without this distinction, consumers can silently disagree about which surface
caused a transition.

## Decisions

1. Keep traversal algorithms separate from result materialization and stop
   policy.
2. Introduce an internal query descriptor with explicit query kind and field
   mask.
3. Make reusable result storage the default for internal high-frequency calls.
4. Define boundary-event semantics before migrating surface provenance.
5. Extend compact batch execution around query intent rather than adding
   viewer-specific batch APIs.
6. Keep ordinary rendering, validation, and provenance result types separate;
   allow them to share traces or traversal infrastructure.
7. Preserve existing public APIs through adapters during migration.
8. Do not replace working consumer paths until equivalence and performance are
   demonstrated for lattice and non-lattice models.
9. Treat interruption, budgets, deterministic ordering, and transactional
   result publication as core query semantics.

## Goals

1. Eliminate repeated ray-result allocation in high-frequency internal loops.
2. Allow traversal to stop as soon as a query's answer is known.
3. Avoid computing normals, paths, surface IDs, or complete segments when a
   consumer did not request them.
4. Provide exact, documented boundary-event semantics for provenance.
5. Reuse compact parallel ray execution across slice and rendering workloads.
6. Preserve deterministic results across scalar/batch execution and thread
   counts.
7. Unify lattice and non-lattice observable semantics even when their internal
   traversal remains different.
8. Improve 3D solid, X-ray, lattice, boundary-map, and validation performance
   without weakening correctness.

## Non-goals

- Do not merge every traversal implementation in the first phase.
- Do not expose internal event/vector layouts as ABI-stable public structs.
- Do not make the geometry validator a rendering or provenance API.
- Do not require full-segment batch materialization for first-hit rendering.
- Do not change public raycast memory-retention behavior without separate API
  review.
- Do not treat forward/reverse agreement as proof of valid geometry.
- Do not migrate physics/transport consumers before their precision and event
  requirements are explicitly audited.

## Phase 0: semantic and performance audit

Create a maintained matrix for every existing ray entry point:

| Property | Questions |
|---|---|
| Traversal | global surface BVH, cell-aware, hierarchical, lattice-aware? |
| Output | hits, segments, first cell, any hit, normals, paths? |
| Surface semantics | complete intersections, ownership boundaries, best effort? |
| Lattice | physical and synthetic boundaries preserved? |
| Allocation | frees, clears, reuses, compact result? |
| Cache | prepares caches, requires prepared caches, no-cache internal? |
| Stop behavior | full ray, first cell, first visible, any hit? |
| Thread safety | safe after cache preparation? mutable TLS/global counters? |
| Budgets | bounded segments, paths, bytes, interruption? |

Cover at least:

- `alea_raycast()`;
- `alea_raycast_cell_aware()`;
- hierarchical segment and hit-producing internal paths;
- `alea_raycast_hier_fast_segments()`;
- `alea_raycast_hier_batch()`;
- `alea_trace_ray_slice_compact()`;
- `alea_ray_first_cell()`;
- `alea_ray_is_occluded()`;
- renderer-specific no-cache calls;
- geometry-validator ray paths.

Instrument representative workloads before changing behavior:

- rays traced and average crossings;
- surfaces and cells tested;
- result allocations/reallocations and retained capacity;
- segments/hits/events materialized;
- early-stop distance and work avoided;
- cache misses/fallbacks;
- time split between traversal, classification, and assembly;
- peak and transient memory.

## Phase 1: internal reusable canonical tracing

**Implementation status (2026-07-22):** complete.  The internal entry point is
`alea_raycast_global_reuse_nocache()`.  It is equivalence-tested against the
public global raycast on a rectangular lattice, including synthetic DDA
crossings and a repeated call on the same result buffer.  It now serves the
boundary-map short-edge trace and lattice solid/X-ray rendering.

Add a private canonical trace entry point whose observable output matches
`alea_raycast()` but which calls `alea_raycast_result_clear()` and retains
vector capacity:

```c
int alea_raycast_global_reuse_nocache(
    alea_system_t* sys,
    const alea_ray_t* ray,
    double t_max,
    alea_raycast_result_t* result);
```

The exact name may change. Required behavior:

- caches are prepared by the caller;
- the result was initialized by the caller;
- prior logical contents are cleared;
- allocated vector capacity is retained;
- failure leaves the result empty and reusable;
- lattice and non-lattice output matches canonical public raycast behavior.

Initial migrations:

1. slice boundary-map short-edge fallback;
2. lattice solid rendering;
3. lattice X-ray rendering.

Do not change public `alea_raycast()` allocation behavior in this phase.

## Phase 2: query descriptor and result policy

Define an internal query descriptor:

```c
typedef enum {
    ALEA_RAY_QUERY_ANY_HIT,
    ALEA_RAY_QUERY_FIRST_CELL,
    ALEA_RAY_QUERY_FIRST_VISIBLE,
    ALEA_RAY_QUERY_SEGMENTS,
    ALEA_RAY_QUERY_BOUNDARY_EVENTS
} alea_ray_query_kind_t;

typedef struct {
    alea_ray_query_kind_t kind;
    uint32_t fields;
    double t_min;
    double t_max;
    int material_filter;
    uint64_t max_events;
    uint64_t max_output_bytes;
    const void* clip_policy; /* internal/optional; avoid renderer ABI here */
} alea_ray_query_t;
```

Suggested requested fields:

```text
cell ID
material ID
density
surface ID
surface normal
resolution/coverage flags
projected owner
occurrence key
full hierarchy path
primitive identity (internal diagnostics only)
```

The descriptor is initially internal. Public exposure should wait until at
least two consumers use it successfully and its error/budget semantics are
stable.

Implement result writers/policies rather than one universal result object:

- Boolean any-hit writer;
- first-hit/first-visible record;
- rich hit/segment vectors;
- boundary-event vector;
- compact CSR batch writer.

Traversal calls the writer for ordered events. The writer decides whether to
retain data and whether the answer is complete, allowing traversal to stop.

## Phase 3: boundary-event contract

Define the internal event before optimizing provenance:

```c
typedef enum {
    ALEA_RAY_EVENT_PHYSICAL,
    ALEA_RAY_EVENT_SYNTHETIC_LATTICE,
    ALEA_RAY_EVENT_UNRESOLVED
} alea_ray_event_kind_t;

typedef struct {
    double t;
    alea_ray_event_kind_t kind;
    int surface_id;          /* -1 none, 0 synthetic, >0 physical */
    uint32_t primitive_id;   /* canonical match identity, optional */
    int cell_before;
    int cell_after;
    int material_before;
    int material_after;
    uint32_t resolution_flags;
    double normal[3];        /* present only when requested */
} alea_ray_boundary_event_t;
```

Resolve and document:

1. whether unchanged ownership intersections are emitted;
2. how coincident events are grouped;
3. whether multiple MC surface IDs sharing a primitive are retained;
4. ordering of physical and synthetic hits at equal `t`;
5. surface ID selected at Boolean/complement boundaries;
6. behavior at gaps, overlaps, undefined fills, and exterior void;
7. transform/lattice occurrence identity;
8. tolerance and deduplication rules;
9. forward/reverse normalization and normal orientation.

Recommended default: `BOUNDARY_EVENTS` emits ownership-changing events plus
synthetic lattice transitions. A separate diagnostic flag may request all
mathematical/physical intersections when the traversal can provide them.

Add canonical event fixtures before any consumer migration.

## Phase 4: specialized query policies

### Any hit

Keep and align the existing occlusion path. Stop on the first qualifying
occluder. Do not materialize segments or normals.

### First cell

Return the first resolved cell matching the query filter. Preserve existing
void/material conventions.

### First visible

For solid 3D rendering, stop at the first non-void visible segment after the
query's `t_min`. Return only:

- hit distance;
- cell/material;
- surface ID;
- normal;
- optional occurrence/resolution flags.

Clipping must be expressed as a traversal-independent acceptance policy or as
one/more ray intervals. Do not put renderer clip-plane structs into the public
raycast API.

### Segments

Preserve complete ordered intervals for X-ray, transport, and analysis. Avoid
normals and paths unless requested.

### Boundary events

Return the documented event stream needed by slice provenance and diagnostics.

## Phase 5: compact batch executor

Generalize compact batch orchestration to accept query intent and requested
fields while preserving current batch APIs through adapters.

Required capabilities:

- per-ray origins and directions;
- per-ray or grouped `t_min`/`t_max` ranges;
- query-kind-specific compact outputs;
- deterministic input-ray order;
- per-thread scratch and no shared hot-loop append;
- transactional CSR/SoA publication;
- segment/event/path/byte budgets;
- interruption;
- reusable result capacity when safe;
- no nested OpenMP oversubscription.

Per-ray ranges are important for arbitrary short edges, clipped primary rays,
and mixed workloads. If adding per-ray `t_max` to the existing public batch API
would destabilize it, introduce an internal executor first and retain the
current scalar-`t_max` adapter.

## Phase 6: shared directional slice traces

Create a shareable slice trace cache beneath validation and provenance:

```text
U forward compact trace
U reverse compact trace
V forward compact trace
V reverse compact trace
```

The cache records:

- system identity and geometry generation;
- complete slice view and sampling dimensions;
- direction/orientation;
- requested fields;
- projected depth and occurrence policy;
- trace budgets and completeness.

Consumers:

- directional validation compares forward/reverse ownership;
- boundary provenance maps events to right/down grid edges;
- renderers may reuse compatible forward ownership traces;
- label-only verification can query selected edges without a full map.

Do not make boundary provenance consume validator disagreement intervals. Both
consume the shared compact traces. Extend reverse validation traces to request
surface fields only when a compatible downstream consumer needs them.

Map right edges from U traces and down edges from V traces. If trace semantics
are incomplete or inconsistent for an edge, fall back to canonical reusable
short-edge boundary-event tracing.

## Phase 7: 3D renderer migration

### Solid mode

Benchmark the current per-thread full trace against:

1. scalar reusable first-visible queries;
2. tile-batched first-visible queries.

Retain current shading and shadow behavior. Request normals only for accepted
hits. Validate cell/material/depth/normal images bitwise or within documented
floating tolerance.

### Shadow mode

Keep specialized any-hit behavior. Migrate only if the shared query layer
matches or improves early termination.

### X-ray mode

Evaluate compact all-segment batches per tile. This workload naturally
consumes every material interval and is a stronger batch candidate than solid
shading.

### Lattice models

Prioritize reusable storage first. Lattice semantics, synthetic DDA events,
fill occurrence, and normals must match the current canonical path before
first-visible or batch migration.

## Phase 8: public API consolidation

Only after internal migrations stabilize:

1. document a small public query-options structure;
2. expose first-visible and boundary-event APIs if external demand exists;
3. keep existing APIs as compatibility wrappers;
4. mark APIs deprecated only after equivalent replacements exist in C and Lua;
5. keep opaque results and accessor-based ownership;
6. use `struct_size` parsing that accepts older prefixes safely.

Do not expose renderer-specific policies or internal primitive payloads by
default.

## Correctness matrix

Every query kind and migration must cover:

- simple planes, spheres, cylinders, and boxes;
- Boolean union/intersection/complement;
- void and graveyard transitions;
- exact and near-coincident surfaces;
- duplicate/equivalent primitives with distinct MC IDs;
- transformed cells and surfaces;
- nested fills and projected-depth ownership;
- rectangular and hexagonal lattices;
- synthetic DDA boundaries;
- undefined fills, gaps, overlaps, and ambiguous ownership;
- rays starting on a boundary;
- grazing/tangent rays;
- zero, finite, and effectively infinite `t_max`;
- forward/reverse normalization;
- deterministic event ordering;
- interruption and every budget failure mode.

Compare:

- scalar versus batch;
- reusable versus current public tracing;
- lattice versus canonical global behavior;
- forward versus normalized reverse;
- old versus migrated renderer images;
- old versus migrated boundary maps;
- thread counts and OpenMP schedules.

## Performance gates

Do not migrate a consumer on intuition alone. Record:

- rays/second and pixels/second;
- median and tail ray cost;
- traversal steps and tested surfaces/cells;
- result allocations and bytes retained;
- hits/segments/events materialized;
- work avoided by early termination;
- compact assembly time;
- peak memory;
- scaling across 1, 2, 4, 8, and available maximum threads.

Suggested acceptance gates:

1. No correctness regression in the matrix above.
2. No material slowdown on representative non-lattice scenes.
3. Clear allocation reduction for boundary and lattice render loops.
4. First-visible mode demonstrates meaningful work reduction before replacing
   solid rendering.
5. Batch mode remains within an explicit transient-memory budget.

## Error and ownership semantics

- Query functions return distinct invalid-argument, invalid-state,
  out-of-memory, budget, interruption, and unsupported-mode errors.
- Failed rich/compact results publish transactionally.
- Reusable internal results remain empty and reusable after failure.
- No silent truncation of hits, segments, events, or paths.
- Thread-local failures are collected deterministically.
- Public pointer lifetimes remain tied to opaque result reuse/destruction.

## Source layout

Retain existing files initially. Split only when responsibilities become clear:

```text
src/raycast/
  raycast.c                 traversal implementations and legacy adapters
  raycast.h                 internal query/event contracts
  raycast_api.c             public/compact result APIs
  ray_query.c               query dispatch and result policies (new, later)
  ray_batch.c               generalized batch executor (new, later)
  ray_events.c              event grouping/dedup semantics (new, later)

src/slice/
  slice_api.c
  slice_directional_trace.c shared U/V cache orchestration (new, later)

src/render/
  render3d.c                consumer; no duplicate traversal logic
```

Add all new source files to `Makefile` and `Makefile.msvc` in the same change.

## Rollout sequence

1. Complete the semantic/performance audit.
2. Add reusable canonical internal tracing and migrate allocation-heavy loops.
3. Define and test boundary-event semantics.
4. Introduce internal query descriptors and result writers.
5. Generalize compact batch execution.
6. Add shared U/V directional slice traces.
7. Migrate boundary provenance with canonical fallback.
8. Benchmark and migrate solid first-visible rendering.
9. Benchmark and migrate X-ray batch rendering.
10. Consolidate public APIs only after internal consumers prove the design.

Each step must be independently mergeable and leave existing public behavior
working. Do not combine the semantic event change, batch rewrite, renderer
migration, and public API redesign into one patch series.

## Open questions

1. Should `FIRST_VISIBLE` accept a generic interval/acceptance callback, or
   should the renderer precompute acceptable ray intervals from clips?
2. Is the canonical provenance event ownership-changing only, or must a public
   diagnostic mode expose all mathematical intersections?
3. Can the hierarchical fast stepper guarantee physical surface IDs for every
   lattice/transform transition, or which cases require canonical fallback?
4. Should compact batch results retain capacity across reuse, and how is peak
   retained memory bounded?
5. Is a per-ray range array worth adding to the public batch API, or should it
   remain internal until external demand exists?
6. Which trace-cache layer owns reverse and V-direction results: `raycast`,
   `slice`, or `geo_validator`? Preferred answer: `slice` owns view sampling,
   `raycast` owns rays/results, and `geo_validator` consumes them.
