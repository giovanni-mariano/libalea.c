# Plan: native ray-slice rasterization

**Date:** 2026-08-08

**Status:** implemented — standalone and fused APIs, deterministic row
parallelism, native tests, and a focused C benchmark are available.

**Target branch:** `feature/unified-ray-query`

**Primary implementation areas:** `include/alea_slice.h`, `src/slice/`,
`src/raycast/raycast_api.c`, `src/raycast/raycast.h`, `tests/unit/`, and a
focused C benchmark/example.

**Related plans:**

- `plans/UNIFIED_RAY_QUERY_ARCHITECTURE.md`
- `plans/PLAN_NATIVE_BATCH_LATTICE_DDA.md`

## Executive decision

Add a renderer-independent ray-slice rasterizer to libalea. It will consume the
existing compact row-oriented ray-slice result and fill caller-owned,
contiguous output buffers. Add a fused libalea convenience entry point that
traces and rasterizes a slice for callers that do not need the intermediate
compact segments.

This plan ends at the libalea C API. Bindings and viewer migrations are
downstream work and are not part of its implementation or completion gates.
The continuous-segment API remains available for validation, diagnostics,
provenance, and frontends that need exact intervals.

## Motivation

The native compact slice path currently stops at continuous CSR intervals.
Every frontend that wants a pixel grid must independently implement the same
interval-to-pixel conversion, allocate and initialize its own fields, and fill
each interval run. That duplicates a numerically sensitive ownership rule and
encourages per-interval dispatch and unnecessary field materialization outside
libalea.

Native run filling is a natural extension of compact row tracing: rows are
independent, interval ownership is already resolved, and each interval maps to
a contiguous range of horizontal pixel centers. A fused adapter also avoids
publishing or copying compact arrays when the caller only needs a raster,
although its first implementation still uses a temporary compact result
internally.

## Goals

1. Move compact-interval-to-pixel conversion into libalea.
2. Freeze and preserve the documented pixel-center ownership rule.
3. Fill only caller-requested raster fields.
4. Keep every output buffer caller-owned and independently reusable.
5. Reject generic ray batches or mismatched slice results instead of
   silently interpreting ray distances as view-U coordinates.
6. Keep the existing compact ray-slice and validation APIs intact.
7. Establish native correctness and performance evidence.
8. Make the rasterizer useful to C, Lua, CLI, and future frontends without a
   dependency on any language runtime or array library.

## Non-goals

- Do not produce RGBA pixels in libalea.
- Do not move frontend color, contour, label, filter, cache, or presentation
  policy into C.
- Do not replace grid point classification with ray tracing.
- Do not change ray ownership, lattice, projected-depth, or clipping
  semantics.
- Do not change compact validation or provenance result formats initially.
- Do not introduce a system-owned or process-global reusable framebuffer.
- Do not retain caller output pointers in libalea.
- Do not add particle importance to the generic raster API. MCNP neutron,
  photon, and electron importances are model-specific metadata, are not stored
  in `alea_system_t`, and are not compact ray fields.
- Do not implement or modify Python bindings or frontend consumers in this
  plan.
- Do not add timing thresholds to ordinary correctness tests.

## Current behavior and ownership baseline

`alea_trace_ray_slice_compact()` fills an opaque, libalea-owned
`alea_raycast_batch_result_t`. Its row offsets form a CSR index and its
`t_enter`/`t_exit` values are rewritten from generic ray distances into clipped
view-U coordinates. The result also carries private provenance containing the
slice view, row count, projected depth, system identity, and geometry
generation.

The raster API will not transfer ownership of that result. Standalone
rasterization borrows it for the duration of the call. Fused rasterization
creates and destroys a temporary compact result internally. In both cases the
caller allocates and owns every output buffer.

## Required semantic invariants

### Pixel sampling

- Pixel ownership is evaluated at horizontal pixel centers.
- A compact interval `[u_enter, u_exit)` covers exactly the pixel indices
  defined by this normative formula:

  ```text
  lo = ceil((u_enter - u_min) / du - 0.5)
  hi = ceil((u_exit  - u_min) / du - 0.5)
  clamp lo and hi to [0, nu]
  fill [lo, hi)
  ```

- Rows correspond to increasing V and use the centered V locations generated
  by `alea_trace_ray_slice_compact()`.
- `row_count` must equal the raster height `nv`.
- Degenerate and fully clipped intervals write no pixels.
- Adjacent interval behavior at coincident or tolerance-adjusted endpoints
  must remain deterministic.
- Intervals within a row must be ordered by nondecreasing `u_enter`. Small
  overlaps are permitted; intervals are applied in CSR source order, so a
  later interval deterministically owns pixels also covered by an earlier one.

### Field values

- Cell, material, universe, and fill-universe values come from the compact
  segment's projected owner at the requested hierarchy depth.
- Density is the resolved leaf segment density regardless of projected depth.
  The current hierarchy path and projected-owner arrays do not carry density;
  adding projected density requires a separate raycast contract change.
- The documented defaults are cell `-1`, material `0`, universe `-1`, fill
  universe `-1`, density `0.0`, and resolution flags `0`.
- Material zero remains data, not presentation transparency.
- Resolution flags are copied as native `ALEA_RESOLVE_*` bits. Mapping those
  flags to frontend-specific error codes is outside libalea.
- Unrequested fields may have null output pointers and must incur no storage
  or fill cost.
- Importance is not a raster field. In particular, a single importance value
  cannot represent MCNP's particle-specific importance metadata.

### Ownership and lifetime

- Output buffers are caller-owned.
- Libalea borrows output pointers only for the duration of the call.
- Libalea never allocates, frees, or retains caller raster buffers.
- The standalone function completes all validation before its first output
  write. Validation failure leaves every output buffer unchanged.
- Once filling begins, the initial implementation has no fallible operation or
  interruption point; successful calls initialize every requested pixel.

### Concurrency

- Different rows may be filled concurrently because their output regions do
  not overlap.
- Rasterization must be deterministic for one and multiple OpenMP threads.
- Output-buffer independence does not by itself establish that simultaneous
  calls on the same mutable system are safe; existing system/cache concurrency
  rules remain authoritative.

## Public libalea API design

Place the public types and functions in `include/alea_slice.h`. Use fixed-width
types for public raster fields and explicit dimensions. The final names may be
adjusted to match neighboring APIs, but the contract should resemble the
following.

### Requested fields

```c
typedef uint32_t alea_slice_raster_fields_t;

#define ALEA_SLICE_RASTER_CELL_ID          (1u << 0)
#define ALEA_SLICE_RASTER_MATERIAL_ID      (1u << 1)
#define ALEA_SLICE_RASTER_UNIVERSE_ID      (1u << 2)
#define ALEA_SLICE_RASTER_FILL_UNIVERSE    (1u << 3)
#define ALEA_SLICE_RASTER_DENSITY          (1u << 4)
#define ALEA_SLICE_RASTER_RESOLUTION_FLAGS (1u << 5)
```

Reject unknown bits so future fields cannot be silently ignored by older
libraries.

### Output descriptor

```c
typedef struct {
    size_t struct_size;
    size_t nu;
    size_t nv;
    uint32_t fields;

    int32_t* cell_ids;
    int32_t* material_ids;
    int32_t* universe_ids;
    int32_t* fill_universe_ids;
    double* densities;
    uint8_t* resolution_flags;
} alea_slice_raster_t;
```

Phase 1 should require tightly packed, row-major buffers of `nu * nv`
elements. Require `nu > 0`, `nv > 0`, and at least one known field bit. Every
requested pointer must identify a writable array of the corresponding element
type and length; requested output ranges must not overlap. Unrequested pointers
are ignored. Add byte strides only if a real non-contiguous consumer requires
them; unnecessary stride support complicates overflow validation and parallel
filling.

Provide an initializer that sets `struct_size` and clears all fields:

```c
void alea_slice_raster_init(alea_slice_raster_t* raster);
```

### Standalone rasterizer

```c
int alea_rasterize_ray_slice_compact(
    const alea_slice_view_t* view,
    const alea_raycast_batch_result_t* segments,
    alea_slice_raster_t* output);
```

This function:

- validates the view, compact-slice provenance, row count, CSR offsets,
  interval finiteness/order/range, positive output dimensions, nonempty known
  field mask, pointer presence/non-overlap, source-field availability, and all
  element/byte/address-range arithmetic before writing any output;
- initializes every requested output buffer to its documented default;
- fills pixel runs from the compact result;
- does not mutate or take ownership of the compact result;
- does not allocate proportional-to-pixel temporary storage.

The compact result must expose every source field requested by the raster
descriptor. Missing compact fields are an error, not an implicit fallback.
Generic `alea_raycast_hier_batch()` results must be rejected: their interval
coordinates are ray distances, not view-U coordinates.

Add a narrow internal raycast helper that checks the result's private compact
slice provenance against `view` and `output->nv` and returns its recorded
projected depth. Rasterization does not need a live system or a matching current
geometry generation because it consumes an immutable result snapshot; it does
need exact view-coordinate, row-count, and projection provenance. Keep the
existing stronger system/generation check used for trace reuse unchanged.

### Fused trace-and-raster entry point

Do not reuse `alea_raycast_batch_options_t` directly: its caller-provided
`fields` member conflicts with the fused API's derived compact field mask, and
its `max_output_bytes` is specifically a compact-array limit. Use a dedicated,
extensible options structure:

```c
typedef struct {
    size_t struct_size;
    int projected_depth;             /* -1 = leaf; >= 0 = path depth */
    uint64_t max_segments;           /* 0 = unbounded */
    uint64_t max_trace_output_bytes; /* compact temporary only; 0 = unbounded */
} alea_slice_raster_options_t;

void alea_slice_raster_options_init(alea_slice_raster_options_t* options);
```

The initializer sets `projected_depth = -1` and all limits to zero. A null
options pointer has those same leaf/unbounded defaults. Validate `struct_size`
using the extension convention already used by neighboring public options
structures.

`max_path_entries` is intentionally absent. The raster path never publishes
`ALEA_RAY_BATCH_FULL_PATHS`; projected-owner capture uses internal hierarchy
paths but not the flattened full-path output governed by that limit.

```c
int alea_trace_ray_slice_raster(
    alea_system_t* sys,
    const alea_slice_view_t* view,
    const alea_slice_raster_options_t* options,
    alea_slice_raster_t* output);
```

This function uses `output->nv` as the ray row count and `output->nu` as the
horizontal raster width. It should:

1. derive the complete compact field mask from the raster fields and projected
   depth, ignoring no caller-controlled trace field mask because none exists;
2. copy the relevant limits into an internal `alea_raycast_batch_options_t`;
3. trace into a temporary compact result;
4. call the standalone rasterizer;
5. destroy the temporary result on every exit path.

At leaf depth, cell IDs use the mandatory base cell array and material IDs use
`ALEA_RAY_BATCH_MATERIAL`; these already describe the resolved leaf owner. At a
nonnegative projected depth, cell and material IDs require
`ALEA_RAY_BATCH_PROJECTED_OWNER`. Universe or fill-universe output requires
`ALEA_RAY_BATCH_PROJECTED_OWNER` at every depth because no base arrays exist for
those fields. Density uses `ALEA_RAY_BATCH_DENSITY` and always remains leaf
density. Resolution flags use `ALEA_RAY_BATCH_RESOLUTION_FLAGS`. The standalone
rasterizer applies the same selection using the projected depth recorded in
compact-slice provenance and rejects a missing required source field.

Initially reuse the existing compact tracer rather than creating a second
traversal or direct-to-pixel ownership algorithm. Direct run publication may
be considered later only after the reusable implementation is verified and
profiled.

### Budgets and interruption

- Existing compact limits (`max_segments`, `max_path_entries`, and
  `max_output_bytes`) retain their current meaning wherever the existing batch
  options are used.
- `max_trace_output_bytes` maps to the temporary compact result's existing
  `max_output_bytes`; it does not include caller-owned raster buffers.
- Libalea must overflow-check `nu * nv`, each requested field's byte count, and
  their sum before execution. The caller is responsible for applying any
  allocation budget before supplying buffers.
- Do not reinterpret the existing `max_output_bytes` as a combined
  compact-plus-raster budget.
- Preserve existing signal/interruption polling in compact tracing.
- The initial raster fill has no interruption polling so that validation
  failure is the only error after buffers are supplied and no partial output is
  observable. Reconsider this only together with an explicit partial-output or
  transactional-buffer contract.

## Native implementation layout

Create a focused slice implementation file, for example:

```text
src/slice/slice_ray_raster.c
```

Keep pixel-index conversion and field filling out of `raycast_api.c`. The
raycast module owns compact traversal, result accessors, and the narrow private
provenance predicate; the slice module owns view-to-raster projection.

Recommended internal structure:

1. descriptor, field-mask, pointer, and overflow validation;
2. compact-slice provenance and complete CSR/interval preflight;
3. requested-field default initialization;
4. interval-to-pixel-range conversion;
5. requested-field run fills;
6. fused options translation, trace adapter, and cleanup.

Use explicit typed fill loops or small private helpers. Avoid calling `memset`
for nonzero integer defaults such as `-1` unless the representation assumption
is already guaranteed and documented.

### Parallelization

Start with deterministic parallelism over rows if benchmark evidence supports
it:

```text
parallel for row in [0, nv):
    initialize row defaults
    for each compact interval in row:
        calculate [lo, hi)
        fill requested fields for that run
```

Combining initialization and filling per row improves locality and avoids two
full raster passes. Do not add nested OpenMP parallelism when the fused compact
trace has already entered a parallel region. The simplest first version may
use serial rasterization and add row parallelism as an independently measured
phase.

## Validation strategy

### Libalea unit tests

Add native tests covering:

- pure interval-to-pixel helper cases for one interval spanning the complete
  row, multiple adjacent intervals, and void gaps before, between, and after
  intervals;
- helper cases exactly on pixel centers and pixel edges, clipped at `u_min`
  and `u_max`, and subpixel intervals containing zero or one pixel center;
- one-column and one-row rasters;
- projected root, intermediate, and leaf depth;
- leaf-density behavior at every projected depth;
- material-zero intervals;
- undefined-fill/resolution flags;
- arbitrary slice planes;
- lattice and nested-lattice projected ownership;
- null unrequested output pointers;
- zero dimensions, an empty field mask, unknown field bits, null requested
  pointers, and overlapping requested output ranges;
- missing required compact source fields;
- generic non-slice batch rejection and mismatched view/row-count provenance;
- inconsistent CSR offsets and intervals through a private validation helper or
  focused fault-injection fixture;
- dimension, element-count, individual-field byte-count, and total-byte-count
  overflow;
- validation failures leave pre-seeded caller buffers byte-identical;
- deterministic equality at one and multiple OpenMP thread counts;
- compact input remains immutable and reusable after success and failure;
- fused output is byte-identical to explicit compact-trace-then-raster output;
- ASan/LSan coverage for temporary-result cleanup on every fused failure path.

Because `alea_raycast_batch_result_t` is opaque, public tests should obtain
compact inputs through `alea_trace_ray_slice_compact()`. Exercise adversarial
raw intervals and malformed CSR through small private helpers or explicit test
hooks rather than weakening the public result abstraction.

## Benchmark plan

Add a focused C benchmark, for example
`examples/c/ray_slice_raster_bench.c`, with separately measured stages:

```text
ray_compact_trace
native_ray_raster_from_compact
native_ray_trace_and_raster
```

Benchmark at several raster sizes and geometry complexities, including cases
with few long intervals and many short lattice intervals. Record:

- median and minimum wall time;
- trace time;
- raster time;
- compact temporary segment count and allocated output bytes;
- raster field bytes and peak resident or allocated bytes where practical;
- one-thread and representative multi-thread configurations.

The main acceptance evidence is:

1. exact raster parity;
2. bounded raster fill cost relative to compact tracing;
3. no regression in compact tracing or grid rendering;
4. measured benefit, if any, from row parallelism on representative sizes.

Do not claim downstream end-to-end speedups from the libalea microbenchmark.
The benchmark establishes native cost, allocation behavior, and scaling; each
frontend remains responsible for measuring its own integration.

## Phased implementation

### Phase 0: baseline and contract

- Add the focused C benchmark and record compact trace baselines.
- Freeze the pixel-center formula and documented defaults in focused native
  helper tests.
- Finalize fixed-width public fields, leaf-density semantics, provenance
  requirements, dedicated fused options, and compact-only byte-budget
  semantics.

**Gate:** documented C contract and reproducible native baselines.

### Phase 1: standalone libalea rasterizer

- Add public raster descriptor, field flags, initializer, and standalone
  raster function.
- Add the narrow internal compact-slice provenance predicate.
- Preflight all descriptor, source-field, CSR, interval, and overflow
  conditions before the first output write.
- Implement serial row filling from an existing compact result.
- Add native edge, overflow, field-mask, and parity tests.

**Gate:** native helper tests and exact raster results from real compact slice
fixtures; generic batches and mismatched provenance are rejected.

### Phase 2: fused libalea API

- Add `alea_trace_ray_slice_raster()` as an adapter over compact tracing plus
  the standalone rasterizer.
- Add `alea_slice_raster_options_t`; do not expose a caller-controlled compact
  field mask through the fused API.
- Derive minimal compact trace fields from requested raster fields.
- Preserve compact trace limits, interruption, transactional validation, and
  cleanup behavior.
- Test lattice/non-lattice parity against explicit trace-then-raster calls.

**Gate:** fused and two-step native results are byte-identical.

### Phase 3: demand-driven fields and parallel rasterization

- Verify that every raster field requests only the compact data necessary for
  its documented semantics.
- Benchmark serial versus OpenMP row filling.
- Add row parallelism only where it improves representative end-to-end cases.
- Ensure there is no harmful nested OpenMP behavior.

**Gate:** measured improvement with deterministic outputs and no excessive
small-render overhead.

### Phase 4: cleanup and documentation

- Consolidate duplicate field/default and interval-index mapping helpers.
- Document standalone versus fused use, caller ownership, exact defaults,
  density semantics, provenance rejection, and budget behavior.
- Add or update a C example showing allocation, invocation, and cleanup.
- Record final benchmark results and implementation status in this plan.
- Keep compact APIs as first-class diagnostic interfaces.

## Risks and mitigations

### Boundary pixels change ownership

The half-pixel/`ceil` formula is easy to alter through floating-point
rearrangement.

**Mitigation:** freeze adversarial endpoint fixtures first, keep the arithmetic
in one native helper, and compare its index ranges with the documented formula.

### Field semantics diverge from grid rendering

Ray projected-owner fields and point-grid fields may not represent fills,
void, or errors identically in every hierarchy case.

**Mitigation:** preserve current ray semantics in phase 1; do not redefine
them to resemble grid semantics without separate geometry analysis and tests.
Document density as leaf density even when owner IDs are projected.

### A generic batch is mistaken for a compact slice

Generic ray batches expose distances from their ray origins, while compact
slice batches expose clipped view-U coordinates. Both use the same opaque
result type and accessors.

**Mitigation:** require and validate private compact-slice provenance before
any output write. Test generic-batch and mismatched-view rejection explicitly.

### Fused call performs unnecessary compact work

The existing compact tracer may still generate surfaces, paths, or fields not
needed by the raster.

**Mitigation:** derive the compact batch field mask from requested raster
fields and profile before pursuing direct-to-raster traversal. Record that
`ALEA_RAY_BATCH_PROJECTED_OWNER` currently materializes its complete grouped
field set even if the raster requests only one projected owner field. Split
that group only as a separately justified raycast API optimization.

### Caller buffer size cannot be proven

The descriptor supplies pointers and dimensions but no independently
verifiable allocation capacities. Libalea can prove arithmetic safety, not
that a caller actually allocated the declared number of elements.

**Mitigation:** document required capacities precisely, validate all count and
byte arithmetic, never retain or free the pointers, and exercise correct and
deliberately invalid descriptors under sanitizers where possible.

### Failure after output writes

Interruption or late validation during raster filling could leave partially
initialized caller buffers.

**Mitigation:** perform a complete preflight before initialization and keep the
initial fill loop infallible and non-interruptible. Any future cancellation
support must explicitly revise the output-on-error contract or use temporary
transactional storage.

### Rasterization is not the dominant cost

Compact tracing may dominate the native call, and frontend presentation or
transport may dominate end-to-end rendering.

**Mitigation:** retain native stage-level timing, report the rasterizer's
measured contribution honestly, and make no downstream performance claim from
libalea-only evidence.

## Completion criteria

This program is complete when:

1. libalea exposes documented standalone and fused ray-slice raster APIs;
2. output buffers are caller-owned and never retained by libalea;
3. validation completes before output writes and failures preserve caller
   buffers;
4. generic batches and mismatched compact-slice provenance are rejected;
5. cell/material/universe/fill projection, leaf density, defaults, and
   resolution-flag semantics are documented and tested;
6. fused and standalone paths are byte-identical across ordinary, fill,
   lattice, nested-lattice, arbitrary-plane, and projected-depth fixtures;
7. compact and segment consumers remain supported unchanged;
8. overflow, cleanup, determinism, and memory-tool tests pass;
9. the focused C benchmark records trace, raster, fused, allocation, and thread
   scaling results without imposing correctness-test timing thresholds;
10. C API documentation and examples accurately describe ownership,
    provenance, field availability, and budget semantics.
