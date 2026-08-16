# Ray-query architecture with first-class geometry diagnostics

**Status:** In progress — Phases 3–8 substantially implemented; Phase 9 measurement active

**Priority:** Correct architecture first; preserve and recover performance within
that architecture

**Threading:** OpenMP remains the current implementation. TinyPar is out of
scope.

## 1. Status, ordering, and relationship to other plans

This plan defines the target ray-query architecture and makes geometry-error
detection a first-class consumer. It should become the architectural authority
for ray traversal, ownership, provenance, and diagnostics.

It incorporates the strongest parts of:

- `plans/UNIFIED_RAY_QUERY_ARCHITECTURE.md`;
- `plans/raycast_consolidation_prerequisite.md`;
- the user proposal, *Lean unified ray traversal with optional surface
  provenance*;
- `plans/PLAN_ELITE_RAY_TRACE_LATTICE_ENTRY_PERFORMANCE.md`;
- `plans/PLAN_NATIVE_BATCH_LATTICE_DDA.md`;
- `plans/SURFACE_BOUNDARY_PROVENANCE_LABELS.md`;
- `plans/RAY_EVENT_SURFACE_ID_PLOTTING.md`;
- `plans/NATIVE_BIDIRECTIONAL_RAY_SLICE_VALIDATION.md`.

Those documents remain useful as implementation history, performance evidence,
or focused subplans. Where they conflict on architecture, this plan takes
precedence after review and acceptance.

The implementation order is:

1. Establish the semantic architecture in this plan.
2. Recover known scalar performance within it.
3. Complete renderer and batch migration onto it.
4. Begin `plans/luporini_implementation.md` for executor-owned scheduling,
   worker arenas, and compact publication.

The executor plan must not begin by inventing another traversal abstraction.

### Implementation checkpoint (2026-08-15)

The selected-owner walker and global coverage sweep are now separate concrete
engines.  Scalar first-cell, first-visible, any-hit, selected-segment, and
X-ray consumers use the selected walker; coverage classification and complete
boundary provenance use global breakpoints and occurrence-sensitive all-owner
resolution.

Geometry validation is the active migration area.  It consumes complete
coverage intervals for overlap and undefined-fill findings, groups coincident
mathematical hits into one crossing decision, preserves adaptive boundary
sampling for transition evidence, and refreshes concrete ownership across
synthetic lattice DDA transitions.  Tests now cross-check selected traces
against uniquely owned coverage for both flat and lattice fixtures (including
exterior gaps), assert visible truncation under the crossing budget, and verify
that requesting segment surface IDs does not change selected endpoints or
ownership.

The remaining Phase 5 work is to complete the remaining
production-versus-coverage matrix.  The validator now retains the complete
coverage trace long enough to use exact before/after ownership at grouped
crossings, while retaining adaptive samples for ambiguity evidence and for
legacy small-gap neighbor probing.  The diagnostic limit audit now propagates
the shared nested-lattice recursion budget and reports saturated coverage owner
sets as explicit truncated intervals rather than discarding the ray with an
overflow failure.  Global breakpoint
enumeration respects enclosing fill transforms and transformed nested lattice
occurrences; parity accepts the explicit undefined-fill classification and the
native DDA boundary tolerance, including reverse traversal through the
transformed nested-lattice fixture.  Phase 6 has begun by consolidating AUTO
and FAST_FORWARD first-cell, first-visible, and any-hit dispatch onto their
shared selected-walker implementation.  Scalar X-ray compositing now streams
selected intervals directly and can stop at saturation without publishing a
per-pixel segment vector.  Compact X-ray tiles stream into fixed accumulators
rather than a public batch-result CSR round trip.  X-ray rendering now owns one frame-level OpenMP
tile region with serial worker-local ray consumers, rather than opening a
parallel batch region for every tile.  Performance
recovery, and executor work remain later phases.

The X-ray performance fixture now measures the fixed-output batch visitor used
by compact tiles, rather than the retired CSR publication route.  Its timings
are informational until repeated under a controlled benchmark protocol.

### Local Phase 9 baseline (2026-08-16)

Three warmed, single-thread fixture runs on this checkout produced these median
costs: 20-shell full trace 146.2 us/ray versus early first-visible 12.3
us/ray; packet/scalar first-visible 6.0/6.8 us/ray; packet/scalar any-hit
5.6/6.7 us/ray; scalar/fixed-output X-ray visitor 73.5/72.4 us/ray; and the
96x96 fixed-tile X-ray frame 38.2 us/pixel.  These figures establish a local
comparison point only: compiler flags, affinity, and the optional E-lite input
are not yet controlled or available in this checkout.

The historical E-lite timings in the focused performance subplan were measured
on a different computer.  They remain useful for paired A/B comparisons made
there, but are not numerically comparable to this local fixture baseline and
must not be used as a local acceptance threshold.  Any E-lite recovery claim
requires repeated baseline/current runs on the same host, compiler/build
configuration, thread setting, affinity policy, and warmed-cache protocol.

### Same-host release comparison (2026-08-16)

Three warmed `RELEASE=1 PORTABLE=1` runs compared this branch with
`feature/unified-ray-query` (`3858845`) in isolated worktrees on this host.
The current branch retained identical 20-shell segment and step counts, but
its median full selected trace was 114.4 versus 86.0 us/ray (+33%).  Scalar
first-visible was 9.03 versus 8.14 us/ray (+11%); scalar/compact batch segment
traces were 14.88/16.44 versus 11.91/12.80 us/ray (+25%/+28%); and the shared
X-ray scalar/fixed-output fixture was 56.91/59.67 versus 43.37/45.19 us/pixel
(+31%/+32%).  This is a same-work regression signal, not an E-lite result;
the retained-output semantics must remain intact while Phase 9 attributes and
recovers the added overhead.

The regression was traced to whole-capacity hierarchy-path copies at every
one-interval walker yield.  `alea_ray_walk_t` now owns the authoritative path
in place; the stepper snapshots it only when midpoint verification actually
re-resolves a candidate.  Three portable-release candidate runs produced
medians of 60.8 us/ray for the full selected trace, 11.35/11.29 us/ray for
scalar/compact batch segments, and 34.97/34.02 us/pixel for scalar/fixed-output
X-ray.  Fresh reference runs in the same window measured 86.0, 12.24/14.47,
and 43.62/44.34 respectively.  Host variance remains visible, but the same
outputs and the recovery across both rich and streaming consumers confirm that
the avoidable resumable-state copying—not geometric work or scheduling—was the
dominant local regression.

The performance suite also now prints deterministic lattice-entry attribution
for the two critical local cases.  The transformed inactive occurrence makes
one entry call, tests one TLAS node, and is pruned before a leaf, candidate,
or DDA step; the exact ancestor-support case makes six entry calls, visits six
TLAS nodes and leaves, examines six candidates, takes five DDA steps, and
records 18 ancestor surface tests and eight ownership events.  This is a
reproducible local guard for the entry-search fixes while the optional E-lite
deck remains unavailable.  Its eight-ray compact-batch companion reports the
maximum per-ray values (6 calls, 6 TLAS nodes/leaves, 6 candidates, 5 DDA
steps, 18 ancestor tests, and 8 events), so aggregate throughput cannot hide
an expensive lattice-entry outlier.

A same-walker spot measurement separates selected-segment publication from
traversal: reusable rich segment-vector output took 145.5 us/ray for the
20-shell fixture, versus 137.2 us/ray for streamed fixed output, while
avoiding 3,584 retained result bytes.  This is not yet a controlled regression
threshold, but it confirms the renderer migration removes measurable
publication work without changing selected traversal semantics.

Fine-grained, counter-only attribution now covers each selected-owner
resolution tier (neighbor, coherent path, root universe, and full hierarchy),
selected boundary enrichment, rare retry-path snapshots including the live
entry count, selected intervals yielded, and result-vector capacity growths
including bytes.  Batch reports retain the maximum per ray so outliers are not
hidden by aggregate throughput.  On the 20-shell worst ray, neighbor reuse
resolved 19 of 39 attempts, coherent-path reuse resolved 0 of 20, root lookup
resolved 20 of 22, and the full fallback resolved 0 of 2.  The ray yielded and
enriched 41 selected intervals, took no retry snapshot, and a fresh compact
trace grew its result storage seven times by 3,584 bytes.  The warmed reusable
rich trace and fixed-output visitor both grew no result storage.

Three warmed portable-release runs with this attribution had medians of 53.44
us/ray for the full selected trace, 8.78/11.32 us/ray for scalar/compact batch
segments, 32.36/32.39 us/pixel for scalar/fixed-output X-ray, and 54.87/53.49
us/ray for rich/streamed selected publication.  Semantic query lowering alone
had a 12.95 ns median over one million calls.  Host variance prevents treating
the improvement over the earlier instrument-free candidate as a speedup, but
these samples expose no recovered-path regression from the coarse counters.
The release-only lattice-walk maybe-uninitialized warning was also removed by
initializing the prior lattice location before its guarded first use.

## 2. Decision summary

The architecture has two scalar traversal engines, not one monolithic walker:

1. **Hierarchical selected-owner walker** for production queries such as
   first-visible, any-hit, ordered material intervals, rendering, and transport.
2. **Global diagnostic coverage sweep** for completeness-oriented questions
   such as gaps, overlaps, undefined fills, all crossing participants, and
   canonical boundary evidence.

They share:

- Ray and numerical conventions.
- Primitive intersection kernels.
- Hierarchy occurrence identity.
- Lattice coordinate and transition rules.
- Crossing vocabulary.
- Ownership query services.
- Query validation, budgets, and error semantics.

They do not share a single hot loop. A selected-owner walker deliberately hides
co-claimants; a diagnostic sweep must enumerate them. Forcing both policies
through one loop would either make normal traversal expensive or make
diagnostics incomplete.

The central rule is:

> Share geometric truth and result semantics; specialize traversal strategy and
> materialization for the question being asked.

## 3. Why selected-owner traversal cannot diagnose all geometry errors

A production ray needs one owner at each point. In legal geometry that owner is
unique. In overlapping geometry, a selected-owner policy chooses one claimant
and can continue coherently through it.

That is insufficient for diagnostics:

- A cell entirely contained inside another can form an overlap without changing
  the selected owner.
- Two occurrences of the same cell definition under different transforms can
  overlap while sharing the same cell ID.
- A canonical deck-order owner still hides later co-claimants.
- Forward and reverse traces can agree while both choose the same incomplete or
  incorrect ownership interpretation.
- An ownership-transition stream cannot expose a surface that does not change
  selected ownership.

Reliable overlap and gap detection along a ray therefore requires:

1. A complete ordered breakpoint set for the ray.
2. A complete owner-set query in each open interval.
3. Occurrence-sensitive identity, not cell ID alone.

The existing `alea_ray_classify_intervals()` already demonstrates this model:
it uses global crossings as breakpoints and an uncached recursive all-owner
query in each elementary interval. This behavior should become an explicit
diagnostic engine rather than remain a special function at the end of
`raycast_api.c`.

## 4. Semantic vocabulary

Use these terms consistently.

### Mathematical hit

One primitive/ray intersection before ownership interpretation. It includes a
distance and canonical primitive identity. Multiple MCNP surface cards may map
to the same primitive.

### Crossing group

All mathematical hits coincident at one distance within the documented
tolerance, plus an optional synthetic lattice transition. A diagnostic crossing
group retains all distinct primitive participants subject to budgets.

### Selected owner

One concrete hierarchy occurrence chosen for production tracking. It is more
specific than a cell definition or MCNP cell ID.

### Coverage set

Every concrete owner occurrence containing an interval sample, organized by
hierarchy depth. This is the basis of overlap, gap, and undefined-fill
classification.

### Selected interval

A verified open range `[t_enter, t_exit]` associated with one selected owner,
or with selected void/undefined state.

### Coverage interval

An elementary or merged range associated with a complete coverage set and a
diagnostic classification.

### Ownership transition

A change in selected or complete ownership across a crossing group. A physical
surface hit is not automatically an ownership transition.

### Finding

A structured diagnostic conclusion supported by one or more intervals,
crossing groups, coverage samples, or directional comparisons.

## 5. Three ownership semantics

The architecture must represent three policies explicitly.

### `TRACK_COHERENT`

Purpose: normal rendering, transport, early-stop queries, and fast segment
tracing.

Behavior:

- Select a deterministic initial owner.
- Retain the current concrete occurrence while exact path containment remains
  true.
- Change ownership only at a resolved boundary.
- Rebuild only the invalid path suffix.
- Do not repeatedly run ordered-universe overlap resolution.

This policy may be history-dependent in invalid overlapping geometry. That is
acceptable only for APIs documented as production tracking.

### `SELECT_CANONICAL`

Purpose: deterministic point ownership and compatibility operations that
promise one history-independent owner.

Behavior:

- Apply the documented universe/deck ordering.
- Return one owner occurrence.
- Do not claim that the result proves unique ownership.

Canonical selection is not overlap detection: it still hides co-claimants.

### `COMPLETE_COVERAGE`

Purpose: overlap/gap/undefined-fill diagnostics and oracle comparisons.

Behavior:

- Return all owner occurrences at relevant depths.
- Retain occurrence identity for repeated/transformed/lattice placements.
- Preserve enough hierarchy information to distinguish a valid nested chain
  from multiple claimants at the same depth.
- Report truncation if a configured owner limit is exceeded.

Normal and canonical-selected results must agree wherever complete coverage
contains exactly one valid chain.

## 6. Architecture

```text
                         Semantic ray request
                                  │
                                  ▼
                      Validation and plan lowering
                product + requirements + ownership policy
                                  │
                 ┌────────────────┴────────────────┐
                 ▼                                 ▼
     Hierarchical selected-owner          Global diagnostic sweep
              walker                     breakpoint enumeration
                 │                                 │
       local next crossing                grouped global crossings
       coherent/canonical owner                    │
       verified selected interval          complete owner coverage
                 │                                 │
                 └────────────────┬────────────────┘
                                  ▼
                         Concrete consumers
       first-visible | any-hit | selected segments | boundary evidence
       coverage intervals | geometry findings | rendering | estimators
                                  │
                                  ▼
                 Public compatibility and batch adapters
```

### Shared geometry services

Both engines use the same services for:

- Primitive intersection.
- World/local ray transformation.
- Hierarchy occurrence and transform composition.
- Lattice location, entry, DDA transition, and suffix rebuilding.
- Exact containment.
- Surface-card to canonical-primitive mapping.
- Numerical progress and coincidence tolerances.

No consumer may maintain its own divergent lattice transition rules.

### Independent engine invariants

The hierarchical engine may use locality and selected ownership. It is not
required to enumerate surfaces that cannot affect its selected path.

The diagnostic engine must retain a completeness-first global breakpoint
enumerator independent of hierarchical selected-owner decisions. It must not be
implemented as “run the fast walker with stricter flags.”

This independence makes the diagnostic engine useful as an oracle for the
production engine.

## 7. Semantic request and internal plan

Keep existing public APIs as adapters. Do not introduce a new public ABI during
the initial migration.

Internally, validate a semantic request once and lower it into a small immutable
plan. The plan chooses the engine and geometric requirements; the engine does
not inspect public API field masks directly.

Illustrative request products:

```c
typedef enum {
    ALEA_RAY_PRODUCT_FIRST_CELL,
    ALEA_RAY_PRODUCT_FIRST_VISIBLE,
    ALEA_RAY_PRODUCT_ANY_HIT,
    ALEA_RAY_PRODUCT_SELECTED_INTERVALS,
    ALEA_RAY_PRODUCT_BOUNDARY_TRANSITIONS,
    ALEA_RAY_PRODUCT_CROSSING_GROUPS,
    ALEA_RAY_PRODUCT_COVERAGE_INTERVALS,
    ALEA_RAY_PRODUCT_DIRECTIONAL_COMPARISON
} alea_ray_product_t;
```

Illustrative traversal requirements:

```c
typedef struct {
    uint8_t need_selected_owner;
    uint8_t need_complete_coverage;
    uint8_t need_density;
    uint8_t need_surface_identity;
    uint8_t need_all_coincident_primitives;
    uint8_t need_normal;
    uint8_t need_occurrence_key;
    uint8_t need_projected_owner;
    uint8_t need_full_path;
} alea_ray_requirements_t;
```

The requirements are not the public field mask. They represent work the engine
must perform. Examples:

- A normal requires primitive identity and a valid occurrence transform.
- Neighbor lookup may require one winning surface token even when no surface is
  published.
- Projected owner requires an internal occurrence stack but not necessarily a
  published full path.
- Complete coverage requires all-owner queries but does not automatically
  require normals.
- Coincident labels require a crossing group, not a single segment surface ID.

The lowered plan rejects unsupported combinations before traversal begins and
records limits, tolerances, geometry generation, and required caches.

Diagnostic lowering also records the coverage domain:

- target universe occurrence or hierarchy depth;
- finite ray interval or enclosing validation bounds;
- whether exterior unowned space is allowed;
- whether material-zero cells are reported as owned void;
- how undefined or missing fills are classified.

A material-zero cell is still defined ownership and is not automatically a
gap. Likewise, unowned space outside the configured model domain may be allowed
exterior rather than an error.

## 8. Product-to-engine matrix

| Product | Default engine | Ownership | Materialization |
| --- | --- | --- | --- |
| First cell | Hierarchical walker | coherent selected | one fixed answer |
| First visible | Hierarchical walker | coherent selected | one fixed answer |
| Any hit/shadow | Hierarchical walker | coherent selected | one boolean |
| Material/path intervals | Hierarchical walker | coherent selected | stream or segments |
| Depth/cell/material rendering | Hierarchical walker | coherent selected | direct framebuffer |
| X-ray rendering | Hierarchical walker | coherent selected | streaming accumulator |
| Segment surface IDs | Hierarchical walker | coherent selected | one enter/exit ID |
| Complete boundary provenance | Global sweep | selected plus crossing groups | events |
| Coverage intervals | Global sweep | complete coverage | diagnostic intervals |
| Gap/overlap classification | Global sweep | complete coverage | findings |
| Point ownership compatibility | Point resolver | canonical selected | one owner |
| Forward/reverse comparison | Hierarchical walker twice | coherent selected | mismatch evidence |

An explicit diagnostic caller may compare hierarchical selected intervals with
coverage intervals, but one is not a substitute for the other.

## 9. Core internal data

The following are design sketches, not ABI commitments.

### Concrete owner reference

```c
typedef struct {
    uint32_t cell_index;
    int cell_id;
    int material_id;
    int universe_id;
    int depth;
    uint64_t occurrence_key;
    uint32_t path_slot;
} alea_ray_owner_ref_t;
```

Occurrence identity is mandatory internally. Cell ID alone is insufficient for
repeated and transformed placements.

### Minimal crossing token

```c
typedef enum {
    ALEA_RAY_CROSSING_NONE,
    ALEA_RAY_CROSSING_PHYSICAL,
    ALEA_RAY_CROSSING_SYNTHETIC_LATTICE,
    ALEA_RAY_CROSSING_UNRESOLVED
} alea_ray_crossing_kind_t;

typedef struct {
    double t;
    alea_ray_crossing_kind_t kind;
    int surface_id;
    uint32_t primitive_id;
    uint32_t frame_slot;
} alea_ray_crossing_t;
```

The token retains only continuation/provenance identity. Do not copy a complete
transform matrix into every interval when provenance is disabled. A normal or
public event is enriched from scratch-owned frame state before that state is
reused.

### Selected interval

```c
typedef struct {
    double t_enter;
    double t_exit;
    alea_ray_owner_ref_t owner;
    alea_ray_crossing_t enter;
    alea_ray_crossing_t exit;
    double density;
    uint8_t resolution_flags;
} alea_ray_selected_interval_t;
```

### Coverage interval

```c
typedef enum {
    ALEA_RAY_COVERAGE_UNIQUE,
    ALEA_RAY_COVERAGE_GAP,
    ALEA_RAY_COVERAGE_ALLOWED_EXTERIOR,
    ALEA_RAY_COVERAGE_OVERLAP,
    ALEA_RAY_COVERAGE_UNDEFINED_FILL,
    ALEA_RAY_COVERAGE_TRUNCATED
} alea_ray_coverage_kind_t;

typedef struct {
    double t_enter;
    double t_exit;
    alea_ray_coverage_kind_t kind;
    const alea_ray_owner_ref_t* owners;
    size_t owner_count;
    int diagnostic_depth;
    uint32_t flags;
} alea_ray_coverage_interval_t;
```

Scratch-backed spans are valid only until the next step. Consumers that retain
them must materialize them explicitly.

## 10. Hierarchical selected-owner walker

The production walker exposes `init` and `next selected interval` operations.

```c
int alea_ray_walk_init(...);

alea_ray_walk_status_t alea_ray_walk_next_selected(
    const alea_system_t* sys,
    const alea_ray_plan_t* plan,
    alea_ray_walk_t* walk,
    alea_ray_scratch_t* scratch,
    alea_ray_selected_interval_t* interval);
```

One step performs:

1. Resolve or retain the current concrete occurrence.
2. Find the nearest relevant terminal, ancestor, lattice, or void crossing.
3. Enforce numerical progress.
4. Verify selected ownership inside the interval.
5. Rebuild only invalid hierarchy/lattice suffix state.
6. Return one selected interval and minimal crossing tokens.

It does not:

- Decide whether a material is visible.
- Apply a consumer's stop policy.
- Append public vectors.
- Compute normals unless required.
- Enumerate complete owner sets.
- Perform canonical diagnostic overlap resolution in normal tracking.
- Create threads or prepare shared caches.

Concrete scalar consumers implement stopping and output. Avoid a generic
callback in the hot loop until profiling proves it harmless and useful.

## 11. Global diagnostic coverage sweep

The diagnostic engine is an ordered interval partitioner.

### Step A: enumerate breakpoints

Collect every relevant crossing along the requested finite ray interval:

- Root and terminal physical primitives.
- Filled-universe content in concrete transforms.
- Enclosing ancestor surfaces.
- Lattice wrapper surfaces.
- Synthetic lattice-element transitions where ownership/occurrence partitioning
  requires them.

Sort and group coincident crossings deterministically. Deduplicate duplicate
surface cards by canonical primitive for geometry classification while retaining
surface-card labels when requested for provenance.

Do not build selected ownership segments as an accidental side effect.

### Step B: classify open intervals

For every non-degenerate interval between breakpoints:

1. Choose a robust interior sample away from its boundaries.
2. Execute `COMPLETE_COVERAGE` resolution.
3. Organize claimants by hierarchy depth and concrete occurrence.
4. Classify:
   - no claimant inside the configured validation domain: gap;
   - no claimant outside an allow-exterior domain: allowed exterior;
   - one valid hierarchy chain: unique;
   - multiple claimants at one depth: overlap;
   - selected container with no resolved fill content: undefined fill;
   - owner limit exceeded: truncated.
5. Merge consecutive intervals only when their complete normalized coverage
   sets and classifications match.

### Step C: classify boundaries

When transition diagnostics are requested:

- Sample complete coverage on both sides with an adaptive offset ladder.
- Retain the crossing group and all distinct primitives involved.
- Compare expected adjacency with actual before/after coverage.
- Report missing neighbor, non-adjacent transition, ambiguous coincident
  boundary, or unresolved sampling.

This folds the sound parts of the current geometry validator into the same
coverage semantics without conflating findings with traversal output.

## 12. Geometry-diagnostic query products

Geometry diagnostics should consume ray infrastructure through explicit
products rather than inspecting renderer-oriented segments.

### Coverage findings

Required findings include:

- Gap interval.
- Allowed exterior interval when the caller requests explicit reporting.
- Overlap interval, including every retained claimant.
- Undefined fill interval.
- Truncated owner coverage.
- Unresolved numerical interval.

### Transition findings

Required findings include:

- Undefined after crossing.
- Overlap after crossing.
- Missing neighbor.
- Non-adjacent transition.
- Ambiguous/coincident boundary.
- Synthetic lattice transition with inconsistent occurrence change.

### Consistency findings

Required signals include:

- Forward/reverse selected-owner disagreement.
- Production selected owner absent from complete coverage.
- Hierarchical interval endpoint absent from the global breakpoint set.
- Global unique coverage disagreeing with hierarchical selected ownership.

Consistency findings are evidence of a traversal or ambiguity problem. A
forward/reverse mismatch is not by itself a geometry error, and agreement is
not proof that geometry is valid.

### Finding evidence

Every diagnostic finding should be able to retain, subject to field requests
and budgets:

- Source kind: random ray, explicit ray, slice row, surface sample, or curve.
- Ray origin/direction and `[t_enter, t_exit]` or crossing `t`.
- World crossing and interior sample points.
- Sample offsets and tolerances.
- Concrete owner occurrences and hierarchy depths.
- Surface IDs and canonical primitive IDs.
- Physical, synthetic, coincident, unresolved, and truncated flags.
- Geometry generation and query configuration.

Diagnostics must distinguish “no finding” from “coverage truncated” and
“query failed.”

## 13. Limits of ray-based validation

A coverage sweep can be complete along one ray if breakpoint enumeration and
owner resolution are complete. A finite collection of rays cannot prove an
arbitrary 3D CSG model valid.

The validation system therefore combines:

- Deterministic explicit regression rays.
- Random Cauchy-Crofton or bounded rays.
- Axis- and slice-aligned directional sweeps.
- Surface-normal samples on analytical boundaries.
- Slice-curve sampling for small or nested features.
- Forward/reverse consistency checks.

Reports and documentation must describe sampling coverage. “No errors found”
means no errors were found by the configured probes, not a mathematical proof
of global validity.

Surface-driven sampling remains necessary for hidden or small nested overlaps
that random rays can miss. Tangential and zero-measure contacts require
explicit boundary-oriented tests.

## 14. Surface provenance policy

Preserve these meanings:

- `-1`: no reportable surface.
- `0`: synthetic lattice/DDA boundary.
- `>0`: physical MCNP/user surface ID.

Keep three publication levels:

1. Selected intervals without surface publication.
2. Selected intervals with one reported enter/exit surface ID.
3. Crossing groups/boundary evidence with all requested participants,
   primitive identity, normals, and before/after ownership.

Turning provenance publication off must not:

- Change selected ownership.
- Change interval endpoints.
- Weaken exact boundary selection.
- Remove an internal winning surface token needed for neighbor traversal.

Turning provenance on must add only the requested enrichment and publication.
A single segment surface ID never claims to represent all coincident surfaces.

## 15. Scalar consumers

Build small concrete consumers over the appropriate engine:

```c
int ray_first_cell_serial(...);
int ray_first_visible_serial(...);
int ray_any_hit_serial(...);
int ray_selected_segments_serial(...);
int ray_xray_accumulate_serial(...);
int ray_boundary_transitions_serial(...);
int ray_coverage_intervals_serial(...);
int ray_geometry_findings_serial(...);
```

Each serial consumer:

- Never creates threads.
- Never prepares shared caches.
- Uses caller-owned reusable scratch.
- Does not allocate for fixed output.
- Allocates variable output only through an explicit result builder/arena.
- Preserves transactional failure semantics at public adapter boundaries.

Boundary and coverage consumers may use the global engine. Production consumers
use the hierarchical engine.

## 16. Packet and SIMD policy

Preserve the current measured lane-interleaved packet driver for first-visible
and any-hit batches. It is not SIMD, but existing measurements show substantial
throughput value. Treat it as a separate driver over the same exact selected-
interval transition semantics.

Do not make the scalar full-segment path carry packet lane machinery when a
lean direct scalar driver is faster.

Distinguish future SIMD work from the current packet driver:

- Current: several resumable scalar lanes interleaved and written directly to
  fixed SoA output.
- Future: vectorized primitive/AABB/surface kernels where profiling and lane
  coherence justify them.

Do not vectorize the complete hierarchical state machine before measuring where
coherence survives.

## 17. OpenMP and execution boundary

During this architecture plan:

- OpenMP remains the only parallel backend.
- Scalar engines and consumers contain no OpenMP calls.
- Public batch, renderer, or validator orchestration owns the parallel region.
- Shared caches are prepared before workers start.
- Worker calls use serial/no-cache consumers.
- Use local `num_threads(...)`; do not mutate process-global thread count.
- A build without OpenMP executes the same semantics serially.

Diagnostic batches partition rays across workers. Each worker owns scratch and
a local finding arena; results are merged deterministically by input ray index
and within-ray distance. Random-ray generation must be deterministic by seed
and ray index, independent of thread count.

The more general executor and per-worker CSR arena design remain later work in
`luporini_implementation.md`.

## 18. Performance principles

Architecture has priority, but the architecture must make fast behavior
possible by construction.

### Production queries do not pay diagnostic cost

Normal tracking must not perform:

- Complete owner-set resolution.
- Ordered-universe overlap scans after every coherent restart.
- All-coincident crossing enrichment.
- Public full-path materialization.
- Diagnostic sampling ladders.
- Finding construction.

### Optional work is controlled at its source

- Normals are computed only when requested.
- Full paths are published only when requested.
- An internal compact occurrence stack remains available when continuation or
  projected ownership needs it.
- Surface arrays and event groups are built only when requested.
- Instrumentation is optional or benchmark-only where hot-state cost matters.

### Copying and allocation are architectural concerns

- Copy only live path prefixes.
- Prefer reusable path buffers and swaps over full fixed-capacity assignments.
- Keep scratch separate from public result storage.
- Avoid per-ray allocation in fixed-output queries.
- Avoid rich per-ray temporaries when a consumer streams output.

### Benchmark interpretation

Track both work and wall time. Fewer primitive tests do not imply faster code,
and equal primitive-test counts do not imply equal execution overhead.

The known E-lite scalar regression remains a priority benchmark, but matching a
historical time is a performance target rather than permission to weaken
correctness or architecture.

## 19. Proposed source ownership

Move code only after interfaces stabilize.

```text
src/raycast/
    ray_types.h             shared ray/crossing/owner vocabulary
    ray_request.c           validation and internal plan lowering
    ray_walk.c              hierarchical selected-owner walker
    ray_global.c            complete breakpoint enumeration
    ray_coverage.c          all-owner interval coverage sweep
    ray_owner.c             coherent/canonical/complete ownership services
    ray_lattice.c           lattice occurrence, entry, and transition rules
    ray_query_scalar.c      concrete scalar consumers
    ray_output.c            rich, SoA, CSR, and finding builders
    ray_intersect.c         primitive kernels
    bvh.c                   acceleration structures
    raycast_api.c           public compatibility adapters

src/geo_validator/
    geo_validator.c         sampling orchestration and finding policy
    ray_findings.c          mapping coverage/transition evidence to errors
    ray_slice_validation.c  directional comparison product
```

The raycast module supplies evidence products. The geometry validator decides
which probes to run, which findings to report, and how to aggregate them. The
validator must not become a rendering API, and rendering must not infer errors
from its selected segments.

## 20. Migration phases

### Phase 0: freeze contracts and baselines

- Inventory every ray and validator entry point, direct caller, engine, ownership
  policy, output, allocation, cache, threading, and failure contract.
- Freeze deterministic fixtures for hierarchy, transforms, fills, rectangular
  and hexagonal lattices, coincidence, gaps, isolated overlaps, undefined
  fills, clipping, occurrence paths, and forward/reverse comparison.
- Record scalar, packet, compact-batch, renderer, and validator performance.
- Preserve the E-lite ray measurements and exact build commands.
- Record intermediate backup-branch revisions or archive reproducible patches;
  several cited hashes are not reachable from the current checkout.

**Gate:** maintained semantic matrix, reproducible tests, and benchmark record.

### Phase 1: establish vocabulary and policy boundaries

- Introduce internal names/types for crossing, selected owner, coverage set,
  selected interval, and coverage interval.
- Introduce the three ownership policies explicitly.
- Add a product-to-engine legality table in code/tests.
- Lower existing query descriptors and public options into internal requirements.
- Keep old functions as adapters; no algorithm change yet.
- Correct comments that call selected ownership “canonical” when it is coherent,
  or call directional agreement proof of valid geometry.

**Gate:** no output or performance change beyond measurement noise.

### Phase 2: separate state, scratch, and publication

- Split resumable traversal state from instrumentation and public rich results.
- Retain compact occurrence-sensitive hierarchy/lattice state.
- Stop copying full transforms or fixed-capacity paths when requirements do not
  need them.
- Introduce caller-owned reusable scalar scratch.
- Preserve compatibility adapters around the current implementation.

**Gate:** bit/tolerance-equivalent outputs and no material scalar regression.

### Phase 3: extract the hierarchical selected-interval walker

- Extract initialization and one verified selected interval at a time from
  `raycast_cell_aware_resume()`.
- Preserve numerical progress, retry, enclosing-fill validity, lattice suffix
  rebuild, and occurrence identity.
- Remove first-cell, first-visible, material-filter, output-vector, and public
  field-mask decisions from the geometric step.
- Implement coherent and canonical-selected owner services explicitly.
- Retain the old combined entry point as a temporary comparison adapter.

**Gate:** old/new hierarchical interval equivalence on uniquely owned geometry;
documented coherent/canonical difference on deliberate overlaps.

### Phase 4: promote the global coverage sweep

- Separate global breakpoint enumeration from selected segment construction.
- Group coincident crossings by primitive identity with optional surface-card
  labels.
- Implement complete occurrence-sensitive coverage resolution.
- Produce coverage intervals and merge identical coverage sets.
- Migrate `alea_ray_classify_intervals()` onto this engine without changing its
  public behavior.
- Add internal output capable of retaining more than two overlap claimants,
  subject to budgets.

**Gate:** exact extents for isolated overlap and gap fixtures; selected-owner
traversal is demonstrably unable to replace this oracle.

### Phase 5: migrate geometry validation

- Make the ray-driven validator consume grouped crossings and complete
  before/after coverage.
- Preserve adaptive sampling ladders and ambiguity reporting.
- Map coverage and transition evidence into existing structured errors.
- Retain fast adjacency checks only as an optimization that falls back to
  complete coverage when proof is unavailable.
- Preserve surface/slice-driven validation for hidden nested overlaps.
- Treat forward/reverse mismatch as a consistency diagnostic, not a geometry
  verdict.
- Add production-versus-coverage cross-checks for test and diagnostic modes.

**Gate:** existing validator tests pass, isolated/nested overlaps remain
detectable, and truncated diagnostics cannot report success silently.

### Phase 6: migrate concrete production consumers

In order:

1. First cell.
2. First visible.
3. Any hit/occlusion.
4. Complete selected intervals.
5. Projected owner and full-path adapters.
6. Streaming X-ray accumulation.

Verify each consumer against its old path before removing query-specific
branches from the compatibility adapter.

**Gate:** one scalar semantic implementation per production query.

### Phase 7: provenance and boundary products

- Make requested segment surface IDs an enrichment of selected intervals.
- Make complete coincident provenance consume grouped global crossings.
- Preserve physical/synthetic/unresolved distinctions.
- Compute normals only when requested.
- Keep dense slice ownership separate from sparse boundary labels/evidence.
- Verify surface-off/on ownership and endpoint parity.

**Gate:** optional provenance changes only requested output and cost.

### Phase 8: renderer, slice, and batch consolidation

- Make solid/depth/ID rendering consume fixed selected answers directly.
- Remove reconstruction of one-element rich results.
- Stream X-ray intervals directly into the pixel accumulator.
- Remove per-tile public compact-batch round trips from internal rendering.
- Preserve one frame-level OpenMP tile region with serial ray consumers.
- Make fixed-output batches share scalar semantics while preserving the measured
  lane-packet driver.
- Keep variable CSR publication compatible until executor arenas land.

**Gate:** identical requested image/slice products and simpler execution paths.

### Phase 9: performance recovery within the architecture

- Complete the E-lite attribution work.
- Measure query lowering, state copying, path copying, owner resolution,
  crossing enrichment, and publication separately.
- Apply coherent tracking to normal legal-geometry traversal.
- Reduce live-prefix and transform copying.
- Keep diagnostic work out of production queries.
- Preserve packet throughput.
- Add accelerators only when same-work overhead has been recovered and a
  measured memory/performance case remains.

This phase can begin measurement in Phase 0 and apply small proven fixes during
earlier phases, but it must not dictate incorrect interfaces.

**Gate:** no unexplained major regression; known remaining costs are attributed.

### Phase 10: handoff to executor architecture

Only after the prior gates:

- Implement worker scratch lifetime and capacity.
- Centralize OpenMP region ownership.
- Add fixed-output executor publication.
- Replace per-ray rich temporaries with per-worker variable-output arenas.
- Compact transactionally into public CSR order.
- Consider additional SIMD only from profiling evidence.

**Gate:** executor changes scheduling and publication, not semantic traversal.

## 21. Test matrix

| Area | Required cases |
| --- | --- |
| Primitive crossings | tangent, near-tangent, large `t`, coincident cards, distinct coincident primitives |
| Selected ownership | root, void, graveyard, ordinary fill, transformed fill, deep chain |
| Occurrence identity | repeated cell at different transforms, re-entrant fill, same cell in multiple lattice elements |
| Lattices | rect/hex, finite/repeating, inactive support, neighbor changes, forward/reverse |
| Coverage | unique chain, exterior gap, interior gap, isolated total overlap, partial overlap, more than two claimants |
| Undefined state | missing fill universe, unresolved container, truncated owner set |
| Transitions | expected neighbor, implicit neighbor, missing neighbor, non-adjacent transition |
| Provenance | no surfaces, one surface ID, all coincident participants, primitive, normal, synthetic DDA |
| Query products | first cell, first visible, any hit, selected intervals, coverage intervals, events |
| Directional checks | agreement, mismatch, same wrong selected owner, occurrence-key mismatch |
| Diagnostics | explicit ray, random rays, slice rows, surface samples, curve samples |
| Limits | events, segments, owners, paths, bytes, interruption, allocation failure |
| Threading | serial/OpenMP parity, deterministic ordering, deterministic random seeds |

For every uniquely owned fixture, assert parity among coherent selected,
canonical selected, and the unique chain from complete coverage.

For deliberate overlap fixtures, assert the intended difference rather than
forcing selected-owner outputs to pretend they represent complete coverage.

## 22. Performance acceptance

Correctness and architectural separation are mandatory. Performance is a
release-quality requirement but not a reason to merge diagnostic and production
semantics.

Track:

- Wall time and rays/pixels per second.
- Traversal iterations.
- Terminal, ancestor, lattice, and global primitive tests.
- Coherent path reuse and fallback counts.
- Complete owner queries.
- Live path bytes copied.
- Scratch and published bytes.
- Allocations per ray/tile/frame.
- OpenMP regions per operation.

Targets:

1. Production queries perform zero complete-coverage queries unless an explicit
   diagnostic comparison is requested.
2. Surface-off and path-off queries perform no corresponding publication or
   enrichment work.
3. The E-lite normal segment trace should approach the `9313e4e` same-build
   baseline; within 15% is a target, not permission to weaken semantics.
4. Any remaining material baseline gap must be attributed.
5. Existing packet first-visible/any-hit throughput must not regress by more
   than normal measurement noise without a documented architectural reason.
6. Architecture migrations should normally stay within a 10% scalar regression
   envelope per phase; exceptions require an explicit temporary record and a
   recovery phase before rollout.
7. Diagnostic cost is measured separately and may be higher because it answers
   a strictly richer question.

## 23. Non-goals

- Do not collapse production tracking and complete diagnostics into one hot
  traversal loop.
- Do not infer geometry validity from selected segments alone.
- Do not infer validity from forward/reverse agreement alone.
- Do not remove the independent global breakpoint oracle.
- Do not weaken concrete lattice occurrence identity or enclosing-fill checks.
- Do not remove surface provenance or make it mandatory for every query.
- Do not expose new public query structs until internal semantics are stable.
- Do not replace OpenMP in this plan.
- Do not introduce a general executor before scalar semantics stabilize.
- Do not add a large accelerator solely because a counter is large; require a
  wall-time and memory case.
- Do not claim that finite ray sampling proves a 3D model globally valid.

## 24. Completion criteria

This architecture is ready for executor work only when:

- [ ] The three ownership policies are explicit and tested.
- [ ] The hierarchical selected-owner walker yields verified intervals without
      consumer-specific output branches.
- [ ] The global diagnostic sweep enumerates breakpoints independently and
      yields complete coverage intervals.
- [ ] Isolated overlaps and gaps are detected with correct extents.
- [ ] Occurrence-sensitive overlaps cannot collapse by cell ID.
- [ ] Geometry validators consume explicit coverage/transition evidence.
- [ ] Directional mismatches remain diagnostic signals, not validity verdicts.
- [ ] Production queries perform no complete diagnostic ownership work.
- [ ] First-visible, any-hit, segments, and X-ray use concrete scalar consumers.
- [ ] Surface provenance is optional and ownership-neutral.
- [ ] The measured packet driver is preserved over shared exact semantics.
- [ ] Renderer and fixed-output batch consumers use the same scalar semantics.
- [ ] Serial/OpenMP, lattice, provenance, diagnostic, and failure test matrices
      pass.
- [ ] Known performance regressions are recovered or explicitly attributed.
- [ ] Public APIs remain compatible through adapters.

At that point, the executor plan can optimize how independent rays are
scheduled and how their outputs are published without reopening the meaning of
a ray query or a geometry finding.
