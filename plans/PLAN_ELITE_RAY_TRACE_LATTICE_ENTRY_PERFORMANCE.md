# Plan: Recover E-lite hierarchical ray and ray-slice performance

**Date:** 2026-08-09
**Status:** In progress — support-bounds pruning, lattice-entry exact-ancestor
event acceleration, canonical-owner gating, and coherent ordinary-fill ray
restarts implemented. The supplied non-lattice E-lite ray is now attributed to
the former per-restart canonical owner query. This recovered a substantial
part of the elapsed-time regression; the remaining baseline gap is tracked.
**Target branch:** `feature/unified-ray-query`
**Reference baseline:** `9313e4e` (2026-08-02, repeating `LAT=1` fill fix)
**Primary workload:** `/home/giovanni/projects/test_aleathor/mcnp_files/E-lite_R250630.mcnp`

**Comparability:** The recorded E-lite timings were obtained on a different
computer from the current checkout.  They are historical evidence and may be
used only for paired baseline/current runs under the same host and build
matrix; they must not be compared numerically with the local synthetic-fixture
measurements in the architecture plan.

## 1. Problem statement

On the E-lite model, both of these operations are approximately five times slower than with the 2026-08-02 libalea revision:

- `Model.trace()` for hierarchical rays;
- `plot_ray_slice()` for the plane `Z=60`, with `X=[-1700, 1700]` and `Y=[-1700, 1700]`.

The similar slowdown in a single trace is the important signal. It places the primary regression in the shared native hierarchical ray traversal, rather than in raster assembly, NumPy conversion, plotting, or OpenMP scheduling. The long single-core portion observed during a slice may still be a secondary load-balancing problem, but it should be investigated only after the per-ray cost is corrected.

The leading hypothesis is a bad interaction between the new lattice-entry search and simple repeating `LAT=1` placements. It is a hypothesis to prove with counters and focused tests before changing the algorithm.

## 2. Current evidence

### 2.0 Implementation progress (2026-08-10)

The first safe implementation slice is complete and covered by focused MCNP
integration tests:

- Repeating lattice occurrences now inherit trusted conservative bounds from
  enclosing fill cells when the hierarchical placement index is built. This
  removes inactive occurrences from TLAS candidates without treating a
  fallback/unbounded box as an exact clip.
- `alea_raycast_result_t` now records private per-ray lattice-entry calls,
  candidates, DDA steps, outcomes, and exact-ancestor surface work.
- The internal lattice-entry result is explicit (`NO_ENTRY`, `FUTURE_ENTRY`,
  `ALREADY_INSIDE`, `ERROR`); errors propagate to the caller instead of
  silently becoming `DBL_MAX`.
- While a candidate is outside its exact enclosing-fill CSG, the entry search
  jumps to the next enclosing-cell surface instead of walking repeated lattice
  pitches. Once the occurrence is active, normal lattice DDA still determines
  the selected element.
- Before publishing a synthetic lattice entry, every enclosing fill is now
  checked as the canonical deck-order owner. A shadowed later occurrence is
  rejected and left to the ordinary global surface walker until ownership can
  change.

The shell-support regression uses a placement bound spanning 200 pitches with
an exact support shell beginning at radius 99. Its hierarchical trace now
performs fewer than ten DDA steps across the entire ray; before this change the
first entry search would walk approximately one hundred inactive pitches.

This is deliberately a conservative partial implementation. Competing-owner
events are not yet folded into the event cursor: a shadowed occurrence is
correctly rejected, but the global walker—not a dedicated ownership cursor—
advances to the competitor's boundary. That is correct but may leave
performance on pathological overlapping input.

The large-model probe now reports all new counters. The supplied E-lite lattice
point `(-859.680, -566.265, 60.000)` was probed with +X traces from both
`x=-959.680` (200 units) and `x=-1200` (500 units). Both traces exercised the
lattice walker (`190` lattice surface tests) but recorded
`lattice_entry.calls=0`: they entered the lattice through an already-resolved
hierarchical path, not through the suspected unresolved-entry helper. Their
warmed scalar times were 1.232 s and 2.046 s respectively. They validate the
lattice location, but they are not performance evidence for this entry-search
fix. The plan still requires the actual slow trace origin/direction (or a
placement-chain inspection that reproduces it) before reporting a speedup.

### 2.5 Supplied non-lattice slow ray changes the primary diagnosis

The exact ray from `(-1800, 0, 60)` to `(1800, 0, 60)` was measured as
`origin=(-1800,0,60)`, `direction=(1,0,0)`, `t_max=3600`. It has no lattice
work: every lattice-entry counter and the active lattice-surface counter was
zero. Its warmed fast-segments trace was 2.749 s in the controlled run,
emitting 160 segments in 793 iterations.

The dominant cost is instead the general hierarchical ancestor path:

- 171,345 total primitive surface tests;
- 137,288 (80.1%) ancestor-cell surface tests;
- 34,057 terminal-cell surface tests.

A direct isolated build of the stated baseline (`9313e4e`) now exists for the
same compiler flags and probe command. It produces the **same** 160 segments,
793 iterations, and 171,345 surface tests as the current tree, but its warmed
mean is 0.746 s versus approximately 2.90 s current (3.9x slower). The
regression is therefore not extra lattice DDA, ancestor queries, or primitive
intersection count; it is per-iteration/per-intersection overhead introduced
by `4b1045c`'s rich-path/packet traversal rewrite. This comparison supersedes
the earlier assumption that fewer ancestor surface tests alone would recover
the reported timing.

`raycast_hier_path_ancestor_surfaces()` was introduced after the selected
baseline (`647edbb`). It recomputes the nearest enclosing-cell boundary during
each resolved child-cell step and is therefore the leading regression target
for this ray, independent of lattice entry.

An exact per-ray ancestor-event cache was implemented and benchmarked, then
removed. It reduced tests by 20.3% (171,345 to 136,484; 411 cache reuses), but
made the controlled warmed trace slower (2.962 s versus 2.749 s). Do not
reintroduce that cache without a design that has lower path-key overhead and a
repeated A/B wall-time win.

The follow-up attribution confirms that this ray has only one enclosing path
entry at a time (`max_depth=1`), yet it performs 791 ancestor boundary queries.
The parent changes frequently along the line: after recording the first 16
distinct ancestor cells, the remaining cells still accounted for 644 queries
and 93,395 surface tests. This rules out a single globally stable parent event
as the main solution.

A second, deliberately narrow cache for only an exactly identical depth-one
parent was also tested. It achieved the same 411 reuses and the same 136,484
surface tests, but lost a same-build A/B wall-time comparison (2.913 s cached
versus 2.897 s without it). A plane-only dispatch fast path likewise showed no
repeatable win. Neither experiment is retained in the implementation.

The next performance implementation must therefore reduce the cost of a
*missed* parent query, rather than only reusing an event. The appropriate
design target is a bounded, immutable per-cell surface accelerator for the
large ancestor CSG cells (with an explicit memory budget and a benchmarked
build threshold). It must preserve all candidate surfaces and exact closest-t
ordering; an AABB or owner approximation alone is not valid for this task.
The E-lite cell-surface index contains 9,823,752 references across 450,648
cells, but only 3,866 cells have more than 128 surfaces. That threshold is the
starting memory budget for the accelerator; building a structure for every
cell is not acceptable.

### 2.6 Fine-grained backup-branch bisect: canonical restart validation

The detailed history in `feature/unified-ray-query-backup-20260807` isolates
the actual regression boundary more precisely than the squashed main-branch
commit. With the same isolated build and E-lite fast-segments command:

| Revision | Warmed mean |
| --- | ---: |
| `3a5cfef` (logical lattice occurrences in paths) | 0.750 s |
| `7805626` (rebuild paths across lattice element transitions) | 2.917 s |
| current tree | about 2.90 s |

All three produce the same segments, iteration count, and surface-test count
for this ray. `7805626` introduced `hier_path_check_canonical_owners()` and
calls it from `alea_hier_spatial_find_path_from_parent()` on every restarted
prefix, including ordinary non-lattice fills. That helper invokes
`alea_hier_spatial_find_ordered_cell_in_universe()` for each ancestor, turning
the former cheap CSG containment check into a BLAS candidate query and
deck-order scan at many of the ray's 793 steps.

The implemented fix retains canonical ownership for root resolution and lattice
element transitions, but uses coherent ownership when ray traversal restarts
below an already-containing ordinary fill. That path validates the full
ancestor containment chain, then descends without re-running the costly
deck-order universe query. In illegal overlaps it intentionally preserves the
tracked owner; the existing canonical entry point remains available for strict
point-query or diagnostic callers. The non-overlap hierarchy regression and
the existing competing-owner coherence tests cover those two contracts.

On the supplied non-lattice E-lite ray, this retains the exact prior output
(160 segments, 793 steps, and 171,345 surface tests) and changes the warmed
mean from about 2.90 s to **1.409 s**: a 2.1x speedup. The `9313e4e` reference
is still about 0.746 s, so this commit isolates and removes the known ordered-
universe restart overhead but does not claim to have recovered every cost added
by the later rich-path traversal work.

### 2.1 Shared execution path

The relevant public operations converge on the cell-aware ray implementation:

```text
Model.trace()
  -> alea_raycast_hier_fast_segments()
  -> raycast_cell_aware_impl()
  -> raycast_cell_aware_resume()

compact/native ray-slice batch
  -> one raycast_cell_aware_resume() state per active ray
```

Consequences:

- `include_paths=False` does not avoid the regression; path result materialization is not the main suspect.
- A fivefold slowdown in one ray cannot be explained by the raster row scheduler.
- The scalar trace must be measured first. Slice measurements are useful afterward for throughput and load balance.

### 2.2 Regression window

Use `9313e4e` as the initial known-fast baseline. The first commits to bisect and inspect are the hierarchical/lattice traversal changes beginning with `4b1045c`, followed by `507e81e`, `4b17f31`, and `ab86055`. The later native raster commits (`62407d6` and `0e4ec28`) should remain in the comparison matrix, but the scalar-trace symptom makes them less likely to be the root cause.

This list is a starting point, not a substitute for a measured commit bisect.

### 2.3 E-lite geometry relevant to `Z=60`

The input is large (about 410 MB and 6.9 million lines), with six `LAT=1` cell cards and thousands of fill placements. In particular:

- cell `586239`: `U=20043 FILL=20046 LAT=1`, whose fundamental element is described by surface `589099`;
- surface `589099`: `RPP 844.29 853.441 -77.96 -68.806 -33.71 62.011`;
- corresponding cell `616259` and surface `620720` have the same bounds.

The simple `LAT=1 FILL=N` form is an infinite repeating lattice. The `RPP` is
the fundamental element used to infer the pitch and element-local coordinate
mapping; it is **not** the world-space outer bound of the lattice occurrence.
In the current index such a lattice placement is deliberately assigned a very
large bounding box. Its effective support is the intersection of the enclosing
fill-placement cells that instantiate universe `20043` or `20053`. Those parent
cells can be transformed and non-convex.

Consequently, `Z=60` and Y values near `[-77.96, -68.806]` are useful only
after a loaded-model inspection proves that a particular world-space occurrence
crosses them. Phase 0 must report the placement identifiers, transforms,
ancestor cells, and conservative world support bounds of the expensive
occurrences before freezing the focused E-lite rays.

### 2.4 Suspected failure mode

The current `lattice_first_valid_entry()` initializes `previous_valid` to false, starts at `t_enter`, and returns only for an invalid-to-valid transition strictly after `t_min + epsilon`.

For a repeating lattice, `lattice_raycast_interval()` begins at the caller's
clipped `t_min`, and every position maps to the single repeated fill. If the
first element is already valid:

1. the strict first-step return condition fails because `t_current == t_min`;
2. `previous_valid` becomes true;
3. every subsequent repeating element is also valid;
4. the DDA can walk the whole candidate interval before returning no future entry.

That search is reached through `system_first_lattice_entry()`, which visits
hierarchy placements and filters lattice placements while traversing. Because
repeating placements currently have effectively unbounded boxes, unrelated
occurrences can survive TLAS pruning. If invoked at multiple void or
unresolved-container intervals, it can repeatedly traverse placement
acceleration data and repeatedly perform an unnecessary lattice DDA for each
surviving occurrence.

The plan must verify all parts of that chain. A large DDA count alone is not enough to establish that every call is semantically unnecessary.

## 3. Goals

1. Identify the first regressing commit and quantify where the extra work occurs.
2. Restore correct first-entry semantics for repeating and finite lattices.
3. Avoid rediscovering the same lattice placement candidates at every ray interval.
4. Preserve exact segment ownership, path identity, boundary ordering, and forward/reverse behavior.
5. Recover E-lite scalar-trace performance close to the 2026-08-02 baseline.
6. Recover slice throughput and eliminate a pathological single-thread tail where practical.

## 4. Non-goals

- Do not redesign all hierarchical point location or all lattice traversal in the first fix.
- Do not optimize Python plotting before the native scalar regression is resolved.
- Do not weaken geometric validation merely to reduce the benchmark time.
- Do not add permanent wall-clock calls, logging, or environment-variable polling to the hot loop.
- Do not add the external E-lite input to the repository or to ordinary CI.
- Do not combine this work with unrelated ownership projection or public API changes.

## 5. Required semantic and geometry model

The first-entry helper previously overloaded a floating-point result and
`DBL_MAX` with too many meanings. Its internal result is now explicit:

```text
NO_ENTRY       candidate cannot contribute in the remaining ray interval
FUTURE_ENTRY   candidate becomes valid at a returned t > current_t
ALREADY_INSIDE the current ray point is already in a valid lattice element
ERROR          malformed state or a failed required geometric query
```

The result must also carry the candidate placement identifier, boundary `t`, a
strictly interior sample `t`, and, for `ALREADY_INSIDE`, enough resolved
hit/path state for the caller to adopt the canonical lattice path without a
second ambiguous point query. `ERROR` must propagate through
`system_first_lattice_entry()`, `raycast_cell_aware_resume()`, and the public
return code; it must never be converted to `NO_ENTRY` or `DBL_MAX`.

The entry query must also retain three different quantities:

- the current traversal position (`search_t_min`);
- the conservative occurrence-support interval entry (`support_t_enter`);
- the conservative occurrence-support interval exit (`support_t_exit`).

Clipping all three into one `t_enter` loses the information needed to distinguish “the ray will enter this valid repeating lattice later” from “the ray is already inside this valid repeating lattice.”

Use these geometry terms consistently:

- **fundamental element:** the lattice cell card's CSG and inferred pitch in
  element-local coordinates;
- **lattice occurrence:** one placement of that lattice through the hierarchy;
- **occurrence support:** the exact, possibly disconnected region where every
  enclosing placement cell is the canonical owner;
- **support bounds:** a conservative world-space AABB enclosing the occurrence
  support, used only for pruning;
- **active candidate:** a candidate for which canonical hierarchy resolution at
  the sampled point selects the same placement chain and a valid lattice
  location.

For a simple repeating lattice, the fundamental element is not an outer
container. A CSG check against it, if needed, must be performed after reducing
the point into the selected element frame. Future entry is controlled by the
ancestor occurrence support, not by the repeating lattice's own huge AABB.

## 6. Phase 0: Freeze the benchmark and add low-overhead observability

### 6.1 Reproducible build matrix

For every compared commit, record and hold constant:

- release/optimization mode;
- compiler and version;
- `PORTABLE` and architecture flags;
- OpenMP enabled/disabled state and runtime;
- thread count and affinity;
- Python environment and libalea consumer revision;
- cache warm-up policy.

First compare scalar native calls with one thread. Then compare native batches with 1, 2, 4, and the normal production thread count. Report minimum and median of repeated warmed runs; retain individual samples to expose outliers.

### 6.2 Counter block

Extend the existing internal `alea_raycast_result_t` work counters where that
keeps the data naturally per-ray; do not create a second competing
instrumentation mechanism. Add a compile-time or test-only counter block only
for hierarchy details that do not belong in the result. Update counters at
coarse branch points rather than timing inside the hot loop. At minimum collect:

- ray step/segment-loop iterations;
- existing-path coherence attempts and hits;
- upward-parent, root-universe, and full point-query fallbacks;
- containment checks, failures, and re-resolution retries;
- calls to `system_first_lattice_entry()`;
- placement TLAS nodes and leaves visited by those calls;
- lattice placement candidates examined;
- repeating candidates rejected by support bounds, exact ancestor validation,
  and canonical-owner validation;
- ancestor and competing-owner surface events examined while seeking an
  occurrence-support entry;
- lattice DDA steps;
- `NO_ENTRY`, `FUTURE_ENTRY`, and `ALREADY_INSIDE` results;
- terminal-cell, lattice-container, and ancestor surface tests;
- emitted segment and path counts;
- maximum counter values per ray in a batch, not only totals.

Implemented so far: per-ray entry calls, TLAS node/leaf visits, candidates,
DDA steps, outcome counts, ancestor surface/event counts, and canonical-owner
rejections, plus compact-batch maxima for each of those per-ray counters.

Expose new counters only through the already-private result or an internal
test/benchmark accessor. Do not change the stable public API for temporary
diagnostics. For the E-lite inspection, also print each expensive lattice
placement identifier, lattice cell, transform, ancestor cell chain, support
bounds, and per-placement work contribution.

### 6.3 Commit bisect

Run the same small ray set and counter report on the baseline and each suspect commit. The bisect result is the gate for Phase 2: if the first regressing commit does not correlate with lattice-entry calls/DDA work, revise the hypothesis before implementing a lattice-specific fix.

## 7. Phase 1: Create a small deterministic reproducer

Build a compact test model that captures the E-lite topology:

- a root void or container region;
- many placements, including non-lattice fills;
- a transformed universe containing a simple repeating `LAT=1` fill whose
  fundamental element is small but whose occurrence is limited by a larger,
  non-trivial ancestor fill cell;
- a ray that starts before the placement;
- a ray that starts inside it;
- a ray that misses it narrowly;
- a non-convex ancestor that produces two separated active intervals for one
  lattice occurrence;
- an overlapping earlier-deck-order fill that makes the lattice occurrence
  geometrically containing but non-canonical;
- both forward and reverse directions.

The fixture should demonstrate the current excessive DDA or repeated placement scan without depending on the 410 MB external file. Check both returned segments and diagnostic work counts.

Keep the E-lite file as an opt-in local benchmark/smoke test configured by a path argument. If absent, the benchmark must skip cleanly.

## 8. Phase 2: Correct first-valid-entry behavior

### 8.1 Build conservative occurrence-support bounds

Before relying on `placement_t_enter`, make lattice occurrence bounds useful and
correct. While collecting placements, carry a conservative world-space support
AABB formed by intersecting the transformed conservative bounds of all
enclosing fill cells. For a finite lattice, also intersect its finite lattice
extent. For a repeating lattice, do not intersect the support with the
fundamental-element `RPP`.

Only use ancestor bounds that are known conservative. The current artificial
fallback used for unbounded or unknown cell bounds must not become a correctness
clip; mark such bounds as untrusted and leave that axis or occurrence unbounded.
The resulting support AABB is a pruning envelope, not proof of CSG containment.
Add tests proving that rotated bounds remain conservative and that unknown or
unbounded ancestors cannot cause false rejection.

Use the support bounds in the existing placement TLAS and expose their ray
entry/exit to the lattice visitor. This removes the huge-box pruning failure
where possible and provides a meaningful candidate interval.

### 8.2 Add exact candidate activation and boundary cursors

Refactor `lattice_first_valid_entry()` and its caller context so the helper knows
the unclipped support interval, the current search position, and the placement
identifier. Candidate validity is the conjunction of:

1. the point lies within the conservative support interval;
2. all enclosing placement cells contain it;
3. those cells are the canonical deck-order owners in their universes;
4. lattice index lookup succeeds; and
5. the resolved path belongs to the same lattice occurrence.

Add an internal hierarchy helper that performs this canonical occurrence
resolution and can return a resolved hit/path. Do not use the current
containment-only ancestor check as proof of canonical ownership. If comparing
occurrences through paths is otherwise ambiguous, add an internal placement
identifier to `alea_hier_ray_path_entry_t`; do not expose it in the stable
public path ABI merely for this fix.

At the candidate start:

1. compute a scale-aware forward sample offset;
2. clamp the sample strictly inside the candidate interval;
3. resolve the candidate through the exact canonical occurrence helper;
4. return an explicit result state or advance its boundary cursor.

For a simple repeating lattice:

- if the canonical candidate is active at the current traversal sample, return
  `ALREADY_INSIDE` immediately with the resolved hit/path;
- if it is active just after a support or ownership boundary greater than
  `search_t_min`, return `FUTURE_ENTRY` at that exact boundary;
- if the first sample is inactive, seek the next possible ownership boundary;
  do not DDA through periodic element boundaries, because every element maps to
  the same repeated fill and those boundaries cannot activate the occurrence;
- return `NO_ENTRY` only after exhausting the support interval's ownership
  events.

Obtain the next possible activation event from an **ownership-event cursor**.
It must include:

- surfaces of the occurrence's ancestor cells in their proper local frames,
  because those change containment; and
- surfaces of higher-precedence competing cells in each ancestor universe,
  because leaving an overlapping earlier-deck-order cell can make the candidate
  canonical without crossing a surface of the candidate itself.

Enumerate competing cells conservatively with the existing per-universe ray
acceleration and support interval; exact canonical sampling remains the final
authority. It is safe to visit extra mathematical surface intersections: after
every event, take a strictly forward sample and accept it only through exact
canonical occurrence resolution. Continue after rejected events so a
non-convex ancestor can be entered later or re-entered, and so a non-canonical
overlap can become canonical after its competitor exits. Use a scale-aware
progress rule and a structural event/step bound; never advance by repeatedly
adding a fixed epsilon.

Do not equate a support-AABB entry with an exact CSG entry. Axis-aligned boxes
may use a proven exact slab entry, but transformed, boolean, or non-convex
ancestors still require the surface-event plus canonical-sample procedure.

### 8.3 Preserve finite-lattice behavior

Finite lattices can contain invalid indices followed by valid elements. They
still require transition scanning. Merge two ordered event sources for each
finite candidate: lattice DDA boundaries and exact ownership events (ancestor
and competing-owner surfaces).
Validate the open interval after the earlier event through canonical occurrence
resolution. Preserve continuation state for:

- invalid-to-valid transitions;
- finite fill extents and holes;
- rays that can leave and re-enter valid islands;
- transformed and nested lattices.

The helper must make forward progress by a scale-aware tolerance and must never
rediscover the same boundary indefinitely. A candidate cannot be discarded
merely because an ancestor is inactive at its first support-AABB interval
sample; a non-convex ancestor or a finite lattice can have later active islands.

### 8.4 Define how the ray walker consumes each result

`raycast_cell_aware_resume()` must consume the explicit result as follows:

- `NO_ENTRY`: keep the normal physical-surface result for the current interval;
- `FUTURE_ENTRY`: if earlier than the current `t_next` by the accepted boundary
  tolerance, emit one synthetic lattice entry at the returned boundary and use
  the returned strictly interior sample on the next iteration;
- `ALREADY_INSIDE`: emit no boundary and do not change `t_current`; adopt the
  returned canonical hit/path and continue the current iteration's terminal,
  lattice, and ancestor surface search from that state. Guard this adoption so
  the same candidate cannot be adopted twice at the same `t_current`;
- `ERROR`: terminate the ray and propagate failure.

This contract is required because the entry search is currently called only
after the ordinary point resolver returned void or an unresolved container.
`ALREADY_INSIDE` is therefore a controlled resolution repair, not a zero-length
future entry. Count such repairs, and assert in tests that they do not duplicate
segments or loop.

### 8.5 Correctness gate

Before doing broader performance work, the focused tests must establish:

- a future repeating occurrence is returned at its exact ancestor entry;
- a ray already inside does not receive a false future entry;
- reverse rays choose the corresponding reverse entry;
- a missed/non-convex ancestor is not accepted from its support AABB alone, and
  a later valid island of that ancestor remains discoverable;
- the fundamental-element `RPP` is never treated as the outer extent of a
  simple repeating occurrence;
- finite invalid-to-valid transitions remain discoverable;
- nested ancestor validation rejects geometrically impossible candidates;
- geometrically containing but non-canonical overlapping placements are
  rejected in favor of the same deck-order owner selected by point resolution;
- `ALREADY_INSIDE` adoption emits no zero-length segment and makes progress;
- internal hierarchy failure propagates as `ERROR`, never as no entry.

## 9. Phase 3: Discover lattice candidates once per ray

The semantic fix and support-bound correction may eliminate most of the gap.
Enter this phase only if post-Phase-2 counters show that repeated candidate
enumeration remains material. Do not make a second index or per-ray dynamic
storage mandatory merely to satisfy an architectural preference.

### 9.1 Add a lattice placement index

During hierarchical spatial-index construction, retain a compact list of lattice
placement indices, conservative support bounds, and its count. This enables:

- immediate return when a model has no lattice placements;
- lattice-only candidate enumeration without visiting unrelated fill leaves;
- stable placement identifiers for diagnostics and canonical-owner checks.

The existing hierarchy stats already contain `lattice_placement_count`; use it
for the no-lattice early return without waiting for a new index. Start with a
direct lattice-only list only when the occurrence count is small. Add a
dedicated lattice-only BVH only if measured candidate counts justify its
construction and memory cost. A lattice-only BVH is allowed only after Phase
8.1 supplies conservative occurrence-support bounds; indexing the current huge
repeating boxes would not improve pruning. Use an explicit threshold selected
from benchmark data rather than an unconditional second BVH.

Likely implementation area: `src/core/alea_spatial_hier.c` and its internal/public header, with lifecycle ownership beside the existing hierarchical index.

### 9.2 Add per-ray candidate state

Extend the internal cell-aware ray state with a candidate cache:

- candidate placement identifier;
- conservative support `t_enter` and `t_exit`;
- deterministic enumeration key `(support_t_enter, placement_id)`;
- current ownership-event cursor, including ancestor and competing-owner
  surfaces;
- cached next active-entry result;
- finite-lattice DDA continuation/next-entry state where needed;
- an exhausted flag that is set only after every possible support island has
  been considered.

Build or lazily initialize this cache once for each scalar ray or packet lane, then consume it as `current_t` advances. Do not rescan the general placement TLAS at each void/container interval.

Important constraints:

- no shared mutable cursor across rays or OpenMP threads;
- no allocation at every segment step;
- storage grows once per ray or comes from thread-local scratch;
- cleanup occurs on success, early return, and error paths;
- a candidate capable of multiple valid islands due to either a finite lattice
  or non-convex ancestors cannot be discarded after its first entry;
- ancestor support bounds, transforms, and exact canonical validation are
  retained in candidate state;
- equal-distance candidates are resolved by canonical deck-order ownership;
  placement ID is only a deterministic secondary key, never the geometry
  precedence rule.

Likely implementation area: `src/raycast/raycast.c` and `src/raycast/raycast.h` internal state declarations.

### 9.3 Candidate-cache gate

On the E-lite ray set, placement traversal for first-lattice discovery should
occur at most once per ray (or once per explicitly documented cache rebuild),
not once per emitted interval. The cache must return the same canonical owner as
a fresh hierarchy query at every accepted entry. Counter assertions may be used
in performance-focused tests, but ordinary correctness tests should avoid
fragile machine-time thresholds.

## 10. Phase 4: Reduce remaining hierarchical re-resolution work

Only enter this phase if counters still show a material gap after Phase 2 and,
when justified by counters, Phase 3.

Investigate in this order:

1. why coherence rejects container cells;
2. whether unchanged placement paths can safely carry a resolved child hint across a boundary step;
3. which interior verification failures trigger the two re-resolution retries;
4. whether terminal, lattice-container, and ancestor surfaces are duplicated within one interval;
5. whether the root/full point lookup is repeated at an unchanged sampled point.

Any cache added here must have a precise invalidation event: crossing a tested surface, changing placement/lattice index, exhausting the candidate interval, or failing containment. Do not suppress validation based only on an epsilon-distance heuristic.

## 11. Phase 5: Address batch load balance separately

After scalar work per ray is near baseline, measure row duration or work counters in the native ray-slice batch.

The current static row scheduling can leave a single busy thread at the end if a small number of Y rows intersect disproportionately expensive geometry. Compare:

- `schedule(static)`;
- static scheduling with smaller chunks;
- `schedule(dynamic, chunk)`;
- `schedule(guided)`.

Choose using full-slice wall time and repeatability, including cheap models where dynamic scheduling overhead can lose. Prefer a fixed, documented policy or a model/work-size threshold; do not expose a tuning option unless real workloads require it.

This phase can improve the visible single-core tail, but it must not be presented as the fix for a fivefold scalar `trace()` regression.

## 12. Benchmark matrix

### 12.1 Focused E-lite rays

Start with X-directed probes using origin `(-1700, y, 60)`, direction
`(+1, 0, 0)`, and `t_max=3400`. The historical probe values are:

```text
y = -100, -80, -77, -73, -69, -60, 0
```

These values are discovery probes, not a claim that `[-77.96, -68.806]` is a
world-space outer lattice band. Freeze the regression set only after the Phase
0 placement report identifies which lattice occurrence and ancestor support
each ray actually crosses. Include at least one ray for every observed expensive
class: active at the origin, future ancestor entry, AABB-only miss, non-convex
re-entry, and unrelated occurrence rejected by pruning. Add rays through
transformed copies located during model inspection and repeat each ray in the
reverse direction.

Vary ray length/bounds to expose bad interval scaling:

- `[-425, 425]`;
- `[-850, 850]`;
- `[-1700, 1700]`.

A linear or worse rise in repeating-lattice DDA work with empty travel distance
is the expected signature of the current bug. After Phase 2, repeating entry
search must perform zero periodic DDA steps; work may grow only with relevant
support candidates and ownership events.

### 12.2 Slice cases

At `Z=60` and the production bounds, measure multiple resolutions such as 128, 512, and the normal production resolution. Record:

- total native traversal time;
- result projection/assembly time separately;
- per-thread or per-row work distribution;
- maximum ray work versus median ray work;
- segment checksum and output dimensions.

Test the minimal field projection used by plotting and owner/path-enriched variants separately. Run `include_paths=False` and `True` where the API permits it.

### 12.3 Control workloads

Include:

- no lattices;
- one repeating lattice;
- one repeating lattice with a small fundamental `RPP` inside a much larger
  ancestor support;
- one finite lattice with holes;
- deeply nested transformed placements;
- a non-convex ancestor with two active ray intervals;
- overlapping ancestor fills with known deck-order precedence;
- many non-lattice placements and one lattice placement;
- rays entirely in void;
- rays starting inside a lattice.

These controls prevent an E-lite-specific shortcut from becoming a general regression.

## 13. Test plan

Add low-level cases primarily to `tests/unit/test_raycast.c`,
`tests/unit/test_raycast_perf.c`, and `tests/unit/test_spatial_hier.c`. Put MCNP
simple-repeat parsing and deck-order fixtures beside the existing lattice-entry
coverage in `tests/integration/test_mcnp.c`. Reuse existing fixtures where
possible.

Required coverage:

1. repeating lattice, future entry through an ancestor cell;
2. repeating lattice, origin already inside;
3. repeating lattice, reverse ray;
4. transformed repeating placement;
5. nested transformed ancestors;
6. finite lattice with invalid-to-valid and valid-to-invalid transitions;
7. multiple separated valid islands in one finite candidate;
8. non-rectangular ancestor whose support AABB is entered but CSG is missed;
9. non-convex ancestor with a missed first interval and a valid later entry;
10. simple repeat whose fundamental-element bounds are far smaller than its
    ancestor support;
11. trusted, rotated, untrusted, and unbounded ancestor support bounds;
12. overlapping candidates where canonical deck order differs from placement
    ID order;
13. `ALREADY_INSIDE` path adoption without a synthetic boundary;
14. injected internal candidate-resolution error propagation;
15. ray exactly on or within tolerance of a placement/lattice boundary;
16. `t_max` clipping before, at, and after entry;
17. void and container-cell transitions;
18. no-lattice model fast rejection;
19. scalar, generic batch, compact slice, and packet-lane parity;
20. identical results with one and multiple OpenMP threads;
21. cleanup under early stop, callback stop, and error return;
22. ASan/UBSan and leak-check runs for new candidate storage.

Compare complete segment sequences, including `t_enter`, `t_exit`, cell/material/owner identifiers, and path data where requested. Use the repository's accepted geometric tolerance; do not replace structural equality with only an image checksum.

## 14. Acceptance criteria

### Correctness

- All existing raycast, spatial hierarchy, lattice, MCNP integration, compact batch, and bidirectional tests pass.
- New focused tests pass in forward and reverse directions.
- No zero-progress loop or repeated emission at one boundary is possible.
- A simple repeating lattice uses its fundamental cell only for periodic
  element mapping and uses its canonical ancestor occurrence for world support.
- Every accepted candidate matches canonical deck-order point/path ownership;
  equal-distance placement ID ordering cannot override geometry precedence.
- Candidate and hierarchy failures propagate as errors.
- Scalar, batch, compact-slice, and packet paths return equivalent segment semantics.
- E-lite produces the same accepted material/owner segmentation as the corrected reference behavior.

### Performance

- Repeating first-entry search performs no periodic lattice DDA; it scales with
  pruned support candidates and relevant ownership events instead of ray
  length divided by lattice pitch.
- If Phase 3 is justified and implemented, first-lattice placement discovery
  does not repeatedly traverse the general placement TLAS for each interval of
  one ray.
- On the fixed build matrix, representative E-lite scalar traces are at most 1.2x the 2026-08-02 median, or every remaining difference is explained by measured required work and explicitly accepted.
- The production `Z=60` slice shows no approximately fivefold regression against the baseline.
- Multi-thread slice completion has no large pathological single-thread tail attributable to a few statically assigned rows.
- Non-lattice and small-model control benchmarks do not materially regress.

Keep hard wall-time assertions out of normal unit tests. Store timings and counter summaries in benchmark reports; use structural counter bounds only for deterministic synthetic fixtures.

## 15. Likely files affected during implementation

- `src/raycast/raycast.c`: entry result semantics, per-ray cache consumption, counters, and batch scheduling experiment;
- `src/raycast/raycast.h`: internal state/result declarations;
- `src/core/alea_spatial_hier.c`: conservative occurrence-support bounds,
  canonical occurrence validation, and, if justified, a lattice placement list
  and optional dedicated index;
- corresponding spatial hierarchy header: lifecycle/query declarations;
- `src/raycast/raycast_api.c` only if an internal benchmark hook cannot remain below the API layer;
- `tests/unit/test_raycast.c`: semantic and scalar/batch parity cases;
- `tests/unit/test_raycast_perf.c`: deterministic work-count and microbenchmark cases;
- `tests/unit/test_spatial_hier.c`: lattice-only candidate index tests;
- `tests/integration/test_mcnp.c`: simple-repeat, canonical precedence, and
  transformed/non-convex ancestor tests;
- an existing benchmark/example location for the optional external E-lite driver.

Exact headers and benchmark placement should follow the current ownership boundaries discovered during implementation. Avoid promoting diagnostic types into `include/alea_raycast.h` unless a durable public use case is established.

## 16. Risks and mitigations

- **Fundamental-element/support confusion:** the simple-repeat cell CSG is one
  periodic template, not the occurrence's world extent. Keep these types and
  coordinate frames explicit in names and tests.
- **AABB/CSG confusion:** an early repeating-lattice return could create false
  entries. Treat support bounds only as pruning and require exact canonical
  occurrence validation after surface events.
- **Unsafe support clipping:** artificial or non-conservative ancestor bounds
  could remove real geometry. Track bound trust and intersect only proven
  conservative bounds.
- **Boundary tolerance changes:** sampling changes can alter ownership at coincident surfaces. Centralize and test the forward-progress/sample tolerance.
- **Finite/non-convex island regression:** a one-shot candidate cache could miss
  later valid islands. Preserve resumable DDA and ownership-event cursor state.
- **Memory growth:** per-ray candidate arrays can multiply in large packets. Use compact entries, small inline capacity or reusable scratch, and measure peak memory.
- **Wrong deterministic owner:** stable placement ordering can still disagree
  with deck-order geometry precedence. Resolve ties canonically; use placement
  ID only as a secondary deterministic key.
- **Suppressed failures:** converting an internal query failure to `NO_ENTRY`
  can silently erase geometry. Propagate the explicit `ERROR` state end to end.
- **Overfitting to E-lite:** retain synthetic transformed, finite, non-lattice, and miss cases.
- **Misleading scheduler win:** dynamic scheduling can hide expensive rays without fixing them. Gate scheduling work on scalar recovery.

Each phase should be separately reviewable and revertible. In particular, land the semantic fix before the candidate-index optimization, and keep scheduling changes separate from traversal changes.

## 17. Implementation sequence and review gates

1. Record the build matrix and baseline timings/checksums.
2. Add temporary/internal counters and run the commit bisect.
3. Land the synthetic reproducer as a failing regression test.
4. Add explicit first-entry result/error semantics and the internal canonical
   occurrence resolver, with no fast-path behavior change yet.
5. Add trusted conservative occurrence-support bounds and verify that TLAS
   pruning cannot reject real geometry.
6. Implement repeating-candidate ownership-event cursors and the defined
   `ALREADY_INSIDE`/`FUTURE_ENTRY` consumer behavior; assert zero periodic DDA
   work in the synthetic repeat fixture.
7. Merge finite-lattice DDA and ownership events, then validate all
   finite, nested, transformed, non-convex, and canonical-overlap cases.
8. Rerun E-lite scalar rays and the scalar, batch, compact-slice, reverse-ray,
   and sanitizer suites.
9. Add the no-lattice early exit using existing hierarchy stats. Use counters
   to decide whether a lattice-only list/index and per-ray candidate cache are
   warranted; skip them if traversal is no longer material.
10. If warranted, add the lattice-only index and per-ray/lane cursor cache with
    complete cleanup, then repeat correctness, memory, and performance gates.
11. Use counters to decide whether Phase 4 resolution caching is warranted.
12. Measure row imbalance and select a scheduler only if the tail remains material.
13. Remove or compile out temporary diagnostics, retaining only useful benchmark instrumentation.
14. Publish a before/after report containing commit IDs, build flags, timings, counters, and output-parity evidence.

## 18. Questions to answer with measurements

- Which exact commit first increases scalar ray time and which counter grows with it?
- How many times per E-lite ray is `system_first_lattice_entry()` currently called?
- What fraction of its time is general TLAS traversal, transforms/ancestor validation, and lattice DDA?
- Which actual world-space occurrences and ancestor chains dominate at `Z=60`?
- How many repeating candidates survive the current huge bounds, and how many
  survive trusted ancestor-support clipping?
- Does exact canonical validation classify each expensive candidate as active
  now, future ancestor entry, non-canonical overlap, or no entry?
- How many ancestor and competing-owner surface events are examined per
  candidate after periodic DDA is removed?
- After the semantic/support fix, is any lattice-only list or BVH still
  justified, and are its occurrence bounds selective enough to help?
- If per-ray caching is added, are containment retries or duplicated surface
  tests the remaining gap?
- After scalar recovery, how skewed is per-row work and which OpenMP schedule wins across both E-lite and cheap controls?

The next implementation action should be Phase 0 instrumentation and the focused synthetic reproducer—not a scheduler change and not a broad hierarchy rewrite.
