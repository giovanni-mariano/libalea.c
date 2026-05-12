# TLAS/BLAS Spatial Index Implementation Plan

**Date:** 2026-05-12
**Branch:** `feature/tlas-blas-spatial-index`
**Target model:** `/home/giovanni/projects/mcnp_files/E-lite_R250630.mcnp`
**Status:** implementation plan before replacing the flat spatial index

## Decision

Use the hierarchical two-level spatial index approach from
`PLAN_TWO_LEVEL_BVH_QUERY.md` and Phase 4 of
`PLAN_SPATIAL_INDEX_LARGE_MODEL.md`.

The current flat global spatial index is still the scaling limiter for
`E-lite_R250630.mcnp`. The model loads, but spatial build fails while
collecting expanded instances in `collect_instances_recursive()`, before BVH
construction. Further flat-instance compression is therefore unlikely to be
enough unless it changes the instance-count scaling law.

The recent per-universe point-BVH experiment is useful diagnostic work, but it
is not the final answer by itself. On E-lite-like queries it reduces exact cell
tests, but can still produce large candidate sets and extra bbox traversal
work. It should remain guarded by thresholds/environment switches while the
real spatial-index replacement is built.

## Goals

- Build a spatial acceleration structure for E-lite without materializing the
  full recursively expanded flat instance array.
- Preserve point, region, slice, and ray-style query correctness on normal
  models.
- Keep the existing flat spatial path available as a fallback until the
  hierarchical path has parity tests and E-lite validation.
- Keep lattice handling O(1) through existing lattice coordinate math instead
  of expanding lattice elements into flat instances.
- Make the memory/performance tradeoff measurable with dedicated diagnostics.

## Non-Goals

- Do not optimize the recursive all-cells path as the primary solution.
- Do not delete the existing TLS cache or flat spatial path until the
  hierarchical path is proven.
- Do not silently rebuild a flat instance array inside the hierarchical path for
  API compatibility.
- Do not expand every lattice element into a TLAS placement unless a later,
  measured case proves that this is necessary and bounded.

## Architecture

### BLAS: Per-Universe Local Geometry

Build one bottom-level acceleration structure per unique universe, in that
universe's local coordinates.

Each BLAS should contain:

- the universe id and local cell indices
- local cell bounding boxes
- a compact BVH node array
- an index/order array for candidate cells
- build/query statistics

Small universes should keep the current linear scan path. The threshold must be
configurable because the first experiment showed that an aggressive low
threshold can hurt performance on E-lite.

### TLAS: Placements Instead Of Expanded Instances

Represent recursive fill relationships as placements rather than terminal
global cell instances.

Each placement should contain:

- placement id
- parent placement id, or an invalid/root marker
- parent cell index
- universe id being placed
- fill/lattice metadata
- depth
- transform handle or composed transform reference
- conservative world bounding box when available
- flags for root, fill placement, lattice container, and terminal/container use

The TLAS may initially be an explicit placement graph plus traversal stack. A
separate BVH over placements can be added after the first correctness pass if
region/slice queries need it.

### Transform Storage

Transforms should be stored separately from hot traversal records.

Initial implementation:

- store composed local-to-world and world-to-local matrices per placement where
  needed for correctness
- use transform handles in placements and hits
- leave deduplication for a later optimization

The important first property is that transform storage grows with placements,
not with fully expanded terminal instances.

### Lattice Handling

Lattices remain a special fast path.

For point queries, traversal should:

1. transform the point into the lattice cell's local coordinates
2. resolve the lattice element with the existing rectangular/hex lookup code
3. descend into the resolved universe

For spatial build, the first hierarchical implementation should store a lattice
container placement, not one placement per lattice element. This is the key rule
that prevents accidentally recreating the E-lite flat-expansion failure.

## Implementation Phases

### Phase 0: Guardrails And Diagnostics

Keep the current flat path and experimental per-universe point BVH intact.

Add or extend diagnostics so every build mode can report:

- cell count and universe count
- BLAS count, node count, and memory estimate
- placement count and transform count
- lattice container count
- max placement depth
- peak RSS where the probe can report it
- query counts, candidate counts, exact containment tests, and fallback counts

Add a build mode switch:

- `ALEA_SPATIAL_INDEX_MODE=flat`
- `ALEA_SPATIAL_INDEX_MODE=hier`
- `ALEA_SPATIAL_INDEX_MODE=auto`

The default should stay conservative until parity and E-lite validation are
complete.

Exit criteria:

- existing unit and integration tests still pass
- `tools/large_model_probe` can exercise the new diagnostics without changing
  public behavior

### Phase 1: Hierarchical Data Structures

Add a new implementation unit, preferably `src/core/alea_spatial_hier.c` with a
matching internal header, unless the surrounding code strongly favors keeping
the first pass inside `alea_spatial.c`.

Define:

- `alea_universe_blas_t`
- `alea_spatial_placement_t`
- `alea_spatial_transform_ref_t` or equivalent transform handle
- `alea_hier_spatial_index_t`
- statistics structs for build and query diagnostics

Add allocation and cleanup hooks wherever the current spatial index is owned and
destroyed.

Exit criteria:

- the hierarchical index can be allocated and freed without leaks in small tests
- no public query path uses it yet

### Phase 2: Build BLAS Per Unique Universe

Build BLAS data from each universe's local cell list without recursive
expansion.

Implementation notes:

- reuse existing local bounding-box code where possible
- preserve universe cell order for deterministic query behavior
- keep linear scan for small universes
- store build failures as explicit status codes, not partial silent fallback

Tests:

- empty universe
- small universe below threshold
- larger universe above threshold
- overlapping local cell bounding boxes

E-lite gate:

- BLAS construction completes after model load
- memory increase is proportional to unique cells, not expanded placements

### Phase 3: Collect Placements Without Flat Expansion

Replace the flat recursive collection strategy with hierarchical placement
collection for the new mode.

Rules:

- create a root placement for the root universe
- create placement records for filled universes
- stop at fill/lattice boundaries instead of collecting all terminal global
  instances
- store lattice cells as lattice containers
- do not materialize per-lattice-element placements in the first pass

Exit criteria:

- E-lite placement collection completes without OOM
- placement count is explainable and much smaller than the failed flat expanded
  instance count
- diagnostics identify root/fill/lattice placement categories

### Phase 4: Hierarchical Point Query

Implement a new point-query entry point for the hierarchical index before
switching public callers.

Traversal:

1. start from the root placement
2. transform the query point into placement-local coordinates
3. query the universe BLAS or scan the small universe
4. run exact containment tests on candidates
5. record hits with cell id, universe id, depth, placement id, and transform
6. descend into fill cells
7. resolve lattice cells through existing O(1) lookup and descend

Correctness requirements:

- preserve overlap/all-cells semantics
- preserve expected deepest-cell behavior for single-cell callers
- preserve deterministic ordering by depth and universe cell order

Tests:

- simple single universe
- nested fill with transform
- multiple fills of the same universe with different transforms
- overlapping cells at the same depth
- rectangular lattice
- hex lattice, if existing fixtures support it
- comparison against `alea_find_all_cells_at_point_recursive()` on small models

Exit criteria:

- point-query parity tests pass
- query diagnostics show candidate and exact-test counts
- flat path remains selectable

### Phase 5: Region, Slice, And Ray-Style Queries

Port region/slice/raycast-style users only after point-query parity is stable.

Implementation notes:

- transform region boxes conservatively into local placement coordinates
- query local BLAS boxes for conservative candidates
- run existing exact geometry checks after candidate filtering
- make slice hit records carry enough placement/transform data to replace flat
  instance metadata

Tests:

- region query on root cells
- region query through nested fills
- slice query through transformed fill
- overlap grid behavior compared to the recursive reference path

Exit criteria:

- existing slice/render/raycast tests pass
- no hierarchical query path depends on a hidden flat global instance array

### Phase 6: Public Spatial API Migration

Switch public spatial-index construction to mode-aware behavior.

Proposed behavior:

- `flat`: existing implementation
- `hier`: hierarchical implementation, fail explicitly if an unsupported API is
  called
- `auto`: choose hierarchical when flat collection is estimated to be too large

Audit APIs that currently assume flat instances:

- `alea_spatial_index_instance_count()`
- `alea_estimate_instance_volumes()`
- any Lua or example code that expects stable flat instance ids

Do not fake flat instance counts by rematerializing the hierarchy. If a caller
needs per-expanded-instance volume behavior, keep it flat-only until a proper
hierarchical instance-id scheme exists.

Exit criteria:

- normal models can still opt into flat behavior
- E-lite can use hierarchical behavior
- unsupported flat-only APIs fail clearly in hierarchical mode

### Phase 7: E-lite Validation

Use `E-lite_R250630.mcnp` as the acceptance model.

Required checks:

- load completes
- universe/BLAS build completes
- placement collection completes
- hierarchical spatial build completes
- representative point queries complete
- at least one bounded region or slice query completes
- peak memory stays close enough to load RSS to avoid OOM

Run:

```sh
make test-unit
make test-integration
bin/large_model_probe /home/giovanni/projects/mcnp_files/E-lite_R250630.mcnp
```

Add more specific probe flags as they become available, for example:

```sh
bin/large_model_probe /home/giovanni/projects/mcnp_files/E-lite_R250630.mcnp --hier-build
bin/large_model_probe /home/giovanni/projects/mcnp_files/E-lite_R250630.mcnp --hier-query-samples 200
```

Exit criteria:

- no OOM during hierarchical spatial build
- no correctness regression in the automated tests
- no meaningful slowdown on small/normal models when flat mode is selected
- hierarchical mode has clear diagnostics for any remaining slow query cases

### Phase 8: Cleanup After Merge Confidence

Only after the hierarchical path passes the gates:

- make `auto` prefer hierarchical for large models
- retire or delete the old TLS cache if no remaining caller needs it
- delete the flat spatial path only if every public API either has a
  hierarchical equivalent or a deliberate flat-only replacement
- remove or demote the experimental per-universe point-BVH branch if it is fully
  superseded by the new BLAS

## Performance Risks And Mitigations

### Candidate Explosion

Risk: local BLAS queries can still return many candidates when cell bounding
boxes overlap heavily.

Mitigation:

- keep exact containment counters
- keep small-universe linear fallback
- tune thresholds with measured data
- consider better split heuristics only after the hierarchy is correct

### Transform Overhead

Risk: hierarchical traversal may spend too much time composing or applying
transforms.

Mitigation:

- precompute placement transforms where memory allows
- store hot transform handles separately from larger payloads
- avoid recomputing inverse matrices per query

### Lattice Expansion Regression

Risk: adding TLAS placements naively can recreate the same flat expansion that
breaks E-lite.

Mitigation:

- first implementation stores lattice containers
- diagnostics report lattice container count and placement count
- E-lite acceptance explicitly checks that placement count stays bounded

### API Compatibility

Risk: existing APIs and examples may assume flat instance ids/counts.

Mitigation:

- keep flat mode
- fail clearly for unsupported hierarchical calls
- add new hierarchical APIs instead of overloading flat semantics where needed

### Region/Slice Conservatism

Risk: conservative transformed bboxes may over-query and hurt slice performance.

Mitigation:

- implement point-query parity first
- add query counters for region/slice candidate rates
- optimize transformed-region pruning after correctness is locked

## Commit Plan

Use small commits with independent verification:

1. add hierarchical spatial structs and cleanup
2. build per-universe BLAS with diagnostics
3. collect hierarchical placements without flat terminal expansion
4. add hierarchical point query and parity tests
5. port region/slice query paths
6. add mode switch and API behavior
7. validate E-lite and update docs
8. clean up obsolete flat/TLS pieces after confidence

## Immediate Next Coding Step

Start with a scaffold that does not replace public behavior:

1. add `alea_spatial_hier` internal types
2. add build/free functions
3. build BLAS-per-universe diagnostics only
4. extend `large_model_probe` with a hierarchical build/probe mode

This gives an early E-lite memory signal before any risky public query migration.
