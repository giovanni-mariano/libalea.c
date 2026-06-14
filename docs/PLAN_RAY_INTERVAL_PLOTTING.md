# Ray-Interval 2D Plotting Plan

**Date:** 2026-06-14
**Status:** experimental implementation landed
**Scope:** faster 2D plot-time coverage/error detection using raycast-derived
row intervals.

## Where We Are

The current implemented experiment is a new opt-in coverage mode:

```c
ALEA_GRID_COVERAGE_RAY_INTERVALS
```

It is wired into `alea_find_cells_grid_coverage()` and into the Lua
`sys:find_cells_grid()` binding behind:

```text
ALEA_PLOT_RAY_INTERVALS=1
```

The implementation is intentionally not the older scanline fill path
`alea_find_cells_grid_raycast()`. That older path is fast, but it reports
`out_errors` as clean and is documented as fill-only. It also depends on the
hierarchical segment raycaster, which still has known lattice DDA boundary bugs.

The new experiment instead uses raycasting only to find row intervals:

1. Cast one hit-producing ray per slice row along the slice `u_axis`.
2. Split the row into intervals between sorted surface hits.
3. Sample exact coverage once at each interval midpoint.
4. Fill all pixel centers in that interval with:
   - `ALEA_COVERAGE_NONE` for undefined/gap,
   - `ALEA_COVERAGE_ONE` for valid single ownership,
   - `ALEA_COVERAGE_MULTI` for overlap.
5. Run the existing boundary-ambiguity filter to suppress shared-surface false
   positives.

This keeps the semantics of exact all-hit coverage while reducing repeated
point checks when many adjacent pixels share the same CSG ownership.

## Current Code Touch Points

- Public experimental flag:
  - `include/alea_slice.h`
  - `ALEA_GRID_COVERAGE_RAY_INTERVALS`
- C implementation:
  - `src/slice/slice_api.c`
  - `find_cells_grid_coverage_ray_intervals()`
  - interval filling helpers
  - fallback to existing exact coverage when ray intervals are unsupported or
    fail
- Lua opt-in:
  - `src/lua_bind/lua_slice.c`
  - `ALEA_PLOT_RAY_INTERVALS=1`
  - when enabled, returned grid tables include `coverage` and
    `secondary_cell_ids`
- Regression coverage:
  - `tests/unit/test_grid_query.c`

## What Works Now

Implemented and tested:

- partial overlap parity against `ALEA_GRID_COVERAGE_EXACT`
- nested/total overlap parity against exact coverage
- undefined gap parity against exact coverage
- lattice fixture parity through the global hit-producing raycast path
- coplanar split-surface suppression through the existing ambiguity filter
- Lua `find_cells_grid()` still works normally and with
  `ALEA_PLOT_RAY_INTERVALS=1`

Verification run:

```text
./bin/tests/unit/test_grid_query
./bin/alea tests/lua/test_slice.lua
ALEA_PLOT_RAY_INTERVALS=1 ./bin/alea tests/lua/test_slice.lua
make test-unit USE_OPENMP=1 RELEASE=1
```

All passed on 2026-06-14.

## Important Semantics

The ray interval path is a coverage/error detector, not only a color fill.

For each interval, exact coverage decides ownership. The ray itself does not
choose the cell. This matters because a ray segment has only one winning cell,
while plot diagnostics need count-up-to-two semantics:

```text
0 claimants -> undefined
1 claimant  -> valid
2+ claimants -> overlap
```

This is why the path can detect nested overlaps that have no visible cell-id
boundary in the primary grid.

## Fallback Policy

The path is conservative. It falls back to the existing exact coverage path for:

- unsupported `universe_depth` values; currently only `-1` is handled directly;
- invalid slice bounds or degenerate grid spacing;
- raycast setup or row raycast failures;
- exact interval midpoint coverage failures.

The fallback is intentionally exact, not fast, because the experimental path
must not hide geometry errors.

## Known Limitations

- No dedicated public stats yet for rows, intervals, hits, fallbacks, or
  verifier mismatches.
- No sampled verifier yet for production plots.
- No large-model benchmark has been committed for this path.
- The implementation uses the global hit-producing raycast path for correctness
  around lattices. It avoids the known hierarchical DDA bugs, but may be slower
  than the hierarchical segment-only path on some models.
- Interval midpoint classification assumes CSG ownership does not change inside
  a surface-hit interval. That should hold when all relevant surfaces are hit by
  the row ray, but degenerate/tangent/coincident cases need broader stress
  testing.
- The Lua switch is environment-gated rather than a formal options table.

## Next Plan

1. Add diagnostics.
   - rows traced
   - total hits
   - total intervals
   - exact interval coverage queries
   - row/ray fallbacks
   - exact fallback pixels
   - boundary-filter suppress/retain counts

2. Add a sampled verifier.
   - Proposed environment variable:
     `ALEA_RAY_INTERVAL_VERIFY=N`
   - Every Nth pixel, compare ray-interval classification against exact
     per-pixel coverage.
   - Report mismatches and optionally force exact fallback on mismatch.

3. Add a real benchmark executable.
   - Suggested name:
     `examples/c/slice_error_bench.c`
   - Compare:
     - point grid
     - old raycast fill
     - exact coverage
     - tile/path exact coverage
     - new ray-interval coverage

4. Run representative large-model tests.
   - Reuse the C-Model benchmark from `PLAN_FAST_PLOT_ERROR_DETECTION.md`.
   - Test multiple resolutions and slices.
   - Require zero verifier mismatches before promotion.

5. Improve degenerate handling.
   - Stress tangent rays.
   - Stress coincident/coplanar surfaces.
   - Stress repeated surfaces and macrobody-derived surfaces.
   - Consider local exact fallback around dense hit clusters.

6. Decide API shape.
   - Keep C flag experimental until validated.
   - Replace the Lua environment switch with an explicit options table if this
     becomes user-facing plotting behavior.

7. Promotion criteria.
   - Exact parity on unit/integration fixtures.
   - Zero sampled verifier mismatches on representative large models.
   - Faster than full exact coverage on large models.
   - No regression in false-positive suppression for shared boundaries.
   - Clear behavior for lattices before considering use of the faster
     hierarchical segment raycaster.

## Current Recommendation

Keep the feature experimental and opt-in:

```text
ALEA_PLOT_RAY_INTERVALS=1
```

Use it for benchmarking and validation, not as the default plot error path yet.
The implementation has the right structure and passes focused parity tests, but
it still needs diagnostics, verifier support, and large-model evidence.
