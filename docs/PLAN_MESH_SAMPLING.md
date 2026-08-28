<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO

SPDX-License-Identifier: MPL-2.0
-->

# Structured Grid Sampling Plan

## Status and scope

Implementation status (August 2026): phases 0-6 are implemented. The adaptive
grid is intentionally nonconforming and unbalanced; balancing remains a
consumer-driven option rather than a current requirement. OpenMP sampling is
available for lean fixed-grid results, while ordered callbacks, packed sparse
fractions, and adaptive global budgets retain serial execution. Geometry
ownership ambiguity classification remains deferred as agreed.

`alea_mesh` currently samples CSG geometry on a rectilinear hexahedral grid and
exports the result to legacy Gmsh or VTK. It does not generate a
boundary-conforming mesh. The implementation is useful as an experimental grid
sampler, but its data semantics, accuracy reporting, scaling, and exporters need
to be hardened before it is used as analysis infrastructure.

This plan covers:

- truthful sampled-composition semantics
- coherent dominant cell/material results and deterministic ties
- complete material discovery
- scalable time and memory behavior
- reliable and useful exports
- robust automatic bounds
- adaptive composition estimation and, later, adaptive grids
- C and Lua API documentation and tests

Classification of gaps, overlaps, undefined fills, and query failures is
explicitly deferred. Until that work is designed, this plan preserves the
current point-query behavior and does not claim that grid sampling validates
geometry ownership.

Boundary-conforming surface and volume meshing is also out of scope. It should
be a separate module or an integration with a dedicated meshing library.

## Terminology and invariants

Use the following terms consistently in code and documentation:

- **grid**: the rectilinear collection of hexahedral voxels
- **voxel**: one grid element
- **sampled fraction**: sample count for a material divided by total samples;
  it is an estimator, not an exact volume fraction
- **dominant material**: material with the largest sampled fraction
- **dominant cell**: the most frequently sampled cell among samples belonging
  to the dominant material
- **mixed voxel**: non-dominant sampled fraction exceeds the configured
  threshold
- **tie**: two or more materials have the same largest sample count

The following invariants should hold for every successful result:

1. Every voxel has at least one sample and one fraction entry.
2. Sampled fractions for a voxel are finite, non-negative, and sum to one
   within a documented floating-point tolerance.
3. `material_ids[v]` and `cell_ids[v]` describe a material/cell combination
   observed together. The cell is selected only from samples of the dominant
   material.
4. A tie is recorded explicitly and resolved using a stable policy independent
   of sample traversal order. The initial policy is the lowest material ID,
   followed by the lowest cell ID.
5. `unique_materials` contains every material present in any fraction entry,
   including materials that are never dominant.
6. Grid indexing is public and fixed: X varies fastest, followed by Y and Z.
7. Export never reports success after a failed or incomplete write.

## API direction

Keep the existing `alea_mesh_*` names for compatibility in the first pass, but
describe them as structured-grid APIs. Consider `alea_grid_*` names only at a
planned API-version boundary.

### Result metadata

Add enough metadata to make the estimator interpretable:

```c
typedef enum {
    ALEA_MESH_TIE_NONE = 0,
    ALEA_MESH_TIE_MATERIAL = 1u << 0,
    ALEA_MESH_TIE_CELL = 1u << 1
} alea_mesh_tie_flag_t;

/* One value per voxel when requested. */
uint32_t *sample_counts;
uint8_t  *tie_flags;
double   *estimated_errors;     /* added with adaptive estimation */
```

The existing fraction structs may remain source-compatible, but their
documentation and Lua fields must say `sampled_fraction`. If a breaking API
revision is acceptable, rename the C types and members in one coordinated
change rather than maintaining two permanent names.

### Optional result fields

Large grids should not pay for data they do not request. Add a field mask:

```c
ALEA_MESH_FIELD_MATERIAL_ID
ALEA_MESH_FIELD_CELL_ID
ALEA_MESH_FIELD_MIXED_FLAG
ALEA_MESH_FIELD_DOMINANT_FRACTION
ALEA_MESH_FIELD_SAMPLED_FRACTIONS
ALEA_MESH_FIELD_SAMPLE_COUNT
ALEA_MESH_FIELD_TIE_FLAG
ALEA_MESH_FIELD_ESTIMATED_ERROR
```

Preserve the current complete result as the compatibility default initially.
After downstream callers are migrated, reconsider a smaller default.

Sampling configuration and export configuration should be separated. In
particular, `alea_mesh_sample()` should not reject a sampling request because an
unused output-format field is invalid.

## Implementation phases

### Phase 0: establish the contract

- Update `alea_mesh.h`, the API reference, README, C example, and Lua example to
  call the operation structured-grid sampling.
- Document sampled fractions as estimates based on a stated sample pattern.
- Document X-fastest indexing and the exact mixed-threshold rule.
- Document that dominant-only Gmsh/VTK output loses composition unless extended
  data export is enabled.
- Mark geometry-ownership classification as deferred, not solved.
- Add a short compatibility note for every public field that will be renamed or
  reinterpreted.

Exit criterion: a caller can determine exactly what every returned value means
without reading the implementation.

### Phase 1: correct result semantics

- Sample `(cell_id, material_id)` together instead of tallying material alone.
- Select the dominant material first, then select its dominant cell from the
  same observations.
- Detect material and cell ties and store tie flags.
- Resolve ties using the documented stable lowest-ID rule.
- Accumulate `unique_materials` from all fraction entries, not just dominant
  material IDs.
- Define `mixed_threshold` precisely as the maximum tolerated total
  non-dominant sampled fraction:

  ```text
  mixed = (1 - dominant_sampled_fraction) > mixed_threshold
  ```

- Validate fraction normalization and result invariants before returning from
  debug/test builds.
- Avoid querying the voxel center twice in center-sampling mode.

Exit criterion: material/cell pairs are coherent, ties are visible and stable,
and minority-only materials appear in the global material table.

### Phase 2: harden export

- Add an internal result validator used by both exporters:
  - positive dimensions
  - increasing finite nodes
  - required arrays present
  - checked node/cell counts and index arithmetic
  - valid sparse spans and finite normalized fractions
- Check every write path using a small checked-output helper.
- Check `ferror()` and flushing for stream export.
- Check both write and close errors for filename export.
- Write filename exports to a sibling temporary file and rename only after a
  complete successful close, so an existing output is not replaced by a
  truncated file.
- Set detailed Alea errors for null arguments, invalid formats, malformed
  results, and I/O failures.

Add opt-in composition export:

- VTK: emit one cell-data fraction array per sampled material when the material
  count is below a configurable limit.
- Gmsh 2.2: emit `$ElementData` for mixed flag, dominant sampled fraction, tie
  flag, and opt-in per-material sampled fractions.
- Preserve physical groups as dominant-material groups and state that they do
  not encode mixtures.
- Consider VTK XML/XDMF as a later scalable format instead of extending the
  legacy writers indefinitely.

Validate generated files with an independent reader in integration tests, not
only substring checks.

Exit criterion: injected write failures return failure, successful files parse
externally, and composition can survive an opt-in export/import round trip.

### Phase 3: remove scaling bottlenecks

- Replace the insertion sort in `collect_unique()` with incremental set
  accumulation or an `O(N log N)` sort. Prefer collecting materials while
  voxel results are produced.
- Add checked byte-size helpers to every allocation, including fraction-vector
  growth.
- Allocate optional arrays only when requested by the field mask.
- Replace fixed 512-entry stack arrays with scratch sized from the actual sample
  count or a bounded per-worker scratch object.
- Build all query caches before entering the sampling loop, then keep the hot
  loop read-only.
- Process voxels in slabs/tiles to improve locality and make bounded-memory
  operation possible.
- Add parallel sampling after confirming that prepared point queries are safe
  for concurrent readers. Use per-worker fraction/material buffers and a
  deterministic merge.
- Add cancellation checks at least once per tile and an optional progress
  callback outside parallel critical paths.
- Provide a streaming/callback sampling API for callers that do not need the
  complete result resident in memory.

Add benchmarks for:

- uniform single-material grids
- two-material interfaces
- many-material grids
- nested fills and lattices
- center versus `n x n x n` sampling
- memory per voxel for each field-mask combination

Exit criterion: no quadratic stage remains; time scales approximately with
`voxel_count * samples_per_voxel`, and optional fields have measured memory
costs.

### Phase 4: robust domain bounds

- Prefer an explicit world-space AABB supplied by the caller for production
  analysis.
- Replace the hard-coded `1.0` auto-bound tolerance with scale-aware tolerances
  derived from geometry configuration or explicit mesh options.
- Compute an AABB of placed root-universe occurrences rather than a sphere over
  every stored cell, including unplaced universe definitions.
- Do not convert the AABB to a cube; retain independent bounds per axis.
- Make auto selection explicit with a bounds mode instead of overloading six
  zero values indefinitely.
- Define behavior for mixed custom/automatic axes.
- Record whether bounds were explicit or inferred and which tolerance/padding
  was used in the result metadata.

Exit criterion: translated, elongated, small-scale, large-scale, filled, and
lattice models receive reproducible world-space bounds without magic absolute
lengths.

### Phase 5: improve sampled-composition accuracy

Keep center and regular subcell sampling as predictable baseline modes, but do
not present them as exact volume integration.

- Add deterministic stratified or low-discrepancy sampling to reduce alignment
  aliasing. Store the pattern and seed in result metadata.
- Replace the fixed `subsamples_per_axis <= 8` policy with explicit work limits:
  maximum samples per voxel, total query budget, and optional time/cancellation
  limits.
- Add adaptive composition estimation within a fixed voxel:
  1. sample a coarse level
  2. refine voxels that appear mixed or fail stability checks
  3. stop on requested tolerance, maximum depth, or work budget
  4. report sample count, estimated error, and limit-reached status
- Distinguish empirical convergence estimates from rigorous geometric error
  bounds.
- Investigate CSG box classification or boundary-assisted subdivision for
  rigorous lower/upper material-volume bounds. This is a separate milestone and
  must work through universe placements, not only flat cell trees.

Exit criterion: callers can request an accuracy/work target and can tell when a
voxel failed to meet it.

### Phase 6: adaptive grid refinement

Only after fixed-grid composition estimation is trustworthy:

- Introduce a separate adaptive-grid result representation.
- Refine mixed or high-error voxels into octree children.
- Balance neighboring refinement levels if required by consumers.
- Define stable parent/child IDs and material/cell metadata.
- Export adaptive cells as unstructured hexahedra.
- Keep the existing rectilinear result unchanged for non-adaptive sampling.

This phase improves spatial resolution but still does not create a
boundary-conforming mesh.

## Test plan

### Semantics

- exact 50/50 plane split records a tie and follows the stable tie policy
- reversing sample traversal does not change dominant IDs
- dominant cell belongs to the dominant material
- two cells with the same material select a coherent dominant cell
- minority-only materials appear in `unique_materials`
- sampled fractions are normalized for every sampling mode
- thresholds immediately below, at, and above a sampled minority fraction obey
  the documented strict comparison
- center, near-corner, and subcell modes report their actual sample counts

### Sampling quality

Use analytic or high-resolution references for:

- axis-aligned two-material planes
- offset thin slabs
- small spheres away from regular sample locations
- spherical shells
- repeated periodic structures that expose sampling aliasing
- non-uniform rectilinear cells with very different aspect ratios

Assert convergence trends rather than exact equality for sampled curved
volumes. Include cases where a coarse pattern misses an inclusion and confirm
that adaptive/stratified modes either find it or report insufficient accuracy.

### Hierarchy and bounds

- translated fills
- nested fills
- rectangular and hexagonal lattices
- repeated universe placements
- unplaced universe definitions
- elongated geometry
- geometry at very small and very large coordinate scales
- mixed explicit/custom/automatic axes

Ownership ambiguity behavior remains outside this plan, but ordinary unique
hierarchical cases must sample and bound correctly.

### Export and failure handling

- uniform VTK parses as a structured dataset
- non-uniform VTK parses as a rectilinear dataset
- Gmsh node ordering produces positive-volume hexes
- physical tags map back to original material IDs
- mixed, dominant-fraction, tie, and per-material arrays round-trip
- invalid result objects fail without out-of-bounds access
- invalid format and null arguments set useful errors
- write, flush, close, and rename failures return failure
- an interrupted filename export does not replace a previous valid file

### Robustness and performance

- allocation overflow at every size calculation
- allocation-failure injection for each result component
- cancellation during a large sampling operation
- serial/parallel deterministic equivalence
- sanitizer runs for sampling and both exporters
- fuzzed configuration validation and malformed-result export
- performance regression thresholds for one-, two-, and many-material grids

## Recommended delivery order

1. Contract and terminology
2. Coherent results, ties, and complete material discovery
3. Export correctness and composition round-trip
4. Scaling and optional-memory work
5. World-space automatic bounds
6. Adaptive composition estimation
7. Adaptive grid representation

Do not begin boundary-conforming meshing inside `alea_mesh`. Once the structured
grid work is stable, surface extraction and conforming volume meshing should be
designed as separate consumers of the CSG kernel.
