# Architecture

This document explains how Alea works internally. Read this if you want to contribute, embed the library, or understand why it makes the decisions it does.

## The Big Picture

Alea is a single-struct library. One `alea_system_t` holds everything: all nodes, all primitives, all surfaces, all cells, all materials, all transforms, all acceleration structures. You create a system, fill it with data (either by loading a file or building programmatically), and then query it.

There is no global state except a single interrupt flag and an error message buffer (both thread-local where it matters). Two systems are completely independent. You can load two MCNP files into two systems and query them concurrently.

## The Data Model

### alea_system_t

The central structure. Here's what it contains, conceptually:

```
alea_system_t
├── nodes[]           Flat array of CSG tree nodes
├── primitives[]      Flat array of deduplicated geometric primitives
├── primitive_index   Hash table for O(1) dedup lookup
├── surfaces[]        Surface entries (MCNP surface ID -> primitive + sense nodes)
├── cells[]           Cell entries (MCNP cell ID -> root node + metadata)
├── cell_index        Hash table for O(1) cell lookup by ID
├── materials[]       Material definitions (nuclide compositions)
├── mixtures[]        Material mixtures
├── transforms[]      Transform matrices (TR cards, TRCL, FILL transforms)
├── universes[]       Universe groupings (built on demand)
├── spatial_index     BVH over cell instances for fast spatial queries
├── surface_bvh       Bounding volume hierarchy for ray tracing
├── config            Tolerance, export options, void parameters
├── stats             Dedup hits, memory usage, conversion counts
└── error state       Last error code and message
```

All arrays are dynamic vectors (`alea_vec`) that grow automatically. They start small and double when full. You never need to pre-size anything.

### Nodes

A CSG tree node (`alea_node_t`) is either a leaf (primitive) or an internal node (boolean operation):

```c
struct alea_node {
    uint32_t type_and_flags;   // operation type in lower 8 bits
    alea_material_id_t material_id;
    alea_bbox_t bbox;

    union {
        struct {               // leaf node
            alea_primitive_id_t primitive_id;
            alea_primitive_type_t prim_type;
            int8_t sense;      // +1 or -1
            int8_t inverted;   // canonicalization flag
            int mc_surface_id;
        } primitive;

        struct {               // internal node
            alea_node_id_t left;
            alea_node_id_t right;
        } operation;
    };
};
```

Nodes are stored in a flat array and referenced by `alea_node_id_t` (a `uint32_t` index). This is critical: because nodes are referenced by index, the array can grow (reallocate) without invalidating any references. This eliminates an entire class of use-after-realloc bugs.

Each cell's geometry is a tree rooted at `cell->root_node_id`. To evaluate "is point P inside cell C?", you recursively evaluate the tree: leaves test the primitive, internal nodes combine results with AND/OR/NOT.

### Primitives

Primitives (`alea_primitive_entry_t`) are the geometric objects: planes, spheres, cylinders, etc. They are stored separately from nodes because multiple nodes can reference the same primitive.

For example, surface 5 in MCNP creates two nodes: one with `sense = -1` (inside) and one with `sense = +1` (outside). Both point to the same primitive. The primitive stores only the geometric parameters (center, radius, coefficients). The sense is on the node.

Primitives carry a `ref_count` — the number of nodes pointing at them. This is used for diagnostics, not for lifecycle management (primitives live as long as the system lives).

### Surfaces

A surface entry (`alea_surface_entry_t`) maps an MCNP surface ID to a primitive and pre-allocates both sense nodes:

```c
struct alea_surface_entry {
    int mc_surface_id;
    alea_primitive_id_t primitive_id;
    alea_node_id_t pos_node;           // +S (sense = +1)
    alea_node_id_t neg_node;           // -S (sense = -1)
    int transform_id;
    alea_boundary_type_t boundary_type;
    // ... macrobody expansion nodes
};
```

When an MCNP cell card references `-5`, the parser looks up surface 5's `neg_node`. This pre-allocation means the parser never needs to check if a node already exists for a given surface + sense combination.

### Cells

A cell entry (`alea_cell_entry_t`) is the richest structure. It holds:

- The MCNP cell ID and its index in the array
- The root node of the CSG tree
- Material, density, universe
- Fill information (universe ID, transform)
- Lattice parameters (type, pitch, fill array, dimensions)
- Cell importance values (IMP:N, IMP:P, etc.)
- A list of surface indices (for cell-aware ray tracing)
- Neighbor list (for adjacency-based traversal)

Cells are looked up by ID through a hash table (`cell_index`), giving O(1) access.

## Primitive Deduplication

Large models contain many duplicate surfaces. A tokamak model might define the same cylindrical surface dozens of times with slightly different floating-point coefficients (because different CAD tools or human authors wrote them independently).

For the full pipeline — canonicalization, hashing, tolerance-based equality, canonical surface maps, and the XOR sense correction formula — see [Surface Deduplication](SURFACE_DEDUP.md).

Alea deduplicates automatically:

1. **Canonicalize**: Normalize the primitive so the first non-zero coefficient is positive. This means `0x + 0y - 1z + 5 = 0` and `0x + 0y + 1z - 5 = 0` become the same canonical form. The `inverted` flag on the node records whether the sign was flipped.

2. **Hash**: Compute a hash of the canonical coefficients (with tolerance rounding).

3. **Lookup**: Check the hash table. If a primitive with the same hash exists, compare coefficients within tolerance (`abs_tol` and `rel_tol` from config).

4. **Reuse or insert**: If a match is found, return the existing `primitive_id`. Otherwise, insert the new primitive and return a fresh ID.

The `inverted` flag is the subtle part. When two primitives differ only by sign (one has `(0, 0, 1, -5)` and the other `(0, 0, -1, 5)`), they map to the same canonical primitive, but one of them has `inverted = 1`. The sense on the node is adjusted accordingly, so the final evaluation is correct:

```
effective_sense = sense * (inverted ? -1 : 1)
```

The emitted surface during export always uses the raw canonical coefficients. The `inverted` flag is only used for evaluation.

## How Point Queries Work

`alea_find_cell(sys, x, y, z)` does the following:

### Simple case (no fills, no lattices)

1. Build the universe index if not already done
2. Get the list of cells in universe 0
3. For each cell, evaluate its CSG tree at (x, y, z)
4. Return the first cell whose tree evaluates to `true` (inside)

Tree evaluation is recursive:

- **Primitive leaf**: evaluate the primitive's signed distance function (SDF) at (x, y, z). If `sense * SDF(x, y, z) < 0`, the point is inside.
- **Intersection**: left is inside AND right is inside
- **Union**: left is inside OR right is inside
- **Difference**: left is inside AND right is NOT inside
- **Complement**: NOT (child is inside)

Bounding boxes provide early rejection: if the point is outside the node's bounding box, skip it.

### With fills

When a cell has `FILL=N`:

1. The point matches the container cell's region
2. Get the fill transform (if any)
3. Apply the inverse transform to the point: `local = inverse_transform(x, y, z)`
4. Search for `local` among the cells of universe N
5. If that cell also has a fill, recurse

This is `alea_find_cell_recursive()`. The depth is bounded by the universe nesting depth (typically 3-5 levels in real models).

### With lattices

When a cell has `lat_type != 0`:

1. The point matches the lattice cell's region
2. Compute which lattice element contains the point:
   - **Rectangular**: `i = floor((x - lower_left_x) / pitch_x)`, same for j, k
   - **Hexagonal**: convert to fractional hex coordinates, round to nearest hex center using cube coordinate constraints
3. Look up the universe ID in `lat_fill[i][j][k]`
4. Translate the point to the element's local frame
5. Search for the point in that universe

Element computation is O(1) — no iteration over elements. A lattice with 100,000 elements takes the same time as one with 4.

### Spatial index

For the first point query, a spatial index (BVH over cell instances) is built lazily. Subsequent queries use it to narrow down candidate cells, avoiding a linear scan of all cells.

The spatial index is particularly important for models with deep universe hierarchies, where the naive approach would test cells at every level. See [Spatial Index](#spatial-index) for the full design.

## How Export Works

Export (`mcnp_export` or `openmc_export`) does:

1. **Build canonical surface map**: assign each primitive a unique export surface ID. If dedup is enabled, identical primitives share an ID.

2. **Handle macrobodies**: depending on `surface_policy`:
   - `ALEA_EMIT_MACROBODY`: emit RCC, BOX, etc. as-is
   - `ALEA_EMIT_SURFACES`: decompose each macrobody into constituent primitives and rewrite the cell's CSG tree

3. **Emit surfaces**: write each unique surface with its coefficients.

4. **Emit cells**: for each cell, walk its CSG tree and emit the boolean expression. Leaf nodes are rewritten using the export surface ID and the effective sense (accounting for dedup and the `inverted` flag).

5. **Emit data cards**: materials, transforms, and other MCNP data cards.

The key invariant: the `inverted` flag is never applied to surface coefficients during output. Surfaces are always emitted with their canonical coefficients. The sense in the cell expression absorbs the inversion.

## Lazy Universe Evaluation

When you load a model with universe fills, Alea does NOT immediately materialize every universe instance. A model with 1000 uses of universe 5 would need 1000 copies of every primitive and node in universe 5 — that's expensive and usually unnecessary.

Instead:

- **At load time**: fill parameters are stored on the cell. Nothing is expanded.
- **At query time**: the point is inverse-transformed and the query descends into the original universe. No copies needed.
- **At flatten time** (explicit `alea_flatten()`): instances are materialized. Each primitive is cloned and transformed to global coordinates.

## Ray Tracing

Ray tracing (`alea_raycast`) finds all cells intersected by a ray, in order:

### Global approach

1. Build a BVH (bounding volume hierarchy) over all surface bounding boxes
2. Traverse the BVH to find candidate ray-surface intersections
3. Sort intersections by distance
4. Deduplicate hits that share both the same distance and surface ID (but keep different surfaces at the same distance — these represent triple junctions)
5. Walk along the ray, at each intersection testing which cell the point is in
6. Build segments: consecutive regions of same-cell occupancy

The BVH traversal uses a bounded stack (128 entries, supporting trees up to ~64 levels deep). A warning is logged if the stack overflows. See [Surface BVH](#surface-bvh) for details on the construction and traversal algorithm.

### Cell-aware approach

`alea_raycast_cell_aware` is semantically equivalent to `alea_raycast`. For
non-lattice models it is faster for large models:

1. Find which cell the ray origin is in
2. For that cell, test only its surfaces for intersection
3. Find the nearest intersection
4. Step past it, find the next cell
5. Repeat until t_max

This avoids testing all surfaces globally. The per-cell surface list is built by `alea_build_cell_surface_index`, which walks each cell's CSG tree and collects referenced surface indices. For lattice models the cell-aware entry point uses the canonical DDA-aware path so synthetic lattice element-boundary hits are preserved.

### Lattice DDA

For rays traversing lattices, a DDA (Digital Differential Analyzer) algorithm steps through lattice elements:

- **Rectangular lattices**: Amanatides-Woo DDA. Compute `t_step` and `t_next` for each axis. At each step, advance to the nearest boundary.
- **Hexagonal lattices**: 3-axis oblique DDA using the hex coordinate system.

At each element boundary, a synthetic hit is emitted (`surface_id = 0`) so the segment builder knows the ray crossed into a new lattice element.

### Ray-Primitive Intersections

Ray-surface intersection routines respect all geometric constraints stored on the primitive:

- **Cylinders**: hits outside `[axis_min, axis_max]` are rejected
- **Cones**: hits outside axis bounds are rejected, and `sheet_selection` restricts to a single nappe (positive or negative side of the apex) when non-zero
- **Normals**: computed analytically for planes, spheres, cylinders, cones (gradient of implicit equation), boxes (closest face), quadrics (gradient), and tori (gradient of quartic implicit). Other types fall back to `(0,0,1)`.

## Spatial Index

**Files:** `src/core/alea_spatial.h`, `src/core/alea_spatial.c`

The spatial index (`alea_spatial_index_t`) is a BVH (bounding volume hierarchy) that stores **cell instances** — combinations of (cell, transform, depth) that represent a specific materialized position of a cell in global coordinates.

**Why a BVH over instances, not over raw geometry:** a model with deep universe hierarchies and thousands of fill placements would require flattening every surface into global coordinates to build a traditional spatial index. That is expensive in both time and memory, and destroys the compact representation of shared universes. Instead, the spatial index flattens only the **bounding boxes**: for each cell instance, it records `(cell_index, transform, global_bbox)` — lightweight metadata. The actual CSG geometry stays in its original universe-local form and is evaluated on demand through the stored transform. This keeps the index small while still enabling O(log N) spatial queries.

### Building the index

1. Walk the universe hierarchy recursively (`collect_instances_recursive`)
2. For each terminal cell (one with a material, not a fill), compute its bounding box in global coordinates by transforming the local bbox through the accumulated transform chain, then clip against the parent cell's bbox
3. For each container cell, record it too (needed for fill descent during point queries)
4. Build a BVH over all instances using **median splits** on centroids (leaf size = 4, max depth = 30)

**Why median splits instead of SAH:** the spatial index is built over cell instances whose bounding boxes often overlap heavily (nested universes, lattice elements). SAH cost estimation assumes non-overlapping primitives and provides little benefit here, while median splits are simpler and faster to build — important because the index is rebuilt whenever the model changes.

The build is **thread-safe**: an atomic CAS on `spatial_build_state` (0=pending → 1=building → 2=done) ensures only one thread builds while others spin-wait. The index is lazy-built on first query if not built explicitly.

### Querying

Four query types are available:

- **Region query** (`alea_spatial_query_region`): find all terminal instances whose bbox overlaps a query bbox. Used by the slice module to find cells intersecting a cutting plane.
- **Z-slice query** (`alea_spatial_query_slice_z`): convenience wrapper that converts a Z-plane slice into a thin bbox region query.
- **Point query** (`alea_spatial_query_point`): find all instances (including non-terminal) at a point, sorted by depth. Used for full hierarchy traversal.
- **Point-in-cell with coherence** (`alea_spatial_find_cells_at_point`): high-level query with a thread-local coherence cache. Caches the last found cell path so consecutive nearby queries (adjacent pixels in a grid) reuse the result without re-traversing the BVH.

The spatial index avoids linear scans over cells. On a model with 10,000 cells, a point query touches maybe 5-10 cells instead of all 10,000.

### Role in slice plotting

The spatial index is central to both slice operations:

1. **Curve computation** (`alea_compute_slice_curves_spatial`): converts the slice plane into a thin bounding box, queries the BVH to find which cell instances intersect, deduplicates hits by `(cell_index, transform)`, then computes analytical curve intersections only for those cells. Without the spatial index, every cell in the model would need to be tested.

2. **Grid-based cell identification** (`alea_find_cells_grid`): uses `alea_spatial_find_cells_at_point` per pixel with the coherence cache. Combined with per-universe adjacency walking, this gives interactive-speed slicing even on models with hundreds of thousands of cells.

## Surface BVH

**Files:** `src/raycast/bvh.h`, `src/raycast/bvh.c`

The surface BVH (`alea_bvh_t`) is a separate acceleration structure from the spatial index. It indexes **primitive surface bounding boxes** and is used exclusively for **ray-surface intersection** queries.

**Why a separate BVH for surfaces:** ray tracing and spatial queries have fundamentally different access patterns. Ray tracing needs to find which surfaces a ray crosses, ordered by distance. Spatial queries need to find which cells contain a point or overlap a region. Trying to use one structure for both would either degrade ray traversal performance (because cell instances can deeply overlap) or fail at spatial queries (because surface BVHs have no concept of cell containment). Two purpose-built structures, each optimal for its query type, are simpler and faster than one compromise.

### Construction

The surface BVH uses **SAH (Surface Area Heuristic)** with 16 bins per axis:

1. For each surface in `sys->surfaces`, compute its bounding box
2. Recursively partition using SAH cost estimation: for each axis, bin the surfaces by centroid position, evaluate `cost = traversal_cost + (SA_left * N_left + SA_right * N_right) / SA_parent * intersect_cost`, and pick the split with minimum cost
3. Stop splitting when leaf size ≤ 4 or depth exceeds 30

**Why SAH here but not for the spatial index:** surface primitives in a ray tracing context have relatively tight, non-overlapping bounding boxes. SAH excels in this regime — it produces trees where rays traverse fewer nodes by putting large surfaces in large nodes and grouping spatially coherent surfaces together. The build cost (O(N log N) with binning) is acceptable because the BVH is built once and reused for many ray queries.

The BVH nodes are **64 bytes** (cache-line-friendly): bbox + child indices + leaf metadata. Surface indices are reordered for cache locality during traversal.

### Traversal

`alea_bvh_traverse` does **stack-based** (not recursive) traversal:

- Precomputes inverse ray direction and sign bits
- Uses a 128-entry stack (supports ~64 tree levels)
- Pushes the far child first so the near child is processed next
- Calls a user callback for each candidate surface whose bbox the ray hits
- Falls back to linear scan (`raycast_surfaces_linear`) if BVH is disabled or the build failed

### Lazy build and invalidation

The surface BVH is cached on `sys->surface_bvh` and lazy-built on first raycast via `alea_raycast_ensure_caches()`. The `bvh_dirty` flag triggers a rebuild when surfaces are added or modified. If the surface count changes since the last build, the BVH is also rebuilt.

## Cell Adjacency

Cell adjacency is an optimization for grid-based queries (slicing). The idea: if you know pixel (i, j) is in cell A, then pixel (i+1, j) is probably in cell A or one of its neighbors.

### Building adjacency

1. For each cell, find which surfaces bound it (from the CSG tree)
2. For each surface, find which other cells use the same surface
3. Two cells sharing a surface are neighbors
4. Deduplicate using a bitset (one 64-bit word array, ~14.5 KB for 116K cells) with dirty-word tracking for O(neighbors) reset per cell

**Why a 128-cell threshold:** surfaces shared by more than 128 cells (`ADJACENCY_MAX_CELLS_PER_SURFACE`) are skipped when building adjacency. Global planes (e.g., a bounding plane used by thousands of cells) create O(n²) cliques — every cell sharing that plane becomes a neighbor of every other. These cliques waste memory without helping local adjacency walking, because a global plane doesn't imply geometric proximity. Skipping high-fanout surfaces keeps the neighbor lists small and focused on actual geometric neighbors.

### Memory layout

All neighbor data uses a **pool allocator** (`sys->neighbor_pool`): a single `malloc` for the entire adjacency structure. This avoids heap fragmentation from thousands of small per-cell neighbor arrays. `alea_destroy` and `alea_reset` check `neighbor_pool` before attempting per-cell `free(neighbors)`.

### Usage in grid queries

During grid queries, after finding the cell for one pixel, the next pixel first checks the same cell and its neighbors before falling back to a full search. This gives ~10x speedup on grid queries because consecutive pixels almost always share a cell or neighbor.

Only cells in the same universe are considered neighbors. This prevents cross-universe false matches.

## Memory Management

Alea uses a simple ownership model:

- `alea_system_t` owns all its arrays (nodes, primitives, surfaces, cells, etc.)
- `alea_destroy()` frees everything. `alea_reset()` frees all per-cell data (surface indices, neighbors, lattice fills) before clearing the vectors.
- Opaque result objects (raycast results, void results, slice curves) must be freed by the caller using the corresponding `_free()` function. Raycast result objects are safe to reuse — `alea_raycast()` frees prior allocations internally before populating new results.

Internally, the MCNP parser uses an arena allocator for temporary parse data. The arena is freed after conversion to the `alea_system_t`.

Dynamic arrays (`alea_vec`) start with a small capacity and double when full. This amortizes reallocation cost to O(1) per insertion. Node IDs (not pointers) ensure that reallocation doesn't invalidate existing references.

## Error Handling

Two mechanisms:

1. **Public API** (`alea_error()`, `alea_error_code()`): Simple, user-facing. Functions return NULL or negative on error; call `alea_error()` for details.

2. **Internal** (`alea_set_error_detail()`): Thread-local error state with printf-style formatting. Used by internal functions to propagate detailed error messages up the call chain.

## Interrupt Support

Long-running operations (loading large files, grid queries, ray tracing) check a global `volatile sig_atomic_t` flag periodically. The `ALEA_CHECK_INTERRUPTED(retval)` macro does this check and returns early if the flag is set.

The Python binding installs a SIGINT handler that sets this flag. After the operation returns `ALEA_ERR_INTERRUPTED`, the binding calls `alea_clear_interrupt()` and raises `KeyboardInterrupt`.

The flag is signal-safe: it's a `sig_atomic_t`, written atomically, checked with a simple comparison. No mutexes, no race conditions.

## Library Structure

The library is split into a core library and optional format modules:

- **Core** (`libalea.a`): The complete geometry engine — CSG evaluation, tree operations, export framework, dedup, primitives, bounding boxes, raycast (BVH, ray-primitive intersection, segment building), slice (grid queries, analytical curve extraction, label positioning), render (3D batch renderer with Phong shading, shadow rays, cutaway views, OpenMP parallelized), and mesh export (structured hexahedral mesh sampling, Gmsh/VTK output).
- **MCNP** (`libalea_mcnp.a`): Optional module. MCNP lexer, parser, cell/surface conversion, MCNP export formatting.
- **OpenMC** (`libalea_openmc.a`): Optional module. OpenMC XML parsing and export.
- **Nuclear Data** (`libalea_nucdata.a`): Optional module. ACE-format nuclear data reader, cross-section lookup, reaction classification, Doppler broadening, multigroup collapse.

Link against `libalea.a` for everything except format/physics modules. Add the format modules only if you need to read or write those formats. Add `libalea_nucdata.a` for nuclear data and transport physics. For convenience, `libalea_full.a` bundles everything (core + MCNP + OpenMC) into a single archive.

The **Lua binding** layer (`src/lua_bind/`) wraps the public C API for the interactive CLI (`bin/alea`). It is compiled into the CLI binary, not shipped as a separate library.

The source is organized by concern (`src/core/`, `src/primitives/`, `src/raycast/`, `src/slice/`, `src/render/`, `src/mesh/`, `src/util/`) but all compile into the single `libalea.a`. Format modules (`src/mcnp/`, `src/openmc/`) and the nuclear data module (`src/nucdata/`) produce separate archives.

## Nuclear Data Module

The nuclear data module (`libalea_nucdata.a`) reads ACE nuclear data files and provides cross-section lookups, reaction classification, and data decoding. It sits between raw evaluated nuclear data and the geometry engine.

```
  xsdir/xsdata          ACE files (.c, .p)
       |                      |
       v                      v
  +---------------+     +--------------+
  | alea_nuc_xsdir|----->| ace_reader  |
  +---------------+     +--------------+
       |                      |
       v                      v
  +-----------------+    +------------+
  | load_nuclide()  |--->| xs_decode  |
  +-----------------+    +------------+
                               |
              +----------+     |
              | nuclide  |<----+
              +----------+
              /    |    \
             v     v     v
         lookup  reaction  material
```

### Source Organization

| File | Role |
|------|------|
| `alea_nucdata.h` | Public API header |
| `alea_nucdata_types.h` | Type definitions — structs, enums, constants |
| `nuclear_internal.h` | Internal helpers shared between source files (not public) |
| `constants.h` | Physical and algorithmic constants |
| `context.c` | ZAID parsing, nuclide free |
| `ace_reader.c` | Read ACE files — Type 1 (ASCII) and Type 2 (binary) |
| `xsdir.c` | Parse xsdir/xsdata directory files |
| `xs_decode.c` | Decode XSS array into nuclide structure: ESZ, SIG, NU, URR, photon blocks |
| `angular.c` | Decode angular distributions from ACE AND block |
| `energy_dist.c` | Decode energy distributions from ACE DLW block (Laws 3, 4, 7, 9, 11, 44, 61, 66) |
| `lookup.c` | Energy grid binary search, cross-section interpolation, URR factors, photon XS |
| `reaction.c` | Reaction classification, neutron yield (TYR), ν̄ evaluation |
| `material.c` | Material composition, macroscopic cross sections, nuclide/reaction sampling |
| `multigroup.c` | Collapse pointwise data to multigroup constants |
| `doppler.c` | On-the-fly Doppler broadening with exact kernel |

### Data Pipeline

**1. Loading: xsdir → ACE → nuclide.** The user first loads an xsdir file, which maps ZAID strings like `"92235.80c"` to ACE file paths and byte offsets. Then `alea_nuc_load_nuclide` reads the ACE file and decodes the physics.

ACE tables use Fortran-style 1-based indexing throughout. The NXS array (16 ints) gives table dimensions, the JXS array (32 ints) gives block locators into the XSS data array. The internal `xss()` helper converts 1-based indices to 0-based with bounds checking.

**2. Decode: XSS blocks.** `xs_decode.c` extracts physics blocks — ESZ (energy grid + principal cross sections), SIG (non-elastic reactions), NU (fission nu-bar), photon data, and URR probability tables. After decoding, an MT-to-reaction lookup table is built for O(1) reaction finding. Angular and energy distributions are decoded lazily on first collision sample.

**3. Runtime: lookups and sampling.** All cross-section lookups start with a binary search on the nuclide's energy grid. Reaction sampling does one O(log N) search then walks all reaction cross sections using direct sub-grid indexing — no per-reaction binary search. This is the hot path in transport.

### Key Design Decisions

- **Flat arrays, not trees**: cross sections are flat `double*` arrays indexed by the energy grid. Cache-friendly sequential access, simple interpolation, no pointer chasing.
- **MT lookup table**: 1000-element table maps MT directly to reaction index. 4 KB cost, O(1) lookup.
- **Log-space photon storage**: pre-stored `ln(energy)` and `ln(sigma)` arrays. Runtime interpolation is linear in log space + one `exp()` — ~3x faster than naive log-log.
- **Bounds-checked XSS access**: all access goes through `xss()`, `xss_int()`, `xss_copy()` helpers that validate indices.
- **User-owned objects**: no hidden context or cache. The user loads an xsdir, loads nuclides from it, and frees them when done. Multiple xsdir files can coexist.

### Performance

- **Single binary search** per collision event for reaction sampling
- **Sub-grid indexing**: each reaction's threshold offset converts main-grid index to sub-array index without additional searches
- **Doppler broadening**: O(N × W) where W is the local integration window width (~50-200 points at room temperature)
- **Memory**: loaded nuclide dominated by XSS array (0.5-5 MB for neutron tables). Decoded structures add ~30-50% overhead.

## See Also

- [Concepts](CONCEPTS.md) — domain concepts (surfaces, sense, cells, universes, lattices, materials)
- [Surface Deduplication](SURFACE_DEDUP.md) — full dedup pipeline: canonicalization, hashing, tolerance equality, canonical surface maps, XOR sense correction
- [API Reference](API.md) — complete public function listing
- [Tutorial](TUTORIAL.md) — C API walkthrough with working examples
- [Lua Tutorial](LUA_TUTORIAL.md) — Lua binding reference and examples
