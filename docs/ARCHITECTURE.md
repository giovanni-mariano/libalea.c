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
├── instance_cache    Cache for materialized universe instances
├── spatial_index     KD-tree for fast cell instance queries
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
            int mcnp_surface_id;
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
    int mcnp_surface_id;
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

For the first point query, a spatial index (KD-tree over cell instances) is built lazily. Subsequent queries use it to narrow down candidate cells, avoiding a linear scan of all cells.

The spatial index is particularly important for models with deep universe hierarchies, where the naive approach would test cells at every level.

## How Export Works

Export (`alea_export_mcnp` or `alea_export_openmc`) does:

1. **Build canonical surface map**: assign each primitive a unique export surface ID. If dedup is enabled, identical primitives share an ID.

2. **Handle macrobodies**: depending on `surface_policy`:
   - `ALEA_EMIT_MACROBODY`: emit RCC, BOX, etc. as-is
   - `ALEA_EMIT_SURFACES`: decompose each macrobody into constituent primitives and rewrite the cell's CSG tree

3. **Emit surfaces**: write each unique surface with its coefficients.

4. **Emit cells**: for each cell, walk its CSG tree and emit the boolean expression. Leaf nodes are rewritten using the export surface ID and the effective sense (accounting for dedup and the `inverted` flag).

5. **Emit data cards**: materials, transforms, and other MCNP data cards.

The key invariant: the `inverted` flag is never applied to surface coefficients during output. Surfaces are always emitted with their canonical coefficients. The sense in the cell expression absorbs the inversion.

## Lazy Universe Instantiation

When you load a model with universe fills, Alea does NOT immediately materialize every universe instance. A model with 1000 uses of universe 5 would need 1000 copies of every primitive and node in universe 5 — that's expensive and usually unnecessary.

Instead, Alea uses lazy instantiation:

- **At load time**: fill parameters are stored on the cell. Nothing is expanded.
- **At query time**: the point is inverse-transformed and the query descends into the original universe. No copies needed.
- **At flatten time** (explicit `alea_flatten()`): instances are materialized. Each primitive is cloned and transformed to global coordinates. The instance cache prevents re-materializing the same universe+transform pair.

The instance cache (`alea_instance_cache_t`) maps `(universe_id, transform_id)` to a materialized instance. Each instance stores:

- A primitive remap table (old primitive IDs to new ones)
- A node remap table
- The set of cloned cell roots

## Ray Tracing

Ray tracing (`alea_raycast`) finds all cells intersected by a ray, in order:

### Global approach

1. Build a BVH (bounding volume hierarchy) over all surface bounding boxes
2. Traverse the BVH to find candidate ray-surface intersections
3. Sort intersections by distance
4. Deduplicate hits that share both the same distance and surface ID (but keep different surfaces at the same distance — these represent triple junctions)
5. Walk along the ray, at each intersection testing which cell the point is in
6. Build segments: consecutive regions of same-cell occupancy

The BVH traversal uses a bounded stack (128 entries, supporting trees up to ~64 levels deep). A warning is logged if the stack overflows.

### Cell-aware approach

`alea_raycast_cell_aware` is faster for large models:

1. Find which cell the ray origin is in
2. For that cell, test only its surfaces for intersection
3. Find the nearest intersection
4. Step past it, find the next cell
5. Repeat until t_max

This avoids testing all surfaces globally. The per-cell surface list is built by `alea_build_cell_surface_index`, which walks each cell's CSG tree and collects referenced surface indices.

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

The spatial index (`alea_spatial_index_t`) is a KD-tree that stores **cell instances** — combinations of (cell, transform, depth) that represent a specific materialized position of a cell in global coordinates.

Building the index:

1. Walk the universe hierarchy recursively
2. For each terminal cell (one with a material, not a fill), record its bounding box in global coordinates
3. For each container cell, record it too (needed for fill descent)
4. Build a KD-tree over all instances

Querying:

1. Find all instances whose bounding box contains the query point
2. For each candidate, evaluate the cell's CSG tree (with inverse transform if needed)
3. Return the best match (innermost, deepest in hierarchy)

The spatial index avoids linear scans over cells. On a model with 10,000 cells, a point query touches maybe 5-10 cells instead of all 10,000.

## Cell Adjacency

Cell adjacency is an optimization for grid-based queries (slicing). The idea: if you know pixel (i, j) is in cell A, then pixel (i+1, j) is probably in cell A or one of its neighbors.

Building adjacency:

1. For each cell, find which surfaces bound it (from the CSG tree)
2. For each surface, find which other cells use the same surface
3. Two cells sharing a surface are neighbors

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

## Module Separation

The library is split into modules that can be linked independently:

- **Core** (`libalea.a`): CSG evaluation, tree operations, export framework, dedup. No format-specific code.
- **MCNP** (`libalea_mcnp.a`): MCNP lexer, parser, cell/surface conversion, MCNP export formatting.
- **OpenMC** (`libalea_openmc.a`): OpenMC XML parsing and export.
- **Raycast** (`libalea_raycast.a`): BVH, ray-primitive intersection, segment building.
- **Slice** (`libalea_slice.a`): Grid queries, analytical curve extraction, label positioning.

The `libalea_full.a` archive combines all modules. Use it unless binary size matters.

Module boundaries follow a dependency rule: MCNP and OpenMC depend on Core, but not on each other. Raycast and Slice depend on Core but not on MCNP or OpenMC. Core depends on nothing except the util/ layer.
