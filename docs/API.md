# API Reference

Every public function in Alea, grouped by what you're trying to do.

**Headers**: `alea.h` (main), `alea_types.h` (types and constants), `alea_raycast.h` (ray tracing), `alea_slice.h` (visualization), `alea_render.h` (3D rendering), `alea_mesh.h` (mesh export), `alea_mcnp.h` (MCNP I/O), `alea_openmc.h` (OpenMC I/O), `alea_nucdata.h` (nuclear data).

**Conventions**:
- Functions returning pointers return `NULL` on error
- Functions returning `int` return 0 on success, negative on error
- After any error, `alea_error()` gives a human-readable message
- Node IDs, cell IDs, and surface IDs are all 0-based internally but MCNP-numbered externally

---

## Version

### alea_version

```c
const char* alea_version(void);
```

Returns a static string like `"0.1.0"`.

---

## Error Handling

### alea_error

```c
const char* alea_error(void);
```

Returns the last error message, or an empty string if no error.

### alea_error_code

```c
int alea_error_code(void);
```

Returns the last error code.

### alea_error_clear

```c
void alea_error_clear(void);
```

Clears the last error.

### alea_error_string

```c
const char* alea_error_string(alea_error_t error);
```

Returns a static description string for an error code. Defined in `alea_types.h`.

### alea_set_error_detail

```c
void alea_set_error_detail(alea_error_t code, const char* fmt, ...);
```

Set a detailed error message (thread-local) with printf-style formatting.

### alea_get_error_detail

```c
const char* alea_get_error_detail(void);
```

Get detailed error message. Returns the detailed message if set, otherwise the static string from the error code.

### alea_get_last_error

```c
alea_error_t alea_get_last_error(void);
```

Get the current error code (thread-local).

### alea_clear_error_detail

```c
void alea_clear_error_detail(void);
```

Clear the error detail.

---

## Interrupt Support

Cooperative interrupt mechanism for cancelling long-running operations. The Python binding's SIGINT handler calls `alea_interrupt()`, then the next iteration of any long-running loop checks the flag and returns `ALEA_ERR_INTERRUPTED`.

### alea_interrupt

```c
void alea_interrupt(void);
```

Request interruption. Signal-safe: can be called from a signal handler.

### alea_clear_interrupt

```c
void alea_clear_interrupt(void);
```

Clear the interrupt flag. Call after handling an interrupted operation.

### alea_interrupted

```c
bool alea_interrupted(void);
```

Check if an interrupt has been requested.

---

## System Lifecycle

### alea_create

```c
alea_system_t* alea_create(void);
```

Create an empty system. Use this when building geometry programmatically. Returns NULL on allocation failure.

### alea_destroy

```c
void alea_destroy(alea_system_t* sys);
```

Free all memory associated with a system. Safe to call with NULL.

### alea_clone

```c
alea_system_t* alea_clone(const alea_system_t* sys);
```

Deep-copy a system. The clone is fully independent.

### alea_reset

```c
void alea_reset(alea_system_t* sys);
```

Clear all data but keep allocated memory. Faster than destroy + create if you're reusing a system.

---

## Configuration

### alea_get_config

```c
alea_config_t alea_get_config(const alea_system_t* sys);
```

Get the current configuration. Start from the returned config or `ALEA_CONFIG_DEFAULT`, tweak fields, then call `alea_set_config()`.

### alea_set_config

```c
void alea_set_config(alea_system_t* sys, const alea_config_t* config);
```

Set the configuration.

**`alea_config_t` fields** (defaults in parentheses):

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `abs_tol` | double | 1e-14 | Absolute tolerance |
| `rel_tol` | double | 1e-12 | Relative tolerance |
| `zero_threshold` | double | 1e-14 | Zero threshold |
| `dedup` | bool | true | Deduplicate primitives |
| `log_level` | int | 2 (WARN) | Log verbosity (0=none..5=trace) |
| `export_materials` | bool | true | Export material data |
| `export_transforms` | bool | true | Export transforms |
| `universe_depth` | int | -1 | Universe depth for export (-1=all) |
| `fill_depth` | int | 0 | Fill expansion depth (0=none) |
| `void_max_depth` | int | 8 | Octree depth for void generation |
| `void_min_size` | double | 0.1 | Minimum void region size |
| `void_probes_per_axis` | int | 3 | Probes per axis for void detection |
| `merge_cell_weight` | double | 1.0 | Weight for cell count in merge cost |
| `merge_surface_weight` | double | 0.1 | Weight for surface count in merge cost |
| `merge_max_surfaces` | int | 24 | Maximum surfaces per merged void cell |
| `merge_min_cells` | int | 1 | Minimum cells before merging stops |
| `merge_use_greedy` | bool | false | Use greedy merge algorithm |
| `void_consolidate` | int | 100 | Consolidation threshold (0=off) |
| `flatten_max_depth` | int | 0 | Max flatten depth (0=unlimited) |

---

## Loading Models

Loading functions live in format-specific headers (`alea_mcnp.h`, `alea_openmc.h`). See [MCNP Module](#mcnp-module-alea_mcnph) and [OpenMC Module](#openmc-module-alea_openmch).

---

## Indexing

### alea_build_universe_index

```c
int alea_build_universe_index(alea_system_t* sys);
```

Build the universe hierarchy index.

### alea_prepare_query_acceleration

```c
int alea_prepare_query_acceleration(alea_system_t* sys);
```

Build query acceleration (the hierarchical spatial index, surface BVH, and
related caches) explicitly before point, spatial, raycast, slice, render, or
mesh queries. Query paths do not lazily build shared caches. Geometry mutation
invalidates prepared caches; call this again after adding/removing/changing
surfaces, transforms, or cells.

### alea_query_acceleration_stats

```c
int alea_query_acceleration_stats(const alea_system_t* sys,
                                  alea_query_acceleration_stats_t* out_stats);
```

Inspect the currently built query-acceleration structure without building
caches. Call `alea_prepare_query_acceleration()` first when a built index is
required. The hierarchical fields report the TLAS/BLAS shape, including
universe, BLAS, placement, transform, and memory counts.

---

## Geometry Queries

### alea_find_cell

```c
int alea_find_cell(const alea_system_t* sys, double x, double y, double z);
```

Legacy convenience wrapper for single-point queries. Prefer `alea_find_cell_at()` for resolved point answers, or `alea_find_all_cells()` when hierarchy detail is needed.

Returns the internal cell index, or `-1` if not found.

### alea_find_all_cells

```c
int alea_find_all_cells(const alea_system_t* sys, double x, double y, double z,
                        alea_cell_hit_t* hits, size_t max_hits);
```

Find all cells containing a point (useful for overlap detection). Returns number of hits, or -1 on error.

### alea_point_inside

```c
bool alea_point_inside(const alea_system_t* sys, alea_node_id_t node,
                       double x, double y, double z);
```

Test if a point is inside a CSG node's region.

### alea_material_at

```c
int alea_material_at(const alea_system_t* sys, double x, double y, double z);
```

Legacy convenience wrapper for single-point queries. Prefer `alea_find_cell_at()`.

Returns material ID, `0` for void, or `-1` on error.

### alea_find_overlaps

```c
int alea_find_overlaps(const alea_system_t* sys, int* pairs, size_t max_pairs);
```

Find overlapping cell pairs by random sampling. `pairs` is filled as `[cell1, cell2, cell1, cell2, ...]`. Returns number of pairs found.

### alea_find_cell_at

```c
int alea_find_cell_at(const alea_system_t* sys, double x, double y, double z,
                      int* out_cell_id, int* out_material);
```

Preferred single-point query API. Finds the resolved cell and returns both cell ID and material in one call.

### alea_set_debug_trace

```c
void alea_set_debug_trace(int enable);
```

Enable/disable debug tracing for point queries. When enabled, prints detailed traversal information.

---

## CSG Construction - Surfaces

All surface creation functions return the index into `sys->surfaces`, or -1 on error. Pass `surface_id=0` for automatic ID assignment. Use `alea_halfspace()` to get the +/- sense node.

```c
// Example
int idx = alea_sphere_surface(sys, 0, 0.0, 0.0, 0.0, 5.0);
alea_node_id_t interior = alea_halfspace(sys, idx, -1);
alea_node_id_t exterior = alea_halfspace(sys, idx, +1);
```

### Primitive Surfaces

```c
int alea_plane_surface(alea_system_t* sys, int surface_id,
                       double a, double b, double c, double d);
```

Plane: `ax + by + cz + d = 0`.

```c
int alea_sphere_surface(alea_system_t* sys, int surface_id,
                        double cx, double cy, double cz, double r);
```

Sphere centered at `(cx, cy, cz)` with radius `r`.

```c
int alea_cylinder_x_surface(alea_system_t* sys, int surface_id,
                            double cy, double cz, double r);
int alea_cylinder_y_surface(alea_system_t* sys, int surface_id,
                            double cx, double cz, double r);
int alea_cylinder_z_surface(alea_system_t* sys, int surface_id,
                            double cx, double cy, double r);
```

Axis-aligned cylinders with given center offsets and radius.

```c
int alea_box_surface(alea_system_t* sys, int surface_id,
                     double xmin, double xmax,
                     double ymin, double ymax,
                     double zmin, double zmax);
```

Axis-aligned box (RPP).

```c
int alea_cone_x_surface(alea_system_t* sys, int surface_id,
                        double cx, double cy, double cz, double t_squared);
int alea_cone_y_surface(alea_system_t* sys, int surface_id,
                        double cx, double cy, double cz, double t_squared);
int alea_cone_z_surface(alea_system_t* sys, int surface_id,
                        double cx, double cy, double cz, double t_squared);
```

Axis-aligned cones. `t_squared` is tan^2(half-angle).

```c
int alea_torus_x_surface(alea_system_t* sys, int surface_id,
                         double cx, double cy, double cz,
                         double major_radius, double minor_radius);
int alea_torus_y_surface(alea_system_t* sys, int surface_id,
                         double cx, double cy, double cz,
                         double major_radius, double minor_radius);
int alea_torus_z_surface(alea_system_t* sys, int surface_id,
                         double cx, double cy, double cz,
                         double major_radius, double minor_radius);
```

Axis-aligned tori. `major_radius` is distance from axis to tube center.

```c
int alea_quadric_surface(alea_system_t* sys, int surface_id,
                         double A, double B, double C,
                         double D, double E, double F,
                         double G, double H, double I, double J);
```

General quadric: `Ax^2 + By^2 + Cz^2 + Dxy + Eyz + Fxz + Gx + Hy + Iz + J = 0`.

### Macrobody Surfaces

```c
int alea_rcc_surface(alea_system_t* sys, int surface_id,
                     double base_x, double base_y, double base_z,
                     double height_x, double height_y, double height_z,
                     double radius);
```

Right Circular Cylinder. `base` is the center of the bottom cap; `height` is the axis vector.

```c
int alea_box_general_surface(alea_system_t* sys, int surface_id,
                             double corner_x, double corner_y, double corner_z,
                             double v1_x, double v1_y, double v1_z,
                             double v2_x, double v2_y, double v2_z,
                             double v3_x, double v3_y, double v3_z);
```

General oriented box. `corner` is one vertex; `v1`, `v2`, `v3` are the three edge vectors.

```c
int alea_sph_surface(alea_system_t* sys, int surface_id,
                     double cx, double cy, double cz, double r);
```

Sphere macrobody (same geometry as `alea_sphere_surface`, stored as MCNP macrobody format).

```c
int alea_trc_surface(alea_system_t* sys, int surface_id,
                     double base_x, double base_y, double base_z,
                     double height_x, double height_y, double height_z,
                     double base_radius, double top_radius);
```

Truncated Right Cone. A cone frustum (or cylinder if both radii equal).

```c
int alea_ell_surface(alea_system_t* sys, int surface_id,
                     double v1_x, double v1_y, double v1_z,
                     double v2_x, double v2_y, double v2_z,
                     double major_axis_len);
```

Ellipsoid defined by two foci (`v1`, `v2`) and major axis length.

```c
int alea_rec_surface(alea_system_t* sys, int surface_id,
                     double base_x, double base_y, double base_z,
                     double height_x, double height_y, double height_z,
                     double axis1_x, double axis1_y, double axis1_z,
                     double axis2_x, double axis2_y, double axis2_z);
```

Right Elliptical Cylinder. `base` is center of bottom, `height` is axis vector, `axis1`/`axis2` are perpendicular semi-axis vectors.

```c
int alea_wed_surface(alea_system_t* sys, int surface_id,
                     double vertex_x, double vertex_y, double vertex_z,
                     double v1_x, double v1_y, double v1_z,
                     double v2_x, double v2_y, double v2_z,
                     double v3_x, double v3_y, double v3_z);
```

Wedge (triangular prism). `vertex` is the right-angle corner; `v1`/`v2` are the triangle legs; `v3` is the extrusion vector.

```c
int alea_rhp_surface(alea_system_t* sys, int surface_id,
                     double base_x, double base_y, double base_z,
                     double height_x, double height_y, double height_z,
                     double r1_x, double r1_y, double r1_z,
                     double r2_x, double r2_y, double r2_z,
                     double r3_x, double r3_y, double r3_z);
```

Right Hexagonal Prism. `base` is bottom center, `height` is axis vector, `r1`/`r2`/`r3` are vectors to vertex pairs.

### alea_halfspace

```c
alea_node_id_t alea_halfspace(const alea_system_t* sys, int surface_index, int sense);
```

Get the sense node from a surface index. `sense`: -1 for negative (inside), +1 for positive (outside). Returns `ALEA_NODE_ID_INVALID` if out of range.

---

## CSG Construction - Operations

### alea_union

```c
alea_node_id_t alea_union(alea_system_t* sys, alea_node_id_t a, alea_node_id_t b);
```

Boolean union of two nodes (A | B).

### alea_intersection

```c
alea_node_id_t alea_intersection(alea_system_t* sys, alea_node_id_t a, alea_node_id_t b);
```

Boolean intersection of two nodes (A & B).

### alea_difference

```c
alea_node_id_t alea_difference(alea_system_t* sys, alea_node_id_t a, alea_node_id_t b);
```

Boolean difference (A - B). Equivalent to `alea_intersection(sys, a, alea_complement(sys, b))`.

### alea_complement

```c
alea_node_id_t alea_complement(alea_system_t* sys, alea_node_id_t a);
```

Boolean complement (NOT A).

### alea_union_n

```c
alea_node_id_t alea_union_n(alea_system_t* sys, const alea_node_id_t* nodes, size_t count);
```

N-ary union. Builds a balanced tree from `count` nodes.

### alea_intersection_n

```c
alea_node_id_t alea_intersection_n(alea_system_t* sys, const alea_node_id_t* nodes, size_t count);
```

N-ary intersection. Builds a balanced tree from `count` nodes.

---

## CSG Construction - Cells

### alea_add_cell

```c
int alea_add_cell(alea_system_t* sys, int cell_id, alea_node_id_t root,
                  int material_index, double density, int universe);
```

Add a cell to the system. `cell_id` is the MCNP-style cell number (0 for auto-assign). `root` is the CSG tree root. `material_index` is the material number (0 for void). Returns cell index or -1 on error.

### alea_set_fill

```c
int alea_set_fill(alea_system_t* sys, int cell_index, int fill_universe, int transform);
```

Set a cell's fill universe and optional transform index.

### alea_cell_set_comment

```c
int alea_cell_set_comment(alea_system_t* sys, int cell_index, const char* comment);
```

Set comment lines (preceding "C" comments in MCNP) for a cell.

### alea_cell_set_inline_comment

```c
int alea_cell_set_inline_comment(alea_system_t* sys, int cell_index, const char* comment);
```

Set inline comment ("$" comment in MCNP) for a cell.

### alea_cell_set_material

```c
int alea_cell_set_material(alea_system_t* sys, int cell_index, int material_index);
```

Set cell material. Pass `ALEA_MATERIAL_VOID` (-1) for void. Derives the MCNP material ID from the material index automatically.

### alea_cell_set_density

```c
int alea_cell_set_density(alea_system_t* sys, int cell_index, double density);
```

Set cell density. Signed convention: negative = g/cm3, positive = atoms/b-cm, 0 = void.

### alea_cell_set_universe

```c
int alea_cell_set_universe(alea_system_t* sys, int cell_index, int universe_id);
```

Set which universe the cell belongs to. Invalidates the universe index.

### alea_cell_remove

```c
int alea_cell_remove(alea_system_t* sys, int cell_index);
```

Remove a cell by index. Frees per-cell allocations, compacts the array, rebuilds the cell ID hashmap, and triggers the `on_cell_removed` hook. Invalidates the universe index.

---

## Universe Operations

### alea_flatten

```c
int alea_flatten(alea_system_t* sys, int universe_id);
```

Flatten a universe hierarchy: resolve fill references and merge cells into a single level.

### alea_simplify_and_prune_cells

```c
void alea_simplify_and_prune_cells(alea_system_t* sys,
                                   alea_simplify_stats_t* stats);
```

Full optimization pass on all cells: NNF conversion, associative Boolean normalization, balancing, contradiction detection, and empty-cell removal. This does not flatten universe hierarchy. `stats` can be NULL.

**`alea_simplify_stats_t` fields**: `nodes_before`, `nodes_after`, `complements_eliminated`, `double_negations`, `idempotent_reductions`, `absorption_reductions`, `subtrees_deduplicated`, `cell_complements_expanded`, `contradictions_found`, `tautologies_found`, `empty_cells_removed`, `union_branches_absorbed`, `union_common_factors`, `union_branches_subsumed`.

### alea_split_union_cells

```c
int alea_split_union_cells(alea_system_t* sys);
```

Split cells whose root is a union into multiple simpler cells. Each branch becomes a separate cell inheriting the original's material, density, and universe. Returns number of new cells created. Call `alea_build_universe_index()` after.

### alea_extract_universe

```c
alea_system_t* alea_extract_universe(const alea_system_t* sys, int universe_id);
```

Extract a universe into a standalone system.

### alea_merge

```c
int alea_merge(alea_system_t* target, const alea_system_t* source, int id_offset);
```

Merge cells/surfaces from `source` into `target`, offsetting IDs by `id_offset`.

---

## Information

### alea_cell_count

```c
size_t alea_cell_count(const alea_system_t* sys);
```

### alea_surface_count

```c
size_t alea_surface_count(const alea_system_t* sys);
```

### alea_universe_count

```c
size_t alea_universe_count(const alea_system_t* sys);
```

### alea_stats

```c
void alea_stats(const alea_system_t* sys, alea_stats_t* stats);
```

Fill a stats structure with system information.

### alea_print_summary

```c
void alea_print_summary(const alea_system_t* sys);
```

Print a one-line summary to stdout.

### alea_tree_print

```c
void alea_tree_print(const alea_system_t* sys, alea_node_id_t node_id);
```

Pretty-print a CSG tree to stdout.

---

## Cell Information

### alea_cell_get

```c
int alea_cell_get(const alea_system_t* sys, size_t index,
                  int* cell_id, int* material_id, double* density,
                  int* universe_id, int* fill_universe, alea_node_id_t* root);
```

Get cell info by index. Any output pointer can be NULL.

### alea_cell_find

```c
int alea_cell_find(const alea_system_t* sys, int cell_id);
```

Find cell index by MCNP cell ID. Returns -1 if not found.

### alea_cell_get_info

```c
int alea_cell_get_info(const alea_system_t* sys, size_t index, alea_cell_info_t* info);
```

Get comprehensive cell info (struct-based). See `alea_cell_info_t` in `alea_types.h` for all fields including temperature, lattice data, and comments.

### alea_cell_find_info

```c
int alea_cell_find_info(const alea_system_t* sys, int cell_id, alea_cell_info_t* info);
```

Find cell by MCNP ID and get info.

### alea_cell_expr

```c
char* alea_cell_expr(const alea_system_t* sys, size_t cell_index,
                     const char* union_op, const char* inter_op,
                     const char* compl_op);
```

Get CSG expression string for a cell. Caller must `free()` the returned string.

```c
// MCNP-style
char* expr = alea_cell_expr(sys, i, ":", " ", "#");
// OpenMC-style
char* expr = alea_cell_expr(sys, i, " | ", " ", "~");
```

### alea_get_cell_id

```c
int alea_get_cell_id(const alea_system_t* sys, int cell_index);
```

Get MCNP cell ID from cell index.

### alea_cells_in_universe

```c
int alea_cells_in_universe(const alea_system_t* sys, int universe_id,
                           int* out_indices, size_t max_count);
```

Get cell indices in a universe. Returns number of cells found.

---

## Surface Information

### alea_surface_get

```c
int alea_surface_get(const alea_system_t* sys, size_t index,
                     int* surface_id, alea_primitive_type_t* type,
                     alea_node_id_t* pos_node, alea_node_id_t* neg_node,
                     alea_boundary_type_t* boundary_type);
```

Get surface info by index. Any output pointer can be NULL.

### alea_surface_find

```c
int alea_surface_find(const alea_system_t* sys, int surface_id);
```

Find surface index by MCNP surface ID. Returns -1 if not found.

---

## Universe Information

### alea_universe_get

```c
int alea_universe_get(const alea_system_t* sys, size_t index,
                      int* universe_id, size_t* cell_count, alea_bbox_t* bbox);
```

Get universe info by index.

### alea_universe_find

```c
int alea_universe_find(const alea_system_t* sys, int universe_id);
```

Find universe index by ID. Returns -1 if not found.

---

## CSG Node Inspection

Walk the CSG tree starting from a cell's root node. Enables reconstruction of region descriptions.

```c
// Example: walk a cell's tree
alea_node_id_t root;
alea_cell_get(sys, idx, NULL, NULL, NULL, NULL, NULL, &root);
alea_operation_t op = alea_node_operation(sys, root);
if (op == ALEA_OP_PRIMITIVE) {
    int sid = alea_node_surface_id(sys, root);
    int sense = alea_node_sense(sys, root);
} else {
    alea_node_id_t left = alea_node_left(sys, root);
    alea_node_id_t right = alea_node_right(sys, root);
}
```

### alea_node_operation

```c
alea_operation_t alea_node_operation(const alea_system_t* sys, alea_node_id_t node);
```

Get operation type: `ALEA_OP_PRIMITIVE`, `ALEA_OP_UNION`, `ALEA_OP_INTERSECTION`, `ALEA_OP_DIFFERENCE`, or `ALEA_OP_COMPLEMENT`.

### alea_node_left / alea_node_right

```c
alea_node_id_t alea_node_left(const alea_system_t* sys, alea_node_id_t node);
alea_node_id_t alea_node_right(const alea_system_t* sys, alea_node_id_t node);
```

Get children of a boolean node. Returns `ALEA_NODE_ID_INVALID` for leaf nodes.

### alea_node_primitive_type

```c
alea_primitive_type_t alea_node_primitive_type(const alea_system_t* sys, alea_node_id_t node);
```

Get primitive type for a leaf node (e.g., `ALEA_PRIMITIVE_SPHERE`).

### alea_node_primitive_id

```c
alea_primitive_id_t alea_node_primitive_id(const alea_system_t* sys, alea_node_id_t node);
```

Get primitive ID. Two nodes with the same primitive_id reference the same geometric surface (possibly with different sense).

### alea_node_primitive_data

```c
int alea_node_primitive_data(const alea_system_t* sys, alea_node_id_t node,
                             alea_primitive_data_t* out);
```

Get primitive geometric data (the union type `alea_primitive_data_t` — access the field matching the primitive type).

### alea_node_sense

```c
int alea_node_sense(const alea_system_t* sys, alea_node_id_t node);
```

Get sense: +1 (positive/outside) or -1 (negative/inside), 0 if not a primitive.

### alea_node_surface_id

```c
int alea_node_surface_id(const alea_system_t* sys, alea_node_id_t node);
```

Get MCNP surface ID associated with a primitive node.

---

## Filtering and Extraction

### alea_extract_region

```c
alea_system_t* alea_extract_region(const alea_system_t* sys, const alea_bbox_t* bbox);
```

Create a new system containing cells whose bounding boxes intersect the given region.

### alea_get_cells_in_bbox

```c
size_t alea_get_cells_in_bbox(const alea_system_t* sys, const alea_bbox_t* bbox,
                              int* out_indices, size_t max_count);
```

Get cell indices whose bounding box intersects a region.

### alea_get_cells_by_material

```c
size_t alea_get_cells_by_material(const alea_system_t* sys, int material_id,
                                  int* out_indices, size_t max_count);
```

Get cell indices with the given material.

### alea_get_cells_by_universe

```c
size_t alea_get_cells_by_universe(const alea_system_t* sys, int universe_id,
                                  int* out_indices, size_t max_count);
```

Get cell indices in a universe.

### alea_get_cells_filling_universe

```c
size_t alea_get_cells_filling_universe(const alea_system_t* sys, int universe_id,
                                       int* out_indices, size_t max_count);
```

Get cell indices that fill a given universe.

---

## Renumbering

### alea_renumber_cells

```c
int alea_renumber_cells(alea_system_t* sys, int start_id);
```

Renumber all cells sequentially starting from `start_id`.

### alea_renumber_surfaces

```c
int alea_renumber_surfaces(alea_system_t* sys, int start_id);
```

Renumber all surfaces sequentially starting from `start_id`.

### alea_offset_cell_ids

```c
int alea_offset_cell_ids(alea_system_t* sys, int offset);
```

Add `offset` to all cell IDs.

### alea_offset_surface_ids

```c
int alea_offset_surface_ids(alea_system_t* sys, int offset);
```

Add `offset` to all surface IDs.

### alea_offset_material_ids

```c
int alea_offset_material_ids(alea_system_t* sys, int offset);
```

Add `offset` to all material IDs.

---

## Materials

### alea_add_material

```c
int alea_add_material(alea_system_t* sys, int material_id);
```

Register a material. Returns material index or -1 on error. Use 0 or negative for auto-assign.

### alea_find_material_by_id

```c
int alea_find_material_by_id(const alea_system_t* sys, int material_id);
```

Find material index by MCNP ID. Returns -1 if not found.

### alea_create_mixture

```c
int alea_create_mixture(alea_system_t* sys, const int* mat_ids,
                        const double* fractions, size_t count, int new_mat_id);
```

Create a new material as a weighted mixture. Returns the assigned material ID or -1 on error.

---

## Volume Estimation

### alea_compute_bounding_sphere

```c
int alea_compute_bounding_sphere(const alea_system_t* sys, double tol,
                                 double* cx, double* cy, double* cz,
                                 double* radius);
```

Compute a tight bounding sphere for the entire model using interval arithmetic.

### alea_estimate_cell_volumes

```c
int alea_estimate_cell_volumes(const alea_system_t* sys,
                               double ox, double oy, double oz,
                               double radius, int n_rays,
                               double* volumes, double* rel_errors);
```

Estimate cell volumes using Cauchy-Crofton random ray tracing. `volumes` and `rel_errors` must be arrays of size `alea_cell_count(sys)`. `rel_errors` can be NULL.

### alea_volume_path_count

```c
size_t alea_volume_path_count(alea_system_t* sys);
```

Return the number of concrete hierarchical volume paths. The count is stable
until geometry mutation invalidates query caches.

### alea_volume_paths_get

```c
size_t alea_volume_paths_get(alea_system_t* sys,
                             alea_volume_path_t* out_paths,
                             size_t max_paths);
```

Enumerate concrete hierarchical volume paths. Each path identifies one terminal
cell under a concrete fill/lattice ancestry and has a dense `path_id` for volume
arrays. If the return value is greater than `max_paths`, the output was
truncated.

### alea_volume_path_at_point

```c
int alea_volume_path_at_point(alea_system_t* sys,
                              double x, double y, double z,
                              alea_volume_path_t* out_path);
```

Resolve a world-space point to a concrete hierarchical volume path. Returns `1`
when a path is found, `0` when no terminal path contains the point, and `-1` on
error.

### alea_estimate_path_volumes

```c
int alea_estimate_path_volumes(alea_system_t* sys, int n_rays,
                               double* volumes, double* rel_errors);
```

Estimate volumes per concrete path. Arrays must be sized to
`alea_volume_path_count(sys)`. This provides repeated-cell separation using the
hierarchical spatial index.

### alea_remove_cells_by_volume

```c
int alea_remove_cells_by_volume(alea_system_t* sys,
                                const double* volumes, double threshold);
```

Remove cells with estimated volume at or below `threshold`. Returns number of cells removed.

---

## Bounding Box Tightening

### alea_tighten_cell_bbox

```c
int alea_tighten_cell_bbox(const alea_system_t* sys, size_t cell_index,
                           double tol, alea_bbox_t* out);
```

Tighten a single cell's bounding box via interval arithmetic.

### alea_tighten_all_bboxes

```c
int alea_tighten_all_bboxes(alea_system_t* sys, double tol);
```

Tighten all cell bboxes. Returns number of cells tightened.

### alea_tighten_cell_bbox_numerical

```c
int alea_tighten_cell_bbox_numerical(alea_system_t* sys, int cell_index);
```

Numerically tighten a cell's bbox using LP vertex enumeration (plane-only cells) or octree recursive bisection (general cells). Useful for cells with infinite analytical bboxes.

---

## Macrobody Expansion

### alea_is_macrobody

```c
bool alea_is_macrobody(alea_primitive_type_t type);
```

Check if a primitive type is a macrobody (RCC, BOX, SPH, TRC, ELL, REC, WED, RHP, ARB).

### alea_expand_macrobody

```c
alea_node_id_t alea_expand_macrobody(alea_system_t* sys, alea_node_id_t node_id);
```

Expand a single macrobody node into primitive surfaces.

### alea_expand_all_macrobodies

```c
alea_node_id_t alea_expand_all_macrobodies(alea_system_t* sys, alea_node_id_t root_id);
```

Expand all macrobodies in a subtree.

### alea_expand_macrobodies_in_system

```c
int alea_expand_macrobodies_in_system(alea_system_t* sys);
```

Expand all macrobodies in all cells.

---

## Void Generation

Void options are read from `sys->config` (see [Configuration](#configuration)).

### alea_void_generate

```c
void_result_t* alea_void_generate(alea_system_t* sys, const alea_bbox_t* bounds);
```

Generate void regions via octree subdivision. Returns opaque result (must be freed with `alea_void_free`).

### alea_void_add_cells

```c
int alea_void_add_cells(alea_system_t* sys, void_result_t* result);
```

Add generated void regions as cells to the system.

### alea_void_add_graveyard

```c
int alea_void_add_graveyard(alea_system_t* sys, void_result_t* result);
```

Add a graveyard cell enclosing all void regions.

### alea_void_count

```c
size_t alea_void_count(const void_result_t* result);
```

Get number of void regions.

### alea_void_get

```c
int alea_void_get(const void_result_t* result, size_t index, alea_bbox_t* box);
```

Get the bounding box of a void region.

### alea_void_to_node

```c
alea_node_id_t alea_void_to_node(alea_system_t* sys, const void_result_t* result);
```

Convert void result to a CSG node (union of boxes).

### alea_void_merge

```c
int alea_void_merge(alea_system_t* sys, void_result_t* result);
```

Merge adjacent void cells to reduce count while balancing complexity.

### alea_void_free

```c
void alea_void_free(void_result_t* result);
```

Free a void result.

---

## Validation

### alea_validate

```c
int alea_validate(const alea_system_t* sys);
```

Validate system integrity. Returns 0 if valid, otherwise the number of issues found.

---

## Logging

### alea_log_set_level

```c
void alea_log_set_level(alea_log_level_t level);
```

Set the log level. Values: `ALEA_LOG_LEVEL_NONE` (0), `ALEA_LOG_LEVEL_ERROR` (1), `ALEA_LOG_LEVEL_WARN` (2, default), `ALEA_LOG_LEVEL_INFO` (3), `ALEA_LOG_LEVEL_DEBUG` (4), `ALEA_LOG_LEVEL_TRACE` (5).

### alea_log_get_level

```c
alea_log_level_t alea_log_get_level(void);
```

Get the current log level.

---

## Raycast (alea_raycast.h)

Ray tracing through CSG geometry.

### alea_raycast_result_create

```c
alea_raycast_result_t* alea_raycast_result_create(void);
```

Create a heap-allocated raycast result. Destroy with `alea_raycast_result_destroy()`.

### alea_raycast

```c
int alea_raycast(const alea_system_t* sys,
                 double ox, double oy, double oz,
                 double dx, double dy, double dz,
                 double t_max,
                 alea_raycast_result_t* result);
```

Cast a ray and find all cell intersections. Direction is normalized internally. `t_max=0` for infinite.

### alea_raycast_cell_aware

```c
int alea_raycast_cell_aware(const alea_system_t* sys,
                            double ox, double oy, double oz,
                            double dx, double dy, double dz,
                            double t_max,
                            alea_raycast_result_t* result);
```

Cell-aware raycast using per-cell surface index. It is semantically equivalent to
`alea_raycast`. For non-lattice models it may track through cells one at a time,
testing only surfaces belonging to each cell. For lattice models it uses the
canonical lattice-aware path so DDA element-boundary hits are included.

### alea_ray_first_cell

```c
int alea_ray_first_cell(const alea_system_t* sys,
                        double ox, double oy, double oz,
                        double dx, double dy, double dz,
                        double t_max, double* out_t);
```

Find first cell along ray. Returns cell ID or -1 if none. `out_t` receives distance (can be NULL).

### alea_raycast_segment_count

```c
size_t alea_raycast_segment_count(const alea_raycast_result_t* result);
```

Get number of segments in a raycast result.

### alea_raycast_segment_get

```c
int alea_raycast_segment_get(const alea_raycast_result_t* result, size_t index,
                             double* t_enter, double* t_exit,
                             int* cell_id, int* material_id, double* density,
                             int* enter_surface_id, int* exit_surface_id);
```

Get data for a segment. Any output pointer can be NULL. Surface IDs use `-1` for no boundary, `0` for synthetic lattice/DDA boundaries, and positive IDs for physical model surfaces.

### alea_raycast_path_length

```c
double alea_raycast_path_length(const alea_raycast_result_t* result, int material_id);
```

Total path length through a material. Pass -1 to sum all materials.

### alea_raycast_result_free

```c
void alea_raycast_result_free(alea_raycast_result_t* result);
```

Free internal buffers of a stack-allocated result.

### alea_raycast_result_destroy

```c
void alea_raycast_result_destroy(alea_raycast_result_t* result);
```

Destroy a heap-allocated result (from `alea_raycast_result_create`).

---

## Slice / Visualization (alea_slice.h)

2D slicing for geometry visualization: grid-based cell queries and analytical curve extraction.

### Slice View Setup

#### alea_slice_view_axis

```c
void alea_slice_view_axis(alea_slice_view_t* view,
                          int axis, double value,
                          double u_min, double u_max,
                          double v_min, double v_max);
```

Initialize a slice view for an axis-aligned slice. `axis`: 0=X (YZ plane), 1=Y (XZ plane), 2=Z (XY plane).

#### alea_slice_view_init

```c
void alea_slice_view_init(alea_slice_view_t* view,
                          double ox, double oy, double oz,
                          double nx, double ny, double nz,
                          double ux, double uy, double uz,
                          double u_min, double u_max,
                          double v_min, double v_max);
```

Initialize a slice view with arbitrary orientation. Normal vector is normalized; up vector is orthogonalized.

### Grid Queries

#### alea_find_cells_grid

```c
int alea_find_cells_grid(const alea_system_t* sys,
                         const alea_slice_view_t* view,
                         int nu, int nv,
                         int universe_depth,
                         int* out_cell_ids,
                         int* out_material_ids,
                         uint8_t* out_errors);
```

Find cells on a 2D grid. `universe_depth`: -1=innermost, 0=root, N=level N. `out_material_ids` and `out_errors` can be NULL. Error codes: 0=ok, 1=overlap, 2=undefined.

#### alea_check_grid_overlaps

```c
int alea_check_grid_overlaps(const alea_system_t* sys,
                             const alea_slice_view_t* view,
                             int nu, int nv,
                             int universe_depth,
                             const int* cell_ids,
                             uint8_t* errors);
```

Comprehensive overlap detection. Re-queries every non-void pixel to detect fully-nested overlaps.

#### alea_check_grid_overlaps_curves

```c
int alea_check_grid_overlaps_curves(const alea_system_t* sys,
                                    const alea_slice_view_t* view,
                                    const alea_slice_curves_t* curves,
                                    int nu, int nv,
                                    int universe_depth,
                                    const int* cell_ids,
                                    uint8_t* errors);
```

Efficient nested overlap detection. Only probes pixels where a surface curve crosses without a cell-ID transition. Returns number of new overlap pixels found.

### Analytical Curves

#### alea_get_slice_curves

```c
alea_slice_curves_t* alea_get_slice_curves(const alea_system_t* sys,
                                           const alea_slice_view_t* view);
```

Get exact surface boundaries as parametric curves. Returns opaque result (free with `alea_slice_curves_free`).

#### alea_slice_curves_count

```c
size_t alea_slice_curves_count(const alea_slice_curves_t* curves);
```

#### alea_slice_curves_get

```c
int alea_slice_curves_get(const alea_slice_curves_t* curves, size_t index, alea_curve_t* out);
```

Get curve data. Curve types include `ALEA_CURVE_LINE`, `ALEA_CURVE_CIRCLE`, `ALEA_CURVE_ELLIPSE`, `ALEA_CURVE_POLYGON`, `ALEA_CURVE_QUARTIC`, and more (see `alea_curve_type_t`).

#### alea_slice_curves_bounds

```c
void alea_slice_curves_bounds(const alea_slice_curves_t* curves,
                              double* u_min, double* u_max,
                              double* v_min, double* v_max);
```

Get bounding box of all curves.

#### alea_slice_curves_free

```c
void alea_slice_curves_free(alea_slice_curves_t* curves);
```

### Label Positions

#### alea_find_label_positions

```c
int alea_find_label_positions(const int* ids, int width, int height,
                              int min_pixels,
                              alea_label_position_t** out_labels,
                              int* out_count);
```

Find optimal label positions for cell/material regions. Handles non-convex regions. Caller frees `*out_labels` with `free()`.

#### alea_find_surface_label_positions

```c
int alea_find_surface_label_positions(const alea_slice_curves_t* curves,
                                      double x_min, double x_max,
                                      double y_min, double y_max,
                                      int width, int height, int margin,
                                      alea_label_position_t** out_labels,
                                      int* out_count);
```

Find label positions for surface curves. Caller frees `*out_labels` with `free()`.

```c
int alea_find_surface_label_positions_on_boundaries(
    const alea_slice_curves_t* curves, const int* boundary_ids,
    double x_min, double x_max, double y_min, double y_max,
    int width, int height, int margin,
    alea_label_position_t** out_labels, int* out_count);
```

Boundary-aware variant. Pass the same `cell_ids` or `material_ids` grid used for contour drawing so labels are only placed on drawn contours.

### Error Checking

#### alea_check_slice_errors

```c
alea_slice_error_result_t* alea_check_slice_errors(
    const alea_system_t* sys,
    const alea_slice_view_t* view,
    const alea_slice_curves_t* curves,
    int universe_depth);
```

Check surface curves for geometry errors (overlaps/gaps) using CSG point-containment queries. Free with `alea_slice_errors_free()`.

#### alea_check_slice_errors_grid

```c
alea_slice_error_result_t* alea_check_slice_errors_grid(
    const alea_system_t* sys,
    const alea_slice_view_t* view,
    const alea_slice_curves_t* curves,
    const int* cell_ids,
    const uint8_t* grid_errors,
    int nu, int nv);
```

Fast error checking using a pre-computed cell grid. Uses O(1) pixel lookups with CSG fallback for ambiguous cases. Free with `alea_slice_errors_free()`.

#### alea_slice_errors_free

```c
void alea_slice_errors_free(alea_slice_error_result_t* result);
```

### Debug

#### alea_slice_curve_set_debug

```c
void alea_slice_curve_set_debug(int enable);
```

Enable/disable debug output for slice curve generation.

#### alea_slice_point_trace_set_debug

```c
void alea_slice_point_trace_set_debug(int enable);
```

Enable/disable debug tracing for cell lookup in slices.

---

## 3D Rendering (alea_render.h)

Publication-quality 3D rendering with Phong shading, cutaway views, and shadow rays. Pure CPU, OpenMP parallelized.

### Configuration

#### render_config_init

```c
void render_config_init(render_config_t* cfg);
```

Initialize with defaults (1920x1080, 45 FOV, material coloring, no clipping).

#### render_config_free

```c
void render_config_free(render_config_t* cfg);
```

Free config resources (custom color table).

Key `render_config_t` fields:

| Field | Default | Description |
|-------|---------|-------------|
| `width`, `height` | 1920, 1080 | Image size |
| `fov` | 45.0 | Vertical FOV in degrees (0 = orthographic) |
| `color_mode` | `RENDER_COLOR_MATERIAL` | Color by material/cell/universe/density |
| `render_mode` | `RENDER_MODE_SOLID` | Solid/X-ray/depth/cellid/matid |
| `shadows` | 0 | Enable shadow rays |
| `edges` | 0 | Enable edge darkening |
| `aa_samples` | 1 | NxN supersampling (1=off) |
| `clips[N]` | - | Up to 16 clipping planes |
| `universe_depth` | -1 | Universe hierarchy level |

### Camera

#### render_camera_setup

```c
int render_camera_setup(render_camera_t* cam,
                        const render_config_t* cfg,
                        const alea_system_t* sys);
```

Set up camera from config. Auto-fits from model bounding sphere if eye/target not set.

#### render_camera_ray

```c
void render_camera_ray(const render_camera_t* cam,
                       int width, int height,
                       double px, double py,
                       double* ox, double* oy, double* oz,
                       double* dx, double* dy, double* dz);
```

Generate camera ray for pixel `(px, py)`.

### Rendering

#### render_scene

```c
int render_scene(const alea_system_t* sys,
                 const render_config_t* cfg,
                 const render_camera_t* cam,
                 render_framebuffer_t* fb);
```

Main render function. The system must be loaded and indexed.

#### render_framebuffer_create

```c
render_framebuffer_t* render_framebuffer_create(int width, int height, int aux);
```

Create a framebuffer. `aux=1` to allocate auxiliary buffers (depth, cell_id, material_id, normal).

#### render_framebuffer_free

```c
void render_framebuffer_free(render_framebuffer_t* fb);
```

### Post-Processing

#### render_edge_darken

```c
void render_edge_darken(render_framebuffer_t* fb);
```

Darken edges where cell ID changes.

#### render_tonemap

```c
void render_tonemap(const render_framebuffer_t* fb, uint8_t* pixels);
```

Convert float framebuffer to 8-bit RGB. `pixels` must be `width*height*3` bytes.

#### render_get_color

```c
void render_get_color(int id, render_color_mode_t mode,
                      const render_config_t* cfg,
                      float* r, float* g, float* b);
```

Get palette color for a material/cell/universe ID.

### Image Output

```c
int render_write_png(const char* filename, const uint8_t* pixels, int width, int height);
int render_write_bmp(const char* filename, const uint8_t* pixels, int width, int height);
int render_write_ppm(const char* filename, const uint8_t* pixels, int width, int height);
int render_write_image(const char* filename, const uint8_t* pixels, int width, int height);
```

Write images. `render_write_image` auto-detects format from extension (.png, .bmp, .ppm).

#### render_write_aux

```c
int render_write_aux(const char* base_filename, const render_framebuffer_t* fb);
```

Write auxiliary maps (depth, cellid, matid, normal).

#### render_load_colors

```c
int render_load_colors(const char* filename, render_config_t* cfg);
```

Load custom color palette from file.

---

## Mesh Export (alea_mesh.h)

Sample CSG geometry onto a structured rectilinear grid and export to Gmsh or VTK format. This is grid sampling, not boundary-conforming mesh generation.

### alea_mesh_config_init

```c
void alea_mesh_config_init(alea_mesh_config_t* cfg);
```

Initialize with defaults: nx=ny=nz=10, format=GMSH, auto bounds with 1% padding.

Key `alea_mesh_config_t` fields:

| Field | Default | Description |
|-------|---------|-------------|
| `x_min`..`z_max` | 0 | Bounds (0,0 = auto-detect from model) |
| `nx`, `ny`, `nz` | 10 | Elements per axis |
| `x_nodes`, `y_nodes`, `z_nodes` | NULL | Custom node positions (NULL=uniform) |
| `format` | `ALEA_MESH_GMSH` | `ALEA_MESH_GMSH` or `ALEA_MESH_VTK` |
| `void_material_id` | 0 | Material ID for void regions |
| `auto_pad` | 0.01 | Fractional padding for auto-bounds |
| `sampling_mode` | `ALEA_MESH_SAMPLE_SUBCELL` | Center, corners, regular, stratified, or adaptive sampling |
| `subsamples_per_axis` | 2 | Initial per-axis sample count |
| `mixed_threshold` | 0 | Tolerance before a voxel is flagged as mixed |
| `target_error` | 0.05 | Adaptive total-variation convergence target |
| `max_refine_depth` | 3 | Maximum adaptive refinements |
| `max_samples_per_voxel` | 32768 | Hard cumulative query limit per voxel |
| `max_total_samples` | 0 | Whole-run query budget; zero is unlimited |
| `sampling_seed` | fixed | Reproducible stratified-sampling seed |
| `workers` | 1 | Sampling workers; 0 selects the OpenMP runtime default |
| `bounds_mode` | `ALEA_MESH_BOUNDS_LEGACY` | Compatibility, explicit, or root-AABB inference |
| `fields` | all current result fields | Arrays retained in the result |
| `progress` | NULL | Optional callback after each completed Z slab; nonzero cancels |

### alea_mesh_sample

```c
alea_mesh_result_t* alea_mesh_sample(const alea_system_t* sys,
                                     const alea_mesh_config_t* cfg);
```

Sample CSG geometry onto a structured grid. Auto-detects bounds if all zero. Returns result (free with `alea_mesh_result_free`), or NULL on error.

By default, each voxel is sampled with a 2x2x2 subcell lattice. `material_ids` stores the dominant sampled material for each voxel; `cell_ids` stores the most frequently observed cell among samples of that material. Exact ties select the lowest material ID and then the lowest cell ID, with `tie_flags` recording the ambiguity. Arrays use X-fastest, then Y, then Z ordering.

Material fractions are point-count estimates, not exact volume fractions. A voxel is mixed when `(1 - dominant_fraction) > mixed_threshold`. In adaptive mode, fractions come from the final refinement level while `sample_counts` records all point queries spent across levels.

`alea_mesh_result_t` also stores sparse per-voxel material fractions:

| Field | Description |
|-------|-------------|
| `mixed_flags` | `nx*ny*nz` flags, one for each voxel |
| `dominant_fractions` | Largest sampled material fraction per voxel |
| `sample_counts` | Number of point samples used for each voxel |
| `tie_flags` | `ALEA_MESH_TIE_MATERIAL` / `ALEA_MESH_TIE_CELL` diagnostics |
| `estimated_errors` | Empirical difference between the last two adaptive estimates |
| `refinement_flags` | Indicates that an adaptive work/depth limit was reached |
| `mixed_count` | Number of voxels flagged as mixed |
| `fraction_spans` | `nx*ny*nz` spans into the packed fraction array |
| `fractions` | Packed `(material_id, fraction)` entries |
| `fraction_count` | Number of packed entries |
| `fields` | Materialized `ALEA_MESH_FIELD_*` arrays |
| `bounds_source` | Whether coordinates were explicit, custom, or inferred |
| `bounds_padding` | Fractional automatic padding actually applied |
| `sampling_mode`, `sampling_seed`, `target_error` | Estimator provenance |

Automatic bounds now use the union of bounded cells placed directly in the
root universe. Unplaced universe definitions are ignored, each axis retains its
own extent, and padding is relative to that axis. Production workflows should
still prefer explicit bounds.

Set `cfg.fields` to retain only arrays required by the caller. Material
discovery remains complete even when sparse per-voxel fractions are omitted.

### alea_mesh_export

```c
int alea_mesh_export(const alea_mesh_result_t* mesh,
                     alea_mesh_format_t fmt, const char* filename);
```

Export a sampled mesh to file.

### alea_mesh_export_stream

```c
int alea_mesh_export_stream(const alea_mesh_result_t* mesh,
                            alea_mesh_format_t fmt, FILE* out);
```

Export a sampled mesh to an open stream.

Both convenience exporters write mixed flag, dominant sampled fraction, tie
flag, sample count, estimated error, and refinement status when those arrays
are available. Filename export uses a
sibling temporary file and replaces the destination only after all writes and
the close succeed.

For explicit output fields, initialize `alea_mesh_export_options_t` and use the
`_ex` forms:

```c
alea_mesh_export_options_t options;
alea_mesh_export_options_init(&options);
options.fields |= ALEA_MESH_EXPORT_MATERIAL_FRACTIONS;
options.max_fraction_materials = 32;

int alea_mesh_export_ex(const alea_mesh_result_t* mesh,
                        alea_mesh_format_t fmt, const char* filename,
                        const alea_mesh_export_options_t* options);
int alea_mesh_export_stream_ex(const alea_mesh_result_t* mesh,
                               alea_mesh_format_t fmt, FILE* out,
                               const alea_mesh_export_options_t* options);
```

Per-material arrays are opt-in because they make sparse in-memory composition
dense in the output. Gmsh keeps dominant-material physical groups and writes
diagnostics as `$ElementData`; VTK writes cell scalar arrays.

`alea_mesh_visit()` runs the same sampler but invokes a callback for each voxel
instead of retaining per-voxel result arrays. Fraction pointers passed to the
callback are valid only for that call.

When built with `USE_OPENMP=1`, `workers > 1` parallelizes fixed, non-adaptive
sampling when sparse fractions and callbacks are disabled. The implementation
uses per-worker scratch and a deterministic material-table merge. Requests that
need ordered callbacks, packed sparse fractions, adaptive whole-run budgets, or
progress callbacks use the serial path.

### Adaptive grids

`alea_adaptive_grid_sample()` returns a separate octree representation; it does
not alter `alea_mesh_result_t`. Start with `alea_adaptive_grid_config_init()`,
then configure its `sampling` member and the spatial limits:

```c
alea_adaptive_grid_config_t cfg;
alea_adaptive_grid_config_init(&cfg);
cfg.sampling.nx = cfg.sampling.ny = cfg.sampling.nz = 4;
cfg.max_grid_depth = 3;
cfg.max_cells = 250000;

alea_adaptive_grid_result_t *grid =
    alea_adaptive_grid_sample(sys, &cfg);
alea_adaptive_grid_export(grid, ALEA_MESH_VTK, "adaptive.vtk");
alea_adaptive_grid_result_free(grid);
```

Mixed or empirically unstable voxels are split into eight children. Cell IDs
are stable and 1-based; every node records its parent and child IDs. Depth and
cell-count limits are visible in each leaf's flags. The result is currently an
unbalanced, nonconforming octree (`balanced == 0`). Export writes leaf voxels as
unstructured hexahedra in Gmsh 2.2 or VTK legacy format. It improves spatial
resolution but is not boundary-conforming meshing.

### alea_mesh_export_system

```c
int alea_mesh_export_system(const alea_system_t* sys,
                            const alea_mesh_config_t* cfg,
                            const char* filename);
```

One-shot: sample and export in one call.

### alea_mesh_result_free

```c
void alea_mesh_result_free(alea_mesh_result_t* mesh);
```

---

## MCNP Module (alea_mcnp.h)

Link `libalea_mcnp.a` (or use `libalea_full.a` which bundles everything).

### Loading

#### mcnp_load

```c
mcnp_model_t* mcnp_load(const char* filename);
```

Load an MCNP input file. Returns model with `model->sys` ready for queries.

#### mcnp_load_string

```c
mcnp_model_t* mcnp_load_string(const char* input, size_t len);
```

Load MCNP from a string buffer. Pass `len=0` to use `strlen`.

### Export

#### mcnp_export

```c
int mcnp_export(const mcnp_model_t* model, const char* filename);
```

Export MCNP model to file.

#### mcnp_export_stream

```c
int mcnp_export_stream(const mcnp_model_t* model, FILE* out);
```

Export to an open stream.

#### mcnp_export_system

```c
int mcnp_export_system(const alea_system_t* sys, const char* filename);
```

Export a bare system (no MCNP params) with default config.

#### mcnp_export_system_stream

```c
int mcnp_export_system_stream(const alea_system_t* sys, FILE* out);
```

Export a bare system to a stream.

### Model Management

#### mcnp_model_destroy

```c
void mcnp_model_destroy(mcnp_model_t* model);
```

Destroy model, freeing sys (if owned) and params.

#### mcnp_model_wrap

```c
mcnp_model_t* mcnp_model_wrap(alea_system_t* sys);
```

Create a non-owning model wrapper. `mcnp_model_destroy()` will NOT free the system. Cell params are populated with defaults.

#### mcnp_cell_params / mcnp_cell_params_const

```c
mcnp_cell_params_t* mcnp_cell_params(mcnp_model_t* m, size_t idx);
const mcnp_cell_params_t* mcnp_cell_params_const(const mcnp_model_t* m, size_t idx);
```

Get cell params by index (bounds-checked). Returns NULL if out of range.

#### mcnp_model_reserve_params

```c
int mcnp_model_reserve_params(mcnp_model_t* model, size_t cap);
```

Ensure the params array has room for at least `cap` entries.

#### mcnp_model_add_params

```c
int mcnp_model_add_params(mcnp_model_t* model);
```

Append a new cell params entry with defaults (imp=1.0). Returns index or -1.

#### mcnp_model_register_hooks

```c
void mcnp_model_register_hooks(mcnp_model_t* model);
```

Register cell callbacks so that `alea_add_cell()`, split, merge, etc. automatically grow/copy the parallel params array.

### MCNP Export Configuration

The `mcnp_export_config_t` struct controls export behavior:

| Field | Default | Description |
|-------|---------|-------------|
| `surface_policy` | 0 | 0=emit macrobodies, 1=force primitives |
| `trcl_mode` | 0 | 0=preserve, 1=bake transforms |
| `transform_mode` | 0 | 0=original, 1=inline, 2=cards |
| `mcnp_max_col` | 80 | Maximum line width |
| `mcnp_cont_indent` | 5 | Continuation indent |

Access via `model->export_config`. Default available as `MCNP_EXPORT_CONFIG_DEFAULT`.

---

## OpenMC Module (alea_openmc.h)

Link `libalea_openmc.a` (or use `libalea_full.a`).

### Loading

#### openmc_load

```c
openmc_model_t* openmc_load(const char* filename);
```

Load an OpenMC XML file (geometry.xml or model.xml).

#### openmc_load_string

```c
openmc_model_t* openmc_load_string(const char* input, size_t len);
```

Load from a string buffer.

### Export

#### openmc_export

```c
int openmc_export(const openmc_model_t* model, const char* filename);
```

Export OpenMC model to file.

#### openmc_export_stream

```c
int openmc_export_stream(const openmc_model_t* model, FILE* out);
```

Export to stream.

#### openmc_export_system

```c
int openmc_export_system(const alea_system_t* sys, const char* filename);
```

Export a bare system with default config.

#### openmc_export_system_stream

```c
int openmc_export_system_stream(const alea_system_t* sys, FILE* out);
```

### Model Management

#### openmc_model_destroy

```c
void openmc_model_destroy(openmc_model_t* model);
```

Destroy model. Frees system if owned.

#### openmc_model_wrap

```c
openmc_model_t* openmc_model_wrap(alea_system_t* sys);
```

Create a non-owning wrapper around an existing system.

---

## Nuclear Data (`alea_nucdata.h`)

Nuclear data module for ACE-format cross sections, reaction classification, and data decoding. All functions are declared in `alea_nucdata.h` with types in `alea_nucdata_types.h`. Energies are in MeV, cross sections in barns (microscopic) or cm⁻¹ (macroscopic). All functions are thread-safe for read-only access to loaded nuclides.

All objects are user-owned — no hidden state or caches. The user loads an xsdir, loads nuclides from it, and frees them when done.

```c
alea_nuc_xsdir_t* xsdir = alea_nuc_xsdir_load("/path/to/xsdir");
alea_nuc_nuclide_t* u235 = alea_nuc_load_nuclide(xsdir, "92235.80c");
double sigma = alea_nuc_xs_total(u235, 1.0);
alea_nuc_nuclide_free(u235);
alea_nuc_xsdir_free(xsdir);
```

### Cross-Section Directory (xsdir)

#### alea_nuc_xsdir_load

```c
alea_nuc_xsdir_t* alea_nuc_xsdir_load(const char* path);
```

Load an xsdir or xsdata file. Returns user-owned xsdir, or NULL on error.

#### alea_nuc_xsdir_load_dir

```c
alea_nuc_xsdir_t* alea_nuc_xsdir_load_dir(const char* dirpath);
```

Load all `.xsd` files from a directory (FENDL-style per-nuclide xsdir). Returns user-owned xsdir, or NULL on error.

#### alea_nuc_xsdir_free

```c
void alea_nuc_xsdir_free(alea_nuc_xsdir_t* xsdir);
```

Free an xsdir and all its entries.

#### alea_nuc_xsdir_find

```c
const alea_nuc_xsdir_entry_t* alea_nuc_xsdir_find(const alea_nuc_xsdir_t* xsdir, const char* zaid);
```

Find an xsdir entry by ZAID string (e.g. `"92235.80c"`). Returns NULL if not found.

#### alea_nuc_xsdir_count

```c
size_t alea_nuc_xsdir_count(const alea_nuc_xsdir_t* xsdir);
```

Number of entries in the loaded xsdir.

### ACE Table I/O

#### alea_nuc_ace_read

```c
alea_error_t alea_nuc_ace_read(const char* path, int address, int file_type, alea_nuc_ace_table_t* table);
```

Read raw ACE table from file. `file_type` is 1 (ASCII) or 2 (binary). `address` is the start line or byte offset.

#### alea_nuc_ace_free

```c
void alea_nuc_ace_free(alea_nuc_ace_table_t* table);
```

Free ACE table internals. Does not free the struct itself.

### Nuclide Loading

#### alea_nuc_load_nuclide

```c
alea_nuc_nuclide_t* alea_nuc_load_nuclide(const alea_nuc_xsdir_t* xsdir, const char* zaid);
```

Load and decode a nuclide by ZAID. Caller owns the returned nuclide and must free it with `alea_nuc_nuclide_free()`.

#### alea_nuc_nuclide_free

```c
void alea_nuc_nuclide_free(alea_nuc_nuclide_t* nuc);
```

Free a nuclide and all its data.

### Microscopic Cross Sections

All take a nuclide pointer and energy in MeV, return barns.

| Function | Description |
|----------|-------------|
| `alea_nuc_xs_total(nuc, energy)` | Total cross section |
| `alea_nuc_xs_absorption(nuc, energy)` | Absorption cross section (capture + fission) |
| `alea_nuc_xs_elastic(nuc, energy)` | Elastic scattering cross section |
| `alea_nuc_xs_reaction(nuc, mt, energy)` | Cross section for a specific MT reaction number |
| `alea_nuc_xs_heating(nuc, energy)` | Heating number (MeV-barn). Works for both neutron and photon nuclides |
| `alea_nuc_heating_per_collision(nuc, energy)` | Average energy deposited per collision (MeV) = heating / sigma_total |

### Photon Cross Sections

Log-log interpolation on photoatomic data. Takes energy in MeV, returns barns.

| Function | Description |
|----------|-------------|
| `alea_nuc_photon_xs_incoherent(nuc, energy)` | Compton (incoherent) scattering |
| `alea_nuc_photon_xs_coherent(nuc, energy)` | Rayleigh (coherent) scattering |
| `alea_nuc_photon_xs_photoelectric(nuc, energy)` | Photoelectric absorption |
| `alea_nuc_photon_xs_pair(nuc, energy)` | Pair production (including triplet) |

### Energy Grid Utilities

#### alea_nuc_energy_lookup

```c
int alea_nuc_energy_lookup(const double* energy_grid, int n, double E, double* frac);
```

Binary search on an ascending energy grid. Returns index `i` where `grid[i] <= E < grid[i+1]` and interpolation fraction in `frac`. Returns -1 if out of range.

#### alea_nuc_interp_loglog

```c
double alea_nuc_interp_loglog(const double* grid, const double* values, int n, double x);
```

Log-log interpolation on a grid. Falls back to lin-lin for zero or negative values.

### URR Probability Tables

#### alea_nuc_urr_factors

```c
int alea_nuc_urr_factors(const alea_nuc_nuclide_t* nuc, double energy, double xi, double factors[5]);
```

Sample URR probability band using random number `xi`. Writes 5 multiplicative factors to `factors`: [total, elastic, fission, capture, heating]. Returns 1 if URR applies, 0 otherwise.

### Material Composition

#### alea_nuc_material_create

```c
alea_nuc_material_t* alea_nuc_material_create(void);
```

Create an empty material.

#### alea_nuc_material_destroy

```c
void alea_nuc_material_destroy(alea_nuc_material_t* mat);
```

Free the material (does not free its nuclides).

#### alea_nuc_material_add

```c
alea_error_t alea_nuc_material_add(alea_nuc_material_t* mat, alea_nuc_nuclide_t* nuclide, double number_density);
```

Add a nuclide with number density in atoms/barn-cm.

### Macroscopic Cross Sections

Σ_macro = Σᵢ Nᵢ · σᵢ(E). Returns cm⁻¹.

| Function | Description |
|----------|-------------|
| `alea_nuc_mat_xs_total(mat, energy)` | Macroscopic total cross section |
| `alea_nuc_mat_xs_absorption(mat, energy)` | Macroscopic absorption cross section |
| `alea_nuc_mat_xs_elastic(mat, energy)` | Macroscopic elastic cross section |
| `alea_nuc_mean_free_path(mat, energy)` | Mean free path (cm) = 1 / Σ_total |
| `alea_nuc_sample_distance(mat, energy, xi)` | Sample distance to next collision: −ln(1−ξ) / Σ_total |

### Nuclide and Reaction Selection

#### alea_nuc_sample_nuclide

```c
int alea_nuc_sample_nuclide(const alea_nuc_material_t* mat, double energy, double xi, alea_nuc_nuclide_t** out_nuclide);
```

Sample which nuclide is hit in a material. Returns component index.

#### alea_nuc_sample_reaction

```c
int alea_nuc_sample_reaction(const alea_nuc_nuclide_t* nuc, double energy, double xi, int* out_mt);
```

Sample which reaction (MT) occurs on a nuclide. Returns reaction index (0 = elastic).

### Reaction Classification

| Function | Description |
|----------|-------------|
| `alea_nuc_reaction_classify(mt)` | Classify by MT number. Returns `ALEA_NUC_RXN_ABSORPTION`, `ALEA_NUC_RXN_SCATTER`, or `ALEA_NUC_RXN_MULTIPLY` |
| `alea_nuc_reaction_yield(nuc, mt, energy)` | Average neutron yield for a reaction |

### Fission

#### alea_nuc_nu_bar

```c
double alea_nuc_nu_bar(const alea_nuc_nuclide_t* nuc, double energy);
```

Average number of neutrons per fission, ν̄(E). Supports polynomial and tabular representations.

### Multigroup

Collapse continuous-energy data into group constants.

#### alea_nuc_mg_create

```c
alea_nuc_multigroup_t* alea_nuc_mg_create(int n_groups, const double* bounds);
```

Create multigroup structure. `bounds` has `n_groups+1` entries in descending order (MeV).

#### alea_nuc_mg_destroy

```c
void alea_nuc_mg_destroy(alea_nuc_multigroup_t* mg);
```

Free multigroup data.

#### alea_nuc_mg_collapse

```c
alea_error_t alea_nuc_mg_collapse(alea_nuc_multigroup_t* mg, const alea_nuc_nuclide_t* nuc);
```

Collapse pointwise cross sections into group constants using 1/E flux weighting. Fills sigma_t, sigma_a, sigma_s, sigma_f, nu_sigma_f, chi, and the scattering matrix.

#### alea_nuc_mg_scatter / alea_nuc_mg_scatter_adjoint

```c
double alea_nuc_mg_scatter(const alea_nuc_multigroup_t* mg, int g_from, int g_to);
double alea_nuc_mg_scatter_adjoint(const alea_nuc_multigroup_t* mg, int g_from, int g_to);
```

Forward and adjoint (transposed) scattering transfer cross sections.

#### alea_nuc_mg_sample_scatter

```c
int alea_nuc_mg_sample_scatter(const alea_nuc_multigroup_t* mg, int g_from, double xi, int adjoint);
```

Sample outgoing group from scattering CDF. Set `adjoint=1` for transposed matrix.

### Doppler Broadening

#### alea_nuc_doppler_broaden

```c
alea_error_t alea_nuc_doppler_broaden(alea_nuc_nuclide_t* nuc, double kT_target);
```

Broaden cross sections in-place to temperature `kT_target` (MeV). Can only broaden to higher temperatures. Modifies total, absorption, elastic, heating, and per-reaction cross sections.

### Utility

#### alea_nuc_parse_zaid

```c
alea_error_t alea_nuc_parse_zaid(const char* zaid, int* Z, int* A, int* meta, alea_nuc_table_type_t* type);
```

Parse ZAID string (e.g. `"92235.80c"`) into atomic number, mass number, metastable state, and table type.

#### alea_error_string

```c
const char* alea_error_string(alea_error_t err);
```

Human-readable error message for an error code.

### Nuclear Data Error Codes

| Code | Meaning |
|------|---------|
| `ALEA_OK` | Success |
| `ALEA_ERR_NULL_ARG` | NULL pointer argument |
| `ALEA_ERR_FILE_NOT_FOUND` | File does not exist |
| `ALEA_ERR_FILE_READ` | I/O error reading file |
| `ALEA_ERR_PARSE_ERROR` | Parse error in file format |
| `ALEA_ERR_OUT_OF_MEMORY` | Memory allocation failure |
| `ALEA_ERR_NOT_FOUND` | ZAID not found in xsdir |
| `ALEA_ERR_INVALID_ARG` | Invalid data format |
| `ALEA_ERR_UNSUPPORTED` | Unsupported operation (e.g. broadening to lower T) |

## Ray-query convenience APIs

`alea_ray_first_visible_query()` and `alea_ray_boundary_event_query()` use
caller-owned opaque reusable results. Create a result once, reuse it for
successive rays, then destroy it. A failed query clears its result.

Both option structures accept an older prefix when `struct_size` is set to
that prefix's size. Zero `t_max` means an unbounded ray.

### First visible

`alea_ray_first_visible_query()` returns the first non-void interval in the
requested range. Request `ALEA_RAY_FIRST_VISIBLE_SURFACE_ID` and/or
`ALEA_RAY_FIRST_VISIBLE_SURFACE_NORMAL` when those fields are needed.
If `t_min` cuts through an existing interval, the result is a cross-section
and has no reportable surface ID or normal.

### Boundary events

`alea_ray_boundary_event_query()` returns ordered ownership boundaries.
Kinds are physical, synthetic lattice, or unresolved. By default coincident
physical surfaces use the lowest positive surface ID; set
`include_all_coincident_physical` for every participant. Event order and
ownership fields are deterministic. `max_events` and `max_output_bytes` reject
oversized results rather than truncating them.

### Complete coverage slices

`alea_ray_coverage_query()` and `alea_ray_coverage_slice_query()` answer the
diagnostic question that selected-owner tracing cannot: which concrete owner
occurrences claim each elementary ray interval. They use the global
breakpoint/complete-owner engine and publish input-order CSR data:

```text
row_offsets:   rows      -> intervals
owner_offsets: intervals -> concrete owner occurrences
```

Create one `alea_ray_coverage_slice_result_t`, reuse it for many queries, and
destroy it when finished. Every accessor returns a borrowed read-only array;
it remains valid until the next successful query using that result or until
the result is destroyed. A failed query leaves the previous publication
unchanged.

```c
const double origins[] = {-2, 0, 0, -2, 2, 0};
const double directions[] = {1, 0, 0, 1, 0, 0};
const uint8_t direction_tags[] = {0, 0};
const double row_coordinates[] = {0, 2};

alea_ray_coverage_slice_options_t options;
alea_ray_coverage_slice_options_init(&options);
options.t_max = 4.0;
options.flags = ALEA_RAY_COVERAGE_DOMAIN;
options.domain_t_min = 0.5;
options.domain_t_max = 3.5;

alea_ray_coverage_slice_result_t* coverage =
    alea_ray_coverage_slice_result_create();
if (alea_ray_coverage_slice_query(sys, origins, directions, 2,
                                  direction_tags, row_coordinates,
                                  &options, coverage) != 0) {
    /* coverage still contains its previous successful result, if any */
}

const size_t* rows = alea_ray_coverage_slice_row_offsets(coverage);
const double* t0 = alea_ray_coverage_slice_t_enter(coverage);
const double* t1 = alea_ray_coverage_slice_t_exit(coverage);
const uint8_t* kinds = alea_ray_coverage_slice_kinds(coverage);
const size_t* owners = alea_ray_coverage_slice_owner_offsets(coverage);
const int* owner_cells = alea_ray_coverage_slice_owner_cell_ids(coverage);

for (size_t row = 0; row < alea_ray_coverage_slice_row_count(coverage); row++)
    for (size_t interval = rows[row]; interval < rows[row + 1]; interval++)
        for (size_t owner = owners[interval]; owner < owners[interval + 1]; owner++)
            /* [t0[interval], t1[interval]] has kind kinds[interval],
               claimed by owner_cells[owner] */;

alea_ray_coverage_slice_result_destroy(coverage);
```

`ALEA_RAY_COVERAGE_UNIQUE` means one complete owner chain; `GAP`, `OVERLAP`,
`UNDEFINED_FILL`, and `UNRESOLVED` are geometry findings or indeterminate
coverage. `TRUNCATED` means the owner budget retained only a prefix, so its
owners must not be treated as complete; use `owner_count_lower_bounds` to see
the minimum known count. `ALLOWED_EXTERIOR` is reported only with both
`ALEA_RAY_COVERAGE_DOMAIN` and `ALEA_RAY_COVERAGE_REPORT_EXTERIOR`.

The owner arrays expose cell, material, universe, fill-universe, depth,
occurrence key, immediate parent occurrence key, and resolution flags. The
occurrence-key arrays distinguish concrete placements such as repeated lattice
elements; cell ID alone is not owner identity.

`max_rows`, `max_intervals`, `max_owners`, and `max_output_bytes` bound the
published CSR output; zero means unlimited. Set `max_refinement_depth` to
enable deterministic midpoint refinement. The optional signature,
displacement, density, and finding signals select adjacent rows to split;
`alea_ray_coverage_slice_refinement_status()` reports whether the published
sample completed or stopped at a depth, row, or spacing limit.
`alea_ray_coverage_query()` is the scalar one-row adapter over the same
contract.

## Ray-slice rasterization

`alea_rasterize_ray_slice_compact()` fills caller-owned, tightly packed
row-major raster arrays from a result created by
`alea_trace_ray_slice_compact()`. The compact result must match the exact view
and row count; generic ray batches are rejected because their intervals are ray
distances rather than slice-U coordinates. Validation failures leave requested
output buffers unchanged.

`alea_trace_ray_slice_raster()` is the fused alternative. It derives the
compact fields required by the requested raster fields, creates a temporary
compact result, rasterizes it, and releases it before returning. Its
`max_trace_output_bytes` limit applies only to that temporary compact result;
callers own and budget the raster buffers themselves.

Raster ownership uses horizontal pixel centers. For an interval
`[u_enter, u_exit)`, a row with width `nu`, and pixel pitch
`du = (u_max-u_min)/nu`, the filled range is `[ceil((u_enter-u_min)/du-0.5),
ceil((u_exit-u_min)/du-0.5))`, clamped to `[0, nu]`. Void defaults are cell
`-1`, material `0`, universe `-1`, fill universe `-1`, density `0.0`, and
resolution flags `0`. Density is always the resolved leaf density; projected
owner depth applies to cell, material, universe, and fill-universe fields.

```c
size_t pixels = nu * nv;
int32_t* cells = malloc(pixels * sizeof(*cells));
int32_t* materials = malloc(pixels * sizeof(*materials));
alea_slice_raster_t raster;
alea_slice_raster_init(&raster);
raster.nu = nu;
raster.nv = nv;
raster.fields = ALEA_SLICE_RASTER_CELL_ID |
                ALEA_SLICE_RASTER_MATERIAL_ID;
raster.cell_ids = cells;
raster.material_ids = materials;

alea_slice_raster_options_t options;
alea_slice_raster_options_init(&options);  /* leaf ownership, no trace limit */
if (alea_trace_ray_slice_raster(sys, &view, &options, &raster) != 0) {
    /* cells and materials are unchanged on validation/trace failure */
}
free(materials);
free(cells);
```

Requested output buffers must not overlap. The rasterizer may fill independent
rows in parallel for large rasters, but produces the same bytes at every OpenMP
thread count. See `examples/c/ray_slice_raster_bench.c` for compact-trace,
standalone-raster, and fused timings.

## Directional slice trace caches

`alea_slice_directional_trace_cache_create()` creates an opaque reusable cache
for a slice view and sampling dimensions. It contains canonical U+/U-/V+/V-
boundary events plus ownership traces, and is invalid after a geometry change.
Pass it to `alea_validate_ray_slice_compact_with_directional_cache()` to reuse
the U ownership traces and obtain endpoint provenance arrays through the
`alea_ray_slice_validation_u_*_surface_ids()` and provenance-flag accessors.
The cache must be destroyed before its system.
It is valid only for the exact system generation, slice basis and bounds, and
sampling dimensions used at construction. Cache-aware calls reject a mismatch
without modifying an already published validation result. The cache is
read-only after construction and may be shared by independent consumers that
only read it; do not mutate the system while it is in use.

`alea_slice_surface_boundary_map_create_with_directional_cache()` is the
public provenance-map consumer for the same cache. It retains the canonical
short-edge fallback when directional stream evidence is incomplete or
contradictory.

Lua exposes the same workflow as
`sys:directional_trace_cache(view, width, height)` and
`sys:validate_ray_slice(view, rows [, options [, cache]])`. The returned table
contains intervals, trace reuse/execution masks, and (when a cache is supplied)
canonical forward endpoint surface IDs.
Lua cache userdata defers `sys:destroy()` until the cache is collected, so Lua
callers cannot free the system before a live cache; C callers must still destroy
their cache before destroying the system.

The currently published cache owns leaf-level ownership traces
(`projected_depth = -1`) with occurrence keys. A validation request for another
projected depth still receives canonical provenance from a matching cache, but
executes fresh ownership traces for that depth. Lua options accept
`projected_depth`, `include_agreements`, `absolute_tolerance`, and
`relative_tolerance`; C callers use `alea_ray_slice_validation_options_t` and
its documented interval/path/output budgets.
