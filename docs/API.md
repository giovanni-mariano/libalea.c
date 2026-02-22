# API Reference

Every public function in Alea, grouped by what you're trying to do.

**Headers**: `alea.h` (main), `alea_raycast.h` (ray tracing), `alea_slice.h` (visualization), `alea_types.h` (types and constants).

**Conventions**:
- Functions returning pointers return `NULL` on error
- Functions returning `int` return 0 on success, negative on error
- After any error, `alea_error()` gives a human-readable message
- Node IDs, cell IDs, and surface IDs are all 0-based internally but MCNP-numbered externally

---

## System Lifecycle

Create, destroy, copy, or reset a system.

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

Deep-copy a system. The clone is fully independent — modifying one does not affect the other.

### alea_reset

```c
void alea_reset(alea_system_t* sys);
```

Clear all data but keep allocated memory. Faster than destroy + create if you're reusing a system.

---

## Loading Models

### alea_load_mcnp

```c
alea_system_t* alea_load_mcnp(const char* filename);
```

Parse an MCNP input file and build a system. Handles cell cards, surface cards, data cards (materials, transforms), `LIKE BUT`, cell complements (`#cell`), macrobodies, universe fills, and lattices.

Returns NULL on parse error. Call `alea_error()` for details including line number.

### alea_load_mcnp_string

```c
alea_system_t* alea_load_mcnp_string(const char* input, size_t length);
```

Same as `alea_load_mcnp` but reads from a string buffer instead of a file.

### alea_load_openmc

```c
alea_system_t* alea_load_openmc(const char* filename);
```

Parse an OpenMC XML geometry file. Expects a `geometry.xml` (or similar) containing `<surface>`, `<cell>`, and `<lattice>` elements.

---

## Indexing

Call these after loading or building geometry, before doing queries.

### alea_build_universe_index

```c
int alea_build_universe_index(alea_system_t* sys);
```

Group cells by universe, compute bounding boxes, resolve cell complements. **Required** before point queries, overlap detection, slicing, or ray tracing. Returns 0 on success.

### alea_build_spatial_index

```c
int alea_build_spatial_index(alea_system_t* sys);
```

Build a KD-tree over cell instances for fast point queries. Called automatically on first query if not called explicitly. Useful to control when the (potentially slow) build happens.

**Thread safety**: Call this before any concurrent ray tracing or point queries. The lazy build on first use is not thread-safe — concurrent calls may corrupt shared caches. A warning is logged if lazy building occurs.

---

## Point Queries

Ask what's at a specific point in space.

### alea_find_cell

```c
int alea_find_cell(const alea_system_t* sys, double x, double y, double z);
```

Find the innermost cell containing the point in the root universe. Returns the cell index (into `sys->cells`), or -1 if no cell claims the point (void or undefined). Use `alea_get_cell_id()` to convert to the MCNP cell ID, or use `alea_find_cell_at()` to get both cell ID and material in one call.

### alea_material_at

```c
int alea_material_at(const alea_system_t* sys, double x, double y, double z);
```

Find the material number at a point. Returns 0 for void cells.

### alea_find_cell_at

```c
int alea_find_cell_at(const alea_system_t* sys, double x, double y, double z,
                           int* out_cell_id, int* out_material);
```

Find both cell ID and material at a point in one call. Either output pointer can be NULL.

### alea_find_all_cells

```c
int alea_find_all_cells(const alea_system_t* sys, double x, double y, double z,
                            alea_cell_hit_t* hits, size_t max_hits);
```

Find every cell in the hierarchy that contains the point, from outermost to innermost. Returns the number of hits. Each hit includes:

- `cell_id`, `material_id`, `universe_id`
- `fill_universe` (>0 if this cell fills with another universe)
- `depth` (0 for root, increases with nesting)
- `local_x`, `local_y`, `local_z` (point in this universe's coordinate frame)

Use this for debugging fill/transform issues.

### alea_point_inside

```c
bool alea_point_inside(const alea_system_t* sys, alea_node_id_t node,
                           double x, double y, double z);
```

Test if a point is inside a specific CSG region. Use when you have a node ID (e.g., from `alea_cell_get`) and want to check containment directly.

---

## Finding Problems

### alea_find_overlaps

```c
int alea_find_overlaps(const alea_system_t* sys, int* pairs, size_t max_pairs);
```

Find overlapping cell pairs by statistical sampling. Returns the number of overlap pairs found. The `pairs` array is filled with alternating cell IDs: `[cell_a1, cell_b1, cell_a2, cell_b2, ...]`. Provide `max_pairs` as half the array size.

This is statistical, not exhaustive. It samples random points within cell bounding boxes and checks for multiple occupancy.

### alea_validate

```c
int alea_validate(const alea_system_t* sys);
```

Check system integrity: dangling node references, invalid primitive IDs, etc. Returns 0 if valid, or the number of issues found.

---

## Bounding & Volume Estimation

### alea_compute_bounding_sphere

```c
int alea_compute_bounding_sphere(const alea_system_t* sys, double tol,
                                 double* cx, double* cy, double* cz, double* radius);
```

Compute the bounding sphere for the entire model by sampling rays from multiple directions. The sphere center and radius are returned through the output pointers. Returns 0 on success.

### alea_estimate_cell_volumes

```c
int alea_estimate_cell_volumes(const alea_system_t* sys, int n_rays,
                               double cx, double cy, double cz, double radius,
                               alea_cell_volume_t* volumes);
```

Estimate cell volumes using Monte Carlo ray tracing. Casts `n_rays` random rays through a sphere of the given center and radius. The `volumes` array must be pre-allocated with `alea_cell_count(sys)` entries. Each entry contains `volume` and `rel_error` (relative error estimate).

### alea_tighten_cell_bbox

```c
int alea_tighten_cell_bbox(const alea_system_t* sys, int cell_index,
                           double tol, alea_bbox_t* bbox);
```

Compute a tighter bounding box for a single cell using ray sampling. The `tol` parameter controls the precision.

### alea_tighten_all_bboxes

```c
int alea_tighten_all_bboxes(alea_system_t* sys, double tol);
```

Tighten bounding boxes for all cells in the system. Updates the bounding boxes in place.

### alea_tighten_cell_bbox_numerical

```c
int alea_tighten_cell_bbox_numerical(alea_system_t* sys, int cell_index);
```

Tighten a single cell's bounding box using numerical optimization (gradient descent on the CSG region).

---

## Building Geometry

### Creating Surfaces

Create surfaces with these functions. Each returns a surface index (into `sys->surfaces`), or -1 on error. Pass `surface_id=0` for automatic ID assignment.

Access the halfspace nodes via:

```c
alea_surface_entry_t* s = alea_surface_at(sys, idx);
alea_node_id_t inside  = s->neg_node;   // inside the surface (-S in MCNP)
alea_node_id_t outside = s->pos_node;   // outside the surface (+S in MCNP)
```

Or use the helper function:

```c
alea_node_id_t alea_halfspace(const alea_system_t* sys, int surface_index, int sense);
```

Where `sense = -1` returns the interior node and `sense = +1` returns the exterior node.

### Primitive Surfaces

```c
int alea_plane_surface(alea_system_t* sys, int surface_id,
                           double a, double b, double c, double d);
int alea_sphere_surface(alea_system_t* sys, int surface_id,
                            double cx, double cy, double cz, double r);
int alea_cylinder_z_surface(alea_system_t* sys, int surface_id,
                                double cx, double cy, double r);
int alea_cylinder_x_surface(alea_system_t* sys, int surface_id,
                                double cy, double cz, double r);
int alea_cylinder_y_surface(alea_system_t* sys, int surface_id,
                                double cx, double cz, double r);
int alea_box_surface(alea_system_t* sys, int surface_id,
                         double xmin, double xmax,
                         double ymin, double ymax,
                         double zmin, double zmax);
int alea_cone_z_surface(alea_system_t* sys, int surface_id,
                            double cx, double cy, double cz, double t_squared);
int alea_cone_x_surface(alea_system_t* sys, int surface_id,
                            double cx, double cy, double cz, double t_squared);
int alea_cone_y_surface(alea_system_t* sys, int surface_id,
                            double cx, double cy, double cz, double t_squared);
int alea_torus_z_surface(alea_system_t* sys, int surface_id,
                             double cx, double cy, double cz,
                             double major_radius, double minor_radius);
int alea_torus_x_surface(alea_system_t* sys, int surface_id,
                             double cx, double cy, double cz,
                             double major_radius, double minor_radius);
int alea_torus_y_surface(alea_system_t* sys, int surface_id,
                             double cx, double cy, double cz,
                             double major_radius, double minor_radius);
int alea_quadric_surface(alea_system_t* sys, int surface_id,
                             double A, double B, double C,
                             double D, double E, double F,
                             double G, double H, double I, double J);
```

All return the surface index (into `sys->surfaces`), or -1 on error. Access the entry with:

```c
alea_surface_entry_t* s = alea_surface_at(sys, idx);
alea_node_id_t inside  = s->neg_node;   // -S
alea_node_id_t outside = s->pos_node;   // +S
```

### Macrobody Surfaces

```c
int alea_rcc_surface(alea_system_t* sys, int surface_id,
                         double base_x, double base_y, double base_z,
                         double height_x, double height_y, double height_z,
                         double radius);
int alea_box_general_surface(alea_system_t* sys, int surface_id,
                                 double corner_x, double corner_y, double corner_z,
                                 double v1_x, double v1_y, double v1_z,
                                 double v2_x, double v2_y, double v2_z,
                                 double v3_x, double v3_y, double v3_z);
int alea_sph_surface(alea_system_t* sys, int surface_id,
                         double cx, double cy, double cz, double r);
int alea_trc_surface(alea_system_t* sys, int surface_id,
                         double base_x, double base_y, double base_z,
                         double height_x, double height_y, double height_z,
                         double base_radius, double top_radius);
int alea_ell_surface(alea_system_t* sys, int surface_id,
                         double v1_x, double v1_y, double v1_z,
                         double v2_x, double v2_y, double v2_z,
                         double major_axis_len);
int alea_rec_surface(alea_system_t* sys, int surface_id,
                         double base_x, double base_y, double base_z,
                         double height_x, double height_y, double height_z,
                         double axis1_x, double axis1_y, double axis1_z,
                         double axis2_x, double axis2_y, double axis2_z);
int alea_wed_surface(alea_system_t* sys, int surface_id,
                         double vertex_x, double vertex_y, double vertex_z,
                         double v1_x, double v1_y, double v1_z,
                         double v2_x, double v2_y, double v2_z,
                         double v3_x, double v3_y, double v3_z);
int alea_rhp_surface(alea_system_t* sys, int surface_id,
                         double base_x, double base_y, double base_z,
                         double height_x, double height_y, double height_z,
                         double r1_x, double r1_y, double r1_z,
                         double r2_x, double r2_y, double r2_z,
                         double r3_x, double r3_y, double r3_z);
```

### Boolean Operations

```c
alea_node_id_t alea_union(alea_system_t* sys, alea_node_id_t a, alea_node_id_t b);
alea_node_id_t alea_intersection(alea_system_t* sys, alea_node_id_t a, alea_node_id_t b);
alea_node_id_t alea_difference(alea_system_t* sys, alea_node_id_t a, alea_node_id_t b);
alea_node_id_t alea_complement(alea_system_t* sys, alea_node_id_t a);
```

Binary boolean operations. `difference(a, b)` is `a AND NOT b`.

```c
alea_node_id_t alea_union_n(alea_system_t* sys, const alea_node_id_t* nodes, size_t count);
alea_node_id_t alea_intersection_n(alea_system_t* sys, const alea_node_id_t* nodes, size_t count);
```

N-ary versions. Builds a balanced tree internally for better evaluation performance.

### Cells

```c
int alea_add_cell(alea_system_t* sys, int cell_id, alea_node_id_t root,
                      int material, double density, int universe);
```

Add a cell to the system. `cell_id` is the MCNP cell number. `root` is the CSG tree defining its region. `material` is the material number (0 for void). `universe` is the universe this cell belongs to (0 for default). Returns the cell index, or -1 on error.

```c
int alea_set_fill(alea_system_t* sys, int cell_index, int fill_universe, int transform);
```

Set a universe fill on a cell. `cell_index` is the return value from `alea_add_cell`. `fill_universe` is the universe ID to fill with. `transform` is a transform ID (0 for none).

### Surface Lookup

```c
alea_node_id_t alea_halfspace(const alea_system_t* sys, int surface_index, int sense);
```

Get the halfspace node for a surface at the given index. `sense = -1` returns the interior node, `sense = +1` returns the exterior node. Returns `ALEA_NODE_ID_INVALID` if the index is out of range.

To look up by surface ID rather than index:

```c
int alea_surface_find(const alea_system_t* sys, int surface_id);
```

Returns the surface index, or -1 if not found. Then use `alea_halfspace()` or `alea_surface_at()` to access the nodes.

---

## Visualization (alea_slice.h)

### Slice View

All slice operations use an `alea_slice_view_t` that bundles the plane definition with viewport bounds:

```c
typedef struct {
    alea_slice_plane_t plane;   /* Slice plane (orientation + position) */
    double u_min, u_max;            /* Horizontal viewport bounds */
    double v_min, v_max;            /* Vertical viewport bounds */
} alea_slice_view_t;
```

### View Setup

```c
void alea_slice_view_axis(alea_slice_view_t* view,
                               int axis, double value,
                               double u_min, double u_max,
                               double v_min, double v_max);
```

Initialize a view for axis-aligned slicing. `axis`: 0=X (YZ plane), 1=Y (XZ plane), 2=Z (XY plane). `value` is the coordinate along that axis. Bounds define the viewport in the slice plane.

```c
void alea_slice_view_init(alea_slice_view_t* view,
                               double ox, double oy, double oz,
                               double nx, double ny, double nz,
                               double ux, double uy, double uz,
                               double u_min, double u_max,
                               double v_min, double v_max);
```

Initialize a view with arbitrary orientation. `(ox, oy, oz)` is a point on the plane. `(nx, ny, nz)` is the normal (will be normalized). `(ux, uy, uz)` is an "up" hint for the horizontal direction (will be orthogonalized).

### Grid Queries

Sample cell/material IDs on a 2D pixel grid.

```c
int alea_find_cells_grid(const alea_system_t* sys,
                              const alea_slice_view_t* view,
                              int nu, int nv,
                              int universe_depth,
                              int* out_cell_ids, int* out_material_ids,
                              uint8_t* out_errors);
```

`nu` and `nv` are pixel counts. `universe_depth`: -1 = innermost cell, 0 = root level only, N = cell at depth N. Output arrays must be `nu * nv` elements. `out_material_ids` and `out_errors` can be NULL.

Error values: 0 = OK, 1 = overlap, 2 = undefined region.

When `out_errors` is provided, the grid query automatically rechecks boundary pixels (where adjacent pixels have different cell IDs) for overlaps. This catches overlapping geometry at cell transitions with minimal cost.

### Comprehensive Overlap Check

```c
int alea_check_grid_overlaps(const alea_system_t* sys,
                                  const alea_slice_view_t* view,
                                  int nu, int nv,
                                  int universe_depth,
                                  const int* cell_ids,
                                  uint8_t* errors);
```

Re-query every non-void pixel with a full hierarchy search to detect all overlaps, including fully-nested ones (e.g., concentric spheres where one cell is entirely inside another). Call after `alea_find_cells_grid()`. Pixels already flagged in `errors` are skipped. Returns 0 on success.

This is O(area) and significantly slower than the boundary recheck — use it only when comprehensive geometry validation is needed.

### Analytical Curves

Extract exact surface boundaries as parametric curves.

```c
alea_slice_curves_t* alea_get_slice_curves(const alea_system_t* sys,
                                                    const alea_slice_view_t* view);
```

Returns an opaque curves object. Must be freed with `alea_slice_curves_free`.

```c
size_t alea_slice_curves_count(const alea_slice_curves_t* curves);
```

Number of curves in the result.

```c
int alea_slice_curves_get(const alea_slice_curves_t* curves,
                              size_t index, alea_curve_t* out);
```

Get curve data at index. The `alea_curve_t` struct contains:

- `type`: one of `ALEA_CURVE_LINE`, `ALEA_CURVE_CIRCLE`, `ALEA_CURVE_ELLIPSE`, `ALEA_CURVE_LINE_SEGMENT`, `ALEA_CURVE_ARC`, `ALEA_CURVE_ELLIPSE_ARC`, `ALEA_CURVE_POLYGON`, `ALEA_CURVE_PARALLEL_LINES`
- `surface_id`: which surface generated this curve
- `data`: union with geometry-specific fields (point + direction for lines, center + radius for circles, etc.)
- `t_min`, `t_max`: parameter range

```c
void alea_slice_curves_bounds(const alea_slice_curves_t* curves,
                                  double* u_min, double* u_max,
                                  double* v_min, double* v_max);
```

Get bounding box of all curves.

```c
void alea_slice_curves_free(alea_slice_curves_t* curves);
```

Free curves result.

### Label Positioning

```c
int alea_find_label_positions(const int* ids, int width, int height,
                                  int min_pixels,
                                  alea_label_position_t** out_labels,
                                  int* out_count);
```

Find optimal label positions for regions in a grid. `ids` is a cell ID or material ID grid (from grid queries). `min_pixels` filters out small regions. Output: array of `alea_label_position_t` (caller frees with `free()`). Each label has `id`, `px`, `py` (pixel coordinates), and `pixel_count`.

```c
int alea_find_surface_label_positions(const alea_slice_curves_t* curves,
                                          double x_min, double x_max,
                                          double y_min, double y_max,
                                          int width, int height, int margin,
                                          alea_label_position_t** out_labels,
                                          int* out_count);
```

Find label positions for surfaces along curve boundaries. `margin` excludes labels too close to image edges.

### Debug

```c
void alea_slice_curve_set_debug(int enable);
void alea_slice_point_trace_set_debug(int enable);
```

Enable verbose debug output for curve generation or point lookups. Prints to stdout.

---

## Ray Tracing (alea_raycast.h)

### alea_raycast

```c
alea_raycast_result_t* alea_raycast_result_create(void);

int alea_raycast(const alea_system_t* sys,
                     double ox, double oy, double oz,
                     double dx, double dy, double dz,
                     double t_max, alea_raycast_result_t* result);
```

Cast a ray from `(ox, oy, oz)` in direction `(dx, dy, dz)` up to distance `t_max` (0 = infinite). The direction is normalized internally. Results are stored in the provided result object. The result object is safe to reuse across calls — prior allocations are freed internally.

### alea_raycast_cell_aware

```c
int alea_raycast_cell_aware(const alea_system_t* sys,
                                double ox, double oy, double oz,
                                double dx, double dy, double dz,
                                double t_max, alea_raycast_result_t* result);
```

Cell-aware ray tracing. Faster for large models — tests only surfaces belonging to the current cell at each step.

### alea_ray_first_cell

```c
int alea_ray_first_cell(const alea_system_t* sys,
                            double ox, double oy, double oz,
                            double dx, double dy, double dz,
                            double t_max, double* out_t);
```

Find the first cell hit by a ray. Returns the cell ID, or -1 if none. `out_t` receives the distance (can be NULL).

### Result Queries

```c
size_t alea_raycast_segment_count(const alea_raycast_result_t* result);

int alea_raycast_segment_get(const alea_raycast_result_t* result, size_t index,
                                 double* t_enter, double* t_exit,
                                 int* cell_id, int* material_id, double* density);
```

Get segment data. Each segment is a contiguous region where the ray is inside one cell. `cell_id = -1` means void.

```c
double alea_raycast_path_length(const alea_raycast_result_t* result, int material_id);
```

Total path length through a given material. Pass `material_id = -1` for total path through all materials.

```c
void alea_raycast_result_free(alea_raycast_result_t* result);
```

Free result memory.

---

## Export

### alea_export_mcnp

```c
int alea_export_mcnp(const alea_system_t* sys, const char* filename);
int alea_export_mcnp_stream(const alea_system_t* sys, FILE* out);
```

Write the geometry to an MCNP input file. Export options come from `sys->config`:

- `surface_policy`: `ALEA_EMIT_MACROBODY` or `ALEA_EMIT_SURFACES`
- `export_materials`: include `Mn` cards
- `export_transforms`: include `TRn` cards
- `universe_depth`: -1 = all, 0 = root only
- `fill_depth`: 0 = don't expand fills

### alea_export_openmc

```c
int alea_export_openmc(const alea_system_t* sys, const char* filename);
```

Write the geometry to an OpenMC XML file.

### alea_export_stream

```c
int alea_export_stream(const alea_system_t* sys, alea_export_format_t format, FILE* out);
```

Write to an arbitrary stream. `format`: `ALEA_EXPORT_FORMAT_MCNP` or `ALEA_EXPORT_FORMAT_OPENMC`.

---

## Universe Operations

### alea_flatten

```c
int alea_flatten(alea_system_t* sys, int universe_id);
```

Expand all fills in a universe, materializing every instance with transforms applied. After flattening, all cells are in the specified universe with no fills. Modifies the system in place.

### alea_extract_universe

```c
alea_system_t* alea_extract_universe(const alea_system_t* sys, int universe_id);
```

Create a new system containing only the cells from one universe. The original system is not modified.

### alea_merge

```c
int alea_merge(alea_system_t* target, const alea_system_t* source, int id_offset);
```

Copy all cells, surfaces, materials, and transforms from `source` into `target`, adding `id_offset` to all IDs to avoid collisions.

---

## Void Generation

### alea_void_generate

```c
void_result_t* alea_void_generate(alea_system_t* sys, const alea_bbox_t* bounds);
```

Find void regions within the bounding box using octree subdivision. Returns an opaque result. Void generation parameters come from `sys->config` (`void_max_depth`, `void_min_size`, `void_samples`).

### alea_void_add_cells

```c
int alea_void_add_cells(alea_system_t* sys, void_result_t* result);
```

Add the detected void regions as void cells (material 0) to the system. Returns the number of cells added.

### alea_void_count / alea_void_get

```c
size_t alea_void_count(const void_result_t* result);
int alea_void_get(const void_result_t* result, size_t index, alea_bbox_t* box);
```

Query individual void regions before adding them.

### alea_void_to_node

```c
alea_node_id_t alea_void_to_node(alea_system_t* sys, const void_result_t* result);
```

Convert the void result to a single CSG node (union of boxes). Useful for creating a custom void cell.

### alea_void_merge

```c
int alea_void_merge(alea_system_t* sys, void_result_t* result);
```

Reduce void cell count by merging adjacent regions. Configured by `merge_cell_weight`, `merge_surface_weight`, `merge_max_surfaces`, `merge_min_cells` in config.

### alea_void_free

```c
void alea_void_free(void_result_t* result);
```

Free void result memory.

---

## Model Information

### Counts

```c
size_t alea_cell_count(const alea_system_t* sys);
size_t alea_surface_count(const alea_system_t* sys);
size_t alea_universe_count(const alea_system_t* sys);
```

### Summary

```c
void alea_print_summary(const alea_system_t* sys);
```

Print cell/surface/material counts and other summary info to stdout.

### Statistics

```c
void alea_stats(const alea_system_t* sys, alea_stats_t* stats);
```

Fill a stats struct with dedup hits, memory usage, conversion counts, etc.

---

## Cell Information

### alea_cell_get

```c
int alea_cell_get(const alea_system_t* sys, size_t index,
                      int* cell_id, int* material_id, double* density,
                      int* universe_id, int* fill_universe, alea_node_id_t* root);
```

Get cell properties by index (0 to `alea_cell_count - 1`). Any output pointer can be NULL.

### alea_cell_find

```c
int alea_cell_find(const alea_system_t* sys, int cell_id);
```

Find a cell's index by MCNP cell ID. Returns -1 if not found. O(1) via hash table.

### alea_cell_get_info / alea_cell_find_info

```c
int alea_cell_get_info(const alea_system_t* sys, size_t index, alea_cell_info_t* info);
int alea_cell_find_info(const alea_system_t* sys, int cell_id, alea_cell_info_t* info);
```

Get a `alea_cell_info_t` struct with all cell properties including lattice parameters, fill info, and bounding box.

### alea_cells_in_universe

```c
int alea_cells_in_universe(const alea_system_t* sys, int universe_id,
                                int* out_indices, size_t max_count);
```

Get indices of all cells in a universe.

### alea_get_cell_id

```c
int alea_get_cell_id(const alea_system_t* sys, int cell_index);
```

Get the MCNP cell ID from a cell index.

---

## Surface Information

### alea_surface_get

```c
int alea_surface_get(const alea_system_t* sys, size_t index,
                          int* surface_id, alea_primitive_type_t* type,
                          alea_node_id_t* pos_node, alea_node_id_t* neg_node,
                          alea_boundary_type_t* boundary_type);
```

Get surface properties by index. Any output pointer can be NULL.

### alea_surface_find

```c
int alea_surface_find(const alea_system_t* sys, int surface_id);
```

Find a surface's index by MCNP surface ID. Returns -1 if not found.

---

## Universe Information

```c
int alea_universe_get(const alea_system_t* sys, size_t index,
                           int* universe_id, size_t* cell_count, alea_bbox_t* bbox);
int alea_universe_find(const alea_system_t* sys, int universe_id);
```

---

## Renumbering

```c
int alea_renumber_cells(alea_system_t* sys, int start_id);
int alea_renumber_surfaces(alea_system_t* sys, int start_id);
```

Reassign cell or surface IDs starting from `start_id`, preserving order.

```c
int alea_offset_cell_ids(alea_system_t* sys, int offset);
int alea_offset_surface_ids(alea_system_t* sys, int offset);
int alea_offset_material_ids(alea_system_t* sys, int offset);
```

Add a fixed offset to all IDs. Useful when merging models.

---

## Filtering

```c
size_t alea_get_cells_by_material(const alea_system_t* sys, int material_id,
                                       int* out_indices, size_t max_count);
size_t alea_get_cells_by_universe(const alea_system_t* sys, int universe_id,
                                       int* out_indices, size_t max_count);
size_t alea_get_cells_filling_universe(const alea_system_t* sys, int universe_id,
                                            int* out_indices, size_t max_count);
```

Find cells matching criteria. Returns the number found.

```c
alea_system_t* alea_extract_region(const alea_system_t* sys, const alea_bbox_t* bbox);
```

Create a new system containing only cells whose bounding boxes intersect the given region.

```c
size_t alea_get_cells_in_bbox(const alea_system_t* sys, const alea_bbox_t* bbox,
                                  int* out_indices, size_t max_count);
```

Get cell indices intersecting a bounding box.

```c
size_t alea_spatial_index_instance_count(const alea_system_t* sys);
```

Number of cell instances in the spatial index (reflects the expanded universe hierarchy).

---

## CSG Tree Inspection

Walk the CSG tree of a loaded or constructed model.

### alea_node_operation

```c
alea_operation_t alea_node_operation(const alea_system_t* sys, alea_node_id_t node);
```

Get the operation type: `ALEA_OP_PRIMITIVE`, `ALEA_OP_UNION`, `ALEA_OP_INTERSECTION`, `ALEA_OP_DIFFERENCE`, or `ALEA_OP_COMPLEMENT`.

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

Get primitive type for a leaf node. Returns 0 if not a primitive.

### alea_node_primitive_id

```c
alea_primitive_id_t alea_node_primitive_id(const alea_system_t* sys, alea_node_id_t node);
```

Get the deduplicated primitive ID. Two nodes with the same primitive ID reference the same geometric surface.

### alea_node_primitive_data

```c
int alea_node_primitive_data(const alea_system_t* sys, alea_node_id_t node,
                                 alea_primitive_data_t* out);
```

Get the raw primitive parameters (coefficients, center, radius, etc.).

### alea_node_sense

```c
int alea_node_sense(const alea_system_t* sys, alea_node_id_t node);
```

Get the sense of a primitive node: +1 (positive/outside) or -1 (negative/inside). Returns 0 if not a primitive.

### alea_node_surface_id

```c
int alea_node_surface_id(const alea_system_t* sys, alea_node_id_t node);
```

Get the MCNP surface ID associated with a primitive node. Returns 0 if not registered.

---

## Macrobody Expansion

```c
bool alea_is_macrobody(alea_primitive_type_t type);
```

Check if a primitive type is a macrobody (RCC, BOX, TRC, etc.).

```c
alea_node_id_t alea_expand_macrobody(alea_system_t* sys, alea_node_id_t node_id);
```

Expand a single macrobody node into its constituent primitives. Returns the new root node.

```c
alea_node_id_t alea_expand_all_macrobodies(alea_system_t* sys, alea_node_id_t root_id);
```

Expand all macrobodies in a tree.

```c
int alea_expand_macrobodies_in_system(alea_system_t* sys);
```

Expand all macrobodies in all cells.

---

## Materials

```c
int alea_create_mixture(alea_system_t* sys, const int* mat_ids,
                            const double* fractions, size_t count, int new_mat_id);
```

Create a mixture of existing materials. `new_mat_id = 0` for auto-assignment. Returns the assigned material ID.

---

## Configuration

```c
alea_config_t alea_get_config(const alea_system_t* sys);
void alea_set_config(alea_system_t* sys, const alea_config_t* config);
```

Get or set the system configuration. See [Concepts: Configuration](CONCEPTS.md#configuration) for field descriptions.

The default configuration is available as `ALEA_CONFIG_DEFAULT`.

---

## Logging

```c
void alea_log_set_level(int level);
int alea_log_get_level(void);
```

Set/get log level. Levels: `ALEA_LOG_NONE` (0), `ALEA_LOG_ERROR` (1), `ALEA_LOG_WARN` (2), `ALEA_LOG_INFO` (3), `ALEA_LOG_DEBUG` (4), `ALEA_LOG_TRACE` (5).

```c
void alea_log_set_callback(alea_log_callback_t callback, void* user_data);
```

Set a custom log handler:

```c
typedef void (*alea_log_callback_t)(alea_log_level_t level, const char* file,
                                        int line, const char* message, void* user_data);
```

---

## Error Handling

```c
const char* alea_error(void);
int alea_error_code(void);
void alea_error_clear(void);
```

Get the last error message, error code, or clear the error state.

---

## Interrupt Support

```c
void alea_interrupt(void);
void alea_clear_interrupt(void);
bool alea_interrupted(void);
```

Cooperative interruption for long-running operations. Signal-safe — can be called from a SIGINT handler. After the interrupted operation returns, call `alea_clear_interrupt()` before starting a new operation.

---

## Version

```c
const char* alea_version(void);
```

Returns version string (e.g., "0.2.0").

```c
#define ALEA_VERSION_MAJOR 0
#define ALEA_VERSION_MINOR 2
#define ALEA_VERSION_PATCH 0
```

Compile-time version macros.

---

## Types Quick Reference

### ID Types

| Type | Underlying | Meaning |
|------|-----------|---------|
| `alea_node_id_t` | `uint32_t` | Index into node array |
| `alea_primitive_id_t` | `uint32_t` | Index into primitive array |
| `alea_material_id_t` | `uint32_t` | Material number |
| `alea_surface_id_t` | `uint32_t` | Surface number |

### Special Values

| Constant | Value | Meaning |
|----------|-------|---------|
| `ALEA_NODE_ID_INVALID` | `UINT32_MAX` | Invalid/missing node |
| `ALEA_PRIMITIVE_ID_INVALID` | `UINT32_MAX` | Invalid/missing primitive |
| `ALEA_MATERIAL_NONE` | `0` | Void (no material) |

### Error Codes

| Code | Name | Meaning |
|------|------|---------|
| 0 | `ALEA_OK` | Success |
| 1 | `ALEA_ERR_NULL_ARG` | NULL pointer argument |
| 2 | `ALEA_ERR_INVALID_ID` | ID out of range |
| 3 | `ALEA_ERR_INVALID_ARG` | Invalid argument value |
| 4 | `ALEA_ERR_INVALID_STATE` | Invalid operation for current state |
| 5 | `ALEA_ERR_OUT_OF_MEMORY` | Allocation failure |
| 6 | `ALEA_ERR_FILE_NOT_FOUND` | File not found |
| 7 | `ALEA_ERR_FILE_READ` | File read error |
| 8 | `ALEA_ERR_FILE_WRITE` | File write error |
| 9 | `ALEA_ERR_PARSE_ERROR` | Parse/syntax error |
| 10 | `ALEA_ERR_UNSUPPORTED` | Unsupported feature |
| 13 | `ALEA_ERR_INTERRUPTED` | Interrupted by user |
| 14 | `ALEA_ERR_NOT_FOUND` | Item not found |
| 16 | `ALEA_ERR_OVERFLOW` | Buffer too small |

### Primitive Types

| Enum | MCNP |
|------|------|
| `ALEA_PRIMITIVE_PLANE` | P, PX, PY, PZ |
| `ALEA_PRIMITIVE_SPHERE` | S, SO, SX, SY, SZ |
| `ALEA_PRIMITIVE_CYLINDER_X` | CX, C/X |
| `ALEA_PRIMITIVE_CYLINDER_Y` | CY, C/Y |
| `ALEA_PRIMITIVE_CYLINDER_Z` | CZ, C/Z |
| `ALEA_PRIMITIVE_CONE_X` | KX, K/X |
| `ALEA_PRIMITIVE_CONE_Y` | KY, K/Y |
| `ALEA_PRIMITIVE_CONE_Z` | KZ, K/Z |
| `ALEA_PRIMITIVE_RPP` | RPP |
| `ALEA_PRIMITIVE_QUADRIC` | GQ, SQ |
| `ALEA_PRIMITIVE_TORUS_X` | TX |
| `ALEA_PRIMITIVE_TORUS_Y` | TY |
| `ALEA_PRIMITIVE_TORUS_Z` | TZ |
| `ALEA_PRIMITIVE_RCC` | RCC |
| `ALEA_PRIMITIVE_BOX` | BOX |
| `ALEA_PRIMITIVE_SPH` | SPH |
| `ALEA_PRIMITIVE_TRC` | TRC |
| `ALEA_PRIMITIVE_ELL` | ELL |
| `ALEA_PRIMITIVE_REC` | REC |
| `ALEA_PRIMITIVE_WED` | WED |
| `ALEA_PRIMITIVE_RHP` | RHP/HEX |
| `ALEA_PRIMITIVE_ARB` | ARB |

### Boundary Types

| Enum | MCNP syntax |
|------|-------------|
| `ALEA_BOUNDARY_TRANSMISSIVE` | (default) |
| `ALEA_BOUNDARY_REFLECTIVE` | `*S` |
| `ALEA_BOUNDARY_WHITE` | `+S` |
| `ALEA_BOUNDARY_PERIODIC` | paired surfaces |
| `ALEA_BOUNDARY_VACUUM` | (graveyard) |
