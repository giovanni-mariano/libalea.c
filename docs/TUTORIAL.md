# Tutorial

This tutorial walks through the main things you can do with Alea. Each section is self-contained: real code, real results, no hand-waving.

Build the library first, then compile examples from the repository root with:

```bash
make full
cc -std=c11 -Wall -Wextra -Iinclude example.c bin/libalea_full.a -lm -o example
```

`libalea_full.a` bundles the core library with the MCNP and OpenMC modules. Build it with `make full`. If you only need the core (no MCNP/OpenMC I/O), link against `libalea.a` instead.

## 1. Loading a Model

The most common starting point is an existing MCNP input file:

```c
#include <alea.h>
#include <alea_mcnp.h>
#include <stdio.h>

int main(void) {
    mcnp_model_t* model = mcnp_load("iter_blanket.inp");
    if (!model) {
        fprintf(stderr, "Load failed: %s\n", alea_error());
        return 1;
    }
    alea_system_t* sys = mcnp_model_system(model);  /* borrowed */

    // Build all shared query caches before concurrent or repeated queries.
    if (alea_prepare_query_acceleration(sys) != 0) {
        fprintf(stderr, "Query setup failed: %s\n", alea_error());
        mcnp_model_destroy(model);
        return 1;
    }

    alea_print_summary(sys);
    mcnp_model_destroy(model);
    return 0;
}
```

`mcnp_load` parses the cell cards, surface cards, data cards (materials, transforms), and builds the internal CSG tree. It handles `LIKE BUT`, cell complements (`#cell`), macrobodies, and universe fills.

For OpenMC (requires `alea_openmc.h`):

```c
openmc_model_t* model = openmc_load("geometry.xml");
alea_system_t* sys = openmc_model_system(model);  /* borrowed */
```

You can also load from a string instead of a file:

```c
const char* input = "1 1 -10.0 -1\n2 0 1\n\n1 SO 5.0\n\n";
mcnp_model_t* model = mcnp_load_string(input, strlen(input));
alea_system_t* sys = mcnp_model_system(model);
```

The model owns the borrowed system pointer. Destroying the model invalidates it.
Use `mcnp_model_take_system()` or `openmc_model_take_system()` when the system
must outlive the format wrapper; after detaching it, the caller must eventually
call `alea_destroy()`.

**Important**: prepare query caches after loading or after a geometry mutation.
`alea_prepare_query_acceleration()` is the general entry point for point,
raycast, slice, mesh, render, and validation workloads. The narrower
`alea_build_universe_index()` remains useful when only hierarchy lookup is
needed. Cache preparation is idempotent, and geometry-changing APIs invalidate
affected caches.

## 2. Asking Questions About the Geometry

Once loaded, the most useful thing is asking "what's at this point?"

### Which cell and material are at a point?

```c
int cell_id, material;
if (alea_find_cell_at(sys, 650.0, 0.0, 0.0, &cell_id, &material) == 0) {
    printf("Point is in cell %d\n", cell_id);
    printf("Material is %d\n", material);
} else {
    printf("Point is in void or undefined\n");
}
```

`alea_find_cell_at()` is the preferred single-point query API. It traverses the full universe hierarchy — if the point is in a cell that has `FILL=5`, it descends into universe 5, applies the inverse transform, and continues until it finds a terminal cell (one with a material or void, not another fill).

Returns -1 if no cell claims the point. This means either void or a geometry error.

### Convenience wrappers

```c
int cell_idx = alea_find_cell(sys, 650.0, 0.0, 0.0);  /* internal cell index */
int mat = alea_material_at(sys, 650.0, 0.0, 0.0);
```

`alea_find_cell()` and `alea_material_at()` are still available as convenience wrappers, but for new code prefer `alea_find_cell_at()`.

### Both at once

```c
int cell_id, material;
if (alea_find_cell_at(sys, 650.0, 0.0, 0.0, &cell_id, &material) == 0) {
    printf("cell=%d material=%d\n", cell_id, material);
}
```

### The full hierarchy

For debugging universe fills, you often want to see every cell the point passes through:

```c
alea_cell_hit_t hits[32];
int nhits = alea_find_all_cells(sys, 650.0, 0.0, 0.0, hits, 32);

for (int i = 0; i < nhits; i++) {
    printf("  depth %d: cell %d, universe %d, material %d",
           hits[i].depth, hits[i].cell_id,
           hits[i].universe_id, hits[i].material_id);
    if (hits[i].fill_universe > 0)
        printf(" -> FILL=%d", hits[i].fill_universe);
    printf("\n");
}
```

Each hit includes the local coordinates in that universe's frame (`local_x`, `local_y`, `local_z`), which is essential for debugging transform issues.

### Finding overlaps

```c
int pairs[200];  // pairs of cell IDs: [a1, b1, a2, b2, ...]
int noverlaps = alea_find_overlaps(sys, pairs, 100);

for (int i = 0; i < noverlaps; i++) {
    printf("Overlap: cell %d and cell %d\n", pairs[2*i], pairs[2*i+1]);
}
```

This is a cheap root-universe screen: it intersects each candidate pair's
bounding boxes and probes the eight corners plus the center. Returned values
are zero-based cell indices, not external cell IDs. It is deliberately bounded
and can miss overlaps; use it for a quick hint, not as evidence that a model is
clean.

### Transport-style validation

For actionable geometry diagnostics, use the occurrence-aware validator rather
than treating `alea_find_overlaps()` as a proof that a model is valid:

```c
#include <alea_geo_validator.h>

alea_geom_validator_options_t options;
alea_geom_validator_options_init(&options);
options.ray_count = 20000;
options.seed = 12345;
options.max_errors = 1000;

alea_geom_validator_result_t validation;
alea_geom_validator_result_init(&validation);

if (alea_validate_geometry(sys, &options, &validation) != 0) {
    fprintf(stderr, "validation failed: %s\n", alea_error());
} else {
    for (size_t i = 0; i < alea_geom_validator_error_count(&validation); i++) {
        alea_geom_error_t error;
        if (alea_geom_validator_error_get(&validation, i, &error) == 0) {
            printf("%s near cell %d at (%.6g, %.6g, %.6g)\n",
                   alea_geom_error_type_name(error.type),
                   error.found_cell_id,
                   error.crossing_point[0], error.crossing_point[1],
                   error.crossing_point[2]);
        }
    }
}

alea_geom_validator_result_free(&validation);
```

`alea_validate_geometry()` runs bounded randomized ray validation. For a known
ray use `alea_validate_geometry_ray()`. For a plotted plane, combine
`alea_get_slice_curves()` with `alea_validate_geometry_slice()`; analytical
curve sampling can expose fully nested overlaps that random rays miss. Set
`ALEA_GEOM_VALIDATE_DOMAIN_BOUNDS` plus `validation_bounds` when unowned space
inside a closed world box must be reported as an interior gap. Results may be
truncated when a configured error or work budget is reached, so inspect
`validation.truncated` and its counters.

## 3. Visualizing the Geometry

Alea provides two complementary approaches to 2D visualization:

1. **Grid queries**: sample cell/material IDs on a pixel grid (fast, gives you colors)
2. **Analytical curves**: extract exact surface boundaries (lines, circles, ellipses — gives you contours)

### Grid queries

```c
#include <alea_slice.h>

int width = 800, height = 800;
int* cell_ids = malloc(width * height * sizeof(int));
int* mat_ids  = malloc(width * height * sizeof(int));
uint8_t* errors = malloc(width * height);

// Set up a slice view: XY plane at z=0, from -100 to +100
alea_slice_view_t view;
alea_slice_view_axis(&view, 2, 0.0,
    -100.0, 100.0, -100.0, 100.0);

alea_find_cells_grid(sys, &view, width, height,
    -1,                       // universe_depth: -1 = innermost
    cell_ids, mat_ids, errors);
```

Each pixel gets a cell ID, a material ID, and an error flag:

- `errors[i] == 0`: normal, valid cell
- `errors[i] == 1`: overlap — multiple cells claim this point
- `errors[i] == 2`: undefined — no cell claims this point (geometry error)

The grid query automatically rechecks boundary pixels for overlaps (where adjacent cells differ). This catches overlapping geometry at cell transitions — like offset spheres — with minimal overhead. For comprehensive validation that also catches fully-nested overlaps (e.g., concentric spheres), call `alea_check_grid_overlaps()` afterward:

```c
alea_check_grid_overlaps(sys, &view, width, height, -1, cell_ids, errors);
```

This re-queries every non-void pixel and is O(area), so use it only when thorough validation is needed.

The `universe_depth` parameter controls which level of the hierarchy you see:

- `-1`: innermost cell (follow all fills) — what the transport code sees
- `0`: root-level cells only — useful for seeing the container structure
- `N`: cells at depth N

For other axes, change the axis parameter: 0=X (YZ plane), 1=Y (XZ plane), 2=Z (XY plane).

For an arbitrary plane:

```c
alea_slice_view_t view;
alea_slice_view_init(&view,
    0.0, 0.0, 0.0,     // origin
    0.0, 0.0, 1.0,     // normal (Z direction)
    1.0, 0.0, 0.0,     // up hint (X direction)
    -100.0, 100.0, -100.0, 100.0);  // viewport bounds

alea_find_cells_grid(sys, &view, width, height,
    -1, cell_ids, mat_ids, errors);
```

### Analytical curves

Grid queries give you pixel data. For sharp surface boundaries, extract the exact curves:

```c
alea_slice_view_t view;
alea_slice_view_axis(&view, 2, 0.0,
    -100.0, 100.0, -100.0, 100.0);
alea_slice_curves_t* curves = alea_get_slice_curves(sys, &view);

size_t ncurves = alea_slice_curves_count(curves);
for (size_t i = 0; i < ncurves; i++) {
    alea_curve_t c;
    alea_slice_curves_get(curves, i, &c);

    switch (c.type) {
    case ALEA_CURVE_LINE:
        // c.data.line.point, c.data.line.direction
        break;
    case ALEA_CURVE_CIRCLE:
        // c.data.circle.center, c.data.circle.radius
        break;
    case ALEA_CURVE_ELLIPSE:
        // c.data.ellipse.center, semi_a, semi_b, angle
        break;
    // ...
    }

    printf("Surface %d: %s\n", c.surface_id,
           c.type == ALEA_CURVE_LINE ? "line" :
           c.type == ALEA_CURVE_CIRCLE ? "circle" : "other");
}

alea_slice_curves_free(curves);
```

The typical workflow combines both: use the grid for pixel coloring, and overlay the curves for crisp surface boundaries. The `tools/mc_plotter.c` program does exactly this.

For new diagnostic tooling, the grid's three-state error byte is only a fast
overlay. Use `alea_find_cells_grid_coverage()` when you need explicit
none/unique/multiple coverage, or the APIs in `alea_geo_validator.h` when you
need occurrence-aware findings and provenance.

### Label positioning

To place cell or surface labels on a slice image:

```c
alea_label_position_t* labels;
int nlabels;
alea_find_label_positions(cell_ids, width, height, 100, &labels, &nlabels);

for (int i = 0; i < nlabels; i++) {
    printf("Cell %d: place label at pixel (%d, %d)\n",
           labels[i].id, labels[i].px, labels[i].py);
}
free(labels);
```

The `min_pixels` parameter (100 above) filters out tiny regions that are too small for a readable label. The algorithm finds a point guaranteed to be inside the region, close to its centroid — it handles non-convex shapes correctly.

### Sampling and exporting a 3D mesh

`alea_mesh_sample()` produces a rectilinear hexahedral grid. Initialize the
configuration first so newly added fields receive safe defaults:

```c
#include <alea_mesh.h>

alea_mesh_config_t mesh_options;
alea_mesh_config_init(&mesh_options);
mesh_options.nx = 80;
mesh_options.ny = 80;
mesh_options.nz = 80;
mesh_options.bounds_mode = ALEA_MESH_BOUNDS_AUTO;
mesh_options.sampling_mode = ALEA_MESH_SAMPLE_ADAPTIVE;
mesh_options.target_error = 0.02;
mesh_options.workers = 0;  /* parallel-backend default */

alea_mesh_result_t* mesh = alea_mesh_sample(sys, &mesh_options);
if (!mesh) {
    fprintf(stderr, "mesh sampling failed: %s\n", alea_error());
} else {
    alea_mesh_export(mesh, ALEA_MESH_VTK, "geometry.vtk");
    alea_mesh_result_free(mesh);
}
```

The default result field mask retains all current fields, including packed
per-material and per-concrete-cell sampled fractions. These are sampling
estimates, not exact volume fractions. Use `alea_mesh_export_ex()` to choose
diagnostic arrays written to VTK/Gmsh, `alea_mesh_visit()` to stream voxels
without retaining all arrays, or `alea_adaptive_grid_sample()` for a separate
nonconforming octree grid.

## 4. Tracing Rays

Ray tracing reports every cell the ray passes through, in order:

```c
#include <alea_raycast.h>

alea_raycast_result_t* result = alea_raycast_result_create();

// Ray from origin, going in +X direction, up to 500 cm
alea_raycast(sys, 0.0, 0.0, 0.0,   // origin
                      1.0, 0.0, 0.0,    // direction
                      500.0,             // max distance
                      result);

size_t nseg = alea_raycast_segment_count(result);
for (size_t i = 0; i < nseg; i++) {
    double t_enter, t_exit;
    int cell_id, material_id, enter_surface_id, exit_surface_id;
    double density;
    alea_raycast_segment_get(result, i,
        &t_enter, &t_exit, &cell_id, &material_id, &density,
        &enter_surface_id, &exit_surface_id);

    double thickness = t_exit - t_enter;
    printf("%8.2f - %8.2f cm: cell %d, mat %d, density %.3f, surfaces %d -> %d (%.2f cm)\n",
           t_enter, t_exit, cell_id, material_id, density,
           enter_surface_id, exit_surface_id, thickness);
}

alea_raycast_result_destroy(result);
```

### Finding the first cell along a ray

If you only need the first hit:

```c
double t;
int cell = alea_ray_first_cell(sys,
    0.0, 0.0, 0.0,    // origin
    1.0, 0.0, 0.0,    // direction
    500.0, &t);        // max distance, output distance

if (cell > 0)
    printf("First cell hit: %d at distance %.2f cm\n", cell, t);
```

### Path length through a material

```c
double steel_path = alea_raycast_path_length(result, 5);  // material 5
double total_path = alea_raycast_path_length(result, -1);  // all materials
```

### Cell-aware ray tracing

For large non-lattice models, cell-aware tracing can be faster. Instead of testing every surface globally, it tracks through cells one at a time:

```c
alea_raycast_cell_aware(sys, ox, oy, oz, dx, dy, dz, t_max, result);
```

Same interface and same result format. On lattice models this entry point uses the canonical DDA-aware path so lattice element-boundary hits are preserved.

For high-throughput work, use `alea_raycast_hier_batch()` and request only the
fields you consume. For visibility picking, `alea_ray_first_visible_query()`
avoids constructing a full segment list. For boundary inspection,
`alea_ray_boundary_event_query()` reports physical, synthetic lattice, and
unresolved events. The result objects for these queries are reusable: create
once, execute many queries, then destroy.

### Estimating volumes

Alea provides two different estimators:

- `alea_estimate_volumes_ex()` uses reproducible Cauchy–Crofton rays and
  reports volumes for concrete hierarchy paths, so repeated fill and lattice
  occurrences remain distinct.
- `alea_cell_estimate_volume()` uses a deterministic interval/octree method for
  one cell definition in its universe-local frame and returns lower/upper
  bounds plus convergence and resource-limit information.

```c
size_t path_count = alea_volume_path_count(sys);
double* volumes = calloc(path_count, sizeof(*volumes));
double* errors = calloc(path_count, sizeof(*errors));

alea_volume_estimate_options_t volume_options;
alea_volume_estimate_options_init(&volume_options);
volume_options.max_rays = 1000000;
volume_options.seed = 12345;
volume_options.requested_workers = 0;
volume_options.target_rel_error = 0.01;

alea_volume_estimate_stats_t volume_stats;
if (alea_estimate_volumes_ex(sys, &volume_options, volumes, errors,
                             &volume_stats) == 0) {
    printf("used %zu rays; converged=%d\n",
           volume_stats.rays_completed, volume_stats.converged);
}

free(errors);
free(volumes);
```

The output order is the order returned by `alea_volume_paths_get()`. The
production RNG is counter-based Philox, so a seed is reproducible independently
of worker count. Treat `converged == false` as a reported statistical result,
not necessarily an API failure.

## 5. Building Geometry from Scratch

You don't have to load from a file. You can build geometry programmatically:

```c
alea_system_t* sys = alea_create();

// Create surfaces and get their interior (neg_node) or exterior (pos_node)
int s1 = alea_sphere_surface(sys, 0, 0, 0, 0, 10.0);
int s2 = alea_cylinder_z_surface(sys, 0, 0, 0, 3.0);

alea_node_id_t outer = alea_halfspace(sys, s1, -1);  // inside sphere
alea_node_id_t hole  = alea_halfspace(sys, s2, -1);  // inside cylinder

// Boolean difference: sphere minus cylinder
alea_node_id_t region = alea_difference(sys, outer, hole);

// Add material 1, then cell 1 at 10.0 g/cm3 in universe 0.
int material_index = alea_add_material(sys, 1);
alea_add_cell(sys, 1, region, material_index, -10.0, 0);

alea_prepare_query_acceleration(sys);
```

### Creating surfaces

Each surface function returns a surface index. Convert it to a CSG halfspace
with `alea_halfspace()`; applications must not access internal surface storage.

```c
int idx = alea_plane_surface(sys, id, a, b, c, d);      // ax + by + cz + d = 0
int idx = alea_sphere_surface(sys, id, cx, cy, cz, r);
int idx = alea_cylinder_z_surface(sys, id, cx, cy, r);  // infinite along Z
int idx = alea_cylinder_x_surface(sys, id, cy, cz, r);  // infinite along X
int idx = alea_cylinder_y_surface(sys, id, cx, cz, r);  // infinite along Y
int idx = alea_box_surface(sys, id, xmin, xmax, ymin, ymax, zmin, zmax);
int idx = alea_cone_z_surface(sys, id, cx, cy, cz, t2); // t2 = tan^2(half-angle)

// Negative and positive sides of the surface equation.
alea_node_id_t inside  = alea_halfspace(sys, idx, -1);
alea_node_id_t outside = alea_halfspace(sys, idx, +1);
```

Pass `id=0` for automatic surface ID assignment.

### Boolean operations

```c
alea_node_id_t u = alea_union(sys, a, b);          // a OR b
alea_node_id_t i = alea_intersection(sys, a, b);   // a AND b
alea_node_id_t d = alea_difference(sys, a, b);     // a AND NOT b
alea_node_id_t c = alea_complement(sys, a);        // NOT a
```

For more than two operands:

```c
alea_node_id_t parts[] = {a, b, c, d};
alea_node_id_t all = alea_union_n(sys, parts, 4);
```

### Universe fills

To create nested geometry (like a fuel pin inside a lattice cell):

```c
// Universe 1: the fuel pin
int s_fuel = alea_sphere_surface(sys, 0, 0, 0, 0, 0.5);
alea_node_id_t fuel_r = alea_halfspace(sys, s_fuel, -1);
int fuel_material = alea_add_material(sys, 1);
int clad_material = alea_add_material(sys, 2);
alea_add_cell(sys, 10, fuel_r, fuel_material, -10.0, 1);  // universe 1

int s_clad = alea_sphere_surface(sys, 0, 0, 0, 0, 0.6);
alea_node_id_t clad_r = alea_difference(sys,
    alea_halfspace(sys, s_clad, -1),
    alea_halfspace(sys, s_fuel, -1));
alea_add_cell(sys, 11, clad_r, clad_material, -8.0, 1);  // universe 1

// Universe 0: container that fills with universe 1
int s_box = alea_box_surface(sys, 0, -5, 5, -5, 5, -5, 5);
alea_node_id_t box = alea_halfspace(sys, s_box, -1);
int cell_idx = alea_add_cell(sys, 1, box, ALEA_MATERIAL_VOID, 0.0, 0);
alea_set_fill(sys, cell_idx, 1, 0);  // fill with universe 1, no transform
```

Named transforms accept the normalized MCNP representation: three translation
values or twelve translation-and-rotation values. Angles are interpreted in
degrees when the final argument is nonzero.

```c
double translation[3] = {10.0, 0.0, 0.0};
if (alea_add_transform(sys, 7, translation, 3, 0) != 0)
    fprintf(stderr, "transform failed: %s\n", alea_error());
else
    alea_set_fill(sys, cell_idx, 1, 7);
```

Use `alea_add_inline_transform()` when importing or constructing an anonymous
cell-local transform. It deduplicates the normalized transform and returns the
assigned transform ID.

## 6. Exporting

### To MCNP (alea_mcnp.h)

```c
mcnp_export_system(sys, "output.inp");
```

Or with a model for full control over export settings:

```c
mcnp_export(model, "output.inp");
```

### To OpenMC (alea_openmc.h)

```c
openmc_export_system(sys, "geometry.xml");
```

### To Serpent (`alea_serpent.h`)

```c
#include <alea_serpent.h>

serpent_export_system(sys, "geometry.serpent");
```

### To a file stream

```c
FILE* f = fopen("output.inp", "w");
mcnp_export_system_stream(sys, f);
fclose(f);
```

### Export configuration

Export behavior is controlled by the system configuration:

```c
alea_config_t cfg = alea_get_config(sys);
cfg.surface_policy = ALEA_EMIT_SURFACES;   // decompose macrobodies to primitives
cfg.export_materials = true;           // include material cards
cfg.export_transforms = true;          // include TR cards
cfg.universe_depth = -1;               // export all universes
cfg.fill_depth = 0;                    // don't expand fills
alea_set_config(sys, &cfg);

mcnp_export_system(sys, "output.inp");
```

The `surface_policy` setting matters most:

- `ALEA_EMIT_MACROBODY` (default): keep RCC, BOX, etc. as macrobody surfaces in the output
- `ALEA_EMIT_SURFACES`: decompose macrobodies into their constituent planes, cylinders, etc.

## 7. Void Generation

Void generation finds regions within a bounding box that no cell covers. This is essential for creating the "graveyard" cell in MCNP:

```c
alea_bbox_t bounds = {
    .min_x = -200, .max_x = 200,
    .min_y = -200, .max_y = 200,
    .min_z = -200, .max_z = 200
};

void_result_t* voids = alea_void_generate_in_bbox(sys, &bounds);
if (voids) {
    size_t n = alea_void_count(voids);
    printf("Found %zu void regions\n", n);

    // Add them as cells to the geometry
    int added = alea_void_add_cells(sys, voids);
    printf("Added %d void cells\n", added);

    alea_void_free(voids);
}
```

If too many void cells are generated, merge them:

```c
alea_void_merge(sys, voids);  // reduce count by merging adjacent regions
```

## 8. Format Conversion

Converting between MCNP and OpenMC:

```c
// MCNP to OpenMC
mcnp_model_t* model = mcnp_load("input.inp");
openmc_export_system(mcnp_model_system(model), "geometry.xml");
mcnp_model_destroy(model);

// OpenMC to MCNP
openmc_model_t* omc = openmc_load("geometry.xml");
mcnp_export_system(openmc_model_system(omc), "output.inp");
openmc_model_destroy(omc);
```

For merging multiple models:

```c
mcnp_model_t* a = mcnp_load("model_a.inp");
mcnp_model_t* b = mcnp_load("model_b.inp");

// Merge b into a, offsetting all IDs by 100000 to avoid collisions
alea_merge(mcnp_model_system(a), mcnp_model_system(b), 100000);

mcnp_export_system(mcnp_model_system(a), "combined.inp");
mcnp_model_destroy(a);
mcnp_model_destroy(b);
```

## 9. Flattening and Manipulation

### Flattening universes

Flattening expands the universe hierarchy into a single flat universe. Every filled cell is replaced by the actual geometry it references, with transforms applied:

```c
alea_flatten(sys, 0);  // flatten universe 0
```

After flattening, all cells are in universe 0 with no fills. This is useful for exporting to codes that don't support universe hierarchies, or for simplifying a model.

### Extracting a universe

Pull one universe out into its own system:

```c
alea_system_t* sub = alea_extract_universe(sys, 5);  // universe 5
mcnp_export_system(sub, "universe_5.inp");
alea_destroy(sub);
```

### Renumbering

```c
alea_renumber_cells(sys, 1);       // cells start at 1
alea_renumber_surfaces(sys, 1);    // surfaces start at 1
alea_offset_cell_ids(sys, 10000);  // add 10000 to all cell IDs
```

## Next Steps

- Read [Concepts](CONCEPTS.md) to understand sense, universes, lattices, and other domain concepts
- Read the [API Reference](API.md) for supported C workflows, contracts, and
  advanced query families; installed public headers remain canonical.
- Look at `tools/mc_plotter.c` for a complete visualization example, and the `examples/c/` directory for other usage patterns
