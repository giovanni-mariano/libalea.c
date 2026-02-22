# Tutorial

This tutorial walks through the main things you can do with Alea. Each section is self-contained: real code, real results, no hand-waving.

All code compiles with:

```bash
gcc -o example example.c -I../include ../bin/libalea_full.a -lm
```

## 1. Loading a Model

The most common starting point is an existing MCNP input file:

```c
#include <alea.h>
#include <stdio.h>

int main(void) {
    alea_system_t* sys = alea_load_mcnp("iter_blanket.inp");
    if (!sys) {
        fprintf(stderr, "Load failed: %s\n", alea_error());
        return 1;
    }

    // Build the index — required before any queries
    alea_build_universe_index(sys);

    alea_print_summary(sys);
    alea_destroy(sys);
    return 0;
}
```

`alea_load_mcnp` parses the cell cards, surface cards, data cards (materials, transforms), and builds the internal CSG tree. It handles `LIKE BUT`, cell complements (`#cell`), macrobodies, and universe fills.

For OpenMC:

```c
alea_system_t* sys = alea_load_openmc("geometry.xml");
```

You can also load from a string instead of a file:

```c
const char* input = "1 1 -10.0 -1\n2 0 1\n\n1 SO 5.0\n\n";
alea_system_t* sys = alea_load_mcnp_string(input, strlen(input));
```

**Important**: Always call `alea_build_universe_index(sys)` after loading. Without it, point queries, overlap detection, and slicing will not work correctly.

## 2. Asking Questions About the Geometry

Once loaded, the most useful thing is asking "what's at this point?"

### Which cell contains a point?

```c
int cell = alea_find_cell(sys, 650.0, 0.0, 0.0);
if (cell > 0) {
    printf("Point is in cell %d\n", cell);
} else {
    printf("Point is in void or undefined\n");
}
```

`alea_find_cell` returns the MCNP cell ID of the innermost cell containing the point. It traverses the full universe hierarchy — if the point is in a cell that has `FILL=5`, it descends into universe 5, applies the inverse transform, and continues until it finds a terminal cell (one with a material or void, not another fill).

Returns 0 if no cell claims the point. This means either void or a geometry error.

### What material is there?

```c
int mat = alea_material_at(sys, 650.0, 0.0, 0.0);
```

Returns the material number (MCNP `m` card number). Returns 0 for void cells.

### Both at once

```c
int cell_id, material;
alea_find_cell_at(sys, 650.0, 0.0, 0.0, &cell_id, &material);
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

This samples random points within each cell's bounding box and checks if multiple cells claim the same point. It's statistical — not exhaustive — but it catches the vast majority of real overlaps.

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

The typical workflow combines both: use the grid for pixel coloring, and overlay the curves for crisp surface boundaries. The `examples/plotter.c` program does exactly this.

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
    int cell_id, material_id;
    double density;
    alea_raycast_segment_get(result, i,
        &t_enter, &t_exit, &cell_id, &material_id, &density);

    double thickness = t_exit - t_enter;
    printf("%8.2f - %8.2f cm: cell %d, mat %d, density %.3f (%.2f cm)\n",
           t_enter, t_exit, cell_id, material_id, density, thickness);
}

alea_raycast_result_free(result);
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

For large models, cell-aware tracing is faster. Instead of testing every surface globally, it tracks through cells one at a time:

```c
alea_raycast_cell_aware(sys, ox, oy, oz, dx, dy, dz, t_max, result);
```

Same interface, same result format, better performance on models with many surfaces.

## 5. Building Geometry from Scratch

You don't have to load from a file. You can build geometry programmatically:

```c
alea_system_t* sys = alea_create();

// Create surfaces and get their interior (neg_node) or exterior (pos_node)
int s1 = alea_sphere_surface(sys, 0, 0, 0, 0, 10.0);
int s2 = alea_cylinder_z_surface(sys, 0, 0, 0, 3.0);

alea_node_id_t outer = alea_surface_at(sys, s1)->neg_node;  // inside sphere
alea_node_id_t hole  = alea_surface_at(sys, s2)->neg_node;  // inside cylinder

// Boolean difference: sphere minus cylinder
alea_node_id_t region = alea_difference(sys, outer, hole);

// Add as cell 1, material 1, density 10.0, universe 0
alea_add_cell(sys, 1, region, 1, 10.0, 0);

alea_build_universe_index(sys);
```

### Creating surfaces

Each surface function returns an index. Access `neg_node` (interior) or `pos_node` (exterior):

```c
int idx = alea_plane_surface(sys, id, a, b, c, d);      // ax + by + cz + d = 0
int idx = alea_sphere_surface(sys, id, cx, cy, cz, r);
int idx = alea_cylinder_z_surface(sys, id, cx, cy, r);  // infinite along Z
int idx = alea_cylinder_x_surface(sys, id, cy, cz, r);  // infinite along X
int idx = alea_cylinder_y_surface(sys, id, cx, cz, r);  // infinite along Y
int idx = alea_box_surface(sys, id, xmin, xmax, ymin, ymax, zmin, zmax);
int idx = alea_cone_z_surface(sys, id, cx, cy, cz, t2); // t2 = tan^2(half-angle)

// Get the halfspace nodes
alea_node_id_t inside  = alea_surface_at(sys, idx)->neg_node;
alea_node_id_t outside = alea_surface_at(sys, idx)->pos_node;
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
alea_node_id_t fuel_r = alea_surface_at(sys, s_fuel)->neg_node;
alea_add_cell(sys, 10, fuel_r, 1, 10.0, 1);  // universe 1

int s_clad = alea_sphere_surface(sys, 0, 0, 0, 0, 0.6);
alea_node_id_t clad_r = alea_difference(sys,
    alea_surface_at(sys, s_clad)->neg_node,
    alea_surface_at(sys, s_fuel)->neg_node);
alea_add_cell(sys, 11, clad_r, 2, 8.0, 1);   // universe 1

// Universe 0: container that fills with universe 1
int s_box = alea_box_surface(sys, 0, -5, 5, -5, 5, -5, 5);
alea_node_id_t box = alea_surface_at(sys, s_box)->neg_node;
int cell_idx = alea_add_cell(sys, 1, box, 0, 0.0, 0);  // universe 0
alea_set_fill(sys, cell_idx, 1, 0);  // fill with universe 1, no transform
```

## 6. Exporting

### To MCNP

```c
alea_export_mcnp(sys, "output.inp");
```

### To OpenMC

```c
alea_export_openmc(sys, "geometry.xml");
```

### To a file stream

```c
FILE* f = fopen("output.inp", "w");
alea_export_mcnp_stream(sys, f);
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

alea_export_mcnp(sys, "output.inp");
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

void_result_t* voids = alea_void_generate(sys, &bounds);
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

Converting between MCNP and OpenMC is a two-liner:

```c
// MCNP to OpenMC
alea_system_t* sys = alea_load_mcnp("input.inp");
alea_export_openmc(sys, "geometry.xml");
alea_destroy(sys);

// OpenMC to MCNP
alea_system_t* sys2 = alea_load_openmc("geometry.xml");
alea_export_mcnp(sys2, "output.inp");
alea_destroy(sys2);
```

For merging multiple models:

```c
alea_system_t* a = alea_load_mcnp("model_a.inp");
alea_system_t* b = alea_load_mcnp("model_b.inp");

// Merge b into a, offsetting all IDs by 100000 to avoid collisions
alea_merge(a, b, 100000);

alea_export_mcnp(a, "combined.inp");
alea_destroy(a);
alea_destroy(b);
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
alea_export_mcnp(sub, "universe_5.inp");
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
- Read the [API Reference](API.md) for the complete function listing
- Look at the examples in the `examples/` directory, especially `plotter.c` for visualization
