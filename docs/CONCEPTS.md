# Concepts

This document explains the domain concepts that Alea works with. If you know MCNP, most of this will be familiar, but there are a few places where Alea's representation differs from what you might expect.

## Surfaces

A surface is a mathematical object that divides all of 3D space into two halves. A plane divides space into "above" and "below." A sphere divides space into "inside" and "outside." A cylinder divides space into "inside the tube" and "outside the tube."

In Alea, a surface is always one of these primitive types:

| Type | Equation | MCNP equivalent |
|------|----------|-----------------|
| Plane | ax + by + cz + d = 0 | P, PX, PY, PZ |
| Sphere | (x-cx)^2 + (y-cy)^2 + (z-cz)^2 = r^2 | S, SO, SX, SY, SZ |
| Cylinder X | (y-cy)^2 + (z-cz)^2 = r^2 | CX, C/X |
| Cylinder Y | (x-cx)^2 + (z-cz)^2 = r^2 | CY, C/Y |
| Cylinder Z | (x-cx)^2 + (y-cy)^2 = r^2 | CZ, C/Z |
| Cone X/Y/Z | ... = t^2 * (x-cx)^2 | KX, KY, KZ, K/X, K/Y, K/Z |
| Box (RPP) | axis-aligned box | RPP |
| Quadric | Ax^2 + By^2 + ... + J = 0 | GQ, SQ |
| Torus X/Y/Z | fourth-degree surface of revolution | TX, TY, TZ |

Surfaces have IDs. In MCNP, these are the numbers you put on surface cards (like `1 PZ 5.0`). Alea preserves these IDs through loading and export.

Surface creation functions follow the pattern `alea_*_surface()`:

```c
int idx = alea_sphere_surface(sys, /*surface_id=*/10, /*cx=*/0, /*cy=*/0, /*cz=*/0, /*r=*/5.0);
int idx2 = alea_plane_surface(sys, /*surface_id=*/20, /*a=*/0, /*b=*/0, /*c=*/1, /*d=*/-10.0);
```

These return a surface index (position in the surfaces array), not the surface ID.

## Sense

Every surface divides space into two sides. In Alea (and MCNP), these sides are called **positive** and **negative** sense.

- **Negative sense** (`sense = -1`): the inside. For a sphere, points closer to the center than the radius. For a plane, points on the side opposite to the normal vector. Written as `-S` in MCNP.
- **Positive sense** (`sense = +1`): the outside. For a sphere, points farther from the center than the radius. Written as `+S` or just `S` in MCNP.

To build CSG trees, you create half-space nodes from surfaces:

```c
int idx = alea_sphere_surface(sys, 10, 0, 0, 0, 5.0);

// "inside the sphere" — the region where (x-cx)^2 + ... < r^2
alea_node_id_t inside = alea_halfspace(sys, idx, -1);

// "outside the sphere" — the region where (x-cx)^2 + ... > r^2
alea_node_id_t outside = alea_halfspace(sys, idx, +1);
```

For planes, the convention is: the normal vector `(a, b, c)` points toward the positive side. So `sense = -1` means the half-space on the opposite side of the normal.

## Cells

A cell is a region of space defined by a boolean combination of surface half-spaces. In MCNP terms, a cell card like:

```
1  1  -7.8  -1  2  -3
```

defines cell 1, material 1, density 7.8 g/cc, as the intersection of: inside surface 1, outside surface 2, and inside surface 3. Alea represents this as a CSG tree node.

Every cell has:

- **Cell ID**: the MCNP cell number (unique identifier)
- **Material ID**: which material fills the cell (0 = void)
- **Density**: material density (always stored as a positive value internally, with a separate `is_mass_density` flag to distinguish g/cm^3 from atoms/barn-cm; MCNP's sign convention is applied only at export time)
- **Universe ID**: which universe the cell belongs to (default 0)
- **Region**: a CSG tree node defining the shape
- **Fill**: optionally, a universe that fills this cell instead of a material

A cell either has a material or a fill, never both. A cell with a fill is a "container" — it defines a region of space and says "look in universe N to find what's actually here."

## Boolean Operations (CSG Trees)

Regions are defined by combining surface half-spaces with boolean operations:

- **Intersection** (`a AND b`): points that are inside both a and b. In MCNP, this is the implicit operation when you write `-1 -2` (inside surface 1 AND inside surface 2).
- **Union** (`a OR b`): points that are inside a or b or both. In MCNP, written with `:` — e.g., `-1:-2`.
- **Difference** (`a AND NOT b`): points inside a but not inside b. Equivalent to `intersection(a, complement(b))`.
- **Complement** (`NOT a`): everything that is not inside a. In MCNP, `#` applied to a cell.

Alea stores these as a binary tree. Each leaf is a surface half-space (primitive + sense). Each internal node is a boolean operation with left and right children.

```
     INTERSECTION
      /        \
  -sphere    +plane
  (inside)   (above)
```

This tree represents "inside the sphere AND above the plane" — a hemisphere.

The construction API:

```c
alea_node_id_t region = alea_intersection(sys, inside_sphere, above_plane);
alea_node_id_t region2 = alea_union(sys, region_a, region_b);
alea_node_id_t region3 = alea_difference(sys, region_a, region_b);
alea_node_id_t region4 = alea_complement(sys, region_a);

// N-ary variants for combining many nodes at once
alea_node_id_t nodes[] = {a, b, c, d};
alea_node_id_t all = alea_intersection_n(sys, nodes, 4);
```

## Universes

A universe is a collection of cells that share the same coordinate system. Every cell belongs to exactly one universe. By default, everything is in universe 0.

Universes become important when you have repeated geometry. Instead of defining the same fuel pin 1000 times, you define it once in universe 5, and then fill 1000 container cells with universe 5.

### Fills

A cell with `FILL=N` is a container. It says: "the geometry inside me is defined by the cells in universe N." When Alea evaluates a point query, it:

1. Finds the container cell in the current universe
2. Sees it has `FILL=N`
3. Optionally applies the inverse of the fill transform
4. Searches for the point among the cells of universe N
5. Repeats if it hits another fill

This is recursive. Universe 5 might contain a cell that fills with universe 10, which might contain cells that fill with universe 20. Alea handles arbitrary nesting depth.

### Fill transforms

A fill can include a coordinate transform. In MCNP:

```
10 0 -1 FILL=5 (10.0 0.0 0.0)
```

This means: universe 5 is placed inside cell 10, shifted by (10, 0, 0). The same universe can be filled at different positions with different transforms, creating repeated geometry at different locations.

Transforms can include rotations (3x3 matrix) in addition to translation.

## Lattices

A lattice is an efficient way to tile space with a repeating pattern. Instead of creating thousands of container cells manually, you define one lattice cell and a fill pattern.

Alea supports two lattice types:

- **Rectangular** (`lat=1` in MCNP): grid of axis-aligned boxes
- **Hexagonal** (`lat=2` in MCNP): hexagonal prism tiling

### How lattices work

A lattice cell defines:

- The **pitch**: spacing between elements in each dimension
- The **fill pattern**: a 2D or 3D array of universe IDs
- The **range**: index bounds (e.g., -3 to 3 in X and Y)

When Alea evaluates a point inside a lattice cell, it:

1. Computes which lattice element contains the point (from the coordinates and pitch)
2. Looks up the universe ID for that element in the fill array
3. Translates coordinates to the element's local frame
4. Searches for the point in that universe

Lattice elements are never materialized in memory. The element index is computed on the fly for each query. A lattice with 10,000 elements uses the same memory as one with 4.

### Hexagonal coordinates

For hexagonal lattices, element lookup uses cube coordinate rounding:

1. Convert (x, y) to fractional hex coordinates
2. Round to nearest hex center using cube coordinate constraints
3. Map to the fill array index

## Materials

A material is a composition of nuclides or elements. Alea stores:

- **Material ID**: the MCNP material number
- **Nuclides**: list of (ZAID, fraction) pairs
- **Density**: stored on the cell, not the material

Materials are preserved through loading and export. Alea doesn't do physics with materials — it just carries them along so you can query "what material is at this point" and so exports produce valid input files.

### Mixtures

A mixture is a weighted combination of existing materials. You can create them programmatically:

```c
int mat_ids[] = {1, 2};
double fracs[] = {0.7, 0.3};
int new_id = alea_create_mixture(sys, mat_ids, fracs, 2, 100);
```

## Macrobodies

Macrobodies are composite primitives that MCNP defines as single surface cards but that are really combinations of simpler surfaces. Alea supports:

| Macrobody | Description | Decomposition |
|-----------|-------------|---------------|
| RCC | Right circular cylinder | Cylinder + 2 planes |
| BOX | General oriented box | 6 planes |
| RPP | Axis-aligned box | 6 planes (or single box primitive) |
| SPH | Sphere | Single sphere |
| TRC | Truncated right cone | Cone + 2 planes |
| ELL | Ellipsoid | Quadric surface |
| REC | Right elliptical cylinder | Quadric + 2 planes |
| WED | Wedge | 5 planes |
| RHP/HEX | Hexagonal prism | 8 planes |
| ARB | Arbitrary polyhedron | Up to 6 faces from 8 vertices |

When loading an MCNP file, macrobodies are stored as-is. You can expand them:

```c
alea_expand_macrobodies_in_system(sys);
```

This replaces each macrobody reference in the CSG tree with the equivalent boolean combination of primitive surfaces. Export can be configured to either keep macrobodies (`ALEA_EMIT_MACROBODY`) or always decompose them (`ALEA_EMIT_SURFACES`).

## Transforms

Transforms are rotation + translation operations applied to surfaces or cells.

- **Surface transforms** (MCNP `TRn` cards): modify the position/orientation of a surface. The surface is stored in global coordinates after the transform is applied.
- **Cell transforms** (MCNP `TRCL` parameter): transform the entire cell's geometry. Can be an inline translation `TRCL=(dx dy dz)` or a reference to a TR card `TRCL=n`.
- **Fill transforms**: applied when instantiating a universe into a cell. In MCNP: `FILL=5 (dx dy dz ...)`.

Alea applies surface and cell transforms during loading, so the internal representation is always in global coordinates. Fill transforms are handled at query time by inverse-transforming the point before descending into the filled universe.

## Loading and Exporting

Alea can load and export geometry in both MCNP and OpenMC formats:

```c
// Load
alea_system_t* sys = alea_load_mcnp("model.inp");
alea_system_t* sys = alea_load_openmc("geometry.xml");

// Export
alea_export_mcnp(sys, "output.inp");
alea_export_openmc(sys, "model.xml");
```

The `mc_convert` tool wraps this for command-line format conversion. See [API Reference](API.md) for export options.

## Spatial Index and Point Queries

After loading a model, you build a spatial index (BVH — bounding volume hierarchy) that accelerates all geometric queries:

```c
alea_build_spatial_index(sys);
```

Once the index is built, you can query any point in space:

```c
int cell = alea_find_cell(sys, x, y, z);       // which cell contains this point?
int mat  = alea_material_at(sys, x, y, z);      // what material is here?
```

The spatial index handles universe nesting, fills, and lattices transparently. A single query can descend through multiple universe levels to find the innermost cell.

## Ray Tracing

Alea includes a ray tracing module that casts rays through the geometry and reports every cell crossing:

```c
#include "alea_raycast.h"

alea_raycast_result_t* result = alea_raycast_result_create();
alea_raycast(sys,
    ox, oy, oz,       // ray origin
    dx, dy, dz,       // ray direction
    max_distance,      // maximum distance to trace
    result);
```

Each intersection records the cell entered, the material, the distance along the ray, and the surface crossed.

Ray tracing is also used internally for Monte Carlo volume estimation:

```c
double volumes[ncells], errors[ncells];
alea_estimate_cell_volumes(sys, ox, oy, oz, radius, n_rays, volumes, errors);
```

## 2D Slice Visualization

The slice module computes 2D cross-sections of the geometry. Given a cutting plane (axis-aligned or arbitrary), it queries a grid of points and produces cell/material ID arrays suitable for rendering:

```c
#include "alea_slice.h"

alea_slice_view_t view;
alea_slice_view_axis(&view, 2, z_value, x_min, x_max, y_min, y_max);  // Z plane

int* cell_ids = malloc(width * height * sizeof(int));
int* mat_ids  = malloc(width * height * sizeof(int));
alea_find_cells_grid(sys, &view, width, height, -1, cell_ids, mat_ids, NULL);
```

The `mc_plotter` tool wraps this into a complete plotting application with PNG/BMP output, labels, contours, and batch mode.

## 3D Rendering

The render module produces 3D images of the geometry using ray casting with Phong shading:

```c
#include "alea_render.h"

alea_render3d_params_t params = alea_render3d_default_params();
params.width = 800;
params.height = 600;
alea_render3d(sys, &params, pixels);
```

Multiple color modes are supported: by material, cell, universe, or density.

## Mesh Export

The mesh module exports the geometry as a structured hexahedral mesh for use in external tools:

```c
#include "alea_mesh.h"

alea_mesh_export_vtk(sys, &bounds, nx, ny, nz, "output.vtk");
alea_mesh_export_gmsh(sys, &bounds, nx, ny, nz, "output.msh");
```

Supported formats: Gmsh (.msh v2.2) and VTK (.vtk).

## Surface Deduplication

Large models often have many identical surfaces (e.g., the same plane defined multiple times with slightly different coefficients due to floating-point roundoff). Alea automatically deduplicates surfaces during loading.

Deduplication works by:

1. Canonicalizing each primitive (normalizing coefficients so the first non-zero coefficient is positive)
2. Hashing the canonical form
3. Checking if an equivalent primitive already exists (within tolerance)
4. If so, reusing the existing primitive instead of creating a new one

This is transparent — it doesn't change the behavior of any query. It reduces memory usage and speeds up ray tracing (fewer unique surfaces to test).

You can disable it:

```c
alea_config_t cfg = alea_get_config(sys);
cfg.dedup = false;
alea_set_config(sys, &cfg);
```

## Void Regions

A void region is a part of 3D space that no cell claims. In a well-defined MCNP model, every point should be inside exactly one cell. Void regions indicate geometry errors — or they indicate you haven't defined the outer boundary yet.

Alea's void generation uses an octree. It recursively subdivides a bounding box, probing points on a regular grid in each sub-box to determine if the region is solid (inside a cell), void (no cell claims it), or mixed (partially covered). Void sub-boxes are collected and can be added as void cells to the model.

The octree parameters are configurable:

```c
alea_config_t cfg = alea_get_config(sys);
cfg.void_max_depth = 8;    // maximum octree depth (higher = finer resolution)
cfg.void_min_size = 0.1;   // minimum sub-box size in cm
cfg.void_probes_per_axis = 3; // probe points per axis (3 → 3x3x3 grid)
alea_set_config(sys, &cfg);
```

## Boundary Conditions

Surfaces can have boundary conditions that affect particle behavior in transport simulations:

| Type | MCNP syntax | Meaning |
|------|-------------|---------|
| Transmissive | (default) | Particles cross normally |
| Reflective | `*1` | Particles reflect at the surface |
| White | `+1` | Particles re-enter with random angle |
| Periodic | paired surfaces | Particles teleport to partner surface |
| Vacuum | (graveyard) | Particles are killed |

Alea preserves boundary conditions through loading and export. They don't affect point queries or ray tracing (those are purely geometric operations), but they are important metadata for producing valid transport input files.

## Configuration

All behavior is controlled through a single `alea_config_t` struct:

```c
alea_config_t cfg = alea_get_config(sys);

// Tolerances
cfg.abs_tol = 1e-6;          // absolute tolerance for surface matching
cfg.rel_tol = 1e-9;          // relative tolerance
cfg.zero_threshold = 1e-10;  // values smaller than this are treated as zero

// Surface dedup
cfg.dedup = true;

// Export
cfg.surface_policy = ALEA_EMIT_MACROBODY;  // or ALEA_EMIT_SURFACES
cfg.export_materials = true;
cfg.export_transforms = true;
cfg.universe_depth = -1;     // -1 = all, 0 = root only
cfg.fill_depth = 0;          // 0 = don't expand fills
cfg.trcl_mode = 0;           // 0 = preserve TRCLs, 1 = bake into geometry
cfg.transform_mode = 0;      // 0 = original, 1 = inline, 2 = TR cards

// MCNP formatting
cfg.mcnp_max_col = 80;       // maximum column width
cfg.mcnp_cont_indent = 5;    // continuation line indent

// Void generation
cfg.void_max_depth = 8;
cfg.void_min_size = 0.1;
cfg.void_probes_per_axis = 3;

// Void merge (controls how void boxes are merged into cells)
cfg.merge_cell_weight = 1.0;
cfg.merge_surface_weight = 0.1;
cfg.merge_max_surfaces = 24;
cfg.merge_min_cells = 1;

// Flatten
cfg.flatten_max_depth = 0;   // 0 = unlimited

alea_set_config(sys, &cfg);
```

See the [API Reference](API.md) for the complete list of configuration fields.

## Error Handling

Most Alea functions use simple return conventions:

- Functions returning pointers return `NULL` on error
- Functions returning `int` return 0 on success, negative on error
- After an error, call `alea_error()` to get a human-readable message
- Call `alea_error_code()` to get a machine-readable error code

```c
alea_system_t* sys = alea_load_mcnp("missing.inp");
if (!sys) {
    printf("Error %d: %s\n", alea_error_code(), alea_error());
    // prints: "Error 6: File not found: missing.inp"
}
```

Error codes are defined in `alea_types.h` as the `alea_error_t` enum:

| Code | Name | Meaning |
|------|------|---------|
| 0 | `ALEA_OK` | Success |
| 1 | `ALEA_ERR_NULL_ARG` | NULL pointer where not allowed |
| 2 | `ALEA_ERR_INVALID_ID` | Node/primitive/cell ID out of range |
| 3 | `ALEA_ERR_INVALID_ARG` | Argument value out of valid range |
| 4 | `ALEA_ERR_INVALID_STATE` | Operation not valid in current state |
| 5 | `ALEA_ERR_OUT_OF_MEMORY` | Memory allocation failed |
| 6 | `ALEA_ERR_FILE_NOT_FOUND` | File does not exist |
| 7 | `ALEA_ERR_FILE_READ` | Error reading file |
| 8 | `ALEA_ERR_FILE_WRITE` | Error writing file |
| 9 | `ALEA_ERR_PARSE_ERROR` | Syntax/parse error in input file |
| 10 | `ALEA_ERR_UNSUPPORTED` | Feature not supported |
| 11 | `ALEA_ERR_UNSUPPORTED_SURFACE` | Surface type not supported |
| 12 | `ALEA_ERR_EXPORT_FAILED` | Export operation failed |
| 13 | `ALEA_ERR_NOT_IMPLEMENTED` | Feature not yet implemented |
| 14 | `ALEA_ERR_INTERRUPTED` | Operation interrupted (Ctrl+C) |
| 15 | `ALEA_ERR_NOT_FOUND` | Item not found |
| 16 | `ALEA_ERR_EMPTY` | Collection is empty |
| 17 | `ALEA_ERR_OVERFLOW` | Buffer too small, result truncated |
