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
// Load (alea_mcnp.h / alea_openmc.h)
mcnp_model_t* model = mcnp_load("model.inp");
openmc_model_t* omc = openmc_load("geometry.xml");

// Export
mcnp_export(model, "output.inp");
openmc_export(omc, "model.xml");

// Or export a bare system with defaults
mcnp_export_system(sys, "output.inp");
openmc_export_system(sys, "model.xml");
```

The `mc_convert` tool wraps this for command-line format conversion. See [API Reference](API.md) for export options.

## Spatial Index and Point Queries

After loading a model, prepare query acceleration before geometric queries:

```c
alea_prepare_query_acceleration(sys);
```

This builds the acceleration structures (the hierarchical TLAS/BLAS spatial
index, surface BVH, and related caches). Once the acceleration is prepared, you
can query any point in space:

```c
int cell_id, material_id;
alea_find_cell_at(sys, x, y, z, &cell_id, &material_id);
```

The spatial index handles universe nesting, fills, and lattices transparently. A single query can descend through multiple universe levels to find the innermost cell.

## Ray Tracing

Alea includes ray tracing that casts rays through the geometry and reports every cell crossing:

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

Alea computes 2D cross-sections of the geometry. Given a cutting plane (axis-aligned or arbitrary), it queries a grid of points and produces cell/material ID arrays suitable for rendering:

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

Alea produces 3D images of the geometry using ray casting with Phong shading:

```c
#include "alea_render.h"

render_config_t cfg;
render_config_init(&cfg);
cfg.width = 800;
cfg.height = 600;

render_camera_t cam;
render_camera_setup(&cam, &cfg, sys);

render_framebuffer_t* fb = render_framebuffer_create(cfg.width, cfg.height, 0);
render_scene(sys, &cfg, &cam, fb);
render_write_png("output.png", fb);
render_framebuffer_free(fb);
```

Multiple color modes are supported: by material, cell, universe, or density.

## Mesh Export

Alea exports the geometry as a structured hexahedral mesh for use in external tools:

```c
#include "alea_mesh.h"

alea_mesh_config_t cfg;
alea_mesh_config_init(&cfg);
cfg.nx = cfg.ny = cfg.nz = 100;

alea_mesh_export_system(sys, &cfg, "output.msh");  // format auto-detected from extension
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
cfg.abs_tol = 1e-6;            // absolute tolerance for surface matching
cfg.rel_tol = 1e-9;            // relative tolerance
cfg.zero_threshold = 1e-10;    // values smaller than this are treated as zero

// Behavior
cfg.dedup = true;              // surface deduplication
cfg.log_level = 2;             // 0=none, 1=error, 2=warn, 3=info, 4=debug

// Export
cfg.export_materials = true;
cfg.export_transforms = true;
cfg.universe_depth = -1;       // -1 = all, 0 = root only
cfg.fill_depth = 0;            // 0 = don't expand fills

// Void generation
cfg.void_max_depth = 8;
cfg.void_min_size = 0.1;
cfg.void_probes_per_axis = 3;

// Void merge (controls how void boxes are merged into cells)
cfg.merge_cell_weight = 1.0;
cfg.merge_surface_weight = 0.1;
cfg.merge_max_surfaces = 24;
cfg.merge_min_cells = 1;
cfg.merge_use_greedy = false;
cfg.void_consolidate = 100;    // 0 = off

// Flatten
cfg.flatten_max_depth = 0;     // 0 = unlimited

alea_set_config(sys, &cfg);
```

MCNP-specific export options (surface policy, column width, continuation indent, TRCL mode, transform mode) are configured separately via `mcnp_export_config_t` in `alea_mcnp.h`.

See the [API Reference](API.md) for the complete list of configuration fields.

## Error Handling

Most Alea functions use simple return conventions:

- Functions returning pointers return `NULL` on error
- Functions returning `int` return 0 on success, negative on error
- After an error, call `alea_error()` to get a human-readable message
- Call `alea_error_code()` to get a machine-readable error code

```c
mcnp_model_t* model = mcnp_load("missing.inp");
if (!model) {
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

---

## Nuclear Data Concepts

This section explains the nuclear physics and data format concepts behind the nuclear data module (`libalea_nucdata.a`).

### The ACE Format

ACE (A Compact ENDF) is the standard data format for Monte Carlo transport codes. It originated at Los Alamos for MCNP and is now used by OpenMC, Serpent, and others. An ACE file encodes everything a transport code needs about a nuclide: energy grids, cross sections, angular distributions, energy distributions, fission data, and more.

Each ACE table is identified by a ZAID string like `92235.80c`:

| Field | Example | Meaning |
|-------|---------|---------|
| `92235` | ZZAAA | Z=92 (uranium), A=235 |
| `.80c` | .LLt | Library 80, continuous-energy neutron |

The suffix letter identifies the table type:

| Suffix | Type | Description |
|--------|------|-------------|
| `.c` | Continuous-energy neutron | Full energy-dependent neutron interactions |
| `.p` | Photoatomic | Photon interactions with atoms |
| `.u` | Photonuclear | Photon-induced nuclear reactions |
| `.t` | Thermal S(α,β) | Bound-atom scattering at low energies |
| `.e` | Electron | Electron/positron interactions |

The module supports `.c` and `.p` tables.

Each ACE table has three index structures:

- **NXS[16]**: Table dimensions — energy grid size, number of reactions, etc.
- **JXS[32]**: Block locators — byte offsets into the XSS data array.
- **XSS[N]**: The data itself — a flat array of doubles containing all physics data.

All indices are **1-based** (Fortran convention). The library converts to 0-based internally.

Key blocks for neutron tables:

| Block | JXS index | Contents |
|-------|-----------|----------|
| ESZ | JXS[1] | Energy grid, σ_total, σ_abs, σ_elastic, heating |
| NU | JXS[2] | ν̄ (neutrons per fission) |
| MTR | JXS[3] | MT reaction numbers for non-elastic reactions |
| SIG | JXS[7] | Non-elastic reaction cross sections |
| AND | JXS[8] | Angular distributions |
| DLW | JXS[11] | Energy distributions (outgoing spectra) |
| URR | JXS[23] | Unresolved resonance probability tables |

### The xsdir File

A transport code doesn't read ACE files directly by name. Instead, an **xsdir** (cross-section directory) file maps ZAIDs to file locations:

```
92235.80c  235.043930  endf80/092_U_235.ace  0  1  1  85947  0  0  2.5301E-08
```

Fields: ZAID, atomic weight ratio, filename, access route, file type (1=ASCII, 2=binary), address, table length, record length, entries per record, temperature (kT in MeV).

### Cross Sections

A **cross section** σ(E) gives the probability of a nuclear reaction per unit path length per target atom. Units: barns (1 barn = 10⁻²⁴ cm²).

The ACE ESZ block provides four principal cross sections on a common energy grid:

- **σ_total**: Sum of all interactions. Determines collision rate.
- **σ_absorption**: Reactions that remove the neutron (capture, fission, etc.).
- **σ_elastic**: Elastic scattering — neutron bounces off nucleus.
- **heating**: Energy deposited per collision (MeV-barn).

Each non-elastic reaction has an MT number from the ENDF convention:

| MT | Reaction | Class |
|----|----------|-------|
| 2 | Elastic scattering | Scatter |
| 16 | (n,2n) | Multiply |
| 17 | (n,3n) | Multiply |
| 18 | Total fission | Multiply |
| 51-91 | Inelastic levels and continuum | Scatter |
| 102 | Radiative capture (n,γ) | Absorption |
| 103 | (n,p) proton emission | Absorption |
| 107 | (n,α) | Absorption |

Reaction cross sections are stored on **sub-grids**: each reaction's XS array starts at its energy threshold and shares indices with the main grid above that threshold.

#### Interpolation

Neutron cross sections use **linear-linear** interpolation between grid points. Photon cross sections use **log-log** interpolation:

```
Linear:   σ(E) = σᵢ + f · (σᵢ₊₁ - σᵢ),      f = (E - Eᵢ) / (Eᵢ₊₁ - Eᵢ)
Log-log:  σ(E) = σᵢ · (σᵢ₊₁/σᵢ)^g,           g = ln(E/Eᵢ) / ln(Eᵢ₊₁/Eᵢ)
```

### Elastic Scattering Kinematics

When a neutron of energy E scatters elastically off a nucleus with atomic weight ratio A:

- **Minimum outgoing energy**: E_out = α·E, where α = ((A−1)/(A+1))²
- **Maximum outgoing energy**: E_out = E (glancing collision)

For hydrogen (A=1), α = 0 — a neutron can lose all its energy in one collision. For uranium (A=238), α = 0.983 — very little energy loss per collision.

### Angular Distributions

The AND block provides angular distributions at discrete incident energies in three formats:

1. **Isotropic**: μ_cm uniform on [−1, 1]
2. **Equiprobable bins**: 32 cosine values dividing [−1, 1] into equal-probability bins
3. **Tabular**: Full PDF and CDF on a cosine grid with CDF inversion

At runtime, the library uses stochastic interpolation between the two bracketing incident energies.

### Energy Distributions

The DLW block encodes outgoing energy distributions using several "laws":

| Law | Name | Use |
|-----|------|-----|
| 3 | Level scattering | Discrete inelastic levels |
| 4 | Continuous tabular | General-purpose tabulated spectra |
| 7 | Maxwell fission | Fission spectrum (simple) |
| 9 | Evaporation | Neutron evaporation from compound nucleus |
| 11 | Watt fission | Watt spectrum: χ(E) ∼ exp(−E/a) · sinh(√(bE)) |
| 44 | Kalbach-Mann | Tabular energy with Kalbach angular |
| 61 | Correlated energy-angle | Tabular with per-bin angular distributions |
| 66 | N-body phase space | Multi-body breakup |

### Fission

Fissile nuclides have special data:

- **ν̄(E)**: average neutrons per fission. ~2.4 for U-235 thermal, rising to ~4 at 10 MeV. Stored as polynomial or tabular.
- **Fission spectrum χ(E)**: energy distribution of emitted neutrons. Common model is the Watt spectrum: χ(E) ∼ exp(−E/a) · sinh(√(bE)).

### Photon Interactions

Photoatomic tables (`.p` suffix) contain four interaction types:

- **Compton (incoherent) scattering — MT 504**: photon scatters off electron, losing energy. Klein-Nishina cross section with composition-rejection sampling.
- **Rayleigh (coherent) scattering — MT 502**: elastic scattering off atom, no energy change. Angular distribution from atomic form factor F(q,Z).
- **Photoelectric absorption — MT 522**: photon fully absorbed. Dominant below ~100 keV.
- **Pair production — MT 516**: above 1.022 MeV, photon converts to electron-positron pair. Two 0.511 MeV annihilation photons created.

### Unresolved Resonance Region (URR)

In the unresolved resonance region (typically 1 keV to 1 MeV for heavy nuclides), individual resonances are too dense to resolve. **Probability tables** encode the statistical distribution of cross-section values. A random number selects a band giving multiplicative factors for total, elastic, fission, capture, and heating cross sections, preserving self-shielding effects.

### Doppler Broadening

Cross sections at temperature T₀ must be adjusted for temperature T > T₀. The thermal motion of target nuclei smears resonance peaks:

```
σ_D(E) = 1/(y√π) · ∫ y'·σ(E')·[exp(−(y'−y)²) − exp(−(y'+y)²)] dy'
```

where y = √(AWR·E/ΔkT) and ΔkT = kT_new − kT_old. Broadening can only increase temperature.

### Multigroup Cross Sections

Multigroup methods discretize the energy variable into G groups with boundaries in descending order. Group-averaged cross sections:

```
σ_g = ∫σ(E)·φ(E)dE / ∫φ(E)dE
```

The module uses 1/E weighting (φ(E) ∼ 1/E). The scattering transfer matrix, fission spectrum χ[g], and adjoint scattering matrix are all computed during collapse.

### Nuclear Materials

A material is a mixture of nuclides, each with a number density Nᵢ in atoms/barn-cm. Macroscopic cross section: Σ(E) = Σᵢ Nᵢ · σᵢ(E). Mean free path: λ = 1/Σ. Distance to next collision: s = −ln(1−ξ)/Σ_total(E).
