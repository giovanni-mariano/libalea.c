# Lua Tutorial

This tutorial covers the Lua bindings for Alea. The Lua API provides a concise, scriptable interface for geometry manipulation, analysis, and format conversion.

Run scripts with:

```bash
bin/alea script.lua
```

Pass arguments after the script name; access them via `alea.arg[1]`, `alea.arg[2]`, etc.

## 1. Getting Started

### Creating and inspecting a system

```lua
local sys = alea.create()
print("Cells:    " .. sys:cell_count())
print("Surfaces: " .. sys:surface_count())
```

### Loading from files

```lua
-- MCNP input
local sys = alea.load_mcnp("model.inp")

-- OpenMC XML
local sys = alea.load_openmc("geometry.xml")

-- From a string
local mcnp_input = [[
Test model
1  1  -10.0  -1
2  0          1

1 SO 5.0
]]
local sys = alea.load_mcnp_string(mcnp_input)
```

### Configuration

```lua
local cfg = sys:get_config()
print("abs_tol: " .. cfg.abs_tol)
print("dedup:   " .. tostring(cfg.dedup))

sys:set_config({
    abs_tol   = 1e-8,
    dedup     = false,
    log_level = 0,
})
```

### Cleanup

```lua
sys:destroy()
```

## 2. Building Geometry

Create surfaces, get halfspaces, combine with boolean operators, and add cells.

```lua
local sys = alea.create()

-- Create surfaces (returns surface index)
local s1 = sys:sphere(1, 0, 0, 0, 5)    -- id=1, centered at origin, r=5
local s2 = sys:sphere(2, 0, 0, 0, 6)    -- id=2, r=6
local s3 = sys:sphere(3, 0, 0, 0, 15)   -- id=3, r=15

-- Get halfspace nodes
local in_s1  = sys:inside(s1)   -- interior of sphere 1
local in_s2  = sys:inside(s2)
local in_s3  = sys:inside(s3)
local out_s1 = sys:outside(s1)  -- exterior of sphere 1
local out_s2 = sys:outside(s2)
local out_s3 = sys:outside(s3)

-- Boolean operations use Lua operators:
--   *  intersection (AND)
--   +  union (OR)
--   -  difference (A AND NOT B)
--   ~  complement (NOT)
local fuel_region = in_s1                -- r < 5
local clad_region = out_s1 * in_s2       -- 5 <= r < 6
local mod_region  = out_s2 * in_s3       -- 6 <= r < 15
local void_region = out_s3               -- r >= 15

-- Create cells
sys:cell{id = 1, region = fuel_region, material = 1, density = 10.97}
sys:cell{id = 2, region = clad_region, material = 2, density = 6.56}
sys:cell{id = 3, region = mod_region,  material = 3, density = 1.0}
sys:cell{id = 4, region = void_region, material = 0, density = 0.0}

sys:print_summary()
```

### Available surfaces

```lua
sys:plane(id, a, b, c, d)                      -- ax + by + cz + d = 0
sys:sphere(id, cx, cy, cz, r)
sys:cylinder_z(id, cx, cy, r)                  -- infinite along Z
sys:cylinder_x(id, cy, cz, r)                  -- infinite along X
sys:cylinder_y(id, cx, cz, r)                  -- infinite along Y
sys:box(id, xmin, xmax, ymin, ymax, zmin, zmax)
sys:cone_z(id, cx, cy, cz, t2)                 -- t2 = tan^2(half-angle)
```

Pass `id=0` for automatic surface ID assignment.

### Universe fills

```lua
-- Universe 1: fuel pin
local s_fuel = sys:sphere(0, 0, 0, 0, 0.5)
local s_clad = sys:sphere(0, 0, 0, 0, 0.6)

sys:cell{id = 10, region = sys:inside(s_fuel), material = 1, density = 10.0, universe = 1}
sys:cell{id = 11, region = sys:outside(s_fuel) * sys:inside(s_clad), material = 2, density = 8.0, universe = 1}

-- Universe 0: container filled with universe 1
local s_box = sys:box(0, -5, 5, -5, 5, -5, 5)
sys:cell{id = 1, region = sys:inside(s_box), material = 0, density = 0.0, universe = 0, fill = 1}
```

## 3. Point Queries

Before querying, build the universe index:

```lua
sys:build_universe_index()
```

### Find cell at a point

```lua
local cell_id, material_id = sys:find_cell_at(x, y, z)
if cell_id then
    local info = sys:cell_find_info(cell_id)
    print("Cell ID: " .. info.cell_id)
    print("Material: " .. material_id)
else
    print("No cell found (void or undefined)")
end
```

### Convenience wrappers

```lua
local idx = sys:find_cell(x, y, z)      -- legacy convenience wrapper
local mat = sys:material_at(0, 0, 0)
```

For new Lua code, prefer `sys:find_cell_at(x, y, z)`.

### Full hierarchy traversal

```lua
local hits = sys:find_all_cells(0, 0, 0)
for i, h in ipairs(hits) do
    print(string.format("depth %d: cell %d, universe %d, mat %d",
        h.depth, h.cell_id, h.universe_id, h.material_id))
    if h.fill_universe > 0 then
        print("  -> fills universe " .. h.fill_universe)
    end
end
```

### Test point against a region

```lua
local info = sys:cell_info(0)
local inside = sys:point_inside(info.root, x, y, z)
```

## 4. Model Inspection

### Cell information

```lua
local nc = sys:cell_count()
for i = 0, nc - 1 do
    local info = sys:cell_info(i)
    print(string.format("Cell %d: mat=%d universe=%d",
        info.cell_id, info.material_id, info.universe_id))
end
```

The `cell_info` table contains:
- `cell_id` - MCNP cell ID
- `material_id` - material number
- `density` - cell density
- `universe_id` - universe this cell belongs to
- `fill_universe` - universe this cell fills (0 if none)
- `root` - CSG tree root node
- `bbox` - bounding box (`min_x`, `max_x`, `min_y`, `max_y`, `min_z`, `max_z`)

### CSG tree traversal

```lua
local function walk_tree(sys, node, depth)
    local indent = string.rep("  ", depth)
    local op = sys:node_operation(node)

    if op == "primitive" then
        local sense = sys:node_sense(node)
        local sid = sys:node_surface_id(node)
        print(indent .. (sense > 0 and "+" or "-") .. sid)
    elseif op == "complement" then
        print(indent .. "COMPLEMENT")
        walk_tree(sys, sys:node_left(node), depth + 1)
    else
        print(indent .. string.upper(op))
        walk_tree(sys, sys:node_left(node), depth + 1)
        walk_tree(sys, sys:node_right(node), depth + 1)
    end
end

local root = sys:cell_info(0).root
walk_tree(sys, root, 0)
```

### Find cell by ID

```lua
local idx = sys:cell_find(cell_id)
```

## 5. Overlap Detection

Build both indices first:

```lua
sys:prepare_query_acceleration()
```

Find overlapping cell pairs:

```lua
local overlaps = sys:find_overlaps()
print("Found " .. #overlaps .. " overlapping pairs")

for i, pair in ipairs(overlaps) do
    local a, b = pair[1], pair[2]
    local info_a = sys:cell_info(a)
    local info_b = sys:cell_info(b)
    print(string.format("Overlap: cell %d <-> cell %d",
        info_a.cell_id, info_b.cell_id))
end
```

## 6. Volume Estimation

Monte Carlo volume estimation using ray tracing:

```lua
sys:build_universe_index()

-- Get bounding sphere
local cx, cy, cz, r = sys:bounding_sphere(1e-3)
print(string.format("Bounding sphere: center (%.2f, %.2f, %.2f), radius %.2f",
    cx, cy, cz, r))

-- Estimate volumes
local n_rays = 10000
local volumes = sys:estimate_volumes(n_rays)

for i, v in ipairs(volumes) do
    if v.volume > 0 then
        print(string.format("Cell %d, path %d: volume %.2f (error %.1f%%)",
            v.cell_id, v.path_id, v.volume, v.rel_error * 100))
    end
end
```

## 7. Format Conversion

### MCNP to OpenMC

```lua
local sys = alea.load_mcnp("input.inp")
sys:export_openmc("geometry.xml")
```

### OpenMC to MCNP

```lua
local sys = alea.load_openmc("geometry.xml")
sys:export_mcnp("output.inp")
```

### Round-trip verification

```lua
local sys1 = alea.load_mcnp("original.inp")
sys1:export_openmc("/tmp/converted.xml")

local sys2 = alea.load_openmc("/tmp/converted.xml")
sys2:export_mcnp("/tmp/roundtrip.inp")

print(string.format("Original: %d cells", sys1:cell_count()))
print(string.format("After round-trip: %d cells", sys2:cell_count()))
```

## 8. Flattening and Simplification

### Flatten universes

Expands all fill hierarchies into a single flat universe:

```lua
print("Before: " .. sys:universe_count() .. " universes")
sys:flatten(0)
print("After: " .. sys:universe_count() .. " universes")
```

### Tighten bounding boxes

```lua
sys:tighten_all_bboxes(1e-3)
```

### Expand macrobodies

Decomposes BOX, RCC, etc. into primitive surfaces:

```lua
local count = sys:expand_macrobodies()
print("Expanded " .. count .. " macrobodies")
```

### Split union cells

Cells with top-level unions become multiple cells:

```lua
local before = sys:cell_count()
sys:split_union_cells()
print("Cells: " .. before .. " -> " .. sys:cell_count())
```

## 9. Universe Operations

### Extract a universe

Pull one universe into its own system:

```lua
local extracted = sys:extract_universe(5)
extracted:export_mcnp("universe_5.inp")
```

### Clone a system

```lua
local copy = sys:clone()
```

### Merge systems

```lua
local sys_a = alea.load_mcnp("model_a.inp")
local sys_b = alea.load_mcnp("model_b.inp")

-- Merge b into a (offset=0 means no ID renumbering)
sys_a:merge(sys_b, 0)
```

### Renumbering

```lua
sys:renumber_cells(1)      -- cells start at 1
sys:renumber_surfaces(1)   -- surfaces start at 1
```

### ID offsets

```lua
sys:offset_cell_ids(10000)
sys:offset_surface_ids(10000)
sys:offset_material_ids(100)
```

## 10. Parametric Geometry

Lua's scripting capabilities make it easy to generate parametric geometry:

```lua
local params = {
    fuel_radius  = 0.4,
    clad_radius  = 0.5,
    pitch        = 1.26,
    nx           = 3,
    ny           = 3,
}

local sys = alea.create()
local cell_id = 1
local surf_id = 1

local x0 = -((params.nx - 1) / 2) * params.pitch
local y0 = -((params.ny - 1) / 2) * params.pitch

for iy = 0, params.ny - 1 do
    for ix = 0, params.nx - 1 do
        local cx = x0 + ix * params.pitch
        local cy = y0 + iy * params.pitch

        local s_fuel = sys:cylinder_z(surf_id, cx, cy, params.fuel_radius)
        surf_id = surf_id + 1
        local s_clad = sys:cylinder_z(surf_id, cx, cy, params.clad_radius)
        surf_id = surf_id + 1

        sys:cell{id = cell_id, region = sys:inside(s_fuel), material = 1, density = 10.97}
        cell_id = cell_id + 1

        sys:cell{id = cell_id, region = sys:outside(s_fuel) * sys:inside(s_clad), material = 2, density = 6.56}
        cell_id = cell_id + 1
    end
end

sys:export_mcnp("/tmp/pin_array.inp")
```

## 11. Validation

Check geometry for issues:

```lua
local issues = sys:validate()
if issues == 0 then
    print("Geometry is valid")
else
    print("Found " .. issues .. " issues")
end
```

## 12. Ray Tracing

Cast rays through geometry to find cell intersections:

```lua
sys:build_universe_index()

-- Cast a ray from (-10,0,0) in direction (1,0,0) with max distance 100
local result = sys:raycast(-10, 0, 0, 1, 0, 0, 100)

-- Inspect segments
print("Segments: " .. result:segment_count())
for i = 1, result:segment_count() do
    local seg = result:segment(i)  -- 1-based indexing
    print(string.format("  [%.2f, %.2f] cell=%d mat=%d density=%.2f",
        seg.t_enter, seg.t_exit, seg.cell_id, seg.material_id, seg.density))
end

-- Or get all segments at once
local segs = result:segments()
for i, seg in ipairs(segs) do
    -- same fields: t_enter, t_exit, cell_id, material_id, density
end

-- Path length through a specific material (-1 = all)
local pl = result:path_length(1)     -- total path in material 1
local total = result:path_length(-1) -- total path through all materials
```

### Cell-aware raycasting

Uses per-cell surface indices for better efficiency:

```lua
local result = sys:raycast_cell_aware(-10, 0, 0, 1, 0, 0, 100)
```

### Finding the first cell

```lua
local cell_id, t = sys:ray_first_cell(-10, 0, 0, 1, 0, 0, 100)
if cell_id then
    print("First cell: " .. cell_id .. " at distance " .. t)
end
```

## 13. 2D Slicing

Create 2D cross-sections of your geometry for visualization.

### Slice views

```lua
-- Axis-aligned: axis (0=X, 1=Y, 2=Z), value, u_min, u_max, v_min, v_max
local view = alea.slice_view_axis(2, 0.0, -10, 10, -10, 10)

-- Arbitrary orientation: origin, normal, up-hint, viewport bounds
local view = alea.slice_view(0,0,0,  0,0,1,  0,1,0,  -10,10, -10,10)
```

### Grid cell queries

```lua
sys:build_universe_index()

local grid = sys:find_cells_grid(view, 256, 256)
-- grid.cell_ids      -- flat table (256*256 elements)
-- grid.material_ids  -- flat table
-- grid.errors        -- flat table (0=ok, 1=overlap, 2=undefined)

-- With depth control: -1=innermost, 0=root level, N=level N
local grid = sys:find_cells_grid(view, 256, 256, 0)
```

### Analytical curves

```lua
local curves = sys:get_slice_curves(view)
print("Curves: " .. curves:count())
for i = 1, curves:count() do
    local c = curves:get(i)  -- 1-based
    print(c.type .. " surface=" .. c.surface_id)
    -- c.data contains type-specific fields:
    --   line: px, py, dx, dy
    --   circle/arc: cx, cy, radius
    --   ellipse: cx, cy, semi_a, semi_b, angle
    --   polygon: count, closed, vertices
end

local u0, u1, v0, v1 = curves:bounds()
```

### Label positions

```lua
-- For cell/material regions
local labels = alea.find_label_positions(grid.cell_ids, 256, 256, 100)
for _, lbl in ipairs(labels) do
    print(string.format("ID %d at pixel (%d, %d), size=%d px",
        lbl.id, lbl.px, lbl.py, lbl.pixel_count))
end

-- For surface labels along curves
local slabels = alea.find_surface_label_positions(curves,
    -10, 10, -10, 10, 256, 256, 10)
```

### Debug tracing

```lua
alea.slice_curve_set_debug(true)       -- verbose curve generation
alea.slice_point_trace_set_debug(true) -- verbose point queries
```

## 14. 3D Rendering

Produce publication-quality 3D images of CSG models.

### Basic rendering

```lua
sys:build_universe_index()

local fb = sys:render{
    width  = 800,
    height = 600,
    fov    = 45,
    shadows = 1,
    edges   = 1,
    aa_samples = 2,
}

fb:edge_darken()
fb:write("/tmp/render.png")
```

### Render configuration

```lua
-- Get default config
local cfg = alea.render_config()

-- All available options:
local fb = sys:render{
    width  = 1920,
    height = 1080,

    -- Camera
    eye    = {50, 50, 50},
    target = {0, 0, 0},
    up     = {0, 0, 1},
    fov    = 45,         -- 0 = orthographic
    ortho_height = 20,   -- used when fov=0

    -- Appearance
    color_mode  = 0,     -- 0=material, 1=cell, 2=universe, 3=density
    render_mode = 0,     -- 0=solid, 1=xray, 2=depth, 3=cellid, 4=matid
    shadows = 1,
    edges   = 1,
    ambient = 0.2,
    diffuse = 0.6,

    -- Clipping planes: visible where n.x + d > 0
    clips = {{0, 0, 1, 0}},  -- clip below z=0

    -- Quality
    aa_samples = 2,      -- NxN supersampling
    universe_depth = -1,  -- -1=innermost

    -- Output
    aux_output = 1,      -- enable depth/cellid/matid/normal maps
}
```

### Framebuffer methods

```lua
fb:width()               -- image width
fb:height()              -- image height
fb:edge_darken()         -- post-process edge darkening
fb:write("file.png")     -- auto-detect format from extension
fb:write_png("file.png")
fb:write_bmp("file.bmp")
fb:write_ppm("file.ppm")
fb:write_aux("base")     -- write auxiliary maps (depth, cellid, matid, normal)
```

### Camera setup

```lua
local cam = sys:render_camera_setup{
    eye = {50, 50, 50},
    target = {0, 0, 0},
}
print("eye:", cam.eye[1], cam.eye[2], cam.eye[3])
print("auto_fit:", cam.auto_fit)
```

### Custom colors

```lua
local cfg = alea.render_load_colors("palette.txt", {width=800, height=600})
local fb = sys:render(cfg)
```

## 15. Mesh Export

Sample CSG geometry onto a structured grid and export for mesh-based codes.

### One-shot export

```lua
sys:build_universe_index()
sys:mesh_export({
    nx = 100, ny = 100, nz = 100,
    format = 0,  -- 0=Gmsh, 1=VTK
}, "output.msh")
```

### Step-by-step with inspection

```lua
local mesh = sys:mesh_sample{
    nx = 50, ny = 50, nz = 50,
    x_min = -10, x_max = 10,
    y_min = -10, y_max = 10,
    z_min = -10, z_max = 10,
    void_material_id = 0,
    sampling_mode = 4,        -- 0=center, 1=corners, 2=regular, 3=stratified, 4=adaptive
    subsamples_per_axis = 2,
    target_error = 0.05,
    max_refine_depth = 3,
    max_samples_per_voxel = 32768,
    sampling_seed = 12345,
    workers = 1,              -- USE_TINYPAR build: 0=backend default
    bounds_mode = 0,          -- 0=legacy, 1=auto root AABB, 2=explicit
    fields = 511,             -- ALEA_MESH_FIELD_* mask; 511=current complete result
}

-- Inspect results
local info = mesh:info()
print(string.format("Grid: %dx%dx%d", info.nx, info.ny, info.nz))
print("Materials found: " .. info.num_materials)
print("Mixed voxels: " .. info.mixed_count)
for _, mid in ipairs(info.unique_materials) do
    print("  material " .. mid)
end

-- Access raw arrays
local mids = mesh:material_ids()  -- flat table (nx*ny*nz)
local cids = mesh:cell_ids()      -- flat table
local counts = mesh:sample_counts()
local ties = mesh:tie_flags()
local errors = mesh:estimated_errors()
local refinement = mesh:refinement_flags()
local fractions = mesh:material_fractions(1) -- fractions for voxel 1

-- Export to different formats
mesh:export(0, "output.msh")  -- Gmsh
mesh:export(1, "output.vtk")  -- VTK

-- Opt in to per-material sampled-fraction arrays (bit 16).
-- The default diagnostic field mask is 111; 127 includes every current field.
mesh:export(1, "composition.vtk", {
    export_fields = 127,
    max_fraction_materials = 32,
})

-- A separate, nonconforming adaptive octree.
local grid = sys:adaptive_grid_sample{
    nx = 2, ny = 2, nz = 2,
    max_grid_depth = 3,
    max_cells = 100000,
    target_error = 0.05,
}
print(grid:info().leaf_count)
local cells = grid:cells() -- stable id/parent_id/child_ids and sampled metadata
grid:export(1, "adaptive.vtk")
```

## 16. Error Handling and Logging

### Error state

```lua
-- Check last error
local msg = alea.error()       -- error message string
local code = alea.error_code() -- numeric error code
alea.error_clear()             -- clear error state
```

### Log levels

```lua
-- 0=none, 1=error, 2=warn (default), 3=info, 4=debug, 5=trace
alea.log_set_level(0)    -- suppress all output
alea.log_set_level(3)    -- show info and above

local level = alea.log_get_level()
```

### Debug tracing

```lua
alea.set_debug_trace(true)   -- enable point query tracing
```

### System reset

Clear a system and start over:

```lua
sys:reset()
assert(sys:cell_count() == 0)
```

### Simplify CSG cells

```lua
local stats = sys:simplify_all()
print("Nodes: " .. stats.nodes_before .. " -> " .. stats.nodes_after)
print("Empty cells removed: " .. stats.empty_cells_removed)
```

## Quick Reference

| Operation | Lua |
|-----------|-----|
| **Lifecycle** | |
| Create system | `alea.create()` |
| Destroy | `sys:destroy()` |
| Reset | `sys:reset()` |
| Clone | `sys:clone()` |
| Version | `alea.version()` |
| **Loading** | |
| Load MCNP | `alea.load_mcnp(file)` |
| Load OpenMC | `alea.load_openmc(file)` |
| Load from string | `alea.load_mcnp_string(str)` |
| **Export** | |
| Export MCNP | `sys:export_mcnp(file)` |
| Export OpenMC | `sys:export_openmc(file)` |
| **CSG Construction** | |
| Inside surface | `sys:inside(surf_idx)` |
| Outside surface | `sys:outside(surf_idx)` |
| Intersection | `a * b` |
| Union | `a + b` |
| Difference | `a - b` |
| Complement | `~a` |
| Add cell | `sys:cell{id=N, region=node, ...}` |
| **Indexing** | |
| Prepare query acceleration | `sys:prepare_query_acceleration()` |
| Query acceleration stats | `sys:query_acceleration_stats()` |
| Build universe index | `sys:build_universe_index()` |
| **Point Queries** | |
| Find cell + material | `sys:find_cell_at(x, y, z)` |
| Find cell | `sys:find_cell(x, y, z)` |
| Material at point | `sys:material_at(x, y, z)` |
| All cells at point | `sys:find_all_cells(x, y, z)` |
| Point inside node | `sys:point_inside(node, x, y, z)` |
| **Information** | |
| Cell count | `sys:cell_count()` |
| Surface count | `sys:surface_count()` |
| Universe count | `sys:universe_count()` |
| Cell info by index | `sys:cell_info(idx)` |
| Cell info by ID | `sys:cell_find_info(cell_id)` |
| Cell find | `sys:cell_find(cell_id)` |
| Surface info | `sys:surface_info(idx)` |
| Surface find | `sys:surface_find(surface_id)` |
| Universe info | `sys:universe_info(idx)` |
| Universe find | `sys:universe_find(universe_id)` |
| Statistics | `sys:stats()` |
| Print summary | `sys:print_summary()` |
| **Node Inspection** | |
| Node operation | `sys:node_operation(node)` |
| Node left/right | `sys:node_left(node)`, `sys:node_right(node)` |
| Node sense | `sys:node_sense(node)` |
| Node surface ID | `sys:node_surface_id(node)` |
| Primitive type | `sys:node_primitive_type(node)` |
| Primitive ID | `sys:node_primitive_id(node)` |
| Primitive data | `sys:node_primitive_data(node)` |
| Tree print | `sys:tree_print(node)` |
| **Transforms** | |
| Flatten | `sys:flatten(universe)` |
| Simplify all cells | `sys:simplify_all()` |
| Split union cells | `sys:split_union_cells()` |
| Expand macrobodies | `sys:expand_macrobodies()` |
| Extract universe | `sys:extract_universe(uid)` |
| Extract region | `sys:extract_region(bbox)` |
| Merge | `sys:merge(other, offset)` |
| Create mixture | `sys:create_mixture(ids, fracs, new_id)` |
| **Renumbering** | |
| Renumber cells | `sys:renumber_cells(start)` |
| Renumber surfaces | `sys:renumber_surfaces(start)` |
| Offset cell IDs | `sys:offset_cell_ids(offset)` |
| Offset surface IDs | `sys:offset_surface_ids(offset)` |
| Offset material IDs | `sys:offset_material_ids(offset)` |
| **Volume & BBox** | |
| Bounding sphere | `sys:bounding_sphere(tol)` |
| Estimate physical volumes | `sys:estimate_volumes(n)` |
| Volume paths | `sys:volume_paths()` |
| Tighten bbox | `sys:tighten_cell_bbox(idx, tol)` |
| Tighten all | `sys:tighten_all_bboxes(tol)` |
| Cells in bbox | `sys:cells_in_bbox(bbox)` |
| **Validation** | |
| Validate | `sys:validate()` |
| Find overlaps | `sys:find_overlaps()` |
| **Raycast** | |
| Cast ray | `sys:raycast(ox,oy,oz,dx,dy,dz,t_max)` |
| Cell-aware raycast | `sys:raycast_cell_aware(...)` |
| First cell on ray | `sys:ray_first_cell(...)` |
| First visible material | `sys:first_visible(...[, options])` |
| Ownership boundary events | `sys:boundary_events(...[, options])` |
| Segment count | `result:segment_count()` |
| Get segment | `result:segment(i)` |
| All segments | `result:segments()` |
| Path length | `result:path_length(mat_id)` |
| **2D Slicing** | |
| Axis-aligned view | `alea.slice_view_axis(axis, val, ...)` |
| Arbitrary view | `alea.slice_view(...)` |
| Grid query | `sys:find_cells_grid(view, nu, nv, depth)` |
| Grid overlaps | `sys:check_grid_overlaps(...)` |
| Get curves | `sys:get_slice_curves(view)` |
| Curve count | `curves:count()` |
| Get curve | `curves:get(i)` |
| Curve bounds | `curves:bounds()` |
| Label positions | `alea.find_label_positions(ids, w, h, min)` |
| Surface labels | `alea.find_surface_label_positions(...)` |
| **3D Rendering** | |
| Default config | `alea.render_config()` |
| Render scene | `sys:render(config)` |
| Camera setup | `sys:render_camera_setup(config)` |
| Write image | `fb:write(filename)` |
| Edge darken | `fb:edge_darken()` |
| Write aux maps | `fb:write_aux(base)` |
| Load colors | `alea.render_load_colors(file, cfg)` |
| **Mesh Export** | |
| Sample mesh | `sys:mesh_sample(config)` |
| One-shot export | `sys:mesh_export(config, filename)` |
| Export result | `mesh:export(format, filename)` |
| Mesh info | `mesh:info()` |
| Material IDs | `mesh:material_ids()` |
| Cell IDs | `mesh:cell_ids()` |
| **Error/Logging** | |
| Error message | `alea.error()` |
| Error code | `alea.error_code()` |
| Clear error | `alea.error_clear()` |
| Set log level | `alea.log_set_level(level)` |
| Get log level | `alea.log_get_level()` |
| Debug trace | `alea.set_debug_trace(bool)` |

## Next Steps

- See the example scripts in `examples/lua/` for complete working programs
- Read [Concepts](CONCEPTS.md) for domain concepts (sense, universes, lattices)
- Read [API Reference](API.md) for the complete C API (Lua bindings mirror this closely)
