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
local idx = sys:find_cell(x, y, z)
if idx then
    local info = sys:cell_info(idx)
    print("Cell ID: " .. info.cell_id)
    print("Material: " .. info.material_id)
else
    print("No cell found (void or undefined)")
end
```

### Material at a point

```lua
local mat = sys:material_at(0, 0, 0)
```

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
sys:build_universe_index()
sys:build_spatial_index()
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
local volumes = sys:estimate_volumes(n_rays, cx, cy, cz, r)

for i, v in ipairs(volumes) do
    if v.volume > 0 then
        local info = sys:cell_info(i - 1)
        print(string.format("Cell %d: volume %.2f (error %.1f%%)",
            info.cell_id, v.volume, v.rel_error * 100))
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

## Quick Reference

| Operation | Lua |
|-----------|-----|
| Create system | `alea.create()` |
| Load MCNP | `alea.load_mcnp(file)` |
| Load OpenMC | `alea.load_openmc(file)` |
| Load from string | `alea.load_mcnp_string(str)` |
| Export MCNP | `sys:export_mcnp(file)` |
| Export OpenMC | `sys:export_openmc(file)` |
| Build index | `sys:build_universe_index()` |
| Find cell | `sys:find_cell(x, y, z)` |
| Material at point | `sys:material_at(x, y, z)` |
| Cell info | `sys:cell_info(idx)` |
| Intersection | `a * b` |
| Union | `a + b` |
| Difference | `a - b` |
| Complement | `~a` |
| Inside surface | `sys:inside(surf_idx)` |
| Outside surface | `sys:outside(surf_idx)` |
| Flatten | `sys:flatten(universe)` |
| Extract universe | `sys:extract_universe(uid)` |
| Merge | `sys:merge(other, offset)` |
| Clone | `sys:clone()` |
| Validate | `sys:validate()` |

## Next Steps

- See the example scripts in `examples/lua/` for complete working programs
- Read [Concepts](CONCEPTS.md) for domain concepts (sense, universes, lattices)
- Read [API Reference](API.md) for the complete C API (Lua bindings mirror this closely)
