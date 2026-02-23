-- 12_slice.lua: 2D cross-section slicing and curve extraction
--
-- Standalone: yes (no data files needed)
-- Usage: bin/alea examples/lua/12_slice.lua
--
-- Geometry: concentric spheres (fuel / moderator / void)
--   Cell 1: fuel      (mat 1)  r < 5
--   Cell 2: moderator (mat 2)  5 <= r < 10
--   Cell 3: void      (mat 0)  r >= 10

print("=== 2D Slicing ===\n")

local sys = alea.create()

-- Build concentric spheres
local s1 = sys:sphere(1, 0, 0, 0, 5)
local s2 = sys:sphere(2, 0, 0, 0, 10)

sys:cell{id = 1, region = sys:inside(s1), material = 1, density = 10.0}
sys:cell{id = 2, region = sys:outside(s1) * sys:inside(s2), material = 2, density = 1.0}
sys:cell{id = 3, region = sys:outside(s2), material = 0, density = 0.0}

sys:build_universe_index()

-- Create an axis-aligned Z slice at z=0
print("--- Axis-Aligned Slice (Z=0) ---")
local view = alea.slice_view_axis(2, 0.0, -15, 15, -15, 15)
print(string.format("  Viewport: u=[%.1f, %.1f] v=[%.1f, %.1f]",
    view.u_min, view.u_max, view.v_min, view.v_max))

-- Sample the slice on a grid
local nu, nv = 16, 16
local grid = sys:find_cells_grid(view, nu, nv)
print(string.format("  Grid: %dx%d = %d points", nu, nv, #grid.cell_ids))

-- Count cells and materials found
local cell_counts = {}
local mat_counts = {}
for i = 1, #grid.cell_ids do
    local cid = grid.cell_ids[i]
    local mid = grid.material_ids[i]
    cell_counts[cid] = (cell_counts[cid] or 0) + 1
    mat_counts[mid] = (mat_counts[mid] or 0) + 1
end

print("  Cells found:")
local ckeys = {}
for k in pairs(cell_counts) do ckeys[#ckeys + 1] = k end
table.sort(ckeys)
for _, cid in ipairs(ckeys) do
    print(string.format("    cell %d: %d samples", cid, cell_counts[cid]))
end

-- Count errors
local n_ok, n_overlap, n_undef = 0, 0, 0
for _, e in ipairs(grid.errors) do
    if e == 0 then n_ok = n_ok + 1
    elseif e == 1 then n_overlap = n_overlap + 1
    else n_undef = n_undef + 1
    end
end
print(string.format("  Errors: %d ok, %d overlap, %d undefined", n_ok, n_overlap, n_undef))

-- Analytical curves
print("\n--- Slice Curves ---")
local curves = sys:get_slice_curves(view)
print(string.format("  Curve count: %d", curves:count()))

for i = 1, curves:count() do
    local c = curves:get(i)
    if c.type == "circle" then
        print(string.format("    [%d] %s (surface %d): center=(%.2f, %.2f) radius=%.2f",
            i, c.type, c.surface_id, c.data.cx, c.data.cy, c.data.radius))
    elseif c.type == "line" then
        print(string.format("    [%d] %s (surface %d): p=(%.2f, %.2f) d=(%.2f, %.2f)",
            i, c.type, c.surface_id, c.data.px, c.data.py, c.data.dx, c.data.dy))
    else
        print(string.format("    [%d] %s (surface %d)", i, c.type, c.surface_id))
    end
end

local u0, u1, v0, v1 = curves:bounds()
print(string.format("  Curve bounds: u=[%.2f, %.2f] v=[%.2f, %.2f]", u0, u1, v0, v1))

-- Label positions for cell IDs
print("\n--- Label Positions ---")
local labels = alea.find_label_positions(grid.cell_ids, nu, nv, 1)
print(string.format("  Found %d labels", #labels))
for _, lbl in ipairs(labels) do
    print(string.format("    ID %d at pixel (%d, %d), size=%d px",
        lbl.id, lbl.px, lbl.py, lbl.pixel_count))
end

-- Surface label positions
local slabels = alea.find_surface_label_positions(curves,
    view.u_min, view.u_max, view.v_min, view.v_max, nu, nv, 0)
print(string.format("  Surface labels: %d", #slabels))

-- Arbitrary-orientation slice (diagonal through model)
print("\n--- Arbitrary Slice ---")
local view2 = alea.slice_view(0, 0, 0,  0, 0, 1,  0, 1, 0,  -15, 15, -15, 15)
local grid2 = sys:find_cells_grid(view2, 8, 8)
print(string.format("  Grid: 8x8 = %d points", #grid2.cell_ids))

-- Grid with depth control
print("\n--- Depth Control ---")
local grid_root = sys:find_cells_grid(view, 8, 8, 0)
local grid_inner = sys:find_cells_grid(view, 8, 8, -1)
print(string.format("  Root level:     %d points", #grid_root.cell_ids))
print(string.format("  Innermost (-1): %d points", #grid_inner.cell_ids))

print("\n12_slice: OK")
