-- 11_raycast.lua: Ray tracing through CSG geometry
--
-- Standalone: yes (no data files needed)
-- Usage: bin/alea examples/lua/11_raycast.lua
--
-- Geometry: concentric spheres (fuel / moderator / void)
--   Cell 1: fuel      (mat 1)  r < 5
--   Cell 2: moderator (mat 2)  5 <= r < 10
--   Cell 3: void      (mat 0)  r >= 10

print("=== Ray Tracing ===\n")

local sys = alea.create()

-- Build concentric spheres
local s1 = sys:sphere(1, 0, 0, 0, 5)
local s2 = sys:sphere(2, 0, 0, 0, 10)

sys:cell{id = 1, region = sys:inside(s1), material = 1, density = 10.0}
sys:cell{id = 2, region = sys:outside(s1) * sys:inside(s2), material = 2, density = 1.0}
sys:cell{id = 3, region = sys:outside(s2), material = 0, density = 0.0}

sys:build_universe_index()

-- Cast a ray along the X axis through the center
print("--- Basic Raycast ---")
local result = sys:raycast(-20, 0, 0, 1, 0, 0, 100)

print(string.format("  Segments: %d", result:segment_count()))
print(string.format("  %-12s  %-12s  %-8s  %-6s  %-8s", "t_enter", "t_exit", "cell_id", "mat", "density"))
print(string.rep("-", 54))

for i = 1, result:segment_count() do
    local seg = result:segment(i)
    print(string.format("  %12.4f  %12.4f  %-8d  %-6d  %8.2f",
        seg.t_enter, seg.t_exit, seg.cell_id, seg.material_id, seg.density))
end

-- Path lengths by material
print("\n--- Path Lengths ---")
local total = result:path_length(-1)
local fuel  = result:path_length(1)
local mod   = result:path_length(2)
print(string.format("  Total path:     %8.4f", total))
print(string.format("  In fuel (mat 1):%8.4f", fuel))
print(string.format("  In mod  (mat 2):%8.4f", mod))

-- Get all segments at once
print("\n--- All Segments (batch) ---")
local segs = result:segments()
print(string.format("  Retrieved %d segments at once", #segs))

-- Cast rays from multiple directions
print("\n--- Multi-Direction Rays ---")
local directions = {
    { 1,  0,  0, "along +X"},
    { 0,  1,  0, "along +Y"},
    { 0,  0,  1, "along +Z"},
    { 1,  1,  1, "along diagonal"},
}

for _, dir in ipairs(directions) do
    local dx, dy, dz = dir[1], dir[2], dir[3]
    -- Normalize direction
    local len = math.sqrt(dx*dx + dy*dy + dz*dz)
    dx, dy, dz = dx/len, dy/len, dz/len

    local r = sys:raycast(-20*dx, -20*dy, -20*dz, dx, dy, dz, 100)
    print(string.format("  %-16s  segments=%d  fuel_path=%.4f  total_path=%.4f",
        dir[4], r:segment_count(), r:path_length(1), r:path_length(-1)))
end

-- Cell-aware raycast (uses per-cell surface indices)
print("\n--- Cell-Aware Raycast ---")
local result2 = sys:raycast_cell_aware(-20, 0, 0, 1, 0, 0, 100)
print(string.format("  Segments: %d", result2:segment_count()))
print(string.format("  Fuel path: %.4f (should match basic raycast)", result2:path_length(1)))

-- Find first cell hit by a ray
print("\n--- First Cell Hit ---")
local origins = {
    {-20, 0, 0, "from -X"},
    { 20, 0, 0, "from +X"},
    {  0, 0, 20, "from +Z"},
}

for _, orig in ipairs(origins) do
    local cell_id, t = sys:ray_first_cell(orig[1], orig[2], orig[3], -orig[1], -orig[2], -orig[3], 100)
    if cell_id then
        print(string.format("  %-12s  first cell=%d at t=%.4f", orig[4], cell_id, t))
    else
        print(string.format("  %-12s  no cell hit", orig[4]))
    end
end

print("\n11_raycast: OK")
