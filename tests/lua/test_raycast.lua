-- test_raycast.lua: Raycast module tests

-- Build a sphere geometry
local sys = alea.create()
local s1 = sys:sphere(1, 0, 0, 0, 5)
local m1 = sys:material(1)
sys:cell{id = 1, region = sys:inside(s1),  material = m1, density = 10.0}
sys:cell{id = 2, region = sys:outside(s1), material = -1, density = 0.0}
sys:build_universe_index()

-- Basic raycast through sphere
local result = sys:raycast(-10, 0, 0, 1, 0, 0, 100)
assert(result:segment_count() >= 2, "should have at least 2 segments")

-- Check segments
local seg1 = result:segment(1)
assert(seg1.cell_id >= 0, "first segment should have a valid cell")

-- Get all segments
local segs = result:segments()
assert(#segs == result:segment_count(), "segments() count should match segment_count()")
for _, seg in ipairs(segs) do
    assert(seg.t_enter ~= nil, "segment should have t_enter")
    assert(seg.t_exit ~= nil, "segment should have t_exit")
    assert(seg.t_exit >= seg.t_enter, "t_exit >= t_enter")
end

-- Find the material-1 segment (through the sphere)
local mat1_found = false
for _, seg in ipairs(segs) do
    if seg.material_id == 1 then
        mat1_found = true
        -- Sphere radius is 5, ray passes through center
        -- Expected path length ~10 (diameter)
        local len = seg.t_exit - seg.t_enter
        assert(math.abs(len - 10.0) < 0.1, "path through sphere should be ~10")
    end
end
assert(mat1_found, "should find material 1 segment")

-- Path length
local pl = result:path_length(1)
assert(math.abs(pl - 10.0) < 0.1, "path_length(mat=1) should be ~10, got " .. pl)

local total = result:path_length(-1)
assert(total > pl, "total path should be > material path")

-- Cell-aware raycast
local result2 = sys:raycast_cell_aware(-10, 0, 0, 1, 0, 0, 100)
assert(result2:segment_count() >= 2, "cell-aware should also find segments")
local pl2 = result2:path_length(1)
assert(math.abs(pl2 - 10.0) < 0.1, "cell-aware path_length should match")

-- ray_first_cell
local fc, ft = sys:ray_first_cell(-10, 0, 0, 1, 0, 0, 100)
assert(fc ~= nil, "should find a first cell")
assert(ft ~= nil, "should have a distance")

-- Ray that misses everything meaningful
local fc2 = sys:ray_first_cell(100, 100, 100, 0, 0, 1, 1)
-- may or may not find void cell

sys:destroy()
print("test_raycast: OK")
