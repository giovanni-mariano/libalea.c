-- test_util_extended.lua: Extended utility function tests

-- Build geometry
local sys = alea.create()
local s1 = sys:sphere(1, 0, 0, 0, 5)
local s2 = sys:sphere(2, 0, 0, 0, 8)
local m1 = sys:material(1)
local m2 = sys:material(2)
sys:cell{id = 1, region = sys:inside(s1), material = m1, density = 10.0}
sys:cell{id = 2, region = sys:outside(s1) * sys:inside(s2), material = m2, density = 5.0}
sys:cell{id = 3, region = sys:outside(s2), material = -1, density = 0.0}
sys:build_universe_index()

-- simplify_all
local sys2 = sys:clone()
local stats = sys2:simplify_all()
assert(type(stats) == "table", "simplify_all should return stats")
assert(stats.nodes_before >= 0, "should have nodes_before")
assert(stats.nodes_after >= 0, "should have nodes_after")
sys2:destroy()

-- tighten_cell_bbox
local bbox = sys:tighten_cell_bbox(0, 1.0)
assert(type(bbox) == "table", "should return bbox table")
assert(bbox.min_x ~= nil, "should have min_x")
assert(bbox.max_x ~= nil, "should have max_x")
assert(bbox.max_x > bbox.min_x, "max_x > min_x")

-- extract_region
local region = sys:extract_region{
    min_x = -3, max_x = 3,
    min_y = -3, max_y = 3,
    min_z = -3, max_z = 3,
}
assert(region:cell_count() >= 1, "extracted region should have cells")
region:destroy()

-- cells_in_bbox
local cells = sys:cells_in_bbox{
    min_x = -10, max_x = 10,
    min_y = -10, max_y = 10,
    min_z = -10, max_z = 10,
}
assert(#cells >= 1, "should find cells in large bbox")

-- volume_paths / estimate_volumes preserve physical path identity
sys:prepare_query_acceleration()
local paths = sys:volume_paths()
assert(type(paths) == "table", "volume_paths should return table")
assert(#paths >= 1, "volume_paths should include flat instances")
assert(paths[1].path_id == 1, "Lua path_id is 1-based")
assert(paths[1].cell_id ~= nil, "volume path should include cell_id")

local path_volumes = sys:estimate_volumes(8)
assert(type(path_volumes) == "table", "estimate_volumes should return table")
assert(#path_volumes == #paths, "path volume count should match path count")
assert(path_volumes[1].volume ~= nil, "path volume record should include volume")

-- expand_macrobody_node / expand_all_macrobodies_node
-- These work on nodes; sphere is not a macrobody so result should be same node
local info = sys:cell_info(0)
local expanded = sys:expand_macrobody_node(info.root)
-- Just verify it doesn't crash

-- is_macrobody
assert(alea.is_macrobody(14) == true, "RCC (14) should be a macrobody")
assert(alea.is_macrobody(2) == false, "SPHERE (2) should not be a macrobody")

-- void_get and void_to_node
-- Use a simple box geometry for void generation (sphere geometry is too expensive)
local vsys = alea.create()
local box_s = vsys:box(0, -1, 1, -1, 1, -1, 1)
local vm1 = vsys:material(1)
vsys:cell{id = 1, region = vsys:inside(box_s), material = vm1, density = 1.0}
vsys:cell{id = 2, region = vsys:outside(box_s), material = -1, density = 0.0}
vsys:build_universe_index()

local vr = vsys:void_generate{
    min_x = -5, max_x = 5,
    min_y = -5, max_y = 5,
    min_z = -5, max_z = 5,
}
local vc = vr:count()
assert(vc >= 0, "void count should be >= 0")

if vc > 0 then
    -- void:get
    local vbox = vr:get(1)
    assert(type(vbox) == "table", "void:get should return table")
    assert(vbox.min_x ~= nil, "void box should have min_x")

    -- void:to_node
    local vnode = vr:to_node()
    -- may be nil if no voids
end
vsys:destroy()

-- create_mixture
local mid = sys:create_mixture({1, 2}, {0.5, 0.5}, 100)
assert(mid == 100, "create_mixture should return 100, got " .. tostring(mid))

sys:destroy()
print("test_util_extended: OK")
