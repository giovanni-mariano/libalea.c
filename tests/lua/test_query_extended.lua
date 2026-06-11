-- test_query_extended.lua: Extended query function tests

-- Build geometry with multiple cells/universes
local sys = alea.create()
local s1 = sys:sphere(1, 0, 0, 0, 5)
local s2 = sys:sphere(2, 0, 0, 0, 8)
local m1 = sys:material(1)
local m2 = sys:material(2)
sys:cell{id = 10, region = sys:inside(s1), material = m1, density = 10.0, universe = 0}
sys:cell{id = 20, region = sys:outside(s1) * sys:inside(s2), material = m2, density = 5.0, universe = 0}
sys:cell{id = 30, region = sys:outside(s2), material = -1, density = 0.0, universe = 0}
sys:build_universe_index()

-- find_cell_at
local cid, mat = sys:find_cell_at(0, 0, 0)
assert(cid == 10, "find_cell_at(0,0,0) should return cell 10, got " .. tostring(cid))
assert(mat == 1, "material should be 1, got " .. tostring(mat))

local cid2, mat2 = sys:find_cell_at(6, 0, 0)
assert(cid2 == 20, "find_cell_at(6,0,0) should return cell 20, got " .. tostring(cid2))
assert(mat2 == 2, "material should be 2, got " .. tostring(mat2))

-- cell_find_info
local info = sys:cell_find_info(10)
assert(info ~= nil, "cell_find_info(10) should succeed")
assert(info.cell_id == 10, "cell_id should be 10")
assert(info.material_id == 1, "material_id should be 1")
assert(info.root ~= nil, "should have root node")

local nil_info = sys:cell_find_info(999)
assert(nil_info == nil, "cell_find_info(999) should return nil")

-- cells_in_universe
local cells = sys:cells_in_universe(0)
assert(#cells == 3, "universe 0 should have 3 cells, got " .. #cells)

-- surface_info
local si = sys:surface_info(0)
assert(si.surface_id == 1, "first surface should have id 1")
assert(si.type > 0, "should have a valid type")
assert(si.pos_node ~= nil, "should have pos_node")
assert(si.neg_node ~= nil, "should have neg_node")

-- surface_find
local sf = sys:surface_find(1)
assert(sf == 0, "surface 1 should be at index 0")
local sf2 = sys:surface_find(999)
assert(sf2 == nil, "surface 999 should not be found")

-- universe_info
local ui = sys:universe_info(0)
assert(ui.universe_id == 0, "universe_id should be 0")
assert(ui.cell_count == 3, "should have 3 cells in universe 0")
assert(type(ui.bbox) == "table", "should have bbox")

-- universe_find
local uf = sys:universe_find(0)
assert(uf == 0, "universe 0 should be at index 0")
local uf2 = sys:universe_find(999)
assert(uf2 == nil, "universe 999 should not be found")

-- node_primitive_type
local cell_info = sys:cell_info(0)
local root = cell_info.root
local ptype = sys:node_primitive_type(root)
-- root is inside(s1), which is a primitive (sphere)
assert(ptype == "sphere", "node_primitive_type should be sphere, got " .. ptype)

-- node_primitive_id
local pid = sys:node_primitive_id(root)
assert(pid ~= nil, "primitive_id should not be nil")

-- node_primitive_data
local pdata = sys:node_primitive_data(root)
assert(pdata ~= nil, "primitive_data should not be nil")
assert(pdata.r ~= nil, "sphere data should have radius")
assert(math.abs(pdata.r - 5.0) < 0.01, "sphere radius should be 5")

-- cells_filling_universe (none in this model)
local filling = sys:cells_filling_universe(0)
assert(type(filling) == "table", "should return a table")

-- query acceleration stats (hierarchical spatial index)
local qstats = sys:query_acceleration_stats()
assert(type(qstats) == "table", "query_acceleration_stats should return a table")
assert(qstats.built == false, "query acceleration should not be built yet")
assert(type(qstats.hier_universe_count) == "number", "hier_universe_count should be a number")

-- stats
local st = sys:stats()
assert(type(st) == "table", "stats should be a table")
assert(st.unique_primitives >= 1, "should have primitives")

-- tree_print (just make sure it doesn't crash)
sys:tree_print(root)

-- set_debug_trace (just verify it doesn't error)
alea.set_debug_trace(false)

sys:destroy()
print("test_query_extended: OK")
