-- test_spatial_mode.lua: spatial mode API and hier-mode regression tests

-- ============================================================================
-- set_spatial_mode / get_spatial_mode API
-- ============================================================================

local sys = alea.create()

-- default is "flat"
local mode = sys:get_spatial_mode()
assert(mode == "flat" or mode == "hier" or mode == "auto",
    "get_spatial_mode returns a valid mode string, got: " .. tostring(mode))

sys:set_spatial_mode("flat")
assert(sys:get_spatial_mode() == "flat", "set flat")

sys:set_spatial_mode("hier")
assert(sys:get_spatial_mode() == "hier", "set hier")

sys:set_spatial_mode("auto")
assert(sys:get_spatial_mode() == "auto", "set auto")

-- invalid mode must error
local ok, err = pcall(function() sys:set_spatial_mode("bogus") end)
assert(not ok, "set_spatial_mode('bogus') should raise an error")
assert(err:find("bogus"), "error message should mention the bad value")

sys:destroy()

-- ============================================================================
-- find_cells_grid in hier mode
-- ============================================================================

local function make_sphere_sys()
    local s = alea.create()
    local si = s:sphere(1, 0, 0, 0, 5)
    local m1 = s:material(1)
    s:cell{id=1, region=s:inside(si),  material=m1, density=2.7}
    s:cell{id=2, region=s:outside(si), material=-1, density=0.0}
    s:build_universe_index()
    return s
end

local view = alea.slice_view_axis(2, 0.0, -10, 10, -10, 10)

-- flat reference
local sys_flat = make_sphere_sys()
sys_flat:set_spatial_mode("flat")
local g_flat = sys_flat:find_cells_grid(view, 16, 16)

-- hier
local sys_hier = make_sphere_sys()
sys_hier:set_spatial_mode("hier")
local g_hier = sys_hier:find_cells_grid(view, 16, 16)

local flat_stats = sys_flat:query_acceleration_stats()
assert(flat_stats.resolved_mode == "flat", "flat stats resolved mode")
assert(flat_stats.built, "flat stats built")
assert(flat_stats.flat_instance_count >= 1, "flat stats instance count")

local hier_stats = sys_hier:query_acceleration_stats()
assert(hier_stats.resolved_mode == "hier", "hier stats resolved mode")
assert(hier_stats.built, "hier stats built")
assert(hier_stats.hier_blas_count >= 1, "hier stats BLAS count")
assert(hier_stats.hier_placement_count >= 1, "hier stats placement count")
assert(hier_stats.flat_instance_count == 0, "hier stats flat instance count is zero")

local hier_paths = sys_hier:volume_paths()
assert(type(hier_paths) == "table", "hier volume_paths returns table")
assert(#hier_paths >= 1, "hier volume_paths has at least one path")
assert(hier_paths[1].cell_id ~= nil, "hier path has cell_id")

assert(#g_flat.cell_ids == 16*16, "flat grid size")
assert(#g_hier.cell_ids == 16*16, "hier grid size")

-- material layouts must be identical
for i = 1, #g_flat.material_ids do
    assert(g_flat.material_ids[i] == g_hier.material_ids[i],
        "material mismatch at pixel " .. i ..
        " flat=" .. g_flat.material_ids[i] ..
        " hier=" .. g_hier.material_ids[i])
end

sys_flat:destroy()
sys_hier:destroy()

-- ============================================================================
-- get_slice_curves in hier mode — parity with flat
-- ============================================================================

local sys_fc = make_sphere_sys()
sys_fc:set_spatial_mode("flat")
local curves_flat = sys_fc:get_slice_curves(view)
assert(curves_flat:count() >= 1, "flat: at least one curve")

local sys_hc = make_sphere_sys()
sys_hc:set_spatial_mode("hier")
local curves_hier = sys_hc:get_slice_curves(view)
assert(curves_hier:count() >= 1, "hier: at least one curve")

assert(curves_flat:count() == curves_hier:count(),
    "flat and hier curve counts must match: flat=" .. curves_flat:count() ..
    " hier=" .. curves_hier:count())

-- types must match
for i = 1, curves_flat:count() do
    local cf = curves_flat:get(i)
    local ch = curves_hier:get(i)
    assert(cf.type == ch.type,
        "curve " .. i .. " type mismatch: flat=" .. cf.type .. " hier=" .. ch.type)
    if cf.type == "circle" then
        assert(math.abs(cf.data.radius - ch.data.radius) < 0.01,
            "circle radius mismatch")
    end
end

sys_fc:destroy()
sys_hc:destroy()

-- ============================================================================
-- find_cells_grid in hier mode on fill model
-- ============================================================================

local fill_model = alea.load_mcnp("tests/data/simple_fill.mcnp")
if fill_model then
    fill_model:set_spatial_mode("hier")
    fill_model:prepare_query_acceleration()

    local fview = alea.slice_view_axis(2, 0.0, -110, 110, -110, 110)
    local fg = fill_model:find_cells_grid(fview, 32, 32)
    assert(#fg.material_ids == 32*32, "fill model grid size")

    -- center pixel should be inside the fill sphere (mat 1)
    local cx = 16 * 32 + 16 + 1  -- (row 16, col 16), 1-indexed
    assert(fg.material_ids[cx] == 1,
        "center of fill model should be mat 1, got " .. fg.material_ids[cx])

    fill_model:destroy()
end

print("test_spatial_mode: OK")
