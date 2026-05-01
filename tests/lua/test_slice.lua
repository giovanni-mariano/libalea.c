-- test_slice.lua: Slice module tests

-- Build a sphere geometry
local sys = alea.create()
local s1 = sys:sphere(1, 0, 0, 0, 5)
local m1 = sys:material(1)
sys:cell{id = 1, region = sys:inside(s1),  material = m1, density = 10.0}
sys:cell{id = 2, region = sys:outside(s1), material = -1, density = 0.0}
sys:build_universe_index()

-- Slice view creation (axis-aligned)
local view = alea.slice_view_axis(2, 0.0, -10, 10, -10, 10)
assert(type(view) == "table", "view should be a table")
assert(view.u_min == -10, "u_min should be -10")
assert(view.u_max == 10, "u_max should be 10")
assert(view.v_min == -10, "v_min should be -10")
assert(view.v_max == 10, "v_max should be 10")

-- Arbitrary view
local view2 = alea.slice_view(0, 0, 0, 0, 0, 1, 0, 1, 0, -10, 10, -10, 10)
assert(type(view2) == "table", "arbitrary view should be a table")
assert(view2.u_min == -10, "u_min should be -10")

-- Grid query
local grid = sys:find_cells_grid(view, 8, 8)
assert(type(grid) == "table", "grid should be a table")
assert(#grid.cell_ids == 64, "should have 8*8=64 cell IDs")
assert(#grid.material_ids == 64, "should have 64 material IDs")
assert(#grid.errors == 64, "should have 64 error codes")

-- Center point should be inside sphere (cell 1)
-- Grid is 8x8 over -10..10, center is at index ~(4*8+4) = 36
-- Due to pixel sampling this might vary, just check we have both cells
local has_cell1 = false
local has_cell2 = false
for _, cid in ipairs(grid.cell_ids) do
    if cid == 1 then has_cell1 = true end
    if cid == 2 then has_cell2 = true end
end
assert(has_cell1, "should find cell 1 in grid")
assert(has_cell2, "should find cell 2 in grid")

-- Grid with depth parameter
local grid2 = sys:find_cells_grid(view, 4, 4, -1)
assert(#grid2.cell_ids == 16, "should have 4*4=16 cell IDs")

-- Get slice curves
local curves = sys:get_slice_curves(view)
assert(curves:count() >= 1, "sphere should produce at least 1 curve")

-- Get first curve
local c = curves:get(1)
assert(c.type == "circle", "sphere slice should be a circle, got " .. c.type)
assert(c.surface_id == 1, "curve surface_id should be 1")
assert(c.data ~= nil, "circle should have data")
assert(math.abs(c.data.radius - 5.0) < 0.1, "circle radius should be ~5")

-- Bounds
local u0, u1, v0, v1 = curves:bounds()
assert(type(u0) == "number", "bounds should return numbers")

-- Label positions
local labels = alea.find_label_positions(grid.cell_ids, 8, 8, 1)
assert(type(labels) == "table", "labels should be a table")
assert(#labels >= 1, "should find at least 1 label")
for _, lbl in ipairs(labels) do
    assert(lbl.id ~= nil, "label should have id")
    assert(lbl.px ~= nil, "label should have px")
    assert(lbl.py ~= nil, "label should have py")
end

-- Surface label positions
local slabels = alea.find_surface_label_positions(curves, -10, 10, -10, 10, 8, 8, 0)
assert(type(slabels) == "table", "surface labels should be a table")

-- Boundary-aware surface labels should follow the chosen contour ID grid.
local sys2 = alea.create()
local bs1 = sys2:sphere(1, 0, 0, 0, 3)
local bs2 = sys2:sphere(2, 0, 0, 0, 6)
local bm1 = sys2:material(1)
local bm2 = sys2:material(2)
sys2:cell{id = 1, region = sys2:inside(bs1), material = bm1, density = 1.0}
sys2:cell{id = 2, region = sys2:outside(bs1) * sys2:inside(bs2), material = bm1, density = 1.0}
sys2:cell{id = 3, region = sys2:outside(bs2), material = bm2, density = 1.0}
sys2:build_universe_index()

local bview = alea.slice_view_axis(2, 0.0, -8, 8, -8, 8)
local bgrid = sys2:find_cells_grid(bview, 64, 64)
local bcurves = sys2:get_slice_curves(bview)

local cell_surface_labels = alea.find_surface_label_positions(
    bcurves, -8, 8, -8, 8, 64, 64, 2, bgrid.cell_ids)
local material_surface_labels = alea.find_surface_label_positions(
    bcurves, -8, 8, -8, 8, 64, 64, 2, bgrid.material_ids)

local cell_has_s1 = false
local material_has_s1 = false
local material_has_s2 = false
for _, lbl in ipairs(cell_surface_labels) do
    if lbl.id == 1 then cell_has_s1 = true end
end
for _, lbl in ipairs(material_surface_labels) do
    if lbl.id == 1 then material_has_s1 = true end
    if lbl.id == 2 then material_has_s2 = true end
end
assert(cell_has_s1, "cell contours should label the inner cell boundary")
assert(not material_has_s1, "material contours should not label same-material cell boundary")
assert(material_has_s2, "material contours should label material boundary")

sys2:destroy()

sys:destroy()
print("test_slice: OK")
