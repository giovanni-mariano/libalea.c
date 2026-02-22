-- test_geometry.lua: Surface creation and cell construction

local sys = alea.create()

-- Create surfaces
local s1 = sys:sphere(1, 0, 0, 0, 10)
assert(type(s1) == "number", "sphere returns index")
assert(s1 >= 0, "valid surface index")

local s2 = sys:cylinder_z(2, 0, 0, 3)
assert(s2 >= 0, "valid cylinder index")

local p1 = sys:plane(3, 0, 0, 1, 5)   -- z = -5
local p2 = sys:plane(4, 0, 0, 1, -5)  -- z = 5

-- Halfspaces
local inside_sphere = sys:inside(s1)
assert(inside_sphere, "inside returns a Node")
local outside_cyl = sys:outside(s2)
assert(outside_cyl, "outside returns a Node")

-- Boolean operations via operators
local region = inside_sphere - sys:inside(s2)  -- sphere minus cylinder
assert(region, "difference returns a Node")

-- Verify tostring
local str = tostring(region)
assert(str:find("Node"), "Node tostring contains 'Node'")

-- Create cell
local idx = sys:cell{id=1, region=region, material=1, density=10.0}
assert(idx >= 0, "cell returns valid index")
assert(sys:cell_count() == 1, "one cell created")

-- Cell info
local info = sys:cell_info(idx)
assert(info.cell_id == 1, "cell_id matches")
assert(info.material_id == 1, "material matches")
assert(info.density == 10.0, "density matches")

-- Multiple surfaces
local b = sys:box(5, -20, 20, -20, 20, -20, 20)
assert(b >= 0, "box surface created")

print("test_geometry: OK")
