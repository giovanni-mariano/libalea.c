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

-- Selected cone sheets preserve the existing constructor form while accepting
-- -1 or +1 as an optional final argument.
local cone = sys:cone_z(6, 0, 0, 0, 1, 1)
assert(cone >= 0, "selected cone sheet created")
local bad_sheet = pcall(function()
    sys:cone_z(7, 0, 0, 0, 1, 2)
end)
assert(not bad_sheet, "invalid cone sheet rejected")

-- Boundary conditions use surface IDs rather than storage indices.
assert(sys:surface_get_boundary(1) == "transmissive", "default boundary")
sys:surface_set_boundary(1, "reflective")
assert(sys:surface_get_boundary(1) == "reflective", "updated boundary")
local missing_boundary = pcall(function()
    sys:surface_set_boundary(9999, "vacuum")
end)
assert(not missing_boundary, "missing surface ID rejected")
sys:surface_set_periodic_pair(3, 4)
assert(sys:surface_get_boundary(3) == "periodic" and
       sys:surface_get_boundary(4) == "periodic", "periodic pair")
assert(not pcall(function() sys:surface_set_periodic_pair(1, 3) end),
       "nonplane periodic pair rejected")

-- Standalone primitives use the same checked evaluator as the public C and
-- Python APIs.
local sphere_value = alea.primitive_evaluate(
    alea.PRIMITIVE_SPHERE, {0, 0, 0, 2}, {1, 0, 0})
assert(math.abs(sphere_value + 3) < 1e-12, "sphere implicit value")
local on_sphere = alea.primitive_evaluate(
    alea.PRIMITIVE_SPHERE, {0, 0, 0, 2}, {2, 0, 0})
assert(math.abs(on_sphere) < 1e-12, "sphere boundary value")
local bad_parameters = pcall(function()
    alea.primitive_evaluate(alea.PRIMITIVE_SPHERE, {0, 0, 2}, {0, 0, 0})
end)
assert(not bad_parameters, "wrong primitive parameter count rejected")

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

-- Add material
local mat = sys:material(1)  -- register material with MCNP ID 1, returns index
assert(mat >= 0, "material returns valid index")

-- Create cell
local idx = sys:cell{id=1, region=region, material=mat, density=10.0}
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

-- Programmatic transforms and cell editing
assert(sys:add_transform(10, {1, 2, 3}) == 10, "named transform added")
local inline_transform = sys:add_inline_transform(
    {0, 0, 1}, {cell_id = 1, role = "fill"})
assert(inline_transform > 0, "inline transform assigned an ID")
sys:set_fill(idx, 2, 10)
sys:set_comment(idx, "cell comment")
sys:set_inline_comment(idx, "inline comment")
sys:cell_set_density(idx, -9.5)
sys:cell_set_temperature(idx, 600)
sys:cell_set_universe(idx, 3)
info = sys:cell_info(idx)
assert(info.fill_universe == 2 and info.fill_transform == 10, "fill updated")
assert(info.comments == "cell comment", "cell comment updated")
assert(info.inline_comment == "inline comment", "inline comment updated")
assert(info.density == 9.5 and info.is_mass_density, "cell density updated")
assert(info.temperature == 600, "cell temperature updated")
assert(info.universe_id == 3, "cell universe updated")
sys:cell_clear_temperature(idx)
assert(sys:cell_info(idx).temperature == nil, "cell temperature cleared")

local bad_transform = pcall(function()
    sys:add_transform(11, {1, 2})
end)
assert(not bad_transform, "short transform rejected")

sys:cell_remove(idx)
assert(sys:cell_count() == 0, "cell removed")

print("test_geometry: OK")
