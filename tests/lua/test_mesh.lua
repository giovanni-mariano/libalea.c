-- test_mesh.lua: Mesh module tests

-- Build a sphere geometry
local sys = alea.create()
local s1 = sys:sphere(1, 0, 0, 0, 5)
local m1 = sys:material(1)
sys:cell{id = 1, region = sys:inside(s1),  material = m1, density = 10.0}
sys:cell{id = 2, region = sys:outside(s1), material = -1, density = 0.0}
sys:build_universe_index()

-- Mesh sample
local mesh = sys:mesh_sample{nx = 3, ny = 3, nz = 3}
local info = mesh:info()
assert(info.nx == 3, "nx should be 3")
assert(info.ny == 3, "ny should be 3")
assert(info.nz == 3, "nz should be 3")
assert(info.num_materials >= 1, "should have at least 1 material")
assert(type(info.unique_materials) == "table", "should have unique_materials")

-- Material IDs
local mids = mesh:material_ids()
assert(#mids == 27, "should have 3*3*3=27 material IDs")

-- Cell IDs
local cids = mesh:cell_ids()
assert(#cids == 27, "should have 3*3*3=27 cell IDs")

-- Export Gmsh
mesh:export(0, "/tmp/alea_test_mesh.msh")

-- Export VTK
mesh:export(1, "/tmp/alea_test_mesh.vtk")

-- One-shot export
sys:mesh_export({nx = 3, ny = 3, nz = 3, format = 0}, "/tmp/alea_test_mesh2.msh")
sys:mesh_export({nx = 3, ny = 3, nz = 3, format = 1}, "/tmp/alea_test_mesh2.vtk")

-- With explicit bounds
local mesh2 = sys:mesh_sample{
    nx = 2, ny = 2, nz = 2,
    x_min = -6, x_max = 6,
    y_min = -6, y_max = 6,
    z_min = -6, z_max = 6,
}
local info2 = mesh2:info()
assert(info2.nx == 2, "bounded mesh nx should be 2")

sys:destroy()
print("test_mesh: OK")
