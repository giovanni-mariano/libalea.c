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

-- Mixed diagnostics: small inclusion missed by center sampling
local sys2 = alea.create()
local inc = sys2:sphere(1, 0.25, 0.25, 0.25, 0.20)
local mat = sys2:material(1)
sys2:cell{id = 1, region = sys2:inside(inc), material = mat, density = 1.0}
sys2:build_universe_index()

local mixed_mesh = sys2:mesh_sample{
    nx = 1, ny = 1, nz = 1,
    x_min = 0, x_max = 1,
    y_min = 0, y_max = 1,
    z_min = 0, z_max = 1,
    sampling_mode = 2,
    subsamples_per_axis = 2,
}
local mixed_info = mixed_mesh:info()
assert(mixed_info.mixed_count == 1, "small inclusion should flag one mixed voxel")
assert(mixed_info.fraction_count == 2, "one mixed voxel should have two fraction entries")

local fractions = mixed_mesh:material_fractions(1)
assert(#fractions == 2, "voxel should have two material fractions")
local by_mat = {}
for _, entry in ipairs(fractions) do
    by_mat[entry.material_id] = entry.fraction
end
assert(math.abs(by_mat[1] - 0.125) < 1e-12, "material 1 fraction should be 1/8")
assert(math.abs(by_mat[0] - 0.875) < 1e-12, "void fraction should be 7/8")

sys2:destroy()
sys:destroy()
print("test_mesh: OK")
