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

local sample_counts = mesh:sample_counts()
local tie_flags = mesh:tie_flags()
local estimated_errors = mesh:estimated_errors()
local refinement_flags = mesh:refinement_flags()
assert(#sample_counts == 27, "should have one sample count per voxel")
assert(#tie_flags == 27, "should have one tie flag per voxel")
assert(#estimated_errors == 27, "should have one estimated error per voxel")
assert(#refinement_flags == 27, "should have one refinement flag per voxel")
assert(sample_counts[1] == 8, "default mode should take 2^3 samples")

local function join_path(dir, name)
    local sep = package.config and package.config:sub(1, 1) or "/"
    if not dir or dir == "" then
        dir = "."
    end
    local last = dir:sub(-1)
    if last == "/" or last == "\\" then
        return dir .. name
    end
    return dir .. sep .. name
end

local tmpdir = os.getenv("TMPDIR") or os.getenv("TEMP") or os.getenv("TMP") or "."
local gmsh_file = join_path(tmpdir, "alea_test_mesh.msh")
local vtk_file = join_path(tmpdir, "alea_test_mesh.vtk")
local gmsh_file2 = join_path(tmpdir, "alea_test_mesh2.msh")
local vtk_file2 = join_path(tmpdir, "alea_test_mesh2.vtk")

-- Export Gmsh
mesh:export(0, gmsh_file)

-- Export VTK
mesh:export(1, vtk_file, {export_fields = 31, max_fraction_materials = 16})
local vtk = assert(io.open(vtk_file, "r"))
local vtk_text = vtk:read("*a")
vtk:close()
assert(vtk_text:find("sampled_fraction_material_", 1, true),
       "opt-in VTK export should contain material fraction arrays")

-- One-shot export
sys:mesh_export({nx = 3, ny = 3, nz = 3, format = 0}, gmsh_file2)
sys:mesh_export({nx = 3, ny = 3, nz = 3, format = 1}, vtk_file2)
os.remove(gmsh_file)
os.remove(vtk_file)
os.remove(gmsh_file2)
os.remove(vtk_file2)

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
    assert(entry.sampled_fraction == entry.fraction,
           "sampled_fraction compatibility field should match fraction")
end
assert(math.abs(by_mat[1] - 0.125) < 1e-12, "material 1 fraction should be 1/8")
assert(math.abs(by_mat[0] - 0.875) < 1e-12, "void fraction should be 7/8")

local adaptive = sys2:mesh_sample{
    nx = 1, ny = 1, nz = 1,
    x_min = 0, x_max = 1,
    y_min = 0, y_max = 1,
    z_min = 0, z_max = 1,
    sampling_mode = 4,
    subsamples_per_axis = 2,
    target_error = 0.1,
    max_refine_depth = 1,
    max_samples_per_voxel = 72,
    sampling_seed = 7,
}
assert(adaptive:info().sampling_mode == 4)
assert(adaptive:sample_counts()[1] == 72)
assert(adaptive:estimated_errors()[1] >= 0)

sys2:destroy()
sys:destroy()
print("test_mesh: OK")
