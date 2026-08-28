-- 14_mesh.lua: Mesh sampling and export (Gmsh, VTK)
--
-- Standalone: yes (no data files needed)
-- Usage: bin/alea examples/lua/14_mesh.lua
--
-- Geometry: concentric spheres (fuel / moderator / void)
--   Cell 1: fuel      (mat 1)  r < 5
--   Cell 2: moderator (mat 2)  5 <= r < 10
--   Cell 3: void      (mat 0)  r >= 10
--
-- Writes mesh files to /tmp/alea_example_mesh_*

print("=== Mesh Export ===\n")

local sys = alea.create()

-- Build concentric spheres
local s1 = sys:sphere(1, 0, 0, 0, 5)
local s2 = sys:sphere(2, 0, 0, 0, 10)

sys:cell{id = 1, region = sys:inside(s1), material = 1, density = 10.0}
sys:cell{id = 2, region = sys:outside(s1) * sys:inside(s2), material = 2, density = 1.0}
sys:cell{id = 3, region = sys:outside(s2), material = 0, density = 0.0}

sys:build_universe_index()

-- Step-by-step mesh sampling with inspection
print("--- Mesh Sample (10x10x10) ---")
local mesh = sys:mesh_sample{nx = 10, ny = 10, nz = 10}

local info = mesh:info()
print(string.format("  Grid: %dx%dx%d = %d voxels", info.nx, info.ny, info.nz, info.nx * info.ny * info.nz))
print(string.format("  Materials found: %d", info.num_materials))
print("  Unique materials:")
for _, mid in ipairs(info.unique_materials) do
    print(string.format("    material %d", mid))
end

-- Inspect material distribution
print("\n--- Material Distribution ---")
local mids = mesh:material_ids()
local cids = mesh:cell_ids()
print(string.format("  Material IDs: %d entries", #mids))
print(string.format("  Cell IDs:     %d entries", #cids))

local mat_dist = {}
for _, mid in ipairs(mids) do
    mat_dist[mid] = (mat_dist[mid] or 0) + 1
end

local mkeys = {}
for k in pairs(mat_dist) do mkeys[#mkeys + 1] = k end
table.sort(mkeys)

print(string.format("  %-10s  %-8s  %-8s", "Material", "Voxels", "Fraction"))
print(string.rep("-", 30))
local total = #mids
for _, mid in ipairs(mkeys) do
    local n = mat_dist[mid]
    print(string.format("  %-10d  %-8d  %7.1f%%", mid, n, 100 * n / total))
end

-- Export to Gmsh format
print("\n--- Export Gmsh ---")
mesh:export(0, "/tmp/alea_example_mesh.msh")
print("  Written: /tmp/alea_example_mesh.msh")

-- Export to VTK format
print("\n--- Export VTK ---")
mesh:export(1, "/tmp/alea_example_mesh.vtk")
print("  Written: /tmp/alea_example_mesh.vtk")

-- One-shot export (sample + write in one call)
print("\n--- One-Shot Export ---")
sys:mesh_export({nx = 8, ny = 8, nz = 8, format = 0}, "/tmp/alea_example_mesh_quick.msh")
sys:mesh_export({nx = 8, ny = 8, nz = 8, format = 1}, "/tmp/alea_example_mesh_quick.vtk")
print("  Written: /tmp/alea_example_mesh_quick.msh")
print("  Written: /tmp/alea_example_mesh_quick.vtk")

-- Mesh with explicit bounds (smaller region around the fuel)
print("\n--- Bounded Mesh (fuel region only) ---")
local mesh2 = sys:mesh_sample{
    nx = 8, ny = 8, nz = 8,
    x_min = -6, x_max = 6,
    y_min = -6, y_max = 6,
    z_min = -6, z_max = 6,
}
local info2 = mesh2:info()
print(string.format("  Grid: %dx%dx%d", info2.nx, info2.ny, info2.nz))
print(string.format("  Materials: %d", info2.num_materials))

-- Separate adaptive octree (unstructured export, not boundary conforming)
print("\n--- Adaptive Grid ---")
local adaptive = sys:adaptive_grid_sample{
    nx = 2, ny = 2, nz = 2,
    max_grid_depth = 2,
    max_cells = 10000,
}
local adaptive_info = adaptive:info()
print(string.format("  Cells: %d total, %d leaves",
                    adaptive_info.cell_count, adaptive_info.leaf_count))
adaptive:export(1, "/tmp/alea_example_mesh_adaptive.vtk")

print("\n14_mesh: OK")
