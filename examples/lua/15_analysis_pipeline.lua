-- 15_analysis_pipeline.lua: End-to-end analysis combining all modules
--
-- Standalone: yes (no data files needed)
-- Usage: bin/alea examples/lua/15_analysis_pipeline.lua
--
-- Builds a shielded reactor pin model, then demonstrates a complete analysis
-- pipeline: inspection, ray tracing, 2D slicing, 3D rendering, and mesh export.
--
-- Geometry:
--   Cell 1: fuel pin    (mat 1)  r < 0.4 cm
--   Cell 2: gap         (mat 2)  0.4 <= r < 0.42 cm
--   Cell 3: clad        (mat 3)  0.42 <= r < 0.48 cm
--   Cell 4: moderator   (mat 4)  0.48 <= r < 0.63 cm (pitch/2)
--   Cell 5: void        (mat 0)  r >= 0.63 cm
--
-- Output files: /tmp/alea_pipeline_*

print("=== Analysis Pipeline ===\n")

-- ---------------------------------------------------------------
-- 1. Build the geometry
-- ---------------------------------------------------------------
print("--- 1. Build Geometry ---")
local sys = alea.create()

local r_fuel = 0.4
local r_gap  = 0.42
local r_clad = 0.48
local pitch  = 1.26

local s1 = sys:sphere(1, 0, 0, 0, r_fuel)
local s2 = sys:sphere(2, 0, 0, 0, r_gap)
local s3 = sys:sphere(3, 0, 0, 0, r_clad)
local s4 = sys:sphere(4, 0, 0, 0, pitch / 2)

sys:cell{id = 1, region = sys:inside(s1),                            material = 1, density = 10.97}
sys:cell{id = 2, region = sys:outside(s1) * sys:inside(s2),          material = 2, density = 0.001}
sys:cell{id = 3, region = sys:outside(s2) * sys:inside(s3),          material = 3, density = 6.56}
sys:cell{id = 4, region = sys:outside(s3) * sys:inside(s4),          material = 4, density = 1.0}
sys:cell{id = 5, region = sys:outside(s4),                           material = 0, density = 0.0}

sys:build_universe_index()

local nc = sys:cell_count()
local ns = sys:surface_count()
print(string.format("  Built: %d cells, %d surfaces", nc, ns))

-- ---------------------------------------------------------------
-- 2. Model inspection
-- ---------------------------------------------------------------
print("\n--- 2. Inspect Model ---")

-- Cell summary
print(string.format("  %-6s  %-6s  %-8s  %-10s", "Cell", "Mat", "Density", "Universe"))
print(string.rep("-", 36))
for i = 0, nc - 1 do
    local ci = sys:cell_info(i)
    print(string.format("  %-6d  %-6d  %-8.3f  %-10d",
        ci.cell_id, ci.material_id, ci.density, ci.universe_id))
end

-- Stats
local st = sys:stats()
print(string.format("\n  Unique primitives: %d", st.unique_primitives))

-- Universe listing
local ui = sys:universe_info(0)
print(string.format("  Universe 0: %d cells", ui.cell_count))

-- ---------------------------------------------------------------
-- 3. Ray tracing
-- ---------------------------------------------------------------
print("\n--- 3. Ray Tracing ---")

-- Shoot a ray along X through the center
local result = sys:raycast(-2, 0, 0, 1, 0, 0, 10)
print(string.format("  Ray along X: %d segments", result:segment_count()))

print(string.format("  %-10s  %-10s  %-6s  %-6s  %-8s  %-10s",
    "t_enter", "t_exit", "cell", "mat", "density", "thickness"))
print(string.rep("-", 58))

for i = 1, result:segment_count() do
    local seg = result:segment(i)
    local thick = seg.t_exit - seg.t_enter
    print(string.format("  %10.6f  %10.6f  %-6d  %-6d  %8.3f  %10.6f",
        seg.t_enter, seg.t_exit, seg.cell_id, seg.material_id, seg.density, thick))
end

-- Material path lengths
print("\n  Path lengths:")
print(string.format("    Fuel (mat 1): %.6f cm", result:path_length(1)))
print(string.format("    Gap  (mat 2): %.6f cm", result:path_length(2)))
print(string.format("    Clad (mat 3): %.6f cm", result:path_length(3)))
print(string.format("    Mod  (mat 4): %.6f cm", result:path_length(4)))
print(string.format("    Total:        %.6f cm", result:path_length(-1)))

-- ---------------------------------------------------------------
-- 4. 2D slice
-- ---------------------------------------------------------------
print("\n--- 4. 2D Slice (Z=0 plane) ---")

local extent = pitch
local view = alea.slice_view_axis(2, 0.0, -extent, extent, -extent, extent)
local nu, nv = 32, 32
local grid = sys:find_cells_grid(view, nu, nv)

-- Material map as ASCII art
print(string.format("  %dx%d material map:", nu, nv))
local mat_chars = {[0] = ".", [1] = "F", [2] = "g", [3] = "C", [4] = "M"}
local row = {}
for j = 0, nv - 1 do
    row = {}
    for i = 0, nu - 1 do
        local mid = grid.material_ids[j * nu + i + 1]
        row[#row + 1] = mat_chars[mid] or "?"
    end
    print("  " .. table.concat(row))
end
print("  Legend: F=fuel g=gap C=clad M=moderator .=void")

-- Analytical curves
local curves = sys:get_slice_curves(view)
print(string.format("\n  Slice curves: %d", curves:count()))
for i = 1, curves:count() do
    local c = curves:get(i)
    if c.type == "circle" then
        print(string.format("    surface %d: circle r=%.4f", c.surface_id, c.data.radius))
    else
        print(string.format("    surface %d: %s", c.surface_id, c.type))
    end
end

-- ---------------------------------------------------------------
-- 5. 3D rendering
-- ---------------------------------------------------------------
print("\n--- 5. 3D Render ---")
local fb = sys:render{
    width  = 64,
    height = 64,
    eye    = {2, 2, 1.5},
    target = {0, 0, 0},
    fov    = 45,
    shadows = 1,
    edges   = 1,
}
fb:edge_darken()
fb:write_ppm("/tmp/alea_pipeline_render.ppm")
print(string.format("  Rendered: %dx%d -> /tmp/alea_pipeline_render.ppm", fb:width(), fb:height()))

-- Clipped render showing internal structure
local fb2 = sys:render{
    width  = 64,
    height = 64,
    eye    = {2, 2, 1.5},
    target = {0, 0, 0},
    fov    = 45,
    clips  = {{0, -1, 0, 0}},  -- clip y > 0 to show cross-section
}
fb2:edge_darken()
fb2:write_ppm("/tmp/alea_pipeline_render_cut.ppm")
print("  Clipped:  /tmp/alea_pipeline_render_cut.ppm")

-- ---------------------------------------------------------------
-- 6. Mesh export
-- ---------------------------------------------------------------
print("\n--- 6. Mesh Export ---")
local mesh = sys:mesh_sample{nx = 16, ny = 16, nz = 16}
local mi = mesh:info()
print(string.format("  Mesh: %dx%dx%d = %d voxels", mi.nx, mi.ny, mi.nz, mi.nx * mi.ny * mi.nz))
print(string.format("  Materials: %d unique", mi.num_materials))

mesh:export(0, "/tmp/alea_pipeline_mesh.msh")
mesh:export(1, "/tmp/alea_pipeline_mesh.vtk")
print("  Gmsh: /tmp/alea_pipeline_mesh.msh")
print("  VTK:  /tmp/alea_pipeline_mesh.vtk")

-- Material distribution
local mids = mesh:material_ids()
local mat_dist = {}
for _, mid in ipairs(mids) do
    mat_dist[mid] = (mat_dist[mid] or 0) + 1
end
print("\n  Voxel distribution:")
local total = #mids
local mkeys = {}
for k in pairs(mat_dist) do mkeys[#mkeys + 1] = k end
table.sort(mkeys)
for _, mid in ipairs(mkeys) do
    local name = ({[0]="void", [1]="fuel", [2]="gap", [3]="clad", [4]="moderator"})[mid] or "?"
    print(string.format("    mat %-2d (%-10s): %4d voxels (%5.1f%%)",
        mid, name, mat_dist[mid], 100 * mat_dist[mid] / total))
end

print("\n15_analysis_pipeline: OK")
