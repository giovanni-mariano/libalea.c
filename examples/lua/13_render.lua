-- 13_render.lua: 3D rendering of CSG geometry
--
-- Standalone: yes (no data files needed)
-- Usage: bin/alea examples/lua/13_render.lua
--
-- Geometry: concentric spheres (fuel / moderator / void)
--   Cell 1: fuel      (mat 1)  r < 5
--   Cell 2: moderator (mat 2)  5 <= r < 10
--   Cell 3: void      (mat 0)  r >= 10
--
-- Writes images to /tmp/alea_example_render_*.{ppm,bmp}

print("=== 3D Rendering ===\n")

local sys = alea.create()

-- Build concentric spheres
local s1 = sys:sphere(1, 0, 0, 0, 5)
local s2 = sys:sphere(2, 0, 0, 0, 10)

sys:cell{id = 1, region = sys:inside(s1), material = 1, density = 10.0}
sys:cell{id = 2, region = sys:outside(s1) * sys:inside(s2), material = 2, density = 1.0}
sys:cell{id = 3, region = sys:outside(s2), material = 0, density = 0.0}

sys:build_universe_index()

-- Inspect default config
print("--- Default Config ---")
local cfg = alea.render_config()
print(string.format("  Default size: %dx%d", cfg.width, cfg.height))
print(string.format("  FOV: %.1f  Shadows: %d  Edges: %d", cfg.fov, cfg.shadows, cfg.edges))
print(string.format("  Color mode: %d  Render mode: %d", cfg.color_mode, cfg.render_mode))

-- Auto camera setup
print("\n--- Camera Setup ---")
local cam = sys:render_camera_setup{width = 128, height = 128}
print(string.format("  Eye:     (%.2f, %.2f, %.2f)", cam.eye[1], cam.eye[2], cam.eye[3]))
print(string.format("  Target:  (%.2f, %.2f, %.2f)", cam.target[1], cam.target[2], cam.target[3]))
print(string.format("  Forward: (%.2f, %.2f, %.2f)", cam.forward[1], cam.forward[2], cam.forward[3]))
print(string.format("  Auto-fit: %s", tostring(cam.auto_fit)))

-- Basic render with auto camera
print("\n--- Basic Render (auto camera) ---")
local fb = sys:render{width = 64, height = 64}
print(string.format("  Framebuffer: %dx%d", fb:width(), fb:height()))

fb:edge_darken()
fb:write_ppm("/tmp/alea_example_render_basic.ppm")
print("  Written: /tmp/alea_example_render_basic.ppm")

-- Render with explicit camera position
print("\n--- Custom Camera ---")
local fb2 = sys:render{
    width  = 64,
    height = 64,
    eye    = {25, 25, 15},
    target = {0, 0, 0},
    fov    = 45,
    shadows = 1,
    edges   = 1,
}
fb2:edge_darken()
fb2:write_bmp("/tmp/alea_example_render_custom.bmp")
print(string.format("  Framebuffer: %dx%d", fb2:width(), fb2:height()))
print("  Written: /tmp/alea_example_render_custom.bmp")

-- Render with clipping plane (cut away top half)
print("\n--- Clipping Plane (z < 0 visible) ---")
local fb3 = sys:render{
    width  = 64,
    height = 64,
    eye    = {25, 25, 15},
    target = {0, 0, 0},
    fov    = 45,
    clips  = {{0, 0, -1, 0}},  -- visible where -z + 0 > 0, i.e. z < 0
}
fb3:edge_darken()
fb3:write_ppm("/tmp/alea_example_render_clipped.ppm")
print("  Written: /tmp/alea_example_render_clipped.ppm")

-- Multiple output formats
print("\n--- Output Formats ---")
fb:write_ppm("/tmp/alea_example_render.ppm")
fb:write_bmp("/tmp/alea_example_render.bmp")
print("  PPM: /tmp/alea_example_render.ppm")
print("  BMP: /tmp/alea_example_render.bmp")

print("\n13_render: OK")
