-- test_render.lua: Render module tests

-- Build a sphere geometry
local sys = alea.create()
local s1 = sys:sphere(1, 0, 0, 0, 5)
sys:cell{id = 1, region = sys:inside(s1),  material = 1, density = 10.0}
sys:cell{id = 2, region = sys:outside(s1), material = 0, density = 0.0}
sys:build_universe_index()

-- Default config
local cfg = alea.render_config()
assert(type(cfg) == "table", "render_config should return a table")
assert(cfg.width > 0, "default width should be positive")
assert(cfg.height > 0, "default height should be positive")

-- Render with small image
local fb = sys:render{width = 32, height = 32}
assert(fb:width() == 32, "fb width should be 32")
assert(fb:height() == 32, "fb height should be 32")

-- Edge darken
fb:edge_darken()

-- Write to temp files
local tmpdir = "/tmp"
fb:write(tmpdir .. "/alea_test_render.ppm")
fb:write_ppm(tmpdir .. "/alea_test_render2.ppm")
fb:write_bmp(tmpdir .. "/alea_test_render.bmp")

-- Camera setup
local cam = sys:render_camera_setup{width = 64, height = 64}
assert(type(cam) == "table", "camera should be a table")
assert(type(cam.eye) == "table", "camera should have eye")
assert(#cam.eye == 3, "eye should have 3 components")
assert(type(cam.target) == "table", "camera should have target")
assert(type(cam.forward) == "table", "camera should have forward")

-- Render with custom eye/target
local fb2 = sys:render{
    width = 16,
    height = 16,
    eye = {20, 20, 20},
    target = {0, 0, 0},
    fov = 45,
}
assert(fb2:width() == 16, "custom render width should be 16")

sys:destroy()
print("test_render: OK")
