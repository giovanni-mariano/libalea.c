-- test_render.lua: Render module tests

-- Build a sphere geometry
local sys = alea.create()
local s1 = sys:sphere(1, 0, 0, 0, 5)
local m1 = sys:material(1)
sys:cell{id = 1, region = sys:inside(s1),  material = m1, density = 10.0}
sys:cell{id = 2, region = sys:outside(s1), material = -1, density = 0.0}
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
local ppm_file = join_path(tmpdir, "alea_test_render.ppm")
local ppm_file2 = join_path(tmpdir, "alea_test_render2.ppm")
local bmp_file = join_path(tmpdir, "alea_test_render.bmp")
fb:write(ppm_file)
fb:write_ppm(ppm_file2)
fb:write_bmp(bmp_file)
os.remove(ppm_file)
os.remove(ppm_file2)
os.remove(bmp_file)

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
