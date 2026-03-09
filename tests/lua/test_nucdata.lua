-- test_nucdata.lua: Nuclear data module tests
--
-- Tests the nucdata Lua bindings without requiring ACE data files.
-- Verifies type registration, factory functions, and API surface.

-- Factory functions exist
assert(type(alea.nuc_load_xsdir) == "function", "nuc_load_xsdir should be a function")
assert(type(alea.nuc_load_xsdir_dir) == "function", "nuc_load_xsdir_dir should be a function")
assert(type(alea.nuc_material) == "function", "nuc_material should be a function")
assert(type(alea.nuc_mg_create) == "function", "nuc_mg_create should be a function")

-- SVG plot factory
assert(type(alea.svg_plot) == "function", "svg_plot should be a function")

-- NucMaterial creation and basic API
local mat = alea.nuc_material()
assert(mat, "nuc_material should return a material")
assert(mat:n_components() == 0, "new material has 0 components")

-- Multigroup creation
local bounds = {20.0, 1.0, 0.1, 1e-5, 1e-11}
local mg = alea.nuc_mg_create(4, bounds)
assert(mg, "nuc_mg_create should return a multigroup")
assert(mg:n_groups() == 4, "should have 4 groups")

-- Group data starts at zero before collapse
assert(mg:sigma_t(1) == 0, "sigma_t should be 0 before collapse")
assert(mg:sigma_t(4) == 0, "sigma_t should be 0 before collapse")
assert(mg:scatter(1, 1) == 0, "scatter should be 0 before collapse")

-- Bounds checking on group indices
local ok, err = pcall(function() mg:sigma_t(0) end)
assert(not ok, "group index 0 should error (1-based)")
local ok2, err2 = pcall(function() mg:sigma_t(5) end)
assert(not ok2, "group index 5 should error (only 4 groups)")

-- SVG plot creation and API
local plot = alea.svg_plot()
assert(plot, "svg_plot should return a plot")
assert(plot:n_curves() == 0, "new plot has 0 curves")

plot:set_title("Test Plot")
plot:set_x_label("X")
plot:set_y_label("Y")
plot:set_scale("log", "linear")
plot:set_size(800, 400)

-- Add a curve
plot:add({1, 2, 3, 4}, {10, 20, 30, 40}, "Line 1", "#ff0000")
assert(plot:n_curves() == 1, "should have 1 curve after add")

-- Add with options
plot:add({1, 2, 3}, {5, 15, 25}, "Dashed", "#0000ff", {dashed = true, width = 2.0})
assert(plot:n_curves() == 2, "should have 2 curves")

-- Mismatched array lengths should error
local ok3, _ = pcall(function()
    plot:add({1, 2}, {1, 2, 3}, "Bad", "#000")
end)
assert(not ok3, "mismatched x/y lengths should error")

-- Write to temp file
local tmpfile = os.tmpname() .. ".svg"
plot:write(tmpfile)
local f = io.open(tmpfile, "r")
assert(f, "SVG file should exist")
local content = f:read("*a")
f:close()
os.remove(tmpfile)
assert(content:find("<svg"), "file should contain SVG tag")
assert(content:find("Test Plot"), "file should contain title")

-- xsdir load with bad path should error
local ok4, _ = pcall(function() alea.nuc_load_xsdir("/nonexistent/path") end)
assert(not ok4, "loading nonexistent xsdir should error")

print("test_nucdata: OK")
