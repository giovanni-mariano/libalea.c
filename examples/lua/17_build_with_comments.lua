-- 17_build_with_comments.lua: Build geometry with cell comments, export to MCNP
--
-- Demonstrates programmatic cell construction with "C" comment blocks
-- and inline "$" annotations that are preserved in the exported MCNP file.
--
-- Standalone: yes (no data files needed)
-- Usage: bin/alea examples/lua/17_build_with_comments.lua

print("=== Build Geometry with Comments ===\n")

local sys = alea.create()

-- Materials (returns material index for use in cell definitions)
local mat_uo2 = sys:material(1)    -- UO2 fuel
local mat_zr  = sys:material(2)    -- Zircaloy-4
local mat_h2o = sys:material(3)    -- Light water

-- Surfaces
local s1 = sys:sphere(1, 0, 0, 0, 0.41)   -- fuel pellet
local s2 = sys:sphere(2, 0, 0, 0, 0.47)   -- clad outer
local s3 = sys:sphere(3, 0, 0, 0, 0.55)   -- coolant boundary
local s4 = sys:sphere(4, 0, 0, 0, 100)    -- world boundary

-- Regions
local fuel_r = sys:inside(s1)
local clad_r = sys:outside(s1) * sys:inside(s2)
local cool_r = sys:outside(s2) * sys:inside(s3)
local void_r = sys:outside(s3) * sys:inside(s4)
local grav_r = sys:outside(s4)

-- Cells with comments
sys:cell{
    id       = 1,
    region   = fuel_r,
    material = mat_uo2,
    density  = -10.97,
    comment  = "c Fuel pellet\nc Material: UO2 enriched 4.0%\nc Density: 10.97 g/cm3\n",
    inline_comment = "UO2 fuel"
}

sys:cell{
    id       = 2,
    region   = clad_r,
    material = mat_zr,
    density  = -6.56,
    comment  = "c Cladding\nc Material: Zircaloy-4\n",
    inline_comment = "Zr clad"
}

sys:cell{
    id       = 3,
    region   = cool_r,
    material = mat_h2o,
    density  = -1.0,
    comment  = "c Coolant channel\nc Material: Light water\n",
    inline_comment = "H2O coolant"
}

sys:cell{
    id       = 4,
    region   = void_r,
    comment  = "c Void region outside coolant\n",
}

sys:cell{
    id       = 5,
    region   = grav_r,
    comment  = "c Graveyard\n",
    inline_comment = "kill zone"
}

sys:print_summary()

-- Verify comments were set
print("\n--- Verify Comments ---")
for i = 0, sys:cell_count() - 1 do
    local info = sys:cell_info(i)
    local c_mark = info.comments and "C" or "-"
    local d_mark = info.inline_comment and "$" or "-"
    local d_text = info.inline_comment or ""
    print(string.format("  Cell %2d: [%s%s] %s", info.cell_id, c_mark, d_mark, d_text))
end

-- Export to MCNP
local outfile = "/tmp/alea_commented_model.i"
sys:export_mcnp(outfile)
print("\nExported to: " .. outfile)

-- Show the exported file
print("\n--- Exported MCNP (cell block) ---")
local f = io.open(outfile, "r")
if f then
    for line in f:lines() do
        print(line)
        -- Stop after blank line (end of cell block)
        if line:match("^%s*$") then break end
    end
    f:close()
end

print("\n17_build_with_comments: OK")
