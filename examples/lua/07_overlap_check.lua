-- 07_overlap_check.lua: Overlap detection using the geometry validator
--
-- Standalone: partially (has inline geometry, also loads data file)
-- Usage: bin/alea examples/lua/07_overlap_check.lua <mcnp_file>

print("=== Overlap Check ===\n")

-- Part 1: Check for overlaps in an inline geometry with a known overlap
print("--- Part 1: Inline Geometry with Known Overlap ---")

local mcnp_input = [[
Test overlapping geometry
1  1  -1.0  -1
2  1  -1.0  -2
3  0         1  2

1 so 10.0
2 so  8.0
]]

local sys1 = alea.load_mcnp_string(mcnp_input)
sys1:build_universe_index()
sys1:prepare_query_acceleration()
print(string.format("Loaded inline model: %d cells, %d surfaces",
    sys1:cell_count(), sys1:surface_count()))

local report1 = sys1:validate_geometry{
    ray_count = 20000,
    seed = 12345,
    max_errors = 1000,
}
local summary1 = report1:summary()
print(string.format("Overlap findings: %d",
    summary1.overlap_after_crossing or 0))

-- Part 2: Check overlaps in a data file
print("\n--- Part 2: Data File Overlap Check ---")
local filename = alea.arg and alea.arg[1]
if not filename then
    error("Pass an MCNP input file; try the model exported by 02_build_geometry.lua")
end
print("Loading: " .. filename)

local ok, sys2 = pcall(alea.load_mcnp, filename)
if ok then
    sys2:build_universe_index()
    sys2:prepare_query_acceleration()

    local report2 = sys2:validate_geometry{
        ray_count = 20000,
        seed = 12345,
        max_errors = 1000,
    }
    local summary2 = report2:summary()
    print(string.format("Overlap findings: %d",
        summary2.overlap_after_crossing or 0))
else
    print("  Could not load file (skipping): " .. tostring(sys2))
end

print("\n07_overlap_check: OK")
