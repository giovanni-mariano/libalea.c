-- 07_overlap_check.lua: Overlap detection using find_overlaps and load_mcnp_string
--
-- Standalone: partially (has inline geometry, also loads data file)
-- Usage: bin/alea examples/lua/07_overlap_check.lua [mcnp_file]

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
sys1:build_spatial_index()
print(string.format("Loaded inline model: %d cells, %d surfaces",
    sys1:cell_count(), sys1:surface_count()))

local overlaps1 = sys1:find_overlaps()
print(string.format("Overlapping pairs found: %d", #overlaps1))
for i, pair in ipairs(overlaps1) do
    local a, b = pair[1], pair[2]
    local info_a = sys1:cell_info(a)
    local info_b = sys1:cell_info(b)
    print(string.format("  Pair %d: cell %d (idx %d) <-> cell %d (idx %d)",
        i, info_a.cell_id, a, info_b.cell_id, b))
end

-- Part 2: Check overlaps in a data file
print("\n--- Part 2: Data File Overlap Check ---")
local filename = alea.arg and alea.arg[1]
print("Loading: " .. filename)

local ok, sys2 = pcall(alea.load_mcnp, filename)
if ok then
    sys2:build_universe_index()
    sys2:build_spatial_index()

    local overlaps2 = sys2:find_overlaps()
    print(string.format("Overlapping pairs found: %d", #overlaps2))
    for i, pair in ipairs(overlaps2) do
        if i > 10 then
            print(string.format("  ... and %d more", #overlaps2 - 10))
            break
        end
        local a, b = pair[1], pair[2]
        local info_a = sys2:cell_info(a)
        local info_b = sys2:cell_info(b)
        print(string.format("  Pair %d: cell %d <-> cell %d",
            i, info_a.cell_id, info_b.cell_id))
    end
else
    print("  Could not load file (skipping): " .. tostring(sys2))
end

print("\n07_overlap_check: OK")
