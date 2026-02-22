-- 09_universe_extract_merge.lua: Extract universe, merge systems, renumber IDs
--
-- Standalone: no (loads MCNP data file)
-- Usage: bin/alea -l 0 examples/lua/09_universe_extract_merge.lua [mcnp_file]
-- Note: -l 0 suppresses warnings about unused surface refs during extract

print("=== Universe Extract & Merge ===\n")

local filename = alea.arg and alea.arg[1]
print("Loading: " .. filename)
local sys = alea.load_mcnp(filename)
sys:print_summary()

-- Find available universes
local nu = sys:universe_count()
print(string.format("\nUniverses in model: %d", nu))

-- Extract the first non-zero universe we can find
print("\n--- Extracting a Universe ---")
local nc = sys:cell_count()
local target_uid = nil
for i = 0, nc - 1 do
    local info = sys:cell_info(i)
    if info.universe_id > 0 then
        target_uid = info.universe_id
        break
    end
end

if not target_uid then
    print("No non-zero universes found, using universe 0")
    target_uid = 0
end

print(string.format("Extracting universe %d...", target_uid))
local extracted = sys:extract_universe(target_uid)
local nc_ext = extracted:cell_count()
local ns_ext = extracted:surface_count()
print(string.format("  Extracted: %d cells, %d surfaces", nc_ext, ns_ext))

-- Clone and renumber the extracted system
print("\n--- Cloning and Renumbering ---")
local copy = extracted:clone()
copy:renumber_cells(10000)
copy:renumber_surfaces(10000)
print(string.format("  Renumbered clone: %d cells, %d surfaces",
    copy:cell_count(), copy:surface_count()))

-- Show renumbered cell IDs
print("  Cell IDs after renumber:")
for i = 0, math.min(4, copy:cell_count() - 1) do
    local info = copy:cell_info(i)
    print(string.format("    [%d] id=%d mat=%d", i, info.cell_id, info.material_id))
end

-- Merge the copy back into the extracted system
print("\n--- Merging ---")
local nc_pre = extracted:cell_count()
extracted:merge(copy, 0)
local nc_post = extracted:cell_count()
print(string.format("  Before merge: %d cells", nc_pre))
print(string.format("  After merge:  %d cells", nc_post))

-- ID offset operations
print("\n--- ID Offset Operations ---")
local sys2 = alea.create()
local s1 = sys2:sphere(1, 0, 0, 0, 5)
sys2:cell{id = 1, region = sys2:inside(s1), material = 1, density = 1.0}
sys2:cell{id = 2, region = sys2:outside(s1), material = 0, density = 0.0}

print("  Before offset:")
for i = 0, sys2:cell_count() - 1 do
    local info = sys2:cell_info(i)
    print(string.format("    cell_id=%d mat=%d", info.cell_id, info.material_id))
end

sys2:offset_cell_ids(100)
sys2:offset_surface_ids(200)
sys2:offset_material_ids(50)

print("  After offset (cells+100, surfaces+200, materials+50):")
for i = 0, sys2:cell_count() - 1 do
    local info = sys2:cell_info(i)
    print(string.format("    cell_id=%d mat=%d", info.cell_id, info.material_id))
end

print("\n09_universe_extract_merge: OK")
