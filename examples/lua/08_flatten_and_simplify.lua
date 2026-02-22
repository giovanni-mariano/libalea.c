-- 08_flatten_and_simplify.lua: Flatten universes, tighten bboxes, remove empty cells
--
-- Standalone: no (loads MCNP data file)
-- Usage: bin/alea examples/lua/08_flatten_and_simplify.lua [mcnp_file]

print("=== Flatten and Simplify ===\n")

local filename = alea.arg and alea.arg[1]
print("Loading: " .. filename)
local sys = alea.load_mcnp(filename)

local nc_before = sys:cell_count()
local ns_before = sys:surface_count()
local nu_before = sys:universe_count()
print(string.format("Before flatten: %d cells, %d surfaces, %d universes",
    nc_before, ns_before, nu_before))

-- Flatten universe 0 (merges all fill hierarchies into a single universe)
print("\n--- Flattening Universe 0 ---")
sys:flatten(0)

local nc_flat = sys:cell_count()
local ns_flat = sys:surface_count()
local nu_flat = sys:universe_count()
print(string.format("After flatten: %d cells, %d surfaces, %d universes",
    nc_flat, ns_flat, nu_flat))

-- Tighten bounding boxes
print("\n--- Tightening Bounding Boxes ---")
sys:tighten_all_bboxes(1e-3)
print("Bounding boxes tightened.")

-- Show bbox for first few cells
print("\nSample cell bounding boxes:")
for i = 0, math.min(4, nc_flat - 1) do
    local info = sys:cell_info(i)
    local bb = info.bbox
    print(string.format("  Cell %d: [%.1f, %.1f] x [%.1f, %.1f] x [%.1f, %.1f]",
        info.cell_id, bb.min_x, bb.max_x, bb.min_y, bb.max_y, bb.min_z, bb.max_z))
end

-- Expand macrobodies (if any)
print("\n--- Expanding Macrobodies ---")
local expanded = sys:expand_macrobodies()
print(string.format("Expanded %d macrobodies", expanded))

-- Split union cells (cells with top-level unions become multiple cells)
print("\n--- Splitting Union Cells ---")
sys:split_union_cells()
local nc_split = sys:cell_count()
print(string.format("After split: %d cells (was %d)", nc_split, nc_flat))

-- Export the simplified model
local outfile = "/tmp/alea_flattened.i"
sys:export_mcnp(outfile)
print("\nExported flattened model to: " .. outfile)

print("\n08_flatten_and_simplify: OK")
