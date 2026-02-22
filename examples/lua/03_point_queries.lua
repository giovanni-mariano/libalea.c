-- 03_point_queries.lua: Point-in-cell queries, material lookup, overlap detection
--
-- Standalone: no (loads MCNP data file)
-- Usage: bin/alea examples/lua/03_point_queries.lua [mcnp_file]

print("=== Point Queries ===\n")

local filename = alea.arg and alea.arg[1]
print("Loading: " .. filename)
local sys = alea.load_mcnp(filename)
sys:print_summary()

-- Build indices for queries
sys:build_universe_index()

-- Single-cell query at the origin
print("\n--- find_cell ---")
local test_points = {
    {0, 0, 0},
    {10, 0, 0},
    {0, 50, 0},
    {100, 100, 100},
}

for _, pt in ipairs(test_points) do
    local x, y, z = pt[1], pt[2], pt[3]
    local idx = sys:find_cell(x, y, z)
    if idx then
        local info = sys:cell_info(idx)
        print(string.format("  (%6.1f, %6.1f, %6.1f) -> cell %d (mat %d, universe %d)",
            x, y, z, info.cell_id, info.material_id, info.universe_id))
    else
        print(string.format("  (%6.1f, %6.1f, %6.1f) -> no cell found", x, y, z))
    end
end

-- Material query
print("\n--- material_at ---")
local mat = sys:material_at(0, 0, 0)
if mat then
    print(string.format("  Material at origin: %d", mat))
else
    print("  No material at origin")
end

-- Multi-hit query (finds all cells containing the point, including filled)
print("\n--- find_all_cells ---")
local hits = sys:find_all_cells(0, 0, 0)
print(string.format("  Cells containing origin: %d hits", #hits))
for i, h in ipairs(hits) do
    print(string.format("    [%d] cell_id=%d mat=%d universe=%d fill=%d depth=%d",
        i, h.cell_id, h.material_id, h.universe_id, h.fill_universe, h.depth))
end

-- Point-inside test with a node
print("\n--- point_inside ---")
local info0 = sys:cell_info(0)
local root = info0.root
local inside = sys:point_inside(root, 0, 0, 0)
print(string.format("  Origin inside cell 0's region: %s", tostring(inside)))

print("\n03_point_queries: OK")
