-- 04_model_inspection.lua: Cell/material/universe census, CSG tree walk
--
-- Standalone: no (loads MCNP data file)
-- Usage: bin/alea examples/lua/04_model_inspection.lua [mcnp_file]

print("=== Model Inspection ===\n")

local filename = alea.arg and alea.arg[1]
print("Loading: " .. filename)
local sys = alea.load_mcnp(filename)
sys:print_summary()

-- Cell census
local nc = sys:cell_count()
print(string.format("\n--- Cell Census (%d cells) ---", nc))

local mat_count = {}    -- material_id -> count
local univ_count = {}   -- universe_id -> count
local fill_count = 0

for i = 0, nc - 1 do
    local info = sys:cell_info(i)
    mat_count[info.material_id] = (mat_count[info.material_id] or 0) + 1
    univ_count[info.universe_id] = (univ_count[info.universe_id] or 0) + 1
    if info.fill_universe > 0 then
        fill_count = fill_count + 1
    end
end

-- Materials summary
local mats = {}
for k in pairs(mat_count) do mats[#mats + 1] = k end
table.sort(mats)
print(string.format("\nMaterials used: %d", #mats))
for _, m in ipairs(mats) do
    print(string.format("  mat %4d: %d cells", m, mat_count[m]))
end

-- Universes summary
local univs = {}
for k in pairs(univ_count) do univs[#univs + 1] = k end
table.sort(univs)
print(string.format("\nUniverses: %d", #univs))
for _, u in ipairs(univs) do
    print(string.format("  universe %4d: %d cells", u, univ_count[u]))
end

print(string.format("\nFill cells: %d", fill_count))

-- CSG tree walk for a single cell
print("\n--- CSG Tree Walk (cell 0) ---")
local function walk_tree(sys, node, depth)
    local indent = string.rep("  ", depth)
    local op = sys:node_operation(node)

    if op == "primitive" then
        local sense = sys:node_sense(node)
        local sid = sys:node_surface_id(node)
        local sense_str = sense > 0 and "+" or "-"
        print(string.format("%s%s%d", indent, sense_str, sid))
    elseif op == "complement" then
        print(indent .. "COMPLEMENT")
        walk_tree(sys, sys:node_left(node), depth + 1)
    else
        print(indent .. string.upper(op))
        walk_tree(sys, sys:node_left(node), depth + 1)
        walk_tree(sys, sys:node_right(node), depth + 1)
    end
end

local root = sys:cell_info(0).root
walk_tree(sys, root, 1)

-- Cell lookup by ID
print("\n--- Cell Lookup ---")
local cell_id = sys:get_cell_id(0)
local found = sys:cell_find(cell_id)
print(string.format("  cell_find(%d) -> index %s", cell_id, tostring(found)))

print("\n04_model_inspection: OK")
