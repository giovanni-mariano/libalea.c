-- 16_cell_comments.lua: Inspect and report cell comments from MCNP models
--
-- Demonstrates accessing "C" comment blocks and inline "$" comments
-- preserved during MCNP parsing. Useful for auditing large models
-- where cells carry descriptive annotations.
--
-- Standalone: no (loads MCNP data file)
-- Usage: bin/alea examples/lua/16_cell_comments.lua <mcnp_file>

local filename = alea.arg and alea.arg[1]
if not filename then
    print("Usage: bin/alea 16_cell_comments.lua <mcnp_file>")
    os.exit(1)
end

print("=== Cell Comment Report ===\n")
print("Loading: " .. filename)
local sys = alea.load_mcnp(filename)

local nc = sys:cell_count()
local nu = sys:universe_count()
print(string.format("Cells: %d, Universes: %d\n", nc, nu))

-- Collect cells with comments
local commented = 0
local inline_commented = 0
local both = 0

for i = 0, nc - 1 do
    local info = sys:cell_info(i)
    local has_c = info.comments ~= nil
    local has_d = info.inline_comment ~= nil
    if has_c then commented = commented + 1 end
    if has_d then inline_commented = inline_commented + 1 end
    if has_c and has_d then both = both + 1 end
end

print(string.format("Cells with C comments:  %d / %d", commented, nc))
print(string.format("Cells with $ comments:  %d / %d", inline_commented, nc))
print(string.format("Cells with both:        %d / %d", both, nc))

-- Per-universe comment summary
print("\n--- Comments by Universe ---")
for u = 0, nu - 1 do
    local uinfo = sys:universe_info(u)
    local uid = uinfo.universe_id
    local cells = sys:cells_in_universe(uid)
    local c_count, d_count = 0, 0

    for _, idx in ipairs(cells) do
        local info = sys:cell_info(idx)
        if info.comments then c_count = c_count + 1 end
        if info.inline_comment then d_count = d_count + 1 end
    end

    if c_count > 0 or d_count > 0 then
        print(string.format("  Universe %4d: %d cells, %d C-comments, %d $-comments",
                            uid, #cells, c_count, d_count))
    end
end

-- Show detailed comments for first few commented cells
print("\n--- Sample Cell Comments ---")
local shown = 0
for i = 0, nc - 1 do
    if shown >= 5 then break end
    local info = sys:cell_info(i)
    if info.comments or info.inline_comment then
        print(string.format("\nCell %d (mat=%d, U=%d):",
                            info.cell_id, info.material_id, info.universe_id))
        if info.comments then
            -- Print each comment line with indentation
            for line in info.comments:gmatch("[^\n]+") do
                print("  " .. line)
            end
        end
        if info.inline_comment then
            print("  $ " .. info.inline_comment)
        end
        shown = shown + 1
    end
end

if shown == 0 then
    print("  (no comments found in this model)")
end

print("\n=== Done ===")
