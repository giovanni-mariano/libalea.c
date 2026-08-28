-- 06_volume_estimation.lua: Bounding sphere, Monte Carlo volume estimation
--
-- Standalone: no (loads MCNP data file)
-- Usage: bin/alea examples/lua/06_volume_estimation.lua [mcnp_file]

print("=== Volume Estimation ===\n")

local filename = alea.arg and alea.arg[1]
print("Loading: " .. filename)
local sys = alea.load_mcnp(filename)
sys:build_universe_index()

local nc = sys:cell_count()
print(string.format("Cells: %d\n", nc))

-- Compute bounding sphere
print("--- Bounding Sphere ---")
local cx, cy, cz, r = sys:bounding_sphere(1e-3)
print(string.format("  Center: (%.2f, %.2f, %.2f)", cx, cy, cz))
print(string.format("  Radius: %.2f", r))

-- Estimate volumes with Monte Carlo ray tracing
local n_rays = 10000
print(string.format("\n--- Volume Estimation (%d rays) ---", n_rays))
local volumes = sys:estimate_volumes(n_rays)

-- Print top 20 cells by volume
print(string.format("\nTop 20 cells by volume:"))
print(string.format("  %-6s  %-8s  %-6s  %12s  %10s", "Index", "Cell ID", "Mat", "Volume", "Rel Error"))
print(string.rep("-", 52))

-- Collect and sort by volume
local sorted = {}
for i = 1, #volumes do
    sorted[i] = {
        index = volumes[i].path_id,
        cell_id = volumes[i].cell_id,
        mat = volumes[i].material_id,
        vol = volumes[i].volume,
        err = volumes[i].rel_error,
    }
end
table.sort(sorted, function(a, b) return a.vol > b.vol end)

for i = 1, math.min(20, #sorted) do
    local e = sorted[i]
    if e.vol > 0 then
        print(string.format("  %-6d  %-8d  %-6d  %12.2f  %9.1f%%",
            e.index, e.cell_id, e.mat, e.vol, e.err * 100))
    end
end

-- Total volume
local total = 0
for _, e in ipairs(sorted) do total = total + e.vol end
print(string.format("\nTotal estimated volume: %.2f", total))

print("\n06_volume_estimation: OK")
