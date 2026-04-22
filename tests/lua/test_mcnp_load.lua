-- test_mcnp_load.lua: Load a real MCNP model and exercise the API

local filename = alea.arg and alea.arg[1] or "tests/data/mcnp_flatten_ids.mcnp"
print("Loading: " .. filename)

local sys = alea.load_mcnp(filename)
print(tostring(sys))
sys:print_summary()

-- Basic counts
local nc = sys:cell_count()
local ns = sys:surface_count()
local nu = sys:universe_count()
print(string.format("  cells=%d  surfaces=%d  universes=%d", nc, ns, nu))
assert(nc > 0, "should have cells")
assert(ns > 0, "should have surfaces")

-- Cell info for first few cells
print("\nFirst 5 cells:")
for i = 0, math.min(4, nc - 1) do
    local info = sys:cell_info(i)
    print(string.format("  [%d] id=%d mat=%d density=%.4g universe=%d fill=%d",
        i, info.cell_id, info.material_id, info.density,
        info.universe_id, info.fill_universe))
end

-- Build index and query
sys:build_universe_index()
print("\nPoint query at origin:")
local cell_id, material = sys:find_cell_at(0, 0, 0)
if cell_id then
    print(string.format("  cell_id=%d material=%d", cell_id, material))
else
    print("  no cell found at origin")
end

-- Config
local cfg = sys:get_config()
print(string.format("\nConfig: dedup=%s abs_tol=%g", tostring(cfg.dedup), cfg.abs_tol))

-- Export round-trip
local outfile = "/tmp/alea_lua_test_output.i"
sys:export_mcnp(outfile)
print("\nExported to: " .. outfile)

-- Reload and compare
local sys2 = alea.load_mcnp(outfile)
local nc2 = sys2:cell_count()
local ns2 = sys2:surface_count()
print(string.format("Reloaded: cells=%d surfaces=%d", nc2, ns2))

print("\ntest_mcnp_load: OK")
