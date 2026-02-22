-- test_basic.lua: Basic alea module tests

-- Version
local v = alea.version()
assert(type(v) == "string", "version should be a string")
print("  version: " .. v)

-- Create empty system
local sys = alea.create()
assert(sys, "create should return a system")
assert(sys:cell_count() == 0, "new system has 0 cells")
assert(sys:surface_count() == 0, "new system has 0 surfaces")

-- tostring
local s = tostring(sys)
assert(s:find("System"), "tostring should contain 'System'")

-- Config round-trip
local cfg = sys:get_config()
assert(type(cfg) == "table", "get_config returns table")
assert(cfg.dedup == true, "default dedup is true")
sys:set_config({dedup = false})
cfg = sys:get_config()
assert(cfg.dedup == false, "dedup was set to false")

-- Destroy
sys:destroy()
local ok, err = pcall(function() sys:cell_count() end)
assert(not ok, "should error after destroy")

print("test_basic: OK")
