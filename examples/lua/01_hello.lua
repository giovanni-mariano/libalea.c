-- 01_hello.lua: Version info, system creation, configuration
--
-- Standalone: yes (no data files needed)
-- Usage: bin/alea examples/lua/01_hello.lua

print("=== Alea Lua Bindings - Hello World ===\n")

-- Version
local v = alea.version()
print("Alea version: " .. v)

-- Create an empty system
local sys = alea.create()
print("Created system: " .. tostring(sys))
print("  cells:     " .. sys:cell_count())
print("  surfaces:  " .. sys:surface_count())
print("  universes: " .. sys:universe_count())

-- Inspect default configuration
local cfg = sys:get_config()
print("\nDefault configuration:")
print(string.format("  abs_tol        = %g", cfg.abs_tol))
print(string.format("  rel_tol        = %g", cfg.rel_tol))
print(string.format("  zero_threshold = %g", cfg.zero_threshold))
print(string.format("  dedup          = %s", tostring(cfg.dedup)))
print(string.format("  log_level      = %d", cfg.log_level))

-- Modify configuration
sys:set_config({
    abs_tol  = 1e-8,
    dedup    = false,
    log_level = 0,
})

local cfg2 = sys:get_config()
print("\nAfter set_config:")
print(string.format("  abs_tol   = %g", cfg2.abs_tol))
print(string.format("  dedup     = %s", tostring(cfg2.dedup)))
print(string.format("  log_level = %d", cfg2.log_level))

-- Clean up
sys:destroy()
print("\nSystem destroyed.")
print("\n01_hello: OK")
