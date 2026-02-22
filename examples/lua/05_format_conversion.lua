-- 05_format_conversion.lua: MCNP -> OpenMC -> MCNP round-trip
--
-- Standalone: no (loads MCNP data file)
-- Usage: bin/alea examples/lua/05_format_conversion.lua [mcnp_file]

print("=== Format Conversion (MCNP <-> OpenMC) ===\n")

local filename = alea.arg and alea.arg[1]
print("Loading MCNP: " .. filename)
local sys = alea.load_mcnp(filename)

local nc1 = sys:cell_count()
local ns1 = sys:surface_count()
print(string.format("  Original: %d cells, %d surfaces", nc1, ns1))

-- Step 1: Export to OpenMC XML
local openmc_file = "/tmp/alea_convert.xml"
sys:export_openmc(openmc_file)
print("\nExported to OpenMC: " .. openmc_file)

-- Step 2: Reload from OpenMC
local sys2 = alea.load_openmc(openmc_file)
local nc2 = sys2:cell_count()
local ns2 = sys2:surface_count()
print(string.format("  OpenMC reload: %d cells, %d surfaces", nc2, ns2))

-- Step 3: Export back to MCNP
local mcnp_file = "/tmp/alea_roundtrip.i"
sys2:export_mcnp(mcnp_file)
print("\nExported back to MCNP: " .. mcnp_file)

-- Step 4: Reload the MCNP output to verify
local sys3 = alea.load_mcnp(mcnp_file)
local nc3 = sys3:cell_count()
local ns3 = sys3:surface_count()
print(string.format("  MCNP reload: %d cells, %d surfaces", nc3, ns3))

-- Summary
print("\n--- Round-trip Summary ---")
print(string.format("  Step           Cells  Surfaces"))
print(string.format("  Original       %5d  %8d", nc1, ns1))
print(string.format("  After OpenMC   %5d  %8d", nc2, ns2))
print(string.format("  After MCNP     %5d  %8d", nc3, ns3))

print("\n05_format_conversion: OK")
