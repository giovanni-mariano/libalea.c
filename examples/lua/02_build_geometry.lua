-- 02_build_geometry.lua: Build a shielded sphere from scratch, export to MCNP
--
-- Standalone: yes (no data files needed)
-- Usage: bin/alea examples/lua/02_build_geometry.lua
--
-- Geometry: concentric spheres (fuel / clad / moderator / void)
--   Cell 1: fuel      (mat 1)  r < 5
--   Cell 2: clad      (mat 2)  5 <= r < 6
--   Cell 3: moderator (mat 3)  6 <= r < 15
--   Cell 4: void      (mat 0)  r >= 15

print("=== Build Geometry from Scratch ===\n")

local sys = alea.create()

-- Create spherical surfaces
local s1 = sys:sphere(1, 0, 0, 0, 5)   -- inner fuel boundary
local s2 = sys:sphere(2, 0, 0, 0, 6)   -- clad outer boundary
local s3 = sys:sphere(3, 0, 0, 0, 15)  -- moderator outer boundary

-- Build halfspace nodes
local in_s1  = sys:inside(s1)
local in_s2  = sys:inside(s2)
local in_s3  = sys:inside(s3)
local out_s1 = sys:outside(s1)
local out_s2 = sys:outside(s2)
local out_s3 = sys:outside(s3)

-- Define cell regions using operator overloads
local fuel_region = in_s1                -- r < 5
local clad_region = out_s1 * in_s2       -- 5 <= r < 6
local mod_region  = out_s2 * in_s3       -- 6 <= r < 15
local void_region = out_s3              -- r >= 15

-- Create cells
sys:cell{id = 1, region = fuel_region, material = 1, density = 10.97}
sys:cell{id = 2, region = clad_region, material = 2, density = 6.56}
sys:cell{id = 3, region = mod_region,  material = 3, density = 1.0}
sys:cell{id = 4, region = void_region, material = 0, density = 0.0}

sys:print_summary()

-- Validate
local issues = sys:validate()
print(string.format("Validation issues: %d", issues))

-- Export to MCNP
local outfile = "/tmp/alea_shielded_sphere.i"
sys:export_mcnp(outfile)
print("Exported MCNP to: " .. outfile)

-- Also export to OpenMC
local outxml = "/tmp/alea_shielded_sphere.xml"
sys:export_openmc(outxml)
print("Exported OpenMC to: " .. outxml)

print("\n02_build_geometry: OK")
