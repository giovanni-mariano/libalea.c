-- 10_parametric_geometry.lua: Parameterized pin-cell array with void generation
--
-- Standalone: yes (no data files needed)
-- Usage: bin/alea examples/lua/10_parametric_geometry.lua
--
-- Builds a 3x3 array of fuel pins in a moderator box.
-- Parameters control pin radius, pitch, and height.

print("=== Parametric Pin-Cell Array ===\n")

-- Parameters (easily adjustable)
local params = {
    fuel_radius   = 0.4,    -- cm
    clad_radius   = 0.5,    -- cm
    pitch         = 1.26,   -- cm (center-to-center)
    height        = 10.0,   -- cm (half-height)
    nx            = 3,      -- pins in X
    ny            = 3,      -- pins in Y
    fuel_mat      = 1,
    fuel_density  = 10.97,
    clad_mat      = 2,
    clad_density  = 6.56,
    mod_mat       = 3,
    mod_density   = 1.0,
}

print("Parameters:")
for k, v in pairs(params) do
    print(string.format("  %-14s = %s", k, tostring(v)))
end

local sys = alea.create()
local cell_id = 1
local surf_id = 1

-- Create top and bottom planes for the active region
local pz_bot = sys:plane(surf_id, 0, 0, 1, -params.height); surf_id = surf_id + 1
local pz_top = sys:plane(surf_id, 0, 0, 1,  params.height); surf_id = surf_id + 1
local axial_region = sys:inside(pz_bot) * sys:outside(pz_top)
-- Note: plane(id, a, b, c, d) => inside means a*x+b*y+c*z+d >= 0
-- pz_bot: z + (-height) >= 0 => z >= height... we need to think about sense
-- Actually: inside(pz_bot) = z - height >= 0 for "z + d >= 0"
-- Let's use half() for explicit control:
-- For z >= -height: plane coeffs (0,0,1, height), inside = 0*x + 0*y + 1*z + height >= 0
-- For z <= +height: plane coeffs (0,0,1, -height), outside = NOT(z - height >= 0) = z < height

-- Rebuild with correct sense
sys:destroy()
sys = alea.create()
surf_id = 1
cell_id = 1

local pbot = sys:plane(surf_id, 0, 0, 1, params.height); surf_id = surf_id + 1
local ptop = sys:plane(surf_id, 0, 0, 1, -params.height); surf_id = surf_id + 1
local in_axial = sys:inside(pbot) * sys:outside(ptop)

-- Create pins
print(string.format("\nGenerating %dx%d pin array...", params.nx, params.ny))
local all_pin_nodes = {}

local x0 = -((params.nx - 1) / 2) * params.pitch
local y0 = -((params.ny - 1) / 2) * params.pitch

for iy = 0, params.ny - 1 do
    for ix = 0, params.nx - 1 do
        local cx = x0 + ix * params.pitch
        local cy = y0 + iy * params.pitch

        -- Fuel cylinder
        local s_fuel = sys:cylinder_z(surf_id, cx, cy, params.fuel_radius); surf_id = surf_id + 1
        -- Clad cylinder
        local s_clad = sys:cylinder_z(surf_id, cx, cy, params.clad_radius); surf_id = surf_id + 1

        -- Fuel cell: inside fuel cylinder AND in axial range
        local fuel_node = sys:inside(s_fuel) * in_axial
        sys:cell{
            id = cell_id, region = fuel_node,
            material = params.fuel_mat, density = params.fuel_density,
        }
        cell_id = cell_id + 1

        -- Clad cell: between fuel and clad cylinders AND in axial range
        local clad_node = sys:outside(s_fuel) * sys:inside(s_clad) * in_axial
        sys:cell{
            id = cell_id, region = clad_node,
            material = params.clad_mat, density = params.clad_density,
        }
        cell_id = cell_id + 1

        -- Track outer clad surface for moderator exclusion
        all_pin_nodes[#all_pin_nodes + 1] = sys:inside(s_clad)
    end
end

-- Moderator region: inside the box, outside all pins, within axial range
local half_x = (params.nx * params.pitch) / 2
local half_y = (params.ny * params.pitch) / 2

local box_surf = sys:box(surf_id, -half_x, half_x, -half_y, half_y,
    -params.height, params.height)
surf_id = surf_id + 1

local mod_region = sys:inside(box_surf)
for _, pin_node in ipairs(all_pin_nodes) do
    mod_region = mod_region - pin_node
end

sys:cell{
    id = cell_id, region = mod_region,
    material = params.mod_mat, density = params.mod_density,
}
cell_id = cell_id + 1

-- Void outside the box
sys:cell{
    id = cell_id, region = sys:outside(box_surf),
    material = 0, density = 0.0,
}
cell_id = cell_id + 1

sys:print_summary()

-- Validate
local issues = sys:validate()
print(string.format("Validation issues: %d", issues))

-- Export to both formats
local mcnp_out = "/tmp/alea_pin_array.i"
local openmc_out = "/tmp/alea_pin_array.xml"
sys:export_mcnp(mcnp_out)
sys:export_openmc(openmc_out)
print(string.format("\nExported MCNP:  %s", mcnp_out))
print(string.format("Exported OpenMC: %s", openmc_out))

print(string.format("\nGenerated %d cells, %d surfaces",
    sys:cell_count(), sys:surface_count()))

print("\n10_parametric_geometry: OK")
