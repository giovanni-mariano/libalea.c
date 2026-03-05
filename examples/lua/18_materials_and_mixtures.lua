-- 18_materials_and_mixtures.lua: Define materials with nuclides/elements, create mixtures
--
-- Standalone: yes (no data files needed)
-- Usage: bin/alea examples/lua/18_materials_and_mixtures.lua
--
-- Demonstrates:
--   - Creating materials and populating with nuclides/elements
--   - Setting density and fraction type
--   - Querying material composition
--   - Creating mixtures from base materials
--   - Assigning mixtures to cells
--   - Exporting to MCNP and OpenMC

print("=== Materials and Mixtures ===\n")

local sys = alea.create()

-- ============================================================================
-- Define materials
-- ============================================================================

-- Material 1: UO2 fuel (weight fractions)
local m1 = sys:material(1)
sys:material_set_weight_fraction(m1, true)
sys:material_add_nuclide(m1, 92235, ".80c", 0.0333)   -- U-235
sys:material_add_nuclide(m1, 92238, ".80c", 0.8815)   -- U-238
sys:material_add_nuclide(m1, 8016,  ".80c", 0.0852)   -- O-16
sys:material_set_density(m1, 10.97)

-- Material 2: Zircaloy cladding (by element)
local m2 = sys:material(2)
sys:material_set_weight_fraction(m2, true)
sys:material_add_element(m2, 40, ".80c", 0.98)  -- Zr
sys:material_add_element(m2, 50, ".80c", 0.02)  -- Sn
sys:material_set_density(m2, 6.56)

-- Material 3: Water moderator (atom fractions)
local m3 = sys:material(3)
sys:material_set_weight_fraction(m3, false)
sys:material_add_nuclide(m3, 1001, ".80c", 2.0)   -- H-1
sys:material_add_nuclide(m3, 8016, ".80c", 1.0)   -- O-16
sys:material_set_density(m3, 1.0)

-- ============================================================================
-- Query material properties
-- ============================================================================

print(string.format("Materials in system: %d", sys:material_count()))

for i = 0, sys:material_count() - 1 do
    local id = sys:material_get_id(i)
    local nn = sys:material_nuclide_count(i)
    local ne = sys:material_element_count(i)
    local density, has_density = sys:material_get_density(i)
    local wf = sys:material_is_weight_fraction(i)
    print(string.format("  M%d: %d nuclides, %d elements, density=%.2f, weight_frac=%s",
                        id, nn, ne, density, tostring(wf)))

    for j = 1, nn do
        local nuc = sys:material_nuclide_get(i, j)
        print(string.format("    nuclide %d: ZAID=%d lib=%s frac=%.4f",
                            j, nuc.zaid, nuc.library or "none", nuc.fraction))
    end
    for j = 1, ne do
        local elem = sys:material_element_get(i, j)
        print(string.format("    element %d: Z=%d lib=%s frac=%.4f",
                            j, elem.Z, elem.library or "none", elem.fraction))
    end
end

-- Find material by ID
local idx = sys:find_material(2)
print(string.format("\nfind_material(2) -> index %d", idx))

-- ============================================================================
-- Create a mixture
-- ============================================================================

-- Mix material 1 (70%) and material 3 (30%) -> mixture ID 100
local mix_id = sys:create_mixture({1, 3}, {0.7, 0.3}, 100)
print(string.format("\nCreated mixture ID %d", mix_id))

print(string.format("Mixtures in system: %d", sys:mixture_count()))
local mix_idx = sys:find_mixture(100)
print(string.format("find_mixture(100) -> index %d", mix_idx))
print(string.format("  mixture_get_id(0) = %d", sys:mixture_get_id(mix_idx)))
print(string.format("  components: %d", sys:mixture_component_count(mix_idx)))

for j = 1, sys:mixture_component_count(mix_idx) do
    local comp = sys:mixture_component_get(mix_idx, j)
    print(string.format("    component %d: material_id=%d fraction=%.2f",
                        j, comp.material_id, comp.fraction))
end

-- ============================================================================
-- Build geometry and assign mixture to a cell
-- ============================================================================

local s1 = sys:sphere(1, 0, 0, 0, 5)
local s2 = sys:sphere(2, 0, 0, 0, 10)

local cell_fuel = sys:cell{id = 1, region = sys:inside(s1), material = m1, density = 10.97}
local cell_mix  = sys:cell{id = 2, region = sys:outside(s1) * sys:inside(s2),
                           material = m3, density = 1.0}
local cell_void = sys:cell{id = 3, region = sys:outside(s2), material = -1}

-- Assign mixture to the middle cell
sys:cell_set_mixture(cell_mix, 100)
print(string.format("\nAssigned mixture %d to cell index %d", 100, cell_mix))

-- ============================================================================
-- Export
-- ============================================================================

local mcnp_out = "/tmp/alea_materials_demo.i"
sys:export_mcnp(mcnp_out)
print("\nExported MCNP to: " .. mcnp_out)

local openmc_out = "/tmp/alea_materials_demo.xml"
sys:export_openmc(openmc_out)
print("Exported OpenMC to: " .. openmc_out)

print("\n18_materials_and_mixtures: OK")
