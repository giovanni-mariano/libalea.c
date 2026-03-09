-- 19_nucdata.lua: Nuclear data — load ACE cross sections, query, build materials
--
-- Standalone: no (requires ACE data files and xsdir)
-- Usage: bin/alea examples/lua/19_nucdata.lua <xsdir_path> [zaid]
--
-- Demonstrates:
--   - Loading an xsdir directory file
--   - Loading a nuclide from ACE data
--   - Querying microscopic cross sections
--   - Building a material and computing macroscopic XS
--   - Multigroup collapse

print("=== Nuclear Data Module ===\n")

-- Parse arguments
local xsdir_path = arg and arg[1]
if not xsdir_path then
    print("Usage: bin/alea examples/lua/19_nucdata.lua <xsdir_path> [zaid]")
    print("       xsdir_path: path to xsdir file or Nuclear Data directory")
    print("       zaid:       e.g. 92235.80c (default)")
    print("\n19_nucdata: SKIPPED (no xsdir path provided)")
    return
end

local zaid = (arg and arg[2]) or "92235.80c"

-- ---------------------------------------------------------------
-- 1. Load xsdir
-- ---------------------------------------------------------------
print("--- Loading xsdir ---")

local xsdir
-- Try as file first, fall back to directory
local ok, err = pcall(function()
    xsdir = alea.nuc_load_xsdir(xsdir_path)
end)
if not ok then
    ok, err = pcall(function()
        xsdir = alea.nuc_load_xsdir_dir(xsdir_path)
    end)
end
if not ok then
    print("Failed to load xsdir: " .. tostring(err))
    return
end

print(string.format("  Loaded %d entries", xsdir:count()))

-- Look up entry info
local entry = xsdir:find(zaid)
if not entry then
    print("  ZAID " .. zaid .. " not found in xsdir")
    return
end
print(string.format("  Found: %s  AWR=%.4f  kT=%.3e MeV  file=%s",
      entry.zaid, entry.awr, entry.temperature, entry.filename))

-- ---------------------------------------------------------------
-- 2. Load nuclide
-- ---------------------------------------------------------------
print("\n--- Loading nuclide ---")

local nuc = xsdir:load_nuclide(zaid)
print("  " .. tostring(nuc))

local E_min, E_max = nuc:energy_range()
print(string.format("  Energy range: %.3e — %.3e MeV", E_min, E_max))
print(string.format("  Temperature:  %.3e MeV (%.1f K)",
      nuc:temperature(), nuc:temperature() / 8.617333e-11))

-- ---------------------------------------------------------------
-- 3. Cross-section queries
-- ---------------------------------------------------------------
print("\n--- Cross sections at selected energies ---")
print(string.format("  %-12s  %12s  %12s  %12s  %12s",
      "Energy (MeV)", "Total (b)", "Elastic (b)", "Absorp (b)", "Heating"))

local energies = {1e-8, 1e-6, 1e-4, 0.0253e-6, 1e-2, 1.0, 10.0}
for _, E in ipairs(energies) do
    if E >= E_min and E <= E_max then
        print(string.format("  %-12.3e  %12.4f  %12.4f  %12.4f  %12.4e",
              E,
              nuc:xs_total(E),
              nuc:xs_elastic(E),
              nuc:xs_absorption(E),
              nuc:xs_heating(E)))
    end
end

-- Nu-bar (if fissile)
local nu = nuc:nu_bar(1.0)
if nu > 0 then
    print(string.format("\n  nu-bar at 1 MeV:       %.4f", nu))
    print(string.format("  nu-bar at 0.0253 eV:   %.4f", nuc:nu_bar(0.0253e-6)))
end

-- ---------------------------------------------------------------
-- 4. Reaction list
-- ---------------------------------------------------------------
print("\n--- Reactions ---")
local rxns = nuc:reactions()
print(string.format("  %d non-elastic reactions:", #rxns))
for i, r in ipairs(rxns) do
    if i <= 20 then
        print(string.format("    MT %-4d  Q = %+.4f MeV", r.mt, r.q_value))
    elseif i == 21 then
        print(string.format("    ... and %d more", #rxns - 20))
    end
end

-- ---------------------------------------------------------------
-- 5. Build a material
-- ---------------------------------------------------------------
print("\n--- Material (pure " .. zaid .. ") ---")

local mat = alea.nuc_material()
local number_density = 0.048  -- atoms/barn-cm (typical for U metal)
mat:add(nuc, number_density)
print(string.format("  Components: %d", mat:n_components()))

print(string.format("  %-12s  %12s  %12s  %12s",
      "Energy (MeV)", "Sigma_t (1/cm)", "MFP (cm)", "Sigma_a (1/cm)"))

for _, E in ipairs({0.0253e-6, 1e-4, 1.0, 10.0}) do
    if E >= E_min and E <= E_max then
        print(string.format("  %-12.3e  %12.4f  %12.4f  %12.4f",
              E,
              mat:xs_total(E),
              mat:mean_free_path(E),
              mat:xs_absorption(E)))
    end
end

-- ---------------------------------------------------------------
-- 6. Multigroup collapse
-- ---------------------------------------------------------------
print("\n--- 4-group collapse ---")

local bounds = {20.0, 0.1, 1e-3, 1e-6, 1e-11}
local mg = alea.nuc_mg_create(4, bounds)
mg:collapse(nuc)

print(string.format("  %-6s  %12s  %12s  %12s  %12s  %12s  %8s",
      "Group", "E_hi (MeV)", "sigma_t", "sigma_a", "sigma_s", "nu*sigma_f", "chi"))
for g = 1, mg:n_groups() do
    print(string.format("  %-6d  %12.3e  %12.4f  %12.4f  %12.4f  %12.4f  %8.5f",
          g, bounds[g],
          mg:sigma_t(g), mg:sigma_a(g), mg:sigma_s(g),
          mg:nu_sigma_f(g), mg:chi(g)))
end

print("\n  Scattering matrix (g -> g'):")
print(string.format("  %8s  %12s  %12s  %12s  %12s", "", "g'=1", "g'=2", "g'=3", "g'=4"))
for gf = 1, 4 do
    print(string.format("  g=%-5d  %12.4f  %12.4f  %12.4f  %12.4f",
          gf,
          mg:scatter(gf, 1), mg:scatter(gf, 2),
          mg:scatter(gf, 3), mg:scatter(gf, 4)))
end

print("\n19_nucdata: OK")
