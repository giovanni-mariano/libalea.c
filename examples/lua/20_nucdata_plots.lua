-- 20_nucdata_plots.lua: Nuclear data SVG plotting
--
-- Standalone: no (requires ACE data files and xsdir)
-- Usage: bin/alea examples/lua/20_nucdata_plots.lua <xsdir_path> [zaid] [outdir]
--
-- Generates SVG plots:
--   1. Cross sections vs energy (log-log)
--   2. Per-reaction breakdown
--   3. Nu-bar vs energy (if fissile)
--   4. Doppler broadening comparison
--   5. Multigroup bar chart

print("=== Nuclear Data Plotting ===\n")

local xsdir_path = arg and arg[1]
if not xsdir_path then
    print("Usage: bin/alea examples/lua/20_nucdata_plots.lua <xsdir_path> [zaid] [outdir]")
    print("\n20_nucdata_plots: SKIPPED (no xsdir path provided)")
    return
end

local zaid = (arg and arg[2]) or "92235.80c"
local outdir = (arg and arg[3]) or "."

-- Load data
local xsdir
local ok, err = pcall(function() xsdir = alea.nuc_load_xsdir(xsdir_path) end)
if not ok then
    ok, err = pcall(function() xsdir = alea.nuc_load_xsdir_dir(xsdir_path) end)
end
if not ok then
    print("Failed to load xsdir: " .. tostring(err))
    return
end

local nuc = xsdir:load_nuclide(zaid)
print("Loaded: " .. tostring(nuc))

local E_min, E_max = nuc:energy_range()

-- Helper: sample a cross-section function over the energy grid
local function sample_xs(fn, n)
    n = n or 500
    local xs, es = {}, {}
    local log_min = math.log(E_min) / math.log(10)
    local log_max = math.log(E_max) / math.log(10)
    for i = 1, n do
        local log_E = log_min + (log_max - log_min) * (i - 1) / (n - 1)
        local E = 10 ^ log_E
        es[i] = E
        xs[i] = fn(E)
    end
    return es, xs
end

-- ---------------------------------------------------------------
-- Plot 1: Main cross sections
-- ---------------------------------------------------------------
print("\n--- Plot 1: Cross sections ---")

local plot = alea.svg_plot()
plot:set_title(nuc:zaid() .. " Cross Sections")
plot:set_x_label("Energy (MeV)")
plot:set_y_label("Cross section (barn)")
plot:set_scale("log", "log")

local E, sig_t = sample_xs(function(e) return nuc:xs_total(e) end)
local _, sig_el = sample_xs(function(e) return nuc:xs_elastic(e) end)
local _, sig_abs = sample_xs(function(e) return nuc:xs_absorption(e) end)

plot:add(E, sig_t,   "Total",      "#2266cc")
plot:add(E, sig_el,  "Elastic",    "#22aa44")
plot:add(E, sig_abs, "Absorption", "#cc2222")

local fname = outdir .. "/" .. nuc:zaid():gsub("%.", "_") .. "_xs.svg"
plot:write(fname)
print("  Wrote: " .. fname)

-- ---------------------------------------------------------------
-- Plot 2: Per-reaction breakdown
-- ---------------------------------------------------------------
print("\n--- Plot 2: Reaction cross sections ---")

local rxns = nuc:reactions()
local colors = {"#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4",
                "#42d4f4", "#f032e6", "#bfef45", "#fabebe", "#469990"}

-- Pick the top reactions by peak cross section
local rxn_peaks = {}
for _, r in ipairs(rxns) do
    local peak = 0
    for _, e in ipairs({1e-8, 1e-6, 1e-4, 1e-2, 1.0, 5.0}) do
        if e >= E_min and e <= E_max then
            local v = nuc:xs_reaction(r.mt, e)
            if v > peak then peak = v end
        end
    end
    table.insert(rxn_peaks, {mt = r.mt, peak = peak})
end
table.sort(rxn_peaks, function(a, b) return a.peak > b.peak end)

local plot2 = alea.svg_plot()
plot2:set_title(nuc:zaid() .. " Reaction Cross Sections")
plot2:set_x_label("Energy (MeV)")
plot2:set_y_label("Cross section (barn)")
plot2:set_scale("log", "log")

local n_show = math.min(#rxn_peaks, 8)
for i = 1, n_show do
    local mt = rxn_peaks[i].mt
    local _, sig = sample_xs(function(e) return nuc:xs_reaction(mt, e) end)
    plot2:add(E, sig, "MT " .. mt, colors[i] or "#888888")
end

local fname2 = outdir .. "/" .. nuc:zaid():gsub("%.", "_") .. "_reactions.svg"
plot2:write(fname2)
print("  Wrote: " .. fname2 .. " (" .. n_show .. " reactions)")

-- ---------------------------------------------------------------
-- Plot 3: Nu-bar (if fissile)
-- ---------------------------------------------------------------
local nu_thermal = nuc:nu_bar(0.0253e-6)
if nu_thermal > 0 then
    print("\n--- Plot 3: Nu-bar ---")

    local plot3 = alea.svg_plot()
    plot3:set_title(nuc:zaid() .. " Average Neutrons per Fission")
    plot3:set_x_label("Energy (MeV)")
    plot3:set_y_label("nu-bar")
    plot3:set_scale("log", "linear")

    local E_nu, nu_vals = sample_xs(function(e) return nuc:nu_bar(e) end)
    plot3:add(E_nu, nu_vals, "nu-bar(E)", "#2266cc")

    local fname3 = outdir .. "/" .. nuc:zaid():gsub("%.", "_") .. "_nubar.svg"
    plot3:write(fname3)
    print("  Wrote: " .. fname3)
end

-- ---------------------------------------------------------------
-- Plot 4: Doppler broadening comparison
-- ---------------------------------------------------------------
print("\n--- Plot 4: Doppler broadening ---")

-- Sample original cross section in the resolved resonance region
local E_lo_doppler = math.max(E_min, 1e-4)  -- 0.1 eV
local E_hi_doppler = math.min(E_max, 1e-1)  -- 100 keV
local n_pts = 800

local function sample_range(fn, elo, ehi, n)
    local es, ys = {}, {}
    local log_lo = math.log(elo) / math.log(10)
    local log_hi = math.log(ehi) / math.log(10)
    for i = 1, n do
        local log_E = log_lo + (log_hi - log_lo) * (i - 1) / (n - 1)
        local e = 10 ^ log_E
        es[i] = e
        ys[i] = fn(e)
    end
    return es, ys
end

-- Original (cold) cross section
local kT_cold = nuc:temperature()
local E_dop, sig_cold = sample_range(
    function(e) return nuc:xs_total(e) end,
    E_lo_doppler, E_hi_doppler, n_pts)

-- Broaden to 600 K (returns a new independent nuclide)
local kT_600K = 600.0 * 8.617333e-11  -- 600 K in MeV
local nuc_hot, berr
local bok, result = pcall(function() return nuc:broadened(kT_600K) end)
if bok then nuc_hot = result else berr = result end

if bok then
    local _, sig_hot = sample_range(
        function(e) return nuc_hot:xs_total(e) end,
        E_lo_doppler, E_hi_doppler, n_pts)

    local plot4 = alea.svg_plot()
    plot4:set_title(nuc:zaid() .. " Doppler Broadening")
    plot4:set_x_label("Energy (MeV)")
    plot4:set_y_label("Total cross section (barn)")
    plot4:set_scale("log", "log")

    local T_cold_K = kT_cold / 8.617333e-11
    plot4:add(E_dop, sig_cold,
        string.format("%.0f K", T_cold_K), "#2266cc")
    plot4:add(E_dop, sig_hot, "600 K", "#cc2222")

    local fname4 = outdir .. "/" .. nuc:zaid():gsub("%.", "_") .. "_doppler.svg"
    plot4:write(fname4)
    print("  Wrote: " .. fname4)
else
    print("  Skipped (broaden failed: " .. tostring(berr) .. ")")
end

-- ---------------------------------------------------------------
-- Plot 5: Multigroup constants
-- ---------------------------------------------------------------
print("\n--- Plot 5: Multigroup constants ---")

-- CASMO-like 8-group structure (descending boundaries in MeV)
local mg_bounds = {20.0, 2.231, 0.8208, 5.53e-3, 4.0e-6, 6.25e-7, 1.3e-7, 5.0e-9, 1e-11}
local n_groups = #mg_bounds - 1

local mg = alea.nuc_mg_create(n_groups, mg_bounds)
mg:collapse(nuc)

-- Build bar-chart-style step data for group constants
local bar_E, bar_st, bar_sa, bar_ss = {}, {}, {}, {}
local idx = 1
for g = 1, n_groups do
    -- Left edge
    bar_E[idx] = mg_bounds[g]
    bar_st[idx] = mg:sigma_t(g)
    bar_sa[idx] = mg:sigma_a(g)
    bar_ss[idx] = mg:sigma_s(g)
    idx = idx + 1
    -- Right edge (same height)
    bar_E[idx] = mg_bounds[g + 1]
    bar_st[idx] = mg:sigma_t(g)
    bar_sa[idx] = mg:sigma_a(g)
    bar_ss[idx] = mg:sigma_s(g)
    idx = idx + 1
end

local plot5 = alea.svg_plot()
plot5:set_title(nuc:zaid() .. " " .. n_groups .. "-Group Cross Sections")
plot5:set_x_label("Energy (MeV)")
plot5:set_y_label("Group-averaged cross section (barn)")
plot5:set_scale("log", "log")

-- Add pointwise in background
plot5:add(E, sig_t, "Pointwise total", "#cccccc", {width = 0.8})
plot5:add(bar_E, bar_st, "Group total",      "#2266cc", {width = 2.0})
plot5:add(bar_E, bar_sa, "Group absorption",  "#cc2222", {width = 2.0})
plot5:add(bar_E, bar_ss, "Group elastic",     "#22aa44", {width = 2.0})

local fname5 = outdir .. "/" .. nuc:zaid():gsub("%.", "_") .. "_multigroup.svg"
plot5:write(fname5)
print("  Wrote: " .. fname5)

print("\n20_nucdata_plots: OK")
