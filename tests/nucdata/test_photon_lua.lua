-- SPDX-FileCopyrightText: 2026 Giovanni MARIANO
--
-- SPDX-License-Identifier: MPL-2.0

-- Real-data parity checks for the Lua nuclear-data sampling surface.

local xs07 = alea.nuc_load_xsdir("njoy-test07/xsdir")
local u235 = xs07:load_nuclide("92235.00c")
local channels07 = u235:photon_productions()
assert(#channels07 == 33, "unexpected U-235 photon-production count")
assert(channels07[1].mt == 4001)
assert(channels07[1].parent_mt == 4)
assert(channels07[1].mf == 12)
assert(channels07[1].energy_law == 2)

local line_a = u235:sample_photon(4001, 1.0, 1234, 9, 2, 7)
local line_b = u235:sample_photon(4001, 1.0, 1234, 9, 2, 7)
assert(math.abs(line_a.energy - 0.771) < 1e-14)
assert(line_a.energy == line_b.energy and line_a.mu == line_b.mu,
       "identical RNG addresses must replay")
assert(#line_a.direction == 3)
local norm = line_a.direction[1]^2 + line_a.direction[2]^2 +
             line_a.direction[3]^2
assert(math.abs(norm - 1.0) < 1e-14, "photon direction must be normalized")
assert(line_a.mu == line_a.direction[3])
assert(line_a.weight == 1.0 and line_a.time == 0.0)

local xs08 = alea.nuc_load_xsdir("njoy-test08/xsdir")
local ni61 = xs08:load_nuclide("28061.00c")
local channels08 = ni61:photon_productions()
assert(#channels08 == 36, "unexpected Ni-61 photon-production count")

local continuous_a = ni61:sample_photon(16001, 13.0, 1234, 9, 2, 7)
local continuous_b = ni61:sample_photon(16001, 13.0, 1234, 9, 2, 8)
assert(continuous_a.energy >= 0.0 and continuous_a.energy <= 13.0)
assert(continuous_a.mu >= -1.0 and continuous_a.mu <= 1.0)
assert(continuous_a.energy ~= continuous_b.energy or
       continuous_a.mu ~= continuous_b.mu,
       "different event addresses should select different samples")

local outgoing_a = ni61:sample_reaction_energy(51, 13.0, 8128, 61, 0, 3)
local outgoing_b = ni61:sample_reaction_energy(51, 13.0, 8128, 61, 0, 3)
assert(outgoing_a == outgoing_b, "reaction-energy sample must replay")
assert(outgoing_a >= 0.0 and outgoing_a <= 13.0)

local ok = pcall(function()
    ni61:sample_photon(999999, 13.0, 1, 0, 0, 0)
end)
assert(not ok, "unknown photon-production MT must fail")

print("test_photon_lua: OK")
