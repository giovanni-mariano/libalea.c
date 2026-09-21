-- SPDX-FileCopyrightText: 2026 Giovanni MARIANO
--
-- SPDX-License-Identifier: MPL-2.0

local geometry = [[Void transport smoke
c
1 0 -1 imp:n=1 imp:p=1
99 0 1 imp:n=0 imp:p=0

1 s 0 0 0 2
]]

local path = os.tmpname()
local file = assert(io.open(path, "w"))
file:write("directory\n1001.80c 1.0 dummy.ace 0 1 1 1 0 0 2.5301e-8\n")
file:close()
local xsdir = alea.nuc_load_xsdir(path)
os.remove(path)
local system = alea.load_mcnp_string(geometry)

local prepared = alea.source_prepare({
    particle = "neutron",
    space = {type = "point", position = {0, 0, 0}},
    angle = {type = "isotropic"},
    energy = {type = "mono", value = 14.1},
    time = {type = "constant", value = 0.5},
})
assert(not pcall(alea.source_prepare, {
    kind = "box_isotropic", lower = {0, 0, 0}, upper = {1, 1, 1},
    energy = 14.1,
}))
local preview = alea.sample_source(prepared, 4, 17, 6)
assert(#preview.position == 4 and preview.history_id[1] == 6)
assert(preview.energy[4] == 14.1 and preview.time[1] == 0.5)
assert(preview.particle[1] == 0)
for i = 1, 4 do
    assert(preview.position[i][1] == 0 and preview.position[i][2] == 0)
    local d = preview.direction[i]
    assert(math.abs(d[1]*d[1] + d[2]*d[2] + d[3]*d[3] - 1) < 1e-12)
end

local lines = alea.source_prepare({
    space = {type = "point", position = {0, 0, 0}},
    angle = {type = "isotropic"},
    energy = {type = "lines", values = {2.45, 14.1}, weights = {1, 3}},
})
local line_samples = alea.sample_source(lines, 64, 37)
for _, energy in ipairs(line_samples.energy) do
    assert(energy == 2.45 or energy == 14.1)
end
for _, angle in ipairs({
    {type = "cone", direction = {0, 0, 1}, half_angle = math.pi / 3},
    {type = "cosine", direction = {0, 0, 1}},
    {type = "tabulated_mu", direction = {0, 0, 1},
     mu = {0, 1}, pdf = {0, 2}, interpolation = "linear"},
}) do
    local source = alea.source_prepare({
        space = {type = "point", position = {1, 0, 0}},
        angle = angle, energy = 14.1,
    })
    local samples = alea.sample_source(source, 10, 19)
    for _, direction in ipairs(samples.direction) do
        assert(direction[3] >= 0 and direction[3] <= 1)
    end
end
local radial = alea.source_prepare({
    space = {type = "line", start = {1, 0, 0}, ["end"] = {2, 0, 0}},
    angle = {type = "radial", origin = {0, 0, 0}, inward = true},
    energy = 14.1,
})
assert(alea.sample_source(radial, 1, 19).direction[1][1] == -1)
local tokamak = alea.source_prepare({
    space = {type = "tokamak_rz", r_edges = {0, 1, 2}, z_edges = {0, 1, 2},
             emissivity = {{1, 0}, {1, 1}}, phi_min = 0,
             phi_max = math.pi / 2},
    angle = {type = "isotropic"}, energy = 14.1,
})
assert(math.abs(alea.source_integrated_emissivity(tokamak) - 1.75*math.pi) < 1e-12)
local plasma = alea.sample_source(tokamak, 1000, 73)
local outer = 0
for _, point in ipairs(plasma.position) do
    local r2 = point[1]^2 + point[2]^2
    assert(point[1] >= -1e-14 and point[2] >= -1e-14)
    assert(not (r2 < 1 and point[3] >= 1))
    if r2 >= 1 then outer = outer + 1 end
end
assert(math.abs(outer / 1000 - 6 / 7) < 0.05)
assert(not pcall(alea.source_prepare, {
    space = {type = "tokamak_rz", r_edges = {0, 1}, z_edges = {0, 1, 2},
             emissivity = {{1}}},
    angle = {type = "isotropic"}, energy = 14.1,
}))
collectgarbage("collect")
local energy_pdf = alea.source_prepare({
    space = {type = "point", position = {0, 0, 0}},
    angle = {type = "isotropic"},
    energy = {type = "tabulated", values = {1, 3},
              pdf = {1, 1}, interpolation = "linear"},
})
local sampled_energy = alea.sample_source(energy_pdf, 100, 42).energy
for _, value in ipairs(sampled_energy) do
    assert(value >= 1 and value <= 3)
end
local mixture = alea.source_prepare({
    type = "mixture", components = {
        {strength = 1, source = {
            particle = "neutron", space = {type = "point", position = {0, 0, 0}},
            angle = {type = "isotropic"}, energy = 2}},
        {strength = 3, source = {
            particle = "photon", space = {type = "point", position = {1, 0, 0}},
            angle = {type = "isotropic"}, energy = 3}},
    },
})
local mixed = alea.sample_source(mixture, 1000, 17)
local photons = 0
for i, particle in ipairs(mixed.particle) do
    if particle == 1 then
        photons = photons + 1
        assert(mixed.position[i][1] == 1)
    end
    assert(mixed.weight[i] == 1)
end
assert(math.abs(photons / 1000 - 0.75) < 0.05)
local mixed_result = alea.transport_run(system, xsdir, {
    histories = 8, source = mixture,
    tallies = {{score = "track_length", particle = "all"}},
})
assert(mixed_result.leaked == 8)
local mesh = alea.source_prepare({
    space = {type = "cartesian_mesh", x_edges = {0, 1, 3},
             y_edges = {0, 1}, z_edges = {0, 1, 2},
             values = {{{1, 0}}, {{1, 1}}}, value_mode = "density"},
    angle = {type = "isotropic"}, energy = 14.1,
})
assert(math.abs(alea.source_integrated_emissivity(mesh) - 5) < 1e-12)
local mesh_points = alea.sample_source(mesh, 1000, 63).position
local mesh_outer = 0
for _, point in ipairs(mesh_points) do
    assert(not (point[1] < 1 and point[3] >= 1))
    if point[1] >= 1 then mesh_outer = mesh_outer + 1 end
end
assert(math.abs(mesh_outer / 1000 - 0.8) < 0.05)
local mesh_strength = alea.source_prepare({
    space = {type = "cartesian_mesh", x_edges = {0, 1, 3},
             y_edges = {0, 1}, z_edges = {0, 1, 2},
             values = {{{1, 0}}, {{1, 1}}}, value_mode = "strength"},
    angle = {type = "isotropic"}, energy = 14.1,
})
assert(math.abs(alea.source_integrated_emissivity(mesh_strength) - 3) < 1e-12)
assert(not pcall(alea.source_prepare, {
    space = {type = "cartesian_mesh", x_edges = {0, 1},
             y_edges = {0, 1}, z_edges = {0, 1, 2},
             values = {{{1}}}, value_mode = "density"},
    angle = {type = "isotropic"}, energy = 14.1,
}))
collectgarbage("collect")
local replay = alea.sample_source(prepared, 1, 17, 8)
assert(replay.direction[1][1] == preview.direction[3][1])
local prepared_result = alea.transport_run(system, xsdir, {
    histories = 4, seed = 17, history_offset = 6, source = prepared,
    tallies = {{score = "track_length"}},
})
assert(prepared_result.leaked == 4)

for _, space in ipairs({
    {type = "line", start = {0, 0, 0}, ["end"] = {2, 0, 0}},
    {type = "sphere", center = {0, 0, 0}, inner_radius = 0.5, outer_radius = 1},
    {type = "cylinder", base = {0, 0, 0}, axis = {0, 0, 2},
     inner_radius = 0.5, outer_radius = 1},
}) do
    local source = alea.source_prepare({
        space = space, angle = {type = "isotropic"}, energy = 14.1,
    })
    local samples = alea.sample_source(source, 32, 13)
    for _, point in ipairs(samples.position) do
        if space.type == "line" then
            assert(point[1] >= 0 and point[1] <= 2)
            assert(point[2] == 0 and point[3] == 0)
        elseif space.type == "sphere" then
            local radius = math.sqrt(point[1]^2 + point[2]^2 + point[3]^2)
            assert(radius >= 0.5 - 1e-12 and radius <= 1 + 1e-12)
        else
            local radius = math.sqrt(point[1]^2 + point[2]^2)
            assert(radius >= 0.5 - 1e-12 and radius <= 1 + 1e-12)
            assert(point[3] >= 0 and point[3] <= 2)
        end
    end
end

for _, particle in ipairs({"neutron", "photon"}) do
    local config = {
        histories = 4, seed = 7,
        source = {particle = particle, energy = 2,
                  space = {type = "point", position = {0, 0, 0}},
                  angle = {type = "monodirectional", direction = {1, 0, 0}}},
        tallies = {
            {score = "track_length", particle = particle},
            {score = "track_length", domain = "mesh", particle = particle,
             lower = {-2, -2, -2}, upper = {2, 2, 2},
             dimensions = {2, 1, 1}, energy_edges = {0, 3}},
        },
    }
    local result = alea.transport_run(system, xsdir, config)
    assert(result.histories == 4 and result.leaked == 4)
    assert(math.abs(result.tallies[1].mean[1] - 2) < 1e-10)
    assert(result.tallies[1].standard_error[1] == 0)
    assert(math.abs(result.tallies[2].mean[1]) < 1e-10)
    assert(math.abs(result.tallies[2].mean[2] - 2) < 1e-10)
    local replay = alea.transport_run(system, xsdir, config)
    assert(replay.tallies[1].sum[1] == result.tallies[1].sum[1])
end

local box = alea.transport_run(system, xsdir, {
    histories = 8,
    source = {particle = "neutron", energy = 14.1,
              space = {type = "box", lower = {0, 0, 0}, upper = {0, 0, 0}},
              angle = {type = "isotropic"}},
    tallies = {{score = "track_length"}},
})
assert(box.leaked == 8)
local ok = pcall(function()
    alea.transport_run(system, xsdir, {
        source = {energy = 1,
                  space = {type = "point", position = {0, 0, 0}},
                  angle = {type = "monodirectional", direction = {1, 0, 0}}},
        tallies = {{score = "invalid"}},
    })
end)
assert(not ok)
