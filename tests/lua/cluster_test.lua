-- SPDX-FileCopyrightText: 2026 Giovanni MARIANO
-- SPDX-License-Identifier: MPL-2.0

local cluster = require("alea_cluster")
local alea = require("alea")

cluster.initialize()
local context = cluster.create()
assert(context:size() >= 1)
assert(context:rank() >= 0 and context:rank() < context:size())
assert(type(context:backend()) == "string")
assert(context:is_root() == (context:rank() == 0))

local status, message = context:agree(0)
assert(status == 0 and type(message) == "string")

local invalid_read_ok = pcall(function()
    context:read_file(nil)
end)
assert(not invalid_read_ok, "an invalid root path should fail collectively")
assert(context:agree(0) == 0, "the context should recover after an agreed error")

local text = context:read_file(context:is_root() and "README.md" or nil)
assert(text:match("libalea"), "collective file read should broadcast README content")

local system = alea.create()
local sphere = system:sphere(1, 0, 0, 0, 1)
local material = system:material(1)
system:cell{id = 1, region = system:inside(sphere), material = material, density = 1}
system:build_universe_index()
local result = context:estimate_volumes(system, {max_rays = 2000, seed = 17})
assert(#result.volumes == system:volume_path_count())
assert(result.rank_count == context:size())
assert(result.rays_completed > 0)

local rays = context:is_root() and {
    {origin = {0, 0, -2}, direction = {0, 0, 1}},
    {origin = {2, 0, -2}, direction = {0, 0, 1}},
} or nil
local hits = context:raycast_first_segments(system, rays, 4)
if context:is_root() then
    assert(type(hits[1]) == "table" and hits[1].cell_id == 1)
    assert(hits[2] == false, "a missing first concrete segment should be false")
else
    assert(hits == nil)
end

local bad_rays = context:is_root() and {{origin = {0, 0}, direction = {0, 0, 1}}} or nil
local invalid_rays_ok = pcall(function()
    context:raycast_first_segments(system, bad_rays, 4)
end)
assert(not invalid_rays_ok, "invalid root ray data should fail collectively")
assert(context:agree(0) == 0, "the context should remain usable after ray validation")

context:close()
system:destroy()
cluster.finalize()
print("test_cluster: OK")
