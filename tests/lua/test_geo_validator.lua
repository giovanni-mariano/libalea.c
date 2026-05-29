-- test_geo_validator.lua: Geometry validator Lua binding tests

local function build_clean_system()
    local sys = alea.create()
    sys:set_spatial_mode("hier")
    local s1 = sys:sphere(1, 0, 0, 0, 5)
    local m1 = sys:material(1)
    sys:cell{id = 1, region = sys:inside(s1), material = m1, density = 10.0}
    sys:cell{id = 2, region = sys:outside(s1), material = -1, density = 0.0}
    return sys
end

local clean = build_clean_system()
local clean_ray = clean:validate_geometry_ray(
    -10, 0, 0, 1, 0, 0, 25,
    {hier = true, max_errors = 8}
)
assert(clean_ray:error_count() == 0, "clean adjacent sphere should validate")
local clean_stats = clean_ray:stats()
assert(clean_stats.crossings_checked >= 2, "clean ray should check crossings")

local clean_all = clean:validate_geometry{
    hier = true,
    ray_count = 4,
    max_errors = 8
}
assert(clean_all:error_count() == 0, "clean whole-geometry validation should pass")
clean:destroy()

local undefined = alea.create()
undefined:set_spatial_mode("hier")
local sphere = undefined:sphere(1, 0, 0, 0, 1)
local mat = undefined:material(1)
undefined:cell{id = 1, region = undefined:inside(sphere), material = mat, density = 1.0}

local bad = undefined:validate_geometry_ray(
    0, 0, 0, 1, 0, 0, 3,
    {hier = true, max_errors = 8, sample_offset = 0.01}
)
assert(bad:error_count() > 0, "undefined exterior should produce diagnostics")

local first = bad:error(1)
assert(first.type ~= nil, "error should expose a type")
assert(first.crossing_point ~= nil and #first.crossing_point == 3,
       "error should expose crossing point")

local errors = bad:errors()
assert(#errors == bad:error_count(), "errors() should match error_count()")

local summary = bad:summary()
assert(summary.total == bad:error_count(), "summary total should match error_count()")
assert((summary.undefined_after_crossing or 0) > 0,
       "summary should include undefined-region diagnostics")

undefined:destroy()
print("test_geo_validator: OK")
