-- test_geo_validator.lua: Geometry validator Lua binding tests

local function build_clean_system()
    local sys = alea.create()
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

local gaps = undefined:validate_geometry_ray(
    -2, 0, 0, 1, 0, 0, 4,
    {allow_exterior_void = true,
     validation_bounds = {-1.5, 1.5, -0.5, 0.5, -0.5, 0.5}}
)
assert(gaps:error_count() == 2,
       "an explicit validation domain should report both interior gaps")
assert((gaps:summary().interior_gap or 0) == 2,
       "domain gaps should use the interior_gap finding kind")

undefined:destroy()

-- Surface/slice-driven validation (Phase 5)
local nested = alea.create()
local big = nested:sphere(10, 0, 0, 0, 3)
local small = nested:sphere(20, 0, 0, 0, 1)
local mb = nested:material(1)
local ms = nested:material(2)
nested:cell{id = 1, region = nested:inside(big), material = mb, density = 1.0}
nested:cell{id = 2, region = nested:inside(small), material = ms, density = 1.0}

local sview = alea.slice_view_axis(2, 0.0, -4, 4, -4, 4)
local directional_cache = nested:directional_trace_cache(sview, 16, 16)
local directional = nested:validate_ray_slice(
    sview, 16, {include_agreements = true}, directional_cache)
assert(directional.reused_trace_mask == 3,
       "directional cache should reuse U forward/reverse traces")
assert(directional.executed_trace_mask == 0,
       "matching directional cache should avoid fresh ownership traces")
local mismatched_view = alea.slice_view_axis(2, 0.0, -3, 4, -4, 4)
local ok = pcall(function()
    nested:validate_ray_slice(mismatched_view, 16, nil, directional_cache)
end)
assert(not ok, "a cache must reject a different slice view")
local scurves = nested:get_slice_curves(sview)
local sres = nested:validate_geometry_slice(sview, scurves, {allow_exterior_void = true})
local ssum = sres:summary()
assert((ssum.overlap_after_crossing or 0) > 0,
       "slice validation should detect the hidden nested overlap")
local sstats = sres:stats()
assert(sstats.crossings_checked > 0, "slice validation should sample boundaries")
local serr = sres:error(1)
assert(serr.curve_index ~= nil, "slice events should expose curve index")
nested:sphere(99, 10, 0, 0, 1)
ok = pcall(function()
    nested:validate_ray_slice(sview, 16, nil, directional_cache)
end)
assert(not ok, "a cache must reject geometry changed after construction")
nested:destroy()

print("test_geo_validator: OK")
