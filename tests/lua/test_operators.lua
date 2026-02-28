-- test_operators.lua: CSG boolean operators

local sys = alea.create()

local s1 = sys:sphere(0, 0, 0, 0, 10)
local s2 = sys:sphere(0, 5, 0, 0, 8)
local s3 = sys:sphere(0, -5, 0, 0, 8)

local a = sys:inside(s1)
local b = sys:inside(s2)
local c = sys:inside(s3)

-- Union: a + b
local u = a + b
assert(u, "union works")
assert(tostring(u):find("UNION"), "union node type")

-- Intersection: a * b
local i = a * b
assert(i, "intersection works")
assert(tostring(i):find("INTERSECTION"), "intersection node type")

-- Difference: a - b
local d = a - b
assert(d, "difference works")
assert(tostring(d):find("DIFFERENCE"), "difference node type")

-- Complement: ~a
local comp = ~a
assert(comp, "complement works")
assert(tostring(comp):find("COMPLEMENT"), "complement node type")

-- Chained: (a * b) + c
local chain = (a * b) + c
assert(chain, "chaining works")

-- union_n / intersection_n
local un = alea.union_n(sys, {a, b, c})
assert(un, "union_n works")

local in_ = alea.intersection_n(sys, {a, b, c})
assert(in_, "intersection_n works")

-- Create cell with complex region
local m1 = sys:material(1)
local idx = sys:cell{id=1, region=chain, material=m1, density=1.0}
assert(idx >= 0, "cell from chained ops")

-- Node inspection
local root = sys:cell_info(idx).root
assert(root, "cell has root node")
local op = sys:node_operation(root)
assert(type(op) == "string", "node_operation returns string")

print("test_operators: OK")
