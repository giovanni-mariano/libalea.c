-- SPDX-FileCopyrightText: 2026 Giovanni MARIANO
-- SPDX-License-Identifier: MPL-2.0

local first = require("alea")
local second = require("alea")

assert(type(first) == "table", "require('alea') should return the module table")
assert(first == second, "require should cache and return the same module table")
assert(first == _G.alea, "the CLI-compatible alea global should reference the module")
assert(type(first.version()) == "string", "the loadable module should expose the runtime API")

local system = first.create()
assert(tostring(system):match("^System%(") ~= nil, "the module should construct systems")
system:destroy()

print("test_require: OK")
