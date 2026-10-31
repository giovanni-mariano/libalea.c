# Packaged Metal CSG group evaluation

The compute package provides ordered signed-leaf evaluation for pure intersection
subtrees in Metal point queries. It is available through the regular
`alea_compute_classify()` and `alea_compute_select()` APIs and the normal
`mc_plotter` tool. No benchmark code or separate execution metallib is required.

## Build and enable

```sh
make compute USE_METAL=1 RELEASE=1
make -C tools ../bin/mc_plotter ../bin/render3d USE_METAL=1 USE_SDL=1 RELEASE=1

ALEA_METAL_CSG_GROUPS=1 bin/mc_plotter model.mcnp Z 0 -10 10 -10 10 512 slice.bmp --backend=metal
```

Set `ALEA_METAL_CSG_GROUPS=1` before creating a Metal context or launching the
tool. The same setting applies to all-owner and ordered-selection queries.
Its default is `0`; values other than `0` and `1` produce an invalid-argument
error. The option is supported with wide/compact instructions, jump threading,
point spatial programs and candidate filtering.

It is off by default because preparation and additional storage can improve large
CSG workloads while short expressions can be slower. The previous controlled
[measurements](NATIVE_GPU_GROUP_EXECUTION.md) describe those tradeoffs. Packaging
does not assert a universal speedup or change the declared FP32 precision contract.

## Ownership and fallback

Each context owns its group maps and signed references. Preparation reads the
completed immutable geometry after proof construction, before jump threading and
compact encoding. Source models may subsequently change or be freed. Settings
are captured at context creation; changing the environment does not alter existing
contexts. Separate contexts have independent buffers and pipeline constants.

The host recognizes maximal pure intersection subtrees without changing operand
order or Boolean structure. A direct native loop replaces their stack operations;
the surrounding expression retains its generic bytecode. Signed-zero boundary
rules, repeated leaf occurrences, first nonfinite errors and short-circuit leaf
visitation are preserved. Exact code-range aliases share descriptors; partial
aliases keep their own program maps.

Group preparation and retained group data have a 64-MiB requested-allocation
budget. Construction accounts for host plan arrays and the device map together,
then releases host staging before allocating the leaf pool. Device buffer limits
are checked. If no eligible groups exist, optional allocation fails, or the data
cannot fit, the context retains generic Metal evaluation. It does not change
backend or precision. Interruption and malformed code remain errors, and partial
group buffers are released before returning.

The standard context destructor releases both buffers, including on later
pipeline/upload failures. There is no process-global group state or capture hook
in the package. Generic bytecode remains available for point fallback, host proofs
and the ray evaluator.

## Diagnostics

An enabled context reports aggregate group count and storage bytes at preparation.
A context using generic evaluation reports that fallback instead. These messages
do not include cell IDs, coordinates or model paths.

Point and selection profiles include:

- `csg_groups`: whether this context actually uses group evaluation.
- `instruction_unit`: `bytecode` or `bytecode_and_group_ops`.

In grouped mode, instruction counters count remaining bytecode plus a group
dispatch and visited list entries. They are logical operation counts, not GPU
machine instructions, and cannot be compared directly to bytecode-only totals.
Primitive-evaluation counts retain their original meaning. Group buffer lengths
are included in selection profiles' context storage total.

## Validation of the integration

The packaged shader, not a separate prototype implementation, passes 1,728
trace/value/error cases in wide and compact formats. They cover mixed operators,
signed zeros, nonfinite leaves, exact visitation, lookup misses and partial code
aliases, with an independent CPU evaluator.

API tests construct actual intersection groups in flat, transformed/nested-fill,
finite/repeating-lattice, overlap, void and unbounded models. They compare enabled
and generic Metal contexts with the CPU API, cross the 4,096-point submission
boundary, free the source, change the environment and reuse both query APIs.
Compact mode is also tested with jump threading. No benchmark callbacks supply
the enabled execution path.

Allocation-failure tests exercise map and reference buffer failures, cleanup,
generic fallback and a successful retry. Tests also cover interruption, invalid
settings, profiler units, empty plans, budget exhaustion and source preservation.
The 75-test Metal suite passes in default, wide grouped and compact/threaded
grouped configurations. The CPU-only suite passes 16 tests, with ten builder
tests and thirteen Metal resource tests also passing. Existing render/resume
tests continue to pass with the point option enabled.

A bounded private-model check through the public APIs returns identical all-owner
and selected output buffers with groups disabled/enabled, after the source model
is freed. The rebuilt normal plotter also produces byte-identical public-fixture
images in both ownership modes with the option disabled/enabled. Private models,
views and raw outputs remain under ignored build storage.

The ray evaluator, transport navigation and hardware-gated CUDA work retain their
separate implementations and validation requirements. This option accelerates
point queries; it does not install a ray group evaluator.


Ordinary rays now have an independent supported
[ALEA_METAL_RAY_GROUPS option](NATIVE_GPU_RAY_GROUPS_PACKAGE.md). Both options
share one bounded prepared data set when enabled together. The point option
alone retains its existing scope.
