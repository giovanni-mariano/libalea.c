# Native Metal ray-group contract and coverage

R4 preparation, 2026-10-31, Apple M1. This increment establishes a ray-expression
reference and measures intersection-group opportunity in actual ordinary renders.
The supported point option remains documented in
[the package guide](NATIVE_GPU_CSG_GROUPS.md). Native ray-group execution is the
next implementation candidate; the diagnostic still executes ordinary bytecode.

## Contract for the native evaluator

A ray-expression result is a tuple: signed truth, end of its validity interval,
boundary primitive, and primitive observed exactly at the cursor. Each leaf solves
its original primitive, selects the first root strictly beyond the cursor, and
probes the same finite interior point. Preserve half-space sense, signed boundary
classification, primitive/leaf occurrence order and the first observable error.

Combine signed truth using the existing operator. Combine validity ends using a
strict minimum: equal distances retain the left operand's boundary. The at-cursor
identity follows the last evaluated matching root. Short-circuiting skips the same
leaves, but already visited leaves still constrain the interval. An intersection
loop that returns only the failed leaf's crossing can therefore be incorrect.
Surrounding union, complement and difference operations retain their structure.

Pure intersections may occupy several intervals. The reference slab with a sphere
removed has two occupied intervals along the test ray. Recognition grants no
convexity certificate. Retain interval stepping, spatial proof selection and
fallback, primitive solvers, hierarchy clipping and surface-normal construction.

Resumable evaluation must retain the group position, reduced signed value,
validity end, boundary and at-cursor identity. Charge work per visited leaf and
bounded transition; an entire long group must not become one uninterruptible work
unit. Preserve snapshot/generation ownership, cancellation and slot refill.
Ordinary and resumable variants require separate complete-render validation.

## Implemented reference and diagnostic

`benchmarks/ray_contract_oracle.h` uses recursive source expressions and analytic
one-dimensional plane/sphere roots. It compares 203 cursor/expression cases in
both instruction formats against the actual ordinary evaluator and one-opcode
resumptions. Cases include a slab, a slab with a spherical hole, mixed union,
difference, repeated leaves, distinct coincident-surface ties and explicit
complement. Checks cover truth, exact
next distance, boundary, at-cursor identity and ordinary visited-leaf count/order.
The resumed checks cover returned expression state; they do not independently
trace every resumed leaf or certify full progressive-render scheduling.

`benchmarks/ray_group_work.metal` instruments actual original/spatial program
ranges selected by ordinary rendering, including deferred-parent validation.
It uses the package's immutable group maps and generic bytecode. Per-pixel
saturating counters count expression calls, instructions, primitives, grouped
instructions/primitives and group entries. Totals reconcile with the existing
ray profiler; saturation or changed raw render pixels fails the diagnostic.
Pixel comparison includes IDs, depth, normals, colors and error status.

The benchmark bridge binds diagnostic buffers only in its separate executable.
Shader hooks expand to nothing in normal builds. No ray execution setting or new
public API is introduced. The existing package group option supplies preparation
metadata; it still accelerates only point evaluation.

## Reproduction and limits

```sh
make -j4 RELEASE=1 USE_METAL=1 \
  BUILD_DIR=build/native-gpu/release BIN_DIR=build/native-gpu/bin \
  build/native-gpu/bin/compute/metal/benchmark_metal_ray_groups
mkdir -p build/native-gpu/ray-groups
python3 benchmarks/native_gpu_compare.py \
  --config benchmarks/metal_ray_groups.example.json \
  --output build/native-gpu/ray-groups/public-check --blocks 1 --timeout 90
```

The public fixture renders 32x24 AA1 with translated fills and mixed CSG. Wide
and compact configurations pass Metal API validation. This is an attribution
run, not a timing experiment. Two alternating timing blocks are required for a
future performance decision.

For local private validation, the executable also accepts `OUTPUT.csv MODEL
CAMERA.txt`. The camera file contains ten finite numbers: eye XYZ, target XYZ,
clip normal XYZ and clip offset; a zero normal disables clipping. These inputs,
model paths, raw per-cell profiles and rendered data belong in ignored build
storage. The private diagnostic is bounded to 4x3 AA1. Preparation is checked
against 15 seconds; the runner enforces 90 seconds per process. The additional
group-map/reference/counter storage is checked against 64 MiB, total device
allocation against 1 GiB, and reported process peak against 4 GiB. These are
boundary checks, not a continuous memory sampler. Allocation failure aborts the
benchmark; ordinary package fallback retains its documented policy.

The independent reference currently covers exactly representable axial
plane/sphere cases. Tangencies, nearly coincident roots, nonfinite errors,
transformed/lattice occurrences, clipping and full first-visible decisions need
additional independent cases before ray-group promotion. Existing renderer tests
remain regression evidence. Private CPU/Metal hit/depth disagreements remain open;
byte-identical diagnostic output preserves a baseline without resolving them.

## Coverage result and decision

One bounded diagnostic pass, serialized on Apple M1 with four CPU threads,
compact code and jump threading off. Private views are 4x3 AA1 samples; the public
scene is 32x24 AA1. The generic ordinary renderer supplies fresh comparison pixels.

| Workload | Instructions within groups | Primitives within groups |
| --- | ---: | ---: |
| Public translated-fill scene | 51.857% | 59.582% |
| Private exterior sample | 73.336% | 94.943% |
| Private cutaway sample | 73.827% | 95.145% |
| Private close-up sample | 71.340% | 93.043% |

All outputs match byte-for-byte, profiler totals reconcile and no counter
saturates. The analytic probe passes 203 cases in each of wide/compact formats.
The affected 75-test Metal compute suite passes with point groups, compact code
and jump threading enabled; thirteen resource tests pass with Metal API validation.
Private context preparation is 9.67–10.26 seconds, maximum reported process peak
is approximately 2.70 GiB and device allocation is approximately 491 MiB. All
processes finish within the 90-second guard. The initial sandbox-only attempt
could not access the device; the successful runs use actual Metal device access.

**Decision: proceed to one bounded ordinary-ray group-evaluator prototype.**
The large executed share justifies this experiment. It does not predict a speedup:
primitive solves remain substantial, and grouping initially preserves their work.
Reusing the existing maps avoids designing another complete geometry format.

First expand the independent contract cases identified above, then replace only
eligible group interpretation with an ordered tuple reduction. Keep generic ray
execution available for comparison/fallback. After ordinary-ray correctness and
performance pass, implement bounded resumptions and measure preview/detail
behavior. The current increment establishes opportunity and a reference; it
contains no native ray-group speedup or transport-accuracy claim.

Private configurations, cameras, profiles, manifests and raw outputs remain under
ignored `build/native-gpu/ray-groups/`. Tracked reports contain aggregate findings
and independently constructed fixtures only.


The subsequent [ordinary-ray prototype evaluation](NATIVE_GPU_RAY_EXECUTION.md)
is complete. Native execution passes the expanded reference and render checks,
but misses the performance/regression gates; package/default and resumable
promotion are deferred. The R4a coverage decision above is historical opportunity
evidence, not a claim that native ray groups have been accepted into the package.
