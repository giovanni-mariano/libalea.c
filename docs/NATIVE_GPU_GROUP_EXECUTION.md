# Native Metal intersection-group execution experiment

Package integration is now available as the supported opt-in
`ALEA_METAL_CSG_GROUPS=1` option. See [the package guide](NATIVE_GPU_CSG_GROUPS.md).
The measurements and decision below describe the preceding benchmark experiment;
they remain the reason the package option is off by default.

Completed 2026-10-31 on Apple M1. Native point evaluation now has a bounded
benchmark prototype that executes recognized pure intersection subtrees as ordered
signed-leaf lists inside the existing general CSG interpreter. The
[preceding audit](NATIVE_GPU_INTERSECTION_GROUPS.md) selected this experiment.

**Decision: retain the working benchmark experiment; do not promote a library
default.** Private throughput gates pass, including improvement beyond jump
threading. Public timing regressions and variable preparation overhead prevent
promotion. The normal library and viewers do not use this prototype.

## Implementation and contract

The benchmark captures immutable group descriptors and leaf lists after ordinary
geometry/spatial preparation and proof construction, before optional jump
threading and compact encoding. It leaves all original bytecode available. The
native point evaluator substitutes a direct ordered loop at an admitted group
start, pushes one reduced signed value, then continues at the original group end.
All other instructions use the existing interpreter.

The loop preserves signed primitive references, duplicate occurrences, leaf order,
finite checks, early rejection and surrounding complement/difference semantics.
It computes a signed max reduction, rather than Boolean negation, preserving
surface-zero behavior. It does not eliminate primitive evaluations or change
surface coefficients, coordinate transforms, hierarchy/lattice traversal,
ownership order, overlap counts or spatial certificates. Rays keep their existing
programs and evaluators.

Program descriptors are keyed by exact `(root, count)` ranges, so partial aliases
retain their own groups. The first layout used binary search. The one measured
revision adds a bounded open-addressed hash table and bypasses lookup for programs
shorter than four instructions, which cannot contain an eligible group. Hash
lookup failure falls back to generic interpretation. There is no expression
expansion, full second instruction buffer or new GPU framework.

The implementation lives in `benchmarks/metal_group_execute.m` and
`benchmarks/group_execute.metal`, with conditional hooks in the existing point
evaluator. Normal shader builds compile the hooks out. The experiment owns one
context per process; its capture hook and group buffers are benchmark state, not a
new public backend or a thread-safe library lifecycle implementation.

`ALEA_BENCH_GROUP_EXECUTE=0|1` selects generic/grouped execution in this benchmark;
the default is 0. `ALEA_BENCH_GROUP_HASH=0|1` selects the lookup revision, default 0.
Both settings are strictly parsed. Disabled execution builds no group/list
buffers. Preparation failures clean up unpublished buffers. An unavailable initial
optional map or an exceeded allocation budget retains generic evaluation and
reports fallback; later allocation/validation failures abort context construction.
Interruption remains an error. These settings do not enable the feature in an
ordinary library application.

## Measurements

Isolated release build, serialized baseline/candidate processes in forward and
reverse order, two blocks. Each process performs one initial call and five warm
calls on the same 16,384-point private slice, with 4,096 points per submission and
four CPU threads. Compact code is fixed. Profiling is disabled for all acceptance
timings. The harness includes conversion, submissions, waits and output checks;
these are complete batched point operations, not viewer or image-rendering times.

The following are ranges of per-process warm medians. Reductions compare matched
baseline/candidate processes within each ordering block:

| Operation | Threading | Generic seconds | Grouped seconds | Time reduction |
| --- | --- | ---: | ---: | ---: |
| All-owner classification | Off | 2.623–2.631 | 1.763–1.773 | 32.6–32.8% |
| All-owner classification | On | 1.767–1.771 | 1.375–1.379 | 22.1–22.2% |
| Ordered selection | Off | 1.176–1.245 | 0.785–0.813 | 30.9–37.0% |
| Ordered selection | On | 0.766–0.777 | 0.563–0.574 | 26.1–26.5% |

The threading-on rows measure incremental value beyond the existing optimization;
they do not attribute its earlier gain to grouping. Initial binary-layout runs
showed about 25% and 13% all-owner gains with threading off/on respectively,
motivating the bounded lookup revision. Only the final layout received the full
alternating comparison above.

Every timed output matches a fresh ordinary GPU query byte-for-byte. Retained
outputs also match across grouped/generic and threaded/unthreaded configurations
for each API. This establishes preservation of existing GPU results; it does not
close independent private CPU/GPU accuracy discrepancies.

## Why it remains experimental

The public simple/mixed/long-chain/selection cases pass correctness checks, but
their millisecond timings vary. One bounded repeat collects ten warm observations
per process in two rounds on the same context, again in both ordering blocks.
The small mixed-expression case rises from 2.273 to 2.4525 ms and from 1.857 to
2.090 ms, approximately 7.9% and 12.5%. These exceed the declared 5% regression
gate. Simple and selected results are order-sensitive. Long-chain cases improve,
but do not waive the short-case gate. These observations do not establish a
hardware stall or occupancy diagnosis.

Total snapshot preparation also varies. One all-owner pair rises from 9.409 to
10.697 s (13.7%); one selected threaded pair rises from 9.716 to 12.105 s (24.6%).
Other pairs have smaller or negative changes. The isolated group preparation
portion is approximately 0.39–0.81 s in the final timing runs. This is not evidence
that grouping intrinsically causes the entire paired increase, but the 10% total
preparation gate is unresolved and cannot be reported as passed.

The two planned layouts have been measured and the experiment ends with this
decision. Do not add a third lookup mechanism to rescue the current promotion.
A later integration proposal should examine direct handles from selected cell/
spatial descriptors and cost-based admission of short groups, with fresh gates
and comparisons. Broader native feature/lifecycle coverage is also required
before library integration.

## Correctness, diagnostics and resources

The final native probe compares generic, binary-map and hash-map execution over
1,728 cases: six mixed operator trees plus single-leaf and nested-subtree aliases,
with all combinations of negative, positive, signed-zero, NaN and infinite leaf
inputs. Result bits and exact visited signed-leaf sequences match, including
first nonfinite failures. An independent CPU evaluator checks classification,
errors and visitation. Wide/compact modes pass Metal API validation.

Public analytic sphere/hierarchy ownership checks pass. The fixtures exercise
transformed fills, terminal void, mixed unions and long intersection groups.
Invalid experimental settings, preparation interruption and a successful retry
are checked. Existing ordinary compute tests pass in both formats (74 each),
nine builder tests and ten Metal resource tests pass, and the prior audit
executables remain compatible. These ordinary-library regressions are not proof
that the new driver covers the entire native geometry feature matrix. Direct
opt-in tests for all supported lattice/unbounded/fill cases remain integration
work.

A separate 512-point private diagnostic returns exactly **16,664,952 logical
primitive evaluations** for both generic and grouped execution, with identical
results and no saturation. This is a fixed subset, not a full-slice counter claim.
An initial full-slice profiling attempt hit the 90-second process guard after its
warm timings; the runner stopped it and the failed evidence remains local.
Acceptance timings were then collected in fresh processes without profiling.
Group dispatches/list iterations and bytecode instructions are different units;
their instruction counters are not used as a speedup claim.

The final private hash layout uses 34,614,428 bytes of maps/descriptors and
32,452,304 bytes of signed references: 67,066,732 bytes, approximately 63.96 MiB,
within the 64-MiB extra-allocation guard. This is close to the limit. Construction
releases intermediate maps/host metadata before allocating the leaf pool, and
fills owned shared Metal storage directly without a full second leaf copy.
Larger unsupported allocations fall back rather than raising the allowance.
The small sampled profiling buffers fit within the remaining allowance.

All final timing processes stay below 15 s preparation, 90 s per process, 4 GiB
peak process memory and 1 GiB of Metal buffers. Maximum observed process peak in
the timing comparisons is approximately 3.12 GiB. The prototype's allocation
guards count its buffers in addition to the existing context. No production
memory policy or future CUDA staging strategy is implied by Metal shared storage.

## Trying the public experiment

```sh
make -j4 RELEASE=1 USE_METAL=1 \
  BUILD_DIR=build/native-gpu/release BIN_DIR=build/native-gpu/bin \
  build/native-gpu/bin/compute/metal/benchmark_metal_group_execute

python3 benchmarks/native_gpu_compare.py \
  --config benchmarks/metal_group_execute.example.json \
  --output build/native-gpu/group-execute/example-runs --blocks 2 --timeout 90
```

Use a new output directory. The example compares public fixtures with groups
disabled/enabled and the hash layout selected. Each process checks output parity.
The benchmark emits binary point results and timings; it does not render images.

The benchmark builds with only its execution metallib and the normal compute
metallib; it does not require an unrelated audit metallib. Private models, paths,
coordinates, raw profiles, outputs, configurations and timed/final source snapshots
remain under ignored `build/native-gpu/group-execute/`. Tracked documentation
contains aggregate results and independent public fixtures only.

Ray group evaluation, progressive rendering, full-resolution goals, transport
navigation and hardware-gated CUDA validation remain separate milestones.
