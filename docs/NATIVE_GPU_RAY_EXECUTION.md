# Native Metal ordinary-ray group execution

Package integration is now available as the supported opt-in
`ALEA_METAL_RAY_GROUPS=1` option, following explicit user direction. See
[the package guide](NATIVE_GPU_RAY_GROUPS_PACKAGE.md). The evaluation and
non-promotion decision below describe the preceding prototype and retain its
measured tradeoffs. Packaging does not change those results or the off default.

Completed R4b prototype evaluation, 2026-10-31, Apple M1. This is the R4b ordinary-ray
prototype following [the contract and coverage audit](NATIVE_GPU_RAY_GROUPS.md).
The compute package's supported point option remains independent. The native ray
prototype runs in a dedicated comparison executable pending its acceptance
decision; normal rendering retains the generic evaluator.

## Representation and semantics

Reuse the package's context-owned group maps and signed references, captured
before jump threading and compact encoding. Original/spatial program ranges use
exact lookup keys. Unknown ranges retain bytecode execution. No new model format,
expanded Boolean expression, coefficient copy or primitive solver is introduced.

At an eligible group start, a native ordered leaf loop reduces signed truth,
validity end and boundary identity. It preserves the first rejecting leaf and the
earlier visited leaves' validity limits. A strict minimum retains the left
boundary on equal distances. At-cursor identity follows the last visited matching
root. Both primitive-root order and the finite interior probe match the generic
ray evaluator. The loop checks errors before advancing and resumes generic
interpretation at the original group's end, retaining surrounding operators.

The initial implementation preserves primitive work. Profiled instruction counts
use remaining bytecode plus a group dispatch and visited leaves; they have a
different unit from the generic interpreter. Primitive counts remain comparable.
The experiment disables the separate root-reuse diagnostic.

`benchmarks/ray_group_execute.metal` owns this native evaluator. Its private
function constant 10 selects generic/grouped execution. Constant 11 removes audit
tracing from timed render pipelines; independent probes enable it to check leaf
order. Optional execution hooks in the ordinary renderer compile away in normal
builds. Buffer bindings 19/20 reuse owned maps/references; 17/18 are diagnostic
bindings. No native resumable evaluator or package ray setting is added here.

## Independent validation

The recursive analytic reference now checks 406 cursor/expression cases, including
axial and tangent plane/sphere rays, in wide/compact formats with grouping both
on/off. Checks cover truth, exact validity end, boundary, at-cursor identity and
visited leaf order/count. Generic one-opcode resumptions retain the same expected
expression state. Separate NaN/infinity plane cases verify evaluated errors and
correctly skipped errors in all four execution/format combinations.

Complete public renders compare raw pixels with fresh ordinary production output.
The driver also checks analytic depth, material, cell and source-oriented normal
expectations for translated fills and a separately constructed translated finite
lattice containing 27 spheres. Both are tested with and without a clipping plane.
Wide and compact native variants pass Metal API validation. These checks extend
the R4a expression reference to actual first-visible rendering decisions.

This is an initial ray feature matrix. Repeating/outer lattices, additional
primitive families, complex ancestor boundaries and near-coincident ill-conditioned
rays require broader coverage before package or progressive-render promotion.
Private CPU/Metal accuracy disagreements remain separate open cases. Fresh generic
pixel parity establishes preservation of those baselines rather than independent
private geometry accuracy.

## Comparison protocol

```sh
make -j4 RELEASE=1 USE_METAL=1 \
  BUILD_DIR=build/native-gpu/release BIN_DIR=build/native-gpu/bin \
  build/native-gpu/bin/compute/metal/benchmark_metal_ray_execute
mkdir -p build/native-gpu/ray-execute
python3 benchmarks/native_gpu_compare.py \
  --config benchmarks/metal_ray_execute.example.json \
  --output build/native-gpu/ray-execute/public-check --blocks 1 --timeout 90
```

`ALEA_BENCH_RAY_GROUPS=0|1` strictly selects the experimental renderer; default 0.
`ALEA_BENCH_RAY_VALIDATE=1` runs independent probes/public reference checks.
`ALEA_BENCH_RAY_LATTICE=1` chooses the finite lattice fixture;
`ALEA_BENCH_RAY_CLIP=1` clips the public fixture. These controls belong to the
benchmark executable. Validation settings use presence semantics. A separate
`ALEA_BENCH_RAY_PROFILE=1` process enables profiling, never acceptance timing.

Timed comparisons use one first call plus five warm complete renders per process,
fresh generic reference pixels, compact storage, jump threading off, AA1 and four
CPU threads. Each baseline/candidate pair shares point-group preparation because
that representation already belongs to the supported point snapshot. Thus this
measures incremental ray execution, not the preparation cost of enabling point
groups on an otherwise ungrouped snapshot. The executable currently requires the
prepared representation and reports failure if it is unavailable.

`experiment_setup_seconds` includes the reference render and optional validation;
it must not be interpreted as isolated pipeline or geometry preparation time.
Context preparation is reported separately. First/timed render values include
submission, synchronization and normal output handling. Profiling and leaf tracing
are disabled. Processes run serially in forward/reverse order, with a 90-second
per-process timeout and existing requested-allocation/resource checks.
`ALEA_BENCH_RAY_WARM=1..5` permits explicitly smaller diagnostic runs; five is
the default. Runs with fewer observations do not satisfy the timing gate.

Private tests use the camera-file interface from the R4a benchmark and fixed
4x3 AA1 views. These tiny complete renders are bounded performance probes, not
preview/full-resolution acceptance. Paths, cameras, model data, raw profiles and
outputs remain under ignored build storage. Tracked reports contain aggregate
measurements and independent synthetic fixtures only.

## Measurements and decision

Serial unprofiled baseline/candidate processes, then reversed order, five warm
complete calls each. Compact code, jump threading off, AA1 and four CPU threads.
The native shader is unchanged throughout. After the timeout, the host driver
only gains a bounded warm-count setting; its default five-call timing path is
unchanged. Source snapshots and executable identities are retained locally.

| Private 4x3 sample | Generic warm seconds | Native warm seconds | Time reduction in the two orders |
| --- | ---: | ---: | ---: |
| Cutaway | 8.223–8.294 | 6.639–6.668 | 19.6%, 19.3% |
| Close-up | 4.513–4.514 | 3.701–3.720 | 17.6%, 18.0% |

Every completed frame matches fresh generic production pixels byte-for-byte.
The exterior baseline exceeded the 90-second process guard during its planned
five-warm-call comparison. The runner killed the process group and stopped that
run. The failed evidence is retained. A bounded reverse-order/close-up follow-up
passes. Exterior has no accepted throughput comparison; its reduced-observation
profile is solely a work/output diagnostic. No timeout or resource limit is raised.

Public mixed-expression warm medians rise from approximately 0.730 to 0.769 ms
and from 0.662 to 0.732 ms, 5.3% and 10.6%. Both exceed the declared 5% regression
gate. The longer-chain case improves in both orders (about 26–43%), while the
lattice case changes within 5%. Tiny timings vary; these are complete render
observations rather than a diagnosis of GPU stalls or occupancy.

All ten final validation/regression processes pass: generic/grouped probes in
wide/compact formats, analytic fill/lattice and clipped images, eight material-
interval cases per mode/format, visited/skipped NaN/infinity errors, the 75-test
Metal compute suite and thirteen resource tests. The earlier ray coverage audit
also passes unchanged-output and counter reconciliation checks with the extended
reference. Resource tests and public native rendering use Metal API validation.

The maximum reported peak among the completed unprofiled private comparisons is
approximately 2.90 GiB; device allocation is approximately 476 MiB. Context
preparation is 9.61–10.12 seconds. These snapshots already contain point-group
metadata. Isolated pipeline preparation and total incremental ray setup have not
been established as a separate preparation gate. Complete-frame preview/detail
latency, ordinary cancellation and broader private accuracy remain open.

**Decision: retain the native ordinary-ray prototype; defer package/default and
resumable promotion.** The primary cutaway comparison misses the declared 20%
gate in both orders; close-up gains are smaller, the public regression gate fails,
and the exterior timing comparison remains incomplete. Correctness progress and
real local gains do not waive those conditions. The supported point-group package
option is unaffected by this ray decision.

The next useful experiment should target measured residual ray work, particularly
repeated interval/primitive evaluations, with a bounded attribution pass before
changing candidate reduction or interval reuse. Lookup/admission costs remain
possible contributors but are not isolated by these timings. Do not revive the
previous parked root-cache design without new evidence addressing its failure.
Native resumptions remain a later package requiring bounded leaf progress,
resident-state budgeting and camera/cancellation validation.


The separate twelve-ray exterior diagnostic completes with one warm observation
per mode and profiling enabled. Both generic and native paths perform exactly
**2,672,665 primitive evaluations**, without saturation, and return identical
pixels. These diagnostic times are excluded from the unprofiled table. This is a
sampled work-preservation check, not a complete-view performance result.
