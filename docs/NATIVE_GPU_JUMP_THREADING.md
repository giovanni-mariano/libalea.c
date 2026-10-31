# Native GPU jump-threading experiment

The opt-in preparation pass collapses identical short-circuit branch chains in
existing bytecode. Enable it with `ALEA_METAL_JUMP_THREADING=1`; the default is 0.
The pass runs after geometry/spatial preparation, point certification and BVH
construction, before compact encoding and device upload. It changes jump arguments
only. GPU arithmetic, evaluators and layouts are unchanged.

## Preparation contract

The pass enumerates original cell programs and executable spatial leaves, sorts
and deduplicates their ranges, then validates every range and branch target before
mutation. Nonidentical overlapping ranges are retained unchanged. Backward traversal
redirects union-to-union and intersection-to-intersection taken branches through
their already processed target. Other predicates and sign-changing instructions
are excluded. New arguments remain within the 28-bit compact limit.

A taken branch has already established the second test's predicate on the same
unchanged stack value. Skipping that redundant test preserves primitive visitation,
errors, ray crossing bounds and surface provenance. Program addresses, lengths,
non-taken paths and original instructions remain available. The operation is
idempotent. Range scratch is capped at 64 MiB; interruption aborts preparation and
the unpublished snapshot is discarded.

Consumer review: `compute_local.c` reconstructs expressions before this boundary;
point certification and spatial preparation also finish before it. The optional
host point-candidate proof scans every instruction and ignores skip arguments.
Device classification, selection, ordinary rays and resumed rays already interpret
jump offsets using the same program-counter convention. Compact encoding preserves
logical addresses. No consumer that needs the original structured jump layout runs
on the mutated snapshot.

## Validation

New builder tests compare original and rewritten primitive traces, result bits and
finite failures over all combinations of negative, zero, signed-zero, positive,
NaN and infinite input values for union/intersection chains. They also exercise
idempotence, exact compact decoding, duplicate and overlapping ranges, invalid
ranges/targets, interruption, mixed predicates and excluded reverse differences.

Existing compute tests pass in default, opt-in wide and opt-in compact modes (74
each); all ten Metal resource tests and seven builder tests pass. The optional
`ALEA_BENCH_CSG_CHAIN=1` fixture adds 64 enclosing sphere intersections without
changing the original sphere geometry, allowing independent ownership/intersection
oracles to exercise actual long chains. Point classification, ordered selection,
ordinary and resumed rays, and progressive camera cancellation pass Metal API
validation on that fixture.

The pass reports aggregate rewritten-site count and pass time only. Private raw
profiles, model paths, cameras, images and binary results remain ignored. CUDA
execution and existing CPU/private accuracy issues remain separate work.

## Measurements and decision

**Retained as opt-in; default remains disabled.** The private throughput gates
pass, but the public regression and preparation-overhead gates are unresolved.
No further compiler mechanism was added to compensate for them.

Apple M1 release build; sequential baseline/candidate processes in forward and
reverse order, two blocks, one initial and three warm calls per process. The
private ray views contain only 12 pixels (4x3), matching the bounded diagnostic
cameras; these are not detailed-frame performance results. Normal warm medians:

| Private workload | Baseline seconds | Threaded seconds | Elapsed-time reduction |
| --- | ---: | ---: | ---: |
| 16,384-point slice | 2.609–2.667 | 1.719–1.759 | 32.6–35.6% |
| Exterior, ordinary rays | 11.043–11.154 | 8.621–8.698 | 21.9–22.0% |
| Cutaway, ordinary rays | 8.397–8.415 | 6.340–6.365 | 24.4–24.5% |
| Close-up, ordinary rays | 4.491–4.614 | 3.667–3.750 | 18.3–18.7% |
| Exterior, resumable rays | 10.547–11.061 | 8.758–8.905 | 15.6–20.8% |

All retained point, exterior/cutaway/close-up and resumable output buffers match
their baseline byte-for-byte. The pass rewrites 22,696,276 sites over the complete
prepared scene, including spatial programs. This differs from the audit's static
count over visited original cell programs.

A separate instrumented fixed-slice run reduces executed instructions from
3,671,369,563 to **2,188,393,453** (40.4%), while primitive evaluations remain
**651,832,142**. Its query results also match baseline, with no saturated counters.
Those counters establish the mechanism; instrumented elapsed time was excluded
from the normal-timing table.

The pass costs approximately 0.30–0.54 s in private preparations, with about
4.74 MiB of range scratch and unchanged device-buffer lengths. All processes
remain within the 15 s preparation, 90 s process, 4 GiB peak process and 1 GiB
Metal-buffer guards. Peak observed process memory is approximately 3.28 GiB.
However, one paired total preparation comparison rises from 9.139 s to 10.192 s
(11.5%), above the plan's 10% gate; baseline preparation variability prevents a
clean overhead conclusion.

For cancellation, each baseline/threaded resumable process receives 20 simulated requests
100 ms after work starts. The measured latency is the overrun beyond that request
time, not total work duration. Empirical per-process p95 overrun is 12.8–15.4 ms
baseline and 13.4–14.2 ms threaded. Maximum observed continuation-command times
are 27.6–32.0 ms baseline and 26.4–29.4 ms threaded. The samples show no consistent
responsiveness regression, but do not establish production UI latency or behavior
under concurrent GPU load. Full resumed frames remain slow.

Public simple/long-chain cases pass correctness checks. Their short timings are
variable. A bounded repeat with nine warm observations per process still shows
long-chain ordinary rays regressing by 10.1% in forward order (10.146 to 11.171 ms)
and improving by 13.4% in reverse order (11.977 to 10.370 ms). This does not satisfy
the plan's requirement to resolve regression/noise before promotion. It is not
reported as a proven intrinsic regression either.

The implementation, tests, bounded comparison and promotion decision are complete.
Future promotion work should isolate public timing variability and establish
preparation overhead consistently, while retaining these exact semantic checks.
This result does not close existing accuracy blockers, full-resolution targets,
Metal/CPU parity or hardware-gated CUDA validation.

## Trying the experiment

Rebuild the Metal library/tools and set `ALEA_METAL_JUMP_THREADING=1` for the
process using Metal. Set it to `0` or omit it for the current default. Values other
than 0/1 are rejected; resource tests cover that error path and retry cleanup.
The switch affects newly prepared contexts.

Public performance reproduction:

```sh
make -j4 RELEASE=1 USE_METAL=1 \
  BUILD_DIR=build/native-gpu/release BIN_DIR=build/native-gpu/bin \
  build/native-gpu/bin/compute/metal/benchmark_metal_query_scaling

python3 benchmarks/native_gpu_compare.py \
  --config benchmarks/metal_jump_threading.example.json \
  --output build/native-gpu/jump-threading/example-runs --blocks 2 --timeout 90
```

Use a new output directory. The example alternates the simple and long-chain
public workloads with threading disabled/enabled. All private manifests, source
probes, camera data, binary outputs and raw measurements remain under ignored
`build/native-gpu/jump-threading/`.
