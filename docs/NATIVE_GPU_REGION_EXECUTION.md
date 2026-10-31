# Native GPU structural region audit

Completed 2026-10-31 on Apple M1. This implements the first decision stage of the
dated region-execution plan: a strict recognizer, native diagnostic pipelines,
public regressions and weighted eligibility measurements for both point APIs.
Production geometry evaluation remains unchanged.

## Decision

**Park the whole-program intersection-only point specialization.** On the fixed
private slice, it covers only 3.95% of all-owner instructions and 3.93% of ordered
selection instructions. Static eligibility considerably overstates its dynamic
coverage. The next proposed experiment is recognition of intersection groups
inside mixed expressions, beginning with a work-weighted audit of those groups.
No simple-region production switch or device evaluator was added.

These instruction percentages are not elapsed-time fractions or speedup bounds.
The decision reflects limited measured coverage for the declared point workload;
it does not establish the opportunity for other models or rays. Ray eligibility
has not been dynamically measured by this increment. Ray specialization,
navigation and CUDA remain later packages.

## Contract and implementation

`benchmarks/region_program_audit.h` recognizes only validated packed programs
containing signed primitive leaves, intersection combines and their structured
short-circuit skips. It checks program ranges, referenced primitives, stack shape,
subtree adjacency and exact branch destinations. Explicit complements, constants,
unions and differences remain general. Half-space flips already encoded in leaves
are supported. Single leaves are reported separately. The recognizer expects
unthreaded code; the audit rejects enabled jump threading. It supports both wide
and compact storage and does not rewrite geometry.

The host diagnostic enumerates original and executable spatial programs, sorts
and deduplicates exact ranges, and builds a bounded lookup table. Each instrumented
query attributes work using the actual selected `(root, count)` program, including
original fallback outside a spatial domain. It never infers spatial eligibility
from the source cell. Unmatched lookups fail the audit.

Counters distinguish general/constant, single-leaf and multiple-leaf intersection
programs, root versus child traversal and the five existing program-selection
reasons. They record calls, executed instructions, primitive evaluations, rejecting
instructions and rejecting calls. Sharding across 256 cell-index buckets reduces
diagnostic contention; each uint32 shard saturates and is checked, and host totals
use uint64. Metadata and counters together occupy approximately 8.40 MiB on the
private model. The combined host construction allocations are capped at 64 MiB.

The diagnostic compares every result to a fresh ordinary GPU query. All-owner
counters reconcile with existing per-point profiles for all five fields. Selection
counters independently reconcile calls, instructions, primitives and rejected
calls against the existing cell profiler. That profiler does not record rejected
instruction subtotals; this field is not independently reconciled for selection.
Matching diagnostic outputs establishes preservation, not CPU/GPU accuracy.

Hooks in `classify.metal` and `select.metal` compile only into diagnostic shader
libraries. Normal kernels acquire no new buffers or execution paths. The older
point work audit remains compatible with the extended hook arguments.

## Results

The same 16,384-point private slice was used for both APIs. Compact code was
enabled, jump threading disabled, and the comparison runner set four CPU threads.
Each diagnostic ran after ordinary reference classification in the same process.

| Measurement | All-owner classification | Ordered selection |
| --- | ---: | ---: |
| Evaluated program calls | 107,584,636 | 66,461,011 |
| Executed instructions | 3,671,369,563 | 2,050,865,241 |
| Primitive evaluations | 651,832,142 | 362,923,394 |
| Intersection-only calls | 17,465,200 | 9,513,068 |
| Intersection-only instructions | 144,908,805 (3.947%) | 80,525,078 (3.926%) |
| Intersection-only primitive evaluations | 31,680,990 (4.860%) | 17,836,353 (4.915%) |

The prepared snapshot has 537,380 unique executable ranges: 218,888 strict
multiple-leaf intersections, 101 single-leaf programs and 318,391 general/constant
programs. A deduplicated signed-reference pool for all eligible intersection ranges
would occupy 4,329,800 bytes before descriptors. Storage is feasible; dynamic
coverage is the limiting evidence. This count preserves distinct program ranges
and does not claim content-based deduplication across different addresses.

All-owner output is identical and counters reproduce the earlier work audit
exactly, with no saturation. Selected output is also identical and its four
independently available counter fields reconcile, with no saturation.

All-owner diagnostic submissions used 4,096 points. The first selected diagnostic
at that size failed with Metal's `Impacting Interactivity` command error. One
bounded retry used 256-point submissions with the same full point set and passed.
Its process completed in 56.7 seconds; the instrumented query took 38.9 seconds.
Those times include heavy profiling and are not production performance evidence.
The failed run remains preserved locally.

Successful processes stayed within the existing 15-second snapshot preparation,
90-second process, 4-GiB peak-process and 1-GiB Metal-buffer guards. The selected
retry reported 65.5 ms of additional recognition/table preparation, approximately
3.14 GiB peak process memory, and approximately 377 MiB of context Metal buffers
including selection profiling. Temporary diagnostic buffers are accounted for by
the separate allocation cap and peak process observation. No production gain,
occupancy diagnosis or full-frame result is claimed.

## Validation

- Eight builder tests pass, including a new recognizer test covering structured
  nesting, leaf orientation, compact decoding, unsupported operators, malformed
  branch layout, invalid ranges and primitive references.
- Public ordinary-versus-diagnostic outputs match in wide and compact modes,
  including Metal API validation, long intersection chains and ordered selection.
  The existing independent analytic ownership oracle also passes for all-owner
  cases. Selected diagnostic comparisons use the ordinary selected implementation.
- Existing compute tests pass in both code formats (74 each), ten resource tests
  pass, and the preceding point work audit still passes its oracle and totals.

## Reproduction and remaining work

```sh
make -j4 RELEASE=1 USE_METAL=1 \
  BUILD_DIR=build/native-gpu/release BIN_DIR=build/native-gpu/bin \
  build/native-gpu/bin/compute/metal/benchmark_metal_region_work

python3 benchmarks/native_gpu_compare.py \
  --config benchmarks/metal_region_work.example.json \
  --output build/native-gpu/region-execution/example-runs --blocks 1 --timeout 90

python3 benchmarks/summarize_region_work.py \
  build/native-gpu/region-execution/chain.csv \
  --output build/native-gpu/region-execution/chain-summary.json
```

Use a new runner output directory. Set `ALEA_BENCH_REGION_SELECT=1` for selected
queries. The public executable uses the existing independent sphere/hierarchy
fixture, optionally with `ALEA_BENCH_CSG_CHAIN=1`.

The next preparation experiment should measure short-circuit execution within
mixed programs and determine how much work lies in ordered intersection groups.
It must preserve signed-zero and error behavior, avoid Boolean expansion, and
account for descriptor/list duplication. Keep candidate reduction separate so
any gain can be attributed. Do not infer ray coverage from these point results.

Private probes, configurations, raw logs, model paths, views and source snapshots
remain under ignored `build/native-gpu/region-execution/`. Only aggregate results
and independently constructed public fixtures belong in tracked documentation.


## Follow-up checkpoint (2026-10-31)

The [intersection-group audit](NATIVE_GPU_INTERSECTION_GROUPS.md) is now complete.
Groups inside mixed expressions cover 72.427% of all-owner instructions and
70.869% of selected instructions when whole simple programs are included. Most
coverage is in mixed programs themselves. The whole-program-only experiment
remains parked; a native point group evaluator is the selected next experiment.
No production performance gain or ray coverage has been established.
