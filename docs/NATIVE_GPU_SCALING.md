<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO
SPDX-License-Identifier: MPL-2.0
-->
# Native GPU query scaling and camera-motion evaluation

2026-10-30, Apple M1. This evaluates parallel point queries and rays with the
existing Metal evaluator. It does not change production scheduling limits,
geometry semantics or viewer quality policy. CUDA remains hardware-gated.

## Reproducing the public workload

```sh
make RELEASE=1 USE_METAL=1 bin/compute/metal/benchmark_metal_query_scaling
python3 benchmarks/native_gpu_compare.py \
  --config benchmarks/metal_query_scaling.example.json \
  --output build/native-gpu/scaling-new --blocks 2
```

`benchmark_metal_query_scaling OPERATION BATCH PASSES` supports `all-owner`,
`selected`, `rays`, `resume`, and `interactive`. The example selects representative
endpoints; the evaluation also measured intermediate batch sizes. All fixtures
are public synthetic geometry, sharing the existing analytic oracle with
`gpu_geometry.c` through `gpu_geometry_fixture.h`.

Point modes execute the same 16,384 points in the same order. Ray modes execute
the same 192x128 camera rays; only ordinary submission size or the progressive
resident-slot limit changes. Increasing image resolution would change the rays
and confound this comparison. Progressive storage remains allocated at its normal
capacity; the benchmark restricts how many slots are active. Point submissions
still obey the library's 4,096-query native cap.

The benchmark uses the native bridge, so timings omit public API validation but
include bridge staging, submissions, synchronization and readback. A benchmark-only
command-queue proxy records device start/end timestamps, submission counts and
maximum command time. It adds host overhead; compare within this harness rather
than equating its wall times with release-tool timings. These are throughput and
dispatch measurements, not hardware occupancy measurements. Scene setup, analytic
oracles and output comparisons are outside timing. First calls are reported
separately from five warm calls. Ray modes first execute an untimed ordinary
reference, so their pass zero is not a cold-shader measurement. Runs are serialized,
then reversed by block.

## Public fixed-work results

Ranges are the two per-process warm medians, not percentile estimates.

| Workload | Batch / active slots | Wall time, ms | Sum of GPU command times, ms |
| --- | ---: | ---: | ---: |
| All-owner, 16,384 points | 64 | 142.23–150.28 | 54.32–59.76 |
| All-owner, same points | 4,096 | 2.31–2.37 | 0.92–0.95 |
| Selected, same points | 64 | 122.65–134.21 | 33.00–34.31 |
| Selected, same points | 4,096 | 1.46–2.12 | 0.36–0.60 |
| Ordinary, 24,576 rays | 64 | 199.27–209.77 | 54.53–62.40 |
| Ordinary, same rays | 4,096 | 4.25–4.45 | 2.39–2.45 |
| Ordinary, same rays | 24,576 | 2.07–2.49 | 1.48–1.64 |
| Progressive, same rays | 64 | 236.99–242.27 | 93.28–95.59 |
| Progressive, same rays | 4,096 | 10.07–10.57 | 7.67–8.23 |

The same queries complete much faster when enough work is submitted together.
This includes amortizing command overhead and changing available parallel work;
it does not identify their individual contributions or prove maximum occupancy.
Every timed public call passed the analytic checks. Every chunked/progressive ray
image also matched its ordinary Metal reference bitwise in all output fields.

## Private workload and limits

Only aggregate findings are tracked. Models, views, coordinates, images, binary
outputs, hashes, generated probes and raw logs remain ignored locally. The private
observations below have fewer repetitions than the public comparison and do not
establish a production speedup or percentile latency.

For the same 16,384 all-owner points, 256-point submissions took 27.229 s on the
first call and 27.204 s warm. At 4,096 points, the first call took 2.646 s and two
warm calls took 2.636 and 2.613 s. Results matched exactly. GPU time accounted for
nearly all of the private frame time. The smaller-batch process stopped before a
third frame because its predicted runtime exceeded the process budget.

For the fixed 16x12 exterior ray set, 64-ray submissions did not finish within the
90-second process guard. The runner killed that process and stopped its comparison;
it is a failed observation, not a timing sample. A separate 192-ray submission
completed in 35.103 s, including 35.085 s on the GPU. No higher private image
resolution was attempted. A 32-slot progressive run was skipped after the timeout.
The 192-slot progressive run completed in 27.207 s, with 26.073 s GPU time over
1,650 commands; the longest command was 27.261 ms. It matched the ordinary image
bitwise after cancellation and camera reset.

A separate, unretained host-bridge diagnostic raised the native point cap to
16,384. It improved the simple public workload, but failed on the private workload.
A retry captured one successful first call at 1.754 s, followed by:

```
Metal compute failed: Impacting Interactivity
(0000000e:kIOGPUCommandBufferCallbackErrorImpactingInteractivity)
```

The comparison stopped and did not complete its reversed block. This is not an
accepted speedup, nor evidence of a universal command-duration threshold. The
production cap remains 4,096. Increasing concurrency must be evaluated together
with bounded work per submission and device interactivity.

## Query-cost variation and work attribution

A separate instrumented pass over the same 16,384 private points recorded executed
CSG instructions per point. Instrumented timings are excluded from the performance
comparison; results still matched the uninstrumented outputs exactly.

| Per-point executed instructions | Count |
| --- | ---: |
| Minimum | 135,474 |
| Median | 208,704.5 |
| Empirical 95th percentile over this fixed grid | 330,981 |
| Maximum | 820,960 |
| Total | 3,671,369,563 |

No counters saturated. Original programs account for 71.89% of instructions and
spatial programs for 28.11%. Tests rejecting the point account for **99.82%** of
instructions. This is not a claim that all rejection work can be removed. It
motivates measuring conservative candidate rejection and specialization coverage
by executed work rather than raw cell counts. The current aggregate profile does
not by itself identify which proof or preparation algorithm would succeed.

The earlier, separately sampled twelve-ray exterior profile ranges from 229,217
to 1,974,830 instructions per ray, median 1,403,044.5. It is useful evidence of
uneven work, but is not the same ray set as the 192-ray scaling test and cannot be
used as its throughput denominator.

## Camera motion: preview first, detailed rays after idle

The existing viewer already implements the control flow:

1. A camera change invalidates the current revision and requests cancellation.
2. Render a reduced-resolution preview with AA1 and shadows disabled.
3. After 150 ms without input and completion of the current preview, start detail
   at the requested resolution and supported quality settings.
4. Seed detail only from a preview of the same revision. Reject stale partial and
   completed frames. Saving renders the requested full quality.

Metal previews adapt between 16 and 128 pixels on the longest image edge; large
models start at 16. This policy is present, but its latency is not guaranteed.

The new `interactive` benchmark measures moving-camera 96x64 previews followed by
384x256 detail, AA1, on the public fixture. Warm preview medians were 3.42–3.51 ms;
warm detail medians were 31.63–32.55 ms. Initial previews took about 9.9 ms. It also
cancels a detailed render at a callback threshold and checks that a subsequent
camera renders identically to an ordinary reference. The benchmark excludes
postprocessing, SDL presentation and the idle debounce; it is not a complete UI
latency measurement.

The private 16x12 progressive preview took **27.207 s**, with first completed pixels
after **3.343 s**. A stale-camera cancellation returned 16.838 ms after its simulated
100 ms threshold, and the largest callback gap during the completed preview was
28.638 ms. These are individual observations, not p95 cancellation guarantees.
Fast cancellation and fast fresh-image production are different requirements.

The five existing SDL shell regressions pass: responsive input, close during work,
preview/idle refinement/full-quality saves, progressive publication with stale-frame
rejection, and progressive close cancellation. No viewer policy change is claimed
by this evaluation.

### Next implementation boundary

Keep parallel queries/rays as the basic native GPU work units. Evaluate throughput
and responsiveness separately: enough resident queries to use the device, bounded
work per command to let new input take effect. The point-cap failure supports
investigating bounded point continuations before simply increasing native limits.
It does not establish that continuations will improve total execution time.

For 3D, resolution reduction alone has not produced a fast private motion preview.
Evaluate a separate cheap display preview—such as cached-depth reprojection or a
prepared geometric proxy—while retaining analytic rays for settled detail. Such
previews must be explicitly approximate, handle disoccluded/missing regions, never
supply physics/navigation answers, and never be exported as a completed detailed
render. Reprojection alone cannot solve first-image latency. The representation
and its startup/memory cost need a bounded experiment before adoption.

The useful-preview target remains 200 ms, the detailed-frame target remains 2 s
at the declared final size, and the cancellation target remains 100 ms p95. The
private preview/detail targets are not met. Existing private accuracy discrepancies,
public CPU oracle failures and CUDA hardware validation remain separate open work.

### Cached-depth experiment outcome

The [native Metal reprojection experiment](NATIVE_GPU_REPROJECTION.md) now measures
that candidate separately. Warm 256x192 projection medians are below 0.45 ms, but
small-motion holes and ghost pixels remain, and a 16x12 cached frame cannot supply
a useful dense preview. It is retained as a benchmark, not integrated into the
viewer. A useful first-image representation remains the next prerequisite.
