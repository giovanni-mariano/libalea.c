# Native GPU work-weighted specialization audit

Completed 2026-10-31. This adds a diagnostic pipeline and selects a preparation
optimization for the next experiment. **It does not enable an optimization or
claim a rendering speedup.** The production classifier's behavior is unchanged.

## What was measured

The previous point profiler counted calls by specialization/fallback reason, but
recorded executed work only as original versus specialized. The new benchmark
attributes calls, executed instructions, primitive evaluations, rejecting
instructions and rejected calls to each cell and reason. Its separate Metal
library enables guarded hooks in `classify.metal`; those hooks are absent from the
normal build. No framework abstraction or new geometry evaluator was added.

The same fixed 16,384-point slice is queried by ordinary and instrumented Metal
pipelines. Every output is compared byte-for-byte, and all cell counters reconcile
with the existing per-query counters. Saturation fails the benchmark. Shared
uint32 counters accumulate over bounded dispatches; host totals use uint64.
Counters occupy 42.977 MiB on the private model, below the 64 MiB limit. An initial
attempt stopped before queries because duplicating the table on the host exceeded
that limit; it was replaced by direct shared-buffer reads, not a larger budget.

Raw CSV contains cell identity and program metadata and must remain local for
private models. `benchmarks/summarize_point_work.py` aggregates it without emitting
identifiers or coordinates. Optional existing ray-profile joins check cell index,
cell ID, material and program length before comparing work. This is a metadata
consistency check, not a substitute for an identical-snapshot fingerprint.

## Results on the fixed private slice

Totals exactly reproduce the earlier aggregate profile: **3,671,369,563 executed
instructions**, 107,584,636 cell-program calls and 651,832,142 primitive evaluations.
Rejecting evaluations account for 99.82% of instructions. The new breakdown is:

| Program selection | Executed instructions | Share of all work |
| --- | ---: | ---: |
| Certified container, outside its specialization domain | 1,846,654,752 | 50.30% |
| Specialized container program | 1,032,138,626 | 28.11% |
| Terminal cell without a spatial program | 521,579,773 | 14.21% |
| Certified container, selected region unchanged | 134,871,974 | 3.67% |
| Container with unproven finite evaluation over the domain | 84,612,472 | 2.30% |
| Container without a spatial program | 51,511,966 | 1.40% |

The certification classification follows the current packer's budget/order and
records whether a point certificate is available. No visited uncertified work was
attributed to exhaustion of the proof budget. An unproven-finite result does not
mean the actual sampled evaluation was nonfinite: the interval proof is
conservative. Raising the proof budget is not supported as the priority here.

All 30,782,024 private outside-domain calls rejected. **That is an observation,
not permission to reject by domain.** The public sphere/void hierarchy has 7,168
accepted outside-domain calls on its fixed grid. Specialization domains choose
where cached programs may run; they are not necessarily cell bounds. Replacing
fallback with unconditional rejection would break existing valid ownership.

Work is broad rather than concentrated: the largest cell contributes 0.176%, the
top ten 1.20%, the top 100 7.44%, and the top 1,000 45.65% of slice instructions.
This disfavors a few model-specific special cases.

## Connection to 3D

The existing 12-ray profiles were reused without another private render. Matched
point/ray metadata covers 93.35%, 94.81% and 94.40% of exterior, cutaway and close-up
ray instructions respectively. Unmatched ray cells remain unclassified by this
point audit. The top 1,000 point-cost cells account for only 20.36%, 22.77% and
20.16% of those complete ray profiles. A point-hot-cell optimization therefore
cannot be assumed to cover the ray bottleneck.

The bytecode inspection found a broader candidate: a conditional skip can land
on another identical skip. Among original programs of visited cells, 7,959,955 of
9,700,976 static skip instructions have that directly threadable pattern (82.05%).
170,082 of 170,106 visited cells have at least one such site. These counts are per
cell program, not a deduplicated inventory of underlying stored bytecode.

Almost all point work belongs to cells whose **original** programs have such
chains; matched cells with chains cover 93.34–94.81% of the three ray profiles.
Neither statistic measures how often those particular jumps execute, nor whether
specialized programs still contain them. The ray profiles and point workload are
different samples. No percentage here is a predicted speedup.

## Selected next optimization

Test **short-circuit jump threading during bytecode preparation**: when a taken
union/intersection skip lands on the same predicate, redirect it to the end of
that chain. Preserve non-taken paths, program lengths/addresses, primitive order,
error handling and all numeric operations. Initially exclude reverse-difference
skips, sign-changing operations and mixed predicates.

This is smaller than adding new spatial bounds or translating another CPU
subsystem. The existing point, ordinary ray and resumable ray interpreters use
the same branch conventions. Removing redundant predicate tests can therefore
benefit each without changing their geometry calculations. Actual gains must
still be measured independently for each workload.

The detailed implementation and acceptance plan is saved in the repository's
ignored `plans/native-gpu-jump-threading.md`. It requires fixtures with actual
long chains (the current public sphere fixture has none), trace-equivalent tests,
wide/compact validation, ordinary/resumable ray parity, cancellation measurements,
preparation accounting and controlled performance comparisons before promotion.

## Validation and limits

- Public work audit passed in compact and wide bytecode; the wide run enabled
  Metal API validation. Cell/query totals reconcile; ownership matches the
  independent fixture oracle and ordinary Metal outputs.
- The final private run returned identical outputs, reconciled every counter and
  reported no saturation. It completed in 24.31 s. Preparation and measured peak
  resources remained inside the existing 15 s / 4 GiB / 1 GiB model-buffer guards.
- Existing regressions pass: 74 compute tests in each bytecode format, five builder
  tests and ten Metal resource tests.
- The diagnostic's extra atomics and per-call recording affect scheduling and
  timing. Its elapsed time is not production performance evidence.
- No default viewer, query or renderer behavior changed. CUDA hardware validation,
  existing CPU/private accuracy discrepancies and full-frame targets remain open.

## Reproduction

On a Metal-capable Mac, from the repository root:

```sh
make -j4 RELEASE=1 USE_METAL=1 \
  BUILD_DIR=build/native-gpu/release BIN_DIR=build/native-gpu/bin \
  build/native-gpu/bin/compute/metal/benchmark_metal_point_work

python3 benchmarks/native_gpu_compare.py \
  --config benchmarks/metal_point_work.example.json \
  --output build/native-gpu/work-audit/example-runs --blocks 1 --timeout 60

python3 benchmarks/summarize_point_work.py \
  build/native-gpu/work-audit/public-compact.csv \
  --output build/native-gpu/work-audit/public-summary.json
```

Use a new runner output directory. CSV outputs in the example are overwritten.
Private sources, model paths, cell-level profiles and raw logs remain ignored;
only the aggregate findings above belong in tracked documentation.
