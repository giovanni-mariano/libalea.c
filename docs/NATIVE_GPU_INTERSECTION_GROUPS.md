# Native GPU intersection-group audit

Completed 2026-10-31 on Apple M1. The audit extends the strict whole-program
recognition in [the region report](NATIVE_GPU_REGION_EXECUTION.md) to maximal
intersection-only subtrees inside mixed expressions. It measures a candidate
execution form; production evaluators and defaults are unchanged.

## Result and decision

**Proceed to a bounded native Metal point group-evaluator experiment.** The
eligible subtrees cover much more executed work than whole simple programs:

| Measured fixed-slice work | All-owner classification | Ordered selection |
| --- | ---: | ---: |
| Total instructions | 3,671,369,563 | 2,050,865,241 |
| Instructions within eligible groups | 2,659,059,389 (72.427%) | 1,453,432,653 (70.869%) |
| Group instructions inside mixed programs, as share of total | 68.480% | 66.943% |
| Total primitive evaluations | 651,832,142 | 362,923,394 |
| Primitive evaluations within groups | 622,609,479 (95.517%) | 342,144,202 (94.274%) |
| Group entries | 333,088,234 | 183,067,625 |
| Executed intersection skips within groups | 1,746,928,665 | 952,211,874 |
| Executed intersection combines within groups | 289,521,245 | 159,076,577 |

These are operation counts, not elapsed-time fractions or predicted savings.
The initial specialization preserves primitive evaluation order and therefore
does not eliminate the surface tests represented in this table. It aims to
replace their surrounding interpreter control flow with a direct ordered loop.

Only about 1.87 primitives execute per entered group on average. Long static
groups often reject early. Descriptor access and group dispatch must therefore
be cheap; a large static group is not evidence of a large profitable execution.
The measured skip/combine share also overlaps the earlier jump-threading
opportunity. Compare against both fixed threading configurations before assigning
an improvement to the new mechanism.

The next prototype specification is saved in the ignored dated implementation
plans. No production group evaluator, new default, 3D speedup or CUDA result is
claimed by this audit.

## Recognition and execution attribution

`region_group_audit.h` parses validated unthreaded packed programs with a bounded
stack. It checks branch kind, destination, subtree adjacency, operator arity,
primitive references and complete program coverage. A group contains only signed
primitive leaves and their existing intersection skips/combines, with at least
two leaves. Constants and explicit complements are not group operands. A pure
child beneath a complement, union, difference or mixed intersection can still
form a group. No reassociation, normalization, deduplication of leaf occurrences
or geometry change is performed.

Groups are maximal in the existing expression tree. For example, a mixed
left-associated expression ending with two separate leaf intersections does not
automatically gain a group by regrouping those leaves. That is a different
optimization and remains outside this experiment.

The host builds a separate sorted group list for each exact executable program
range. Exact original/spatial range aliases share a map. Partially overlapping
programs retain their own maps: a single-leaf view inside a larger intersection
must not inherit its parent's eligibility. This refinement left the measured
private totals unchanged, but is required for a valid general diagnostic.

The diagnostic shader counts instructions at the program counters the original
evaluator actually executes. It never evaluates additional leaves to establish
eligibility. A forward cursor follows the selected program's group list across
short-circuit jumps. Per-call totals join the existing region/reason/depth
counters, distinguishing mixed programs from whole simple programs.

Instrumentation is conditional on `ALEA_GROUP_WORK_AUDIT`. Normal shader builds
have no group buffer, lookup or counter work. Earlier region and cell work audits
remain compatible. The public group audit has its own native Metal library and
executable; it does not alter source geometry or prepared bytecode.

## Validation and resources

The same private 16,384-point slice was used for both APIs, with compact code,
jump threading disabled and four CPU threads. All-owner diagnostic submissions
contain 4,096 points; selected submissions contain 256 points, using the previously
established diagnostic limit. Every result matches a fresh ordinary GPU query.
All-owner totals reconcile with existing profiles for all five original fields.
Selected totals reconcile for calls, instructions, primitives and rejected calls;
the existing selector profiler lacks a rejected-instruction subtotal.

For every group-counter row, group instructions equal the sum of its leaf,
intersection-skip and intersection-combine counts. Group instruction/primitive
totals cannot exceed the corresponding complete-program totals. Unknown programs,
overlapping groups within a program and saturation fail the audit.

The snapshot contains 1,468,977 groups over 537,380 unique executable ranges,
covering 21,401,274 static instructions. Program-local maps occupy 20,349,900 bytes;
region metadata and sharded counters occupy 9,007,684 bytes. Construction accounts
for host program records, range-vector capacity, device maps and counters together
before upload, within the 64-MiB additional-allocation guard. The observed combined
accounted peak for this construction is about 52.2 MiB. The range vector has a
separate maximum of 24 MiB and may be capped lower by the remaining allowance.

The final private processes complete in 25.4 s and 68.8 s, with snapshots prepared
in approximately 9.7 s and 9.6 s. Group construction is included in total process
time; isolated group-construction time was not recorded in these private runs.
Peak process memory stays below 4 GiB, and context plus temporary diagnostic Metal
buffers remain below 1 GiB. Instrumented query times are not throughput evidence.

Public validation includes wide/compact code, Metal API validation, long groups,
mixed unions and ordered selection. An independent analytic sphere/hierarchy
oracle checks all-owner geometry. A separate group-work oracle derives expected
leaf/skip/combine/entry counts from the public point positions. A native map probe
checks nested and partial aliases, single-leaf views and missing ranges.

Nine builder tests pass, including mixed-tree recognition, compact decoding,
unsupported/invalid structure, callback failure and interruption. Both existing
compute configurations pass (74 tests each), ten resource tests pass, and both
preceding audit executables still pass their output and counter checks. This does
not close existing independent CPU/private accuracy discrepancies.

## Reproduction

```sh
make -j4 RELEASE=1 USE_METAL=1 \
  BUILD_DIR=build/native-gpu/release BIN_DIR=build/native-gpu/bin \
  build/native-gpu/bin/compute/metal/benchmark_metal_group_work

python3 benchmarks/native_gpu_compare.py \
  --config benchmarks/metal_group_work.example.json \
  --output build/native-gpu/group-execution/example-runs --blocks 1 --timeout 90

python3 benchmarks/summarize_region_work.py \
  build/native-gpu/group-execution/mixed-wide.csv \
  --output build/native-gpu/group-execution/mixed-summary.json
```

Use a new runner output directory. `ALEA_BENCH_CSG_GROUPS=1` selects an independently
constructed public union-of-intersections fixture with unchanged analytic geometry.
`ALEA_BENCH_REGION_SELECT=1` selects ordered queries. The audit requires threading
disabled because recognition checks the original structured branch destinations.

Private paths, models, cameras, source probes and raw measurements remain in
ignored `build/native-gpu/group-execution/`. Tracked material contains aggregate
findings and public fixtures only. Dynamic ray coverage remains unmeasured.


## Native execution follow-up (2026-10-31)

The [native group-evaluator experiment](NATIVE_GPU_GROUP_EXECUTION.md) is now
implemented and evaluated for both point APIs, separately with threading off/on.
It delivers private throughput gains with unchanged outputs and sampled primitive
work. Public timing regressions and preparation variability keep it benchmark-only;
normal library and viewer behavior remain unchanged.
