# Residual work in packaged Metal rays

Completed 2026-10-31 on Apple M1 after ordinary-ray group integration. This is a
bounded attribution pass through the public compute API, not another performance
comparison. The [package guide](NATIVE_GPU_RAY_GROUPS_PACKAGE.md) documents use.

## Protocol and checks

Each of three fixed private views renders 4x3 AA1 once with generic ordinary rays
and once with the packaged ray option. Processes are serial, compact code is on,
point grouping and jump threading are off, and profiling is enabled. The source
model is destroyed after snapshot creation; changing group settings afterward
checks captured context behavior. Profiling times are not acceptance timings.

All three pairs produce identical raw cell/material/depth/normal/color buffers.
Primitive and interval/cell-call totals also match. The summarizer reconciles
per-ray and per-cell totals and rejects saturated or inconsistent counters and
execution metadata. Generic profiles identify bytecode; native profiles identify
bytecode plus group operations, so instruction totals are not compared as equal
units. No private accuracy discrepancy is resolved by this preservation check.

## Aggregate observations

| View | Active packed cells | Cell calls | Primitive evaluations | Intervals per call | Container share of primitives | Top ten cells' primitive share |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Exterior | 25,649 | 98,633 | 2,672,665 | 3.430 | 85.057% | 5.913% |
| Cutaway | 18,909 | 70,565 | 1,981,525 | 2.847 | 85.994% | 5.853% |
| Close-up | 21,353 | 93,419 | 1,113,849 | 1.593 | 44.552% | 4.293% |

A cell call can perform several CSG interval evaluations. In these profiles,
`interval_steps` counts those evaluations, not physical material transitions.
The ratios show repeated expression evaluation inside calls, especially in the
container-heavy views. Active cells identify packed cells; repeated placements
can contribute to the same row, so these aggregates cannot establish that state
can be reused across occurrences or different coordinate frames.

The existing deferred-parent truth-segment cache hits approximately 0.004–0.005%
of parent checks in exterior/cutaway and none in close-up. The denominator
includes all parent checks, rather than only cache-eligible lookups. This describes that
specific cache and workload. It does not contradict the earlier primitive-root
reuse diagnostic or justify reviving its parked implementation. Concentration
is low: accelerating a handful of individual cells would cover little work.

## Selected next investigation

Measure whether conservative expression-validity ends cause avoidable restarts.
The current evaluator takes the earliest crossing of visited operands even when
a rejecting operand can certify the intersection remains false beyond it. A
hypothetical longer truth certificate might eliminate repeated interpretation,
root solving and probing, rather than merely execute the same work more cheaply.
The present counters do not quantify that opportunity or prove it safe.

The next package should compute shadow validity certificates on the same executed
leaf stream without changing stepping or results. Count expressions whose final
truth could remain certified beyond the current next crossing, separating
containers and terminal cells. Require independent analytic checks, identical
outputs and reconciled existing counters. Begin with strict nonzero truth;
uncertain, zero, error or numerical cases retain their original bounds.

Before any optimization, establish how a certificate preserves surface/cursor
identity, surrounding mixed operators, spatial-domain limits and first errors.
Skipping an otherwise evaluated nonfinite leaf is a contract change unless a
finite-domain proof excludes it. An intersection is not necessarily convex.
Measure this one mechanism before changing candidate indexes or introducing a
new cache. The observations here support investigation, not a speedup prediction.

## Reproduction

```sh
python3 benchmarks/summarize_ray_work.py /local/ignored/frame-prefix \
    --output build/native-gpu/ray-package/aggregate.json
```

Inputs are the standard `PREFIX-rays.csv` and `PREFIX-cells.csv` profiles from
`ALEA_METAL_PROFILE`. Output contains aggregates only, without source paths,
coordinates, cell IDs or geometry. An independent small fixture checks known
role totals and rejection of mismatched/saturated counters and mixed units.
Private models, camera files, profiles, images and source/API-run manifests remain
under ignored build storage. Only aggregate findings are included here.
