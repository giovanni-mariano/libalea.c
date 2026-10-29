<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO
SPDX-License-Identifier: MPL-2.0
-->

# Geometry query experiment: integration checkpoint

Status reviewed on 2026-10-29. Latest retained implementation: ordered recovery
child-cell lists following `a8491e2` (ordered point queries), on
`experiment/geometry-recovery-search`; original baseline: `8c585eb`.
This page records the current decision. Linked experiment reports preserve
their original baselines, measurements, and historical conclusions.
The subsequent [corrected-baseline measurement](GEOMETRY_CORRECTED_BASELINE.md)
is complete. The subsequent [context reuse implementation](GEOMETRY_RECOVERY_CONTEXT_REUSE.md)
records that implementation, validation, and measurements. The latest
[ordered child-cell experiment](GEOMETRY_RECOVERY_CHILD_CELLS.md) is complete.

## Decision

Close the CPU query/recovery experiment with its hypothesis narrowed: existing
shared traversal can serve visualization and transport consumers with less
unnecessary work, and bounded recovery optimizations help the measured workloads.
Retain the changes for integration review. Interactive broad-view visualization,
general hardware acceleration, and complete transport geometry validation are
not established by this experiment.

The branch remains local and unmerged. Closing the experiment is a scope
decision, not a merge, release, or certification of arbitrary imported models.
Further performance work starts from the corrected checkpoint and receives its
own measurement and acceptance criteria.

## Changes retained

| Review group | Commits | Result |
| --- | --- | --- |
| Reference workloads | `6221d16`, `b859b99` | Public analytic fixtures, characteristic tracks, complete-slice comparisons, and bounded row profiling |
| Ownership-only queries | `9875c79` | Skip boundary enrichment when callers do not request surface output |
| Navigation correctness | `99106c4` | Resolve ownership inside a short flight when the initial forward sample overruns its next crossing |
| Recovery enumeration | `bbf414d`, `8c8789f`, `5a43f40` | Reuse ordered cell/universe cards and avoid repeated solves within one occurrence |
| Torus rejection | `0319988`, `ba5a3ba` | Profile the cost and reject clear finite-window misses before the quartic solve |
| Elliptical torus correctness | `0e9451e` | Align CPU roots, normals, and analytical slices with the evaluator's axial width |

Review correctness separately from timing. The navigation and ellipse fixes
are supported by independent expected results; equality to the old code would
preserve their defects. The ellipse change also replaces the earlier recovery
envelope with the true axial width. Preserve this dependency if reorganizing
commits. No history rewrite is required for review.

Public APIs are unchanged. Optional universe lists retain a 64 MiB per-model
storage cap and per-cell fallback; temporary seen-card bitsets are capped at
1 MiB. Surface order and aliases are preserved. Solved intersections are never
shared across differently transformed occurrences. Cache publication, index
invalidation, clone ownership, and destruction were inspected in the checkpoint
review and have focused lifecycle coverage. Concurrent geometry mutation is not
made safe by these cache changes.

## Evidence and its scope

| Evidence | Conclusion supported | Limit |
| --- | --- | --- |
| Native unit/integration/Lua and serial unit/integration suites pass at `0e9451e` | Tested CPU consumers remain integrated | Recorded runs, not a new cross-platform or GPU validation |
| 13 geometry tests, 12 ellipse tests, and analytic characteristic tracks | Targeted ownership, recovery, transform, normal, and distance contracts | Not exhaustive coverage of geometric degeneracies |
| Independent elliptical pixel equations at 512 and 1,024, split and fused | Correct material/cell/density/flag fields on the synthetic ellipse view | Does not validate arbitrary private-model ownership |
| 12,288 circular-ray comparisons | Root and normal digests unchanged for those controls | Compatibility evidence rather than an independent root oracle |
| Final 64 private sampled rows complete with unchanged digests | Sampled integration compatibility after the ellipse fix | Those views do not establish correctness of ellipse crossings |
| Four private imported-parameter axial checks: two failures before, zero after | The reproduced ellipse defect is corrected in those checks | Primitive audit only, not full occurrence or overlap validation |
| Final 100 ms cancellation requests return in about 100–105 ms | Observed responsiveness and subsequent query reuse | Cooperative cancellation is not a hard deadline |

The checkpoint review inspected the production diff, cache lifecycle paths,
focused tests, and existing validation logs. It identified no additional blocker
within this scope. That is a bounded review, not a substitute for independent
integration review. Consolidation changes documentation only; the CPU suites
need not be repeated solely for these edits.

## Performance status

On the 32 sampled private rows at resolution 1,024, prepared universe lists
reduced median total tracing time from 17.148 s to 10.281 s in their experiment.
The subsequent torus-rejection experiment measured 10.126 s to 6.893 s across
three alternating process pairs. These are separate measurements at different
checkpoints; do not multiply them into a general library speedup.

The final context-reuse build completes broad private 512 split and fused
first frames in 44.699 s and 45.267 s of frame work, with matching fields.
Loading and preparation are additional. These are single completion checks;
there is no completed baseline broad frame for field or speed comparison and
no broad warm-frame series at that checkpoint. The one-second objective remains
unachieved. Ordinary-resolution broad visualization is still slow.

For the repeated 32-row workload at resolution 1,024, three alternating process
pairs reduce median total tracing time from 7.169 s to 3.264 s (54.5% less).
All 256 comparison rows and eight detail frames match. Public control timing
changes remain below the 5% regression screen after delaying allocation for
short recoveries. Full native/serial suites and fresh ellipse pixel oracles pass.

The subsequent [ordered point-query fix](GEOMETRY_POINT_QUERIES.md) removes the
public single-point API's global hierarchy scan and uses ordered spatial
candidates. Its synthetic measurements are separate; the broad-view timings
above belong to the preceding recovery-context implementation. The current
checkpoint is measured separately below.

The [current broad-view profile](GEOMETRY_VISUALIZATION_PROFILE.md), measured at
`a8491e2`, now rechecks those outputs: one-worker first/warm frames take about
36.0/35.6 s, and four configured workers take 12.5/11.8 s. These are single-process
observations. The complete-frame diagnostic still attributes about 92% of
tracing to recovery, with 3.34 billion fill-loop cell visits. The subsequent
[ordered child-cell experiment](GEOMETRY_RECOVERY_CHILD_CELLS.md) implements
that hypothesis and measures three alternating complete-frame process pairs.
One-worker median warm time falls from 35.275 s to 24.961 s (29.2% less).
Four-worker median warm time falls from 12.779 s to 9.608 s, but the ranges
overlap and one candidate run is slower than its paired baseline. All 24 paired
broad frames match. These paired results supersede the exploratory timings
for evaluating this particular change; the one-second objective remains unmet.

## Latest retained implementation and next investigation

Nearest recovery now iterates optional ordered fill/lattice cell lists, skipping
terminal entries. Index allocations are capped at 8 MiB per model; unavailable
lists retain the original scan. Empty lists, publication, interruption, rebuilds,
clones, edits, budget exhaustion, and allocation failure have explicit coverage.
Final native unit/integration/Lua and serial unit/integration suites pass, as do
288 public frame comparisons, 256 private comparison rows, eight detail frames,
independent ellipse pixel oracles, and cancellation followed by reuse. The private
cache retains 44,344 bytes of indices plus 6,144 bytes of struct metadata.

Nearest recovery now avoids repeated successful indexed scans with an exact
list/ray/transform key and a containing previously scanned window. Optional
scratch is about 624 KiB on the tested ABI, allocated only after sixteen scans,
bounded below 1 MiB, and freed within the recovery call. Existing scans handle
allocation failure, saturation, or uncertain keys. Recursive children, aliases,
first-hit ties, and cancellation retain their contracts. No record crosses a
recovery call or geometry edit. See the implementation report for limitations.

The final diagnostic confirms that fill-loop cell visits fall from 3.34 billion
to 2.64 million, with all other recorded work counters unchanged. Indexed scans
still consume about 88% of recovery and 76% of total diagnostic tracing time,
with 815 million indexed references. Further work should target that remaining
surface-scan cost and recovery strategy. These are instrumented timing shares,
not release throughput measurements or justification for a GPU backend by
themselves. The one-second target remains unmet.

Transport validation remains a separate track: complete coverage, gaps,
overlaps, close/coincident boundaries, transformed occurrences, and navigator
restart behavior need independent expected outcomes. Existing fixed-source and
characteristic-track consumers are useful; selected-owner tracks alone do not
prove unique physical ownership. Deterministic mesh coupling and material
conservation also need their own contracts.

## Evidence index and confidentiality

- [Initial query experiment](GEOMETRY_QUERY_EXPERIMENT.md)
- [Occurrence-scoped recovery](GEOMETRY_RECOVERY_SEARCH.md)
- [Complete slices and bounded rows](GEOMETRY_SLICE_VALIDATION.md)
- [Prepared universe lists](GEOMETRY_UNIVERSE_RECOVERY.md)
- [Historical remaining-cost profile](GEOMETRY_REMAINING_RECOVERY.md)
- [Torus rejection](GEOMETRY_TORUS_RECOVERY.md)
- [Elliptical torus correction](GEOMETRY_ELLIPTICAL_TORUS.md)
- [Corrected baseline and repeated contexts](GEOMETRY_CORRECTED_BASELINE.md)
- [Bounded recovery context reuse](GEOMETRY_RECOVERY_CONTEXT_REUSE.md)
- [Ordered public point queries](GEOMETRY_POINT_QUERIES.md)
- [Broad-view profile before child lists](GEOMETRY_VISUALIZATION_PROFILE.md)
- [Ordered recovery child cells](GEOMETRY_RECOVERY_CHILD_CELLS.md)

Public fixtures are independently constructed. Private model inputs, paths,
coordinates, extracted geometry, raw hashes, diagnostic sources, and logs stay
outside Git or under ignored `build/`. Only aggregate private observations
belong in tracked reports. No model or diagnostic artifact is published by
this consolidation.
