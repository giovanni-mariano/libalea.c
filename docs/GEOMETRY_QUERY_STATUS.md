<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO
SPDX-License-Identifier: MPL-2.0
-->

# Geometry query experiment: integration checkpoint

Status reviewed on 2026-10-29. Production checkpoint: `0e9451e` on
`experiment/geometry-recovery-search`; original baseline: `8c585eb`.
This page records the current decision. Linked experiment reports preserve
their original baselines, measurements, and historical conclusions.

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

The latest broad private 512-view screen still exceeded a 60-second
whole-process limit, including loading and preparation. That screen was run
before the ellipse correction. It remains the latest broad-view observation;
there is no fresh full-frame measurement at `0e9451e`. Completed material-detail
views and sampled rows do not establish broad-view interactivity.

The earlier profile attributing about 89% of recovery time to indexed-list
scanning predates torus rejection. It motivates a new profile but cannot identify
the dominant cost of the corrected checkpoint. No new GPU code, GPU throughput
claim, or device accuracy validation is part of this branch.

## Next bounded investigation

Measure the corrected checkpoint before choosing another index or backend.
Use complete 512 and 1,024 slices for public controls and the existing private
detail/broad views. Record loading, preparation, first frame, subsequent frames,
trace/raster costs where available, process peak RSS, and cancellation separately.
Run private processes sequentially; preserve the existing trace/storage budgets
and external hard timeouts. A timeout is an incomplete result, not a frame time.

The initial engineering target is a complete 512-by-512 broad slice within one
second after preparation, with 100 ms requested cancellation returning within
150 ms on the measured expensive workload. These are investigation targets,
not achieved results or universal latency guarantees. Record the distance to
the target even if no candidate is implemented. The current slice benchmark
excludes display/upload and cannot establish application first-image latency
or 3D render performance; those require separate consumer measurements.

If the complete frame remains too expensive, use the existing bounded 32-row
sets to profile actual recovery, occurrence traversal, primitive work, and
ownership work on the corrected baseline. Separate diagnostic instrumentation
from release timing. Select at most one candidate from that evidence. Require
independent checks for any changed geometry decision, unchanged completed
control outputs, a repeatable benefit exceeding run variability, no repeatable
control regression above 5%, and explicit preparation/memory accounting.
Failure to meet the target is grounds to narrow or reject a candidate.

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

Public fixtures are independently constructed. Private model inputs, paths,
coordinates, extracted geometry, raw hashes, diagnostic sources, and logs stay
outside Git or under ignored `build/`. Only aggregate private observations
belong in tracked reports. No model or diagnostic artifact is published by
this consolidation.
