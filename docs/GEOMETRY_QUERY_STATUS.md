<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO
SPDX-License-Identifier: MPL-2.0
-->

# Geometry query experiment: integration checkpoint

Status reviewed on 2026-10-29. Production checkpoint: `0e9451e` on
`experiment/geometry-recovery-search`; original baseline: `8c585eb`.
This page records the current decision. Linked experiment reports preserve
their original baselines, measurements, and historical conclusions.
The subsequent [corrected-baseline measurement](GEOMETRY_CORRECTED_BASELINE.md)
is complete and supplies the current performance profile.

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

The corrected-baseline broad private 512-view screen still exceeds a 60-second
whole-process limit, including loading and preparation. Repeated completed
material-detail views take about 2.9–3.0 s warm at 512 and 5.8 s at 1,024.
Neither those views nor sampled rows establish broad-view interactivity.

The fresh coarse profile attributes about 87–90% of sampled row time to recovery
and about 83% of recovery time to indexed scanning. A separate shadow diagnostic
finds that 81.6% of indexed calls repeat an exact list/ray/transform context within
one recovery, over a contained window. It skips no work and establishes an
opportunity, not a speedup. No new GPU code, GPU throughput claim, or device
accuracy validation is part of this branch.

## Next bounded investigation

The corrected-baseline measurement is complete. All public controls and
independent ellipse pixel checks pass, as do 36 private detail frames and
64 private rows against their prior outputs. The sampled 100 ms cancellation
requests return in 103–105 ms with successful subsequent queries. The broad
one-second warm-frame target remains unmet.

Test one candidate: avoid repeating an indexed scan within a single nearest
recovery after an identical list/ray/transform context has completed successfully
for a containing window. Preserve first-hit ties, aliases, world normals,
recursive child traversal, and cancellation. Use exact keys, at most 1 MiB of
optional per-recovery scratch, and the existing scan on allocation failure,
saturation, or uncertain inputs. No cache crosses recovery calls or model edits.

First establish the rule with independent synthetic cases and the uncached
fallback. Require repeated release measurements, unchanged completed controls,
no repeatable control regression above 5%, explicit memory accounting, and a
new broad-frame screen before retaining an implementation. The shadow count of
79% repeated references cannot be converted into an expected speedup. The
slice benchmark excludes display/upload and does not establish 3D render or
application first-image performance.

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

Public fixtures are independently constructed. Private model inputs, paths,
coordinates, extracted geometry, raw hashes, diagnostic sources, and logs stay
outside Git or under ignored `build/`. Only aggregate private observations
belong in tracked reports. No model or diagnostic artifact is published by
this consolidation.
