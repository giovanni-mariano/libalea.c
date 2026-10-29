<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO
SPDX-License-Identifier: MPL-2.0
-->

# Reusing exact contexts within nearest recovery

Follow-up to the [corrected baseline](GEOMETRY_CORRECTED_BASELINE.md), compared
with `abddc6f` (production code at `0e9451e`). All fixtures are independently
constructed. Private models and diagnostics remain outside Git or under ignored
`build/context-reuse/`.

## Implementation and contract

A nearest-owner recovery can encounter the same ordered surface list with the
same local ray and transform through several parent paths. Once a scan has
completed successfully, scanning it again over a contained window cannot add
an earlier hit to that same accumulating nearest-hit result. The implementation
uses this rule only within one recovery call.

The exact key includes list identity/count, all nine ray values (including
reciprocal directions), transform/inverse values, and flags. Hash collisions
require full key comparison. Explicit arrays avoid struct-padding comparisons;
no rounding or tolerance is used. Non-finite geometry inputs take the existing
scan; infinite reciprocals of zero directions remain valid ray data.

Entries record the successfully scanned window. Wider, disjoint, or uncertain
windows execute the scan. Only complete, non-interrupted scans publish entries.
Nearest-hit ties retain the first hit; aliases and world-normal transforms are
preserved. A repeated list skips only its indexed scan: child fill/lattice
recursion still runs. Interruption is checked before a reused scan returns.

The table has 2,048 entries and a maximum of sixteen probes per lookup. Its
allocation is about 624 KiB on the measured ABI, below a checked 1 MiB cap.
Allocation failure or probe saturation retains the existing scan. Each recovery
owns and frees its scratch; nothing persists across queries or geometry edits.
The result's scratch pointer and warm-up counter are internal, borrowed state.
Public APIs are unchanged and ordinary all-hit enumeration does not use reuse.

Allocation waits until the sixteenth indexed scan in a recovery. Initial eager
allocation caused a repeatable small-control regression: that workload made
58 short recovery calls per 512 frame. Delaying allocation avoids paying for
large scratch on those short calls. The larger key-building helper also stays
outside the ordinary indexed scan's stack frame.

## Validation

Three new focused tests cover exact key/window checks, ULP changes, aliases,
repeated parents, a nearer translated occurrence, forward/reverse rays, small
and saturated tables, disabled fallback, delayed allocation, interruption and
subsequent reuse, and independent elliptical-torus roots/normals across queries.
The repeated-parent fixture checks the expected crossing and first alias against
uncached enumeration while confirming fewer primitive solves. Its deliberately
overlapping parents test mathematical breakpoints, not unique physical ownership.
Existing transformed/lattice recovery and close-boundary tests remain active.

The final native unit/integration/Lua suite and serial unit/integration suites
pass. Sixteen focused geometry tests, twelve ellipse tests, and characteristic
tracks pass. These are targeted geometric and lifecycle checks, not exhaustive
proof for every degenerate imported model.

## Release measurements

Baseline and candidate use Apple Clang 21 on macOS ARM64 with
`-O3 -march=native -DNDEBUG`, one worker, and identical existing output budgets.
Processes run sequentially. Three alternating process pairs compare complete
512/1,024 public slices in split/fused modes, with three frames per process.
All 288 public frames match their baseline fields. Warm timing changes range
from approximately -4.7% to +3.0%, within the 5% screening threshold. Longer
repeating-lattice controls, with eighteen measured warm frames per process,
show approximately +0.9% at 512 and -0.4% at 1,024 after delayed allocation.
These small changes are not useful speedup claims.

Three alternating private process pairs at resolution 1,024 produce identical
completed digests and segment/interval, snapshot, and ownership-query counters
for every sampled row. The total time for the same 32 rows has these results:

| Build | Median total | Range across three processes |
| --- | ---: | ---: |
| Corrected baseline | 7.169 s | 7.121–7.576 s |
| Context reuse | 3.264 s | 3.257–3.448 s |

This is approximately 54.5% less tracing time, or 2.20 times the throughput for
this bounded workload. A separate 512 pair falls from 5.828 s to 2.688 s.
All 256 row results in these comparisons match. The measurements do not establish
a universal speedup or predict complete-frame timings from the sampled rows.

Eight private detail frames match the prior completed-view fields at both
resolutions and in both modes. No new preparation pass or persistent model
storage is added. The optional per-recovery allocation is separately bounded;
process RSS variation does not establish a whole-process memory improvement.

The final broad private 512 split frame completes with 44.699 s of frame work,
including 1.170 ms rasterization. It has 88,597 material pixels and no resolution
flags. Loading and preparation are additional, within the 60-second process
screen. The baseline did not complete that screen, so no complete baseline
field comparison or numerical full-frame speedup is available. Completion is
a visualization milestone, not validation of full-model physical ownership.
The one-second warm broad-frame target remains unmet.

The final fused broad frame also completes, at 45.267 s, with identical fields
to the split frame. These are single first-frame completion checks per mode;
no broad warm-frame series was measured. The initial eager-allocation candidate
had a fused timeout, so this should not be presented as a reliable hard
60-second completion guarantee under all system conditions.

Fresh independent 512/1,024 ellipse pixel oracles pass in both modes against
the final build. The final one-second cancellation workload completes seven
rows, cancels nine, reports no query errors, and returns in 1.006 s overall.
Successful queries follow cancelled ones; incomplete rows have no digest.
The last cancellation has only the remaining total budget, so its short return
is not a full 100 ms cancellation measurement. Full-budget requests return in
approximately 101–105 ms. Cancellation remains cooperative.

## Decision

Retain the bounded reuse change: independent regressions and full CPU suites
pass, completed controls match, sampled recovery improves repeatedly, and final
public controls show no repeatable regression above 5%. The storage and lifetime
are bounded without model invalidation machinery. No private geometry, paths,
coordinates, raw digests, logs, or extracted diagnostics enter tracked files.

The next investigation should profile the remaining cost of the final broad
view and its recovery/ownership work. The earlier corrected-baseline profile
predates this reuse and cannot establish the final bottleneck. Keep loading,
preparation, first/warm frame times, and display cost separate. This branch
still makes no new GPU performance or accuracy claim.
