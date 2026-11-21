<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO

SPDX-License-Identifier: MPL-2.0
-->

# CPU particle navigation baseline

The particle replay drives the public FP64 navigator. The box workload has paired
periodic X planes, a specular −Y plane, a white +Y plane, and a vacuum +Z exit.
The operation script fixes collision distances and the direction returned after
the white crossing. Seeded transport tests exercise the physical sampler separately.
Run from the repository root:

```sh
make RELEASE=1 bin/model_lifecycle
out=build/transport-navigation
mkdir -p "$out"
for mode in fast strict interval; do
    bin/model_lifecycle navigator-script benchmarks/workloads/transport/box.i \
        benchmarks/workloads/transport/box.ops "$mode" > "$out/box-$mode.json" || exit 1
done
python3 tools/compare_navigator_replays.py "$out/box-strict.json" "$out/box-interval.json"
make RELEASE=1 check-model-lifecycle
```

Each script line contains an operation, particle slot (`parent` or `child`), and
numeric arguments. A line beginning with `#` is a comment. `restart` takes XYZ
and direction, `advance` takes collision and maximum flight distances (`inf`
allows a collision-free flight), and `direction` takes a direction vector.
`clone` replaces the other slot with an independent copy. `reflect` and
`periodic` apply the pending action. A navigator error produces a trace through
the failed operation and a nonzero process status. The analytic regression checks
the event order, endpoint positions, periodic mappings and which event wins on
either side of an exact tie.

For collision-free model flights, use the same binary with `navigator` and a ray
workload. Each flight records its input, initial owner, successful events and
termination. Events contain kind, distance, surface ID, boundary type, normal,
incident and outgoing ownership, and world position. Location records contain
kind, cell ID, material ID, density and occurrence key. Errors retain the last
completed flight distance and diagnostic. The comparison tool checks the trace
in order and reports the first differing field and operation or ray. Matching
failed traces reproduce a baseline failure. The navigator command itself exits
nonzero when advance or restart fails.

To compare both navigators on a checked-in model, run the box script above
again with the regional accelerator and compare its events.

```sh
bin/model_lifecycle navigator-script benchmarks/workloads/transport/box.i \
    benchmarks/workloads/transport/box.ops strict regional > "$out/box-regional.json"
python3 tools/compare_navigator_replays.py "$out/box-strict.json" "$out/box-regional.json"
```

## Opt-in transport recovery

Replay the same workload with `transport` in place of `strict`. Both
producers retain the selected geometry path and perform sampled ownership
validation. When that validation finds a candidate defect, transport can
charge at most `1e-10` cm to the incident placement after checking raw
surface crossings and the complete ownership chain on both sides. The
replay JSON includes recovered gap and overlap counts and distance travelled.

The C navigator regressions cover a verified two-claimant overlap, a third
claimant, a real intervening void, a mirror, an over-width gap, a distance
limit and a collision inside a gap, and reversal and secondary restart
there. Fixed-source and adjoint tests compare identical seeds against a
reference model with the incident material extended to the far surface.
They check particle balance, track length, and adjoint response.

## End-to-end transport fixtures

The existing unit binaries accept a test name. Each invocation below runs one
fixed-seed fixture including preparation, transport, scoring and assertions:

```sh
make RELEASE=1 bin/tests/unit/test_nucdata bin/tests/unit/test_adjoint
bin/tests/unit/test_nucdata fixed_neutron_transport_preserves_optical_depth_across_cells_and_void
bin/tests/unit/test_nucdata fixed_photon_source_samples_photoatomic_collisions_and_deposition
bin/tests/unit/test_adjoint adjoint_vacuum_and_absorber_match_path_integral_for_both_species
```

| Fixture | Histories across all calls | Whole process (s) | Histories/s | Peak RSS (MB) |
|---|---:|---:|---:|---:|
| Fixed-source neutron | 40,000 | 0.273 | 147,000 | 4.5 |
| Fixed-source photon | 200 | 0.011 | 18,100 | 4.8 |
| Multigroup adjoint | 82,000 | 0.142 | 578,000 | 4.1 |

These single runs include process launch and fixture construction. The short
photon fixture is primarily a correctness check. The full nucdata suite also
checks secondary production and white reflection with sampled directions.
