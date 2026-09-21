# Transport tallies

`alea_tally.h` defines an optional tally plan for a transport run. Create the
plan after geometry is built, add score specifications, and set
`alea_transport_options_t.tally_plan` before calling
`alea_transport_run_fixed_neutron()` or `alea_transport_run_fixed_source()`.
The geometry and plan must remain alive
and unchanged for the run. The result owns its tally arrays, so the plan may be
freed after the run.

```c
alea_tally_plan_t* plan = alea_tally_plan_create(system);
alea_tally_spec_t spec = {
    .domain = ALEA_TALLY_CELL,
    .score = ALEA_TALLY_REACTION_EVENT,
    .particle_mask = ALEA_TALLY_NEUTRON,
    .reaction_mt = 16,
    .nuclide_zaid = 3006
};
alea_tally_plan_add(plan, &spec, NULL);
options.tally_plan = plan;
/* Run transport, then inspect result.tallies through
 * alea_tally_results_view(result.tallies, 0, &view). */
```

The current scores are weighted track length in centimetres, weighted
collision count, weighted sampled-reaction count, a track-length
reaction-rate estimator, photon local deposition, and neutron heating. The neutron heating score adds
`weight * distance * Σ_H`, where `Σ_H` is the macroscopic ACE heating response
in MeV/cm; its sum is in MeV per source history and can be signed. The rate
score adds `weight * distance * Σ_mt`,
where `Σ_mt` is the evaluated macroscopic cross section in cm⁻¹ for the
selected reaction and target. It estimates reactions per source particle;
it is distinct from the count of sampled reaction events. For event and rate
scores, `reaction_mt == 0` includes all reactions. Each tally can filter by
particle type, incident energy, particle time, and material ID. Zero particle mask accepts
all types. Both energy bounds zero accept all energies; otherwise the interval
is `[energy_min, energy_max)`. Both time bounds zero accept all times;
otherwise the interval is `[time_min, time_max)` in seconds. Flight segments
are clipped to the time interval before spatial binning; collision events use
their arrival time. Material ID zero accepts all materials, `-1`
selects void, and a positive ID selects one material. `reaction_mt` applies
only to `ALEA_TALLY_REACTION_EVENT`, `ALEA_TALLY_REACTION_RATE`, and
`ALEA_TALLY_LOCAL_DEPOSITION`. Local deposition scores
`weight * local_energy_deposition` in MeV at photoatomic collisions, including
absorption and scattering. Collision, reaction-event, and local-deposition
tallies may filter the sampled target nuclide by
`nuclide_zaid` (`1000*Z+A`; zero accepts all targets). This distinguishes, for
example, Li-6 and Li-7 reactions in one material. Rate tallies apply the same
nuclide filter to the evaluated macroscopic cross section; they do not sample
a target. Heating tallies use the same target filter on the ACE heating
response and do not accept a reaction MT. If a selected prepared neutron
component has no heating data, a requested heating tally reports
`ALEA_ERR_NOT_FOUND` rather than silently scoring zero. The ACE heating
response is a neutron-only estimate. Its relationship to transported photon
deposition depends on how the ACE heating table was processed; do not add the
two scores as a total-heating estimate without checking that convention.

Spatial bins are either all terminal cells, all terminal universe IDs, or a
uniform Cartesian mesh in world coordinates. Cell bin order is geometry cell
index and `bin_ids` gives the corresponding cell IDs. Universe IDs are sorted
and repeated occurrences of a terminal universe contribute to one bin. Mesh
bins are half-open; the upper face is outside. The flat index is
`(z*ny + y)*nx + x`. Flight segments are clipped and divided at mesh faces,
including void segments. Collision, reaction, and local-deposition events
use their world position.

An optional energy axis can be set on any tally with `energy_edges` and
`energy_group_count`. Edges are strictly increasing, nonnegative MeV values;
the plan copies them when the tally is added. Group `g` accepts incident
energy in `[energy_edges[g], energy_edges[g+1])`. Values outside all groups
do not score. The result copies the edges and reports both
`spatial_bin_count` and total `bin_count`. Spatial bins are fastest:
`flat_index = energy_group * spatial_bin_count + spatial_bin`.
Without an energy axis, `energy_group_count` is zero and `bin_count` equals
`spatial_bin_count`. The ordinary `energy_min`/`energy_max` filter can also
be used as a prefilter on an energy-axis tally.

Each tally view contains `sum` and `sum_squared` for complete source-history
bin scores. All neutron and photon descendants of one source particle
contribute to the same history value before it is squared. For `N` histories,
the mean score per source particle is `sum/N`. The squared sums permit a sample variance and
standard error. Histories with no contribution to a bin contribute zero.

Photon track-length, collision, reaction-event, and local-deposition scores
are active for photon sources and neutron-induced photons. Reaction-rate and
ACE heating scores are currently neutron-only. Coupled runs need both neutron
and photon bindings; neutron-induced photon sampling also requires
`ALEA_NUC_CAP_PHOTON_PRODUCTION` in neutron preparation. A neutron-only
binding that encounters an emitted photon in material fails with
`ALEA_ERR_NOT_FOUND` at that history and position. Photon source histories
can use `alea_transport_run_fixed_source()`. The legacy per-cell path arrays
remain neutron-only; use tally plans for photon paths. Ancestor-universe
filters remain future work.
