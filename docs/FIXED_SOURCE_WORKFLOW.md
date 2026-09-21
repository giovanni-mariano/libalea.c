# Fixed-source neutron and photon workflow

The transport API runs independent neutron or photon source histories through
libalea geometry. Nuclear data is prepared explicitly for transport; geometry
loading, plotting, and validation do not require it. The source, bindings,
geometry, and tally plan remain alive during a run.

For an MCNP geometry and an xsdir file or directory of FENDL-style `.xsd`
files pointing to ACE data, build the libraries and example:

```sh
make full
make -C examples/c fusion_fixed_source
examples/c/fusion_fixed_source geometry.i xsdir-or-xsd-dir 100000 14.1 \
    -10 10 -10 10 -10 10 cell_scores.csv
```

The six coordinates define a uniform source box in centimetres; the box must
lie inside the modeled geometry. Every history starts with a unit-weight,
isotropically directed neutron at the specified energy in MeV. The example
uses a fixed seed and writes the mean neutron track length and its standard
error per source neutron for each cell. Add `--coupled` to request neutron
photon-production and photoatomic data, transport emitted photons, and include
their local-deposition score in the CSV. Data preparation fails explicitly
when required tables are unavailable. The example is a starting point for a
fusion calculation, not a built-in model of a plasma or source spectrum.

Python and Lua expose the same `transport_run(system, xsdir, config)` call for
small experiments. Preparation happens inside the call; reading or plotting a
geometry never needs an xsdir. Both return a dictionary/table with history
counters, neutron-only legacy `track_length` arrays, and a `tallies` sequence.
Each tally has `sum`, `sum_squared`, `mean`, and `standard_error` arrays,
`bin_ids` for cell/universe tallies, and mesh dimensions when relevant.
Python returns independent NumPy arrays for numeric results; Lua returns
tables. Arrays are zero-based in Python and one-based in Lua. Mesh scores are
flat and can be reshaped to `(nz, ny, nx)` for one energy group; with an energy
axis, reshape to `(energy_group_count, nz, ny, nx)`.

For new source configurations, specify the spatial and angular distributions
independently. `Source` in Python and `alea.source_prepare` in Lua create a
reusable prepared source. Preparation and preview do not need geometry or
nuclear data. The preview uses the same seed and global history IDs as a
transport run.

```python
source = pyalea.Source({
    "particle": "neutron",
    "space": {"type": "point", "position": [0, 0, 0]},
    "angle": {"type": "isotropic"},
    "energy": {"type": "mono", "value": 14.1},
    "time": {"type": "constant", "value": 0},
    "weight": 1,
})
samples = pyalea.sample_source(source, histories=1000, seed=1)
# NumPy arrays: position/direction (N,3); energy/time/weight/species/IDs (N,)
```

```lua
local source = alea.source_prepare({
    particle = "neutron",
    space = {type = "box", lower = {-1, -1, -1}, upper = {1, 1, 1}},
    angle = {type = "monodirectional", direction = {1, 0, 0}},
    energy = {type = "mono", value = 14.1},
})
local samples = alea.sample_source(source, 1000, 1)
-- position/direction are tables of triples; history_id starts at zero.
```

Pass `source` as the value of `config["source"]` in Python or
`config.source` in Lua. A source description dictionary/table may also be
passed directly. The current prepared source supports point, line, uniform box,
spherical-volume, and cylindrical-volume space; monodirectional, isotropic,
cone, cosine-hemisphere, tabulated-polar, or radial angle; constant energy,
time and weight; and neutron or photon identity. A line
uses `start` and `end`. A sphere uses `center`, `outer_radius`, and optional
`inner_radius`. A cylinder uses `base`, a base-to-top `axis` vector,
`outer_radius`, and optional `inner_radius`. Both radii describe uniform volume
sampling, including shells. `time` defaults to zero and `weight` to one.
Energy can also be a discrete line distribution, for example
`{"type": "lines", "values": [2.45, 14.1], "weights": [1, 3]}` in
Python, or `{type = "lines", values = {2.45, 14.1}, weights = {1, 3}}` in Lua.
Weights are nonnegative probability masses; the sampler normalizes them and
requires a positive total. Source history weight remains separate from line
selection probability.
Angular examples in Python are `{"type": "cone", "direction": [0, 0, 1],
"half_angle": 0.2}` (radians), `{"type": "cosine", "direction": [0, 0, 1]}`,
and `{"type": "tabulated_mu", "direction": [0, 0, 1], "mu": [0, 1],
"pdf": [0, 2], "interpolation": "linear"}`. The cone is uniform per unit
solid angle. The cosine law emits into the hemisphere around its axis.
Tabulated `pdf` is density per unit `mu = cos(theta)`, with uniform azimuth;
`interpolation` is `histogram` (left-bin height) or `linear`. A radial law,
`{"type": "radial", "origin": [0, 0, 0], "inward": False}`, points from the
origin toward each sampled position; `inward=True` reverses it. A sample
exactly at the origin has undefined radial direction and fails. Lua uses the
same field names in tables.
The C interface is in `alea_source.h`. A preview does not start transport or
read any ACE tables. Python returns `particle` as 0 for neutron and 1 for
photon; Lua uses the same codes.

```python
import pyalea

system = pyalea.load_mcnp("geometry.i")
xsdir = pyalea.XsDir("/data/xsdir")  # or XsDir("/data/xsd", directory=True)
result = pyalea.transport_run(system, xsdir, {
    "histories": 10000, "seed": 1, "coupled": True,
    "source": {"particle": "neutron", "energy": 14.1,
               "space": {"type": "box", "lower": [-10, -10, -10],
                         "upper": [10, 10, 10]},
               "angle": {"type": "isotropic"}},
    "tallies": [
        {"score": "track_length", "domain": "cell", "particle": "neutron"},
        {"score": "local_deposition", "domain": "cell", "particle": "photon"},
    ],
})
print(result["tallies"][0]["mean"])
```

```lua
local system = alea.load_mcnp("geometry.i")
local xsdir = alea.nuc_load_xsdir("/data/xsdir") -- or nuc_load_xsdir_dir
local result = alea.transport_run(system, xsdir, {
    histories = 10000, seed = 1, coupled = true,
    source = {particle = "neutron", energy = 14.1,
              space = {type = "box", lower = {-10, -10, -10},
                       upper = {10, 10, 10}},
              angle = {type = "isotropic"}},
    tallies = {
        {score = "track_length", domain = "cell", particle = "neutron"},
        {score = "local_deposition", domain = "cell", particle = "photon"},
    },
})
print(result.tallies[1].mean[1])
```

For a repeated point source, use `space = {type = "point", position = {0,0,0}}`
and `angle = {type = "monodirectional", direction = {1,0,0}}` in Lua (or the
equivalent Python dictionaries). Omitted `particle` defaults to `neutron`.
Tally `score` accepts `track_length`, `collision`,
`reaction_event`, `reaction_rate`, `heating`, and `local_deposition`;
`domain` accepts `cell`, `universe`, and `mesh`. Mesh tallies need `lower`,
`upper`, and `dimensions` arrays of length three. Optional tally filters are
`particle`, `material_id`, `reaction_mt`, `nuclide_zaid`, `energy_min/max`,
`time_min/max`, and `energy_edges`. Run options also accept `history_offset`,
`max_events_per_history`, `max_segment_distance`, and
`max_pending_particles`. Custom source callbacks remain available through C.

The geometry needs a vacuum boundary around its transport domain. If MCNP
graveyard conversion does not mark that boundary, pass
`--vacuum-surface ID` with the enclosing MCNP surface ID; the example sets
that surface to vacuum before preparing transport bindings.

The library entry point `alea_transport_run_sampled_source()` calls a source
sampler once for each history. `alea_source_sample()` adapts any prepared
built-in source to that callback. A custom callback can sample a plasma
shape, beam direction, energy spectrum, starting time, and statistical weight.
It receives the run seed and global history ID, so it can reproduce a given
history independently of batch size. Return an error if a valid source cannot
be sampled; the run reports the failed history rather than counting it as a
zero-score history. A sampled source particle must have a finite position,
positive finite energy and weight, finite time, and a nonzero finite direction.
The driver normalizes the direction.

Set `options.history_offset` to the first global history ID when splitting a
run. With the same seed and sampler, nonoverlapping batches reproduce the
same particles and collision streams as one run. Each tally result contains
the sum and sum of squares of complete source-history scores. Batch sums and
square sums can be added bin by bin; divide the combined sum by the combined
history count. For a view with `N > 1`, its standard error per source particle
is

```text
sqrt(max(0, (sum_squared - sum*sum/N) / (N-1)) / N)
```

Tallies are normalized per sampled source particle. To obtain a rate, multiply
the mean by the physical source rate in particles/s. Cell track length has
units cm per source particle; divide by a cell volume in cm³ to obtain
cell-average fluence in cm⁻² per source particle before multiplying by source
rate to obtain flux in cm⁻² s⁻¹. Photon local deposition is MeV per source
particle, not dose; conversion to heating density or dose requires the
appropriate volume or mass and unit conversion. See
[transport tally semantics](TRANSPORT_TALLIES.md) for filters, energy groups,
and the caveat about combining ACE neutron heating with transported-photon
deposition.

The current driver supports fixed-source neutron/photon transport, sampled
secondaries, vacuum and specular-reflective boundaries, and cell/mesh/terminal-
universe tallies. It does not sample URR probability tables or support
white/periodic boundaries. It does not solve a criticality eigenvalue problem.
