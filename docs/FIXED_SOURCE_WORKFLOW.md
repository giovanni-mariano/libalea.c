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

For an analytic attenuation check with a layered shield, run
[`examples/python/sandwich_slab_transport.py`](../examples/python/sandwich_slab_transport.py)
with a neutron ACE directory containing FENDL-3.2c `26056.32c` and
`13027.32c`:

```sh
python3 examples/python/sandwich_slab_transport.py /path/to/neutron/ace
```

A 14.1 MeV pencil beam crosses vacuum, Fe-56, vacuum, Al-27, vacuum, Fe-56,
and a downstream 1 cm void detector. For each material slab, the script
computes its macroscopic total cross section as atomic density times the ACE
total cross section. If the survival probability on entry is `S`, the
uncollided track length is `S * (1 - exp(-Sigma*d)) / Sigma`; a void slab gives
`S*d`. The final survival probability is `exp(-sum(Sigma*d))`, equal to the
uncollided track length in the downstream 1 cm detector. A narrow source-energy
filter selects this first-flight contribution because collisions almost surely
change the neutron energy. The script first checks exact slab lengths with all
materials replaced by void, then checks each material-run tally within five
reported standard errors of the analytic result. It also prints the unfiltered
track length, which includes scattered neutrons and has no simple exponential
formula. The default run uses 20,000 source histories and a fixed seed.

For lower-energy neutron transport in the unresolved-resonance range, run
[`examples/python/urr_u238_transport.py`](../examples/python/urr_u238_transport.py)
with the Lib80x xsdir containing `92238.00c`:

```sh
python3 examples/python/urr_u238_transport.py /path/to/Lib80x/xsdir
```

This sends 50 keV neutrons through a U-238 sphere, then repeats the same
calculation with the material divided into two adjacent cells. The second
geometry has no physical interface: collision and leakage counts, total
track length, and MT=102 capture-rate tally must match for the fixed seed.
Transport samples one coordinated URR probability-table realization for each
neutron material/flight encounter. It retains that realization when crossing
adjacent cells of the same material or a distance-limit segment, and draws a
new one after a collision or on entering a different material. Nuclear-data
tables without URR data use their ordinary continuous-energy cross sections.

For a reproducible neutron calculation with the Python binding, run
[`examples/python/fusion_mesh_source.py`](../examples/python/fusion_mesh_source.py)
with a neutron ACE directory containing FENDL-3.2c `26056.32c`:

```sh
python3 examples/python/fusion_mesh_source.py /path/to/neutron/ace --histories 5000
```

The script sends 14.1 MeV neutrons from a weighted Cartesian source mesh into
a 2–8 cm Fe-56 spherical shell. It first replaces the iron with void and
checks both cell track-length tallies against independently calculated ray
intersections. It then restores Fe-56 at 7.87 g/cm³ and reports collisions,
leakage, and cell track lengths with standard errors. The check uses a fixed
seed and requires NumPy and an importable `pyalea` package. A leaked-particle
count above the source-history count can occur when neutron reactions produce
additional particles. Track length is in cm per source neutron; divide by
cell volume to obtain volume-averaged flux per source neutron.

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
spherical-volume, cylindrical-volume, tokamak RZ, and Cartesian-mesh space;
monodirectional, isotropic, cone, cosine-hemisphere, tabulated-polar, or radial
angle; constant energy, time and weight; and neutron or photon identity. A line
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
For a continuous spectrum, use `{"type": "tabulated", "values": [1, 3],
"pdf": [1, 1], "interpolation": "linear"}`. `values` are increasing energy
knots in MeV; `pdf` is density per MeV. `histogram` uses each left-bin height,
while `linear` interpolates between adjacent densities. The table is
normalized when the source is prepared.
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
For an axisymmetric plasma, `space` can be `tokamak_rz`. It takes increasing
`r_edges` (nonnegative cm), increasing `z_edges` (cm), and a two-dimensional
`emissivity` array with one row per R bin and one value per Z bin. Each value
is a nonnegative emission density per cm³ within its bin. For example:

```python
source = pyalea.Source({
    "space": {"type": "tokamak_rz", "r_edges": [100, 110, 120],
              "z_edges": [-10, 0, 10],
              "emissivity": [[1, 0], [2, 1]]},
    "angle": {"type": "isotropic"}, "energy": 14.1,
})
rate = source.integrated_emissivity
```

The default toroidal range is a full turn. `phi_min` and `phi_max` can select
a sector in radians, with a span at most 2π. Bin probabilities follow
emissivity times cylindrical volume; radius is sampled uniformly in R² within
the selected bin. Zero-emissivity bins are never selected. The sampler
normalizes probabilities without modifying particle history weight.
`integrated_emissivity` is the sector integral in the input units (particles/s
if emissivity is particles/cm³/s). Multiply a per-source tally mean by this
rate to get a physical rate. In Lua, call
`alea.source_integrated_emissivity(source)` on a prepared source; C uses
`alea_source_integrated_emissivity()`.
For non-axisymmetric imported emission, use a rectilinear Cartesian source
mesh. Supply `x_edges`, `y_edges`, and `z_edges` in cm, plus a three-dimensional
`values` array indexed `[x_bin][y_bin][z_bin]` (Z varies fastest in C):

```python
source = pyalea.Source({
    "space": {"type": "cartesian_mesh",
              "x_edges": [0, 1, 3], "y_edges": [0, 1],
              "z_edges": [0, 1, 2],
              "values": [[[1, 0]], [[1, 1]]],
              "value_mode": "density"},
    "angle": {"type": "isotropic"}, "energy": 14.1,
})
```

`value_mode` is required. With `density`, each voxel's selection mass is its
value times its volume; values have units per cm³. With `strength`, values
are already integrated voxel strengths and are used directly. Within a
selected voxel, position is uniform in X, Y, and Z. Zero-strength voxels are
never sampled. `integrated_emissivity` returns the sum of voxel masses in the
input units. The corresponding Lua tables use the same keys and nesting.
Multiple complete sources can be mixed with physical strengths:

```python
mixed = pyalea.Source({
    "type": "mixture",
    "components": [
        {"strength": 3, "source": plasma_description},
        {"strength": 1, "source": calibration_source_description},
    ],
})
```

Component selection follows the normalized strengths. If strengths are
absolute emission rates, their sum is the mixture's total physical rate;
otherwise only their ratios matter. Selection does not multiply the history
weight. Components may differ in particle species and can themselves be
mixtures. Lua uses the same table layout. In C,
`alea_source_mixture_prepare()` takes ownership of prepared components on
success and leaves ownership with the caller on failure.
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
`max_events_per_history`, `max_segment_distance`, `max_pending_particles`, and
`boundary_distance_tolerance`. The last option defaults to `1e-6` when omitted;
zero disables its user-sized proximity window. It applies only when competing
hierarchy crossings are already proven aliases of the same world surface.
It does not merge distinct close surfaces or define a minimum material-layer
thickness. Custom source callbacks remain available through C.

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
secondaries, vacuum, specular-reflective, white, and translational periodic
plane boundaries, plus cell/mesh/terminal-universe tallies and neutron URR
probability-table sampling. Periodic boundaries require a distinct parallel
plane partner at the root geometry level. Rotational periodic mappings and
criticality eigenvalue calculations remain unsupported.
