# Fixed-source neutron workflow

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

The geometry needs a vacuum boundary around its transport domain. If MCNP
graveyard conversion does not mark that boundary, pass
`--vacuum-surface ID` with the enclosing MCNP surface ID; the example sets
that surface to vacuum before preparing transport bindings.

The library entry point `alea_transport_run_sampled_source()` calls a source
sampler once for each history. The supplied
`alea_transport_sample_box_isotropic()` samples a monoenergetic source
uniformly within an axis-aligned box. A custom callback can sample a plasma
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
