# Nuclear-data transport capabilities

The continuous-energy collision interface is intentionally capability-gated.
Loading a table makes its data available for inspection; it does not imply that
every decoded reaction can be transported. Call `alea_nuc_capabilities()` for a
single nuclide or `alea_nuc_prepare_material()` for a material before transport.

## Current capability matrix

| Data or operation | Read | Evaluate | Sample | Preparation behavior |
|---|---:|---:|---:|---|
| Pointwise neutron total, elastic and absorption cross sections | Yes | Yes | — | Validated for finite, nonnegative and consistent values |
| Stationary-target elastic, isotropic angular data | Yes | Yes | Yes | Accepted |
| Stationary-target elastic, 32-bin equiprobable angular data | Yes | Yes | Yes | Accepted after grid validation |
| Stationary-target elastic, histogram or linear tabular angular data | Yes | Yes | Yes | Accepted after PDF/CDF validation |
| Non-fission absorption with no outgoing neutron | Yes | Yes | Yes | Accepted when event channels reproduce absorption |
| Level, tabular, Maxwell, evaporation and Watt neutron emission | Yes | Yes | Yes | Accepted after interpolation/PDF validation |
| Kalbach-Mann and energy-first correlated emission | Yes | Yes | Yes | Accepted; law 61 nested angular tables are decoded and validated |
| Laboratory angle-first correlated emission | Yes | Yes | Yes | Law 67 equiprobable cosine grids and conditional spectra are decoded and sampled in their defined order |
| N-body phase-space emission | Yes | Yes | Yes | Accepted for three to five products |
| Prompt fission emission | Yes | Yes | Yes | Accepted with caller-owned secondary storage |
| Delayed fission emission | Yes | Yes | Yes | Opt-in capability; group spectra and exponential emission times are sampled |
| Unresolved-resonance probability tables | Yes | Yes | Yes | Accepted through coordinated evaluation workspace |
| Free-gas target motion | Table kT | Yes | Yes | Opt-in constant-cross-section model; every component requires positive kT |
| Exact-temperature table selection | Yes | — | — | Selects a unique closest table within a caller-supplied absolute kT tolerance |
| Bound thermal scattering, discrete IFENG=0/1 | Yes | Yes | Yes | Explicit component association; ZAID/kT validated; replaces free-atom elastic below the table cutoff |
| Bound thermal scattering, continuous IFENG=2 | Yes | Yes | Yes | Variable outgoing grids, linear PDF inversion, correlated cosine sampling, and prepared association |
| Photoatomic coherent and incoherent scattering | Yes | Yes | Yes | EPR Compton profiles Doppler-broaden scattered energy; accepted for homogeneous prepared photon-table mixtures |
| Photoelectric and pair production | Yes | Yes | Yes | Detailed EPR relaxation photons and pair-event 511 keV photons can be banked; charged-particle energy is local |
| EPR subshell and atomic-relaxation data | Yes | Yes | Yes | Subshell cross sections select the initial vacancy; radiative and nonradiative transition cascades are sampled |
| Averaged photoatomic fluorescence metadata | Yes | No | No | Older JFLO edge, phi, cumulative yield, and representative energy arrays are retained for inspection |
| Neutron-induced photon production | Yes | Yes | Yes | Opt-in capability; MF=12/13/16 channels, GPD totals, YP multipliers, aggregate parents, and discrete/continuum law-4 spectra are decoded |

`ALEA_NUC_CAP_RESTRICTED_NEUTRON` combines stationary elastic scattering and
absorption. `ALEA_NUC_CAP_CONTINUOUS_NEUTRON` adds neutron emission, while
`ALEA_NUC_CAP_FISSION`, `ALEA_NUC_CAP_DELAYED_NEUTRON`, and
`ALEA_NUC_CAP_URR` request those features explicitly. Add
`ALEA_NUC_CAP_FREE_GAS` to use target motion for elastic collisions at or below
`400 kT`; higher-energy collisions retain the stationary-target treatment.
Add `ALEA_NUC_CAP_THERMAL_SAB` and provide component/table associations in the
preparation requirements to use bound thermal scattering. Associated thermal
tables take precedence over the free-gas option inside their energy range.
Preparation reports the active MT and law when a distribution is invalid or
unsupported.

EPRDATA12, EPRDATA14, and EPRDATA25 photoatomic layouts are supported. Repeated
energy knots at shell thresholds are interpreted as right-continuous
discontinuities. Evaluated nonradiative transition energies remain inspection
data because their charged particles are deposited locally; photon transition
energies must be finite and nonnegative. The Doppler-shell population, binding,
probability, momentum, PDF, and CDF blocks are retained explicitly. EPRDATA12
differential shell probabilities and the cumulative EPRDATA14/25 representation
are normalized to a common cumulative form during decoding.

The constant-cross-section free-gas kernel has a statistical equilibrium gate.
It checks the expected collision-rate-biased spectrum and recovers the
Maxwellian mean energy after applying the analytic mean-relative-speed
residence weighting.

`ALEA_NUC_CAP_PHOTON_PRODUCTION` enables photons emitted by the selected
neutron reaction. MF=12 and MF=16 multiplicities are evaluated directly;
MF=13 production cross sections are divided by the corresponding parent
reaction cross section. Each yield is stochastically rounded, and its decoded
energy and angular distribution is sampled into the same caller-owned
secondary buffer as emitted neutrons. Discrete primary-photon energies include
the incident-energy recoil term specified by ACE law 2. Aggregate MT=3, MT=4,
and MT=18 parents are applied to their sampled child reactions while their
yield denominator is reconstructed from the aggregate cross section. GPD total
production cross sections are available through
`alea_nuc_xs_photon_production_total()`, and YP multiplier references are
validated while decoding. Charged-particle level parents MT=600--849 are
conditioned on the selected MT=103--107 aggregate with the partial-to-aggregate
cross-section ratio. When an inclusive MT=5 channel replaces an omitted
neutron/charged-particle reaction, its MF=13 photon-production cross section is
conditioned on that MT=5 event. Other parents that cannot be associated with a
selectable event reaction fail preparation.

Law-4 photon spectra accept ACE incident interpolation codes 11/12 and 21/22.
Pure discrete outgoing tables may carry interpolation digit zero because no
continuous interpolation is needed. In mixed spectra, discrete line energies
are interpolated by line index and unit-base scaling uses only the continuous
energy interval.

The required Linux CI gate downloads checksum-pinned U-235 and Ni-61 ACE
fixtures. It checks MF=12, MF=13, and MF=16 decoding, discrete and continuous
law-2/law-4 spectra, aggregate reaction parents, law-67 joint neutron moments,
sampled photon spectrum and angular moments, and deterministic replay through
the core event RNG. Derived MT 201–207 production totals, MT 301 heating, and
MT 444 damage-energy responses remain available as cross sections but are
excluded from the analogue event-channel map.

The production-library anisotropy gate uses natural carbon `6000.71c`. It
integrates each equiprobable cosine bin analytically and compares the resulting
first and second moments with sampled MT=51 photon directions. The gate
requires `ALEA_PRODUCTION_NEUTRON_XSDIR`; an absent production library may skip
the ordinary test suite, while the explicit required target rejects a missing
or invalid xsdir.

`nuc_inventory` reports neutron collision readiness from restricted-neutron
material preparation, independently of optional photon-production readiness.
It still retains the first optional photon-production issue in the row so a
production-library scan distinguishes a usable neutron table from complete
coupled neutron/photon support. ACE main, angular, applicability, and delayed
energy grids may repeat knots to encode right-continuous discontinuities.
Zero-probability trailing spectrum points are ignored after the CDF reaches
one. Small negative roundoff in photon yields and Kalbach precompound fractions
is clamped at use. URR collision cross sections are clamped to zero if a
processed table contains a negative value; signed heating is retained.

Use `alea_nuc_xsdir_find_temperature()` before loading or associating data at a
requested temperature. The supplied ZAID identifies a table family and type;
its numeric library suffix is ignored during the search. Selection succeeds
only inside the caller's absolute kT tolerance and rejects equidistant matches.
`alea_nuc_xsdir_find_temperature_bracket()` instead returns bounding evaluated
tables without extrapolation. Load both and pass them to
`alea_nuc_material_add_temperature_mix()` to obtain linear expected cross
sections while collision selection retains the sampled table's distributions
and URR representation. Coordinated URR evaluation uses one probability
quantile for both bounding tables so their resonance fluctuations remain
correlated. No table is mutated, and exact endpoints use one component.
For unresolved-resonance transport, the ACE ILF and IOA flags determine which
smooth inelastic and other-absorption competition channels are added to the
sampled elastic, fission, and capture values. A negative flag suppresses that
category, zero sums its applicable smooth reactions, and a positive flag names
either a summed competition tabulation or the sole competition reaction.
Aggregate MT=4 rates are distributed over the sampleable inelastic channels
in proportion to their smooth cross sections, as are aggregate charged-particle
absorption rates when level reactions are present. Preparation validates named
tabulations over the URR range. The probability-table total remains diagnostic
because transport reconstructs an exactly normalized total from the selected
partial channels.
In-place Doppler broadening remains a separate caller-owned option;
its numerical gates cover the 1/v invariant and a narrow-resonance reference
quadrature. It rejects a nuclide with an attached URR probability table before
allocation or mutation because broadening only the smooth cross sections would
mix data processed at different temperatures.

`alea_nuc_load_thermal()` decodes the discrete and continuous correlated ACE
representations for incoherent inelastic scattering and coherent, incoherent,
or mixed elastic scattering. `alea_nuc_sample_thermal_collision()` samples
discrete energy-angle quantiles, continuous outgoing-energy PDF/CDF tables,
correlated cosine sets, and Bragg edges without allocation. The table retains
the ACE applicability identifiers. Preparation accepts exact hydrogen and iron
isotopes; for other elements, it applies the ACE identifier to naturally
occurring isotopes, consistent with the established ACE-reader convention. It
checks kT within the requested absolute tolerance. The prepared path
reports thermal cross section separately and replaces the associated
component's free-atom elastic contribution through the thermal table cutoff.

The required Linux continuous-thermal gate uses a checksum-pinned H-in-ZrH
IFENG=2 table. It checks fixed-draw sampling and statistical outgoing-energy,
angular, and upscatter moments against independently integrated reference
values. It also checks prepared-material macroscopic replacement and the
elastic/inelastic event ratio for an H-in-ZrH moderator.
The same required gate covers 300 K light water with a discrete skewed
outgoing-energy law, comparing sampled energy, angle, and upscatter moments
with direct integration of its ACE probabilities and cosine sets.

## Transport contract

Prepared materials borrow the source `alea_nuc_material_t`, its component
array, all referenced nuclides, and associated thermal tables. Keep those
objects alive and unchanged until `alea_nuc_prepared_material_free()` returns.
A prepared object is immutable
after construction and may be shared between workers. Evaluation and collision
results belong to the caller. Emitting collisions append particles to an
`alea_nuc_secondary_buffer_t`; its storage and capacity belong to the caller.
Capacity is checked before RNG consumption and particles are never truncated.
Delayed neutrons are appended to the same buffer with their sampled future
emission time in `particle.time`; the application decides when to transport
them.

Energy is in MeV, distance is in cm, time is in seconds, microscopic cross
sections are in barns, number density is in atoms per barn-cm, and evaluated
macroscopic cross sections are in inverse cm. Directions must be unit vectors.

The RNG callback must return finite values in `[0,1)`. Sampling does not use
global RNG state. The supplied `alea_nuc_rng_uniform` callback uses libalea's
`src/rng` Philox4x32-10 implementation. Initialize its caller-owned
`alea_nuc_rng_t` with `alea_nuc_rng_init`, providing a seed, history ID,
particle ordinal, event index, and flight, collision, or URR domain. Reuse the
state throughout one sampling operation; initialize a new event for the next
operation. Assign particle ordinals deterministically within each history.
Each event allows 2^32 32-bit words (two words per uniform); exhaustion returns
NaN, which the samplers reject. Applications may still supply their own RNG.
An automated gate compares 257 elastic histories run serially and under dynamic
worker scheduling. Each history initializes the same event-addressed RNG state
in both runs, and every published collision field must match bit for bit.
Sampling allocates no memory after preparation; a unit gate arms the internal
allocation hook across evaluation and a multi-neutron collision and requires
zero allocator calls. A failed operation
does not publish an output value, though it may already have consumed values
from the caller's RNG. Flight and collision calls verify that their evaluation
still matches the stored particle state and the prepared material data.

Decoder and construction allocations pass through a module-local layer. Unit
tests inject a failure at each allocation in complete synthetic neutron and
photoatomic ACE loads and verify full cleanup, rejection of incomplete
nuclides, and successful reuse after injection is disabled. Tests also sweep
neutron-induced photon-production yields, channel storage,
non-isotropic angular tables, and nested energy distributions. A failed nested
energy allocation invalidates the decode instead of publishing a partial
distribution. The injection hook is internal and reserved for single-threaded
tests.

Use `alea_nuc_evaluate_urr()` for energies inside an unresolved resonance
range. It samples one probability-table realization per component into a
caller-owned `alea_nuc_evaluation_workspace_t`. Keep that workspace alive and
unchanged through flight and collision sampling. The ordinary
`alea_nuc_evaluate()` call refuses an energy that requires URR state. Sampled
elastic, fission, and capture values are used consistently for the macroscopic
flight rate, target selection, and event selection. Preserve the evaluation
and workspace when an unchanged-energy flight crosses a geometry boundary but
remains in the same material. Evaluate the destination material into its own
workspace after a material boundary, and reevaluate after any collision that
changes the neutron energy. The macroscopic total is the sum of sampled event
partials; the independently rounded probability-table total is retained as a
consistency check with a 0.1% relative tolerance.

Elastic collisions return both center-of-mass and laboratory cosines and a
three-dimensional outgoing direction. The elastic recoil energy is reported as
local deposition. For absorption, the neutron is terminated and eventwise
deposition is marked unavailable because ACE heating values are mean estimators,
not a sampled event energy balance.

Free-gas elastic scattering samples the collision-conditioned Maxwellian target
velocity and performs exact two-body velocity kinematics. This implementation
assumes the elastic cross section is constant over the target-motion energy
shift. It does not implement Doppler broadening rejection correction (DBRC), so
preparation does not claim a resonance-dependent target-motion capability.

`ALEA_NUC_CAP_PHOTON` prepares materials made entirely from photoatomic element
tables. Evaluation computes the mixture attenuation rate, and collision
sampling selects an element by its macroscopic contribution before sampling
the interaction. Neutron and photon tables cannot be mixed in one prepared
material.

`alea_nuc_sample_photon_collision()` also samples photoatomic channels directly
from one element table. Rayleigh angles use the integrated coherent form factor and
Thomson rejection. Compton events use Klein-Nishina rejection followed by the
incoherent scattering function. EPR tables then select a bound-electron shell,
sample its longitudinal momentum from the Compton profile, and Doppler-broaden
the scattered energy subject to its binding energy. For EPR tables,
photoelectric sampling selects an initial vacancy
from the subshell cross sections and follows its radiative, Auger, and
Coster-Kronig transition cascade. Radiative photons are banked through the
secondary-buffer API; photoelectrons and nonradiative relaxation electrons
remain local. Compton-induced vacancies use the same cascade and may also bank
radiative photons. The sampler deposits the incident energy minus transported
photon energy, preserving eventwise energy balance. Pair events bank two
back-to-back annihilation photons through the same API. Older JFLO data retain their averaged
metadata and use full local deposition because they cannot define an individual
shell cascade. Explicit electron and positron transport remains outside this
increment.

## Example

Build the restricted fixed-source example after building the nuclear-data
module:

```sh
make cli
make -C examples/c slab_transport
examples/c/slab_transport /path/to/ace-directory 5 100000 1001.32c 42
```

The application owns slab traversal and history termination. The library owns
material evaluation, flight sampling, target and reaction selection, and the
elastic or absorption collision. The example exits with a capability diagnostic
when its requested restricted capability cannot represent the loaded table.
It remains an elastic/absorption example and does not bank secondaries.
