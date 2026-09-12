# External nuclear-data tests

The broad external-library tests are optional when their evaluated-data
libraries are absent. The pinned moderator fixtures are required by
`make test-thermal-required`.

The thermal fixture is the NNDC ENDF/B-VII.1 processed thermal-scattering ACE
archive:

- Source: `ENDF-B-VII.1-tsl.tar.gz`
- SHA-256: `601db93498b51aa5a44800ee7ae7d439ac2bfb3ad36d7d1eb7e4a7be28968bc7`

Reference values in `test_thermal.c` were independently decoded from the pinned
ACE fixture. The required graphite benchmark checks coherent Bragg-edge and
inelastic channel probabilities plus sampled energy and angular moments against
direct integration of the decoded tables.

Continuous `IFENG=2` decoding and sampling are checked with NJOY2016.79 test
74, an H-in-ZrH table generated specifically to exercise that representation:

- Source: `tests/74/referenceTape71` from the NJOY2016.79 release
- SHA-256: `624147b06fc517cc7a5227ceac834a3febf107d00c1c5e86ecdddcfea51cd2d2`

The fixed-draw continuous collision result was independently calculated from
the ACE probability tables and sampling equations.

Run `make test-thermal-required` to download this table, verify its hash, and
run the non-skipping continuous-thermal CI gate. In addition to the decoded
values and fixed-draw collision, the gate checks the sampled outgoing-energy
mean, angular mean, and upscatter probability against independently integrated
values from the pinned probability and cosine tables. A prepared H-in-ZrH
moderator check also verifies macroscopic thermal replacement and the sampled
elastic/inelastic channel ratio.

The required thermal gate also downloads the 1.0 MB 300 K light-water ACE
table from NJOY2016.79 test 25:

- Source: `tests/25/referenceTape70` from the NJOY2016.79 release
- SHA-256: `0e2bcca8e078e55d135391904a3cf851e3170f6b6d8907fb7e38964f75996e98`

Its discrete skewed outgoing-energy law is checked against analytically
integrated bin probabilities and cosine means. Prepared-material checks verify
replacement of the free-hydrogen elastic cross section by the light-water
inelastic channel.

Neutron-induced photon production is checked with the 1.1 MB U-235 ACE table
from NJOY2016.79 test 07:

- Source: `tests/07/referenceTape26` from the NJOY2016.79 release
- SHA-256: `4459a2c249f064888aa7d38ecffbcae39ebea4170678ae412d0b71a4592d84c8`

Run `make data-njoy-test07 test_photon_production`. The GPD totals, YP
multipliers, repeated-knot MF=12 yields, MF=13 aggregate yield, and a discrete
photon line were independently decoded from the pinned ACE fixture.

MF=16 multiplicities and continuous photon spectra are checked with the 1.8 MB
Ni-61 ACE table from NJOY2016.79 test 08:

- Source: `tests/08/referenceTape25` from the NJOY2016.79 release
- SHA-256: `b969e801a59924f526d4f032540dcb3a9ad912214bf4bf4a7c47ec00bdf4c56b`

Run `make test-photon-production-required` to download both compact photon
fixtures, verify their hashes, and run the non-skipping validation gate. The
Ni-61 checks cover five MF=16 channels, independently decoded multiplicities,
the sampled mean of a continuous law-4 spectrum, isotropic angular statistics,
and exact replay with libalea's event RNG. Unit tests cover code-22 unit-base
interpolation and mixed discrete/continuum handling. The gate also exercises
deterministic Lua access to discrete and continuous photon spectra, angular
samples, and a nonelastic outgoing-energy distribution.

The same Ni-61 table contains neutron law 67. Required tests decode its
angle-first conditional spectra, compare sampled joint energy-angle moments
with direct PDF integration, and require the complete table to pass collision
capability inspection. Its derived MT 444 damage response is retained for
cross-section access and excluded from analogue reaction selection.

Production-library photon anisotropy is checked separately with the natural
carbon `6000.71c` table from the official ENDF/B-VII.1 neutron library. Set
`ALEA_PRODUCTION_NEUTRON_XSDIR` to the extracted library's `xsdir` and run:

```sh
make test-production-photon-angles-required \
  ALEA_PRODUCTION_NEUTRON_XSDIR=/path/to/ENDF-B-VII.1-neutron-293.6K/xsdir
```

The gate verifies the MT=51 MF=12 channel's 72 angular tables and compares
131072 sampled photon directions with independently integrated first and
second cosine moments for its 32-bin equiprobable representation. The source
table SHA-256 is
`1b93295279fa1ea1ea84b120d22e330f7f3f641a3eeb03be1fc880e6c63ddcf3`.

Production unresolved-resonance transport is checked with `92238.00c` from
the official ENDF/B-VIII.0 Lib80 library. Set `ALEA_ENDFB80_XSDIR` to its
extracted `xsdir` and run:

```sh
make test-production-urr-required \
  ALEA_ENDFB80_XSDIR=/path/to/Lib80x/xsdir
```

The gate pins the U-238 probability-table metadata and factors at 50 keV,
checks coordinated macroscopic reaction rates, and compares transmission
through two consecutive same-material geometry segments with the analytic
attenuation. It also verifies that entering a material with a different number
density requires a new evaluation and scales the sampled macroscopic rate.
The source table SHA-256 is
`35f1d10d8248395e8390e0322409a55c76d788197f69134d4335913f06626c27`.

The same pinned U-235 evaluation supplies a resonance-bearing table for the
required temperature-mixing gate. `make test-temperature-mix-required` creates
an independently Doppler-broadened upper table, mixes the two immutable tables,
checks macroscopic cross sections at thermal, resonance, and fast energies, and
verifies the table-selection frequency against the cross-section-weighted
probability.

Photoatomic fluorescence metadata is checked with the 1.45 MB uranium ACE
table from NJOY2016.79 test 59:

- Source: `tests/59/referenceTape27` from the NJOY2016.79 release
- SHA-256: `c3a6e20ad375227a65b861c9bdee6e70fbde34f23b484721d7054e92cd7b4127`

Run `make test-photoatomic-relaxation-required` to download the fixture and
validate all edge, phi, cumulative-yield, and representative-energy values in
its averaged fluorescence block.

The EPRDATA12, EPRDATA14, and EPRDATA25 tests validate explicit form-factor grids,
bound-electron Compton profiles and Doppler-broadened events, detailed subshell
cross sections, atomic-transition networks, sampled relaxation cascades, and
eventwise photon-energy balance. `make data-eprdata12`, `make data-eprdata14`, and
`make data-eprdata25` download the official archives and verify their published
SHA-512 digests before extraction. The newer library also exercises repeated
right-continuous shell-threshold knots and locally deposited nonradiative
transitions whose evaluated electron energy is negative.

`bin/nuc_inventory XSDIR [ZAID ...]` emits a tab-separated capability matrix
for an evaluated library. Each row distinguishes successful decoding and
cross-section access from collision readiness, lists the supported sampling
capabilities, and reports the first blocking MT or distribution law. Build it
with `make -C tools ../bin/nuc_inventory` after building the libraries.
