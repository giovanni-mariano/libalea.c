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
| Inelastic and neutron-multiplying reactions | Yes | Partial | No | Rejected with the active MT |
| Prompt or delayed fission emission | Partial | Mean yield only | No | Rejected |
| Unresolved-resonance probability tables | Yes | Standalone factors | No coordinated collision | Rejected |
| Free-gas target motion | No | No | No | Unsupported capability |
| Bound thermal scattering | No | No | No | Unsupported capability |
| Photon collisions and neutron-induced photon production | Partial | Cross sections only | No | Unsupported capability |

The accepted set is `ALEA_NUC_CAP_RESTRICTED_NEUTRON`, which combines
`ALEA_NUC_CAP_STATIONARY_ELASTIC` and `ALEA_NUC_CAP_ABSORPTION`. Preparation
also rejects active physics outside that set. It never replaces an unsupported
law with an approximate collision model.

## Transport contract

Prepared materials borrow the source `alea_nuc_material_t`, its component
array, and all referenced nuclides. Keep those objects alive and unchanged
until `alea_nuc_prepared_material_free()` returns. A prepared object is immutable
after construction and may be shared between workers. Evaluation and collision
results belong to the caller.

Energy is in MeV, distance is in cm, time is in seconds, microscopic cross
sections are in barns, number density is in atoms per barn-cm, and evaluated
macroscopic cross sections are in inverse cm. Directions must be unit vectors.

The RNG callback must return finite values in `[0,1)`. Sampling does not use
global RNG state and allocates no memory after preparation. A failed operation
does not publish an output value, though it may already have consumed values
from the caller's RNG. Flight and collision calls verify that their evaluation
still matches the stored particle state and the prepared material data.

Elastic collisions return both center-of-mass and laboratory cosines and a
three-dimensional outgoing direction. The elastic recoil energy is reported as
local deposition. For absorption, the neutron is terminated and eventwise
deposition is marked unavailable because ACE heating values are mean estimators,
not a sampled event energy balance.

## Example

Build the restricted fixed-source example after building the nuclear-data
module:

```sh
make cli
make -C examples/c slab_transport
examples/c/slab_transport /path/to/ace-directory 5 100000 1001.32c
```

The application owns slab traversal and history termination. The library owns
material evaluation, flight sampling, target and reaction selection, and the
elastic or absorption collision. The example exits with a capability diagnostic
when its table requires inelastic, fission, URR, thermal, or other unsupported
physics.
