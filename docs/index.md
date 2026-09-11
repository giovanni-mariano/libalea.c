<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO

SPDX-License-Identifier: MPL-2.0
-->

# ALEA

ALEA is a C library for building, debugging, and analyzing Constructive Solid
Geometry (CSG) models used in neutron and gamma transport simulations.

It can load MCNP and OpenMC geometries, query cells and materials, detect
geometry errors, trace rays, generate two- and three-dimensional
visualizations, export meshes, and work with ACE nuclear data.

!!! warning "Active development"

    ALEA is under active development. Public APIs may still change.

## Where to begin

<div class="grid cards" markdown>

-   :material-language-c:{ .lg .middle } **Use the C library**

    ---

    Build and query a geometry, then export the result.

    [:octicons-arrow-right-24: C tutorial](TUTORIAL.md)

-   :material-language-lua:{ .lg .middle } **Automate with Lua**

    ---

    Use the embedded CLI scripting interface for exploratory workflows.

    [:octicons-arrow-right-24: Lua tutorial](LUA_TUTORIAL.md)

-   :material-shape:{ .lg .middle } **Understand the model**

    ---

    Learn about surfaces, senses, cells, universes, lattices, and materials.

    [:octicons-arrow-right-24: Core concepts](CONCEPTS.md)

-   :material-api:{ .lg .middle } **Look up an API**

    ---

    Browse the public C functions grouped by task.

    [:octicons-arrow-right-24: API reference](API.md)

</div>

## Main capabilities

- Load and export MCNP and OpenMC geometry models.
- Build CSG geometry programmatically with boolean operations.
- Query points, overlaps, hierarchy paths, and material assignments.
- Trace rays through cells, fills, universes, and lattices.
- Produce analytical slices, sampled meshes, and 3D renders.
- Generate void regions and validate geometry transitions.
- Read ACE nuclear data and calculate microscopic, macroscopic, and
  multigroup cross sections.

## Minimal C example

```c
#include <alea.h>
#include <alea_mcnp.h>
#include <stdio.h>

int main(void) {
    mcnp_model_t *model = mcnp_load("geometry.inp");
    if (!model) {
        fprintf(stderr, "load failed: %s\n", alea_error());
        return 1;
    }

    int cell = 0;
    int material = 0;
    if (alea_find_cell_at(model->sys, 0.0, 0.0, 0.0,
                          &cell, &material) == 0) {
        printf("cell %d, material %d\n", cell, material);
    }

    mcnp_model_destroy(model);
    return 0;
}
```

The [C tutorial](TUTORIAL.md) covers compilation, ownership, geometry queries,
visualization, ray tracing, and export in more detail.
