# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""Repeat ownership-sensitive binding calls for the Valgrind CI check."""

import gc

import pyalea


def populated_system():
    system = pyalea.System()
    _, _, inside = system.sphere_surface(1, 0.0, 0.0, 0.0, 2.0)
    material = system.add_material(7)
    system.add_cell(101, inside, material, 1.0)
    system.build_universe_index()
    return system


for _ in range(250):
    system = populated_system()
    system.get_config()
    assert system.add_transform(17, (1.0, 2.0, 3.0)) == 17
    inline_id = system.add_inline_transform(
        (1.0, 2.0, 3.0), cell_id=101, role="fill",
    )
    assert system.add_inline_transform((1.0, 2.0, 3.0)) == inline_id
    system.export_mcnp_string()
    system.export_openmc_string()
    system.export_serpent_string()
    system.find_cells_grid_z(
        0.0, -3.0, 3.0, -3.0, 3.0, 16, 12,
        error_mode="fast", _as_buffers=True,
    )
    system.__init__()
    with system.slice_error_query(
        (0, 0, 0), (0, 0, 1), (0, 1, 0),
        (-1, 1, -1, 1), (-1, 1, -1, 1),
    ) as query:
        assert query.run_page(0)["receipt"]["scope_classified"]

    material = pyalea.NucMaterial()
    material.__init__()
    multigroup = pyalea.Multigroup([10.0, 1.0, 0.1])
    multigroup.__init__([20.0, 2.0])

gc.collect()
