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
    system.find_cells_grid_z(
        0.0, -3.0, 3.0, -3.0, 3.0, 16, 12,
        error_mode="fast", _as_buffers=True,
    )
    system.__init__()

    material = pyalea.NucMaterial()
    material.__init__()
    multigroup = pyalea.Multigroup([10.0, 1.0, 0.1])
    multigroup.__init__([20.0, 2.0])

gc.collect()
