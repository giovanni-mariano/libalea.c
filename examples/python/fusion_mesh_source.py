#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""Check mesh-source transport against a void reference, then run an Fe shell.

Run with a neutron ACE directory containing 26056.32c, for example FENDL-3.2c:

    python3 examples/python/fusion_mesh_source.py /path/to/neutron/ace --histories 5000

The source is 14.1 MeV, isotropic, and distributed over two unequal mesh
voxels. Tallies are track length in cm per source neutron; they are not yet
divided by cell volume or multiplied by a physical source rate.
"""

import argparse
from pathlib import Path

import numpy as np
import pyalea


def geometry(iron: bool) -> str:
    shell = "2 1 -7.87 1 -2 imp:n=1" if iron else "2 0 1 -2 imp:n=1"
    material = "\nm1 26056.32c 1.0\n" if iron else ""
    return (
        "14.1 MeV mesh source in a spherical iron shell\n"
        "c\n"
        "1 0 -1 imp:n=1\n"
        f"{shell}\n"
        "99 0 2 imp:n=0\n"
        "\n"
        "1 so 2\n"
        "2 so 8\n"
        f"{material}"
    )


def exit_distance(position: np.ndarray, direction: np.ndarray, radius: float) -> np.ndarray:
    """Distance from each point inside a sphere to its forward intersection."""
    projection = np.einsum("ij,ij->i", position, direction)
    squared_radius = np.einsum("ij,ij->i", position, position)
    return -projection + np.sqrt(projection**2 + radius**2 - squared_radius)


def cell_means(result: dict) -> dict[int, tuple[float, float]]:
    tally = result["tallies"][0]
    return {
        int(cell_id): (float(mean), float(error))
        for cell_id, mean, error in zip(
            tally["bin_ids"], tally["mean"], tally["standard_error"]
        )
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("xsdir", type=Path, help="ACE directory or xsdir file containing 26056.32c")
    parser.add_argument("--histories", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=17)
    args = parser.parse_args()
    if args.histories < 1:
        parser.error("--histories must be positive")

    source = pyalea.Source({
        "particle": "neutron",
        "space": {
            "type": "cartesian_mesh",
            "x_edges": [-0.6, -0.2, 0.6],
            "y_edges": [-0.5, 0.5],
            "z_edges": [-0.5, 0.5],
            "values": [[[1.0]], [[2.0]]],
            "value_mode": "density",
        },
        "angle": {"type": "isotropic"},
        "energy": 14.1,
    })
    # The two voxel volumes are 0.4 and 0.8 cm3, respectively.
    np.testing.assert_allclose(source.integrated_emissivity, 0.4 * 1.0 + 0.8 * 2.0)

    xsdir = pyalea.XsDir(str(args.xsdir), directory=args.xsdir.is_dir())
    config = {
        "histories": args.histories,
        "seed": args.seed,
        "source": source,
        "tallies": [{"domain": "cell", "score": "track_length", "particle": "neutron"}],
    }

    # In void, each sampled neutron crosses both spherical surfaces once.
    samples = pyalea.sample_source(source, args.histories, args.seed)
    inner_exit = exit_distance(samples["position"], samples["direction"], 2.0)
    outer_exit = exit_distance(samples["position"], samples["direction"], 8.0)
    expected = {1: float(inner_exit.mean()), 2: float((outer_exit - inner_exit).mean())}
    void = pyalea.transport_run(pyalea.load_mcnp_string(geometry(False)), xsdir, config)
    void_means = cell_means(void)
    if void["collisions"] != 0 or void["leaked"] != args.histories:
        raise RuntimeError("void reference did not leak every source neutron")
    for cell_id in (1, 2):
        np.testing.assert_allclose(void_means[cell_id][0], expected[cell_id], rtol=1e-11, atol=1e-11)
        print(f"void cell {cell_id}: {void_means[cell_id][0]:.8f} cm/source "
              f"(ray calculation {expected[cell_id]:.8f})")

    iron = pyalea.transport_run(pyalea.load_mcnp_string(geometry(True)), xsdir, config)
    iron_means = cell_means(iron)
    if iron["histories"] != args.histories:
        raise RuntimeError("iron run returned an unexpected history count")
    if not all(np.isfinite(iron_means[cell_id]).all() for cell_id in (1, 2)):
        raise RuntimeError("iron run returned a non-finite tally")
    print(f"Fe-56: {iron['histories']} source histories, "
          f"{iron['collisions']} collisions, {iron['leaked']} leaked particles")
    for cell_id in (1, 2):
        mean, error = iron_means[cell_id]
        print(f"Fe-56 cell {cell_id}: {mean:.8f} +/- {error:.8f} cm/source")


if __name__ == "__main__":
    main()
