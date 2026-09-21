#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""Compare transport through vacuum/Fe/vacuum/Al/vacuum/Fe/vacuum slabs.

The 14.1 MeV pencil beam starts at x=-0.5 cm and travels along +x. Positive
MCNP densities are atomic densities in atoms/barn-cm. The narrow energy tally
selects the uncollided component: a neutron collision changes its energy
almost surely. The last 1 cm void slab measures uncollided transmission;
the full-energy tally also reports all scattered contributions.

Run with a neutron ACE directory containing FENDL-3.2c 26056.32c and
13027.32c:

    python3 examples/python/sandwich_slab_transport.py /path/to/neutron/ace
"""

import argparse
import math
from pathlib import Path

import numpy as np
import pyalea


ENERGY = 14.1
FE_DENSITY = 0.08  # atoms/barn-cm
AL_DENSITY = 0.06  # atoms/barn-cm
LAYERS = ((1, 0.5, None), (2, 1.0, "fe"), (3, 1.0, None),
          (4, 2.0, "al"), (5, 1.0, None), (6, 1.0, "fe"),
          (7, 1.0, None))


def geometry(materials: bool) -> str:
    fe = f"1 {FE_DENSITY:.17g}" if materials else "0"
    al = f"2 {AL_DENSITY:.17g}" if materials else "0"
    cards = "\nm1 26056.32c 1\nm2 13027.32c 1\n" if materials else ""
    return (
        "Sandwich slab neutron benchmark\n"
        "c\n"
        "1 0 -90 -1 imp:n=1\n"
        f"2 {fe} -90 1 -2 imp:n=1\n"
        "3 0 -90 2 -3 imp:n=1\n"
        f"4 {al} -90 3 -4 imp:n=1\n"
        "5 0 -90 4 -5 imp:n=1\n"
        f"6 {fe} -90 5 -6 imp:n=1\n"
        "7 0 -90 6 imp:n=1\n"
        "99 0 90 imp:n=0\n"
        "\n"
        "1 px 0\n"
        "2 px 1\n"
        "3 px 2\n"
        "4 px 4\n"
        "5 px 5\n"
        "6 px 6\n"
        "90 rpp -1 7 -5 5 -5 5\n"
        f"{cards}"
    )


def tally_cells(result: dict, index: int) -> dict[int, tuple[float, float]]:
    tally = result["tallies"][index]
    return {
        int(cell): (float(mean), float(error))
        for cell, mean, error in zip(
            tally["bin_ids"], tally["mean"], tally["standard_error"]
        )
    }


def analytic_uncollided(sigmas: dict[str, float]) -> tuple[dict[int, float], float]:
    """Integrate survival exp(-Sigma*x) over each slab's beam path."""
    survival = 1.0
    expected = {}
    for cell, width, material in LAYERS:
        sigma = sigmas.get(material, 0.0)
        integrated_path = -math.expm1(-sigma * width) / sigma if sigma else width
        expected[cell] = survival * integrated_path
        survival *= math.exp(-sigma * width)
    return expected, survival


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("xsdir", type=Path, help="ACE directory or xsdir file")
    parser.add_argument("--histories", type=int, default=20000)
    parser.add_argument("--seed", type=int, default=17)
    args = parser.parse_args()
    if args.histories < 1000:
        parser.error("--histories must be at least 1000 for the statistical check")

    xsdir = pyalea.XsDir(str(args.xsdir), directory=args.xsdir.is_dir())
    sigmas = {
        "fe": FE_DENSITY * pyalea.Nuclide(xsdir, "26056.32c").xs_total(ENERGY),
        "al": AL_DENSITY * pyalea.Nuclide(xsdir, "13027.32c").xs_total(ENERGY),
    }
    if not all(math.isfinite(sigma) and sigma > 0 for sigma in sigmas.values()):
        raise RuntimeError("invalid macroscopic total cross section")
    source = {
        "particle": "neutron",
        "space": {"type": "point", "position": [-0.5, 0, 0]},
        "angle": {"type": "monodirectional", "direction": [1, 0, 0]},
        "energy": ENERGY,
    }
    full = {"domain": "cell", "score": "track_length", "particle": "neutron"}
    uncollided = {**full, "energy_min": ENERGY,
                  "energy_max": float(np.nextafter(ENERGY, math.inf))}
    config = {"histories": args.histories, "seed": args.seed, "source": source}

    # With every material removed, exact cell paths are the slab thicknesses.
    void = pyalea.transport_run(pyalea.load_mcnp_string(geometry(False)), xsdir,
                                {**config, "tallies": [full]})
    void_cells = tally_cells(void, 0)
    if void["collisions"] != 0 or void["leaked"] != args.histories:
        raise RuntimeError("void slab did not leak every source neutron")
    for cell, width, _ in LAYERS:
        np.testing.assert_allclose(void_cells[cell][0], width, rtol=0, atol=1e-10)
    print("Void geometry: all seven slab track lengths match their exact widths.")

    result = pyalea.transport_run(pyalea.load_mcnp_string(geometry(True)), xsdir,
                                   {**config, "tallies": [uncollided, full]})
    expected, transmission = analytic_uncollided(sigmas)
    first_flight = tally_cells(result, 0)
    all_flights = tally_cells(result, 1)
    print(f"Macroscopic total XS: Fe-56 {sigmas['fe']:.8f} cm^-1, "
          f"Al-27 {sigmas['al']:.8f} cm^-1")
    print(f"Analytic uncollided transmission after the sandwich: {transmission:.6f}")
    print("cell  material  analytic    sampled +/- stderr  all-energy track length")
    for cell, _, material in LAYERS:
        mean, error = first_flight[cell]
        if not math.isfinite(mean) or not math.isfinite(error):
            raise RuntimeError(f"non-finite tally in cell {cell}")
        if abs(mean - expected[cell]) > 5 * error + 1e-3:
            raise RuntimeError(f"cell {cell} differs from analytic attenuation by >5 stderr")
        print(f"{cell:4d}  {(material or 'void'):>8}  {expected[cell]:8.5f}  "
              f"{mean:8.5f} +/- {error:.5f}     {all_flights[cell][0]:8.5f}")
    print(f"Uncollided transmission in the downstream 1 cm void detector: "
          f"{first_flight[7][0]:.6f} +/- {first_flight[7][1]:.6f}")
    print(f"PASS: {args.histories} histories, {result['collisions']} collisions, "
          f"{result['leaked']} leaked particles")


if __name__ == "__main__":
    main()
