#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

"""Exercise 50 keV U-238 URR transport across an artificial cell split.

Supply an xsdir containing Lib80x 92238.00c:

    python3 examples/python/urr_u238_transport.py /path/to/Lib80x/xsdir

The two geometries have the same material, density, and outer sphere. Splitting
the sphere into adjacent cells must not alter the sampled resonance band,
collision history, or total material score for a fixed seed.
"""

import argparse
from pathlib import Path

import numpy as np
import pyalea


def geometry(split: bool) -> str:
    if split:
        cells = "1 1 0.05 -2 imp:n=1\n2 1 0.05 2 -1 imp:n=1\n"
        surfaces = "1 so 3\n2 so 1\n"
    else:
        cells = "1 1 0.05 -1 imp:n=1\n"
        surfaces = "1 so 3\n"
    return ("U-238 unresolved-resonance transport\nc\n" + cells +
            "99 0 1 imp:n=0\n\n" + surfaces +
            "\nm1 92238.00c 1\n")


def material_score(result: dict, tally_index: int) -> float:
    tally = result["tallies"][tally_index]
    return float(sum(mean for cell, mean in zip(tally["bin_ids"], tally["mean"])
                     if cell != 99))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("xsdir", type=Path)
    parser.add_argument("--histories", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()
    if args.histories < 1:
        parser.error("--histories must be positive")

    xsdir = pyalea.XsDir(str(args.xsdir), directory=args.xsdir.is_dir())
    config = {
        "histories": args.histories, "seed": args.seed,
        "source": {
            "particle": "neutron", "energy": 0.05,
            "space": {"type": "point", "position": [0, 0, 0]},
            "angle": {"type": "isotropic"},
        },
        "tallies": [
            {"domain": "cell", "score": "track_length", "particle": "neutron"},
            {"domain": "cell", "score": "reaction_rate", "particle": "neutron",
             "reaction_mt": 102, "nuclide_zaid": 92238},
        ],
    }
    whole = pyalea.transport_run(pyalea.load_mcnp_string(geometry(False)),
                                 xsdir, config)
    split = pyalea.transport_run(pyalea.load_mcnp_string(geometry(True)),
                                 xsdir, config)
    for key in ("histories", "collisions", "absorbed", "leaked",
                "emitted_neutrons"):
        if whole[key] != split[key]:
            raise RuntimeError(f"cell subdivision changed {key}: "
                               f"{whole[key]} versus {split[key]}")
    for index, name in enumerate(("track length", "MT=102 reaction rate")):
        a, b = material_score(whole, index), material_score(split, index)
        np.testing.assert_allclose(a, b, rtol=1e-11, atol=1e-11)
        print(f"{name}: {a:.8f} per source neutron in both geometries")
    print(f"PASS: {args.histories} histories, {whole['collisions']} collisions, "
          f"{whole['leaked']} leaked particles")


if __name__ == "__main__":
    main()
