# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

import numpy as np
import pytest

import pyalea


GEOMETRY = """Void transport smoke
c
1 0 -1 imp:n=1 imp:p=1
99 0 1 imp:n=0 imp:p=0

1 s 0 0 0 2
"""


@pytest.mark.parametrize("particle", ["neutron", "photon"])
def test_fixed_source_and_tallies(tmp_path, particle):
    xsdir_path = tmp_path / "xsdir"
    xsdir_path.write_text(
        "directory\n1001.80c 1.0 dummy.ace 0 1 1 1 0 0 2.5301e-8\n"
    )
    system = pyalea.load_mcnp_string(GEOMETRY)
    xsdir = pyalea.XsDir(str(xsdir_path))
    config = {
        "histories": 4,
        "seed": 7,
        "boundary_distance_tolerance": 0.0,
        "source": {
            "particle": particle,
            "energy": 2.0,
            "space": {"type": "point", "position": [0, 0, 0]},
            "angle": {"type": "monodirectional", "direction": [1, 0, 0]},
        },
        "tallies": [
            {"score": "track_length", "particle": particle},
            {"score": "track_length", "domain": "mesh", "particle": particle,
             "lower": [-2, -2, -2], "upper": [2, 2, 2],
             "dimensions": [2, 1, 1], "energy_edges": [0, 3]},
        ],
    }
    result = pyalea.transport_run(system, xsdir, config)
    assert result["histories"] == 4
    assert result["leaked"] == 4
    assert isinstance(result["track_length"], np.ndarray)
    assert result["track_length"].dtype == np.float64
    assert isinstance(result["track_length_squared"], np.ndarray)
    assert result["tallies"][0]["bin_ids"].dtype == np.dtype(np.intc)
    assert result["tallies"][0]["bin_ids"].tolist() == [1, 99]
    for tally in result["tallies"]:
        for key in ("sum", "sum_squared", "mean", "standard_error"):
            assert isinstance(tally[key], np.ndarray)
            assert tally[key].dtype == np.float64
            assert tally[key].flags.owndata
    assert result["tallies"][1]["bin_ids"] is None
    assert result["tallies"][1]["energy_edges"].tolist() == [0, 3]
    assert result["tallies"][0]["mean"][0] == pytest.approx(2.0)
    assert result["tallies"][0]["standard_error"][0] == 0
    np.testing.assert_allclose(result["tallies"][1]["mean"], [0.0, 2.0])
    replay = pyalea.transport_run(system, xsdir, config)
    for first, second in zip(result["tallies"], replay["tallies"]):
        np.testing.assert_array_equal(first["sum"], second["sum"])
        np.testing.assert_array_equal(first["mean"], second["mean"])


def test_box_source_and_invalid_tally(tmp_path):
    xsdir_path = tmp_path / "xsdir"
    xsdir_path.write_text(
        "directory\n1001.80c 1.0 dummy.ace 0 1 1 1 0 0 2.5301e-8\n"
    )
    system = pyalea.load_mcnp_string(GEOMETRY)
    xsdir = pyalea.XsDir(str(xsdir_path))
    config = {
        "histories": 8,
        "source": {"particle": "neutron", "energy": 14.1,
                   "space": {"type": "box", "lower": [0, 0, 0],
                             "upper": [0, 0, 0]},
                   "angle": {"type": "isotropic"}},
        "tallies": [{"score": "track_length", "domain": "cell"}],
    }
    assert pyalea.transport_run(system, xsdir, config)["leaked"] == 8
    config["tallies"] = [{"score": "invalid"}]
    with pytest.raises(ValueError, match="unknown score"):
        pyalea.transport_run(system, xsdir, config)


def test_prepared_source_preview_and_reuse_without_nuclear_data(tmp_path):
    spec = {
        "particle": "photon",
        "space": {"type": "point", "position": [0, 0, 0]},
        "angle": {"type": "isotropic"},
        "energy": {"type": "mono", "value": 2.0},
        "time": {"type": "constant", "value": 0.25},
    }
    source = pyalea.Source(spec)
    preview = pyalea.sample_source(source, histories=5, seed=77, history_offset=4)
    assert preview["position"].shape == (5, 3)
    assert preview["direction"].shape == (5, 3)
    assert preview["position"].dtype == np.float64
    np.testing.assert_array_equal(preview["history_id"], [4, 5, 6, 7, 8])
    np.testing.assert_allclose(np.linalg.norm(preview["direction"], axis=1), 1)
    np.testing.assert_array_equal(preview["energy"], [2] * 5)
    np.testing.assert_array_equal(preview["particle"], [1] * 5)
    np.testing.assert_array_equal(preview["time"], [0.25] * 5)
    replay = pyalea.sample_source(source, 1, 77, 6)
    np.testing.assert_array_equal(replay["direction"][0], preview["direction"][2])
    with pytest.raises(ValueError, match="inside space and angle"):
        pyalea.Source({**spec, "kind": "fixed"})

    xsdir_path = tmp_path / "xsdir"
    xsdir_path.write_text(
        "directory\n1001.80c 1.0 dummy.ace 0 1 1 1 0 0 2.5301e-8\n"
    )
    system = pyalea.load_mcnp_string(GEOMETRY)
    result = pyalea.transport_run(system, pyalea.XsDir(str(xsdir_path)), {
        "histories": 5, "seed": 77, "history_offset": 4, "source": source,
        "tallies": [{"score": "track_length", "particle": "photon"}],
    })
    assert result["leaked"] == 5
    assert result["tallies"][0]["mean"][0] == pytest.approx(2)


@pytest.mark.parametrize("space", [
    {"type": "line", "start": [0, 0, 0], "end": [2, 0, 0]},
    {"type": "sphere", "center": [0, 0, 0],
     "inner_radius": 0.5, "outer_radius": 1},
    {"type": "cylinder", "base": [0, 0, 0], "axis": [0, 0, 2],
     "inner_radius": 0.5, "outer_radius": 1},
])
def test_analytic_spatial_source_preview(space):
    source = pyalea.Source({
        "space": space, "angle": {"type": "isotropic"}, "energy": 14.1,
    })
    samples = pyalea.sample_source(source, 100, 13)
    points = samples["position"]
    assert points.shape == (100, 3)
    if space["type"] == "line":
        assert np.all((points[:, 0] >= 0) & (points[:, 0] <= 2))
        assert np.all(points[:, 1:] == 0)
    elif space["type"] == "sphere":
        radius = np.linalg.norm(points, axis=1)
        assert np.all((radius >= 0.5) & (radius <= 1))
    else:
        radius = np.linalg.norm(points[:, :2], axis=1)
        assert np.all((radius >= 0.5) & (radius <= 1))
        assert np.all((points[:, 2] >= 0) & (points[:, 2] <= 2))


def test_discrete_energy_lines_preview():
    source = pyalea.Source({
        "space": {"type": "point", "position": [0, 0, 0]},
        "angle": {"type": "isotropic"},
        "energy": {"type": "lines", "values": [2.45, 14.1],
                   "weights": [1, 3]},
    })
    samples = pyalea.sample_source(source, 1000, 37)
    assert set(np.unique(samples["energy"])) == {2.45, 14.1}
    assert 0.7 < np.mean(samples["energy"] == 14.1) < 0.8
    replay = pyalea.sample_source(source, 1, 37, 17)
    assert replay["energy"][0] == samples["energy"][17]


@pytest.mark.parametrize("angle, expected", [
    ({"type": "cone", "direction": [0, 0, 1],
      "half_angle": np.pi / 3}, 0.75),
    ({"type": "cosine", "direction": [0, 0, 1]}, 2 / 3),
    ({"type": "tabulated_mu", "direction": [0, 0, 1],
      "mu": [0, 1], "pdf": [0, 2], "interpolation": "linear"}, 2 / 3),
])
def test_angular_distributions(angle, expected):
    source = pyalea.Source({
        "space": {"type": "point", "position": [1, 0, 0]},
        "angle": angle, "energy": 14.1,
    })
    samples = pyalea.sample_source(source, 10000, 19)
    assert np.mean(samples["direction"][:, 2]) == pytest.approx(expected, abs=0.02)


def test_radial_direction_uses_sampled_position():
    source = pyalea.Source({
        "space": {"type": "line", "start": [1, 0, 0], "end": [2, 0, 0]},
        "angle": {"type": "radial", "origin": [0, 0, 0], "inward": True},
        "energy": 14.1,
    })
    samples = pyalea.sample_source(source, 16, 19)
    np.testing.assert_allclose(samples["direction"], [[-1, 0, 0]] * 16)


def test_tokamak_rz_emissivity_preview(tmp_path):
    r_edges = np.array([0.0, 1.0, 2.0])
    emissivity = np.array([[1.0, 0.0], [1.0, 1.0]])
    source = pyalea.Source({
        "space": {"type": "tokamak_rz", "r_edges": r_edges,
                  "z_edges": [0, 1, 2], "emissivity": emissivity,
                  "phi_min": 0, "phi_max": np.pi / 2},
        "angle": {"type": "isotropic"}, "energy": 14.1,
    })
    assert source.integrated_emissivity == pytest.approx(1.75 * np.pi)
    r_edges[1] = 100
    emissivity[:] = 0
    samples = pyalea.sample_source(source, 20000, 73)
    points = samples["position"]
    r2 = points[:, 0]**2 + points[:, 1]**2
    assert np.all(points[:, 0] >= -1e-14)
    assert np.all(points[:, 1] >= -1e-14)
    assert not np.any((r2 < 1) & (points[:, 2] >= 1))
    assert np.mean(r2 >= 1) == pytest.approx(6 / 7, abs=0.02)
    np.testing.assert_array_equal(
        pyalea.sample_source(source, 1, 73, 123)["position"][0],
        points[123],
    )
    xsdir_path = tmp_path / "xsdir"
    xsdir_path.write_text(
        "directory\n1001.80c 1.0 dummy.ace 0 1 1 1 0 0 2.5301e-8\n"
    )
    inside = pyalea.Source({
        "space": {"type": "tokamak_rz", "r_edges": [0, 0.5, 1],
                  "z_edges": [0, 1], "emissivity": [[1], [1]]},
        "angle": {"type": "isotropic"}, "energy": 14.1,
    })
    result = pyalea.transport_run(
        pyalea.load_mcnp_string(GEOMETRY), pyalea.XsDir(str(xsdir_path)),
        {"histories": 8, "source": inside,
         "tallies": [{"score": "track_length"}]},
    )
    assert result["leaked"] == 8
    with pytest.raises(ValueError, match="one value per Z bin"):
        pyalea.Source({
            "space": {"type": "tokamak_rz", "r_edges": [0, 1],
                      "z_edges": [0, 1, 2], "emissivity": [[1]]},
            "angle": {"type": "isotropic"}, "energy": 14.1,
        })


def test_tabulated_energy_preview():
    source = pyalea.Source({
        "space": {"type": "point", "position": [0, 0, 0]},
        "angle": {"type": "isotropic"},
        "energy": {"type": "tabulated", "values": [1, 3],
                   "pdf": [1, 1], "interpolation": "linear"},
    })
    energy = pyalea.sample_source(source, 10000, 42)["energy"]
    assert np.all((energy >= 1) & (energy <= 3))
    assert np.mean(energy) == pytest.approx(2, abs=0.03)


def test_weighted_source_mixture_preview():
    source = pyalea.Source({
        "type": "mixture",
        "components": [
            {"strength": 1, "source": {
                "particle": "neutron", "space": {"type": "point", "position": [0, 0, 0]},
                "angle": {"type": "isotropic"}, "energy": 2,
            }},
            {"strength": 3, "source": {
                "particle": "photon", "space": {"type": "point", "position": [1, 0, 0]},
                "angle": {"type": "isotropic"}, "energy": 3,
            }},
        ],
    })
    samples = pyalea.sample_source(source, 10000, 17)
    photons = samples["particle"] == 1
    assert np.mean(photons) == pytest.approx(.75, abs=.02)
    assert np.all(samples["position"][photons, 0] == 1)
    assert np.all(samples["weight"] == 1)


def test_cartesian_source_mesh_density_and_strength(tmp_path):
    x_edges = np.array([0.0, 1.0, 3.0])
    values = np.array([[[1.0, 0.0]], [[1.0, 1.0]]])
    space = {"type": "cartesian_mesh", "x_edges": x_edges,
             "y_edges": [0, 1], "z_edges": [0, 1, 2],
             "values": values, "value_mode": "density"}
    source = pyalea.Source({"space": space, "angle": {"type": "isotropic"},
                           "energy": 14.1})
    assert source.integrated_emissivity == pytest.approx(5)
    x_edges[1] = 100
    values[:] = 0
    points = pyalea.sample_source(source, 20000, 63)["position"]
    assert not np.any((points[:, 0] < 1) & (points[:, 2] >= 1))
    assert np.mean(points[:, 0] >= 1) == pytest.approx(.8, abs=.02)
    strength = pyalea.Source({
        "space": {**space, "x_edges": [0, 1, 3],
                  "values": [[[1, 0]], [[1, 1]]], "value_mode": "strength"},
        "angle": {"type": "isotropic"}, "energy": 14.1,
    })
    assert strength.integrated_emissivity == pytest.approx(3)
    assert np.mean(pyalea.sample_source(strength, 10000, 64)["position"][:, 0] >= 1) == pytest.approx(2 / 3, abs=.02)
    with pytest.raises(ValueError, match="one value per Z bin"):
        pyalea.Source({
            "space": {**space, "x_edges": [0, 1, 3],
                      "values": [[[1]], [[1, 1]]]},
            "angle": {"type": "isotropic"}, "energy": 14.1,
        })
    xsdir_path = tmp_path / "xsdir"
    xsdir_path.write_text(
        "directory\n1001.80c 1.0 dummy.ace 0 1 1 1 0 0 2.5301e-8\n"
    )
    inside = pyalea.Source({
        "space": {"type": "cartesian_mesh", "x_edges": [0, 1],
                  "y_edges": [0, 1], "z_edges": [0, 1],
                  "values": [[[1]]], "value_mode": "density"},
        "angle": {"type": "isotropic"}, "energy": 14.1,
    })
    result = pyalea.transport_run(
        pyalea.load_mcnp_string(GEOMETRY), pyalea.XsDir(str(xsdir_path)),
        {"histories": 8, "source": inside,
         "tallies": [{"score": "track_length"}]},
    )
    assert result["leaked"] == 8
