# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0
import json
import pyalea

FIRST = "identity\n1 1 -1 -10\n2 0 10\n\n10 so 5\n\nm1 1001 1\n"
SECOND = "identity\n9 7 -1 -20\n8 0 20\n\n20 so 5\n\nm7 1001 1\n"


def test_fingerprints_and_revision_reports():
    a = pyalea.load_mcnp_string(FIRST)
    b = pyalea.load_mcnp_string(SECOND)
    assert a.fingerprints()["cells"][1] == b.fingerprints()["cells"][9]
    assert a.fingerprints()["version"] == 2
    changes = a.compare(b)
    assert {c["change"] for c in changes if c["kind"] == "cell"} == {"renumbered"}
    assert json.loads(json.dumps(changes)) == changes
    changed = pyalea.load_mcnp_string(FIRST.replace("-1 -10", "-2 -10"))
    assert any(c["change"] == "content" for c in a.compare(changed))


def test_explicit_cache_roundtrip_and_disable(tmp_path, monkeypatch):
    source = tmp_path / "model.inp"
    cache = tmp_path / "prepared.bin"
    source.write_text(FIRST)
    a = pyalea.load_mcnp(str(source), cache_path=str(cache))
    assert cache.exists()
    b = pyalea.load_mcnp(str(source), cache_path=str(cache))
    assert a.export_mcnp_string() == b.export_mcnp_string()
    assert a.fingerprints() == b.fingerprints()
    cache.write_bytes(b"broken")
    monkeypatch.setenv("ALEA_CACHE", "off")
    c = pyalea.load_mcnp(str(source), cache_path=str(cache))
    assert cache.read_bytes() == b"broken"
    assert a.fingerprints() == c.fingerprints()


def test_directory_cache_budget(tmp_path, monkeypatch):
    directory = tmp_path / "cache"
    source = tmp_path / "model.inp"
    source.write_text(FIRST)
    monkeypatch.setenv("ALEA_CACHE_DIR", str(directory))
    pyalea.load_mcnp(str(source))
    assert len(list(directory.glob("*.abin"))) == 1
    source.write_text(FIRST.replace("so 5", "so 6"))
    pyalea.load_mcnp(str(source))
    assert len(list(directory.glob("*.abin"))) == 2
    monkeypatch.setenv("ALEA_CACHE_GB", "0")
    pyalea.load_mcnp(str(source))
    assert not list(directory.glob("*.abin"))


def test_semantic_material_profile():
    import pytest
    a = pyalea.load_mcnp_string(FIRST.replace("m1 1001 1", "m1 1001 1 8016 2"))
    b = pyalea.load_mcnp_string(FIRST.replace("m1 1001 1", "m1 1001 2 8016 4"))
    assert any(c["change"] == "content" for c in a.compare(b))
    assert all(c["change"] == "unchanged" for c in a.compare(b, profile="semantic"))
    assert a.fingerprints(profile="semantic")["profile"] == "semantic"
    assert a.fingerprints()["version"] == 2
    with pytest.raises(ValueError):
        a.fingerprints(profile="unknown")
