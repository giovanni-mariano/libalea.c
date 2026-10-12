# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0

"""Check Windows response files without requiring an MSVC toolchain."""

import importlib.util
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch


def load_config():
    path = Path(__file__).parents[1] / "scripts" / "config.py"
    spec = importlib.util.spec_from_file_location("pyalea_build_config", path)
    module = importlib.util.module_from_spec(spec)
    with patch.dict("sys.modules", {
        "numpy": SimpleNamespace(get_include=lambda: "numpy-include"),
    }):
        spec.loader.exec_module(module)
    return module


class MsvcResponseTests(unittest.TestCase):
    def test_python_import_library(self):
        config = load_config()
        cases = (
            ("python312.dll", "python312.lib"),
            ("python313.dll", "python313.lib"),
            ("python314.dll", "python314.lib"),
            ("python313_d.dll", "python313_d.lib"),
            ("python313t.dll", "python313t.lib"),
            ("python313.lib", "python313.lib"),
            (None, "python313.lib"),
        )
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory).resolve()
            for runtime, expected in cases:
                for configured_dir in (None, str(root / "custom libs")):
                    with self.subTest(runtime=runtime, libdir=configured_dir):
                        variables = {"LDLIBRARY": runtime, "LIBDIR": configured_dir}
                        with patch.object(config.sysconfig, "get_config_var", variables.get), \
                             patch.object(config.sysconfig, "get_path", return_value=str(root / "include")), \
                             patch.object(config.sys, "base_prefix", str(root / "Python")), \
                             patch.object(config.sys, "version_info", SimpleNamespace(major=3, minor=13)):
                            compile_path = root / "build" / "compile.rsp"
                            link_path = root / "build" / "link.rsp"
                            config.write_msvc_responses(compile_path, link_path, root)
                        library_dir = Path(configured_dir) if configured_dir else root / "Python" / "libs"
                        self.assertEqual(link_path.read_text(),
                                         f'/LIBPATH:"{library_dir.resolve()}"\n{expected}\n')
                        self.assertIn(f'/I"{root / "include"}"', compile_path.read_text())


if __name__ == "__main__":
    unittest.main()
