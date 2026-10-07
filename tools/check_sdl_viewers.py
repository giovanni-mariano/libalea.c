#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0
"""Exercise real viewer event loops with SDL's dummy driver (Unix dev build).

Run after make modules: python3 tools/check_sdl_viewers.py
Builds disposable test entry points; no test hooks are added to shipped tools.
"""
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile


def run(command, **kwargs):
    result = subprocess.run(command, capture_output=True, text=True, timeout=60, **kwargs)
    if result.returncode:
        raise RuntimeError(f"{command}\n{result.stdout}\n{result.stderr}")
    return result


def main():
    root = Path(__file__).resolve().parents[1]
    sdl_config = os.environ.get("SDL_CONFIG")
    if sdl_config:
        flags = shlex.split(run([sdl_config, "--cflags", "--static-libs"]).stdout)
    else:
        flags = shlex.split(run(["pkg-config", "--cflags", "--libs", "sdl2"]).stdout)
    compiler = shlex.split(os.environ.get("CC", "cc"))
    libraries = []
    for name in ("mcnp", "openmc", "xml", "serpent"):
        path = str(root / "bin" / f"libalea_{name}.a")
        if sys.platform == "darwin":
            libraries.append(f"-Wl,-force_load,{path}")
        else:
            libraries += ["-Wl,--whole-archive", path, "-Wl,--no-whole-archive"]
    libraries += [str(root / "bin/libalea.a"), "-lm", "-pthread"]
    with tempfile.TemporaryDirectory(prefix="alea-sdl-") as directory:
        tmp = Path(directory)
        shell_test = tmp / "shell-test"
        run(compiler + ["-std=c11", "-O1", "-Wall", "-Wextra",
            str(root / "tests/tools/test_sdl_view.c"), "-o", str(shell_test)] + flags)
        for args in [[], ["close-during-render"]]:
            result = run([str(shell_test), *args], env=dict(os.environ, SDL_VIDEODRIVER="dummy"))
            print(result.stdout.strip())
        model = tmp / "model.i"
        model.write_text("Sphere\n1 1 -1 -1\n2 0 1\n\n1 so 5\n\n")
        for name, source, options in [
            ("slice", "tools/mc_plotter.c", ["Z", "0", "-10", "10", "-10", "10", "64"]),
            ("render", "examples/c/render3d.c", ["--width", "64", "--height", "64"]),
        ]:
            wrapper = tmp / f"{name}.c"
            wrapper.write_text(f'''#define main viewer_main
#include "{root / source}"
#undef main
int main(int argc, char** argv) {{
    if (SDL_Init(SDL_INIT_VIDEO) != 0) return 2;
    const char* keys = getenv("ALEA_TEST_KEYS");
    for (; keys && *keys; ++keys) {{
        SDL_Event event = {{0}};
        event.type = SDL_KEYDOWN;
        event.key.keysym.sym = *keys;
        if (SDL_PushEvent(&event) != 1) return 3;
    }}
    return viewer_main(argc, argv);
}}
''')
            exe = tmp / name
            run(compiler + ["-std=c11", "-O1", "-Wall", "-Wextra", "-DALEA_USE_SDL",
                f"-I{root / 'include'}", str(wrapper), "-o", str(exe)] + libraries + flags)
            images = {}
            for mode, keys in [("initial", "sq"), ("zoom", "+sq"), ("reset", "+rsq")]:
                output = tmp / f"{name}-{mode}.bmp"
                args = [str(exe), str(model), *options]
                if name == "render": args += ["-o"]
                args += [str(output), "--interactive"]
                env = dict(os.environ, SDL_VIDEODRIVER="dummy", ALEA_TEST_KEYS=keys)
                run(args, env=env)
                images[mode] = output.read_bytes()
                assert images[mode][:2] == b"BM"
            assert images["initial"] != images["zoom"], f"{name}: zoom did not change image"
            assert images["initial"] == images["reset"], f"{name}: reset did not restore image"
            print(f"{name}: save, zoom, reset, quit PASS")


if __name__ == "__main__":
    main()
