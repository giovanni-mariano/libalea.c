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
    metal = os.environ.get("USE_METAL") == "1"
    with tempfile.TemporaryDirectory(prefix="alea-sdl-") as directory:
        tmp = Path(directory)
        shell_test = tmp / "shell-test"
        run(compiler + ["-std=c11", "-O1", "-Wall", "-Wextra",
            str(root / "tests/tools/test_sdl_view.c"), "-o", str(shell_test)] + flags)
        for args in [[], ["close-during-render"], ["adaptive"], ["progressive"], ["progressive-close"]]:
            result = run([str(shell_test), *args], env=dict(os.environ, SDL_VIDEODRIVER="dummy"))
            print(result.stdout.strip())
        model = tmp / "model.i"
        model.write_text("Sphere\n1 1 -1 -1\n2 0 1\n\n1 so 5\n\n")
        inspector_obj = str(tmp / "inspector.o")
        run(compiler + ["-std=c11", "-O1", f"-I{root / 'include'}", *flags,
            "-c", str(root / "tools/render_inspector.c"), "-o", inspector_obj])
        build_test = tmp / "render-build"
        run(compiler + ["-std=c11", "-O1", "-Wall", "-Wextra", "-DALEA_USE_SDL",
            f"-I{root / 'include'}", str(root / "tests/tools/test_render_view_build.c"),
            inspector_obj, "-o", str(build_test)] + libraries + flags)
        for scenario in ("ready", "fail", "renderfail", "close", "cpu", "xray"):
            output = tmp / f"build-{scenario}.bmp"
            options = {"cpu": ["--backend", "cpu"], "xray": ["--mode", "xray"]}.get(scenario, [])
            run([str(build_test), str(model), "--interactive", "--width", "64",
                 "--height", "64", "-o", str(output), *options],
                env=dict(os.environ, SDL_VIDEODRIVER="dummy", ALEA_TEST_BUILD=scenario))
            if scenario != "close":
                assert output.read_bytes()[:2] == b"BM"
                if scenario in ("fail", "renderfail"):
                    assert output.read_bytes() == (tmp / "build-ready.bmp").read_bytes()
        run([str(build_test), str(model), "--interactive", "--backend", "regional",
             "--width", "65", "--height", "65", "--eye", "20,0,0",
             "--target", "0,0,0", "--clip-nx=0"],
            env=dict(os.environ, SDL_VIDEODRIVER="dummy", ALEA_TEST_BUILD="pick"))
        inspector_test = tmp / "inspector-test"
        run(compiler + ["-std=c11", "-O1", "-DALEA_USE_SDL", f"-I{root / 'include'}",
            str(root / "tests/tools/test_render_inspector.c"), "-o", str(inspector_test)]
            + libraries + flags)
        inspector_model = tmp / "inspector.i"
        inspector_model.write_text("Materials\n1 1 -1 -1\n2 2 -1 -2\n3 0 1 2\n\n"
            "1 so 5\n2 s 0 7 1 1.6\n\nm1 1001 1\nm2 8016 1\n")
        result = run([str(inspector_test), str(inspector_model), "--interactive", "--backend", "regional",
            "--width", "65", "--height", "65", "--eye", "20,0,0", "--target", "0,0,0",
            "--clip-nx=0"], env=dict(os.environ, SDL_VIDEODRIVER="dummy"))
        print(result.stdout.strip())
        print("Full-resolution cut-face pick from preview and off-center orbit: PASS")
        print("Regional viewer publication, fallback, save, quit, CPU override and X-ray: PASS")
        nested_model = tmp / "nested.i"
        nested_model.write_text(
            "Nested elliptical torus\n1 0 -1 fill=1\n2 1 -1 -2 u=1\n"
            "3 0 2 u=1\n4 0 1\n\n1 so 6\n2 tz 0 0 0 3 1 1.5\n\n")
        if metal:
            cancellation_test = tmp / "metal-slice-cancel"
            run(compiler + ["-std=c11", "-O1", "-Wall", "-Wextra", "-DALEA_VIEW_METAL",
                f"-I{root / 'include'}", str(root / "tests/tools/test_metal_slice.c"),
                "-o", str(cancellation_test), str(root / "bin/compute/metal/libalea_compute.a"),
                "-framework", "Foundation", "-framework", "Metal"] + libraries)
            print(run([str(cancellation_test), str(nested_model)]).stdout.strip())
        viewers = [
            ("slice", "tools/mc_plotter.c", ["Z", "0", "-10", "10", "-10", "10", "64"]),
            ("render", "examples/c/render3d.c", ["--width", "64", "--height", "64"]),
            ("regional", "examples/c/render3d.c", ["--width", "64", "--height", "64", "--backend", "regional"]),
        ]
        if metal:
            viewers.append(("metal", "examples/c/render3d.c",
                ["--width", "64", "--height", "64", "--backend", "metal"]))
            viewers.append(("metal_nested", "examples/c/render3d.c",
                ["--width", "64", "--height", "64", "--backend", "metal"]))
            for name in ("metal_slice", "metal_slice_nested"):
                viewers.append((name, "tools/mc_plotter.c",
                    ["Z", "0", "-10", "10", "-10", "10", "96", "--backend", "metal"]))
            viewers.append(("metal_slice_selected", "tools/mc_plotter.c",
                ["Z", "0", "-10", "10", "-10", "10", "96", "--backend=metal", "--ownership=selected"]))
            viewers.append(("metal_slice_lattice", "tools/mc_plotter.c",
                ["Z", "0", "-10", "10", "-10", "10", "96", "--backend=metal"]))
        for name, source, options in viewers:
            is_slice = source == "tools/mc_plotter.c"
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
        event.key.keysym.sym = *keys == '>' ? SDLK_PAGEUP : *keys;
        if (SDL_PushEvent(&event) != 1) return 3;
    }}
    return viewer_main(argc, argv);
}}
''')
            exe = tmp / name
            extra = (["-DALEA_VIEW_METAL", str(root / "bin/compute/metal/libalea_compute.a"),
                "-framework", "Foundation", "-framework", "Metal"] if name.startswith("metal") else [])
            run(compiler + ["-std=c11", "-O1", "-Wall", "-Wextra", "-DALEA_USE_SDL",
                f"-I{root / 'include'}", str(wrapper), "-o", str(exe)] + ([] if is_slice else [inspector_obj]) + extra + libraries + flags)
            images = {}
            actions = [("initial", "sq"), ("zoom", "+sq"), ("reset", "+rsq")]
            if name.startswith("metal_slice"):
                actions += [("cpu_keys", "mesq"), ("slice_move", ">sq")]
                if name.endswith("nested"):
                    actions += [("axis", "xsq")]
            if not is_slice:
                actions += [("clip", "xsq"), ("clip_move", "x]sq"), ("clip_reset", "x]rsq")]
            for mode, keys in actions:
                output = tmp / f"{name}-{mode}.bmp"
                source_model = (root / "tests/data/mcnp_nested_lattice_transformed.mcnp"
                    if name.endswith("lattice") else nested_model if name.endswith("nested") else model)
                args = [str(exe), str(source_model), *options]
                if not is_slice: args += ["-o"]
                args += [str(output), "--interactive"]
                env = dict(os.environ, SDL_VIDEODRIVER="dummy", ALEA_TEST_KEYS=keys)
                run(args, env=env)
                images[mode] = output.read_bytes()
                assert images[mode][:2] == b"BM"
            assert images["initial"] != images["zoom"], f"{name}: zoom did not change image"
            assert images["initial"] == images["reset"], f"{name}: reset did not restore image"
            if not is_slice:
                assert images["clip"] != images["initial"], f"{name}: clip did not change image"
                assert images["clip"] != images["clip_move"], f"{name}: clip did not move"
                assert images["initial"] == images["clip_reset"], f"{name}: clip reset failed"
            if name == "render":
                cpu_images = images
            if name == "regional":
                assert images == cpu_images, "Regional frames differ from CPU frames"
            if name.startswith("metal_slice"):
                assert images["initial"] == images["cpu_keys"], "CPU-only keys changed Metal view"
                assert images["initial"] != images["slice_move"], "Slice did not move"
                if name.endswith("nested"):
                    assert images["initial"] != images["axis"], "Slice axis did not change"
            print(f"{name}: save, zoom, reset, clip (3D), quit PASS")


if __name__ == "__main__":
    main()
