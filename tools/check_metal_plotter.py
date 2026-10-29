#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0
"""Check native Metal slice images and explicit failures using public fixtures.

Run after building bin/mc_plotter with USE_METAL=1. Requires a Metal device.
"""
import os
from pathlib import Path
import struct
import subprocess
import tempfile


def bitmap(data):
    offset = struct.unpack_from('<I', data, 10)[0]
    width, height = struct.unpack_from('<ii', data, 18)
    assert width == height == 96
    stride = (width * 3 + 3) & ~3
    # BMP rows run from bottom to top; convert stored BGR to RGB.
    return [[tuple(reversed(data[offset+y*stride+x*3:offset+y*stride+x*3+3]))
             for x in range(width)] for y in range(height)]


def main():
    root = Path(__file__).resolve().parents[1]
    env = {k: v for k, v in os.environ.items() if not k.startswith(('ALEA_PLOT_', 'ALEA_GRID_'))}
    env['ALEA_NUM_THREADS'] = '2'
    with tempfile.TemporaryDirectory(prefix='alea-metal-plot-') as directory:
        tmp = Path(directory)
        model = tmp / 'nested.i'
        # A displaced sphere makes wrong axis mapping or vertical flips visible.
        model.write_text('Translated fill\n1 0 -1 fill=1 (1 2 3)\n2 0 1\n'
                         '10 1 -1 -2 u=1\n11 0 2 u=1\n\n1 so 9\n2 so 2\n\n')

        def run(name, options=(), source=model, axis='Z', value='3', environment=env, error=None,
                bounds=('-10', '10', '-10', '10')):
            output = tmp / f'{name}.bmp'
            command = [str(root / 'bin/mc_plotter'), str(source), axis, value,
                       *bounds, '96', *options, str(output)]
            result = subprocess.run(command, env=environment, capture_output=True, text=True, timeout=60)
            if error:
                assert result.returncode != 0, command
                assert error.lower() in (result.stdout + result.stderr).lower(), result.stderr
                assert not output.exists(), 'Failed request published an image'
                return None
            assert result.returncode == 0, result.stdout + result.stderr
            return output.read_bytes()

        def metal_image(name, options, **kwargs):
            scan = run(name, options, environment=dict(env, ALEA_METAL_POINT_CANDIDATES='0'), **kwargs)
            tiled = run(name + '-tiles', options,
                        environment=dict(env, ALEA_METAL_POINT_CANDIDATES='1'), **kwargs)
            assert scan == tiled, ('Candidate tiles changed pixels', name)
            return tiled

        for axis, value in [('X', '1'), ('Y', '2'), ('Z', '3')]:
            for style in [[], ['--color=materials', '--contours=materials'],
                          ['--labels=cells,materials', '--ticks']]:
                cpu = run('cpu', ['--backend=cpu', *style], axis=axis, value=value)
                gpu = metal_image('gpu', ['--backend', 'metal', *style], axis=axis, value=value)
                assert cpu == gpu, (axis, style)
                assert cpu == metal_image('gpu-selected', ['--backend=metal', '--ownership=selected', *style], axis=axis, value=value)
        rows = bitmap(metal_image('orientation', ['--backend=metal']))
        assert rows[57][52] != (255, 255, 255)  # Near world (1,2).
        assert rows[38][52] == (255, 255, 255)  # Near world (1,-2).

        ownership = tmp / 'ownership.i'
        ownership.write_text('Ownership\n1 1 -1 -1\n2 2 -1 -2\n3 0 -3\n\n'
                             '1 s -2 0 0 2\n2 s 0 0 0 2\n3 s 5 0 0 1\n\n')
        rows = bitmap(metal_image('ownership', ['--backend=metal'], source=ownership, value='0'))
        assert rows[48][43] == (255, 0, 255), 'Overlap must be magenta'
        assert rows[48][72] == (255, 255, 255), 'Unique void must be white'
        assert rows[10][10] == (255, 192, 128), 'Gap must be orange'
        assert rows[48][33] not in [(255, 0, 255), (255, 192, 128), (255, 255, 255)]

        # Selected mode matches CPU precedence on overlaps; all-owner mode
        # retains its overlap color. Candidate tiles must preserve both views.
        selected = metal_image('selected', ['--backend=metal', '--ownership=selected'],
                               source=ownership, value='0')
        colors = bitmap(selected)
        assert colors[48][43] != (255, 0, 255), 'Selected overlap must choose an owner'
        cpu_selected = bitmap(run('cpu-selected', ['--backend=cpu'], source=ownership, value='0'))
        assert colors[48][43] == cpu_selected[48][43]
        assert colors[10][10] == (255, 192, 128)
        hole = tmp / 'hole.i'
        hole.write_text('Undefined fill\n1 0 -1 fill=1\n2 2 -1 -1\n10 1 -1 -2 u=1\n\n1 so 9\n2 so 2\n\n')
        selected_hole = bitmap(metal_image('undefined', ['--backend=metal', '--ownership', 'selected'],
                                          source=hole, value='0'))
        assert selected_hole[48][65] == (255, 255, 0), 'Undefined fill must be yellow'
        all_hole = bitmap(metal_image('all-hole', ['--backend=metal', '--ownership=all'], source=hole, value='0'))
        assert all_hole[48][65] != (255, 255, 0)
        for option in ['--ownership=wrong', '--ownership']:
            run('bad-ownership', ['--backend=metal', option], error='ownership')
        run('cpu-all', ['--backend=cpu', '--ownership=all'], error='requires')

        for option in ['--errors', '--method=ray', '--labels=surfaces', '--labels=all',
                       '--debug', '--trace=4,4']:
            run('reject-option', ['--backend=metal', option], error='do not support')
        for key in ['ALEA_PLOT_EXACT_COVERAGE', 'ALEA_PLOT_FULL_EXACT_COVERAGE',
                    'ALEA_PLOT_TILE_COVERAGE', 'ALEA_PLOT_PATH_COVERAGE',
                    'ALEA_PLOT_ERROR_STATS', 'ALEA_PLOT_ERROR_VERIFY']:
            run('reject-env', ['--backend=metal'], environment=dict(env, **{key: '1'}),
                error='do not support')
        for fixture in ['mcnp_lattice_repeating.mcnp', 'mcnp_nested_lattice_transformed.mcnp']:
            source = root / 'tests/data' / fixture
            bounds = ('4.125', '7.875', '-.875', '.875') if 'nested' in fixture else ('-10', '10', '-10', '10')
            cpu = run('lattice-cpu', ['--backend=cpu'], source=source, value='0', bounds=bounds)
            gpu = metal_image('lattice-gpu', ['--backend=metal'], source=source, value='0', bounds=bounds)
            assert cpu == gpu, fixture
        run('reject-lattice', ['--backend=metal'],
            source=root / 'tests/data/mcnp_hex_lattice.mcnp', error='rectangular')
        run('reject-backend', ['--backend=cuda'], error='unknown backend')
        run('reject-bounds', ['--backend=metal'], value='nan', error='finite')
        run('reject-batch', ['--backend=metal', '--batch=missing.txt'], error='combined')
        print('Metal slices: axes/orientation, transformed fill, candidate tiles/edges, styles, ownership colors, explicit failures PASS')


if __name__ == '__main__':
    main()
