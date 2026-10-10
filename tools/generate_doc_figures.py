#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0
"""Regenerate documentation figures: make tools && python3 tools/generate_doc_figures.py.

Uses the native render3d/mc_plotter binaries and Python's standard library.
PNG files are actual library output; SVG files are explanatory diagrams.
"""
from pathlib import Path
import subprocess
import tempfile
import struct
import zlib
from html import escape

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'docs/assets/figures'
INK = '#243447'
MUTED = '#526578'
ORANGE, GREEN, PURPLE = '#ff7f0e', '#2ca02c', '#9467bd'


def text(x, y, value, size=18, color=INK, anchor='start'):
    return f'<text x="{x}" y="{y}" font-size="{size}" fill="{color}" text-anchor="{anchor}">{escape(value)}</text>'


def rect(x, y, w, h, fill, radius=0):
    return f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{radius}" fill="{fill}"/>'


def circle(x, y, r, fill, stroke=INK, extra=''):
    return f'<circle cx="{x}" cy="{y}" r="{r}" fill="{fill}" stroke="{stroke}" stroke-width="2" {extra}/>'


def line(x1, y1, x2, y2, extra=''):
    return f'<path d="M{x1} {y1} L{x2} {y2}" fill="none" stroke="{INK}" stroke-width="2" {extra}/>'


def svg(name, title, description, height, body):
    (OUT / name).write_text(
        '<svg xmlns="http://www.w3.org/2000/svg" width="1080" '
        f'height="{height}" viewBox="0 0 1080 {height}" role="img" aria-labelledby="title desc">'
        '<!-- SPDX-FileCopyrightText: 2026 Giovanni MARIANO; SPDX-License-Identifier: MPL-2.0 -->'
        f'<title id="title">{escape(title)}</title><desc id="desc">{escape(description)}</desc>'
        '<defs><marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" '
        'markerWidth="7" markerHeight="7" orient="auto-start-reverse">'
        f'<path d="M0 0 L10 5 L0 10 Z" fill="{INK}"/></marker></defs>'
        f'<g font-family="Arial, Helvetica, sans-serif">{rect(0, 0, 1080, height, "#ffffff", 16)}'
        + body + '</g></svg>\n', encoding='utf-8')


def diagrams():
    body = text(36, 44, 'Surface sense', 26)
    body += text(36, 73, 'A surface is the boundary; its sense selects a side.', 18, MUTED)
    body += rect(28, 96, 500, 300, '#f3f6fa', 12) + rect(548, 96, 504, 300, '#f3f6fa', 12)
    body += text(54, 130, 'Sphere · f = x² + y² + z² − r²', 20)
    body += circle(208, 248, 83, '#e1d4ef')
    body += text(208, 243, 'inside', 20, anchor='middle') + text(208, 270, 'f < 0  ·  sense −1', 17, anchor='middle')
    body += text(344, 219, 'outside', 20) + text(324, 247, 'f > 0  ·  sense +1', 17)
    body += text(208, 364, 'On the sphere: f = 0', 17, anchor='middle')
    body += text(574, 130, 'Plane · f = x', 20)
    body += rect(578, 161, 206, 166, '#e1d4ef') + rect(784, 161, 238, 166, '#ffffff')
    body += line(784, 150, 784, 339) + line(784, 186, 879, 186, 'marker-end="url(#arrow)"')
    body += text(894, 191, 'normal', 16)
    body += text(682, 253, 'x < 0', 22, anchor='middle') + text(682, 282, 'sense −1', 18, anchor='middle')
    body += text(911, 253, 'x > 0', 22, anchor='middle') + text(911, 282, 'sense +1', 18, anchor='middle')
    body += text(784, 364, 'On the plane: x = 0', 17, anchor='middle')
    svg('surface-sense.svg', 'Negative and positive surface sense',
        'Sphere: negative inside, positive outside. Plane x=0: negative to the left, positive along the normal.', 420, body)

    body = text(36, 44, 'Boolean operations', 26) + text(36, 73, 'Colored areas show the resulting region. Dashed outlines retain the input shapes.', 18, MUTED)
    for i, (title, expr, note) in enumerate([
        ('Union', 'A ∪ B', 'Inside either region'),
        ('Intersection', 'A ∩ B', 'Inside both regions'),
        ('Difference', 'A ∖ B', 'Inside A, outside B'),
    ]):
        x = 20 + 354 * i
        a, b, cy, r = x + 139, x + 215, 243, 80
        body += rect(x, 98, 332, 307, '#f3f6fa', 12)
        body += text(x + 166, 135, title, 21, anchor='middle')
        body += f'<defs><clipPath id="clip{i}">{circle(b, cy, r, "white", "none")}</clipPath>'
        body += f'<mask id="mask{i}">{rect(x, 150, 332, 185, "white")}{circle(b, cy, r, "black", "none")}</mask></defs>'
        if i == 0:
            body += circle(a, cy, r, PURPLE, 'none') + circle(b, cy, r, PURPLE, 'none')
        else:
            attr = f'clip-path="url(#clip{i})"' if i == 1 else f'mask="url(#mask{i})"'
            body += circle(a, cy, r, PURPLE, 'none', attr)
        body += circle(a, cy, r, 'none', INK, 'stroke-dasharray="5 5"') + circle(b, cy, r, 'none', INK, 'stroke-dasharray="5 5"')
        body += text(a - 35, cy + 7, 'A', 22, anchor='middle') + text(b + 35, cy + 7, 'B', 22, anchor='middle')
        body += text(x + 166, 355, expr, 23, anchor='middle') + text(x + 166, 383, note, 17, anchor='middle')
    svg('boolean-operations.svg', 'Union, intersection, and difference',
        'Two overlapping circular regions A and B illustrate each Boolean operation.', 425, body)

    body = text(36, 44, 'Following a ray through the shield', 26)
    body += text(36, 75, 'One-dimensional view: origin at the center, unit direction +X, stop at the outer surface.', 18, MUTED)
    for start, end, color, cell, mat, label in [(0, 3, ORANGE, 1, 1, 'Source region'), (3, 7, GREEN, 2, 2, 'Inner shield'), (7, 10, PURPLE, 3, 4, 'Outer shield')]:
        x, w = 95 + start * 88, (end - start) * 88
        body += rect(x, 122, w, 159, color)
        body += text(x + w / 2, 152, label, 20, '#14212f', 'middle')
        body += text(x + w / 2, 179, f'Cell {cell} · material {mat}', 17, '#14212f', 'middle')
        body += text(x + w / 2, 262, f'{end-start} cm', 22, '#14212f', 'middle')
    body += line(95, 210, 1010, 210, 'marker-end="url(#arrow)"')
    for t in [0, 3, 7, 10]:
        x = 95 + t * 88
        body += circle(x, 210, 5, '#ffffff') + line(x, 286, x, 298)
        body += text(x, 325, f't = {t}', 18, anchor='middle')
    body += text(540, 371, 'Segment length = t_exit − t_enter   •   t is in cm for this unit direction', 19, anchor='middle')
    svg('ray-segments.svg', 'Ray intersections and cell segment lengths',
        'The ray starts in cell 1: t=0 to 3 cm. Cell 2 spans 3 to 7 cm; cell 3 spans 7 to 10 cm.', 405, body)


def compress_png(path):
    """Losslessly compress the native writer's uncompressed PNG IDAT stream."""
    data = path.read_bytes()
    assert data[:8] == b"\x89PNG\r\n\x1a\n"
    chunks = []
    offset = 8
    while offset < len(data):
        size = struct.unpack_from('>I', data, offset)[0]
        kind = data[offset + 4:offset + 8]
        payload = data[offset + 8:offset + 8 + size]
        chunks.append((kind, payload))
        offset += size + 12
    raw = zlib.decompress(b''.join(payload for kind, payload in chunks if kind == b'IDAT'))
    compressed = zlib.compress(raw, 9)
    assert zlib.decompress(compressed) == raw
    output = bytearray(data[:8])
    written = False
    for kind, payload in chunks:
        if kind == b'IDAT':
            if written:
                continue
            payload = compressed
            written = True
        output.extend(struct.pack('>I', len(payload)) + kind + payload)
        output.extend(struct.pack('>I', zlib.crc32(kind + payload) & 0xffffffff))
    path.write_bytes(output)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    model = ROOT / 'examples/data/doc_shield.i'
    with tempfile.TemporaryDirectory(prefix='alea-doc-figures-') as tmp:
        palette = Path(tmp) / 'palette.txt'
        palette.write_text('1 1.0 0.498039 0.054902\n2 0.172549 0.627451 0.172549\n4 0.580392 0.403922 0.741176\n')
        for name, clips in [('shield-full', []), ('shield-cutaway', ['--clip', '0,0,-1,0'])]:
            subprocess.run([str(ROOT / 'bin/render3d'), str(model), '--width', '800', '--height', '640',
                '--eye', '24,28,23', '--target', '0,0,-1', '--fov', '38', '--background', '1,1,1',
                '--colors', str(palette), '--shadows', '--edges', '--aa', '2',
                *clips, '-o', str(OUT / f'{name}.png')], check=True, cwd=ROOT)
        subprocess.run([str(ROOT / 'bin/mc_plotter'), str(model), 'Z', '0', '-12', '12', '-12', '12',
            '480', str(OUT / 'shield-slice.png'), '--color=materials', '--labels=materials'], check=True, cwd=ROOT)
    for name in ("shield-full.png", "shield-cutaway.png", "shield-slice.png"):
        compress_png(OUT / name)
    diagrams()


if __name__ == '__main__':
    main()
