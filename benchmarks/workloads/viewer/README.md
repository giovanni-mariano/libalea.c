<!--
SPDX-FileCopyrightText: 2026 Giovanni MARIANO

SPDX-License-Identifier: MPL-2.0
-->

# SDL viewer replay

Run from the repository root with SDL2 installed. This example renders the
checked-in box with both CPU accelerators and compares every displayed frame.

```sh
make RELEASE=1 modules
make -C tools RELEASE=1 USE_SDL=1 ../bin/interactive_render_replay
out=build/viewer-box
for backend in cpu regional; do
    mkdir -p "$out/$backend"
    SDL_VIDEODRIVER=dummy ALEA_REPLAY_DIR="$out/$backend" \
        bin/interactive_render_replay tests/data/simple_box.mcnp \
        --backend "$backend" --interactive --width 96 --height 72 --threads 2 \
        > "$out/$backend/stdout" 2> "$out/$backend/stderr" || exit 1
done
for frame in 00 01 02 03 04 05 06 07 08 09 10; do
    cmp "$out/cpu/$frame.rgb" "$out/regional/$frame.rgb" || exit 1
done
```

Frames 00–04 and 06–09 are previews. Frames 05 and 10 are full-resolution
refinements. Each run records presentation times in its stderr log.
