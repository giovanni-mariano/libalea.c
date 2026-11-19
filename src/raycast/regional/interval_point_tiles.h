// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

/* An octant of a spatial leaf often has fewer possible owners than the
 * whole leaf. A solid winner terminates the ordered search only after the
 * point has failed the membership tests of every earlier possible owner. */
static uint32_t point_region_program(interval_index_t* idx, uint32_t code,
                                     interval_box_t box, unsigned* remaining) {
    if (code <= GUARD_ANY || !*remaining) return code;
    --*remaining;
    if (code & GUARD_BRANCH) {
        bool join = !!(code & GUARD_UNION);
        guard_pair_t pair = idx->guards.pairs.data[code & ~(GUARD_BRANCH | GUARD_UNION)];
        uint32_t left = point_region_program(idx, pair.a, box, remaining);
        if (left == (join ? GUARD_ANY : GUARD_EMPTY)) return left;
        uint32_t right = point_region_program(idx, pair.b, box, remaining);
        if (left == pair.a && right == pair.b) return code;
        return guard_join(&idx->guards, left, right, join);
    }
    if (code & GUARD_REF) return code;
    if (code & GUARD_BOX) {
        uint32_t node = (code & ~GUARD_BOX) >> 1;
        interval_box_t bound = interval_box_from(alea_node_bbox_get(&idx->sys->nodes.data[node].bbox));
        bool covers = true;
        for (unsigned k = 0; k < 3; k++) {
            if (box.hi[k] < bound.lo[k] || box.lo[k] > bound.hi[k]) return code & 1;
            if (box.lo[k] < bound.lo[k] || box.hi[k] > bound.hi[k]) covers = false;
        }
        return covers ? !(code & 1) : code;
    }

    const guard_surface_t* surface = idx->surfaces + ((code - 2) >> 1);
    int sign = surface->kind ? guard_box_sign(surface, box) : 0;
    return sign ? ((sign > 0) == (bool)(code & 1)) : code;
}

static int prepare_point_tiles(interval_index_t* idx, spatial_node_t* leaf, alea_bbox_t region) {
    /* Small leaves already have a short packet search. Reserve the eight
     * additional candidate lists for leaves with more point tests to avoid. */
    if (leaf->count < (idx->selective ? 24u : 3u) || leaf->count > (UINT32_MAX >> POINT_OFFSET_BITS)) return 0;
    if (idx->tiles.count >= UINT32_MAX) return -1;

    point_tiles_t tiles = {0};
    interval_box_t bounds = interval_box_from(region);
    for (unsigned k = 0; k < 3; k++) {
        tiles.center[k] = .5 * bounds.lo[k] + .5 * bounds.hi[k];
    }

    for (unsigned t = 0; t < 8; t++) {
        interval_box_t tile = bounds;
        for (unsigned k = 0; k < 3; k++) {
            if (t & (1u << k)) tile.lo[k] = tiles.center[k];
            else tile.hi[k] = tiles.center[k];
        }
        tiles.first[t] = (uint32_t)idx->point_items.count;
        tiles.solid[t] = UINT32_MAX;

        for (uint32_t i = 0; i < leaf->count; i++) {
            const bound_packet_t* packet = idx->packets.data + leaf->packet + i / BOUND_LANES;
            const spatial_patch_t* patch = idx->patches.data + leaf->first + i;
            unsigned lane = i % BOUND_LANES;
            interval_box_t clipped = tile;
            unsigned axes = 0;
            for (unsigned k = 0; k < 3; k++) {
                double lo = packet->xyz[k].lo[lane], hi = packet->xyz[k].hi[lane];
                if (clipped.lo[k] < lo) {
                    clipped.lo[k] = lo;
                    axes |= 1u << k;
                }
                if (clipped.hi[k] > hi) {
                    clipped.hi[k] = hi;
                    axes |= 1u << k;
                }
            }
            if (box_empty(clipped)) continue;

            bool possible = true;
            size_t pos = packet->first;
            for (unsigned k = 0; k < 6; k++) {
                if (!(packet->mask & (1u << k))) continue;
                const packet_axis_t* axis = idx->extra_axes.data + pos++;
                unsigned a = k / 2, b = (a + 1) % 3;
                double lo = clipped.lo[a] + (k % 2 ? -clipped.hi[b] : clipped.lo[b]);
                double hi = clipped.hi[a] + (k % 2 ? -clipped.lo[b] : clipped.hi[b]);
                if (lo > axis->hi[lane] || hi < axis->lo[lane]) {
                    possible = false;
                    break;
                }
                if (lo < axis->lo[lane] || hi > axis->hi[lane]) axes |= 1u << (k + 3);
            }
            if (!possible) continue;

            unsigned budget = 96;
            uint32_t program = point_region_program(idx, patch->program, clipped, &budget);
            if (idx->guards.failed) return -1;
            if (program == GUARD_EMPTY) continue;
            const alea_cell_entry_t* cell = idx->sys->cells.data + patch->cell;
            if (program == GUARD_ANY && !axes && !(cell->lat_type && cell->lat_fill)) {
                tiles.solid[t] = patch->cell;
                break;
            }
            uint32_t bounds_ref = (i << POINT_OFFSET_BITS) | axes;
            if (cell->lat_type && cell->lat_fill) bounds_ref |= POINT_LATTICE;
            point_item_t item = {program, bounds_ref};
            if (idx->point_items.count == UINT32_MAX ||
                alea_vec_push(&idx->point_items, item, point_item_t)) return -1;
        }
        uint32_t first = tiles.first[t];
        uint32_t count = (uint32_t)idx->point_items.count - first;
        tiles.count[t] = count;
        if (idx->selective && count) {
            for (unsigned sibling = 0; sibling < t; sibling++) {
                if (tiles.count[sibling] == count &&
                    !memcmp(idx->point_items.data + first, idx->point_items.data + tiles.first[sibling],
                            count * sizeof(point_item_t))) {
                    tiles.first[t] = tiles.first[sibling];
                    idx->point_items.count = first;
                    break;
                }
            }
        }
    }
    leaf->axis = (uint32_t)idx->tiles.count + 1;
    return alea_vec_push(&idx->tiles, tiles, point_tiles_t);
}
