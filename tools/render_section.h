// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* Section window owned by the inspector. Only inspector_slice runs on the worker. */
enum { SECTION_FOOTER = 58, SECTION_PREVIEW_SIZE = 201, SECTION_MAX_PIXELS_PER_AXIS = 1024 };

static void close_section(inspector_t* p) {
    SDL_AtomicAdd(&p->revision, 1);
    SDL_DestroyTexture(p->slice);
    SDL_DestroyRenderer(p->renderer);
    SDL_DestroyWindow(p->window);
    p->slice = NULL;
    p->renderer = NULL;
    p->window = NULL;
    p->request = 0;
    p->section_error = NULL;
}

static void request_section(inspector_t* p) {
    SDL_AtomicAdd(&p->revision, 1);
    p->section.revision = SDL_AtomicGet(&p->revision);
    p->preview = SDL_max(p->section.width, p->section.height) > SECTION_PREVIEW_SIZE;
    p->request = p->redraw = 1;
    p->section_error = NULL;
}

static SDL_Rect section_area(inspector_t* p) {
    int w, h;
    SDL_GetWindowSize(p->window, &w, &h);
    return (SDL_Rect){0, 0, w, h > SECTION_FOOTER ? h - SECTION_FOOTER : 1};
}

static void size_section(inspector_t* p) {
    SDL_Rect area = section_area(p);
    double scale = fmin(1, SECTION_MAX_PIXELS_PER_AXIS / (double)SDL_max(area.w, area.h));
    p->section.width = SDL_max(1, (int)lround(scale * area.w));
    p->section.height = SDL_max(1, (int)lround(scale * area.h));
}

/* Horizontal runs from the first pick toward the second, with the selected
 * world axis exactly vertical. A parallel segment leaves the plane ambiguous. */
static const char* section_from_picks(const double a[3], const double b[3], int up,
                                     alea_slice_plane_t* plane, double* across, double* rise) {
    int i = (up + 1) % 3, j = (up + 2) % 3;
    double di = b[i]-a[i], dj = b[j]-a[j];
    *across = hypot(di, dj);
    *rise = b[up]-a[up];
    if (*across == 0 && *rise == 0) return "Pick two distinct points";
    if (*across <= 32 * DBL_EPSILON * fabs(*rise))
        return "Picks align. Choose another up axis.";

    double n[3] = {0}, vertical[3] = {0};
    n[i] = dj / *across;
    n[j] = -di / *across;
    vertical[up] = 1;
    alea_slice_view_t view;
    alea_slice_view_init(&view, .5*a[0]+.5*b[0], .5*a[1]+.5*b[1], .5*a[2]+.5*b[2],
                         n[0], n[1], n[2], vertical[0], vertical[1], vertical[2], 0, 0, 0, 0);
    *plane = view.plane;
    return NULL;
}

static void open_section(inspector_t* p, const render_config_t* cfg, section_kind_t kind, int up) {
    const double* pt = p->pick.point;
    section_view_t next = {.kind=kind, .up_axis=up, .cell=p->pick.cell_id};
    p->pick_error = NULL;
    if (kind == SECTION_TWO_PICKS) {
        if (!p->has_previous) { p->pick_error = "Pick two points first"; return; }
        p->pick_error = section_from_picks(p->previous_point, pt, up, &next.plane, &next.pick_du, &next.pick_dv);
        if (p->pick_error) return;
    } else {
        alea_slice_view_t view;
        alea_slice_view_axis(&view, kind, pt[kind], 0, 0, 0, 0);
        next.plane = view.plane;
        memcpy(next.plane.origin, pt, sizeof(next.plane.origin));
    }
    double distance = hypot(hypot(cfg->eye[0]-pt[0], cfg->eye[1]-pt[1]), cfg->eye[2]-pt[2]);
    next.span = cfg->fov > 0 ? 2 * distance * tan(cfg->fov * 0.008726646259971648) : cfg->ortho_height;
    if (!isfinite(next.span) || next.span <= 0) {
        p->pick_error = "Invalid section width";
        return;
    }
    if (!p->window) {
        p->window = SDL_CreateWindow("Alea section", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
                                     641, 481 + SECTION_FOOTER, SDL_WINDOW_RESIZABLE);
        if (p->window) {
            SDL_SetWindowMinimumSize(p->window, 552, 240);
            p->renderer = SDL_CreateRenderer(p->window, -1, SDL_RENDERER_ACCELERATED);
            if (!p->renderer) p->renderer = SDL_CreateRenderer(p->window, -1, SDL_RENDERER_SOFTWARE);
        }
        if (!p->renderer) {
            fprintf(stderr, "Cannot open section: %s\n", SDL_GetError());
            close_section(p);
            p->pick_error = "Cannot open section";
            return;
        }
    }
    if (SDL_GetWindowFlags(p->window) & SDL_WINDOW_MINIMIZED) SDL_RestoreWindow(p->window);
    SDL_RaiseWindow(p->window);
    p->section = next;
    size_section(p);
    if (kind == SECTION_TWO_PICKS) {
        /* Leave a margin around both defining points, including in wide windows. */
        double aspect = (double)p->section.width / p->section.height;
        p->section.span = fmax(next.span, 1.1 * fmax(next.pick_du, fabs(next.pick_dv) * aspect));
    }
    p->initial = p->section;
    request_section(p);
}

section_view_t inspector_section_view(const inspector_t* p) {
    section_view_t view = p->section;
    if (p->preview) {
        double scale = SECTION_PREVIEW_SIZE / (double)SDL_max(view.width, view.height);
        view.width = SDL_max(1, (int)lround(scale * view.width));
        view.height = SDL_max(1, (int)lround(scale * view.height));
    }
    return view;
}

void inspector_section_event(inspector_t* p, Uint32 id, const SDL_Event* e) {
    if (!p->window || id != SDL_GetWindowID(p->window)) return;
    section_view_t* v = &p->section;
    SDL_Rect area = section_area(p);
    double units = v->span / SDL_max(1, area.w);
    double dy = units * area.h;
    double zoom = 1;
    switch (e->type) {
    case SDL_WINDOWEVENT:
        if (e->window.event == SDL_WINDOWEVENT_CLOSE) { close_section(p); return; }
        if (e->window.event == SDL_WINDOWEVENT_EXPOSED) p->redraw = 1;
        if (e->window.event != SDL_WINDOWEVENT_SIZE_CHANGED) return;
        size_section(p);
        break;
    case SDL_MOUSEMOTION:
        if (!(e->motion.state & SDL_BUTTON_LMASK)) return;
        v->u -= e->motion.xrel * units;
        v->v += e->motion.yrel * units;
        break;
    case SDL_MOUSEWHEEL: {
        int steps = e->wheel.y;
        if (e->wheel.direction == SDL_MOUSEWHEEL_FLIPPED) steps = -steps;
        if (!steps) return;
        zoom = steps > 0 ? 0.8 : 1.25;
        break;
    }
    case SDL_KEYDOWN:
        switch (e->key.keysym.sym) {
        case SDLK_ESCAPE: case SDLK_q: close_section(p); return;
        case SDLK_LEFT: v->u -= v->span * .1; break;
        case SDLK_RIGHT: v->u += v->span * .1; break;
        case SDLK_UP: v->v += dy * .1; break;
        case SDLK_DOWN: v->v -= dy * .1; break;
        case SDLK_PAGEUP: v->offset += dy * .05; break;
        case SDLK_PAGEDOWN: v->offset -= dy * .05; break;
        case SDLK_EQUALS: case SDLK_PLUS: case SDLK_KP_PLUS: zoom = .8; break;
        case SDLK_MINUS: case SDLK_KP_MINUS: zoom = 1.25; break;
        case SDLK_r:
            *v = p->initial;
            size_section(p);
            break;
        default: return;
        }
        break;
    default: return;
    }
    if (zoom != 1) {
        double span = v->span * zoom;
        if (span < 1e-10 || span > 1e20) return;
        v->span = span;
    }
    request_section(p);
}

static void draw_section_picks(SDL_Renderer* renderer, SDL_Rect area, const section_view_t* v) {
    if (v->kind != SECTION_TWO_PICKS || v->offset != 0) return;
    double scale = area.w / v->span;
    SDL_RenderSetClipRect(renderer, &area);
    for (int i = 0; i < 2; ++i) {
        double side = i ? .5 : -.5;
        double x = area.w*.5 + (side*v->pick_du - v->u)*scale;
        double y = area.h*.5 - (side*v->pick_dv - v->v)*scale;
        if (x < 0 || x >= area.w || y < 0 || y >= area.h) continue;
        int px = (int)x, py = (int)y;
        SDL_Rect marker = {px-3, py-3, 7, 7};
        SDL_Rect label = {px+5, py-8, 14, 18};
        SDL_SetRenderDrawColor(renderer, 32, 32, 32, 255);
        SDL_RenderFillRect(renderer, &label);
        SDL_RenderDrawLine(renderer, px-2, py, px+2, py);
        SDL_RenderDrawLine(renderer, px, py-2, px, py+2);
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        SDL_RenderDrawRect(renderer, &marker);
        alea_ui_text(renderer, label.x+2, label.y+2, i ? "B" : "A");
    }
    SDL_RenderSetClipRect(renderer, NULL);
}

static void draw_section(inspector_t* p) {
    if (!p->window || !p->redraw) return;
    p->redraw = 0;
    SDL_Rect area = section_area(p);
    const section_view_t* v = &p->section;
    const section_view_t* old = &p->displayed;
    SDL_SetRenderDrawColor(p->renderer, 32, 32, 32, 255);
    SDL_RenderClear(p->renderer);
    /* Move the displayed slice with the drag while the worker resamples it. */
    if (p->slice && old->cell == v->cell && old->offset == v->offset &&
        !memcmp(&old->plane, &v->plane, sizeof(v->plane))) {
        double scale = area.w / v->span;
        double w = old->span * scale, h = w * old->height / old->width;
        double x = area.w * .5 + (old->u - v->u) * scale - w * .5;
        double y = area.h * .5 - (old->v - v->v) * scale - h * .5;
        if (w < 1e6 && h < 1e6 && x < area.w && y < area.h && x+w > 0 && y+h > 0) {
            SDL_FRect dst = {(float)x, (float)y, (float)w, (float)h};
            SDL_RenderSetClipRect(p->renderer, &area);
            SDL_RenderCopyF(p->renderer, p->slice, NULL, &dst);
            SDL_RenderSetClipRect(p->renderer, NULL);
        }
    }
    draw_section_picks(p->renderer, area, v);
    SDL_Rect footer = {0, area.h, area.w, SECTION_FOOTER};
    SDL_SetRenderDrawColor(p->renderer, 220, 220, 216, 255);
    SDL_RenderFillRect(p->renderer, &footer);
    SDL_SetRenderDrawColor(p->renderer, 25, 28, 30, 255);
    alea_ui_text(p->renderer, 8, area.h+8, "Drag / arrows: pan   Wheel / +/-: zoom");
    alea_ui_text(p->renderer, 8, area.h+33, "PgUp/PgDn: step plane   R: reset   Esc: close");
    char title[192];
    char position[96];
    static const char* const names[] = {"YZ", "XZ", "XY"};
    if (v->kind == SECTION_TWO_PICKS) {
        snprintf(position, sizeof(position), "Two picks, %c up  |  Offset %+.6g cm",
                 "XYZ"[v->up_axis], v->offset);
    } else {
        snprintf(position, sizeof(position), "%s  |  %c %.6g cm", names[v->kind], "XYZ"[v->kind],
                 v->plane.origin[v->kind] + v->offset);
    }
    snprintf(title, sizeof(title), "Alea section  |  %s  |  Width %.6g cm  |  Gold: cell %d%s",
             position, v->span, v->cell, p->section_error ? "  |  Section unavailable" :
             old->revision != v->revision || old->width != v->width || old->height != v->height ? "  |  Updating" : "");
    SDL_SetWindowTitle(p->window, title);
    SDL_RenderPresent(p->renderer);
}

int inspector_slice(alea_system_t* sys, const render_config_t* cfg, inspector_t* p,
                    const section_view_t* v, unsigned char** rgb) {
    enum { BAND = 16 };
    int w = v->width, h = v->height;
    if (SDL_AtomicGet(&p->revision) != v->revision) return -1;
    int* grid = malloc(2 * BAND * w * sizeof(*grid));
    unsigned char* image = malloc((size_t)3*w*h);
    double units = v->span / w;
    double top = v->v + units * h / 2;
    if (!grid || !image) goto failed;
    alea_slice_view_t view = {.plane = v->plane, .u_min = v->u-v->span/2, .u_max = v->u+v->span/2};
    for (int i = 0; i < 3; ++i) view.plane.origin[i] += v->offset * view.plane.normal[i];
    for (int y = 0; y < h; y += BAND) {
        int rows = SDL_min(BAND, h-y);
        view.v_min = top-(y+rows)*units;
        view.v_max = top-y*units;
        if (SDL_AtomicGet(&p->revision) != v->revision ||
            alea_find_cells_grid(sys, &view, w, rows, -1, grid, grid+w*rows, NULL)) goto failed;
        for (int j = 0; j < rows; ++j) for (int x = 0; x < w; ++x) {
            int i = (rows-1-j)*w+x, mat = grid[w*rows+i];
            unsigned char* pixel = image + 3*((size_t)(y+j)*w+x);
            if (grid[i] == v->cell) {
                pixel[0] = 255; pixel[1] = 191; pixel[2] = 50;
            } else if (grid[i] >= 0 && mat > 0) {
                float r,g,b;
                render_get_color(mat, RENDER_COLOR_MATERIAL, cfg, &r, &g, &b);
                pixel[0] = 255*r; pixel[1] = 255*g; pixel[2] = 255*b;
            } else pixel[0] = pixel[1] = pixel[2] = 32;
        }
    }
    free(grid);
    *rgb = image;
    return 0;
failed:
    free(image);
    free(grid);
    return -1;
}

void inspector_show_slice(inspector_t* p, const section_view_t* v, unsigned char* rgb) {
    if (!p->window || v->revision != p->section.revision) return;
    SDL_Texture* tex = rgb ? SDL_CreateTexture(p->renderer, SDL_PIXELFORMAT_RGB24,
                              SDL_TEXTUREACCESS_STATIC, v->width, v->height) : NULL;
    if (tex && !SDL_UpdateTexture(tex, NULL, rgb, 3*v->width)) {
        SDL_DestroyTexture(p->slice);
        p->slice = tex;
        p->displayed = *v;
        if (p->preview) {
            p->preview = 0;
            p->request = 1;
        }
    } else {
        SDL_DestroyTexture(tex);
        p->section_error = "Section unavailable";
    }
    p->redraw = 1;
}
