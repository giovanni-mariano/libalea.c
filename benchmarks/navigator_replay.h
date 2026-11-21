// SPDX-License-Identifier: MPL-2.0
/* Trace finite, collision-free flights with the public particle navigator.
 * Keep event serialization outside the profiler and the replay timer. */
ALEA_VEC_DEFINE(nav_events, alea_nav_event_t);

typedef struct {
    alea_nav_location_t source;
    size_t first, end;
    const char* stop;
    double failure_distance;
    char diagnostic[256];
} nav_flight_t;

static void print_nav_location(const alea_nav_location_t* loc) {
    printf("[%d,%d,%d,%.17g,%llu]", (int)loc->kind, loc->cell_id, loc->material_id,
           loc->density, (unsigned long long)loc->occurrence_key);
}

static void print_nav_diagnostic(const char* message) {
    putchar('"');
    for (const unsigned char* p = (const unsigned char*)message; *p; ++p) {
        if (*p == '"' || *p == '\\') putchar('\\');
        if (*p < 32) printf("\\u%04x", (unsigned)*p);
        else putchar(*p);
    }
    putchar('"');
}

static int nav_validation(const char* mode, alea_nav_validation_mode_t* validation) {
    if (!strcmp(mode, "strict")) *validation = ALEA_NAV_VALIDATE_STRICT;
    else if (!strcmp(mode, "interval")) *validation = ALEA_NAV_VALIDATE_INTERVAL;
    else if (!strcmp(mode, "transport")) *validation = ALEA_NAV_VALIDATE_TRANSPORT;
    else if (!strcmp(mode, "fast")) *validation = ALEA_NAV_VALIDATE_FAST;
    else return -1;
    return 0;
}

static void print_nav_event(const alea_nav_event_t* e) {
    printf("[%d,%.17g,%d,%d,[%.17g,%.17g,%.17g],", (int)e->kind, e->distance,
           e->surface_id, (int)e->boundary_type, e->normal[0], e->normal[1], e->normal[2]);
    print_nav_location(&e->before);
    putchar(',');
    print_nav_location(&e->after);
    printf(",[%.17g,%.17g,%.17g]]", e->position[0], e->position[1], e->position[2]);
}

static alea_ray_navigator_t* replay_create_navigator(alea_system_t* sys, bool regional,
        const alea_bbox_t* domain, alea_nav_validation_mode_t mode) {
    alea_ray_navigator_t* nav = alea_ray_navigator_create(sys);
    if (!nav) return NULL;
    if (alea_ray_navigator_set_validation_mode(nav, mode) ||
        (regional && (alea_prepare_ray_acceleration(sys, domain) ||
                      !sys->ray_batch_index || alea_ray_navigator_use_regional(nav, sys->ray_batch_index)))) {
        fprintf(stderr, "Could not prepare particle replay: %s\n", alea_error());
        alea_ray_navigator_destroy(nav);
        return NULL;
    }
    return nav;
}

static int navigator_replay(const char* path, const char* workload, const char* mode, bool regional, const alea_bbox_t* domain) {
    alea_nav_validation_mode_t validation;
    if (nav_validation(mode, &validation)) return 2;
    size_t count;
    query_t* rays = read_workload(workload, &count);
    if (!rays) return 2;
    int status = 1, failures = 0;
    nav_events_t events = ALEA_VEC_INIT;
    nav_flight_t* flights = calloc(count, sizeof(*flights));
    double source_start = alea_monotonic_seconds();
    double first_event_s = -1;
    mcnp_model_t* model = mcnp_load(path);
    alea_ray_navigator_t* nav = NULL;
    if (!flights || !model || alea_prepare_query_acceleration(model->sys)) goto done;
    nav = replay_create_navigator(model->sys, regional, domain, validation);
    if (!nav) goto done;
#ifdef ALEA_RAY_PROFILE
    if (alea_ray_profile_start(model->sys->universes.count)) goto done;
#endif
    double start = alea_monotonic_seconds();
    for (size_t i = 0; i < count; ++i) {
        nav_flight_t* f = &flights[i];
        const double* v = rays[i].v;
        if (rays[i].kind != 'r') goto done;
        f->first = f->end = events.count;
        f->source = (alea_nav_location_t){.kind = ALEA_NAV_UNRESOLVED, .cell_index = -1, .cell_id = -1};
        if (alea_ray_navigator_restart(nav, v, v+3, &f->source)) {
            fprintf(stderr, "ray %zu restart: %s\n", i, alea_error());
            snprintf(f->diagnostic, sizeof(f->diagnostic), "%s", alea_error());
            f->stop = "restart-error";
            failures++;
            continue;
        }
        double t = 0;
        alea_nav_location_t owner = f->source;
        while (t < v[6]) {
            if (owner.kind != ALEA_NAV_MATERIAL && owner.kind != ALEA_NAV_VOID &&
                owner.kind != ALEA_NAV_GAP) {
                f->stop = "ownership";
                break;
            }
            alea_nav_event_t event;
            if (alea_ray_navigator_advance(nav, INFINITY, v[6]-t, &event)) {
                fprintf(stderr, "ray %zu at %.17g: %s\n", i, t, alea_error());
                snprintf(f->diagnostic, sizeof(f->diagnostic), "%s", alea_error());
                f->failure_distance = t;
                f->stop = "advance-error";
                failures++;
                break;
            }
            if (!isfinite(event.distance) || !(t + event.distance > t)) {
                f->failure_distance = t;
                snprintf(f->diagnostic, sizeof(f->diagnostic), "distance did not advance");
                f->stop = "no-progress";
                failures++;
                break;
            }
            if (first_event_s < 0) first_event_s = alea_monotonic_seconds() - source_start;
            if (alea_vec_push(&events, event, alea_nav_event_t)) goto done;
            t += event.distance;
            owner = event.after;
            if (event.kind == ALEA_NAV_DISTANCE_LIMIT) break;
            if (event.kind != ALEA_NAV_BOUNDARY) {
                f->stop = event.kind == ALEA_NAV_VACUUM ? "vacuum" : "boundary-action";
                break;
            }
        }
        f->end = events.count;
        if (!f->stop) f->stop = "distance-limit";
    }
    double elapsed = alea_monotonic_seconds() - start;
#ifdef ALEA_RAY_PROFILE
    alea_ray_profile_stop();
#endif
    printf("{\"mode\":\"%s\",\"replay_s\":%.9g,\"failures\":%d", mode, elapsed, failures);
#ifdef ALEA_RAY_PROFILE
    const char* names[] = {"other","hierarchy","membership","intersection","recovery","recovery_scan","recovery_cache"};
    printf(",\"stage_s\":{");
    for (int i = 0; i < RAY_STAGES; ++i)
        printf("%s\"%s\":%.9g", i ? "," : "", names[i], alea_ray_profile.seconds[i]);
    double enumeration_s = 0;
    uint64_t enumerations = 0;
    for (size_t i = 0; i < model->sys->universes.count; i++) {
        enumeration_s += alea_ray_profile.universes[i].enumeration_seconds;
        enumerations += alea_ray_profile.universes[i].enumerations;
    }
    printf("},\"owner_enumeration_s\":%.9g,\"owner_enumerations\":%llu",
           enumeration_s, (unsigned long long)enumerations);
#endif
    if (validation == ALEA_NAV_VALIDATE_TRANSPORT) {
        alea_nav_recovery_stats_t stats;
        alea_ray_navigator_recovery_stats(nav, &stats);
        printf(",\"recovered_gaps\":%llu,\"recovered_overlaps\":%llu,\"recovered_cm\":%.17g",
               (unsigned long long)stats.gaps, (unsigned long long)stats.overlaps, stats.distance_cm);
    }
    size_t spans, visits;
    alea_ray_navigator_regional_counts(nav, &spans, &visits);
    printf(",\"producer\":\"%s\",\"regional_spans\":%zu,\"regional_candidate_tree_visits\":%zu",
           regional ? "regional" : "hierarchy", spans, visits);
    printf(",\"first_event_s\":");
    if (first_event_s < 0) printf("null");
    else printf("%.9g", first_event_s);
    printf(",\"peak_rss_bytes\":");
    peak_rss();
    printf(",\"flights\":[");
    for (size_t i = 0; i < count; ++i) {
        const nav_flight_t* f = &flights[i];
        printf("%s{\"stop\":\"%s\",\"source\":", i ? "," : "", f->stop);
        print_nav_location(&f->source);
        if (f->diagnostic[0]) {
            printf(",\"failure_distance\":%.17g,\"diagnostic\":", f->failure_distance);
            print_nav_diagnostic(f->diagnostic);
        }
        printf(",\"ray\":[");
        for (int k = 0; k < 7; ++k) printf("%s%.17g", k ? "," : "", rays[i].v[k]);
        putchar(']');
        printf(",\"events\":[");
        for (size_t j = f->first; j < f->end; ++j) {
            if (j != f->first) putchar(',');
            print_nav_event(&events.data[j]);
        }
        printf("]}");
    }
    puts("]}");
    status = failures != 0;
done:
#ifdef ALEA_RAY_PROFILE
    alea_ray_profile_free();
#endif
    alea_ray_navigator_destroy(nav);
    if (model) mcnp_model_destroy(model);
    alea_vec_free(&events);
    free(flights);
    free(rays);
    return status;
}

/* Replay public calls in order on a parent and an independently movable clone.
 * Emit each input with its result so a failure can be reproduced at its line. */
static int navigator_script(const char* path, const char* script, const char* mode, bool regional, const alea_bbox_t* domain) {
    alea_nav_validation_mode_t validation;
    if (nav_validation(mode, &validation)) return 2;
    FILE* input = fopen(script, "r");
    if (!input) return 2;
    int status = 1;
    mcnp_model_t* model = mcnp_load(path);
    alea_ray_navigator_t* slots[2] = {NULL, NULL};
    if (!model || alea_prepare_query_acceleration(model->sys)) goto done;
    slots[0] = replay_create_navigator(model->sys, regional, domain, validation);
    if (!slots[0]) goto done;
    printf("{\"mode\":\"%s\",\"producer\":\"%s\",\"operations\":[", mode, regional ? "regional" : "hierarchy");
    enum { RESTART, ADVANCE, DIRECTION, CLONE, REFLECT, PERIODIC };
    char line[256];
    unsigned number = 0, records = 0;
    while (fgets(line, sizeof(line), input)) {
        ++number;
        if (!strchr(line, '\n') && !feof(input)) goto malformed;
        const char* delim = " \t\r\n";
        char* name = strtok(line, delim);
        if (!name || name[0] == '#') continue;
        char* who = strtok(NULL, delim);
        if (!who || (strcmp(who, "parent") && strcmp(who, "child"))) goto malformed;
        int slot = !strcmp(who, "child");
        alea_ray_navigator_t* nav = slots[slot];
        if (!nav) goto malformed;
        static const struct { const char* name; int argc; } ops[] = {
            {"restart", 6}, {"advance", 2}, {"direction", 3},
            {"clone", 0}, {"reflect", 0}, {"periodic", 0}
        };
        size_t op = 0;
        while (op < sizeof(ops)/sizeof(*ops) && strcmp(ops[op].name, name)) ++op;
        if (op == sizeof(ops)/sizeof(*ops)) goto malformed;
        double a[6];
        for (int i = 0; i < ops[op].argc; ++i) {
            char* token = strtok(NULL, delim);
            char* end;
            if (!token) goto malformed;
            errno = 0;
            a[i] = strtod(token, &end);
            if (*end || errno || (!isfinite(a[i]) &&
                !(op == ADVANCE && i == 0 && a[i] == INFINITY))) goto malformed;
        }
        if (strtok(NULL, delim)) goto malformed;
        printf("%s{\"op\":\"%s\",\"particle\":\"%s\",\"args\":[",
               records++ ? "," : "", name, who);
        for (int i = 0; i < ops[op].argc; ++i) {
            if (i) putchar(',');
            if (isinf(a[i])) printf("\"inf\"");
            else printf("%.17g", a[i]);
        }
        printf("]");
        int rc;
        alea_nav_location_t location;
        alea_nav_event_t event;
        double position[3];
        switch (op) {
        case RESTART:
            rc = alea_ray_navigator_restart(nav, a, a+3, &location);
            break;
        case ADVANCE:
            rc = alea_ray_navigator_advance(nav, a[0], a[1], &event);
            break;
        case DIRECTION:
            rc = alea_ray_navigator_set_direction(nav, a);
            break;
        case CLONE: {
            alea_ray_navigator_t* copy = alea_ray_navigator_clone(nav);
            rc = copy ? 0 : -1;
            if (copy) {
                alea_ray_navigator_destroy(slots[1-slot]);
                slots[1-slot] = copy;
            }
            break;
        }
        case REFLECT:
            rc = alea_ray_navigator_reflect_specular(nav);
            break;
        default:
            rc = alea_ray_navigator_apply_periodic(nav, position, &location);
            break;
        }
        printf(",\"rc\":%d", rc);
        if (!rc && (op == RESTART || op == PERIODIC)) {
            printf(",\"location\":");
            print_nav_location(&location);
        }
        if (!rc && op == PERIODIC)
            printf(",\"position\":[%.17g,%.17g,%.17g]", position[0], position[1], position[2]);
        if (!rc && op == ADVANCE) {
            printf(",\"event\":");
            print_nav_event(&event);
        }
        if (rc) {
            printf(",\"diagnostic\":");
            print_nav_diagnostic(alea_error());
        }
        putchar('}');
        if (rc) {
            fprintf(stderr, "%s:%u: %s on %s: %s\n", script, number, name, who, alea_error());
            puts("]}");
            goto done;
        }
    }
    puts("]}");
    status = ferror(input) || !records;
    goto done;
malformed:
    fprintf(stderr, "%s:%u: invalid navigator operation\n", script, number);
    status = 2;
done:
    alea_ray_navigator_destroy(slots[1]);
    alea_ray_navigator_destroy(slots[0]);
    if (model) mcnp_model_destroy(model);
    fclose(input);
    return status;
}
