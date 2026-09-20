// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* This file is part of pyalea_binding.c — do NOT compile separately.
 * It is #included from pyalea_binding.c to keep a single translation unit.
 *
 * Contents: Transport-style geometry validator bindings.
 *           Mirrors csrc/libalea/src/lua_bind/lua_geo_validator.c and exposes
 *           validate_geometry / validate_geometry_ray / validate_geometry_slice
 *           returning plain dicts (no Matplotlib / JSON-serializable).
 */

#include "alea_geo_validator.h"

/* ============================================================================
 * Option parsing
 * ============================================================================ */

/* Read a long-valued dict entry into *out. Returns -1 on conversion error. */
static int geom_opt_long(PyObject* opts, const char* key, long* out) {
    PyObject* item = PyDict_GetItemString(opts, key); /* borrowed */
    if (!item || item == Py_None) return 0;
    long v = PyLong_AsLong(item);
    if (v == -1 && PyErr_Occurred()) {
        PyErr_Format(PyExc_TypeError, "validator option '%s' must be an int", key);
        return -1;
    }
    *out = v;
    return 0;
}

static int geom_opt_double(PyObject* opts, const char* key, double* out) {
    PyObject* item = PyDict_GetItemString(opts, key); /* borrowed */
    if (!item || item == Py_None) return 0;
    double v = PyFloat_AsDouble(item);
    if (v == -1.0 && PyErr_Occurred()) {
        PyErr_Format(PyExc_TypeError, "validator option '%s' must be a number", key);
        return -1;
    }
    *out = v;
    return 0;
}

/* Set/clear a flag bit from a boolean dict entry. */
static void geom_opt_flag(PyObject* opts, const char* key, unsigned* flags, unsigned bit) {
    PyObject* item = PyDict_GetItemString(opts, key); /* borrowed */
    if (!item || item == Py_None) return;
    if (PyObject_IsTrue(item)) *flags |= bit;
    else                       *flags &= ~bit;
}

/* Fill options from a Python dict (may be NULL/None for defaults).
 * Returns 0 on success, -1 with a Python exception set on failure. */
static int parse_validator_options(PyObject* opts, alea_geom_validator_options_t* o) {
    alea_geom_validator_options_init(o);
    if (!opts || opts == Py_None) return 0;
    if (!PyDict_Check(opts)) {
        PyErr_SetString(PyExc_TypeError, "options must be a dict or None");
        return -1;
    }
    if (PyDict_GetItemString(opts, "hierarchical")) {
        PyErr_SetString(PyExc_ValueError,
                        "'hierarchical' was removed: geometry validation "
                        "always uses occurrence-aware hierarchy traversal");
        return -1;
    }

    long lv;
    lv = o->ray_count;       if (geom_opt_long(opts, "ray_count", &lv) < 0) return -1; o->ray_count = (int)lv;
    lv = o->universe_depth;  if (geom_opt_long(opts, "universe_depth", &lv) < 0) return -1; o->universe_depth = (int)lv;
    lv = (long)o->seed;          if (geom_opt_long(opts, "seed", &lv) < 0) return -1; o->seed = (uint64_t)lv;
    lv = (long)o->max_errors;    if (geom_opt_long(opts, "max_errors", &lv) < 0) return -1; o->max_errors = (size_t)lv;
    lv = (long)o->max_samples_per_signature;
    if (geom_opt_long(opts, "max_samples_per_signature", &lv) < 0) return -1;
    if (lv < 0) {
        PyErr_SetString(PyExc_ValueError,
                        "validator option 'max_samples_per_signature' must be non-negative");
        return -1;
    }
    o->max_samples_per_signature = (size_t)lv;
    lv = (long)o->max_samples_per_curve;
    if (geom_opt_long(opts, "max_samples_per_curve", &lv) < 0) return -1;
    if (lv < 0) {
        PyErr_SetString(PyExc_ValueError,
                        "validator option 'max_samples_per_curve' must be non-negative");
        return -1;
    }
    o->max_samples_per_curve = (size_t)lv;
    lv = (long)o->max_crossings; if (geom_opt_long(opts, "max_crossings", &lv) < 0) return -1; o->max_crossings = (size_t)lv;

    if (geom_opt_double(opts, "sample_offset", &o->sample_offset) < 0) return -1;
    if (geom_opt_double(opts, "t_max", &o->t_max) < 0) return -1;

    geom_opt_flag(opts, "strict", &o->flags, ALEA_GEOM_VALIDATE_STRICT_ADJACENCY);
    geom_opt_flag(opts, "allow_exterior_void", &o->flags, ALEA_GEOM_VALIDATE_ALLOW_EXTERIOR_VOID);

    /* An explicit world-space AABB is what turns unowned space into an
     * interior-gap finding; without it, unowned space cannot be classified. */
    PyObject* bounds = PyDict_GetItemString(opts, "domain_bounds"); /* borrowed */
    if (bounds && bounds != Py_None) {
        PyObject* tuple = PySequence_Tuple(bounds);
        int ok = tuple && PyArg_ParseTuple(tuple, "dddddd",
                                           &o->validation_bounds[0], &o->validation_bounds[1],
                                           &o->validation_bounds[2], &o->validation_bounds[3],
                                           &o->validation_bounds[4], &o->validation_bounds[5]);
        Py_XDECREF(tuple);
        if (!ok) {
            PyErr_Clear();
            PyErr_SetString(PyExc_TypeError, "domain_bounds must be a sequence of six floats "
                                             "(min_x, max_x, min_y, max_y, min_z, max_z)");
            return -1;
        }
        o->flags |= ALEA_GEOM_VALIDATE_DOMAIN_BOUNDS;
    }
    return 0;
}

/* ============================================================================
 * Result marshalling
 * ============================================================================ */

static PyObject* geom_vec3(const double v[3]) {
    return Py_BuildValue("(ddd)", v[0], v[1], v[2]);
}

/* Decode the event-flags bitfield into a list of set flag names. */
static PyObject* geom_flag_names(uint32_t flags) {
    static const struct { unsigned bit; const char* name; } table[] = {
        {ALEA_GEOM_EVENT_INITIAL_POINT,           "initial_point"},
        {ALEA_GEOM_EVENT_PREVIOUS_NO_SURFACE,     "previous_no_surface"},
        {ALEA_GEOM_EVENT_MISSING_ADJACENCY,       "missing_adjacency"},
        {ALEA_GEOM_EVENT_EXTERIOR_ALLOWED,        "exterior_allowed"},
        {ALEA_GEOM_EVENT_COINCIDENT_SURFACES,     "coincident_surfaces"},
        {ALEA_GEOM_EVENT_TRUNCATED_COVERAGE,      "truncated_coverage"},
        {ALEA_GEOM_EVENT_FOUND_WITHOUT_ADJACENCY, "found_without_adjacency"},
        {ALEA_GEOM_EVENT_VIEWPORT_EDGE,           "viewport_edge"},
    };
    PyObject* names = PyList_New(0);
    if (!names) return NULL;
    for (size_t i = 0; i < sizeof(table) / sizeof(table[0]); i++) {
        if (flags & table[i].bit) {
            PyObject* s = PyUnicode_FromString(table[i].name);
            if (!s || PyList_Append(names, s) < 0) { Py_XDECREF(s); Py_DECREF(names); return NULL; }
            Py_DECREF(s);
        }
    }
    return names;
}

static const char* geom_source_name(alea_geom_event_source_t source) {
    switch (source) {
        case ALEA_GEOM_EVENT_SOURCE_RAY:           return "ray";
        case ALEA_GEOM_EVENT_SOURCE_SLICE_CURVE:   return "slice_curve";
        case ALEA_GEOM_EVENT_SOURCE_INITIAL_POINT: return "initial_point";
        default:                                   return "unknown";
    }
}

static PyObject* build_error_dict(const alea_geom_error_t* e) {
    PyObject* flag_names = geom_flag_names(e->flags);
    if (!flag_names) return NULL;
    PyObject* d = Py_BuildValue(
        "{s:s,s:i,s:s,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:i,s:N,s:N,s:N,s:d,s:d,s:I,s:N,s:k,s:N}",
        "type", alea_geom_error_type_name(e->type),
        "type_id", (int)e->type,
        "source", geom_source_name(e->source),
        "source_id", (int)e->source,
        "previous_cell_id", e->previous_cell_id,
        "found_cell_id", e->found_cell_id,
        "expected_neighbor_cell_id", e->expected_neighbor_cell_id,
        "secondary_cell_id", e->secondary_cell_id,
        "found_cell_count", e->found_cell_count,
        "surface_id", e->surface_id,
        "universe_id", e->universe_id,
        "universe_depth", e->universe_depth,
        "crossing_point", geom_vec3(e->crossing_point),
        "sample_point", geom_vec3(e->sample_point),
        "direction", geom_vec3(e->direction),
        "t", e->t,
        "offset", e->offset,
        "component_index", (unsigned)e->component_index,
        "uv", Py_BuildValue("(dd)", e->uv[0], e->uv[1]),
        "flags", (unsigned long)e->flags,
        "flag_names", flag_names);
    if (!d) return NULL;

    /* Sentinel fields surface as None rather than huge magic numbers. */
    PyObject* curve_index = (e->curve_index == (size_t)SIZE_MAX)
        ? (Py_INCREF(Py_None), Py_None)
        : PyLong_FromSize_t(e->curve_index);
    if (dict_set_new(d, "curve_index", curve_index) < 0) { Py_DECREF(d); return NULL; }

    PyObject* primitive_id = (e->primitive_id == ALEA_PRIMITIVE_ID_INVALID)
        ? (Py_INCREF(Py_None), Py_None)
        : PyLong_FromUnsignedLong(e->primitive_id);
    if (dict_set_new(d, "primitive_id", primitive_id) < 0) { Py_DECREF(d); return NULL; }

    return d;
}

/* Convert a populated C result into the documented dict, then free nothing
 * (caller owns the C result). Returns a new reference or NULL on error. */
static PyObject* build_result_dict(const alea_geom_validator_result_t* r) {
    size_t n = alea_geom_validator_error_count(r);
    PyObject* errors = PyList_New(n);
    if (!errors) return NULL;
    for (size_t i = 0; i < n; i++) {
        alea_geom_error_t e;
        if (alea_geom_validator_error_get(r, i, &e) != 0) {
            Py_DECREF(errors);
            PyErr_SetString(PyExc_RuntimeError, "failed to read validator error");
            return NULL;
        }
        PyObject* item = build_error_dict(&e);
        if (!item) { Py_DECREF(errors); return NULL; }
        PyList_SET_ITEM(errors, i, item); /* steals ref */
    }
    return Py_BuildValue(
        "{s:N,s:k,s:k,s:k,s:k,s:k,s:k,s:O}",
        "errors", errors,
        "crossings_checked", (unsigned long)r->crossings_checked,
        "adjacency_hits", (unsigned long)r->adjacency_hits,
        "exact_queries", (unsigned long)r->exact_queries,
        "ambiguous_crossings", (unsigned long)r->ambiguous_crossings,
        "suppressed_samples", (unsigned long)r->suppressed_samples,
        "sample_limited_curves", (unsigned long)r->sample_limited_curves,
        "truncated", r->truncated ? Py_True : Py_False);
}

/* ============================================================================
 * PyAleaSystem Methods - Geometry Validator
 * ============================================================================ */

static PyObject* PyAleaSystem_validate_geometry(PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    PyObject* opts = NULL;
    static char* kwlist[] = {"options", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &opts)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }
    if (ensure_query_acceleration(self) < 0) return NULL;

    alea_geom_validator_options_t options;
    if (parse_validator_options(opts, &options) < 0) return NULL;

    alea_geom_validator_result_t result;
    alea_geom_validator_result_init(&result);

    int rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    rc = alea_validate_geometry(self->sys, &options, &result);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) { alea_geom_validator_result_free(&result); return NULL; }

    if (rc != 0) {
        alea_geom_validator_result_free(&result);
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    PyObject* out = build_result_dict(&result);
    alea_geom_validator_result_free(&result);
    return out;
}

static PyObject* PyAleaSystem_validate_geometry_ray(PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    double ox, oy, oz, dx, dy, dz, t_max = 0.0;
    PyObject* opts = NULL;
    static char* kwlist[] = {"ox", "oy", "oz", "dx", "dy", "dz", "t_max", "options", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "dddddd|dO", kwlist,
                                     &ox, &oy, &oz, &dx, &dy, &dz, &t_max, &opts)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }
    if (ensure_query_acceleration(self) < 0) return NULL;

    alea_geom_validator_options_t options;
    if (parse_validator_options(opts, &options) < 0) return NULL;

    alea_geom_validator_result_t result;
    alea_geom_validator_result_init(&result);

    int rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    rc = alea_validate_geometry_ray(self->sys, &options, ox, oy, oz, dx, dy, dz, t_max, &result);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) { alea_geom_validator_result_free(&result); return NULL; }

    if (rc != 0) {
        alea_geom_validator_result_free(&result);
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    PyObject* out = build_result_dict(&result);
    alea_geom_validator_result_free(&result);
    return out;
}

static PyObject* PyAleaSystem_validate_geometry_slice(PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    PyObject *origin_obj, *normal_obj, *up_obj;
    double u_min, u_max, v_min, v_max;
    PyObject* opts = NULL;
    static char* kwlist[] = {"origin", "normal", "up", "u_min", "u_max", "v_min", "v_max", "options", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOdddd|O", kwlist,
            &origin_obj, &normal_obj, &up_obj, &u_min, &u_max, &v_min, &v_max, &opts)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }
    if (ensure_query_acceleration(self) < 0) return NULL;

    double ox, oy, oz, nx, ny, nz, ux, uy, uz;
    if (!PyArg_ParseTuple(origin_obj, "ddd", &ox, &oy, &oz)) {
        PyErr_SetString(PyExc_TypeError, "origin must be (x, y, z) tuple"); return NULL;
    }
    if (!PyArg_ParseTuple(normal_obj, "ddd", &nx, &ny, &nz)) {
        PyErr_SetString(PyExc_TypeError, "normal must be (nx, ny, nz) tuple"); return NULL;
    }
    if (!PyArg_ParseTuple(up_obj, "ddd", &ux, &uy, &uz)) {
        PyErr_SetString(PyExc_TypeError, "up must be (ux, uy, uz) tuple"); return NULL;
    }

    alea_slice_view_t view;
    alea_slice_view_init(&view, ox, oy, oz, nx, ny, nz, ux, uy, uz, u_min, u_max, v_min, v_max);

    alea_geom_validator_options_t options;
    if (parse_validator_options(opts, &options) < 0) return NULL;

    /* The validator now consumes caller-provided slice curves so it can carry
     * curve_index/t/uv provenance back in each event. */
    alea_slice_curves_t* curves =
        get_slice_curves_allow_threads(self->sys, &view);
    if (PyErr_Occurred()) return NULL;

    alea_geom_validator_result_t result;
    alea_geom_validator_result_init(&result);

    int rc = 0;
    if (curves) {
        sighandler_func old_sigint = install_sigint();
        Py_BEGIN_ALLOW_THREADS
        rc = alea_validate_geometry_slice(self->sys, &view, curves, &options, &result);
        Py_END_ALLOW_THREADS
        if (restore_sigint(old_sigint)) {
            alea_slice_curves_free(curves);
            alea_geom_validator_result_free(&result);
            return NULL;
        }
        alea_slice_curves_free(curves);
    }
    /* No curves: no boundaries on this plane, nothing to validate (empty result). */

    if (rc != 0) {
        alea_geom_validator_result_free(&result);
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    PyObject* out = build_result_dict(&result);
    alea_geom_validator_result_free(&result);
    return out;
}

static int transition_parse_vec3(PyObject* object, const char* name,
                                 double out[3]) {
    PyObject* seq = PySequence_Fast(object, name);
    if (!seq) return -1;
    if (PySequence_Fast_GET_SIZE(seq) != 3) {
        Py_DECREF(seq);
        PyErr_Format(PyExc_ValueError, "%s must contain exactly 3 values", name);
        return -1;
    }
    for (Py_ssize_t i = 0; i < 3; i++) {
        out[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq, i));
        if (PyErr_Occurred()) {
            Py_DECREF(seq);
            PyErr_Format(PyExc_TypeError, "%s values must be numeric", name);
            return -1;
        }
    }
    Py_DECREF(seq);
    return 0;
}

static PyObject* build_transition_dict(const alea_transition_result_t* result) {
    PyObject* candidate_ids = PyList_New((Py_ssize_t)result->candidate_cell_count);
    PyObject* owner_ids = PyList_New((Py_ssize_t)result->owner_cell_count);
    if (!candidate_ids || !owner_ids) {
        Py_XDECREF(candidate_ids);
        Py_XDECREF(owner_ids);
        return NULL;
    }
    for (size_t i = 0; i < result->candidate_cell_count; i++) {
        PyObject* value = PyLong_FromLong(result->candidate_cell_ids[i]);
        if (!value) { Py_DECREF(candidate_ids); Py_DECREF(owner_ids); return NULL; }
        PyList_SET_ITEM(candidate_ids, (Py_ssize_t)i, value);
    }
    for (size_t i = 0; i < result->owner_cell_count; i++) {
        PyObject* value = PyLong_FromLong(result->owner_cell_ids[i]);
        if (!value) { Py_DECREF(candidate_ids); Py_DECREF(owner_ids); return NULL; }
        PyList_SET_ITEM(owner_ids, (Py_ssize_t)i, value);
    }

    PyObject* out = Py_BuildValue(
        "{s:s,s:s,s:I,s:i,s:i,s:i,s:i,s:i,s:i,s:d,"
        "s:n,s:n,s:n,s:n,s:n,s:N,s:N,s:N,s:N,s:N,s:N}",
        "kind", alea_transition_kind_name(result->kind),
        "after_coverage_kind", point_coverage_kind_name(result->after_coverage_kind),
        "flags", result->flags,
        "universe_id", result->universe_id,
        "current_cell_id", result->current_cell_id,
        "primary_surface_id", result->primary_surface_id,
        "connecting_surface_id", result->connecting_surface_id,
        "after_cell_id", result->after_cell_id,
        "current_sense", result->current_sense,
        "probe_distance", result->probe_distance,
        "offset_attempts", (Py_ssize_t)result->offset_attempts,
        "coverage_fallbacks", (Py_ssize_t)result->coverage_fallbacks,
        "primary_candidate_count", (Py_ssize_t)result->primary_candidate_count,
        "primary_containing_count", (Py_ssize_t)result->primary_containing_count,
        "after_owner_count", (Py_ssize_t)result->after_owner_count,
        "candidate_cell_ids", candidate_ids,
        "owner_cell_ids", owner_ids,
        "point", geom_vec3(result->crossing_point),
        "direction", geom_vec3(result->direction),
        "before_point", geom_vec3(result->before_point),
        "after_point", geom_vec3(result->after_point));
    if (!out) return NULL;

#define SET_TRANSITION_INT(KEY, VALUE) \
    do { \
        if (dict_set_new(out, KEY, PyLong_FromLong((long)(VALUE))) < 0) { \
            Py_DECREF(out); return NULL; \
        } \
    } while (0)
#define SET_TRANSITION_KEY(KEY, VALUE) \
    do { \
        if (dict_set_new(out, KEY, PyLong_FromUnsignedLongLong( \
                (unsigned long long)(VALUE))) < 0) { \
            Py_DECREF(out); return NULL; \
        } \
    } while (0)
    SET_TRANSITION_INT("occurrence_depth", result->occurrence_depth);
    SET_TRANSITION_KEY("current_occurrence_key", result->current_occurrence_key);
    SET_TRANSITION_KEY("current_parent_occurrence_key",
                       result->current_parent_occurrence_key);
    SET_TRANSITION_KEY("before_occurrence_key", result->before_occurrence_key);
    SET_TRANSITION_KEY("before_parent_occurrence_key",
                       result->before_parent_occurrence_key);
    SET_TRANSITION_KEY("selected_after_occurrence_key",
                       result->selected_after_occurrence_key);
    SET_TRANSITION_KEY("selected_after_parent_occurrence_key",
                       result->selected_after_parent_occurrence_key);
#undef SET_TRANSITION_INT
#undef SET_TRANSITION_KEY
    return out;
}

static int transition_slice_opt_u64(PyObject* opts, const char* key,
                                    uint64_t* out) {
    PyObject* item = PyDict_GetItemString(opts, key);
    if (!item || item == Py_None) return 0;
    unsigned long long value = PyLong_AsUnsignedLongLong(item);
    if (PyErr_Occurred()) {
        PyErr_Clear();
        PyErr_Format(PyExc_TypeError,
                     "transition slice option '%s' must be a non-negative integer",
                     key);
        return -1;
    }
    *out = (uint64_t)value;
    return 0;
}

static int parse_transition_slice_options(
    PyObject* opts, alea_transition_slice_options_t* options) {
    alea_transition_slice_options_init(options);
    if (!opts || opts == Py_None) return 0;
    if (!PyDict_Check(opts)) {
        PyErr_SetString(PyExc_TypeError, "options must be a dict or None");
        return -1;
    }
#define PARSE_SCREEN_U64(FIELD) \
    do { \
        uint64_t value = (uint64_t)options->FIELD; \
        if (transition_slice_opt_u64(opts, #FIELD, &value) < 0) return -1; \
        if (sizeof(options->FIELD) < sizeof(value) && value > (uint64_t)SIZE_MAX) { \
            PyErr_Format(PyExc_OverflowError, \
                         "transition slice option '%s' does not fit this platform", \
                         #FIELD); \
            return -1; \
        } \
        options->FIELD = value; \
    } while (0)
    PARSE_SCREEN_U64(horizontal_rays);
    PARSE_SCREEN_U64(vertical_rays);
    PARSE_SCREEN_U64(max_rays);
    PARSE_SCREEN_U64(max_events);
    PARSE_SCREEN_U64(max_events_per_ray);
    PARSE_SCREEN_U64(max_findings);
    PARSE_SCREEN_U64(max_components);
    PARSE_SCREEN_U64(max_output_bytes);
    PARSE_SCREEN_U64(max_scratch_bytes);
    PARSE_SCREEN_U64(max_coverage_fallbacks);
    PARSE_SCREEN_U64(max_coverage_hits);
    PARSE_SCREEN_U64(max_row_scratch_bytes);
    PARSE_SCREEN_U64(coverage_uniform_probes_per_ray);
    PARSE_SCREEN_U64(max_coverage_probes);
    PARSE_SCREEN_U64(max_coverage_findings);
    PARSE_SCREEN_U64(max_coverage_components);
    PARSE_SCREEN_U64(max_component_links);
    PARSE_SCREEN_U64(max_refinement_frontiers);
    PARSE_SCREEN_U64(max_critical_tiles);
    PARSE_SCREEN_U64(max_critical_tile_sources);
    PARSE_SCREEN_U64(max_critical_scratch_bytes);
    PARSE_SCREEN_U64(max_curves_per_tile);
    PARSE_SCREEN_U64(max_critical_points);
    PARSE_SCREEN_U64(max_active_boundary_tests);
    PARSE_SCREEN_U64(max_critical_probes);
    PARSE_SCREEN_U64(max_critical_findings);
    PARSE_SCREEN_U64(max_critical_boundary_evidence);
    PARSE_SCREEN_U64(max_curve_pairs);
    PARSE_SCREEN_U64(max_critical_sector_witnesses);
    PARSE_SCREEN_U64(max_exhaustive_occurrence_hits);
#undef PARSE_SCREEN_U64
    uint64_t u64 = options->max_refinement_depth;
    if (transition_slice_opt_u64(opts, "max_refinement_depth", &u64) < 0)
        return -1;
    if (u64 > UINT32_MAX) {
        PyErr_SetString(PyExc_OverflowError,
                        "max_refinement_depth exceeds uint32 range");
        return -1;
    }
    options->max_refinement_depth = (uint32_t)u64;
    u64 = options->refine_signals;
    if (transition_slice_opt_u64(opts, "refine_signals", &u64) < 0)
        return -1;
    if (u64 > UINT32_MAX) {
        PyErr_SetString(PyExc_OverflowError,
                        "refine_signals exceeds uint32 range");
        return -1;
    }
    options->refine_signals = (uint32_t)u64;
    if (options->max_coverage_hits == 0 ||
        options->max_coverage_hits > 16384) {
        PyErr_SetString(PyExc_ValueError,
                        "max_coverage_hits must be between 1 and 16384");
        return -1;
    }
    if (geom_opt_double(opts, "probe_distance", &options->probe_distance) < 0 ||
        geom_opt_double(opts, "max_probe_distance",
                        &options->max_probe_distance) < 0 ||
        geom_opt_double(opts, "min_transverse_spacing",
                        &options->min_transverse_spacing) < 0 ||
        geom_opt_double(opts, "critical_tile_padding",
                        &options->critical_tile_padding) < 0 ||
        geom_opt_double(opts, "critical_probe_radius",
                        &options->critical_probe_radius) < 0)
        return -1;
    PyObject* include_void =
        PyDict_GetItemString(opts, "include_void_transitions");
    if (include_void && include_void != Py_None) {
        int truth = PyObject_IsTrue(include_void);
        if (truth < 0) return -1;
        options->include_void_transitions = truth;
    }
    PyObject* probe_intervals =
        PyDict_GetItemString(opts, "coverage_probe_selected_intervals");
    if (probe_intervals && probe_intervals != Py_None) {
        int truth = PyObject_IsTrue(probe_intervals);
        if (truth < 0) return -1;
        options->coverage_probe_selected_intervals = truth;
    }
    PyObject* report_unowned =
        PyDict_GetItemString(opts, "report_unowned_coverage");
    if (report_unowned && report_unowned != Py_None) {
        int truth = PyObject_IsTrue(report_unowned);
        if (truth < 0) return -1;
        options->report_unowned_coverage = truth;
    }
    PyObject* critical =
        PyDict_GetItemString(opts, "enable_critical_refinement");
    if (critical && critical != Py_None) {
        int truth = PyObject_IsTrue(critical);
        if (truth < 0) return -1;
        options->enable_critical_refinement = truth;
    }
    PyObject* full_view = PyDict_GetItemString(opts, "critical_full_view");
    if (full_view && full_view != Py_None) {
        int truth = PyObject_IsTrue(full_view);
        if (truth < 0) return -1;
        options->critical_full_view = truth;
    }
    PyObject* occurrence_discovery =
        PyDict_GetItemString(opts, "occurrence_discovery");
    if (occurrence_discovery && occurrence_discovery != Py_None) {
        if (!PyUnicode_Check(occurrence_discovery)) {
            PyErr_SetString(PyExc_TypeError,
                            "occurrence_discovery must be 'sampled' or 'exhaustive'");
            return -1;
        }
        const char* value = PyUnicode_AsUTF8(occurrence_discovery);
        if (!value) return -1;
        if (strcmp(value, "sampled") == 0)
            options->occurrence_discovery =
                ALEA_TRANSITION_SLICE_OCCURRENCE_SAMPLED;
        else if (strcmp(value, "exhaustive") == 0)
            options->occurrence_discovery =
                ALEA_TRANSITION_SLICE_OCCURRENCE_EXHAUSTIVE;
        else {
            PyErr_SetString(PyExc_ValueError,
                            "occurrence_discovery must be 'sampled' or 'exhaustive'");
            return -1;
        }
    }
    return 0;
}

static PyObject* build_transition_slice_coverage_finding_dict(
    const alea_transition_slice_coverage_finding_t* finding) {
    PyObject* owners = PyList_New((Py_ssize_t)finding->owner_count);
    if (!owners) return NULL;
    for (size_t i = 0; i < finding->owner_count; i++) {
        PyObject* owner = Py_BuildValue(
            "{s:i,s:i,s:i,s:K,s:K}",
            "cell_id", finding->owner_cell_ids[i],
            "universe_id", finding->owner_universe_ids[i],
            "depth", finding->owner_depths[i],
            "occurrence_key",
                (unsigned long long)finding->owner_occurrence_keys[i],
            "parent_occurrence_key",
                (unsigned long long)finding->owner_parent_occurrence_keys[i]);
        if (!owner) { Py_DECREF(owners); return NULL; }
        PyList_SET_ITEM(owners, (Py_ssize_t)i, owner);
    }
    PyObject* base_ray_index = finding->base_ray_index == SIZE_MAX
        ? (Py_INCREF(Py_None), Py_None)
        : PyLong_FromSize_t(finding->base_ray_index);
    if (!base_ray_index) { Py_DECREF(owners); return NULL; }
    PyObject* out = Py_BuildValue(
        "{s:s,s:O,s:i,s:n,s:n,s:N,s:s,s:n,s:N,s:I,s:d,s:d,s:d,s:d,s:N,s:N}",
        "kind", point_coverage_kind_name(finding->kind),
        "truncated", finding->truncated ? Py_True : Py_False,
        "target_depth", finding->target_depth,
        "owner_count", (Py_ssize_t)finding->owner_count,
        "owner_count_lower_bound",
            (Py_ssize_t)finding->owner_count_lower_bound,
        "owners", owners,
        "orientation",
            finding->orientation == ALEA_TRANSITION_SLICE_HORIZONTAL
                ? "horizontal" : "vertical",
        "ray_index", (Py_ssize_t)finding->ray_index,
        "base_ray_index", base_ray_index,
        "refinement_depth", finding->refinement_depth,
        "transverse_coordinate", finding->transverse_coordinate,
        "ray_t", finding->ray_t,
        "bracket_t_enter", finding->bracket_t_enter,
        "bracket_t_exit", finding->bracket_t_exit,
        "uv", Py_BuildValue("(dd)", finding->uv[0], finding->uv[1]),
        "world_point", geom_vec3(finding->world_point));
    return out;
}

static const char* transition_slice_tile_source_kind_name(
    alea_transition_slice_tile_source_kind_t kind) {
    switch (kind) {
    case ALEA_TRANSITION_SLICE_TILE_SOURCE_TRANSITION_COMPONENT:
        return "transition_component";
    case ALEA_TRANSITION_SLICE_TILE_SOURCE_COVERAGE_COMPONENT:
        return "coverage_component";
    case ALEA_TRANSITION_SLICE_TILE_SOURCE_REFINEMENT_FRONTIER:
        return "refinement_frontier";
    case ALEA_TRANSITION_SLICE_TILE_SOURCE_FULL_VIEW:
        return "full_view";
    }
    return "unknown";
}

static PyObject* transition_slice_result_to_py(
    alea_transition_slice_result_t* result) {
    size_t count = alea_transition_slice_finding_count(result);
    if (count > (size_t)PY_SSIZE_T_MAX) {
        alea_transition_slice_result_destroy(result);
        return PyErr_NoMemory();
    }
    PyObject* findings = PyList_New((Py_ssize_t)count);
    if (!findings) {
        alea_transition_slice_result_destroy(result);
        return NULL;
    }
    for (size_t i = 0; i < count; i++) {
        alea_transition_slice_finding_t finding;
        if (alea_transition_slice_finding_get(result, i, &finding) != 0) {
            Py_DECREF(findings);
            alea_transition_slice_result_destroy(result);
            PyErr_SetString(PyExc_RuntimeError, "failed to read transition finding");
            return NULL;
        }
        PyObject* item = build_transition_dict(&finding.transition);
        PyObject* base_ray_index = finding.base_ray_index == SIZE_MAX
            ? (Py_INCREF(Py_None), Py_None)
            : PyLong_FromSize_t(finding.base_ray_index);
        if (!item) {
            Py_DECREF(base_ray_index);
            Py_DECREF(findings);
            alea_transition_slice_result_destroy(result);
            return NULL;
        }
        if (dict_set_new(item, "base_ray_index", base_ray_index) < 0 ||
            dict_set_new(item, "orientation", PyUnicode_FromString(
                finding.orientation == ALEA_TRANSITION_SLICE_HORIZONTAL
                    ? "horizontal" : "vertical")) < 0 ||
            dict_set_new(item, "ray_index", PyLong_FromSize_t(finding.ray_index)) < 0 ||
            dict_set_new(item, "event_index", PyLong_FromSize_t(finding.event_index)) < 0 ||
            dict_set_new(item, "refinement_depth",
                         PyLong_FromUnsignedLong(finding.refinement_depth)) < 0 ||
            dict_set_new(item, "transverse_coordinate",
                         PyFloat_FromDouble(finding.transverse_coordinate)) < 0 ||
            dict_set_new(item, "ray_t", PyFloat_FromDouble(finding.ray_t)) < 0 ||
            dict_set_new(item, "uv", Py_BuildValue("(dd)", finding.uv[0], finding.uv[1])) < 0 ||
            dict_set_new(item, "world_point", geom_vec3(finding.world_point)) < 0) {
            Py_XDECREF(item);
            Py_DECREF(findings);
            alea_transition_slice_result_destroy(result);
            return NULL;
        }
        PyList_SET_ITEM(findings, (Py_ssize_t)i, item);
    }
    size_t component_count = alea_transition_slice_component_count(result);
    if (component_count > (size_t)PY_SSIZE_T_MAX) {
        Py_DECREF(findings);
        alea_transition_slice_result_destroy(result);
        return PyErr_NoMemory();
    }
    PyObject* components = PyList_New((Py_ssize_t)component_count);
    if (!components) {
        Py_DECREF(findings);
        alea_transition_slice_result_destroy(result);
        return NULL;
    }
    for (size_t i = 0; i < component_count; i++) {
        alea_transition_slice_component_t component;
        if (alea_transition_slice_component_get(result, i, &component) != 0) {
            Py_DECREF(components);
            Py_DECREF(findings);
            alea_transition_slice_result_destroy(result);
            PyErr_SetString(PyExc_RuntimeError,
                            "failed to read transition component");
            return NULL;
        }
        PyObject* item = Py_BuildValue(
            "{s:s,s:s,s:i,s:i,s:i,s:i,s:i,s:K,s:n,s:n,s:I,s:N,s:N,s:N,s:N}",
            "kind", alea_transition_kind_name(component.kind),
            "orientation",
                component.orientation == ALEA_TRANSITION_SLICE_HORIZONTAL
                    ? "horizontal" : "vertical",
            "universe_id", component.universe_id,
            "current_cell_id", component.current_cell_id,
            "after_cell_id", component.after_cell_id,
            "primary_surface_id", component.primary_surface_id,
            "connecting_surface_id", component.connecting_surface_id,
            "current_occurrence_key",
                (unsigned long long)component.current_occurrence_key,
            "first_finding_index", (Py_ssize_t)component.first_finding_index,
            "finding_count", (Py_ssize_t)component.finding_count,
            "max_refinement_depth", component.max_refinement_depth,
            "uv_min", Py_BuildValue("(dd)", component.uv_min[0],
                                     component.uv_min[1]),
            "uv_max", Py_BuildValue("(dd)", component.uv_max[0],
                                     component.uv_max[1]),
            "world_min", geom_vec3(component.world_min),
            "world_max", geom_vec3(component.world_max));
        if (!item) {
            Py_DECREF(components);
            Py_DECREF(findings);
            alea_transition_slice_result_destroy(result);
            return NULL;
        }
        PyList_SET_ITEM(components, (Py_ssize_t)i, item);
    }
    size_t coverage_count =
        alea_transition_slice_coverage_finding_count(result);
    if (coverage_count > (size_t)PY_SSIZE_T_MAX) {
        Py_DECREF(components); Py_DECREF(findings);
        alea_transition_slice_result_destroy(result);
        return PyErr_NoMemory();
    }
    PyObject* coverage_findings = PyList_New((Py_ssize_t)coverage_count);
    if (!coverage_findings) {
        Py_DECREF(components); Py_DECREF(findings);
        alea_transition_slice_result_destroy(result);
        return NULL;
    }
    for (size_t i = 0; i < coverage_count; i++) {
        alea_transition_slice_coverage_finding_t finding;
        if (alea_transition_slice_coverage_finding_get(
                result, i, &finding) != 0) {
            Py_DECREF(coverage_findings); Py_DECREF(components);
            Py_DECREF(findings);
            alea_transition_slice_result_destroy(result);
            PyErr_SetString(PyExc_RuntimeError,
                            "failed to read slice coverage finding");
            return NULL;
        }
        PyObject* item =
            build_transition_slice_coverage_finding_dict(&finding);
        if (!item) {
            Py_DECREF(coverage_findings); Py_DECREF(components);
            Py_DECREF(findings);
            alea_transition_slice_result_destroy(result);
            return NULL;
        }
        PyList_SET_ITEM(coverage_findings, (Py_ssize_t)i, item);
    }
    size_t coverage_component_count =
        alea_transition_slice_coverage_component_count(result);
    if (coverage_component_count > (size_t)PY_SSIZE_T_MAX) {
        Py_DECREF(coverage_findings); Py_DECREF(components);
        Py_DECREF(findings);
        alea_transition_slice_result_destroy(result);
        return PyErr_NoMemory();
    }
    PyObject* coverage_components =
        PyList_New((Py_ssize_t)coverage_component_count);
    if (!coverage_components) {
        Py_DECREF(coverage_findings); Py_DECREF(components);
        Py_DECREF(findings);
        alea_transition_slice_result_destroy(result);
        return NULL;
    }
    for (size_t i = 0; i < coverage_component_count; i++) {
        alea_transition_slice_coverage_component_t component;
        if (alea_transition_slice_coverage_component_get(
                result, i, &component) != 0) {
            Py_DECREF(coverage_components); Py_DECREF(coverage_findings);
            Py_DECREF(components); Py_DECREF(findings);
            alea_transition_slice_result_destroy(result);
            PyErr_SetString(PyExc_RuntimeError,
                            "failed to read slice coverage component");
            return NULL;
        }
        PyObject* item = Py_BuildValue(
            "{s:s,s:O,s:s,s:n,s:n,s:n,s:I,s:N,s:N,s:N,s:N}",
            "kind", point_coverage_kind_name(component.kind),
            "truncated", component.truncated ? Py_True : Py_False,
            "orientation",
                component.orientation == ALEA_TRANSITION_SLICE_HORIZONTAL
                    ? "horizontal" : "vertical",
            "first_finding_index",
                (Py_ssize_t)component.first_finding_index,
            "finding_count", (Py_ssize_t)component.finding_count,
            "owner_count_lower_bound",
                (Py_ssize_t)component.owner_count_lower_bound,
            "max_refinement_depth", component.max_refinement_depth,
            "uv_min", Py_BuildValue("(dd)", component.uv_min[0],
                                     component.uv_min[1]),
            "uv_max", Py_BuildValue("(dd)", component.uv_max[0],
                                     component.uv_max[1]),
            "world_min", geom_vec3(component.world_min),
            "world_max", geom_vec3(component.world_max));
        if (!item) {
            Py_DECREF(coverage_components); Py_DECREF(coverage_findings);
            Py_DECREF(components); Py_DECREF(findings);
            alea_transition_slice_result_destroy(result);
            return NULL;
        }
        PyList_SET_ITEM(coverage_components, (Py_ssize_t)i, item);
    }
    size_t component_link_count =
        alea_transition_slice_component_link_count(result);
    if (component_link_count > (size_t)PY_SSIZE_T_MAX) {
        Py_DECREF(coverage_components); Py_DECREF(coverage_findings);
        Py_DECREF(components); Py_DECREF(findings);
        alea_transition_slice_result_destroy(result);
        return PyErr_NoMemory();
    }
    PyObject* component_links =
        PyList_New((Py_ssize_t)component_link_count);
    if (!component_links) {
        Py_DECREF(coverage_components); Py_DECREF(coverage_findings);
        Py_DECREF(components); Py_DECREF(findings);
        alea_transition_slice_result_destroy(result);
        return NULL;
    }
    for (size_t i = 0; i < component_link_count; i++) {
        alea_transition_slice_component_link_t link;
        if (alea_transition_slice_component_link_get(result, i, &link) != 0) {
            Py_DECREF(component_links); Py_DECREF(coverage_components);
            Py_DECREF(coverage_findings); Py_DECREF(components);
            Py_DECREF(findings);
            alea_transition_slice_result_destroy(result);
            PyErr_SetString(PyExc_RuntimeError,
                            "failed to read slice component link");
            return NULL;
        }
        PyObject* item = Py_BuildValue(
            "{s:n,s:n,s:I,s:n}",
            "transition_component_index",
                (Py_ssize_t)link.transition_component_index,
            "coverage_component_index",
                (Py_ssize_t)link.coverage_component_index,
            "boundary_side_flags", link.boundary_sides,
            "witness_pair_count", (Py_ssize_t)link.witness_pair_count);
        if (!item) {
            Py_DECREF(component_links); Py_DECREF(coverage_components);
            Py_DECREF(coverage_findings); Py_DECREF(components);
            Py_DECREF(findings);
            alea_transition_slice_result_destroy(result);
            return NULL;
        }
        PyList_SET_ITEM(component_links, (Py_ssize_t)i, item);
    }
    alea_transition_slice_stats_t stats;
    if (alea_transition_slice_stats(result, &stats) != 0) {
        Py_DECREF(component_links);
        Py_DECREF(coverage_components);
        Py_DECREF(coverage_findings);
        Py_DECREF(components);
        Py_DECREF(findings);
        alea_transition_slice_result_destroy(result);
        PyErr_SetString(PyExc_RuntimeError, "failed to read transition slice stats");
        return NULL;
    }
    PyObject* out = Py_BuildValue(
        "{s:N,s:O,s:O,s:s,s:K,s:K,s:K,s:K,s:K,s:K,s:K,s:K,s:K,s:K,s:K,s:K,s:K}",
        "findings", findings,
        "complete", stats.complete ? Py_True : Py_False,
        "converged", stats.converged ? Py_True : Py_False,
        "stop_reason", alea_transition_slice_stop_reason_name(stats.stop_reason),
        "requested_rays", (unsigned long long)stats.requested_rays,
        "executed_rays", (unsigned long long)stats.executed_rays,
        "horizontal_rays_executed", (unsigned long long)stats.horizontal_rays_executed,
        "vertical_rays_executed", (unsigned long long)stats.vertical_rays_executed,
        "events_checked", (unsigned long long)stats.events_checked,
        "physical_events_seen", (unsigned long long)stats.physical_events_seen,
        "valid_transitions", (unsigned long long)stats.valid_transitions,
        "finding_count", (unsigned long long)stats.findings,
        "coverage_fallbacks", (unsigned long long)stats.coverage_fallbacks,
        "skipped_void_transitions", (unsigned long long)stats.skipped_void_transitions,
        "peak_live_events", (unsigned long long)stats.peak_live_events,
        "peak_live_event_bytes", (unsigned long long)stats.peak_live_event_bytes,
        "retained_output_bytes", (unsigned long long)stats.retained_output_bytes);
    if (!out) {
        Py_DECREF(components);
        Py_DECREF(coverage_findings);
        Py_DECREF(coverage_components);
        Py_DECREF(component_links);
    } else if (dict_set_new(out, "components", components) < 0) {
        Py_DECREF(coverage_findings);
        Py_DECREF(coverage_components);
        Py_DECREF(component_links);
        Py_DECREF(out);
        out = NULL;
    } else if (dict_set_new(
                   out, "coverage_findings", coverage_findings) < 0) {
        Py_DECREF(coverage_components);
        Py_DECREF(component_links);
        Py_DECREF(out);
        out = NULL;
    } else if (dict_set_new(
                   out, "coverage_components", coverage_components) < 0) {
        Py_DECREF(component_links);
        Py_DECREF(out);
        out = NULL;
    } else if (dict_set_new(
                   out, "component_links", component_links) < 0) {
        Py_DECREF(out);
        out = NULL;
    }
    if (out && (
        dict_set_new(out, "refinement_status", PyUnicode_FromString(
            alea_transition_slice_refinement_status_name(
                stats.refinement_status))) < 0 ||
        dict_set_new(out, "refined_rays_executed",
                     PyLong_FromSize_t(stats.refined_rays_executed)) < 0 ||
        dict_set_new(out, "component_count",
                     PyLong_FromSize_t(stats.components)) < 0 ||
        dict_set_new(out, "coverage_probe_count",
                     PyLong_FromSize_t(stats.coverage_probes)) < 0 ||
        dict_set_new(out, "unique_coverage_probe_count",
                     PyLong_FromSize_t(stats.unique_coverage_probes)) < 0 ||
        dict_set_new(out, "coverage_finding_count",
                     PyLong_FromSize_t(stats.coverage_findings)) < 0 ||
        dict_set_new(out, "coverage_component_count",
                     PyLong_FromSize_t(stats.coverage_components)) < 0 ||
        dict_set_new(out, "component_link_count",
                     PyLong_FromSize_t(stats.component_links)) < 0 ||
        dict_set_new(out, "critical_enabled",
                     PyBool_FromLong(stats.critical_enabled)) < 0 ||
        dict_set_new(out, "critical_complete",
                     PyBool_FromLong(stats.critical_complete)) < 0 ||
        dict_set_new(out, "critical_scope", PyUnicode_FromString(
                     stats.critical_enabled ? "suspicious_tiles" :
                                              "disabled")) < 0 ||
        dict_set_new(out, "critical_stage", PyUnicode_FromString(
                     stats.critical_enabled ? "sector_coverage" :
                                              "disabled")) < 0 ||
        dict_set_new(out, "critical_stop_reason", PyUnicode_FromString(
            alea_transition_slice_critical_stop_reason_name(
                stats.critical_stop_reason))) < 0 ||
        dict_set_new(out, "refinement_frontier_count",
                     PyLong_FromSize_t(stats.refinement_frontiers)) < 0 ||
        dict_set_new(out, "omitted_refinement_frontier_count",
                     PyLong_FromSize_t(
                         stats.omitted_refinement_frontiers)) < 0 ||
        dict_set_new(out, "critical_tile_seed_count",
                     PyLong_FromSize_t(stats.critical_tile_seeds)) < 0 ||
        dict_set_new(out, "critical_tile_count",
                     PyLong_FromSize_t(stats.critical_tiles)) < 0 ||
        dict_set_new(out, "critical_tile_source_count",
                     PyLong_FromSize_t(stats.critical_tile_sources)) < 0 ||
        dict_set_new(out, "omitted_critical_tile_source_count",
                     PyLong_FromSize_t(
                         stats.omitted_critical_tile_sources)) < 0 ||
        dict_set_new(out, "peak_critical_scratch_bytes",
                     PyLong_FromSize_t(
                         stats.peak_critical_scratch_bytes)) < 0 ||
        dict_set_new(out, "critical_tiles_processed",
                     PyLong_FromSize_t(stats.critical_tiles_processed)) < 0 ||
        dict_set_new(out, "critical_tiles_saturated",
                     PyLong_FromSize_t(stats.critical_tiles_saturated)) < 0 ||
        dict_set_new(out, "critical_region_hits",
                     PyLong_FromSize_t(stats.critical_region_hits)) < 0 ||
        dict_set_new(out, "critical_region_candidates_scanned",
                     PyLong_FromSize_t(
                         stats.critical_region_candidates_scanned)) < 0 ||
        dict_set_new(out, "critical_occurrence_seed_points",
                     PyLong_FromSize_t(
                         stats.critical_occurrence_seed_points)) < 0 ||
        dict_set_new(out, "critical_occurrence_paths",
                     PyLong_FromSize_t(
                         stats.critical_occurrence_paths)) < 0 ||
        dict_set_new(out, "critical_occurrence_universe_queries",
                     PyLong_FromSize_t(
                         stats.critical_occurrence_universe_queries)) < 0 ||
        dict_set_new(out, "critical_exhaustive_chain_hits",
                     PyLong_FromSize_t(
                         stats.critical_exhaustive_chain_hits)) < 0 ||
        dict_set_new(out, "occurrence_enumeration_mode",
                     PyUnicode_FromString(
                         stats.occurrence_discovery ==
                             ALEA_TRANSITION_SLICE_OCCURRENCE_EXHAUSTIVE
                         ? "exhaustive_region_chain_v1"
                         : "sampled_seed_paths")) < 0 ||
        dict_set_new(out, "occurrence_enumeration_complete",
                     PyBool_FromLong(
                         stats.occurrence_enumeration_complete)) < 0 ||
        dict_set_new(out, "critical_root_region_fallbacks",
                     PyLong_FromSize_t(
                         stats.critical_root_region_fallbacks)) < 0 ||
        dict_set_new(out, "critical_chain_truncated_hits",
                     PyLong_FromSize_t(
                         stats.critical_chain_truncated_hits)) < 0 ||
        dict_set_new(out, "critical_surface_references",
                     PyLong_FromSize_t(
                         stats.critical_surface_references)) < 0 ||
        dict_set_new(out, "critical_duplicate_surface_occurrences",
                     PyLong_FromSize_t(
                         stats.critical_duplicate_surface_occurrences)) < 0 ||
        dict_set_new(out, "critical_curve_count",
                     PyLong_FromSize_t(stats.critical_curves)) < 0 ||
        dict_set_new(out, "critical_curves_culled",
                     PyLong_FromSize_t(stats.critical_curves_culled)) < 0 ||
        dict_set_new(out, "critical_ranked_curves_omitted",
                     PyLong_FromSize_t(
                         stats.critical_ranked_curves_omitted)) < 0 ||
        dict_set_new(out, "critical_active_boundary_test_count",
                     PyLong_FromSize_t(
                         stats.critical_active_boundary_tests)) < 0 ||
        dict_set_new(out, "critical_active_boundary_fallback_count",
                     PyLong_FromSize_t(
                         stats.critical_active_boundary_fallbacks)) < 0 ||
        dict_set_new(out, "critical_active_capacity_fallback_count",
                     PyLong_FromSize_t(
                         stats.critical_active_capacity_fallbacks)) < 0 ||
        dict_set_new(out, "critical_active_test_budget_fallback_count",
                     PyLong_FromSize_t(
                         stats.critical_active_test_budget_fallbacks)) < 0 ||
        dict_set_new(out, "critical_active_unsupported_parabola_fallback_count",
                     PyLong_FromSize_t(stats.
                         critical_active_unsupported_parabola_fallbacks)) < 0 ||
        dict_set_new(out, "critical_active_unsupported_hyperbola_fallback_count",
                     PyLong_FromSize_t(stats.
                         critical_active_unsupported_hyperbola_fallbacks)) < 0 ||
        dict_set_new(out, "critical_active_unsupported_quartic_fallback_count",
                     PyLong_FromSize_t(stats.
                         critical_active_unsupported_quartic_fallbacks)) < 0 ||
        dict_set_new(out, "critical_active_unsupported_polygon_fallback_count",
                     PyLong_FromSize_t(stats.
                         critical_active_unsupported_polygon_fallbacks)) < 0 ||
        dict_set_new(out, "critical_active_unsupported_point_fallback_count",
                     PyLong_FromSize_t(stats.
                         critical_active_unsupported_point_fallbacks)) < 0 ||
        dict_set_new(out, "critical_active_unsupported_other_fallback_count",
                     PyLong_FromSize_t(stats.
                         critical_active_unsupported_other_fallbacks)) < 0 ||
        dict_set_new(out, "critical_active_evaluation_line_fallback_count",
                     PyLong_FromSize_t(stats.
                         critical_active_evaluation_line_fallbacks)) < 0 ||
        dict_set_new(out,
                     "critical_active_evaluation_closed_conic_fallback_count",
                     PyLong_FromSize_t(stats.
                         critical_active_evaluation_closed_conic_fallbacks)) < 0 ||
        dict_set_new(out,
                     "critical_active_evaluation_general_conic_fallback_count",
                     PyLong_FromSize_t(stats.
                         critical_active_evaluation_general_conic_fallbacks)) < 0 ||
        dict_set_new(out, "critical_active_evaluation_quartic_fallback_count",
                     PyLong_FromSize_t(stats.
                         critical_active_evaluation_quartic_fallbacks)) < 0 ||
        dict_set_new(out, "critical_active_evaluation_polygon_fallback_count",
                     PyLong_FromSize_t(stats.
                         critical_active_evaluation_polygon_fallbacks)) < 0 ||
        dict_set_new(out, "critical_active_evaluation_other_fallback_count",
                     PyLong_FromSize_t(stats.
                         critical_active_evaluation_other_fallbacks)) < 0 ||
        dict_set_new(out,
                     "critical_active_open_conic_canonical_fallback_count",
                     PyLong_FromSize_t(stats.
                         critical_active_open_conic_canonical_fallbacks)) < 0 ||
        dict_set_new(out,
                     "critical_active_open_conic_breakpoint_fallback_count",
                     PyLong_FromSize_t(stats.
                         critical_active_open_conic_breakpoint_fallbacks)) < 0 ||
        dict_set_new(out, "critical_active_segment_count",
                     PyLong_FromSize_t(stats.critical_active_segments)) < 0 ||
        dict_set_new(out, "critical_whole_curve_fallback_count",
                     PyLong_FromSize_t(
                         stats.critical_whole_curve_fallbacks)) < 0 ||
        dict_set_new(out, "peak_critical_curves",
                     PyLong_FromSize_t(stats.peak_critical_curves)) < 0 ||
        dict_set_new(out, "critical_point_candidate_count",
                     PyLong_FromSize_t(stats.critical_point_candidates)) < 0 ||
        dict_set_new(out, "critical_point_count",
                     PyLong_FromSize_t(stats.critical_points)) < 0 ||
        dict_set_new(out, "critical_duplicate_point_count",
                     PyLong_FromSize_t(stats.critical_duplicate_points)) < 0 ||
        dict_set_new(out, "critical_unsupported_curve_count",
                     PyLong_FromSize_t(stats.critical_unsupported_curves)) < 0 ||
        dict_set_new(out, "critical_unsupported_parabola_curve_count",
                     PyLong_FromSize_t(
                         stats.critical_unsupported_parabola_curves)) < 0 ||
        dict_set_new(out, "critical_unsupported_hyperbola_curve_count",
                     PyLong_FromSize_t(
                         stats.critical_unsupported_hyperbola_curves)) < 0 ||
        dict_set_new(out, "critical_unsupported_quartic_curve_count",
                     PyLong_FromSize_t(
                         stats.critical_unsupported_quartic_curves)) < 0 ||
        dict_set_new(out, "critical_unsupported_other_curve_count",
                     PyLong_FromSize_t(
                         stats.critical_unsupported_other_curves)) < 0 ||
        dict_set_new(out, "critical_probe_count",
                     PyLong_FromSize_t(stats.critical_probes)) < 0 ||
        dict_set_new(out, "critical_probe_event_count",
                     PyLong_FromSize_t(stats.critical_probe_events)) < 0 ||
        dict_set_new(out, "critical_probe_finding_count",
                     PyLong_FromSize_t(stats.critical_probe_findings)) < 0 ||
        dict_set_new(out, "critical_finding_count",
                     PyLong_FromSize_t(stats.critical_findings)) < 0 ||
        dict_set_new(out, "omitted_critical_finding_count",
                     PyLong_FromSize_t(stats.omitted_critical_findings)) < 0 ||
        dict_set_new(out, "critical_boundary_evidence_count",
                     PyLong_FromSize_t(
                         stats.critical_boundary_evidence)) < 0 ||
        dict_set_new(out, "omitted_critical_boundary_evidence_count",
                     PyLong_FromSize_t(
                         stats.omitted_critical_boundary_evidence)) < 0 ||
        dict_set_new(out, "critical_curve_pair_candidate_count",
                     PyLong_FromSize_t(
                         stats.critical_curve_pair_candidates)) < 0 ||
        dict_set_new(out, "critical_curve_pair_tested_count",
                     PyLong_FromSize_t(
                         stats.critical_curve_pairs_tested)) < 0 ||
        dict_set_new(out, "critical_unsupported_curve_pair_count",
                     PyLong_FromSize_t(
                         stats.critical_unsupported_curve_pairs)) < 0 ||
        dict_set_new(out, "critical_unsupported_general_conic_pair_count",
                     PyLong_FromSize_t(
                         stats.critical_unsupported_general_conic_pairs)) < 0 ||
        dict_set_new(out, "critical_unsupported_quartic_pair_count",
                     PyLong_FromSize_t(
                         stats.critical_unsupported_quartic_pairs)) < 0 ||
        dict_set_new(out, "critical_unsupported_quartic_line_pair_count",
                     PyLong_FromSize_t(
                         stats.critical_unsupported_quartic_line_pairs)) < 0 ||
        dict_set_new(out, "critical_unsupported_quartic_closed_conic_pair_count",
                     PyLong_FromSize_t(stats.
                         critical_unsupported_quartic_closed_conic_pairs)) < 0 ||
        dict_set_new(out, "critical_unsupported_quartic_general_conic_pair_count",
                     PyLong_FromSize_t(stats.
                         critical_unsupported_quartic_general_conic_pairs)) < 0 ||
        dict_set_new(out, "critical_unsupported_quartic_quartic_pair_count",
                     PyLong_FromSize_t(stats.
                         critical_unsupported_quartic_quartic_pairs)) < 0 ||
        dict_set_new(out, "critical_unsupported_quartic_other_pair_count",
                     PyLong_FromSize_t(
                         stats.critical_unsupported_quartic_other_pairs)) < 0 ||
        dict_set_new(out, "critical_unsupported_polygon_pair_count",
                     PyLong_FromSize_t(
                         stats.critical_unsupported_polygon_pairs)) < 0 ||
        dict_set_new(out, "critical_unsupported_other_pair_count",
                     PyLong_FromSize_t(
                         stats.critical_unsupported_other_pairs)) < 0 ||
        dict_set_new(out, "critical_pair_intersection_point_count",
                     PyLong_FromSize_t(
                         stats.critical_pair_intersection_points)) < 0 ||
        dict_set_new(out, "critical_pair_algebraic_point_count",
                     PyLong_FromSize_t(
                         stats.critical_pair_algebraic_points)) < 0 ||
        dict_set_new(out, "critical_pair_domain_rejection_count",
                     PyLong_FromSize_t(
                         stats.critical_pair_domain_rejections)) < 0 ||
        dict_set_new(out, "critical_sector_witness_count",
                     PyLong_FromSize_t(stats.critical_sector_witnesses)) < 0 ||
        dict_set_new(out, "critical_sector_gap_witness_count",
                     PyLong_FromSize_t(
                         stats.critical_sector_gap_witnesses)) < 0 ||
        dict_set_new(out, "critical_sector_overlap_witness_count",
                     PyLong_FromSize_t(
                         stats.critical_sector_overlap_witnesses)) < 0 ||
        dict_set_new(out, "critical_sector_unresolved_witness_count",
                     PyLong_FromSize_t(
                         stats.critical_sector_unresolved_witnesses)) < 0 ||
        dict_set_new(out, "truncated_coverage_probe_count",
                     PyLong_FromSize_t(stats.truncated_coverage_probes)) < 0 ||
        dict_set_new(out, "skipped_unowned_coverage_probe_count",
                     PyLong_FromSize_t(
                         stats.skipped_unowned_coverage_probes)) < 0 ||
        dict_set_new(out, "max_refinement_depth_reached",
                     PyLong_FromUnsignedLong(
                         stats.max_refinement_depth_reached)) < 0 ||
        dict_set_new(out, "peak_row_scratch_bytes",
                     PyLong_FromSize_t(stats.peak_row_scratch_bytes)) < 0 ||
        dict_set_new(out, "peak_scratch_bytes",
                     PyLong_FromSize_t(stats.peak_scratch_bytes)) < 0)) {
        Py_DECREF(out);
        out = NULL;
    }
    if (out) {
        size_t frontier_count =
            alea_transition_slice_refinement_frontier_count(result);
        PyObject* frontiers = PyList_New((Py_ssize_t)frontier_count);
        if (!frontiers) {
            Py_DECREF(out);
            out = NULL;
        }
        for (size_t i = 0; out && i < frontier_count; i++) {
            alea_transition_slice_refinement_frontier_t frontier;
            if (alea_transition_slice_refinement_frontier_get(
                    result, i, &frontier) != 0) {
                Py_DECREF(frontiers); Py_DECREF(out); out = NULL;
                PyErr_SetString(PyExc_RuntimeError,
                                "failed to read refinement frontier");
                break;
            }
            PyObject* item = Py_BuildValue(
                "{s:s,s:I,s:d,s:d,s:K,s:K,s:K,s:K,s:n,s:N,s:N}",
                "orientation",
                    frontier.orientation == ALEA_TRANSITION_SLICE_HORIZONTAL
                        ? "horizontal" : "vertical",
                "refinement_depth", frontier.refinement_depth,
                "transverse_min", frontier.transverse_min,
                "transverse_max", frontier.transverse_max,
                "signature_a_min",
                    (unsigned long long)frontier.signature_a[0],
                "signature_a_max",
                    (unsigned long long)frontier.signature_a[1],
                "signature_b_min",
                    (unsigned long long)frontier.signature_b[0],
                "signature_b_max",
                    (unsigned long long)frontier.signature_b[1],
                "max_event_count", (Py_ssize_t)frontier.max_event_count,
                "uv_min", Py_BuildValue("(dd)", frontier.uv_min[0],
                                          frontier.uv_min[1]),
                "uv_max", Py_BuildValue("(dd)", frontier.uv_max[0],
                                          frontier.uv_max[1]));
            if (!item) {
                Py_DECREF(frontiers); Py_DECREF(out); out = NULL;
                break;
            }
            PyList_SET_ITEM(frontiers, (Py_ssize_t)i, item);
        }
        if (out && dict_set_new(out, "refinement_frontiers", frontiers) < 0) {
            Py_DECREF(out);
            out = NULL;
        }
    }
    if (out) {
        size_t tile_count = alea_transition_slice_critical_tile_count(result);
        PyObject* tiles = PyList_New((Py_ssize_t)tile_count);
        if (!tiles) {
            Py_DECREF(out);
            out = NULL;
        }
        for (size_t i = 0; out && i < tile_count; i++) {
            alea_transition_slice_critical_tile_t tile;
            if (alea_transition_slice_critical_tile_get(
                    result, i, &tile) != 0) {
                Py_DECREF(tiles); Py_DECREF(out); out = NULL;
                PyErr_SetString(PyExc_RuntimeError,
                                "failed to read critical tile");
                break;
            }
            PyObject* sources = PyList_New((Py_ssize_t)tile.source_count);
            if (!sources) {
                Py_DECREF(tiles); Py_DECREF(out); out = NULL;
                break;
            }
            for (size_t j = 0; j < tile.source_count; j++) {
                alea_transition_slice_critical_tile_source_t source;
                if (alea_transition_slice_critical_tile_source_get(
                        result, tile.first_source_index + j, &source) != 0) {
                    Py_DECREF(sources); Py_DECREF(tiles); Py_DECREF(out);
                    out = NULL;
                    PyErr_SetString(PyExc_RuntimeError,
                                    "failed to read critical tile source");
                    break;
                }
                PyObject* source_item = Py_BuildValue(
                    "{s:s,s:n}", "kind",
                    transition_slice_tile_source_kind_name(source.kind),
                    "index", (Py_ssize_t)source.source_index);
                if (!source_item) {
                    Py_DECREF(sources); Py_DECREF(tiles); Py_DECREF(out);
                    out = NULL;
                    break;
                }
                PyList_SET_ITEM(sources, (Py_ssize_t)j, source_item);
            }
            if (!out) break;
            PyObject* item = Py_BuildValue(
                "{s:N,s:I,s:N,s:N}",
                "sources", sources,
                "source_flags", tile.source_flags,
                "uv_min", Py_BuildValue("(dd)", tile.uv_min[0],
                                          tile.uv_min[1]),
                "uv_max", Py_BuildValue("(dd)", tile.uv_max[0],
                                          tile.uv_max[1]));
            if (!item) {
                Py_DECREF(tiles); Py_DECREF(out); out = NULL;
                break;
            }
            PyList_SET_ITEM(tiles, (Py_ssize_t)i, item);
        }
        if (out && dict_set_new(out, "critical_tiles", tiles) < 0) {
            Py_DECREF(out);
            out = NULL;
        }
    }
    if (out) {
        const size_t critical_count =
            alea_transition_slice_critical_finding_count(result);
        PyObject* critical_findings =
            PyList_New((Py_ssize_t)critical_count);
        if (!critical_findings) {
            Py_DECREF(out);
            out = NULL;
        }
        for (size_t i = 0; out && i < critical_count; i++) {
            alea_transition_slice_critical_finding_t finding;
            if (alea_transition_slice_critical_finding_get(
                    result, i, &finding) != 0) {
                Py_DECREF(critical_findings); Py_DECREF(out); out = NULL;
                PyErr_SetString(PyExc_RuntimeError,
                                "failed to read critical finding");
                break;
            }
            PyObject* item = build_transition_dict(&finding.transition);
            PyObject* boundary_pieces = PyList_New(
                (Py_ssize_t)finding.boundary_piece_count);
            if (!boundary_pieces) {
                Py_XDECREF(item);
                Py_DECREF(critical_findings); Py_DECREF(out); out = NULL;
                break;
            }
            for (size_t piece_index = 0;
                 piece_index < finding.boundary_piece_count;
                 piece_index++) {
                const alea_transition_slice_boundary_piece_t* piece =
                    &finding.boundary_pieces[piece_index];
                PyObject* points = PyList_New((Py_ssize_t)piece->point_count);
                if (!points) {
                    Py_DECREF(boundary_pieces); Py_XDECREF(item);
                    Py_DECREF(critical_findings); Py_DECREF(out); out = NULL;
                    break;
                }
                for (size_t point_index = 0;
                     point_index < piece->point_count; point_index++) {
                    PyObject* point = Py_BuildValue(
                        "(dd)", piece->uv[point_index][0],
                                piece->uv[point_index][1]);
                    if (!point) {
                        Py_DECREF(points); Py_DECREF(boundary_pieces);
                        Py_XDECREF(item); Py_DECREF(critical_findings);
                        Py_DECREF(out); out = NULL;
                        break;
                    }
                    PyList_SET_ITEM(points, (Py_ssize_t)point_index, point);
                }
                if (!out) break;
                PyObject* piece_item = Py_BuildValue(
                    "{s:i,s:I,s:N}",
                    "surface_id", piece->surface_id,
                    "role_flags", piece->role_flags,
                    "uv_points", points);
                if (!piece_item) {
                    Py_DECREF(boundary_pieces); Py_XDECREF(item);
                    Py_DECREF(critical_findings); Py_DECREF(out); out = NULL;
                    break;
                }
                PyList_SET_ITEM(boundary_pieces,
                                (Py_ssize_t)piece_index, piece_item);
            }
            if (!out) break;
            if (!item ||
                dict_set_new(item, "tile_index",
                             PyLong_FromSize_t(finding.tile_index)) < 0 ||
                dict_set_new(item, "point_index",
                             PyLong_FromSize_t(finding.point_index)) < 0 ||
                dict_set_new(item, "source_cell_id",
                             PyLong_FromLong(finding.source_cell_id)) < 0 ||
                dict_set_new(item, "source_surface_id",
                             PyLong_FromLong(finding.source_surface_id)) < 0 ||
                dict_set_new(item, "source_occurrence_key",
                             PyLong_FromUnsignedLongLong(
                                 finding.source_occurrence_key)) < 0 ||
                dict_set_new(item, "source_universe_occurrence_key",
                             PyLong_FromUnsignedLongLong(
                                 finding.source_universe_occurrence_key)) < 0 ||
                dict_set_new(item, "uv", Py_BuildValue(
                    "(dd)", finding.uv[0], finding.uv[1])) < 0 ||
                dict_set_new(item, "world_point",
                             geom_vec3(finding.world_point)) < 0 ||
                dict_set_new(item, "probe_direction",
                             geom_vec3(finding.direction)) < 0 ||
                dict_set_new(item, "probe_radius",
                             PyFloat_FromDouble(finding.radius)) < 0 ||
                dict_set_new(item, "boundary_evidence_truncated",
                             PyBool_FromLong(
                                 finding.boundary_evidence_truncated)) < 0) {
                Py_XDECREF(item);
                Py_DECREF(boundary_pieces);
                Py_DECREF(critical_findings); Py_DECREF(out); out = NULL;
                break;
            }
            if (dict_set_new(item, "boundary_pieces", boundary_pieces) < 0) {
                Py_DECREF(item);
                Py_DECREF(critical_findings); Py_DECREF(out); out = NULL;
                break;
            }
            PyList_SET_ITEM(critical_findings, (Py_ssize_t)i, item);
        }
        if (out && dict_set_new(
                out, "critical_findings", critical_findings) < 0) {
            Py_DECREF(out);
            out = NULL;
        }
    }
    alea_transition_slice_result_destroy(result);
    return out;
}

static PyObject* PyAleaSystem_transition_slice_screen(
    PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    PyObject *origin_obj, *normal_obj, *up_obj, *opts = NULL;
    double u_min, u_max, v_min, v_max;
    static char* kwlist[] = {
        "origin", "normal", "up", "u_min", "u_max", "v_min", "v_max",
        "options", NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOdddd|O", kwlist,
            &origin_obj, &normal_obj, &up_obj,
            &u_min, &u_max, &v_min, &v_max, &opts))
        return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    double origin[3], normal[3], up[3];
    if (transition_parse_vec3(origin_obj, "origin", origin) < 0 ||
        transition_parse_vec3(normal_obj, "normal", normal) < 0 ||
        transition_parse_vec3(up_obj, "up", up) < 0)
        return NULL;
    alea_transition_slice_options_t options;
    if (parse_transition_slice_options(opts, &options) < 0) return NULL;
    if (ensure_query_acceleration(self) < 0) return NULL;

    alea_slice_view_t view;
    alea_slice_view_init(&view, origin[0], origin[1], origin[2],
                         normal[0], normal[1], normal[2],
                         up[0], up[1], up[2],
                         u_min, u_max, v_min, v_max);
    alea_transition_slice_result_t* result =
        alea_transition_slice_result_create();
    if (!result) return PyErr_NoMemory();
    int rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    rc = alea_transition_slice_screen(self->sys, &view, &options, result);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) {
        alea_transition_slice_result_destroy(result);
        return NULL;
    }
    if (rc != 0) {
        alea_transition_slice_result_destroy(result);
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }
    return transition_slice_result_to_py(result);
}

static PyObject* PyAleaSystem_transition_slice_screen_batch(
    PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    PyObject *origin_obj, *normal_obj, *up_obj, *bounds_obj, *opts = NULL;
    Py_ssize_t requested_workers = 0;
    unsigned long long max_parallel_scratch_bytes = 0;
    static char* kwlist[] = {
        "origin", "normal", "up", "bounds", "options",
        "requested_workers", "max_parallel_scratch_bytes", NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOO|OnK", kwlist,
            &origin_obj, &normal_obj, &up_obj, &bounds_obj, &opts,
            &requested_workers, &max_parallel_scratch_bytes))
        return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    if (requested_workers < 0) {
        PyErr_SetString(PyExc_ValueError,
                        "requested_workers must be non-negative");
        return NULL;
    }
    double origin[3], normal[3], up[3];
    if (transition_parse_vec3(origin_obj, "origin", origin) < 0 ||
        transition_parse_vec3(normal_obj, "normal", normal) < 0 ||
        transition_parse_vec3(up_obj, "up", up) < 0)
        return NULL;
    alea_transition_slice_options_t options;
    if (parse_transition_slice_options(opts, &options) < 0) return NULL;
    if (ensure_query_acceleration(self) < 0) return NULL;

    PyObject* bounds = PySequence_Fast(
        bounds_obj, "bounds must be a sequence of four-number sequences");
    if (!bounds) return NULL;
    Py_ssize_t page_count_py = PySequence_Fast_GET_SIZE(bounds);
    size_t page_count = (size_t)page_count_py;
    alea_slice_view_t* views = page_count
        ? PyMem_Calloc(page_count, sizeof(*views)) : NULL;
    alea_transition_slice_result_t** results = page_count
        ? PyMem_Calloc(page_count, sizeof(*results)) : NULL;
    if (page_count && (!views || !results)) {
        Py_DECREF(bounds); PyMem_Free(views); PyMem_Free(results);
        return PyErr_NoMemory();
    }
    int parse_failed = 0;
    for (size_t page = 0; page < page_count; page++) {
        PyObject* item = PySequence_Fast(
            PySequence_Fast_GET_ITEM(bounds, (Py_ssize_t)page),
            "each page bound must contain four numbers");
        if (!item) { parse_failed = 1; break; }
        if (PySequence_Fast_GET_SIZE(item) != 4) {
            Py_DECREF(item);
            PyErr_SetString(PyExc_ValueError,
                            "each page bound must contain four numbers");
            parse_failed = 1;
            break;
        }
        double values[4];
        for (size_t i = 0; i < 4; i++) {
            values[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(
                item, (Py_ssize_t)i));
            if (PyErr_Occurred()) { parse_failed = 1; break; }
        }
        Py_DECREF(item);
        if (parse_failed) break;
        alea_slice_view_init(&views[page], origin[0], origin[1], origin[2],
                             normal[0], normal[1], normal[2],
                             up[0], up[1], up[2],
                             values[0], values[1], values[2], values[3]);
        results[page] = alea_transition_slice_result_create();
        if (!results[page]) { PyErr_NoMemory(); parse_failed = 1; break; }
    }
    Py_DECREF(bounds);
    if (parse_failed) {
        for (size_t page = 0; page < page_count; page++)
            alea_transition_slice_result_destroy(results[page]);
        PyMem_Free(results); PyMem_Free(views);
        return NULL;
    }

    alea_transition_slice_batch_stats_t stats;
    int rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    rc = alea_transition_slice_screen_batch(
        self->sys, views, page_count, &options,
        (size_t)requested_workers, (uint64_t)max_parallel_scratch_bytes,
        results, &stats);
    Py_END_ALLOW_THREADS
    PyMem_Free(views);
    if (restore_sigint(old_sigint)) rc = -1;
    if (rc != 0) {
        for (size_t page = 0; page < page_count; page++)
            alea_transition_slice_result_destroy(results[page]);
        PyMem_Free(results);
        if (!PyErr_Occurred())
            PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }
    PyObject* pages = PyList_New(page_count_py);
    if (!pages) {
        for (size_t page = 0; page < page_count; page++)
            alea_transition_slice_result_destroy(results[page]);
        PyMem_Free(results);
        return NULL;
    }
    for (size_t page = 0; page < page_count; page++) {
        PyObject* converted = transition_slice_result_to_py(results[page]);
        results[page] = NULL;
        if (!converted) {
            for (size_t later = page + 1; later < page_count; later++)
                alea_transition_slice_result_destroy(results[later]);
            Py_DECREF(pages); PyMem_Free(results);
            return NULL;
        }
        PyList_SET_ITEM(pages, (Py_ssize_t)page, converted);
    }
    PyMem_Free(results);
    return Py_BuildValue(
        "{s:N,s:n,s:n,s:n,s:n,s:K,s:K}",
        "pages", pages,
        "page_count", (Py_ssize_t)stats.page_count,
        "completed_page_count", (Py_ssize_t)stats.completed_page_count,
        "requested_workers", (Py_ssize_t)stats.requested_workers,
        "actual_workers", (Py_ssize_t)stats.actual_workers,
        "reserved_scratch_bytes_per_worker",
            (unsigned long long)stats.reserved_scratch_bytes_per_worker,
        "reserved_parallel_scratch_bytes",
            (unsigned long long)stats.reserved_parallel_scratch_bytes);
}

static PyObject* PyAleaSystem_check_transition(
    PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    int universe_id, current_cell_id, surface_id;
    PyObject *point_obj, *direction_obj, *tied_obj = Py_None;
    double probe_distance = 0.0, max_probe_distance = 0.0;
    Py_ssize_t max_coverage_hits = 256;
    static char* kwlist[] = {
        "universe", "current_cell", "surface", "point", "direction",
        "tied_surfaces", "probe_distance", "max_probe_distance",
        "max_coverage_hits", NULL
    };
    if (!PyArg_ParseTupleAndKeywords(
            args, kwds, "iiiOO|Oddn", kwlist,
            &universe_id, &current_cell_id, &surface_id,
            &point_obj, &direction_obj, &tied_obj,
            &probe_distance, &max_probe_distance, &max_coverage_hits))
        return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    if (max_coverage_hits <= 0 || max_coverage_hits > 16384) {
        PyErr_SetString(PyExc_ValueError,
                        "max_coverage_hits must be between 1 and 16384");
        return NULL;
    }
    double point[3], direction[3];
    if (transition_parse_vec3(point_obj, "point", point) < 0 ||
        transition_parse_vec3(direction_obj, "direction", direction) < 0)
        return NULL;

    int* tied_ids = NULL;
    size_t tied_count = 0;
    PyObject* tied_seq = NULL;
    if (tied_obj != Py_None) {
        tied_seq = PySequence_Fast(tied_obj, "tied_surfaces must be iterable");
        if (!tied_seq) return NULL;
        Py_ssize_t count = PySequence_Fast_GET_SIZE(tied_seq);
        if (count < 0 || count > 256) {
            Py_DECREF(tied_seq);
            PyErr_SetString(PyExc_ValueError,
                            "tied_surfaces may contain at most 256 IDs");
            return NULL;
        }
        tied_count = (size_t)count;
        if (tied_count > 0) {
            tied_ids = PyMem_Malloc(tied_count * sizeof(*tied_ids));
            if (!tied_ids) { Py_DECREF(tied_seq); return PyErr_NoMemory(); }
        }
        for (size_t i = 0; i < tied_count; i++) {
            long id = PyLong_AsLong(PySequence_Fast_GET_ITEM(tied_seq, (Py_ssize_t)i));
            if ((id == -1 && PyErr_Occurred()) || id < INT_MIN || id > INT_MAX) {
                PyMem_Free(tied_ids); Py_DECREF(tied_seq);
                PyErr_SetString(PyExc_TypeError,
                                "tied_surfaces must contain integer IDs");
                return NULL;
            }
            tied_ids[i] = (int)id;
        }
        Py_DECREF(tied_seq);
    }

    if (ensure_query_acceleration(self) < 0) {
        PyMem_Free(tied_ids);
        return NULL;
    }
    alea_transition_options_t options;
    alea_transition_options_init(&options);
    if (probe_distance > 0.0) options.probe_distance = probe_distance;
    if (max_probe_distance > 0.0)
        options.max_probe_distance = max_probe_distance;
    options.max_coverage_hits = (size_t)max_coverage_hits;
    alea_transition_result_t result;
    int rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    rc = alea_check_transition_local(
        self->sys, universe_id, current_cell_id, surface_id,
        tied_ids, tied_count, point, direction, &options, &result);
    Py_END_ALLOW_THREADS
    PyMem_Free(tied_ids);
    if (restore_sigint(old_sigint)) return NULL;
    if (rc != 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    return build_transition_dict(&result);
}
