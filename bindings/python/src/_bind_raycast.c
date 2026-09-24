// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* This file is part of pyalea_binding.c — do NOT compile separately.
 * It is #included from pyalea_binding.c to keep a single translation unit.
 *
 * Contents: Raycast, cell-aware raycast, find overlaps,
 *           volume estimation, bounding sphere.
 */

/* ============================================================================
 * PyAleaSystem Methods - Raycast
 * ============================================================================ */

static int raycast_attach_path(PyObject* segment, const alea_raycast_result_t* result,
                               size_t segment_index) {
    size_t count = alea_raycast_segment_path_count(result, segment_index);
    PyObject* path = PyList_New(count);
    if (!path) return -1;
    for (size_t i = 0; i < count; ++i) {
        alea_raycast_path_entry_t entry;
        if (alea_raycast_segment_path_get(result, segment_index, i, &entry) != 0) {
            Py_DECREF(path);
            PyErr_SetString(PyExc_RuntimeError, "Failed to read ray segment path");
            return -1;
        }
        PyObject* item = Py_BuildValue(
            "{s:I,s:i,s:i,s:i,s:i,s:i,s:i,s:(ddd),s:K}",
            "cell_index", entry.cell_index,
            "cell_id", entry.cell_id,
            "material_id", entry.material_id,
            "universe_id", entry.universe_id,
            "fill_universe", entry.fill_universe,
            "depth", entry.depth,
            "is_lattice", (int)entry.is_lattice,
            "lattice_origin", entry.lattice_origin[0], entry.lattice_origin[1],
            entry.lattice_origin[2],
            "occurrence_key", (unsigned long long)entry.occurrence_key);
        if (!item) { Py_DECREF(path); return -1; }
        PyList_SET_ITEM(path, i, item);
    }
    if (PyDict_SetItemString(segment, "path", path) != 0) {
        Py_DECREF(path);
        return -1;
    }
    Py_DECREF(path);
    return 0;
}

static PyObject* PyAleaSystem_raycast(PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    double ox, oy, oz, dx, dy, dz;
    double t_max = 0.0;
    static char* kwlist[] = {"ox", "oy", "oz", "dx", "dy", "dz", "t_max", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "dddddd|d", kwlist,
                                     &ox, &oy, &oz, &dx, &dy, &dz, &t_max)) {
        return NULL;
    }
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }
    if (ensure_query_acceleration(self) < 0) return NULL;

    alea_raycast_result_t* result = alea_raycast_result_create();
    if (!result) {
        PyErr_SetString(PyExc_MemoryError, "Failed to allocate raycast result");
        return NULL;
    }

    int raycast_result;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    raycast_result = alea_raycast(self->sys, ox, oy, oz, dx, dy, dz, t_max, result);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) {
        alea_raycast_result_free(result);
        return NULL;
    }

    if (raycast_result < 0) {
        alea_raycast_result_free(result);
        PyErr_SetString(PyExc_RuntimeError, "Raycast failed");
        return NULL;
    }

    size_t count = alea_raycast_segment_count(result);
    PyObject* segments = PyList_New(count);
    for (size_t i = 0; i < count; i++) {
        double t_enter, t_exit, density;
        int cell_id, material_id, enter_surface_id, exit_surface_id;
        alea_raycast_segment_get(
            result, i, &t_enter, &t_exit, &cell_id, &material_id, &density,
            &enter_surface_id, &exit_surface_id);
        uint8_t resolution_flags;
        if (alea_raycast_segment_resolution_flags(result, i, &resolution_flags) != 0) {
            Py_DECREF(segments);
            alea_raycast_result_free(result);
            PyErr_SetString(PyExc_RuntimeError, "Failed to read ray segment flags");
            return NULL;
        }
        PyObject* item = Py_BuildValue("{s:d,s:d,s:i,s:i,s:d,s:i,s:i,s:i}",
            "t_enter", t_enter, "t_exit", t_exit,
            "cell_id", cell_id, "material_id", material_id,
            "density", density,
            "enter_surface_id", enter_surface_id,
            "exit_surface_id", exit_surface_id,
            "resolution_flags", resolution_flags);
        PyList_SET_ITEM(segments, i, item);
    }

    /* Expose the global hit list (all surface crossings, sorted, deduped) —
     * the raw breakpoints underneath the precedence-resolved segments. */
    size_t hit_count = alea_raycast_hit_count(result);
    PyObject* hits = PyList_New(hit_count);
    if (!hits) { Py_DECREF(segments); alea_raycast_result_free(result); return NULL; }
    for (size_t i = 0; i < hit_count; i++) {
        double t;
        int surface_id;
        if (alea_raycast_hit_get(result, i, &t, &surface_id) != 0) {
            Py_DECREF(hits);
            Py_DECREF(segments);
            alea_raycast_result_free(result);
            PyErr_SetString(PyExc_RuntimeError, "Failed to read ray hit");
            return NULL;
        }
        PyObject* h = Py_BuildValue("{s:d,s:i}",
            "t", t, "surface_id", surface_id);
        if (!h) { Py_DECREF(hits); Py_DECREF(segments); alea_raycast_result_free(result); return NULL; }
        PyList_SET_ITEM(hits, i, h);
    }
    alea_raycast_result_free(result);
    PyObject* out = Py_BuildValue("{s:N,s:N}", "segments", segments,
                                  "hits", hits);
    return out;
}

/* Interval defect classification along a ray (oracle-grade error detection).
 * Owner sets are complete: catches overlaps in isolation and gaps that the
 * precedence-resolved segments hide. */
static PyObject* PyAleaSystem_classify_ray_intervals(PyAleaSystemObject* self, PyObject* args) {
    double ox, oy, oz, dx, dy, dz, t_max;
    if (!PyArg_ParseTuple(args, "ddddddd", &ox, &oy, &oz, &dx, &dy, &dz,
                          &t_max)) return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    if (ensure_query_acceleration(self) < 0) return NULL;

    enum { STACK_CAP = 512 };
    alea_ray_interval_finding_t stack_buf[STACK_CAP];
    alea_ray_interval_finding_t* buf = stack_buf;
    int n;
    Py_BEGIN_ALLOW_THREADS
    n = alea_ray_classify_intervals(self->sys, ox, oy, oz, dx, dy, dz,
                                    t_max, buf, STACK_CAP);
    Py_END_ALLOW_THREADS
    if (n < 0) {
        PyErr_SetString(PyExc_RuntimeError, "interval classification failed");
        return NULL;
    }
    if (n > STACK_CAP) {
        buf = malloc((size_t)n * sizeof(*buf));
        if (!buf) return PyErr_NoMemory();
        int n2;
        Py_BEGIN_ALLOW_THREADS
        n2 = alea_ray_classify_intervals(self->sys, ox, oy, oz, dx, dy, dz,
                                         t_max, buf, (size_t)n);
        Py_END_ALLOW_THREADS
        if (n2 < 0 || n2 > n) { free(buf);
            PyErr_SetString(PyExc_RuntimeError, "interval classification failed");
            return NULL; }
        n = n2;
    }

    PyObject* list = PyList_New(n);
    if (!list) { if (buf != stack_buf) free(buf); return NULL; }
    for (int i = 0; i < n; i++) {
        PyObject* item = Py_BuildValue("{s:d,s:d,s:i,s:i,s:i,s:i}",
            "t_enter", buf[i].t_enter, "t_exit", buf[i].t_exit,
            "kind", buf[i].kind, "cell_id", buf[i].cell_id,
            "overlap_cell_id", buf[i].overlap_cell_id,
            "depth", buf[i].depth);
        if (!item) { Py_DECREF(list); if (buf != stack_buf) free(buf); return NULL; }
        PyList_SET_ITEM(list, i, item);
    }
    if (buf != stack_buf) free(buf);
    return list;
}

static PyObject* PyAleaSystem_ray_first_cell(PyAleaSystemObject* self, PyObject* args) {
    double ox, oy, oz, dx, dy, dz, t_max = 0.0;
    if (!PyArg_ParseTuple(args, "dddddd|d", &ox, &oy, &oz, &dx, &dy, &dz, &t_max)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }
    if (ensure_query_acceleration(self) < 0) return NULL;

    double t;
    int cell_id = alea_ray_first_cell(self->sys, ox, oy, oz, dx, dy, dz, t_max, &t);
    if (cell_id < 0) Py_RETURN_NONE;
    return Py_BuildValue("(id)", cell_id, t);
}

/* First non-void interval along a ray, with optional boundary details.  This
 * is the public unified-query counterpart of the legacy ray_first_cell(). */
static PyObject* PyAleaSystem_ray_first_visible(PyAleaSystemObject* self,
                                                   PyObject* args,
                                                   PyObject* kwds) {
    double ox, oy, oz, dx, dy, dz, t_min = 0.0, t_max = 0.0;
    int material_filter = -1, include_surface = 0, include_normal = 0;
    static char* kwlist[] = {
        "ox", "oy", "oz", "dx", "dy", "dz", "t_min", "t_max",
        "material_filter", "include_surface", "include_normal", NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "dddddd|ddipp", kwlist,
                                     &ox, &oy, &oz, &dx, &dy, &dz,
                                     &t_min, &t_max, &material_filter,
                                     &include_surface, &include_normal))
        return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    if (ensure_query_acceleration(self) < 0) return NULL;

    alea_ray_first_visible_options_t options;
    alea_ray_first_visible_options_init(&options);
    options.t_min = t_min;
    options.t_max = t_max;
    options.material_filter = material_filter;
    options.fields = (include_surface ? ALEA_RAY_FIRST_VISIBLE_SURFACE_ID : 0) |
                     (include_normal ? ALEA_RAY_FIRST_VISIBLE_SURFACE_NORMAL : 0);
    alea_ray_first_visible_query_result_t* result =
        alea_ray_first_visible_query_result_create();
    if (!result) return PyErr_NoMemory();

    int rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    rc = alea_ray_first_visible_query(self->sys, ox, oy, oz, dx, dy, dz,
                                      &options, result);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) {
        alea_ray_first_visible_query_result_destroy(result);
        return NULL;
    }
    if (rc != 0) {
        const char* detail = alea_error();
        alea_ray_first_visible_query_result_destroy(result);
        PyErr_Format(PyExc_RuntimeError, "first-visible ray query failed: %s",
                     detail ? detail : "unknown libalea error");
        return NULL;
    }
    if (!alea_ray_first_visible_found(result)) {
        alea_ray_first_visible_query_result_destroy(result);
        Py_RETURN_NONE;
    }

    PyObject* out = Py_BuildValue("{s:d,s:i,s:i,s:d}",
                                  "t", alea_ray_first_visible_t(result),
                                  "cell_id", alea_ray_first_visible_cell_id(result),
                                  "material_id", alea_ray_first_visible_material_id(result),
                                  "density", alea_ray_first_visible_density(result));
    if (!out) goto failed;
    if (include_surface) {
        PyObject* surface = PyLong_FromLong(alea_ray_first_visible_surface_id(result));
        if (!surface || PyDict_SetItemString(out, "surface_id", surface) != 0) {
            Py_XDECREF(surface);
            goto failed;
        }
        Py_DECREF(surface);
    }
    if (include_normal) {
        double nx, ny, nz;
        if (alea_ray_first_visible_normal(result, &nx, &ny, &nz) != 0) {
            PyErr_SetString(PyExc_RuntimeError, "first-visible ray query did not return a normal");
            goto failed;
        }
        PyObject* normal = Py_BuildValue("(ddd)", nx, ny, nz);
        if (!normal || PyDict_SetItemString(out, "normal", normal) != 0) {
            Py_XDECREF(normal);
            goto failed;
        }
        Py_DECREF(normal);
    }
    alea_ray_first_visible_query_result_destroy(result);
    return out;

failed:
    Py_XDECREF(out);
    alea_ray_first_visible_query_result_destroy(result);
    return NULL;
}

/* Ordered canonical boundary events, including synthetic lattice boundaries
 * when the geometry traversal emits them. */
static PyObject* PyAleaSystem_ray_boundary_events(PyAleaSystemObject* self,
                                                     PyObject* args,
                                                     PyObject* kwds) {
    double ox, oy, oz, dx, dy, dz, t_min = 0.0, t_max = 0.0;
    unsigned long long max_events = 0, max_output_bytes = 0;
    int include_primitive_id = 0, include_normal = 0;
    int include_all_coincident_physical = 0, include_occurrence_provenance = 0;
    static char* kwlist[] = {
        "ox", "oy", "oz", "dx", "dy", "dz", "t_min", "t_max",
        "max_events", "max_output_bytes", "include_primitive_id",
        "include_normal", "include_all_coincident_physical",
        "include_occurrence_provenance", NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "dddddd|ddKKpppp", kwlist,
                                     &ox, &oy, &oz, &dx, &dy, &dz,
                                     &t_min, &t_max, &max_events, &max_output_bytes,
                                     &include_primitive_id, &include_normal,
                                     &include_all_coincident_physical,
                                     &include_occurrence_provenance))
        return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    if (ensure_query_acceleration(self) < 0) return NULL;

    alea_ray_boundary_event_options_t options;
    alea_ray_boundary_event_options_init(&options);
    options.t_min = t_min;
    options.t_max = t_max;
    options.max_events = (uint64_t)max_events;
    options.max_output_bytes = (uint64_t)max_output_bytes;
    options.include_all_coincident_physical = include_all_coincident_physical;
    options.include_occurrence_provenance = include_occurrence_provenance;
    options.fields = (include_primitive_id ? ALEA_RAY_BOUNDARY_EVENT_PRIMITIVE_ID : 0) |
                     (include_normal ? ALEA_RAY_BOUNDARY_EVENT_NORMAL : 0);
    alea_ray_boundary_event_query_result_t* result =
        alea_ray_boundary_event_query_result_create();
    if (!result) return PyErr_NoMemory();

    int rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    rc = alea_ray_boundary_event_query(self->sys, ox, oy, oz, dx, dy, dz,
                                       &options, result);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) {
        alea_ray_boundary_event_query_result_destroy(result);
        return NULL;
    }
    if (rc != 0) {
        const char* detail = alea_error();
        alea_ray_boundary_event_query_result_destroy(result);
        PyErr_Format(PyExc_RuntimeError, "ray boundary-event query failed: %s",
                     detail ? detail : "unknown libalea error");
        return NULL;
    }

    size_t count = alea_ray_boundary_event_count(result);
    PyObject* out = PyList_New((Py_ssize_t)count);
    if (!out) goto failed;
    for (size_t i = 0; i < count; ++i) {
        double t, nx, ny, nz;
        int kind, surface_id, cell_before, cell_after, material_before, material_after;
        uint32_t resolution_flags, primitive_id;
        if (alea_ray_boundary_event_get(result, i, &t, &kind, &surface_id,
                                        &cell_before, &cell_after, &material_before,
                                        &material_after, &resolution_flags,
                                        &primitive_id, &nx, &ny, &nz) != 0) {
            PyErr_SetString(PyExc_RuntimeError, "failed to read ray boundary event");
            goto failed;
        }
        PyObject* item = Py_BuildValue("{s:d,s:i,s:i,s:i,s:i,s:i,s:i,s:I}",
                                       "t", t, "kind", kind, "surface_id", surface_id,
                                       "cell_before", cell_before, "cell_after", cell_after,
                                       "material_before", material_before,
                                       "material_after", material_after,
                                       "resolution_flags", resolution_flags);
        if (!item) goto failed;
        if (include_primitive_id) {
            PyObject* value = PyLong_FromUnsignedLong(primitive_id);
            if (!value || PyDict_SetItemString(item, "primitive_id", value) != 0) {
                Py_XDECREF(value);
                Py_DECREF(item);
                goto failed;
            }
            Py_DECREF(value);
        }
        if (include_normal) {
            PyObject* normal = Py_BuildValue("(ddd)", nx, ny, nz);
            if (!normal || PyDict_SetItemString(item, "normal", normal) != 0) {
                Py_XDECREF(normal);
                Py_DECREF(item);
                goto failed;
            }
            Py_DECREF(normal);
        }
        if (include_occurrence_provenance) {
            alea_ray_boundary_event_provenance_t provenance;
            if (alea_ray_boundary_event_provenance_get(
                    result, i, &provenance) != 0) {
                Py_DECREF(item);
                PyErr_SetString(PyExc_RuntimeError,
                                "failed to read boundary-event occurrence provenance");
                goto failed;
            }
            PyObject* local_surfaces = PyList_New(
                (Py_ssize_t)provenance.local_surface_count);
            PyObject* local_point = Py_BuildValue(
                "(ddd)", provenance.local_point[0], provenance.local_point[1],
                provenance.local_point[2]);
            PyObject* local_direction = Py_BuildValue(
                "(ddd)", provenance.local_direction[0],
                provenance.local_direction[1], provenance.local_direction[2]);
            if (!local_surfaces || !local_point || !local_direction) {
                Py_XDECREF(local_surfaces); Py_XDECREF(local_point);
                Py_XDECREF(local_direction); Py_DECREF(item);
                goto failed;
            }
            for (size_t surface = 0;
                 surface < provenance.local_surface_count; surface++) {
                PyObject* value = PyLong_FromLong(
                    provenance.local_surface_ids[surface]);
                if (!value) {
                    Py_DECREF(local_surfaces); Py_DECREF(local_point);
                    Py_DECREF(local_direction); Py_DECREF(item);
                    goto failed;
                }
                PyList_SET_ITEM(local_surfaces, (Py_ssize_t)surface, value);
            }
            PyObject* receipt = Py_BuildValue(
                "{s:I,s:i,s:i,s:i,s:K,s:K,s:K,s:K,s:K,s:K,s:O,s:O,s:N,s:O}",
                "provenance_flags", provenance.flags,
                "active_cell_id", provenance.active_cell_id,
                "active_universe_id", provenance.active_universe_id,
                "active_depth", provenance.active_depth,
                "active_occurrence_key",
                    (unsigned long long)provenance.active_occurrence_key,
                "active_parent_occurrence_key",
                    (unsigned long long)provenance.active_parent_occurrence_key,
                "before_occurrence_key",
                    (unsigned long long)provenance.before_occurrence_key,
                "before_parent_occurrence_key",
                    (unsigned long long)provenance.before_parent_occurrence_key,
                "after_occurrence_key",
                    (unsigned long long)provenance.after_occurrence_key,
                "after_parent_occurrence_key",
                    (unsigned long long)provenance.after_parent_occurrence_key,
                "local_point", local_point,
                "local_direction", local_direction,
                "local_surface_ids", local_surfaces,
                "local_surface_complete",
                    provenance.local_surface_complete ? Py_True : Py_False);
            Py_DECREF(local_point);
            Py_DECREF(local_direction);
            if (!receipt || PyDict_Update(item, receipt) != 0) {
                Py_XDECREF(receipt); Py_DECREF(item);
                goto failed;
            }
            Py_DECREF(receipt);
        }
        PyList_SET_ITEM(out, (Py_ssize_t)i, item);
    }
    alea_ray_boundary_event_query_result_destroy(result);
    return out;

failed:
    Py_XDECREF(out);
    alea_ray_boundary_event_query_result_destroy(result);
    return NULL;
}

static PyObject* PyAleaSystem_find_overlaps(PyAleaSystemObject* self, PyObject* args) {
    size_t max_pairs = 100;

    if (!PyArg_ParseTuple(args, "|n", &max_pairs)) {
        return NULL;
    }

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    int* pairs = malloc(max_pairs * 2 * sizeof(int));
    if (!pairs) {
        PyErr_NoMemory();
        return NULL;
    }

    int n;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    n = alea_find_overlaps(self->sys, pairs, max_pairs);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) {
        free(pairs);
        return NULL;
    }

    PyObject* list = PyList_New(n);
    if (!list) {
        free(pairs);
        return NULL;
    }

    for (int i = 0; i < n; i++) {
        PyObject* pair = Py_BuildValue("(ii)", pairs[i*2], pairs[i*2+1]);
        if (!pair) {
            free(pairs);
            Py_DECREF(list);
            return NULL;
        }
        PyList_SET_ITEM(list, i, pair);
    }

    free(pairs);
    return list;
}

/* ============================================================================
 * PyAleaSystem Methods - Volume Estimation
 * ============================================================================ */

static PyObject* PyAleaSystem_compute_bounding_sphere(PyAleaSystemObject* self, PyObject* args) {
    double tol = 1.0;
    if (!PyArg_ParseTuple(args, "|d", &tol)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    double cx, cy, cz, radius;
    int result;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    result = alea_compute_bounding_sphere(self->sys, tol, &cx, &cy, &cz, &radius);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) return NULL;

    if (result < 0) {
        PyErr_SetString(PyExc_RuntimeError, "No bounded cells found");
        return NULL;
    }
    return Py_BuildValue("(dddd)", cx, cy, cz, radius);
}

static const char* cell_volume_bounds_source_name(
        alea_cell_volume_bounds_source_t source) {
    switch (source) {
        case ALEA_CELL_VOLUME_BOUNDS_EXPLICIT: return "explicit";
        case ALEA_CELL_VOLUME_BOUNDS_STORED: return "stored";
        case ALEA_CELL_VOLUME_BOUNDS_PLANE_CONSTRAINTS: return "plane_constraints";
        case ALEA_CELL_VOLUME_BOUNDS_ADAPTIVE_SEARCH: return "adaptive_search";
        default: return "unknown";
    }
}

static int parse_cell_volume_bounds(PyObject* obj, alea_bbox_t* bounds) {
    PyObject* seq = PySequence_Fast(obj, "bounds must be a sequence of six numbers");
    if (!seq) return -1;
    if (PySequence_Fast_GET_SIZE(seq) != 6) {
        Py_DECREF(seq);
        PyErr_SetString(PyExc_ValueError, "bounds must contain exactly six values");
        return -1;
    }
    double values[6];
    for (Py_ssize_t i = 0; i < 6; i++) {
        values[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq, i));
        if (PyErr_Occurred()) { Py_DECREF(seq); return -1; }
    }
    Py_DECREF(seq);
    *bounds = (alea_bbox_t){values[0], values[1], values[2],
                            values[3], values[4], values[5]};
    return 0;
}

static PyObject* PyAleaSystem_estimate_cell_volume(
        PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    Py_ssize_t cell_index, workers = 0;
    PyObject* bounds_obj = Py_None;
    double relative_tolerance = 1e-3, absolute_tolerance = 0.0, min_size = 0.0;
    int max_depth = 10, samples_per_axis = 2;
    unsigned long long max_parallel_scratch_bytes = 64u * 1024u * 1024u;
    static char* kwlist[] = {
        "cell_index", "bounds", "relative_tolerance", "absolute_tolerance",
        "max_depth", "min_size", "samples_per_axis", "workers",
        "max_parallel_scratch_bytes", NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "n|OddidinK", kwlist,
            &cell_index, &bounds_obj, &relative_tolerance, &absolute_tolerance,
            &max_depth, &min_size, &samples_per_axis, &workers,
            &max_parallel_scratch_bytes)) return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    if (cell_index < 0 || workers < 0) {
        PyErr_SetString(PyExc_ValueError, "cell_index and workers must be non-negative");
        return NULL;
    }

    alea_cell_volume_options_t options;
    alea_cell_volume_options_init(&options);
    options.relative_tolerance = relative_tolerance;
    options.absolute_tolerance = absolute_tolerance;
    options.max_depth = max_depth;
    options.min_size = min_size;
    options.samples_per_axis = samples_per_axis;
    options.requested_workers = (size_t)workers;
    options.max_parallel_scratch_bytes = (uint64_t)max_parallel_scratch_bytes;
    if (bounds_obj != Py_None) {
        if (parse_cell_volume_bounds(bounds_obj, &options.bounds) != 0) return NULL;
        options.has_bounds = true;
    }

    alea_cell_volume_result_t result;
    int rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    rc = alea_cell_estimate_volume(self->sys, (size_t)cell_index, &options, &result);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) return NULL;
    if (rc != 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    PyObject* stats = Py_BuildValue(
        "{s:n,s:n,s:n,s:n,s:n,s:n,s:K,s:K}",
        "total_nodes", (Py_ssize_t)result.total_nodes,
        "inside_nodes", (Py_ssize_t)result.inside_nodes,
        "outside_nodes", (Py_ssize_t)result.outside_nodes,
        "unresolved_leaf_nodes", (Py_ssize_t)result.unresolved_leaf_nodes,
        "max_depth_reached", (Py_ssize_t)result.max_depth_reached,
        "frontier_task_count", (Py_ssize_t)result.frontier_task_count,
        "scratch_bytes_per_worker", (unsigned long long)result.scratch_bytes_per_worker,
        "reserved_parallel_scratch_bytes",
            (unsigned long long)result.reserved_parallel_scratch_bytes);
    if (!stats) return NULL;
    return Py_BuildValue(
        "{s:d,s:d,s:d,s:d,s:d,s:(dddddd),s:s,s:O,s:O,s:O,s:n,s:n,s:n,s:N}",
        "volume", result.volume,
        "lower_bound", result.lower_bound,
        "upper_bound", result.upper_bound,
        "unresolved_volume", result.unresolved_volume,
        "relative_uncertainty", result.relative_uncertainty,
        "bounds", result.bounds.min_x, result.bounds.max_x,
            result.bounds.min_y, result.bounds.max_y,
            result.bounds.min_z, result.bounds.max_z,
        "bounds_source", cell_volume_bounds_source_name(result.bounds_source),
        "complete_cell_domain", result.complete_cell_domain ? Py_True : Py_False,
        "converged", result.converged ? Py_True : Py_False,
        "resource_limit_reached", result.resource_limit_reached ? Py_True : Py_False,
        "bounds_search_expansions", (Py_ssize_t)result.bounds_search_expansions,
        "requested_workers", (Py_ssize_t)result.requested_workers,
        "actual_workers", (Py_ssize_t)result.actual_workers,
        "stats", stats);
}

/* Convert a concrete hierarchical volume path to a Python dict describing its
 * identity and placement. */
static PyObject* volume_path_to_dict(const alea_volume_path_t* p) {
    PyObject* ancestor_cells = PyList_New(p->ancestor_count);
    PyObject* ancestor_universes = PyList_New(p->ancestor_count);
    if (!ancestor_cells || !ancestor_universes) {
        Py_XDECREF(ancestor_cells); Py_XDECREF(ancestor_universes);
        return NULL;
    }
    for (int i = 0; i < p->ancestor_count; i++) {
        PyList_SET_ITEM(ancestor_cells, i, PyLong_FromLong(p->ancestor_cell_indices[i]));
        PyList_SET_ITEM(ancestor_universes, i, PyLong_FromLong(p->ancestor_universe_ids[i]));
    }

    PyObject* lattice_steps = PyList_New(p->lattice_step_count);
    if (!lattice_steps) {
        Py_DECREF(ancestor_cells); Py_DECREF(ancestor_universes);
        return NULL;
    }
    for (int i = 0; i < p->lattice_step_count; i++) {
        const alea_volume_lattice_step_t* s = &p->lattice_steps[i];
        PyObject* step = Py_BuildValue("{s:i,s:i,s:(iii),s:i}",
            "lattice_cell_index", s->lattice_cell_index,
            "fill_universe", s->fill_universe,
            "ijk", s->i, s->j, s->k,
            "linear_index", s->linear_index);
        if (!step) {
            Py_DECREF(ancestor_cells); Py_DECREF(ancestor_universes); Py_DECREF(lattice_steps);
            return NULL;
        }
        PyList_SET_ITEM(lattice_steps, i, step);
    }

    PyObject* world_to_local = PyTuple_New(12);
    if (!world_to_local) {
        Py_DECREF(ancestor_cells);
        Py_DECREF(ancestor_universes);
        Py_DECREF(lattice_steps);
        return NULL;
    }
    for (int i = 0; i < 12; i++) {
        PyObject* value = PyFloat_FromDouble(p->world_to_local[i]);
        if (!value) {
            Py_DECREF(ancestor_cells);
            Py_DECREF(ancestor_universes);
            Py_DECREF(lattice_steps);
            Py_DECREF(world_to_local);
            return NULL;
        }
        PyTuple_SET_ITEM(world_to_local, i, value);
    }

    return Py_BuildValue(
        "{s:K,s:i,s:i,s:i,s:i,s:i,s:O,s:N,s:N,s:N,s:N}",
        "path_id", (unsigned long long)p->path_id,
        "cell_index", p->terminal_cell_index,
        "cell_id", p->terminal_cell_id,
        "material_id", p->material_id,
        "universe_id", p->universe_id,
        "depth", p->depth,
        "truncated", p->truncated ? Py_True : Py_False,
        "ancestor_cell_indices", ancestor_cells,
        "ancestor_universe_ids", ancestor_universes,
        "lattice_steps", lattice_steps,
        "world_to_local", world_to_local);
}

typedef struct {
    PyObject* callable;
    int failed;
} volume_progress_context_t;

static int volume_progress_to_python(size_t completed, size_t maximum,
                                     double maximum_error, void* user_data) {
    volume_progress_context_t* context = user_data;
    PyGILState_STATE gil = PyGILState_Ensure();
    if (PyErr_CheckSignals() != 0) {
        context->failed = 1;
        PyGILState_Release(gil);
        return 1;
    }
    PyObject* result = PyObject_CallFunction(
        context->callable, "KKd", (unsigned long long)completed,
        (unsigned long long)maximum, maximum_error);
    if (!result) {
        context->failed = 1;
        PyGILState_Release(gil);
        return 1;
    }
    int cancel = PyObject_IsTrue(result);
    Py_DECREF(result);
    if (cancel < 0) {
        context->failed = 1;
        cancel = 1;
    }
    PyGILState_Release(gil);
    return cancel;
}

static PyObject* PyAleaSystem_estimate_volumes(
        PyAleaSystemObject* self, PyObject* args, PyObject* kwargs) {
    int n_rays = 100000;
    unsigned long long seed = 42;
    const char* rng_name = "philox4x32-10";
    Py_ssize_t workers = 0;
    PyObject* target_obj = Py_None;
    PyObject* max_rays_obj = Py_None;
    Py_ssize_t batch_size = 10000;
    PyObject* progress_obj = Py_None;
    unsigned long long max_parallel_scratch_bytes = 0;
    static char* kwlist[] = {
        "n_rays", "seed", "workers", "target_rel_error", "max_rays",
        "batch_size", "progress", "rng", "max_parallel_scratch_bytes", NULL
    };
    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, "|iKnOOnOsK:estimate_volumes", kwlist,
            &n_rays, &seed, &workers, &target_obj, &max_rays_obj,
            &batch_size, &progress_obj, &rng_name,
            &max_parallel_scratch_bytes)) {
        return NULL;
    }
    if (n_rays <= 0 || workers < 0 || batch_size < 0) {
        PyErr_SetString(PyExc_ValueError,
                        "n_rays must be positive; workers and batch_size cannot be negative");
        return NULL;
    }
    if (max_parallel_scratch_bytes > (unsigned long long)SIZE_MAX) {
        PyErr_SetString(PyExc_OverflowError,
                        "max_parallel_scratch_bytes exceeds size_t");
        return NULL;
    }
    size_t maximum_rays = (size_t)n_rays;
    if (max_rays_obj != Py_None) {
        maximum_rays = PyLong_AsSize_t(max_rays_obj);
        if (PyErr_Occurred()) return NULL;
        if (maximum_rays == 0) {
            PyErr_SetString(PyExc_ValueError, "max_rays must be positive");
            return NULL;
        }
    }
    double target_rel_error = 0.0;
    if (target_obj != Py_None) {
        target_rel_error = PyFloat_AsDouble(target_obj);
        if (PyErr_Occurred()) return NULL;
        if (!(target_rel_error > 0.0 && target_rel_error <= 1.0)) {
            PyErr_SetString(PyExc_ValueError,
                            "target_rel_error must be in (0, 1]");
            return NULL;
        }
    }
    if (progress_obj != Py_None && !PyCallable_Check(progress_obj)) {
        PyErr_SetString(PyExc_TypeError, "progress must be callable or None");
        return NULL;
    }
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }
    if (ensure_query_acceleration(self) < 0) return NULL;

    size_t count = alea_volume_path_count(self->sys);
    if (count == 0) {
        PyErr_SetString(PyExc_RuntimeError, "No hierarchical volume paths available");
        return NULL;
    }

    double* volumes = calloc(count, sizeof(double));
    double* rel_errors = calloc(count, sizeof(double));
    if (!volumes || !rel_errors) {
        free(volumes); free(rel_errors);
        return PyErr_NoMemory();
    }

    alea_volume_estimate_options_t options;
    alea_volume_estimate_options_init(&options);
    if (strcmp(rng_name, "philox4x32-10") == 0) {
        options.rng_algorithm = ALEA_RNG_PHILOX4X32_10;
    } else if (strcmp(rng_name, "legacy-lcg32") == 0) {
        options.rng_algorithm = ALEA_RNG_LEGACY_LCG;
    } else {
        free(volumes); free(rel_errors);
        PyErr_Format(PyExc_ValueError,
                     "unknown RNG '%s'; expected 'philox4x32-10' or 'legacy-lcg32'",
                     rng_name);
        return NULL;
    }
    options.max_rays = maximum_rays;
    options.seed = (uint64_t)seed;
    options.requested_workers = (size_t)workers;
    options.batch_size = (size_t)batch_size;
    if (max_parallel_scratch_bytes)
        options.max_parallel_scratch_bytes =
            (size_t)max_parallel_scratch_bytes;
    options.target_rel_error = target_rel_error;
    volume_progress_context_t progress_context = {
        .callable = progress_obj,
        .failed = 0
    };
    if (progress_obj != Py_None) {
        options.progress = volume_progress_to_python;
        options.progress_user_data = &progress_context;
    }
    alea_volume_estimate_stats_t stats;
    int result;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    result = alea_estimate_volumes_ex(
        self->sys, &options, volumes, rel_errors, &stats);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) {
        free(volumes); free(rel_errors);
        return NULL;
    }

    if (progress_context.failed) {
        free(volumes); free(rel_errors);
        if (!PyErr_Occurred())
            PyErr_SetString(PyExc_RuntimeError, "volume progress callback failed");
        return NULL;
    }

    if (result < 0) {
        free(volumes); free(rel_errors);
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    PyObject* vol_list = PyList_New(count);
    PyObject* err_list = PyList_New(count);
    PyObject* path_list = PyList_New(count);
    if (!vol_list || !err_list || !path_list) {
        Py_XDECREF(vol_list); Py_XDECREF(err_list); Py_XDECREF(path_list);
        free(volumes); free(rel_errors);
        return NULL;
    }
    for (size_t i = 0; i < count; i++) {
        PyList_SET_ITEM(vol_list, i, PyFloat_FromDouble(volumes[i]));
        PyList_SET_ITEM(err_list, i, PyFloat_FromDouble(rel_errors[i]));
        alea_volume_path_t path;
        PyObject* pd = alea_volume_paths_get_range(
            self->sys, i, &path, 1) == 1
            ? volume_path_to_dict(&path) : (Py_INCREF(Py_None), Py_None);
        if (!pd) {
            Py_DECREF(vol_list); Py_DECREF(err_list); Py_DECREF(path_list);
            free(volumes); free(rel_errors);
            return NULL;
        }
        PyList_SET_ITEM(path_list, i, pd);
    }
    free(volumes);
    free(rel_errors);

    return Py_BuildValue(
        "{s:N,s:N,s:N,s:K,s:K,s:K,s:K,s:d,s:O,s:O,s:K,s:s,s:I,"
        "s:K,s:K,s:K}",
        "volumes", vol_list,
        "rel_errors", err_list,
        "paths", path_list,
        "rays_completed", (unsigned long long)stats.rays_completed,
        "requested_workers", (unsigned long long)stats.requested_workers,
        "actual_workers", (unsigned long long)stats.actual_workers,
        "batch_size", (unsigned long long)stats.batch_size,
        "maximum_relative_error", stats.maximum_relative_error,
        "converged", stats.converged ? Py_True : Py_False,
        "cancelled", stats.cancelled ? Py_True : Py_False,
        "seed", (unsigned long long)stats.seed,
        "rng_algorithm", alea_rng_algorithm_name(stats.rng_algorithm),
        "rng_address_version", stats.rng_address_version,
        "parallel_scratch_limit_bytes",
            (unsigned long long)stats.parallel_scratch_limit_bytes,
        "parallel_scratch_bytes",
            (unsigned long long)stats.parallel_scratch_bytes,
        "worker_scratch_bytes",
            (unsigned long long)stats.worker_scratch_bytes);
}

static PyObject* PyAleaSystem_volume_path_at_point(PyAleaSystemObject* self, PyObject* args) {
    double x, y, z;
    if (!PyArg_ParseTuple(args, "ddd", &x, &y, &z)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }
    if (ensure_query_acceleration(self) < 0) return NULL;

    alea_volume_path_t path;
    int rc = alea_volume_path_at_point(self->sys, x, y, z, &path);
    if (rc < 0) Py_RETURN_NONE;
    return volume_path_to_dict(&path);
}

static PyObject* PyAleaSystem_volume_path_resolve_at_point(
        PyAleaSystemObject* self, PyObject* args) {
    double x, y, z;
    if (!PyArg_ParseTuple(args, "ddd", &x, &y, &z)) return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    if (ensure_query_acceleration(self) < 0) return NULL;

    alea_volume_path_t path;
    int rc = alea_volume_path_resolve_at_point(self->sys, x, y, z, &path);
    if (rc < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }
    if (rc == 0) Py_RETURN_NONE;
    return volume_path_to_dict(&path);
}

static PyObject* PyAleaSystem_volume_path_resolve_cell_at_point(
        PyAleaSystemObject* self, PyObject* args) {
    double x, y, z;
    int cell_id, universe_id;
    if (!PyArg_ParseTuple(args, "dddii", &x, &y, &z,
                          &cell_id, &universe_id)) return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    if (ensure_query_acceleration(self) < 0) return NULL;

    alea_volume_path_t path;
    int rc = alea_volume_path_resolve_cell_at_point(
        self->sys, x, y, z, cell_id, universe_id, &path);
    if (rc < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    const char* status = rc == 0 ? "none" : rc == 1 ? "unique" : "multiple";
    PyObject* path_dict = rc == 1 ? volume_path_to_dict(&path) : NULL;
    if (rc == 1 && !path_dict) return NULL;
    if (rc != 1) {
        path_dict = Py_None;
        Py_INCREF(path_dict);
    }
    return Py_BuildValue("{s:s,s:N}", "status", status, "path", path_dict);
}

static PyObject* PyAleaSystem_volume_path_resolve_cell_point_sets(
        PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    PyObject* sets_object;
    Py_ssize_t requested_workers = 0;
    unsigned long long max_parallel_scratch_bytes = 64ULL * 1024ULL * 1024ULL;
    static char* kwlist[] = {
        "point_sets", "requested_workers", "max_parallel_scratch_bytes", NULL};
    if (!PyArg_ParseTupleAndKeywords(
            args, kwds, "O|nK", kwlist, &sets_object, &requested_workers,
            &max_parallel_scratch_bytes)) return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    if (requested_workers < 0) {
        PyErr_SetString(PyExc_ValueError,
                        "requested_workers must be non-negative");
        return NULL;
    }

    PyObject* outer = PySequence_Fast(
        sets_object, "point_sets must be a sequence");
    if (!outer) return NULL;
    Py_ssize_t set_count_py = PySequence_Fast_GET_SIZE(outer);
    size_t set_count = (size_t)set_count_py;
    alea_volume_path_point_set_t* sets =
        set_count ? calloc(set_count, sizeof(*sets)) : NULL;
    alea_volume_path_point_set_result_t* native_results =
        set_count ? calloc(set_count, sizeof(*native_results)) : NULL;
    PyObject** point_sequences =
        set_count ? calloc(set_count, sizeof(*point_sequences)) : NULL;
    if (set_count && (!sets || !native_results || !point_sequences)) {
        Py_DECREF(outer);
        free(sets); free(native_results); free(point_sequences);
        return PyErr_NoMemory();
    }

    size_t point_count = 0;
    for (size_t set = 0; set < set_count; set++) {
        PyObject* descriptor = PySequence_Fast(
            PySequence_Fast_GET_ITEM(outer, (Py_ssize_t)set),
            "each point set must be (points, cell_id, universe_id)");
        if (!descriptor) goto parse_fail;
        if (PySequence_Fast_GET_SIZE(descriptor) != 3) {
            Py_DECREF(descriptor);
            PyErr_SetString(PyExc_ValueError,
                            "each point set must contain three items");
            goto parse_fail;
        }
        PyObject* points = PySequence_Fast(
            PySequence_Fast_GET_ITEM(descriptor, 0),
            "point-set points must be a sequence");
        long cell_id = PyLong_AsLong(PySequence_Fast_GET_ITEM(descriptor, 1));
        long universe_id = PyLong_AsLong(
            PySequence_Fast_GET_ITEM(descriptor, 2));
        Py_DECREF(descriptor);
        if (!points || PyErr_Occurred()) {
            Py_XDECREF(points);
            goto parse_fail;
        }
        if (cell_id < INT_MIN || cell_id > INT_MAX ||
            universe_id < INT_MIN || universe_id > INT_MAX) {
            Py_DECREF(points);
            PyErr_SetString(PyExc_OverflowError,
                            "point-set cell or universe ID is out of range");
            goto parse_fail;
        }
        size_t count = (size_t)PySequence_Fast_GET_SIZE(points);
        if (count > SIZE_MAX - point_count) {
            Py_DECREF(points);
            PyErr_SetString(PyExc_OverflowError,
                            "point-set point count overflows");
            goto parse_fail;
        }
        point_sequences[set] = points;
        sets[set].point_offset = point_count;
        sets[set].point_count = count;
        sets[set].target_cell_id = (int)cell_id;
        sets[set].target_universe_id = (int)universe_id;
        point_count += count;
    }
    if (point_count > SIZE_MAX / (3 * sizeof(double))) {
        PyErr_SetString(PyExc_OverflowError, "point-set coordinate storage overflows");
        goto parse_fail;
    }
    double* points_xyz = point_count
        ? malloc(3 * point_count * sizeof(*points_xyz)) : NULL;
    if (point_count && !points_xyz) {
        PyErr_NoMemory();
        goto parse_fail;
    }
    for (size_t set = 0; set < set_count; set++) {
        PyObject* points = point_sequences[set];
        for (size_t item = 0; item < sets[set].point_count; item++) {
            PyObject* point = PySequence_Fast(
                PySequence_Fast_GET_ITEM(points, (Py_ssize_t)item),
                "each point must contain three coordinates");
            if (!point) {
                free(points_xyz);
                goto parse_fail;
            }
            if (PySequence_Fast_GET_SIZE(point) != 3) {
                Py_DECREF(point);
                free(points_xyz);
                PyErr_SetString(PyExc_ValueError,
                                "each point must contain three coordinates");
                goto parse_fail;
            }
            size_t index = sets[set].point_offset + item;
            for (int axis = 0; axis < 3; axis++) {
                points_xyz[3 * index + (size_t)axis] = PyFloat_AsDouble(
                    PySequence_Fast_GET_ITEM(point, axis));
            }
            Py_DECREF(point);
            if (PyErr_Occurred()) {
                free(points_xyz);
                goto parse_fail;
            }
        }
    }
    for (size_t set = 0; set < set_count; set++) {
        Py_DECREF(point_sequences[set]);
        point_sequences[set] = NULL;
    }
    Py_DECREF(outer);
    outer = NULL;

    alea_volume_path_point_set_batch_stats_t stats;
    int rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    rc = alea_volume_path_resolve_cell_point_sets(
        self->sys, points_xyz, point_count, sets, set_count,
        (size_t)requested_workers, (uint64_t)max_parallel_scratch_bytes,
        native_results, &stats);
    Py_END_ALLOW_THREADS
    free(points_xyz);
    free(sets);
    if (restore_sigint(old_sigint)) {
        free(native_results); free(point_sequences);
        return NULL;
    }
    if (rc != 0) {
        free(native_results); free(point_sequences);
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    PyObject* results = PyList_New(set_count_py);
    if (!results) {
        free(native_results); free(point_sequences);
        return NULL;
    }
    for (size_t set = 0; set < set_count; set++) {
        const alea_volume_path_point_set_result_t* native = &native_results[set];
        const char* status = native->status == ALEA_VOLUME_PATH_NONE
            ? "none" : native->status == ALEA_VOLUME_PATH_UNIQUE
            ? "unique" : "multiple";
        PyObject* path = native->status == ALEA_VOLUME_PATH_UNIQUE
            ? volume_path_to_dict(&native->unique_path) : NULL;
        if (native->status != ALEA_VOLUME_PATH_UNIQUE) {
            path = Py_None;
            Py_INCREF(path);
        }
        if (!path) {
            Py_DECREF(results); free(native_results); free(point_sequences);
            return NULL;
        }
        PyObject* item = Py_BuildValue(
            "{s:s,s:n,s:n,s:N}",
            "status", status,
            "tested_point_count", (Py_ssize_t)native->tested_point_count,
            "matching_occurrence_count",
                (Py_ssize_t)native->matching_occurrence_count,
            "path", path);
        if (!item) {
            Py_DECREF(results); free(native_results); free(point_sequences);
            return NULL;
        }
        PyList_SET_ITEM(results, (Py_ssize_t)set, item);
    }
    free(native_results);
    free(point_sequences);
    PyObject* stats_dict = Py_BuildValue(
        "{s:n,s:n,s:n,s:n,s:K,s:K}",
        "point_set_count", (Py_ssize_t)stats.point_set_count,
        "completed_point_set_count",
            (Py_ssize_t)stats.completed_point_set_count,
        "requested_workers", (Py_ssize_t)stats.requested_workers,
        "actual_workers", (Py_ssize_t)stats.actual_workers,
        "reserved_scratch_bytes_per_worker",
            (unsigned long long)stats.reserved_scratch_bytes_per_worker,
        "reserved_parallel_scratch_bytes",
            (unsigned long long)stats.reserved_parallel_scratch_bytes);
    if (!stats_dict) {
        Py_DECREF(results);
        return NULL;
    }
    return Py_BuildValue("{s:N,s:N}", "results", results, "stats", stats_dict);

parse_fail:
    for (size_t set = 0; set < set_count; set++)
        Py_XDECREF(point_sequences[set]);
    Py_XDECREF(outer);
    free(sets); free(native_results); free(point_sequences);
    return NULL;
}

static int transform_evidence_vector(PyObject* object, const char* name,
                                     double out[3]) {
    PyObject* sequence = PySequence_Fast(object, name);
    if (!sequence) return -1;
    if (PySequence_Fast_GET_SIZE(sequence) != 3) {
        Py_DECREF(sequence);
        PyErr_Format(PyExc_ValueError, "%s must contain exactly 3 values", name);
        return -1;
    }
    for (int axis = 0; axis < 3; axis++) {
        out[axis] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(sequence, axis));
        if (PyErr_Occurred()) {
            Py_DECREF(sequence);
            return -1;
        }
    }
    Py_DECREF(sequence);
    return 0;
}

static PyObject* PyAleaSystem_volume_path_resolve_cell_from_transform_evidence(
        PyAleaSystemObject* self, PyObject* args) {
    int cell_id, universe_id;
    PyObject *objects[8];
    if (!PyArg_ParseTuple(args, "iiOOOOOOOO", &cell_id, &universe_id,
                          &objects[0], &objects[1], &objects[2], &objects[3],
                          &objects[4], &objects[5], &objects[6], &objects[7])) {
        return NULL;
    }
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    if (ensure_query_acceleration(self) < 0) return NULL;

    const char* names[8] = {
        "world_point", "world_direction", "local_point", "local_direction",
        "world_point_precision", "world_direction_precision",
        "local_point_precision", "local_direction_precision",
    };
    double vectors[8][3];
    for (int i = 0; i < 8; i++) {
        if (transform_evidence_vector(objects[i], names[i], vectors[i]) != 0)
            return NULL;
    }

    alea_volume_path_t path;
    int rc = alea_volume_path_resolve_cell_from_transform_evidence(
        self->sys, cell_id, universe_id,
        vectors[0], vectors[1], vectors[2], vectors[3],
        vectors[4], vectors[5], vectors[6], vectors[7], &path);
    if (rc < 0) {
        PyErr_SetString(PyExc_ValueError,
                        "Invalid or unavailable transform evidence");
        return NULL;
    }

    const char* status = rc == 0 ? "none" : rc == 1 ? "unique" : "multiple";
    PyObject* path_dict = rc == 1 ? volume_path_to_dict(&path) : NULL;
    if (rc == 1 && !path_dict) return NULL;
    if (rc != 1) {
        path_dict = Py_None;
        Py_INCREF(path_dict);
    }
    return Py_BuildValue("{s:s,s:N}", "status", status, "path", path_dict);
}

/* ============================================================================
 * PyAleaSystem Methods - Cell-Aware Raycast
 * ============================================================================ */

static PyObject* PyAleaSystem_raycast_cell_aware(PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    double ox, oy, oz, dx, dy, dz;
    double t_max = 0.0;
    static char* kwlist[] = {"ox", "oy", "oz", "dx", "dy", "dz", "t_max", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "dddddd|d", kwlist,
                                     &ox, &oy, &oz, &dx, &dy, &dz, &t_max)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }
    if (ensure_query_acceleration(self) < 0) return NULL;

    alea_raycast_result_t* result = alea_raycast_result_create();
    if (!result) { PyErr_SetString(PyExc_MemoryError, "Failed to allocate raycast result"); return NULL; }

    int rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    rc = alea_raycast_cell_aware(self->sys, ox, oy, oz, dx, dy, dz, t_max, result);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) { alea_raycast_result_free(result); return NULL; }

    if (rc < 0) {
        alea_raycast_result_free(result);
        PyErr_SetString(PyExc_RuntimeError, "Cell-aware raycast failed");
        return NULL;
    }

    size_t count = alea_raycast_segment_count(result);
    PyObject* segments = PyList_New(count);
    for (size_t i = 0; i < count; i++) {
        double t_enter, t_exit, density;
        int cell_id, material_id, enter_surface_id, exit_surface_id;
        alea_raycast_segment_get(
            result, i, &t_enter, &t_exit, &cell_id, &material_id, &density,
            &enter_surface_id, &exit_surface_id);
        uint8_t resolution_flags;
        if (alea_raycast_segment_resolution_flags(result, i, &resolution_flags) != 0) {
            Py_DECREF(segments);
            alea_raycast_result_free(result);
            PyErr_SetString(PyExc_RuntimeError, "Failed to read ray segment flags");
            return NULL;
        }
        PyObject* item = Py_BuildValue("{s:d,s:d,s:i,s:i,s:d,s:i,s:i,s:i}",
            "t_enter", t_enter, "t_exit", t_exit,
            "cell_id", cell_id, "material_id", material_id,
            "density", density,
            "enter_surface_id", enter_surface_id,
            "exit_surface_id", exit_surface_id,
            "resolution_flags", resolution_flags);
        PyList_SET_ITEM(segments, i, item);
    }
    alea_raycast_result_free(result);
    return segments;
}

/* Fast hierarchical segment tracer. Walks the hierarchical spatial index and
 * tests only the surfaces of cells actually along the ray, instead of the
 * global surface BVH used by alea_raycast(). This is dramatically faster on
 * large models (hundreds of thousands of surfaces). It produces segments
 * (with boundary surface IDs where available) but not a full global hit list,
 * which the Python layer does not expose anyway. */
static PyObject* PyAleaSystem_raycast_hier_fast(PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    double ox, oy, oz, dx, dy, dz;
    double t_max = 0.0;
    int include_paths = 0;
    static char* kwlist[] = {"ox", "oy", "oz", "dx", "dy", "dz", "t_max", "include_paths", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "dddddd|dp", kwlist,
                                     &ox, &oy, &oz, &dx, &dy, &dz, &t_max,
                                     &include_paths)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }
    if (ensure_query_acceleration(self) < 0) return NULL;

    alea_raycast_result_t* result = alea_raycast_result_create();
    if (!result) { PyErr_SetString(PyExc_MemoryError, "Failed to allocate raycast result"); return NULL; }
    alea_raycast_result_set_path_capture(result, include_paths);

    int rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    rc = alea_raycast_hier_fast_segments(self->sys, ox, oy, oz, dx, dy, dz, t_max, result);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) { alea_raycast_result_free(result); return NULL; }

    if (rc < 0) {
        alea_raycast_result_free(result);
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    size_t count = alea_raycast_segment_count(result);
    PyObject* segments = PyList_New(count);
    for (size_t i = 0; i < count; i++) {
        double t_enter, t_exit, density;
        int cell_id, material_id, enter_surface_id, exit_surface_id;
        alea_raycast_segment_get(
            result, i, &t_enter, &t_exit, &cell_id, &material_id, &density,
            &enter_surface_id, &exit_surface_id);
        uint8_t resolution_flags;
        if (alea_raycast_segment_resolution_flags(result, i, &resolution_flags) != 0) {
            Py_DECREF(segments);
            alea_raycast_result_free(result);
            PyErr_SetString(PyExc_RuntimeError, "Failed to read ray segment flags");
            return NULL;
        }
        PyObject* item = Py_BuildValue("{s:d,s:d,s:i,s:i,s:d,s:i,s:i,s:i}",
            "t_enter", t_enter, "t_exit", t_exit,
            "cell_id", cell_id, "material_id", material_id,
            "density", density,
            "enter_surface_id", enter_surface_id,
            "exit_surface_id", exit_surface_id,
            "resolution_flags", resolution_flags);
        if (!item || (include_paths && raycast_attach_path(item, result, i) != 0)) {
            Py_XDECREF(item);
            Py_DECREF(segments);
            alea_raycast_result_free(result);
            return NULL;
        }
        PyList_SET_ITEM(segments, i, item);
    }
    alea_raycast_result_free(result);
    return segments;
}

/* Batch form of the segment-only hierarchical tracer.  Traversal happens
 * without the GIL; Python objects are materialised only after all rays finish. */
static PyObject* PyAleaSystem_raycast_hier_fast_batch(PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    PyObject *origins_obj, *directions_obj;
    double t_max = 0.0;
    int include_paths = 0;
    static char* kwlist[] = {"origins", "directions", "t_max", "include_paths", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OO|dp", kwlist,
                                     &origins_obj, &directions_obj, &t_max,
                                     &include_paths)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }
    if (ensure_query_acceleration(self) < 0) return NULL;
    PyArrayObject *origins = (PyArrayObject*)PyArray_FROM_OTF(origins_obj, NPY_DOUBLE, NPY_ARRAY_CARRAY_RO);
    PyArrayObject *directions = (PyArrayObject*)PyArray_FROM_OTF(directions_obj, NPY_DOUBLE, NPY_ARRAY_CARRAY_RO);
    if (!origins || !directions) { Py_XDECREF(origins); Py_XDECREF(directions); return NULL; }
    if (PyArray_NDIM(origins) != 2 || PyArray_NDIM(directions) != 2 ||
        PyArray_DIM(origins, 1) != 3 || PyArray_DIM(directions, 1) != 3 ||
        PyArray_DIM(origins, 0) != PyArray_DIM(directions, 0)) {
        Py_DECREF(origins); Py_DECREF(directions);
        PyErr_SetString(PyExc_ValueError, "origins and directions must have shape (n_rays, 3)"); return NULL;
    }
    Py_ssize_t n = PyArray_DIM(origins, 0);
    double *o = PyArray_DATA(origins), *d = PyArray_DATA(directions);
    alea_raycast_batch_result_t* result = alea_raycast_batch_result_create();
    if (!result) {
        Py_DECREF(origins); Py_DECREF(directions);
        return PyErr_NoMemory();
    }
    alea_raycast_batch_options_t options = {
        .struct_size = sizeof(options),
        .fields = ALEA_RAY_BATCH_MATERIAL | ALEA_RAY_BATCH_DENSITY |
                  ALEA_RAY_BATCH_SURFACES | ALEA_RAY_BATCH_RESOLUTION_FLAGS |
                  (include_paths ? ALEA_RAY_BATCH_FULL_PATHS : 0),
        .projected_depth = -1,
    };
    int rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    rc = alea_raycast_hier_batch(
        self->sys, o, d, (size_t)n, t_max, &options, result);
    Py_END_ALLOW_THREADS
    Py_DECREF(origins); Py_DECREF(directions);
    if (restore_sigint(old_sigint) || rc != 0) {
        alea_raycast_batch_result_destroy(result);
        if (!PyErr_Occurred()) PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    const uint64_t* ray_offsets = alea_raycast_batch_ray_offsets(result);
    const double* enters = alea_raycast_batch_t_enter(result);
    const double* exits = alea_raycast_batch_t_exit(result);
    const int32_t* cells = alea_raycast_batch_cell_ids(result);
    const int32_t* materials = alea_raycast_batch_material_ids(result);
    const double* densities = alea_raycast_batch_densities(result);
    const int32_t* enter_surfaces = alea_raycast_batch_enter_surface_ids(result);
    const int32_t* exit_surfaces = alea_raycast_batch_exit_surface_ids(result);
    const uint8_t* resolution_flags = alea_raycast_batch_resolution_flags(result);
    const uint64_t* path_offsets = include_paths
        ? alea_raycast_batch_segment_path_offsets(result) : NULL;
    const int32_t* path_cells = include_paths
        ? alea_raycast_batch_path_cell_ids(result) : NULL;
    const int32_t* path_materials = include_paths
        ? alea_raycast_batch_path_material_ids(result) : NULL;
    const int32_t* path_universes = include_paths
        ? alea_raycast_batch_path_universe_ids(result) : NULL;
    const int32_t* path_fills = include_paths
        ? alea_raycast_batch_path_fill_universes(result) : NULL;
    const int32_t* path_depths = include_paths
        ? alea_raycast_batch_path_depths(result) : NULL;
    const uint8_t* path_lattice = include_paths
        ? alea_raycast_batch_path_is_lattice(result) : NULL;
    const double* path_origins = include_paths
        ? alea_raycast_batch_path_lattice_origins_xyz(result) : NULL;
    const uint64_t* path_keys = include_paths
        ? alea_raycast_batch_path_occurrence_keys(result) : NULL;

    PyObject *batch = PyList_New(n);
    if (!batch) { alea_raycast_batch_result_destroy(result); return NULL; }
    for (Py_ssize_t i = 0; i < n; ++i) {
        size_t count = (size_t)(ray_offsets[i + 1] - ray_offsets[i]);
        PyObject *segments = PyList_New(count);
        if (!segments) { Py_DECREF(batch); batch = NULL; break; }
        for (size_t j = 0; j < count; ++j) {
            size_t segment_index = (size_t)ray_offsets[i] + j;
            PyObject *item = Py_BuildValue(
                "{s:d,s:d,s:i,s:i,s:d,s:i,s:i,s:i}",
                "t_enter", enters[segment_index],
                "t_exit", exits[segment_index],
                "cell_id", cells[segment_index],
                "material_id", materials[segment_index],
                "density", densities[segment_index],
                "enter_surface_id", enter_surfaces[segment_index],
                "exit_surface_id", exit_surfaces[segment_index],
                "resolution_flags", resolution_flags[segment_index]);
            if (item && include_paths) {
                size_t path_count = (size_t)(
                    path_offsets[segment_index + 1] - path_offsets[segment_index]);
                PyObject* path = PyList_New(path_count);
                if (!path) { Py_DECREF(item); item = NULL; }
                for (size_t k = 0; item && k < path_count; ++k) {
                    size_t path_index = (size_t)path_offsets[segment_index] + k;
                    PyObject* entry = Py_BuildValue(
                        "{s:i,s:i,s:i,s:i,s:i,s:i,s:(ddd),s:K}",
                        "cell_id", path_cells[path_index],
                        "material_id", path_materials[path_index],
                        "universe_id", path_universes[path_index],
                        "fill_universe", path_fills[path_index],
                        "depth", path_depths[path_index],
                        "is_lattice", (int)path_lattice[path_index],
                        "lattice_origin", path_origins[3 * path_index],
                        path_origins[3 * path_index + 1],
                        path_origins[3 * path_index + 2],
                        "occurrence_key", (unsigned long long)path_keys[path_index]);
                    if (!entry) { Py_DECREF(item); item = NULL; break; }
                    PyList_SET_ITEM(path, k, entry);
                }
                if (item && PyDict_SetItemString(item, "path", path) != 0) {
                    Py_DECREF(item); item = NULL;
                }
                Py_XDECREF(path);
            }
            if (!item) {
                Py_XDECREF(item);
                Py_DECREF(segments); Py_DECREF(batch); batch = NULL; break;
            }
            PyList_SET_ITEM(segments, j, item);
        }
        if (!batch) break;
        PyList_SET_ITEM(batch, i, segments);
    }
    alea_raycast_batch_result_destroy(result);
    return batch;
}

/* Copy one compact native field into a NumPy array.  The compact batch result
 * remains owned by libalea and is released immediately after conversion, so
 * Python receives independently owned, contiguous buffers without creating a
 * Python object for every ray segment. */
static PyObject* compact_array_copy(const void* values, size_t count, int typenum) {
    if (count > (size_t)NPY_MAX_INTP) {
        PyErr_SetString(PyExc_OverflowError, "compact ray result is too large for NumPy");
        return NULL;
    }
    npy_intp dims[1] = {(npy_intp)count};
    PyArrayObject* array = (PyArrayObject*)PyArray_SimpleNew(1, dims, typenum);
    if (!array) return NULL;
    if (values && count > 0)
        memcpy(PyArray_DATA(array), values, PyArray_NBYTES(array));
    return (PyObject*)array;
}

static int compact_dict_add_array(PyObject* out, const char* key,
                                  const void* values, size_t count, int typenum) {
    PyObject* array = compact_array_copy(values, count, typenum);
    if (!array) return -1;
    int rc = PyDict_SetItemString(out, key, array);
    Py_DECREF(array);
    return rc;
}

static int compact_dict_add_size(PyObject* out, const char* key, size_t value) {
    PyObject* number = PyLong_FromSize_t(value);
    if (!number) return -1;
    int rc = PyDict_SetItemString(out, key, number);
    Py_DECREF(number);
    return rc;
}

/* Compact arbitrary-plane ray slice.  This is deliberately separate from the
 * legacy raycast_hier_fast_batch() method: it exposes CSR NumPy fields rather
 * than list[ray][dict], and libalea performs row generation and U clipping. */
static PyObject* PyAleaSystem_trace_ray_slice_compact(PyAleaSystemObject* self,
                                                          PyObject* args,
                                                          PyObject* kwds) {
    PyObject *origin_obj, *normal_obj, *up_obj;
    double u_min, u_max, v_min, v_max;
    Py_ssize_t row_count;
    int projected_depth = -1;
    int include_paths = 0;
    unsigned long long max_segments = 0, max_path_entries = 0, max_output_bytes = 0;
    static char* kwlist[] = {
        "origin", "normal", "up", "u_min", "u_max", "v_min", "v_max", "row_count",
        "projected_depth", "include_paths", "max_segments", "max_path_entries",
        "max_output_bytes", NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOddddn|ipKKK", kwlist,
                                     &origin_obj, &normal_obj, &up_obj,
                                     &u_min, &u_max, &v_min, &v_max, &row_count,
                                     &projected_depth, &include_paths,
                                     &max_segments, &max_path_entries,
                                     &max_output_bytes)) return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    if (row_count < 1) {
        PyErr_SetString(PyExc_ValueError, "row_count must be positive");
        return NULL;
    }

    double ox, oy, oz, nx, ny, nz, ux, uy, uz;
    if (!PyArg_ParseTuple(origin_obj, "ddd", &ox, &oy, &oz) ||
        !PyArg_ParseTuple(normal_obj, "ddd", &nx, &ny, &nz) ||
        !PyArg_ParseTuple(up_obj, "ddd", &ux, &uy, &uz)) return NULL;

    alea_slice_view_t view;
    alea_slice_view_init(&view, ox, oy, oz, nx, ny, nz, ux, uy, uz,
                         u_min, u_max, v_min, v_max);
    alea_raycast_batch_options_t options = {
        .struct_size = sizeof(options),
        .fields = ALEA_RAY_BATCH_MATERIAL | ALEA_RAY_BATCH_DENSITY |
                  ALEA_RAY_BATCH_SURFACES | ALEA_RAY_BATCH_RESOLUTION_FLAGS |
                  ALEA_RAY_BATCH_PROJECTED_OWNER |
                  (include_paths ? ALEA_RAY_BATCH_FULL_PATHS : 0),
        .projected_depth = projected_depth,
        .max_segments = (uint64_t)max_segments,
        .max_path_entries = (uint64_t)max_path_entries,
        .max_output_bytes = (uint64_t)max_output_bytes,
    };
    alea_raycast_batch_result_t* result = alea_raycast_batch_result_create();
    if (!result) return PyErr_NoMemory();

    int rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    rc = alea_trace_ray_slice_compact(self->sys, &view, (size_t)row_count,
                                      &options, result);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) {
        alea_raycast_batch_result_destroy(result);
        return NULL;
    }
    if (rc != 0) {
        const char* detail = alea_error();
        alea_raycast_batch_result_destroy(result);
        PyErr_Format(PyExc_RuntimeError, "compact ray slice failed: %s",
                     detail ? detail : "unknown libalea error");
        return NULL;
    }

    size_t rays = alea_raycast_batch_ray_count(result);
    size_t segments = alea_raycast_batch_segment_count(result);
    size_t path_entries = alea_raycast_batch_path_entry_count(result);
    PyObject* out = PyDict_New();
    if (!out) goto failed;
    if (compact_dict_add_array(out, "row_offsets",
            alea_raycast_batch_ray_offsets(result), rays + 1, NPY_UINT64) != 0 ||
        compact_dict_add_array(out, "u_enter",
            alea_raycast_batch_t_enter(result), segments, NPY_DOUBLE) != 0 ||
        compact_dict_add_array(out, "u_exit",
            alea_raycast_batch_t_exit(result), segments, NPY_DOUBLE) != 0 ||
        compact_dict_add_array(out, "cell_ids",
            alea_raycast_batch_cell_ids(result), segments, NPY_INT32) != 0 ||
        compact_dict_add_array(out, "material_ids",
            alea_raycast_batch_material_ids(result), segments, NPY_INT32) != 0 ||
        compact_dict_add_array(out, "densities",
            alea_raycast_batch_densities(result), segments, NPY_DOUBLE) != 0 ||
        compact_dict_add_array(out, "enter_surface_ids",
            alea_raycast_batch_enter_surface_ids(result), segments, NPY_INT32) != 0 ||
        compact_dict_add_array(out, "exit_surface_ids",
            alea_raycast_batch_exit_surface_ids(result), segments, NPY_INT32) != 0 ||
        compact_dict_add_array(out, "resolution_flags",
            alea_raycast_batch_resolution_flags(result), segments, NPY_UBYTE) != 0 ||
        compact_dict_add_array(out, "projected_cell_ids",
            alea_raycast_batch_projected_cell_ids(result), segments, NPY_INT32) != 0 ||
        compact_dict_add_array(out, "projected_material_ids",
            alea_raycast_batch_projected_material_ids(result), segments, NPY_INT32) != 0 ||
        compact_dict_add_array(out, "projected_universe_ids",
            alea_raycast_batch_projected_universe_ids(result), segments, NPY_INT32) != 0 ||
        compact_dict_add_array(out, "projected_fill_universes",
            alea_raycast_batch_projected_fill_universes(result), segments, NPY_INT32) != 0 ||
        compact_dict_add_array(out, "projected_depths",
            alea_raycast_batch_projected_depths(result), segments, NPY_INT32) != 0 ||
        compact_dict_add_array(out, "projected_is_lattice",
            alea_raycast_batch_projected_is_lattice(result), segments, NPY_UBYTE) != 0 ||
        compact_dict_add_array(out, "projected_occurrence_keys",
            alea_raycast_batch_projected_occurrence_keys(result), segments, NPY_UINT64) != 0)
        goto failed;
    if (include_paths &&
        (compact_dict_add_array(out, "segment_path_offsets",
             alea_raycast_batch_segment_path_offsets(result), segments + 1, NPY_UINT64) != 0 ||
         compact_dict_add_array(out, "path_cell_ids",
             alea_raycast_batch_path_cell_ids(result), path_entries, NPY_INT32) != 0 ||
         compact_dict_add_array(out, "path_material_ids",
             alea_raycast_batch_path_material_ids(result), path_entries, NPY_INT32) != 0 ||
         compact_dict_add_array(out, "path_universe_ids",
             alea_raycast_batch_path_universe_ids(result), path_entries, NPY_INT32) != 0 ||
         compact_dict_add_array(out, "path_fill_universes",
             alea_raycast_batch_path_fill_universes(result), path_entries, NPY_INT32) != 0 ||
         compact_dict_add_array(out, "path_depths",
             alea_raycast_batch_path_depths(result), path_entries, NPY_INT32) != 0 ||
         compact_dict_add_array(out, "path_is_lattice",
             alea_raycast_batch_path_is_lattice(result), path_entries, NPY_UBYTE) != 0 ||
         compact_dict_add_array(out, "path_lattice_origins_xyz",
             alea_raycast_batch_path_lattice_origins_xyz(result), path_entries * 3, NPY_DOUBLE) != 0 ||
         compact_dict_add_array(out, "path_occurrence_keys",
             alea_raycast_batch_path_occurrence_keys(result), path_entries, NPY_UINT64) != 0))
        goto failed;
    if (compact_dict_add_size(out, "row_count", rays) != 0 ||
        compact_dict_add_size(out, "segment_count", segments) != 0 ||
        compact_dict_add_size(out, "path_entry_count", path_entries) != 0)
        goto failed;
    alea_raycast_batch_result_destroy(result);
    return out;

failed:
    Py_XDECREF(out);
    alea_raycast_batch_result_destroy(result);
    return NULL;
}

/* Transfer a malloc-owned raster buffer into NumPy without a full-grid copy.
 * This mirrors the grid-query ownership model while keeping this raycast
 * binding independent of the compact-result conversion helpers above. */
static void free_ray_slice_raster_array_data(PyObject* capsule) {
    free(PyCapsule_GetPointer(capsule, "pyalea.ray_slice_raster_array"));
}

static PyObject* ray_slice_raster_array_from_owned_data(void* values,
                                                         Py_ssize_t count,
                                                         int typenum) {
    npy_intp dims[1] = {count};
    PyObject* array = PyArray_SimpleNewFromData(1, dims, typenum, values);
    if (!array) {
        free(values);
        return NULL;
    }
    PyObject* capsule = PyCapsule_New(values, "pyalea.ray_slice_raster_array",
                                      free_ray_slice_raster_array_data);
    if (!capsule) {
        Py_DECREF(array);
        free(values);
        return NULL;
    }
    if (PyArray_SetBaseObject((PyArrayObject*)array, capsule) < 0) {
        /* PyArray_SetBaseObject steals the capsule reference even on error. */
        Py_DECREF(array);
        return NULL;
    }
    return array;
}

/* Fused native trace and raster path.  Every successful field allocation is
 * transferred independently to NumPy, so retained renders never alias. */
static PyObject* PyAleaSystem_trace_ray_slice_grid(PyAleaSystemObject* self,
                                                       PyObject* args,
                                                       PyObject* kwds) {
    PyObject *origin_obj, *normal_obj, *up_obj;
    double u_min, u_max, v_min, v_max;
    Py_ssize_t nu, nv;
    int projected_depth = -1;
    unsigned long long max_segments = 0, max_output_bytes = 0;
    static char* kwlist[] = {
        "origin", "normal", "up", "u_min", "u_max", "v_min", "v_max", "nu", "nv",
        "projected_depth", "max_segments", "max_output_bytes", NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOddddnn|iKK", kwlist,
                                     &origin_obj, &normal_obj, &up_obj,
                                     &u_min, &u_max, &v_min, &v_max, &nu, &nv,
                                     &projected_depth, &max_segments, &max_output_bytes)) return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    if (nu < 1 || nv < 1 || (size_t)nu > SIZE_MAX / (size_t)nv ||
        (size_t)nu * (size_t)nv > (size_t)PY_SSIZE_T_MAX) {
        PyErr_SetString(PyExc_OverflowError, "ray-slice raster dimensions are invalid or too large");
        return NULL;
    }
    const size_t pixels = (size_t)nu * (size_t)nv;
    if (pixels > SIZE_MAX / (4 * sizeof(int32_t) + sizeof(double) + sizeof(uint8_t))) {
        PyErr_SetString(PyExc_OverflowError, "ray-slice raster is too large");
        return NULL;
    }
    double ox, oy, oz, nx, ny, nz, ux, uy, uz;
    if (!PyArg_ParseTuple(origin_obj, "ddd", &ox, &oy, &oz) ||
        !PyArg_ParseTuple(normal_obj, "ddd", &nx, &ny, &nz) ||
        !PyArg_ParseTuple(up_obj, "ddd", &ux, &uy, &uz)) return NULL;

    int32_t* cell_ids = malloc(pixels * sizeof(*cell_ids));
    int32_t* material_ids = malloc(pixels * sizeof(*material_ids));
    int32_t* universe_ids = malloc(pixels * sizeof(*universe_ids));
    int32_t* fill_universe_ids = malloc(pixels * sizeof(*fill_universe_ids));
    double* densities = malloc(pixels * sizeof(*densities));
    uint8_t* flags = malloc(pixels * sizeof(*flags));
    if (!cell_ids || !material_ids || !universe_ids || !fill_universe_ids || !densities || !flags) {
        free(cell_ids); free(material_ids); free(universe_ids); free(fill_universe_ids);
        free(densities); free(flags);
        return PyErr_NoMemory();
    }

    alea_slice_view_t view;
    alea_slice_view_init(&view, ox, oy, oz, nx, ny, nz, ux, uy, uz,
                         u_min, u_max, v_min, v_max);
    alea_slice_raster_t raster;
    alea_slice_raster_init(&raster);
    raster.nu = (size_t)nu;
    raster.nv = (size_t)nv;
    raster.fields = ALEA_SLICE_RASTER_CELL_ID | ALEA_SLICE_RASTER_MATERIAL_ID |
                    ALEA_SLICE_RASTER_UNIVERSE_ID | ALEA_SLICE_RASTER_FILL_UNIVERSE |
                    ALEA_SLICE_RASTER_DENSITY | ALEA_SLICE_RASTER_RESOLUTION_FLAGS;
    raster.cell_ids = cell_ids; raster.material_ids = material_ids;
    raster.universe_ids = universe_ids; raster.fill_universe_ids = fill_universe_ids;
    raster.densities = densities; raster.resolution_flags = flags;
    alea_slice_raster_options_t options;
    alea_slice_raster_options_init(&options);
    options.projected_depth = projected_depth;
    options.max_segments = (uint64_t)max_segments;
    options.max_trace_output_bytes = (uint64_t)max_output_bytes;

    int rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    rc = alea_trace_ray_slice_raster(self->sys, &view, &options, &raster);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint) || rc != 0) {
        const char* detail = rc != 0 ? alea_error() : NULL;
        free(cell_ids); free(material_ids); free(universe_ids); free(fill_universe_ids);
        free(densities); free(flags);
        if (rc != 0) PyErr_Format(PyExc_RuntimeError, "ray slice raster failed: %s",
                                  detail ? detail : "unknown libalea error");
        return NULL;
    }
    /* Keep the viewer's shared error-grid contract: undefined fill is code 3. */
    for (size_t i = 0; i < pixels; i++) flags[i] = (flags[i] & 0x1u) ? 3u : 0u;

    PyObject* out = PyDict_New();
    PyObject *cells = NULL, *materials = NULL, *universes = NULL, *fills = NULL,
             *density_values = NULL, *errors = NULL;
    if (!out) {
        free(cell_ids); free(material_ids); free(universe_ids); free(fill_universe_ids);
        free(densities); free(flags);
        return NULL;
    }
    if (!(cells = ray_slice_raster_array_from_owned_data(cell_ids, (Py_ssize_t)pixels, NPY_INT32)) ||
        !(materials = ray_slice_raster_array_from_owned_data(material_ids, (Py_ssize_t)pixels, NPY_INT32)) ||
        !(universes = ray_slice_raster_array_from_owned_data(universe_ids, (Py_ssize_t)pixels, NPY_INT32)) ||
        !(fills = ray_slice_raster_array_from_owned_data(fill_universe_ids, (Py_ssize_t)pixels, NPY_INT32)) ||
        !(density_values = ray_slice_raster_array_from_owned_data(densities, (Py_ssize_t)pixels, NPY_DOUBLE)) ||
        !(errors = ray_slice_raster_array_from_owned_data(flags, (Py_ssize_t)pixels, NPY_UBYTE))) {
        Py_XDECREF(out); Py_XDECREF(cells); Py_XDECREF(materials); Py_XDECREF(universes);
        Py_XDECREF(fills); Py_XDECREF(density_values); Py_XDECREF(errors);
        return NULL;
    }
    cell_ids = material_ids = universe_ids = fill_universe_ids = NULL;
    densities = NULL; flags = NULL;
    if (PyDict_SetItemString(out, "cell_ids", cells) ||
        PyDict_SetItemString(out, "material_ids", materials) ||
        PyDict_SetItemString(out, "universe_ids", universes) ||
        PyDict_SetItemString(out, "fill_universe_ids", fills) ||
        PyDict_SetItemString(out, "densities", density_values) ||
        PyDict_SetItemString(out, "errors", errors)) {
        Py_DECREF(out); out = NULL;
    }
    if (out) {
        PyObject* metadata = Py_BuildValue("{s:n,s:n,s:d,s:d,s:d,s:d}",
                                           "nu", nu, "nv", nv,
                                           "u_min", u_min, "u_max", u_max,
                                           "v_min", v_min, "v_max", v_max);
        if (!metadata || PyDict_Update(out, metadata) != 0) {
            Py_DECREF(out);
            out = NULL;
        }
        Py_XDECREF(metadata);
    }
    Py_DECREF(cells); Py_DECREF(materials); Py_DECREF(universes);
    Py_DECREF(fills); Py_DECREF(density_values); Py_DECREF(errors);
    return out;
}

/* Parse an (a, b) pair from any sequence, reporting the caller's field name. */
static int parse_double_pair(PyObject* obj, const char* what,
                             double* first, double* second) {
    PyObject* tuple = PySequence_Tuple(obj);
    int ok = tuple && PyArg_ParseTuple(tuple, "dd", first, second);
    Py_XDECREF(tuple);
    if (!ok) {
        PyErr_Clear();
        PyErr_Format(PyExc_TypeError, "%s must be a sequence of two floats", what);
        return -1;
    }
    return 0;
}

/* Adaptive base-row provenance marks refined rows with SIZE_MAX.  Publish it as
 * signed -1 so NumPy consumers test one sentinel rather than a platform width. */
static int compact_dict_add_base_indices(PyObject* out, const char* key,
                                         const size_t* values, size_t count) {
    if (count > (size_t)NPY_MAX_INTP) {
        PyErr_SetString(PyExc_OverflowError, "compact ray result is too large for NumPy");
        return -1;
    }
    npy_intp dims[1] = {(npy_intp)count};
    PyArrayObject* array = (PyArrayObject*)PyArray_SimpleNew(1, dims, NPY_INT64);
    if (!array) return -1;
    int64_t* data = (int64_t*)PyArray_DATA(array);
    for (size_t i = 0; i < count; i++)
        data[i] = (values && values[i] != SIZE_MAX) ? (int64_t)values[i] : -1;
    int rc = PyDict_SetItemString(out, key, (PyObject*)array);
    Py_DECREF(array);
    return rc;
}

static int compact_dict_add_string(PyObject* out, const char* key, const char* value) {
    PyObject* text = PyUnicode_FromString(value);
    if (!text) return -1;
    int rc = PyDict_SetItemString(out, key, text);
    Py_DECREF(text);
    return rc;
}

static const char* slice_refinement_status_name(alea_ray_slice_refinement_status_t status) {
    switch (status) {
        case ALEA_RAY_SLICE_REFINEMENT_NOT_REQUESTED: return "not_requested";
        case ALEA_RAY_SLICE_REFINEMENT_CONVERGED:     return "converged";
        case ALEA_RAY_SLICE_REFINEMENT_MAX_DEPTH:     return "max_depth";
        case ALEA_RAY_SLICE_REFINEMENT_MAX_ROWS:      return "max_rows";
        case ALEA_RAY_SLICE_REFINEMENT_MIN_SPACING:   return "min_spacing";
        default:                                      return "unknown";
    }
}

/* Native fast forward/reverse validation.  Like the compact slice method, this
 * copies CSR buffers directly into NumPy arrays and keeps diagnostics separate
 * from the normal forward render result.
 *
 * The bidirectional check is a trace-consistency signal only; only the coverage
 * check classifies gaps and overlaps, because both directions can agree on the
 * same incomplete ownership interpretation. */
static PyObject* PyAleaSystem_validate_ray_slice_compact(PyAleaSystemObject* self,
                                                             PyObject* args,
                                                             PyObject* kwds) {
    PyObject *origin_obj, *normal_obj, *up_obj;
    double u_min, u_max, v_min, v_max;
    Py_ssize_t row_count;
    int projected_depth = -1, include_agreements = 0, include_provenance = 0;
    unsigned long long max_trace_intervals = 0, max_path_entries = 0;
    unsigned long long max_output_intervals = 0, max_output_bytes = 0;
    Py_ssize_t cache_width = 0;
    int bidirectional = 1, coverage = 0;
    PyObject* coverage_domain_obj = NULL;
    int coverage_unowned_is_exterior = 0, coverage_report_exterior = 0;
    unsigned long long max_coverage_owners = 0, max_coverage_rows = 0;
    unsigned int refinement_depth = 0, refine_signals = 0;
    double refine_displacement = 0.0, min_transverse_spacing = 0.0;
    unsigned long long refine_crossing_density = 0;
    static char* kwlist[] = {
        "origin", "normal", "up", "u_min", "u_max", "v_min", "v_max", "row_count",
        "projected_depth", "include_agreements", "max_trace_intervals",
        "max_path_entries", "max_output_intervals", "max_output_bytes",
        "include_provenance", "cache_width", "bidirectional", "coverage",
        "coverage_domain", "coverage_unowned_is_exterior",
        "coverage_report_exterior", "max_coverage_owners", "max_coverage_rows",
        "refinement_depth", "refine_signals", "refine_displacement",
        "refine_crossing_density", "min_transverse_spacing", NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOddddn|ipKKKKpnppOppKKIIdKd", kwlist,
                                     &origin_obj, &normal_obj, &up_obj,
                                     &u_min, &u_max, &v_min, &v_max, &row_count,
                                     &projected_depth, &include_agreements,
                                     &max_trace_intervals, &max_path_entries,
                                     &max_output_intervals, &max_output_bytes,
                                     &include_provenance, &cache_width,
                                     &bidirectional, &coverage, &coverage_domain_obj,
                                     &coverage_unowned_is_exterior,
                                     &coverage_report_exterior, &max_coverage_owners,
                                     &max_coverage_rows, &refinement_depth,
                                     &refine_signals, &refine_displacement,
                                     &refine_crossing_density,
                                     &min_transverse_spacing)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }
    if (row_count < 1) { PyErr_SetString(PyExc_ValueError, "row_count must be positive"); return NULL; }
    if (!bidirectional && !coverage) {
        PyErr_SetString(PyExc_ValueError,
                        "validation requires the bidirectional or the coverage check");
        return NULL;
    }
    double coverage_domain_u_min = 0.0, coverage_domain_u_max = 0.0;
    int has_coverage_domain = coverage_domain_obj && coverage_domain_obj != Py_None;
    if (has_coverage_domain &&
        parse_double_pair(coverage_domain_obj, "coverage_domain",
                          &coverage_domain_u_min, &coverage_domain_u_max) < 0) return NULL;
    if (cache_width < 0 || (include_provenance &&
        (row_count > INT_MAX || (cache_width && cache_width > INT_MAX)))) {
        PyErr_SetString(PyExc_ValueError, "directional cache dimensions are out of range");
        return NULL;
    }

    double ox, oy, oz, nx, ny, nz, ux, uy, uz;
    if (!PyArg_ParseTuple(origin_obj, "ddd", &ox, &oy, &oz) ||
        !PyArg_ParseTuple(normal_obj, "ddd", &nx, &ny, &nz) ||
        !PyArg_ParseTuple(up_obj, "ddd", &ux, &uy, &uz)) return NULL;
    alea_slice_view_t view;
    alea_slice_view_init(&view, ox, oy, oz, nx, ny, nz, ux, uy, uz,
                         u_min, u_max, v_min, v_max);
    alea_raycast_batch_options_t render_options = {
        .struct_size = sizeof(render_options),
        .fields = ALEA_RAY_BATCH_MATERIAL | ALEA_RAY_BATCH_DENSITY |
                  ALEA_RAY_BATCH_SURFACES | ALEA_RAY_BATCH_RESOLUTION_FLAGS |
                  ALEA_RAY_BATCH_PROJECTED_OWNER,
        .projected_depth = projected_depth,
        .max_segments = (uint64_t)max_trace_intervals,
        .max_path_entries = (uint64_t)max_path_entries,
        .max_output_bytes = (uint64_t)max_output_bytes,
    };
    alea_ray_slice_validation_options_t validation_options;
    alea_ray_slice_validation_options_init(&validation_options);
    validation_options.checks =
        (bidirectional ? (uint32_t)ALEA_RAY_SLICE_VALIDATE_FAST_BIDIRECTIONAL : 0u) |
        (coverage ? (uint32_t)ALEA_RAY_SLICE_VALIDATE_COVERAGE : 0u);
    validation_options.flags = include_agreements ? ALEA_RAY_SLICE_VALIDATION_INCLUDE_AGREEMENTS : 0;
    validation_options.projected_depth = projected_depth;
    validation_options.max_trace_intervals = (uint64_t)max_trace_intervals;
    validation_options.max_path_entries = (uint64_t)max_path_entries;
    validation_options.max_output_intervals = (uint64_t)max_output_intervals;
    validation_options.max_output_bytes = (uint64_t)max_output_bytes;
    /* libalea requires exactly one unowned-space policy for the coverage check;
     * it rejects any other combination with a specific message. */
    validation_options.coverage_flags =
        (has_coverage_domain ? (uint32_t)ALEA_RAY_SLICE_COVERAGE_HAS_DOMAIN : 0u) |
        (coverage_unowned_is_exterior ? (uint32_t)ALEA_RAY_SLICE_COVERAGE_UNOWNED_IS_EXTERIOR : 0u) |
        (coverage_report_exterior ? (uint32_t)ALEA_RAY_SLICE_COVERAGE_REPORT_EXTERIOR : 0u);
    validation_options.coverage_domain_u_min = coverage_domain_u_min;
    validation_options.coverage_domain_u_max = coverage_domain_u_max;
    validation_options.max_coverage_owners = (uint64_t)max_coverage_owners;
    validation_options.max_coverage_rows = (uint64_t)max_coverage_rows;
    validation_options.coverage_max_refinement_depth = refinement_depth;
    validation_options.coverage_refine_signals = refine_signals;
    validation_options.coverage_endpoint_displacement = refine_displacement;
    validation_options.coverage_crossing_density = (uint64_t)refine_crossing_density;
    validation_options.coverage_min_transverse_spacing = min_transverse_spacing;
    alea_raycast_batch_result_t* forward = alea_raycast_batch_result_create();
    alea_ray_slice_validation_result_t* validation = alea_ray_slice_validation_result_create();
    alea_slice_directional_trace_cache_t* cache = NULL;
    if (!forward || !validation) { alea_raycast_batch_result_destroy(forward); alea_ray_slice_validation_result_destroy(validation); return PyErr_NoMemory(); }
    if (include_provenance) {
        int width = cache_width ? (int)cache_width : (int)row_count;
        cache = alea_slice_directional_trace_cache_create(self->sys, &view, width,
                                                          (int)row_count);
        if (!cache) {
            const char* detail = alea_error();
            alea_raycast_batch_result_destroy(forward);
            alea_ray_slice_validation_result_destroy(validation);
            PyErr_Format(PyExc_RuntimeError, "directional ray cache failed: %s",
                         detail ? detail : "unknown libalea error");
            return NULL;
        }
    }
    int rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    if (cache) {
        /* Cached ownership traces stay owned by libalea's directional cache,
         * whereas pyAlea returns an independently owned compact render.
         * Validate from the cache, then materialize that public render once. */
        rc = alea_validate_ray_slice_compact_with_directional_cache(
            self->sys, &view, (size_t)row_count, &validation_options,
            &render_options, NULL, cache, validation);
        if (rc == 0)
            rc = alea_trace_ray_slice_compact(self->sys, &view, (size_t)row_count,
                                              &render_options, forward);
    } else {
        rc = alea_validate_ray_slice_compact(self->sys, &view, (size_t)row_count,
                                             &validation_options, &render_options,
                                             forward, validation);
    }
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) goto failed;
    if (rc != 0) {
        const char* detail = alea_error();
        PyErr_Format(PyExc_RuntimeError, "compact ray-slice validation failed: %s", detail ? detail : "unknown libalea error");
        goto failed;
    }
    size_t rays = alea_raycast_batch_ray_count(forward);
    size_t segments = alea_raycast_batch_segment_count(forward);
    size_t intervals = alea_ray_slice_validation_interval_count(validation);
    /* Adaptive refinement publishes rows between the requested ones, so the
     * validation row count is independent of the render's ray count. */
    size_t validation_rows = alea_ray_slice_validation_row_count(validation);
    uint32_t validation_fields = alea_ray_slice_validation_fields(validation);
    uint32_t executed_trace_mask =
        alea_ray_slice_validation_executed_trace_mask(validation);
    uint32_t reused_trace_mask =
        alea_ray_slice_validation_reused_trace_mask(validation);
    PyObject *out = PyDict_New(), *render = PyDict_New();
    if (!out || !render) { Py_XDECREF(out); Py_XDECREF(render); goto failed; }
    if (compact_dict_add_array(render, "row_offsets", alea_raycast_batch_ray_offsets(forward), rays + 1, NPY_UINT64) ||
        compact_dict_add_array(render, "u_enter", alea_raycast_batch_t_enter(forward), segments, NPY_DOUBLE) ||
        compact_dict_add_array(render, "u_exit", alea_raycast_batch_t_exit(forward), segments, NPY_DOUBLE) ||
        compact_dict_add_array(render, "cell_ids", alea_raycast_batch_cell_ids(forward), segments, NPY_INT32) ||
        compact_dict_add_array(render, "material_ids", alea_raycast_batch_material_ids(forward), segments, NPY_INT32) ||
        compact_dict_add_array(render, "densities", alea_raycast_batch_densities(forward), segments, NPY_DOUBLE) ||
        compact_dict_add_array(render, "enter_surface_ids", alea_raycast_batch_enter_surface_ids(forward), segments, NPY_INT32) ||
        compact_dict_add_array(render, "exit_surface_ids", alea_raycast_batch_exit_surface_ids(forward), segments, NPY_INT32) ||
        compact_dict_add_array(render, "resolution_flags", alea_raycast_batch_resolution_flags(forward), segments, NPY_UBYTE) ||
        compact_dict_add_array(render, "projected_cell_ids", alea_raycast_batch_projected_cell_ids(forward), segments, NPY_INT32) ||
        compact_dict_add_array(render, "projected_material_ids", alea_raycast_batch_projected_material_ids(forward), segments, NPY_INT32) ||
        compact_dict_add_array(render, "projected_universe_ids", alea_raycast_batch_projected_universe_ids(forward), segments, NPY_INT32) ||
        compact_dict_add_array(render, "projected_fill_universes", alea_raycast_batch_projected_fill_universes(forward), segments, NPY_INT32) ||
        compact_dict_add_array(render, "projected_depths", alea_raycast_batch_projected_depths(forward), segments, NPY_INT32) ||
        compact_dict_add_array(render, "projected_is_lattice", alea_raycast_batch_projected_is_lattice(forward), segments, NPY_UBYTE) ||
        compact_dict_add_array(render, "projected_occurrence_keys", alea_raycast_batch_projected_occurrence_keys(forward), segments, NPY_UINT64) ||
        compact_dict_add_size(render, "row_count", rays) || compact_dict_add_size(render, "segment_count", segments) ||
        compact_dict_add_array(out, "row_offsets", alea_ray_slice_validation_row_offsets(validation), validation_rows + 1, NPY_UINT64) ||
        compact_dict_add_array(out, "u_enter", alea_ray_slice_validation_u_enter(validation), intervals, NPY_DOUBLE) ||
        compact_dict_add_array(out, "u_exit", alea_ray_slice_validation_u_exit(validation), intervals, NPY_DOUBLE) ||
        compact_dict_add_array(out, "flags", alea_ray_slice_validation_diagnostic_flags(validation), intervals, NPY_UINT32) ||
        compact_dict_add_array(out, "fast_forward_cell_ids", alea_ray_slice_validation_fast_forward_cell_ids(validation), intervals, NPY_INT32) ||
        compact_dict_add_array(out, "fast_reverse_cell_ids", alea_ray_slice_validation_fast_reverse_cell_ids(validation), intervals, NPY_INT32) ||
        compact_dict_add_array(out, "fast_forward_occurrence_keys", alea_ray_slice_validation_fast_forward_occurrence_keys(validation), intervals, NPY_UINT64) ||
        compact_dict_add_array(out, "fast_reverse_occurrence_keys", alea_ray_slice_validation_fast_reverse_occurrence_keys(validation), intervals, NPY_UINT64) ||
        (cache && (compact_dict_add_array(out, "u_enter_forward_surface_ids", alea_ray_slice_validation_u_enter_forward_surface_ids(validation), intervals, NPY_INT32) ||
                   compact_dict_add_array(out, "u_enter_reverse_surface_ids", alea_ray_slice_validation_u_enter_reverse_surface_ids(validation), intervals, NPY_INT32) ||
                   compact_dict_add_array(out, "u_exit_forward_surface_ids", alea_ray_slice_validation_u_exit_forward_surface_ids(validation), intervals, NPY_INT32) ||
                   compact_dict_add_array(out, "u_exit_reverse_surface_ids", alea_ray_slice_validation_u_exit_reverse_surface_ids(validation), intervals, NPY_INT32) ||
                   compact_dict_add_array(out, "u_enter_provenance_flags", alea_ray_slice_validation_u_enter_provenance_flags(validation), intervals, NPY_UINT32) ||
                   compact_dict_add_array(out, "u_exit_provenance_flags", alea_ray_slice_validation_u_exit_provenance_flags(validation), intervals, NPY_UINT32))) ||
        compact_dict_add_size(out, "interval_count", intervals) ||
        compact_dict_add_size(out, "row_count", validation_rows) ||
        compact_dict_add_size(out, "fields", (size_t)validation_fields) ||
        compact_dict_add_size(out, "executed_trace_mask",
                              (size_t)executed_trace_mask) ||
        compact_dict_add_size(out, "reused_trace_mask",
                              (size_t)reused_trace_mask) ||
        PyDict_SetItemString(out, "render", render) != 0) { Py_DECREF(render); Py_DECREF(out); goto failed; }
    if (PyDict_SetItemString(out, "used_directional_cache", cache ? Py_True : Py_False) != 0) {
        Py_DECREF(render); Py_DECREF(out); goto failed;
    }
    Py_DECREF(render);
    alea_ray_slice_refinement_status_t refinement =
        alea_ray_slice_validation_refinement_status(validation);
    if (compact_dict_add_size(out, "refinement_status", (size_t)refinement) ||
        compact_dict_add_string(out, "refinement_status_name",
                                slice_refinement_status_name(refinement))) {
        Py_DECREF(out); goto failed;
    }
    /* Owner counts are the retained prefix for a TRUNCATED interval, not a
     * complete count; row provenance appears only once refinement is requested. */
    if ((validation_fields & ALEA_RAY_SLICE_VALIDATION_FIELD_COVERAGE) &&
        compact_dict_add_array(out, "coverage_owner_counts",
                               alea_ray_slice_validation_coverage_owner_counts(validation),
                               intervals, NPY_UINT32)) {
        Py_DECREF(out); goto failed;
    }
    if ((validation_fields & ALEA_RAY_SLICE_VALIDATION_FIELD_ADAPTIVE_ROWS) &&
        (compact_dict_add_array(out, "row_transverse_coordinates",
                                alea_ray_slice_validation_row_transverse_coordinates(validation),
                                validation_rows, NPY_DOUBLE) ||
         compact_dict_add_array(out, "row_direction_tags",
                                alea_ray_slice_validation_row_direction_tags(validation),
                                validation_rows, NPY_UBYTE) ||
         compact_dict_add_base_indices(out, "row_base_indices",
                                       alea_ray_slice_validation_row_base_indices(validation),
                                       validation_rows))) {
        Py_DECREF(out); goto failed;
    }
    alea_raycast_batch_result_destroy(forward);
    alea_ray_slice_validation_result_destroy(validation);
    alea_slice_directional_trace_cache_destroy(cache);
    return out;
failed:
    alea_raycast_batch_result_destroy(forward);
    alea_ray_slice_validation_result_destroy(validation);
    alea_slice_directional_trace_cache_destroy(cache);
    return NULL;
}

/* ============================================================================
 * PyAleaSystem Methods - Complete-coverage diagnostics
 * ========================================================================= */

static const char* coverage_refinement_status_name(alea_ray_coverage_refinement_status_t status) {
    switch (status) {
        case ALEA_RAY_COVERAGE_REFINEMENT_COMPLETE:    return "complete";
        case ALEA_RAY_COVERAGE_REFINEMENT_MAX_DEPTH:   return "max_depth";
        case ALEA_RAY_COVERAGE_REFINEMENT_MAX_ROWS:    return "max_rows";
        case ALEA_RAY_COVERAGE_REFINEMENT_MIN_SPACING: return "min_spacing";
        default:                                       return "unknown";
    }
}

/* Borrow an optional per-row array of `expected` elements.  Sets *out to NULL
 * when the argument is absent, and returns -1 with an exception on mismatch. */
static int coverage_optional_row_array(PyObject* obj, const char* name, int typenum,
                                       Py_ssize_t expected, PyArrayObject** out) {
    *out = NULL;
    if (!obj || obj == Py_None) return 0;
    PyArrayObject* array = (PyArrayObject*)PyArray_FROM_OTF(obj, typenum, NPY_ARRAY_CARRAY_RO);
    if (!array) return -1;
    if (PyArray_NDIM(array) != 1 || PyArray_DIM(array, 0) != expected) {
        Py_DECREF(array);
        PyErr_Format(PyExc_ValueError, "%s must have shape (n_rays,)", name);
        return -1;
    }
    *out = array;
    return 0;
}

/* Complete ownership coverage for packed rays.  Unlike the tracing methods,
 * which publish the one selected owner per interval, this publishes every
 * concrete owner occurrence claiming each elementary interval — the only way
 * to see gaps and overlaps rather than one self-consistent interpretation. */
static PyObject* PyAleaSystem_ray_coverage_slice(PyAleaSystemObject* self,
                                                   PyObject* args, PyObject* kwds) {
    PyObject *origins_obj, *directions_obj;
    PyObject *tags_obj = NULL, *transverse_obj = NULL, *domain_obj = NULL;
    double t_max;
    int report_exterior = 0;
    unsigned long long max_rows = 0, max_intervals = 0, max_owners = 0, max_output_bytes = 0;
    Py_ssize_t max_refinement_depth = 0, crossing_density = 0;
    unsigned int refine_signals = 0;
    double endpoint_displacement = 0.0, min_transverse_spacing = 0.0;
    static char* kwlist[] = {
        "origins", "directions", "t_max", "direction_tags",
        "transverse_coordinates", "domain", "report_exterior", "max_rows",
        "max_intervals", "max_owners", "max_output_bytes",
        "max_refinement_depth", "refine_signals", "endpoint_displacement",
        "crossing_density", "min_transverse_spacing", NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOd|OOOpKKKKnIdnd", kwlist,
                                     &origins_obj, &directions_obj, &t_max,
                                     &tags_obj, &transverse_obj, &domain_obj,
                                     &report_exterior, &max_rows, &max_intervals,
                                     &max_owners, &max_output_bytes,
                                     &max_refinement_depth, &refine_signals,
                                     &endpoint_displacement, &crossing_density,
                                     &min_transverse_spacing)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }
    if (max_refinement_depth < 0 || crossing_density < 0) {
        PyErr_SetString(PyExc_ValueError, "coverage refinement limits must be non-negative");
        return NULL;
    }
    double domain_t_min = 0.0, domain_t_max = 0.0;
    int has_domain = domain_obj && domain_obj != Py_None;
    if (has_domain &&
        parse_double_pair(domain_obj, "domain", &domain_t_min, &domain_t_max) < 0) return NULL;
    if (ensure_query_acceleration(self) < 0) return NULL;

    PyArrayObject* origins = (PyArrayObject*)PyArray_FROM_OTF(origins_obj, NPY_DOUBLE, NPY_ARRAY_CARRAY_RO);
    PyArrayObject* directions = (PyArrayObject*)PyArray_FROM_OTF(directions_obj, NPY_DOUBLE, NPY_ARRAY_CARRAY_RO);
    if (!origins || !directions) { Py_XDECREF(origins); Py_XDECREF(directions); return NULL; }
    if (PyArray_NDIM(origins) != 2 || PyArray_NDIM(directions) != 2 ||
        PyArray_DIM(origins, 1) != 3 || PyArray_DIM(directions, 1) != 3 ||
        PyArray_DIM(origins, 0) != PyArray_DIM(directions, 0)) {
        Py_DECREF(origins); Py_DECREF(directions);
        PyErr_SetString(PyExc_ValueError, "origins and directions must have shape (n_rays, 3)");
        return NULL;
    }
    Py_ssize_t rows = PyArray_DIM(origins, 0);
    if (rows < 1) {
        Py_DECREF(origins); Py_DECREF(directions);
        PyErr_SetString(PyExc_ValueError, "coverage query requires at least one ray");
        return NULL;
    }
    PyArrayObject *tags = NULL, *transverse = NULL;
    if (coverage_optional_row_array(tags_obj, "direction_tags", NPY_UINT8, rows, &tags) < 0 ||
        coverage_optional_row_array(transverse_obj, "transverse_coordinates",
                                    NPY_DOUBLE, rows, &transverse) < 0) {
        Py_DECREF(origins); Py_DECREF(directions); Py_XDECREF(tags);
        return NULL;
    }

    alea_ray_coverage_slice_options_t options;
    alea_ray_coverage_slice_options_init(&options);
    options.flags = (has_domain ? (uint32_t)ALEA_RAY_COVERAGE_DOMAIN : 0u) |
                    (report_exterior ? (uint32_t)ALEA_RAY_COVERAGE_REPORT_EXTERIOR : 0u);
    options.t_max = t_max;
    options.domain_t_min = domain_t_min;
    options.domain_t_max = domain_t_max;
    options.max_rows = (uint64_t)max_rows;
    options.max_intervals = (uint64_t)max_intervals;
    options.max_owners = (uint64_t)max_owners;
    options.max_output_bytes = (uint64_t)max_output_bytes;
    options.max_refinement_depth = (size_t)max_refinement_depth;
    options.refinement_signals = refine_signals;
    options.min_transverse_spacing = min_transverse_spacing;
    options.endpoint_displacement = endpoint_displacement;
    options.crossing_density = (size_t)crossing_density;

    alea_ray_coverage_slice_result_t* result = alea_ray_coverage_slice_result_create();
    if (!result) {
        Py_DECREF(origins); Py_DECREF(directions); Py_XDECREF(tags); Py_XDECREF(transverse);
        return PyErr_NoMemory();
    }
    const double* origin_data = (const double*)PyArray_DATA(origins);
    const double* direction_data = (const double*)PyArray_DATA(directions);
    const uint8_t* tag_data = tags ? (const uint8_t*)PyArray_DATA(tags) : NULL;
    const double* transverse_data = transverse ? (const double*)PyArray_DATA(transverse) : NULL;
    int rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    rc = alea_ray_coverage_slice_query(self->sys, origin_data, direction_data,
                                       (size_t)rows, tag_data, transverse_data,
                                       &options, result);
    Py_END_ALLOW_THREADS
    Py_DECREF(origins); Py_DECREF(directions); Py_XDECREF(tags); Py_XDECREF(transverse);
    if (restore_sigint(old_sigint) || rc != 0) {
        if (rc != 0) {
            const char* detail = alea_error();
            PyErr_Format(PyExc_RuntimeError, "ray coverage query failed: %s",
                         detail ? detail : "unknown libalea error");
        }
        alea_ray_coverage_slice_result_destroy(result);
        return NULL;
    }

    size_t row_count = alea_ray_coverage_slice_row_count(result);
    size_t interval_count = alea_ray_coverage_slice_interval_count(result);
    size_t owner_count = alea_ray_coverage_slice_owner_count(result);
    alea_ray_coverage_refinement_status_t refinement =
        (alea_ray_coverage_refinement_status_t)alea_ray_coverage_slice_refinement_status(result);
    PyObject* out = PyDict_New();
    if (!out) { alea_ray_coverage_slice_result_destroy(result); return NULL; }
    if (compact_dict_add_array(out, "row_offsets", alea_ray_coverage_slice_row_offsets(result), row_count + 1, NPY_UINTP) ||
        compact_dict_add_array(out, "row_direction_tags", alea_ray_coverage_slice_row_direction_tags(result), row_count, NPY_UBYTE) ||
        compact_dict_add_array(out, "row_transverse_coordinates", alea_ray_coverage_slice_row_transverse_coordinates(result), row_count, NPY_DOUBLE) ||
        compact_dict_add_array(out, "t_enter", alea_ray_coverage_slice_t_enter(result), interval_count, NPY_DOUBLE) ||
        compact_dict_add_array(out, "t_exit", alea_ray_coverage_slice_t_exit(result), interval_count, NPY_DOUBLE) ||
        compact_dict_add_array(out, "kinds", alea_ray_coverage_slice_kinds(result), interval_count, NPY_UBYTE) ||
        compact_dict_add_array(out, "owner_offsets", alea_ray_coverage_slice_owner_offsets(result), interval_count + 1, NPY_UINTP) ||
        compact_dict_add_array(out, "owner_count_lower_bounds", alea_ray_coverage_slice_owner_count_lower_bounds(result), interval_count, NPY_UINTP) ||
        compact_dict_add_array(out, "owner_cell_ids", alea_ray_coverage_slice_owner_cell_ids(result), owner_count, NPY_INT32) ||
        compact_dict_add_array(out, "owner_material_ids", alea_ray_coverage_slice_owner_material_ids(result), owner_count, NPY_INT32) ||
        compact_dict_add_array(out, "owner_universe_ids", alea_ray_coverage_slice_owner_universe_ids(result), owner_count, NPY_INT32) ||
        compact_dict_add_array(out, "owner_fill_universes", alea_ray_coverage_slice_owner_fill_universes(result), owner_count, NPY_INT32) ||
        compact_dict_add_array(out, "owner_depths", alea_ray_coverage_slice_owner_depths(result), owner_count, NPY_INT32) ||
        compact_dict_add_array(out, "owner_occurrence_keys", alea_ray_coverage_slice_owner_occurrence_keys(result), owner_count, NPY_UINT64) ||
        compact_dict_add_array(out, "owner_parent_occurrence_keys", alea_ray_coverage_slice_owner_parent_occurrence_keys(result), owner_count, NPY_UINT64) ||
        compact_dict_add_array(out, "owner_resolution_flags", alea_ray_coverage_slice_owner_resolution_flags(result), owner_count, NPY_UBYTE) ||
        compact_dict_add_size(out, "row_count", row_count) ||
        compact_dict_add_size(out, "interval_count", interval_count) ||
        compact_dict_add_size(out, "owner_count", owner_count) ||
        compact_dict_add_size(out, "refinement_status", (size_t)refinement) ||
        compact_dict_add_string(out, "refinement_status_name",
                                coverage_refinement_status_name(refinement))) {
        Py_DECREF(out);
        alea_ray_coverage_slice_result_destroy(result);
        return NULL;
    }
    alea_ray_coverage_slice_result_destroy(result);
    return out;
}
