// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* This file is part of pyalea_binding.c — do NOT compile separately.
 * It is #included from pyalea_binding.c to keep a single translation unit.
 *
 * Contents: Slice curves (Z/Y/X/arbitrary), grid cell queries,
 *           label positions, grid overlap check, surface label positions.
 */

/* Inline curve geometry into an error-segment dict so Python can render
 * analytical error lines without re-fetching the curve list. */
static void inline_curve_into_dict(PyObject* dict, const alea_curve_t* c) {
    dict_set_new(dict, "curve_type",
                         PyUnicode_FromString(curve_type_to_string(c->type)));
    switch (c->type) {
        case ALEA_CURVE_LINE:
        case ALEA_CURVE_LINE_SEGMENT:
            dict_set_new(dict, "point",
                Py_BuildValue("(dd)", c->data.line.point[0], c->data.line.point[1]));
            dict_set_new(dict, "direction",
                Py_BuildValue("(dd)", c->data.line.direction[0], c->data.line.direction[1]));
            break;
        case ALEA_CURVE_CIRCLE:
        case ALEA_CURVE_ARC:
            dict_set_new(dict, "center",
                Py_BuildValue("(dd)", c->data.circle.center[0], c->data.circle.center[1]));
            dict_set_new(dict, "radius",
                PyFloat_FromDouble(c->data.circle.radius));
            break;
        case ALEA_CURVE_ELLIPSE:
        case ALEA_CURVE_ELLIPSE_ARC:
            dict_set_new(dict, "center",
                Py_BuildValue("(dd)", c->data.ellipse.center[0], c->data.ellipse.center[1]));
            dict_set_new(dict, "semi_a", PyFloat_FromDouble(c->data.ellipse.semi_a));
            dict_set_new(dict, "semi_b", PyFloat_FromDouble(c->data.ellipse.semi_b));
            dict_set_new(dict, "angle", PyFloat_FromDouble(c->data.ellipse.angle));
            break;
        case ALEA_CURVE_PARALLEL_LINES:
            dict_set_new(dict, "point1",
                Py_BuildValue("(dd)", c->data.parallel_lines.point1[0],
                              c->data.parallel_lines.point1[1]));
            dict_set_new(dict, "point2",
                Py_BuildValue("(dd)", c->data.parallel_lines.point2[0],
                              c->data.parallel_lines.point2[1]));
            dict_set_new(dict, "direction",
                Py_BuildValue("(dd)", c->data.parallel_lines.direction[0],
                              c->data.parallel_lines.direction[1]));
            break;
        case ALEA_CURVE_POLYGON: {
            PyObject* verts = PyList_New(c->data.polygon.count);
            for (int j = 0; j < c->data.polygon.count; j++) {
                PyList_SET_ITEM(verts, j, Py_BuildValue("(dd)",
                    c->data.polygon.vertices[j][0],
                    c->data.polygon.vertices[j][1]));
            }
            dict_set_new(dict, "vertices", verts);
            dict_set_new(dict, "closed",
                PyBool_FromLong(c->data.polygon.closed));
            break;
        }
        default:
            break;
    }
}

static PyObject* build_error_components(const int* cell_ids,
                                        const int* secondary_cell_ids,
                                        const uint8_t* coverage,
                                        int nu, int nv) {
    alea_plot_error_component_result_t* comps =
        alea_classify_plot_error_components(cell_ids, secondary_cell_ids,
                                            coverage, nu, nv);
    if (!comps) return NULL;

    PyObject* list = PyList_New(comps->component_count);
    if (!list) {
        alea_plot_error_components_free(comps);
        return NULL;
    }

    for (size_t i = 0; i < comps->component_count; i++) {
        const alea_plot_error_component_t* c = &comps->components[i];
        PyObject* d = PyDict_New();
        const char* kind = "undefined";
        if (c->kind == ALEA_PLOT_ERR_PARTIAL_OVERLAP) kind = "partial_overlap";
        else if (c->kind == ALEA_PLOT_ERR_TOTAL_OVERLAP) kind = "total_overlap";
        dict_set_new(d, "kind", PyUnicode_FromString(kind));
        dict_set_new(d, "primary_cell_id", PyLong_FromLong(c->primary_cell_id));
        dict_set_new(d, "secondary_cell_id", PyLong_FromLong(c->secondary_cell_id));
        dict_set_new(d, "pixel_count", PyLong_FromLong(c->pixel_count));
        dict_set_new(d, "bounds",
            Py_BuildValue("(iiii)", c->min_i, c->min_j, c->max_i, c->max_j));
        dict_set_new(d, "representative_pixel",
            Py_BuildValue("(ii)", c->representative_i, c->representative_j));
        PyList_SET_ITEM(list, i, d);
    }

    alea_plot_error_components_free(comps);
    return list;
}

static PyObject* build_component_list_from_result(
    alea_plot_error_component_result_t* comps) {
    if (!comps) return NULL;
    PyObject* list = PyList_New(comps->component_count);
    if (!list) {
        alea_plot_error_components_free(comps);
        return NULL;
    }
    for (size_t i = 0; i < comps->component_count; i++) {
        const alea_plot_error_component_t* c = &comps->components[i];
        PyObject* d = PyDict_New();
        const char* kind = "undefined";
        if (c->kind == ALEA_PLOT_ERR_PARTIAL_OVERLAP) kind = "partial_overlap";
        else if (c->kind == ALEA_PLOT_ERR_TOTAL_OVERLAP) kind = "total_overlap";
        if (!d ||
            dict_set_new(d, "kind", PyUnicode_FromString(kind)) < 0 ||
            dict_set_new(d, "primary_cell_id", PyLong_FromLong(c->primary_cell_id)) < 0 ||
            dict_set_new(d, "secondary_cell_id", PyLong_FromLong(c->secondary_cell_id)) < 0 ||
            dict_set_new(d, "pixel_count", PyLong_FromLong(c->pixel_count)) < 0 ||
            dict_set_new(d, "bounds",
                                 Py_BuildValue("(iiii)", c->min_i, c->min_j,
                                               c->max_i, c->max_j)) < 0) {
            Py_XDECREF(d);
            Py_DECREF(list);
            alea_plot_error_components_free(comps);
            return NULL;
        }
        if (dict_set_new(d, "representative_pixel",
                                 Py_BuildValue("(ii)", c->representative_i,
                                               c->representative_j)) < 0) {
            Py_DECREF(d);
            Py_DECREF(list);
            alea_plot_error_components_free(comps);
            return NULL;
        }
        PyList_SET_ITEM(list, i, d);
    }
    alea_plot_error_components_free(comps);
    return list;
}

/* Exact local coverage without returning a pixel-sized raster. This keeps
 * viewer deep scans bounded even when the selected region is dense. */
static PyObject* PyAleaSystem_find_local_coverage_components(
    PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    PyObject *origin_obj, *normal_obj, *up_obj;
    double u_min, u_max, v_min, v_max;
    int nu, nv, universe_depth = -1, max_workers = 0;
    Py_ssize_t max_pixels_arg = 0, max_scratch_arg = 0;
    static char* kwlist[] = {
        "origin", "normal", "up", "u_min", "u_max", "v_min", "v_max",
        "nu", "nv", "universe_depth", "max_pixels", "max_scratch_bytes", "max_workers", NULL};
    if (!PyArg_ParseTupleAndKeywords(
            args, kwds, "OOOddddii|inni", kwlist,
            &origin_obj, &normal_obj, &up_obj, &u_min, &u_max, &v_min, &v_max,
            &nu, &nv, &universe_depth, &max_pixels_arg, &max_scratch_arg, &max_workers))
        return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    if (nu <= 0 || nv <= 0 || max_pixels_arg < 0 || max_scratch_arg < 0 || max_workers < 0) {
        PyErr_SetString(PyExc_ValueError,
                        "dimensions, budgets, and worker cap must be non-negative");
        return NULL;
    }
    double ox, oy, oz, nx, ny, nz, ux, uy, uz;
    if (!PyArg_ParseTuple(origin_obj, "ddd", &ox, &oy, &oz) ||
        !PyArg_ParseTuple(normal_obj, "ddd", &nx, &ny, &nz) ||
        !PyArg_ParseTuple(up_obj, "ddd", &ux, &uy, &uz)) {
        PyErr_SetString(PyExc_TypeError,
                        "origin, normal, and up must be (x, y, z) tuples");
        return NULL;
    }
    alea_slice_view_t view;
    alea_slice_view_init(&view, ox, oy, oz, nx, ny, nz, ux, uy, uz,
                         u_min, u_max, v_min, v_max);
    alea_plot_error_component_result_t* components = NULL;
    alea_local_coverage_stats_t stats;
    int rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    rc = alea_find_local_coverage_components(
        self->sys, &view, nu, nv, universe_depth,
        (size_t)max_pixels_arg, (size_t)max_scratch_arg, max_workers,
        &components, &stats);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) {
        alea_plot_error_components_free(components);
        return NULL;
    }
    if (rc == ALEA_LOCAL_COVERAGE_BUDGET_EXCEEDED) {
        PyErr_SetString(PyExc_ValueError,
                        "local coverage request exceeds its pixel or scratch-memory budget");
        return NULL;
    }
    if (rc != 0) {
        alea_plot_error_components_free(components);
        PyErr_SetString(PyExc_RuntimeError, "Local coverage query failed");
        return NULL;
    }
    PyObject* component_list = build_component_list_from_result(components);
    if (!component_list) return NULL;
    PyObject* result = PyDict_New();
    PyObject* point_stats = Py_BuildValue(
        "{s:n,s:n,s:n,s:n,s:n,s:n,s:n,s:n,s:n}",
        "queries", (Py_ssize_t)stats.point_coverage.queries,
        "spatial_queries", (Py_ssize_t)stats.point_coverage.spatial_queries,
        "spatial_multi_early_exit", (Py_ssize_t)stats.point_coverage.spatial_multi_early_exit,
        "recursive_fallbacks", (Py_ssize_t)stats.point_coverage.recursive_fallbacks,
        "lattice_fallbacks", (Py_ssize_t)stats.point_coverage.lattice_fallbacks,
        "truncated_fallbacks", (Py_ssize_t)stats.point_coverage.truncated_fallbacks,
        "candidate_total", (Py_ssize_t)stats.point_coverage.candidate_total,
        "candidate_max", (Py_ssize_t)stats.point_coverage.candidate_max,
        "contains_tests", (Py_ssize_t)stats.point_coverage.contains_tests);
    if (!result || !point_stats ||
        PyDict_SetItemString(result, "components", component_list) < 0 ||
        dict_set_new(result, "pixels", PyLong_FromSize_t(stats.pixels)) < 0 ||
        dict_set_new(result, "scratch_bytes", PyLong_FromSize_t(stats.scratch_bytes)) < 0 ||
        dict_set_new(result, "incomplete_points",
                            PyLong_FromSize_t(stats.incomplete_points)) < 0 ||
        dict_set_new(result, "worker_limit", PyLong_FromLong(stats.worker_limit)) < 0 ||
        PyDict_SetItemString(result, "point_coverage", point_stats) < 0) {
        Py_XDECREF(result); Py_XDECREF(point_stats); Py_DECREF(component_list);
        return NULL;
    }
    Py_DECREF(component_list);
    Py_DECREF(point_stats);
    return result;
}

static int parse_grid_error_mode(const char* error_mode,
                                 int* detect_errors,
                                 int* full_errors) {
    if (!error_mode || strcmp(error_mode, "none") == 0) {
        *detect_errors = 0;
        *full_errors = 0;
        return 0;
    }
    if (strcmp(error_mode, "fast") == 0) {
        *detect_errors = 1;
        *full_errors = 0;
        return 0;
    }
    if (strcmp(error_mode, "full") == 0) {
        *detect_errors = 1;
        *full_errors = 1;
        return 0;
    }
    PyErr_SetString(PyExc_ValueError,
                    "error_mode must be 'none', 'fast', or 'full'");
    return -1;
}

/* Transfer a malloc-owned result buffer into a one-dimensional NumPy array.
 * The array's base capsule frees the C allocation when the last Python view
 * dies, so no intermediate bytes object or full-grid copy is needed. */
static void free_grid_array_data(PyObject* capsule) {
    free(PyCapsule_GetPointer(capsule, "pyalea.grid_array"));
}

static PyObject* grid_array_from_owned_data(void* values, Py_ssize_t count,
                                            int typenum) {
    if (count < 0) {
        free(values);
        PyErr_SetString(PyExc_OverflowError, "grid is too large");
        return NULL;
    }
    npy_intp dims[1] = {count};
    PyObject* array = PyArray_SimpleNewFromData(1, dims, typenum, values);
    if (!array) {
        free(values);
        return NULL;
    }
    PyObject* capsule = PyCapsule_New(values, "pyalea.grid_array",
                                      free_grid_array_data);
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

static int copy_int_values(PyObject* source, int* target, Py_ssize_t count,
                           const char* name) {
    if (PyList_Check(source)) {
        if (PyList_Size(source) != count) goto bad_size;
        for (Py_ssize_t i = 0; i < count; i++) {
            target[i] = (int)PyLong_AsLong(PyList_GET_ITEM(source, i));
            if (PyErr_Occurred()) return -1;
        }
        return 0;
    }
    Py_buffer view;
    if (PyObject_GetBuffer(source, &view, PyBUF_CONTIG_RO) != 0) {
        PyErr_Clear();
        PyErr_Format(PyExc_TypeError, "%s must be a list or contiguous int buffer", name);
        return -1;
    }
    if (view.len != count * (Py_ssize_t)sizeof(*target)) {
        PyBuffer_Release(&view);
        goto bad_size;
    }
    memcpy(target, view.buf, (size_t)view.len);
    PyBuffer_Release(&view);
    return 0;

bad_size:
    PyErr_Format(PyExc_ValueError, "%s size must equal grid dimensions", name);
    return -1;
}

static int copy_u8_values(PyObject* source, uint8_t* target, Py_ssize_t count,
                          const char* name) {
    if (PyList_Check(source)) {
        if (PyList_Size(source) != count) goto bad_size;
        for (Py_ssize_t i = 0; i < count; i++) {
            target[i] = (uint8_t)PyLong_AsLong(PyList_GET_ITEM(source, i));
            if (PyErr_Occurred()) return -1;
        }
        return 0;
    }
    Py_buffer view;
    if (PyObject_GetBuffer(source, &view, PyBUF_CONTIG_RO) != 0) {
        PyErr_Clear();
        PyErr_Format(PyExc_TypeError, "%s must be a list or contiguous uint8 buffer", name);
        return -1;
    }
    if (view.len != count * (Py_ssize_t)sizeof(*target)) {
        PyBuffer_Release(&view);
        goto bad_size;
    }
    memcpy(target, view.buf, (size_t)view.len);
    PyBuffer_Release(&view);
    return 0;

bad_size:
    PyErr_Format(PyExc_ValueError, "%s size must equal grid dimensions", name);
    return -1;
}

/* After find_cells_grid_coverage has produced cell_ids + coverage + errors,
 * build coverage-driven analytical error segments for plotting.
 *
 * Cost is dominated by alea_get_slice_curves; gate this behind an opt-in
 * flag from the caller. Silently no-ops on curve-computation failure.
 */
static PyObject* compute_error_lines(alea_system_t* sys,
                                     const alea_slice_view_t* view,
                                     int nu, int nv,
                                     int universe_depth,
                                     const int* cell_ids,
                                     int* secondary_cell_ids,
                                     uint8_t* coverage,
                                     uint8_t* errors) {
    alea_slice_curves_t* curves = get_slice_curves_allow_threads(sys, view);
    if (!curves) return NULL;

    alea_slice_error_result_t* res;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    if (coverage) {
        alea_check_grid_overlaps_curves(sys, view, curves, nu, nv,
                                        universe_depth, cell_ids, errors);
        for (int i = 0; i < nu * nv; i++) {
            if (errors[i] == ALEA_GRID_OVERLAP) {
                coverage[i] = ALEA_COVERAGE_MULTI;
            } else if (errors[i] == ALEA_GRID_UNDEFINED || cell_ids[i] < 0) {
                coverage[i] = ALEA_COVERAGE_NONE;
            }
        }
    }

    /* These grid-coverage slice-error checkers are deprecated in libalea in
     * favour of alea_validate_geometry_slice. They remain behind the explicit
     * exact_coverage_grid analysis API while callers migrate to the structured
     * validator pipeline. Suppress the deprecation notice for this intentional
     * legacy use. */
#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#endif
    res = coverage
        ? alea_check_slice_errors_grid_ex(view, curves, cell_ids, coverage,
                                          errors, nu, nv)
        : alea_check_slice_errors_grid(sys, view, curves, cell_ids,
                                       errors, nu, nv);
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) {
        if (res) alea_slice_errors_free(res);
        alea_slice_curves_free(curves);
        return NULL;
    }

    PyObject* py_list = PyList_New(0);
    if (!py_list) {
        if (res) alea_slice_errors_free(res);
        alea_slice_curves_free(curves);
        return NULL;
    }

    if (res) {
        for (size_t i = 0; i < res->error_count; i++) {
            const alea_slice_error_t* e = &res->errors[i];
            alea_curve_t c;
            if (alea_slice_curves_get(curves, e->curve_index, &c) < 0) continue;

            PyObject* d = PyDict_New();
            const char* etype = (e->type == ALEA_SLICE_ERR_OVERLAP)
                                ? "overlap" : "gap";
            dict_set_new(d, "type", PyUnicode_FromString(etype));
            dict_set_new(d, "surface_id", PyLong_FromLong(e->surface_id));
            dict_set_new(d, "t_start", PyFloat_FromDouble(e->t_start));
            dict_set_new(d, "t_end", PyFloat_FromDouble(e->t_end));
            inline_curve_into_dict(d, &c);
            PyList_Append(py_list, d);
            Py_DECREF(d);
        }
        alea_slice_errors_free(res);
    }

    alea_slice_curves_free(curves);
    return py_list;
}

/* Full grid diagnostics already classify every pixel before analytical error
 * curves are requested.  Building the full slice-curve collection for a clean
 * viewport is pure overhead and is particularly expensive on large decks.
 * Keep the public result shape stable by publishing an empty Python list in
 * that case rather than omitting ``error_lines``. */
static int grid_has_errors(const uint8_t* errors, size_t count) {
    if (!errors) return 0;
    for (size_t index = 0; index < count; index++)
        if (errors[index] != 0) return 1;
    return 0;
}

static PyObject* PyAleaSystem_get_slice_curves_z(PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    double z, x_min, x_max, y_min, y_max;

    static char* kwlist[] = {"z", "x_min", "x_max", "y_min", "y_max", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "ddddd", kwlist,
            &z, &x_min, &x_max, &y_min, &y_max)) return NULL;

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, z, x_min, x_max, y_min, y_max);
    alea_slice_curves_t* curves =
        get_slice_curves_allow_threads(self->sys, &view);
    if (!curves) {
        if (!PyErr_Occurred())
            PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    size_t count = alea_slice_curves_count(curves);
    PyObject* list = PyList_New(count);
    if (!list) {
        alea_slice_curves_free(curves);
        return NULL;
    }

    for (size_t i = 0; i < count; i++) {
        alea_curve_t c;
        alea_slice_curves_get(curves, i, &c);

        PyObject* dict = PyDict_New();
        dict_set_new(dict, "type", PyUnicode_FromString(curve_type_to_string(c.type)));
        dict_set_new(dict, "surface_id", PyLong_FromLong(c.surface_id));

        switch (c.type) {
            case ALEA_CURVE_LINE:
            case ALEA_CURVE_LINE_SEGMENT:
                dict_set_new(dict, "point", Py_BuildValue("(dd)", c.data.line.point[0], c.data.line.point[1]));
                dict_set_new(dict, "direction", Py_BuildValue("(dd)", c.data.line.direction[0], c.data.line.direction[1]));
                if (c.type == ALEA_CURVE_LINE_SEGMENT) {
                    dict_set_new(dict, "t_min", PyFloat_FromDouble(c.t_min));
                    dict_set_new(dict, "t_max", PyFloat_FromDouble(c.t_max));
                }
                break;

            case ALEA_CURVE_CIRCLE:
            case ALEA_CURVE_ARC:
                dict_set_new(dict, "center", Py_BuildValue("(dd)", c.data.circle.center[0], c.data.circle.center[1]));
                dict_set_new(dict, "radius", PyFloat_FromDouble(c.data.circle.radius));
                if (c.type == ALEA_CURVE_ARC) {
                    dict_set_new(dict, "theta_start", PyFloat_FromDouble(c.t_min));
                    dict_set_new(dict, "theta_end", PyFloat_FromDouble(c.t_max));
                }
                break;

            case ALEA_CURVE_ELLIPSE:
            case ALEA_CURVE_ELLIPSE_ARC:
                dict_set_new(dict, "center", Py_BuildValue("(dd)", c.data.ellipse.center[0], c.data.ellipse.center[1]));
                dict_set_new(dict, "semi_a", PyFloat_FromDouble(c.data.ellipse.semi_a));
                dict_set_new(dict, "semi_b", PyFloat_FromDouble(c.data.ellipse.semi_b));
                dict_set_new(dict, "angle", PyFloat_FromDouble(c.data.ellipse.angle));
                if (c.type == ALEA_CURVE_ELLIPSE_ARC) {
                    dict_set_new(dict, "theta_start", PyFloat_FromDouble(c.t_min));
                    dict_set_new(dict, "theta_end", PyFloat_FromDouble(c.t_max));
                }
                break;

            case ALEA_CURVE_POLYGON: {
                PyObject* verts = PyList_New(c.data.polygon.count);
                for (int j = 0; j < c.data.polygon.count; j++) {
                    PyList_SET_ITEM(verts, j, Py_BuildValue("(dd)",
                        c.data.polygon.vertices[j][0], c.data.polygon.vertices[j][1]));
                }
                dict_set_new(dict, "vertices", verts);
                dict_set_new(dict, "closed", PyBool_FromLong(c.data.polygon.closed));
                break;
            }

            case ALEA_CURVE_PARALLEL_LINES:
                dict_set_new(dict, "point1", Py_BuildValue("(dd)", c.data.parallel_lines.point1[0], c.data.parallel_lines.point1[1]));
                dict_set_new(dict, "point2", Py_BuildValue("(dd)", c.data.parallel_lines.point2[0], c.data.parallel_lines.point2[1]));
                dict_set_new(dict, "direction", Py_BuildValue("(dd)", c.data.parallel_lines.direction[0], c.data.parallel_lines.direction[1]));
                break;

            default:
                break;
        }

        PyList_SET_ITEM(list, i, dict);
    }

    /* Build result with curves and bounds */
    double u_min, u_max, v_min, v_max;
    alea_slice_curves_bounds(curves, &u_min, &u_max, &v_min, &v_max);

    PyObject* result = PyDict_New();
    dict_set_new(result, "curves", list);
    dict_set_new(result, "u_min", PyFloat_FromDouble(u_min));
    dict_set_new(result, "u_max", PyFloat_FromDouble(u_max));
    dict_set_new(result, "v_min", PyFloat_FromDouble(v_min));
    dict_set_new(result, "v_max", PyFloat_FromDouble(v_max));
    dict_set_new(result, "x_min", PyFloat_FromDouble(x_min));
    dict_set_new(result, "x_max", PyFloat_FromDouble(x_max));
    dict_set_new(result, "y_min", PyFloat_FromDouble(y_min));
    dict_set_new(result, "y_max", PyFloat_FromDouble(y_max));
    alea_slice_curves_free(curves);
    return result;
}

static PyObject* PyAleaSystem_get_slice_curves_y(PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    double y, x_min, x_max, z_min, z_max;

    static char* kwlist[] = {"y", "x_min", "x_max", "z_min", "z_max", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "ddddd", kwlist,
            &y, &x_min, &x_max, &z_min, &z_max)) return NULL;

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 1, y, x_min, x_max, z_min, z_max);
    alea_slice_curves_t* curves =
        get_slice_curves_allow_threads(self->sys, &view);
    if (!curves) {
        if (!PyErr_Occurred())
            PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    /* Same processing as Z slice */
    size_t count = alea_slice_curves_count(curves);
    PyObject* list = PyList_New(count);
    if (!list) { alea_slice_curves_free(curves); return NULL; }

    for (size_t i = 0; i < count; i++) {
        alea_curve_t c;
        alea_slice_curves_get(curves, i, &c);

        PyObject* dict = PyDict_New();
        dict_set_new(dict, "type", PyUnicode_FromString(curve_type_to_string(c.type)));
        dict_set_new(dict, "surface_id", PyLong_FromLong(c.surface_id));

        switch (c.type) {
            case ALEA_CURVE_LINE:
            case ALEA_CURVE_LINE_SEGMENT:
                dict_set_new(dict, "point", Py_BuildValue("(dd)", c.data.line.point[0], c.data.line.point[1]));
                dict_set_new(dict, "direction", Py_BuildValue("(dd)", c.data.line.direction[0], c.data.line.direction[1]));
                break;
            case ALEA_CURVE_CIRCLE:
            case ALEA_CURVE_ARC:
                dict_set_new(dict, "center", Py_BuildValue("(dd)", c.data.circle.center[0], c.data.circle.center[1]));
                dict_set_new(dict, "radius", PyFloat_FromDouble(c.data.circle.radius));
                break;
            case ALEA_CURVE_ELLIPSE:
            case ALEA_CURVE_ELLIPSE_ARC:
                dict_set_new(dict, "center", Py_BuildValue("(dd)", c.data.ellipse.center[0], c.data.ellipse.center[1]));
                dict_set_new(dict, "semi_a", PyFloat_FromDouble(c.data.ellipse.semi_a));
                dict_set_new(dict, "semi_b", PyFloat_FromDouble(c.data.ellipse.semi_b));
                dict_set_new(dict, "angle", PyFloat_FromDouble(c.data.ellipse.angle));
                break;
            case ALEA_CURVE_POLYGON: {
                PyObject* verts = PyList_New(c.data.polygon.count);
                for (int j = 0; j < c.data.polygon.count; j++) {
                    PyList_SET_ITEM(verts, j, Py_BuildValue("(dd)",
                        c.data.polygon.vertices[j][0], c.data.polygon.vertices[j][1]));
                }
                dict_set_new(dict, "vertices", verts);
                dict_set_new(dict, "closed", PyBool_FromLong(c.data.polygon.closed));
                break;
            }
            case ALEA_CURVE_PARALLEL_LINES:
                dict_set_new(dict, "point1", Py_BuildValue("(dd)", c.data.parallel_lines.point1[0], c.data.parallel_lines.point1[1]));
                dict_set_new(dict, "point2", Py_BuildValue("(dd)", c.data.parallel_lines.point2[0], c.data.parallel_lines.point2[1]));
                dict_set_new(dict, "direction", Py_BuildValue("(dd)", c.data.parallel_lines.direction[0], c.data.parallel_lines.direction[1]));
                break;
            default:
                break;
        }
        PyList_SET_ITEM(list, i, dict);
    }

    double u_min, u_max, v_min, v_max;
    alea_slice_curves_bounds(curves, &u_min, &u_max, &v_min, &v_max);

    PyObject* result = PyDict_New();
    dict_set_new(result, "curves", list);
    dict_set_new(result, "u_min", PyFloat_FromDouble(u_min));
    dict_set_new(result, "u_max", PyFloat_FromDouble(u_max));
    dict_set_new(result, "v_min", PyFloat_FromDouble(v_min));
    dict_set_new(result, "v_max", PyFloat_FromDouble(v_max));
    /* Viewport bounds for XZ plane (Y slice) */
    dict_set_new(result, "x_min", PyFloat_FromDouble(x_min));
    dict_set_new(result, "x_max", PyFloat_FromDouble(x_max));
    dict_set_new(result, "z_min", PyFloat_FromDouble(z_min));
    dict_set_new(result, "z_max", PyFloat_FromDouble(z_max));
    alea_slice_curves_free(curves);
    return result;
}

static PyObject* PyAleaSystem_get_slice_curves_x(PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    double x, y_min, y_max, z_min, z_max;

    static char* kwlist[] = {"x", "y_min", "y_max", "z_min", "z_max", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "ddddd", kwlist,
            &x, &y_min, &y_max, &z_min, &z_max)) return NULL;

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 0, x, y_min, y_max, z_min, z_max);
    alea_slice_curves_t* curves =
        get_slice_curves_allow_threads(self->sys, &view);
    if (!curves) {
        if (!PyErr_Occurred())
            PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    size_t count = alea_slice_curves_count(curves);
    PyObject* list = PyList_New(count);
    if (!list) { alea_slice_curves_free(curves); return NULL; }

    for (size_t i = 0; i < count; i++) {
        alea_curve_t c;
        alea_slice_curves_get(curves, i, &c);

        PyObject* dict = PyDict_New();
        dict_set_new(dict, "type", PyUnicode_FromString(curve_type_to_string(c.type)));
        dict_set_new(dict, "surface_id", PyLong_FromLong(c.surface_id));

        switch (c.type) {
            case ALEA_CURVE_LINE:
            case ALEA_CURVE_LINE_SEGMENT:
                dict_set_new(dict, "point", Py_BuildValue("(dd)", c.data.line.point[0], c.data.line.point[1]));
                dict_set_new(dict, "direction", Py_BuildValue("(dd)", c.data.line.direction[0], c.data.line.direction[1]));
                break;
            case ALEA_CURVE_CIRCLE:
            case ALEA_CURVE_ARC:
                dict_set_new(dict, "center", Py_BuildValue("(dd)", c.data.circle.center[0], c.data.circle.center[1]));
                dict_set_new(dict, "radius", PyFloat_FromDouble(c.data.circle.radius));
                break;
            case ALEA_CURVE_ELLIPSE:
            case ALEA_CURVE_ELLIPSE_ARC:
                dict_set_new(dict, "center", Py_BuildValue("(dd)", c.data.ellipse.center[0], c.data.ellipse.center[1]));
                dict_set_new(dict, "semi_a", PyFloat_FromDouble(c.data.ellipse.semi_a));
                dict_set_new(dict, "semi_b", PyFloat_FromDouble(c.data.ellipse.semi_b));
                dict_set_new(dict, "angle", PyFloat_FromDouble(c.data.ellipse.angle));
                break;
            case ALEA_CURVE_POLYGON: {
                PyObject* verts = PyList_New(c.data.polygon.count);
                for (int j = 0; j < c.data.polygon.count; j++) {
                    PyList_SET_ITEM(verts, j, Py_BuildValue("(dd)",
                        c.data.polygon.vertices[j][0], c.data.polygon.vertices[j][1]));
                }
                dict_set_new(dict, "vertices", verts);
                dict_set_new(dict, "closed", PyBool_FromLong(c.data.polygon.closed));
                break;
            }
            case ALEA_CURVE_PARALLEL_LINES:
                dict_set_new(dict, "point1", Py_BuildValue("(dd)", c.data.parallel_lines.point1[0], c.data.parallel_lines.point1[1]));
                dict_set_new(dict, "point2", Py_BuildValue("(dd)", c.data.parallel_lines.point2[0], c.data.parallel_lines.point2[1]));
                dict_set_new(dict, "direction", Py_BuildValue("(dd)", c.data.parallel_lines.direction[0], c.data.parallel_lines.direction[1]));
                break;
            default:
                break;
        }
        PyList_SET_ITEM(list, i, dict);
    }

    double u_min, u_max, v_min, v_max;
    alea_slice_curves_bounds(curves, &u_min, &u_max, &v_min, &v_max);

    PyObject* result = PyDict_New();
    dict_set_new(result, "curves", list);
    dict_set_new(result, "u_min", PyFloat_FromDouble(u_min));
    dict_set_new(result, "u_max", PyFloat_FromDouble(u_max));
    dict_set_new(result, "v_min", PyFloat_FromDouble(v_min));
    dict_set_new(result, "v_max", PyFloat_FromDouble(v_max));
    /* Viewport bounds for YZ plane (X slice) */
    dict_set_new(result, "y_min", PyFloat_FromDouble(y_min));
    dict_set_new(result, "y_max", PyFloat_FromDouble(y_max));
    dict_set_new(result, "z_min", PyFloat_FromDouble(z_min));
    dict_set_new(result, "z_max", PyFloat_FromDouble(z_max));
    alea_slice_curves_free(curves);
    return result;
}

/* Arbitrary plane slice curves */

static PyObject* PyAleaSystem_get_slice_curves(PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    PyObject* origin_obj;
    PyObject* normal_obj;
    PyObject* up_obj;
    double u_min, u_max, v_min, v_max;

    static char* kwlist[] = {"origin", "normal", "up", "u_min", "u_max", "v_min", "v_max", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOdddd", kwlist,
            &origin_obj, &normal_obj, &up_obj, &u_min, &u_max, &v_min, &v_max)) return NULL;

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    /* Parse tuples */
    double ox, oy, oz, nx, ny, nz, ux, uy, uz;
    if (!PyArg_ParseTuple(origin_obj, "ddd", &ox, &oy, &oz)) {
        PyErr_SetString(PyExc_TypeError, "origin must be (x, y, z) tuple");
        return NULL;
    }
    if (!PyArg_ParseTuple(normal_obj, "ddd", &nx, &ny, &nz)) {
        PyErr_SetString(PyExc_TypeError, "normal must be (nx, ny, nz) tuple");
        return NULL;
    }
    if (!PyArg_ParseTuple(up_obj, "ddd", &ux, &uy, &uz)) {
        PyErr_SetString(PyExc_TypeError, "up must be (ux, uy, uz) tuple");
        return NULL;
    }

    /* Initialize view */
    alea_slice_view_t view;
    alea_slice_view_init(&view, ox, oy, oz, nx, ny, nz, ux, uy, uz,
                              u_min, u_max, v_min, v_max);

    /* Get curves */
    alea_slice_curves_t* curves =
        get_slice_curves_allow_threads(self->sys, &view);
    if (!curves) {
        if (!PyErr_Occurred())
            PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    /* Build result (same as axis-aligned versions) */
    size_t count = alea_slice_curves_count(curves);
    PyObject* list = PyList_New(count);
    if (!list) { alea_slice_curves_free(curves); return NULL; }

    for (size_t i = 0; i < count; i++) {
        alea_curve_t c;
        alea_slice_curves_get(curves, i, &c);

        PyObject* dict = PyDict_New();
        dict_set_new(dict, "type", PyUnicode_FromString(curve_type_to_string(c.type)));
        dict_set_new(dict, "surface_id", PyLong_FromLong(c.surface_id));

        switch (c.type) {
            case ALEA_CURVE_LINE:
            case ALEA_CURVE_LINE_SEGMENT:
                dict_set_new(dict, "point", Py_BuildValue("(dd)", c.data.line.point[0], c.data.line.point[1]));
                dict_set_new(dict, "direction", Py_BuildValue("(dd)", c.data.line.direction[0], c.data.line.direction[1]));
                break;
            case ALEA_CURVE_CIRCLE:
            case ALEA_CURVE_ARC:
                dict_set_new(dict, "center", Py_BuildValue("(dd)", c.data.circle.center[0], c.data.circle.center[1]));
                dict_set_new(dict, "radius", PyFloat_FromDouble(c.data.circle.radius));
                if (c.type == ALEA_CURVE_ARC) {
                    dict_set_new(dict, "theta_start", PyFloat_FromDouble(c.t_min));
                    dict_set_new(dict, "theta_end", PyFloat_FromDouble(c.t_max));
                }
                break;
            case ALEA_CURVE_ELLIPSE:
            case ALEA_CURVE_ELLIPSE_ARC:
                dict_set_new(dict, "center", Py_BuildValue("(dd)", c.data.ellipse.center[0], c.data.ellipse.center[1]));
                dict_set_new(dict, "semi_a", PyFloat_FromDouble(c.data.ellipse.semi_a));
                dict_set_new(dict, "semi_b", PyFloat_FromDouble(c.data.ellipse.semi_b));
                dict_set_new(dict, "angle", PyFloat_FromDouble(c.data.ellipse.angle));
                if (c.type == ALEA_CURVE_ELLIPSE_ARC) {
                    dict_set_new(dict, "theta_start", PyFloat_FromDouble(c.t_min));
                    dict_set_new(dict, "theta_end", PyFloat_FromDouble(c.t_max));
                }
                break;
            case ALEA_CURVE_POLYGON: {
                PyObject* verts = PyList_New(c.data.polygon.count);
                for (int j = 0; j < c.data.polygon.count; j++) {
                    PyList_SET_ITEM(verts, j, Py_BuildValue("(dd)",
                        c.data.polygon.vertices[j][0], c.data.polygon.vertices[j][1]));
                }
                dict_set_new(dict, "vertices", verts);
                dict_set_new(dict, "closed", PyBool_FromLong(c.data.polygon.closed));
                break;
            }
            case ALEA_CURVE_PARALLEL_LINES:
                dict_set_new(dict, "point1", Py_BuildValue("(dd)", c.data.parallel_lines.point1[0], c.data.parallel_lines.point1[1]));
                dict_set_new(dict, "point2", Py_BuildValue("(dd)", c.data.parallel_lines.point2[0], c.data.parallel_lines.point2[1]));
                dict_set_new(dict, "direction", Py_BuildValue("(dd)", c.data.parallel_lines.direction[0], c.data.parallel_lines.direction[1]));
                break;
            default:
                break;
        }
        PyList_SET_ITEM(list, i, dict);
    }

    double cu_min, cu_max, cv_min, cv_max;
    alea_slice_curves_bounds(curves, &cu_min, &cu_max, &cv_min, &cv_max);

    PyObject* result = PyDict_New();
    dict_set_new(result, "curves", list);
    dict_set_new(result, "u_min", PyFloat_FromDouble(cu_min));
    dict_set_new(result, "u_max", PyFloat_FromDouble(cu_max));
    dict_set_new(result, "v_min", PyFloat_FromDouble(cv_min));
    dict_set_new(result, "v_max", PyFloat_FromDouble(cv_max));
    /* Viewport bounds */
    dict_set_new(result, "view_u_min", PyFloat_FromDouble(u_min));
    dict_set_new(result, "view_u_max", PyFloat_FromDouble(u_max));
    dict_set_new(result, "view_v_min", PyFloat_FromDouble(v_min));
    dict_set_new(result, "view_v_max", PyFloat_FromDouble(v_max));

    alea_slice_curves_free(curves);
    return result;
}

/* Arbitrary plane grid query */

static PyObject* PyAleaSystem_find_cells_grid(PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    PyObject* origin_obj;
    PyObject* normal_obj;
    PyObject* up_obj;
    double u_min, u_max, v_min, v_max;
    int nu, nv;
    int universe_depth = -1;  /* Default: innermost cell */
    const char* error_mode = "none";
    int as_buffers = 0;
    int detect_errors = 0;
    int full_errors = 0;

    static char* kwlist[] = {"origin", "normal", "up", "u_min", "u_max", "v_min", "v_max", "nu", "nv",
                             "universe_depth", "error_mode", "_as_buffers", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOddddii|izp", kwlist,
            &origin_obj, &normal_obj, &up_obj,
            &u_min, &u_max, &v_min, &v_max, &nu, &nv,
            &universe_depth, &error_mode, &as_buffers)) return NULL;
    if (parse_grid_error_mode(error_mode, &detect_errors, &full_errors) < 0)
        return NULL;

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    /* Parse tuples */
    double ox, oy, oz, nx, ny, nz, ux, uy, uz;
    if (!PyArg_ParseTuple(origin_obj, "ddd", &ox, &oy, &oz)) {
        PyErr_SetString(PyExc_TypeError, "origin must be (x, y, z) tuple");
        return NULL;
    }
    if (!PyArg_ParseTuple(normal_obj, "ddd", &nx, &ny, &nz)) {
        PyErr_SetString(PyExc_TypeError, "normal must be (nx, ny, nz) tuple");
        return NULL;
    }
    if (!PyArg_ParseTuple(up_obj, "ddd", &ux, &uy, &uz)) {
        PyErr_SetString(PyExc_TypeError, "up must be (ux, uy, uz) tuple");
        return NULL;
    }

    if (nu <= 0 || nv <= 0) {
        PyErr_SetString(PyExc_ValueError, "Grid dimensions must be positive");
        return NULL;
    }

    /* Initialize view */
    alea_slice_view_t view;
    alea_slice_view_init(&view, ox, oy, oz, nx, ny, nz, ux, uy, uz,
                              u_min, u_max, v_min, v_max);

    /* Allocate arrays */
    int* cell_ids = malloc(nu * nv * sizeof(int));
    int* mat_ids = malloc(nu * nv * sizeof(int));
    uint8_t* errors = detect_errors ? malloc(nu * nv * sizeof(uint8_t)) : NULL;
    uint8_t* coverage = detect_errors ? malloc(nu * nv * sizeof(uint8_t)) : NULL;
    int* secondary_ids = detect_errors ? malloc(nu * nv * sizeof(int)) : NULL;
    if (!cell_ids || !mat_ids || (detect_errors && (!errors || !coverage || !secondary_ids))) {
        free(cell_ids); free(mat_ids); free(errors); free(coverage); free(secondary_ids);
        return PyErr_NoMemory();
    }

    int grid_rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    grid_rc = detect_errors
        ? alea_find_cells_grid_coverage(self->sys, &view, nu, nv,
                                        universe_depth,
                                        full_errors
                                            ? (ALEA_GRID_COVERAGE_EXACT | ALEA_GRID_SECONDARY_CELL_IDS)
                                            : ALEA_GRID_COVERAGE_FAST,
                                        cell_ids, mat_ids, secondary_ids,
                                        coverage, errors)
        : alea_find_cells_grid(self->sys, &view, nu, nv,
                               universe_depth, cell_ids, mat_ids, errors);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) {
        free(cell_ids); free(mat_ids); free(errors); free(coverage); free(secondary_ids);
        return NULL;
    }
    if (grid_rc < 0) {
        free(cell_ids); free(mat_ids); free(errors); free(coverage); free(secondary_ids);
        PyErr_SetString(PyExc_RuntimeError, "Grid query failed");
        return NULL;
    }

    PyObject* error_lines_list = NULL;
    if (detect_errors && full_errors && coverage) {
        old_sigint = install_sigint();
        Py_BEGIN_ALLOW_THREADS
        alea_filter_grid_boundary_ambiguities(
            self->sys, &view, nu, nv, universe_depth,
            cell_ids, secondary_ids, coverage, errors, NULL);
        Py_END_ALLOW_THREADS
        if (restore_sigint(old_sigint)) {
            free(cell_ids); free(mat_ids); free(errors); free(coverage); free(secondary_ids);
            return NULL;
        }
    }
    if (detect_errors && full_errors && errors) {
        error_lines_list = grid_has_errors(errors, (size_t)nu * (size_t)nv)
            ? compute_error_lines(self->sys, &view, nu, nv, universe_depth,
                                  cell_ids, secondary_ids, coverage, errors)
            : PyList_New(0);
    }
    if (!error_lines_list && PyErr_Occurred()) {
        free(cell_ids); free(mat_ids); free(errors); free(coverage);
        free(secondary_ids);
        return NULL;
    }

    /* Build result lists */
    PyObject* cell_list = as_buffers
        ? grid_array_from_owned_data(cell_ids, nu * nv, NPY_INT)
        : PyList_New(nu * nv);
    PyObject* mat_list = as_buffers
        ? grid_array_from_owned_data(mat_ids, nu * nv, NPY_INT)
        : PyList_New(nu * nv);
    if (!as_buffers) {
        for (int i = 0; i < nu * nv; i++) {
            PyList_SET_ITEM(cell_list, i, PyLong_FromLong(cell_ids[i]));
            PyList_SET_ITEM(mat_list, i, PyLong_FromLong(mat_ids[i]));
        }
    }

    PyObject* error_list = NULL;
    if (detect_errors && errors) {
        error_list = as_buffers
            ? grid_array_from_owned_data(errors, nu * nv, NPY_UBYTE)
            : PyList_New(nu * nv);
        if (!as_buffers) {
            for (int i = 0; i < nu * nv; i++) {
                PyList_SET_ITEM(error_list, i, PyLong_FromLong(errors[i]));
            }
        }
    }
    PyObject* coverage_list = NULL;
    PyObject* secondary_list = NULL;
    PyObject* component_list = NULL;
    if (detect_errors && coverage) {
        coverage_list = as_buffers
            ? grid_array_from_owned_data(coverage, nu * nv, NPY_UBYTE)
            : PyList_New(nu * nv);
        secondary_list = as_buffers
            ? grid_array_from_owned_data(secondary_ids, nu * nv, NPY_INT)
            : PyList_New(nu * nv);
        if (!as_buffers) {
            for (int i = 0; i < nu * nv; i++) {
                PyList_SET_ITEM(coverage_list, i, PyLong_FromLong(coverage[i]));
                PyList_SET_ITEM(secondary_list, i, PyLong_FromLong(secondary_ids[i]));
            }
        }
        component_list = build_error_components(cell_ids, secondary_ids,
                                                coverage, nu, nv);
    }

    if (!as_buffers) {
        free(cell_ids); free(mat_ids); free(errors); free(coverage); free(secondary_ids);
    }

    PyObject* result = PyDict_New();
    if (!result) {
        Py_XDECREF(cell_list); Py_XDECREF(mat_list); Py_XDECREF(error_list);
        Py_XDECREF(coverage_list); Py_XDECREF(secondary_list);
        Py_XDECREF(component_list); Py_XDECREF(error_lines_list);
        return NULL;
    }
    PyDict_SetItemString(result, "cell_ids", cell_list);
    PyDict_SetItemString(result, "material_ids", mat_list);
    if (error_list) {
        PyDict_SetItemString(result, "errors", error_list);
    }
    if (coverage_list) {
        PyDict_SetItemString(result, "coverage", coverage_list);
    }
    if (secondary_list) {
        PyDict_SetItemString(result, "secondary_cell_ids", secondary_list);
    }
    if (component_list) {
        PyDict_SetItemString(result, "error_components", component_list);
    }
    if (error_lines_list) {
        PyDict_SetItemString(result, "error_lines", error_lines_list);
        Py_DECREF(error_lines_list);
    }
    dict_set_new(result, "nu", PyLong_FromLong(nu));
    dict_set_new(result, "nv", PyLong_FromLong(nv));
    dict_set_new(result, "u_min", PyFloat_FromDouble(u_min));
    dict_set_new(result, "u_max", PyFloat_FromDouble(u_max));
    dict_set_new(result, "v_min", PyFloat_FromDouble(v_min));
    dict_set_new(result, "v_max", PyFloat_FromDouble(v_max));
    dict_set_new(result, "universe_depth", PyLong_FromLong(universe_depth));

    Py_DECREF(cell_list); Py_DECREF(mat_list); Py_XDECREF(error_list);
    Py_XDECREF(coverage_list); Py_XDECREF(secondary_list);
    Py_XDECREF(component_list);

    return result;
}

/* Grid cell queries for fills */

static PyObject* PyAleaSystem_find_cells_grid_z(PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    double z, x_min, x_max, y_min, y_max;
    int nx, ny;
    int universe_depth = -1;  /* Default: innermost cell */
    const char* error_mode = "none";
    int as_buffers = 0;
    int detect_errors = 0;
    int full_errors = 0;

    static char* kwlist[] = {"z", "x_min", "x_max", "y_min", "y_max", "nx", "ny",
                             "universe_depth", "error_mode", "_as_buffers", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "dddddii|izp", kwlist,
            &z, &x_min, &x_max, &y_min, &y_max, &nx, &ny,
            &universe_depth, &error_mode, &as_buffers)) return NULL;
    if (parse_grid_error_mode(error_mode, &detect_errors, &full_errors) < 0)
        return NULL;

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    if (nx <= 0 || ny <= 0) {
        PyErr_SetString(PyExc_ValueError, "Grid dimensions must be positive");
        return NULL;
    }

    int* cell_ids = malloc(nx * ny * sizeof(int));
    int* mat_ids = malloc(nx * ny * sizeof(int));
    uint8_t* errors = detect_errors ? malloc(nx * ny * sizeof(uint8_t)) : NULL;
    uint8_t* coverage = detect_errors ? malloc(nx * ny * sizeof(uint8_t)) : NULL;
    int* secondary_ids = detect_errors ? malloc(nx * ny * sizeof(int)) : NULL;
    if (!cell_ids || !mat_ids || (detect_errors && (!errors || !coverage || !secondary_ids))) {
        free(cell_ids); free(mat_ids); free(errors); free(coverage); free(secondary_ids);
        return PyErr_NoMemory();
    }

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 2, z, x_min, x_max, y_min, y_max);
    int grid_rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    grid_rc = detect_errors
        ? alea_find_cells_grid_coverage(self->sys, &view, nx, ny,
                                        universe_depth,
                                        full_errors
                                            ? (ALEA_GRID_COVERAGE_EXACT | ALEA_GRID_SECONDARY_CELL_IDS)
                                            : ALEA_GRID_COVERAGE_FAST,
                                        cell_ids, mat_ids, secondary_ids,
                                        coverage, errors)
        : alea_find_cells_grid(self->sys, &view, nx, ny,
                               universe_depth, cell_ids, mat_ids, errors);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) {
        free(cell_ids); free(mat_ids); free(errors); free(coverage); free(secondary_ids);
        return NULL;
    }
    if (grid_rc < 0) {
        free(cell_ids); free(mat_ids); free(errors); free(coverage); free(secondary_ids);
        PyErr_SetString(PyExc_RuntimeError, "Grid query failed");
        return NULL;
    }

    PyObject* error_lines_list = NULL;
    if (detect_errors && full_errors && coverage) {
        old_sigint = install_sigint();
        Py_BEGIN_ALLOW_THREADS
        alea_filter_grid_boundary_ambiguities(
            self->sys, &view, nx, ny, universe_depth,
            cell_ids, secondary_ids, coverage, errors, NULL);
        Py_END_ALLOW_THREADS
        if (restore_sigint(old_sigint)) {
            free(cell_ids); free(mat_ids); free(errors); free(coverage); free(secondary_ids);
            return NULL;
        }
    }
    if (detect_errors && full_errors && errors) {
        error_lines_list = grid_has_errors(errors, (size_t)nx * (size_t)ny)
            ? compute_error_lines(self->sys, &view, nx, ny, universe_depth,
                                  cell_ids, secondary_ids, coverage, errors)
            : PyList_New(0);
    }
    if (!error_lines_list && PyErr_Occurred()) {
        free(cell_ids); free(mat_ids); free(errors); free(coverage);
        free(secondary_ids);
        return NULL;
    }

    /* Create lists for results */
    PyObject* cell_list = as_buffers
        ? grid_array_from_owned_data(cell_ids, nx * ny, NPY_INT)
        : PyList_New(nx * ny);
    PyObject* mat_list = as_buffers
        ? grid_array_from_owned_data(mat_ids, nx * ny, NPY_INT)
        : PyList_New(nx * ny);
    if (!as_buffers) {
        for (int i = 0; i < nx * ny; i++) {
            PyList_SET_ITEM(cell_list, i, PyLong_FromLong(cell_ids[i]));
            PyList_SET_ITEM(mat_list, i, PyLong_FromLong(mat_ids[i]));
        }
    }

    PyObject* error_list = NULL;
    if (detect_errors && errors) {
        error_list = as_buffers
            ? grid_array_from_owned_data(errors, nx * ny, NPY_UBYTE)
            : PyList_New(nx * ny);
        if (!as_buffers) {
            for (int i = 0; i < nx * ny; i++) {
                PyList_SET_ITEM(error_list, i, PyLong_FromLong(errors[i]));
            }
        }
    }
    PyObject* coverage_list = NULL;
    PyObject* secondary_list = NULL;
    PyObject* component_list = NULL;
    if (detect_errors && coverage) {
        coverage_list = as_buffers
            ? grid_array_from_owned_data(coverage, nx * ny, NPY_UBYTE)
            : PyList_New(nx * ny);
        secondary_list = as_buffers
            ? grid_array_from_owned_data(secondary_ids, nx * ny, NPY_INT)
            : PyList_New(nx * ny);
        if (!as_buffers) {
            for (int i = 0; i < nx * ny; i++) {
                PyList_SET_ITEM(coverage_list, i, PyLong_FromLong(coverage[i]));
                PyList_SET_ITEM(secondary_list, i, PyLong_FromLong(secondary_ids[i]));
            }
        }
        component_list = build_error_components(cell_ids, secondary_ids,
                                                coverage, nx, ny);
    }

    if (!as_buffers) {
        free(cell_ids); free(mat_ids); free(errors); free(coverage); free(secondary_ids);
    }

    PyObject* result = PyDict_New();
    if (!result) {
        Py_XDECREF(cell_list); Py_XDECREF(mat_list); Py_XDECREF(error_list);
        Py_XDECREF(coverage_list); Py_XDECREF(secondary_list);
        Py_XDECREF(component_list); Py_XDECREF(error_lines_list);
        return NULL;
    }
    PyDict_SetItemString(result, "cell_ids", cell_list);
    PyDict_SetItemString(result, "material_ids", mat_list);
    if (error_list) {
        PyDict_SetItemString(result, "errors", error_list);
    }
    if (coverage_list) {
        PyDict_SetItemString(result, "coverage", coverage_list);
    }
    if (secondary_list) {
        PyDict_SetItemString(result, "secondary_cell_ids", secondary_list);
    }
    if (component_list) {
        PyDict_SetItemString(result, "error_components", component_list);
    }
    if (error_lines_list) {
        PyDict_SetItemString(result, "error_lines", error_lines_list);
        Py_DECREF(error_lines_list);
    }
    dict_set_new(result, "nx", PyLong_FromLong(nx));
    dict_set_new(result, "ny", PyLong_FromLong(ny));
    dict_set_new(result, "x_min", PyFloat_FromDouble(x_min));
    dict_set_new(result, "x_max", PyFloat_FromDouble(x_max));
    dict_set_new(result, "y_min", PyFloat_FromDouble(y_min));
    dict_set_new(result, "y_max", PyFloat_FromDouble(y_max));
    dict_set_new(result, "universe_depth", PyLong_FromLong(universe_depth));

    Py_DECREF(cell_list); Py_DECREF(mat_list); Py_XDECREF(error_list);
    Py_XDECREF(coverage_list); Py_XDECREF(secondary_list);
    Py_XDECREF(component_list);

    return result;
}

static PyObject* PyAleaSystem_find_cells_grid_y(PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    double y, x_min, x_max, z_min, z_max;
    int nx, nz;
    int universe_depth = -1;  /* Default: innermost cell */
    const char* error_mode = "none";
    int as_buffers = 0;
    int detect_errors = 0;
    int full_errors = 0;

    static char* kwlist[] = {"y", "x_min", "x_max", "z_min", "z_max", "nx", "nz",
                             "universe_depth", "error_mode", "_as_buffers", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "dddddii|izp", kwlist,
            &y, &x_min, &x_max, &z_min, &z_max, &nx, &nz,
            &universe_depth, &error_mode, &as_buffers)) return NULL;
    if (parse_grid_error_mode(error_mode, &detect_errors, &full_errors) < 0)
        return NULL;

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    if (nx <= 0 || nz <= 0) {
        PyErr_SetString(PyExc_ValueError, "Grid dimensions must be positive");
        return NULL;
    }

    int* cell_ids = malloc(nx * nz * sizeof(int));
    int* mat_ids = malloc(nx * nz * sizeof(int));
    uint8_t* errors = detect_errors ? malloc(nx * nz * sizeof(uint8_t)) : NULL;
    uint8_t* coverage = detect_errors ? malloc(nx * nz * sizeof(uint8_t)) : NULL;
    int* secondary_ids = detect_errors ? malloc(nx * nz * sizeof(int)) : NULL;
    if (!cell_ids || !mat_ids || (detect_errors && (!errors || !coverage || !secondary_ids))) {
        free(cell_ids); free(mat_ids); free(errors); free(coverage); free(secondary_ids);
        return PyErr_NoMemory();
    }

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 1, y, x_min, x_max, z_min, z_max);
    int grid_rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    grid_rc = detect_errors
        ? alea_find_cells_grid_coverage(self->sys, &view, nx, nz,
                                        universe_depth,
                                        full_errors
                                            ? (ALEA_GRID_COVERAGE_EXACT | ALEA_GRID_SECONDARY_CELL_IDS)
                                            : ALEA_GRID_COVERAGE_FAST,
                                        cell_ids, mat_ids, secondary_ids,
                                        coverage, errors)
        : alea_find_cells_grid(self->sys, &view, nx, nz,
                               universe_depth, cell_ids, mat_ids, errors);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) {
        free(cell_ids); free(mat_ids); free(errors); free(coverage); free(secondary_ids);
        return NULL;
    }
    if (grid_rc < 0) {
        free(cell_ids); free(mat_ids); free(errors); free(coverage); free(secondary_ids);
        PyErr_SetString(PyExc_RuntimeError, "Grid query failed");
        return NULL;
    }

    PyObject* error_lines_list = NULL;
    if (detect_errors && full_errors && coverage) {
        old_sigint = install_sigint();
        Py_BEGIN_ALLOW_THREADS
        alea_filter_grid_boundary_ambiguities(
            self->sys, &view, nx, nz, universe_depth,
            cell_ids, secondary_ids, coverage, errors, NULL);
        Py_END_ALLOW_THREADS
        if (restore_sigint(old_sigint)) {
            free(cell_ids); free(mat_ids); free(errors); free(coverage); free(secondary_ids);
            return NULL;
        }
    }
    if (detect_errors && full_errors && errors) {
        error_lines_list = grid_has_errors(errors, (size_t)nx * (size_t)nz)
            ? compute_error_lines(self->sys, &view, nx, nz, universe_depth,
                                  cell_ids, secondary_ids, coverage, errors)
            : PyList_New(0);
    }
    if (!error_lines_list && PyErr_Occurred()) {
        free(cell_ids); free(mat_ids); free(errors); free(coverage);
        free(secondary_ids);
        return NULL;
    }

    PyObject* cell_list = as_buffers ? grid_array_from_owned_data(cell_ids, nx * nz, NPY_INT) : PyList_New(nx * nz);
    PyObject* mat_list = as_buffers ? grid_array_from_owned_data(mat_ids, nx * nz, NPY_INT) : PyList_New(nx * nz);
    if (!as_buffers) for (int i = 0; i < nx * nz; i++) {
        PyList_SET_ITEM(cell_list, i, PyLong_FromLong(cell_ids[i]));
        PyList_SET_ITEM(mat_list, i, PyLong_FromLong(mat_ids[i]));
    }

    PyObject* error_list = NULL;
    if (detect_errors && errors) {
        error_list = as_buffers ? grid_array_from_owned_data(errors, nx * nz, NPY_UBYTE) : PyList_New(nx * nz);
        if (!as_buffers) for (int i = 0; i < nx * nz; i++) {
            PyList_SET_ITEM(error_list, i, PyLong_FromLong(errors[i]));
        }
    }
    PyObject* coverage_list = NULL;
    PyObject* secondary_list = NULL;
    PyObject* component_list = NULL;
    if (detect_errors && coverage) {
        coverage_list = as_buffers ? grid_array_from_owned_data(coverage, nx * nz, NPY_UBYTE) : PyList_New(nx * nz);
        secondary_list = as_buffers ? grid_array_from_owned_data(secondary_ids, nx * nz, NPY_INT) : PyList_New(nx * nz);
        if (!as_buffers) for (int i = 0; i < nx * nz; i++) {
            PyList_SET_ITEM(coverage_list, i, PyLong_FromLong(coverage[i]));
            PyList_SET_ITEM(secondary_list, i, PyLong_FromLong(secondary_ids[i]));
        }
        component_list = build_error_components(cell_ids, secondary_ids,
                                                coverage, nx, nz);
    }

    if (!as_buffers) {
        free(cell_ids); free(mat_ids); free(errors); free(coverage); free(secondary_ids);
    }

    PyObject* result = PyDict_New();
    if (!result) {
        Py_XDECREF(cell_list); Py_XDECREF(mat_list); Py_XDECREF(error_list);
        Py_XDECREF(coverage_list); Py_XDECREF(secondary_list);
        Py_XDECREF(component_list); Py_XDECREF(error_lines_list);
        return NULL;
    }
    PyDict_SetItemString(result, "cell_ids", cell_list);
    PyDict_SetItemString(result, "material_ids", mat_list);
    if (error_list) {
        PyDict_SetItemString(result, "errors", error_list);
    }
    if (coverage_list) {
        PyDict_SetItemString(result, "coverage", coverage_list);
    }
    if (secondary_list) {
        PyDict_SetItemString(result, "secondary_cell_ids", secondary_list);
    }
    if (component_list) {
        PyDict_SetItemString(result, "error_components", component_list);
    }
    if (error_lines_list) {
        PyDict_SetItemString(result, "error_lines", error_lines_list);
        Py_DECREF(error_lines_list);
    }
    dict_set_new(result, "nx", PyLong_FromLong(nx));
    dict_set_new(result, "nz", PyLong_FromLong(nz));
    dict_set_new(result, "x_min", PyFloat_FromDouble(x_min));
    dict_set_new(result, "x_max", PyFloat_FromDouble(x_max));
    dict_set_new(result, "z_min", PyFloat_FromDouble(z_min));
    dict_set_new(result, "z_max", PyFloat_FromDouble(z_max));
    dict_set_new(result, "universe_depth", PyLong_FromLong(universe_depth));

    Py_DECREF(cell_list); Py_DECREF(mat_list); Py_XDECREF(error_list);
    Py_XDECREF(coverage_list); Py_XDECREF(secondary_list);
    Py_XDECREF(component_list);

    return result;
}

static PyObject* PyAleaSystem_find_cells_grid_x(PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    double x, y_min, y_max, z_min, z_max;
    int ny, nz;
    int universe_depth = -1;  /* Default: innermost cell */
    const char* error_mode = "none";
    int as_buffers = 0;
    int detect_errors = 0;
    int full_errors = 0;

    static char* kwlist[] = {"x", "y_min", "y_max", "z_min", "z_max", "ny", "nz",
                             "universe_depth", "error_mode", "_as_buffers", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "dddddii|izp", kwlist,
            &x, &y_min, &y_max, &z_min, &z_max, &ny, &nz,
            &universe_depth, &error_mode, &as_buffers)) return NULL;
    if (parse_grid_error_mode(error_mode, &detect_errors, &full_errors) < 0)
        return NULL;

    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    if (ny <= 0 || nz <= 0) {
        PyErr_SetString(PyExc_ValueError, "Grid dimensions must be positive");
        return NULL;
    }

    int* cell_ids = malloc(ny * nz * sizeof(int));
    int* mat_ids = malloc(ny * nz * sizeof(int));
    uint8_t* errors = detect_errors ? malloc(ny * nz * sizeof(uint8_t)) : NULL;
    uint8_t* coverage = detect_errors ? malloc(ny * nz * sizeof(uint8_t)) : NULL;
    int* secondary_ids = detect_errors ? malloc(ny * nz * sizeof(int)) : NULL;
    if (!cell_ids || !mat_ids || (detect_errors && (!errors || !coverage || !secondary_ids))) {
        free(cell_ids); free(mat_ids); free(errors); free(coverage); free(secondary_ids);
        return PyErr_NoMemory();
    }

    alea_slice_view_t view;
    alea_slice_view_axis(&view, 0, x, y_min, y_max, z_min, z_max);
    int grid_rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    grid_rc = detect_errors
        ? alea_find_cells_grid_coverage(self->sys, &view, ny, nz,
                                        universe_depth,
                                        full_errors
                                            ? (ALEA_GRID_COVERAGE_EXACT | ALEA_GRID_SECONDARY_CELL_IDS)
                                            : ALEA_GRID_COVERAGE_FAST,
                                        cell_ids, mat_ids, secondary_ids,
                                        coverage, errors)
        : alea_find_cells_grid(self->sys, &view, ny, nz,
                               universe_depth, cell_ids, mat_ids, errors);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) {
        free(cell_ids); free(mat_ids); free(errors); free(coverage); free(secondary_ids);
        return NULL;
    }
    if (grid_rc < 0) {
        free(cell_ids); free(mat_ids); free(errors); free(coverage); free(secondary_ids);
        PyErr_SetString(PyExc_RuntimeError, "Grid query failed");
        return NULL;
    }

    PyObject* error_lines_list = NULL;
    if (detect_errors && full_errors && coverage) {
        old_sigint = install_sigint();
        Py_BEGIN_ALLOW_THREADS
        alea_filter_grid_boundary_ambiguities(
            self->sys, &view, ny, nz, universe_depth,
            cell_ids, secondary_ids, coverage, errors, NULL);
        Py_END_ALLOW_THREADS
        if (restore_sigint(old_sigint)) {
            free(cell_ids); free(mat_ids); free(errors); free(coverage); free(secondary_ids);
            return NULL;
        }
    }
    if (detect_errors && full_errors && errors) {
        error_lines_list = grid_has_errors(errors, (size_t)ny * (size_t)nz)
            ? compute_error_lines(self->sys, &view, ny, nz, universe_depth,
                                  cell_ids, secondary_ids, coverage, errors)
            : PyList_New(0);
    }
    if (!error_lines_list && PyErr_Occurred()) {
        free(cell_ids); free(mat_ids); free(errors); free(coverage);
        free(secondary_ids);
        return NULL;
    }

    PyObject* cell_list = as_buffers ? grid_array_from_owned_data(cell_ids, ny * nz, NPY_INT) : PyList_New(ny * nz);
    PyObject* mat_list = as_buffers ? grid_array_from_owned_data(mat_ids, ny * nz, NPY_INT) : PyList_New(ny * nz);
    if (!as_buffers) for (int i = 0; i < ny * nz; i++) {
        PyList_SET_ITEM(cell_list, i, PyLong_FromLong(cell_ids[i]));
        PyList_SET_ITEM(mat_list, i, PyLong_FromLong(mat_ids[i]));
    }

    PyObject* error_list = NULL;
    if (detect_errors && errors) {
        error_list = as_buffers ? grid_array_from_owned_data(errors, ny * nz, NPY_UBYTE) : PyList_New(ny * nz);
        if (!as_buffers) for (int i = 0; i < ny * nz; i++) {
            PyList_SET_ITEM(error_list, i, PyLong_FromLong(errors[i]));
        }
    }
    PyObject* coverage_list = NULL;
    PyObject* secondary_list = NULL;
    PyObject* component_list = NULL;
    if (detect_errors && coverage) {
        coverage_list = as_buffers ? grid_array_from_owned_data(coverage, ny * nz, NPY_UBYTE) : PyList_New(ny * nz);
        secondary_list = as_buffers ? grid_array_from_owned_data(secondary_ids, ny * nz, NPY_INT) : PyList_New(ny * nz);
        if (!as_buffers) for (int i = 0; i < ny * nz; i++) {
            PyList_SET_ITEM(coverage_list, i, PyLong_FromLong(coverage[i]));
            PyList_SET_ITEM(secondary_list, i, PyLong_FromLong(secondary_ids[i]));
        }
        component_list = build_error_components(cell_ids, secondary_ids,
                                                coverage, ny, nz);
    }

    if (!as_buffers) {
        free(cell_ids); free(mat_ids); free(errors); free(coverage); free(secondary_ids);
    }

    PyObject* result = PyDict_New();
    if (!result) {
        Py_XDECREF(cell_list); Py_XDECREF(mat_list); Py_XDECREF(error_list);
        Py_XDECREF(coverage_list); Py_XDECREF(secondary_list);
        Py_XDECREF(component_list); Py_XDECREF(error_lines_list);
        return NULL;
    }
    PyDict_SetItemString(result, "cell_ids", cell_list);
    PyDict_SetItemString(result, "material_ids", mat_list);
    if (error_list) {
        PyDict_SetItemString(result, "errors", error_list);
    }
    if (coverage_list) {
        PyDict_SetItemString(result, "coverage", coverage_list);
    }
    if (secondary_list) {
        PyDict_SetItemString(result, "secondary_cell_ids", secondary_list);
    }
    if (component_list) {
        PyDict_SetItemString(result, "error_components", component_list);
    }
    if (error_lines_list) {
        PyDict_SetItemString(result, "error_lines", error_lines_list);
        Py_DECREF(error_lines_list);
    }
    dict_set_new(result, "ny", PyLong_FromLong(ny));
    dict_set_new(result, "nz", PyLong_FromLong(nz));
    dict_set_new(result, "y_min", PyFloat_FromDouble(y_min));
    dict_set_new(result, "y_max", PyFloat_FromDouble(y_max));
    dict_set_new(result, "z_min", PyFloat_FromDouble(z_min));
    dict_set_new(result, "z_max", PyFloat_FromDouble(z_max));
    dict_set_new(result, "universe_depth", PyLong_FromLong(universe_depth));

    Py_DECREF(cell_list); Py_DECREF(mat_list); Py_XDECREF(error_list);
    Py_XDECREF(coverage_list); Py_XDECREF(secondary_list);
    Py_XDECREF(component_list);

    return result;
}

/* Label positioning functions */

static PyObject* PyAleaSystem_find_label_positions(PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    PyObject* ids_obj;
    int width, height;
    int min_pixels = 100;

    static char* kwlist[] = {"ids", "width", "height", "min_pixels", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "Oii|i", kwlist,
            &ids_obj, &width, &height, &min_pixels)) return NULL;

    if (width <= 0 || height <= 0) {
        PyErr_SetString(PyExc_ValueError, "Width and height must be positive");
        return NULL;
    }

    int* ids = malloc(width * height * sizeof(int));
    if (!ids) return PyErr_NoMemory();
    if (copy_int_values(ids_obj, ids, (Py_ssize_t)width * height, "ids") < 0) {
        free(ids);
        return NULL;
    }

    /* Call C function */
    alea_label_position_t* labels = NULL;
    int count = 0;

    int rc = alea_find_label_positions(ids, width, height, min_pixels, &labels, &count);
    free(ids);

    if (rc != 0) {
        PyErr_SetString(PyExc_RuntimeError, "Label position computation failed");
        return NULL;
    }

    /* Build result list */
    PyObject* result = PyList_New(count);
    for (int i = 0; i < count; i++) {
        PyObject* label_dict = PyDict_New();
        dict_set_new(label_dict, "id", PyLong_FromLong(labels[i].id));
        dict_set_new(label_dict, "px", PyLong_FromLong(labels[i].px));
        dict_set_new(label_dict, "py", PyLong_FromLong(labels[i].py));
        dict_set_new(label_dict, "pixel_count", PyLong_FromLong(labels[i].pixel_count));
        PyList_SET_ITEM(result, i, label_dict);
    }

    free(labels);
    return result;
}

/* ============================================================================
 * PyAleaSystem Methods - Grid Overlap Check
 * ============================================================================ */

static PyObject* PyAleaSystem_check_grid_overlaps(PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    PyObject* origin_obj, *normal_obj, *up_obj;
    double u_min, u_max, v_min, v_max;
    int nu, nv;
    int universe_depth = -1;
    PyObject* cell_ids_obj;
    PyObject* errors_obj;

    static char* kwlist[] = {"origin", "normal", "up", "u_min", "u_max", "v_min", "v_max",
                             "nu", "nv", "cell_ids", "errors", "universe_depth", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOddddiiOO|i", kwlist,
            &origin_obj, &normal_obj, &up_obj, &u_min, &u_max, &v_min, &v_max,
            &nu, &nv, &cell_ids_obj, &errors_obj, &universe_depth)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }

    double ox, oy, oz, nx, ny, nz, ux, uy, uz;
    if (!PyArg_ParseTuple(origin_obj, "ddd", &ox, &oy, &oz)) return NULL;
    if (!PyArg_ParseTuple(normal_obj, "ddd", &nx, &ny, &nz)) return NULL;
    if (!PyArg_ParseTuple(up_obj, "ddd", &ux, &uy, &uz)) return NULL;

    Py_ssize_t size = nu * nv;
    int* cell_ids = malloc(size * sizeof(int));
    uint8_t* errors = malloc(size * sizeof(uint8_t));
    if (!cell_ids || !errors) { free(cell_ids); free(errors); return PyErr_NoMemory(); }

    if (copy_int_values(cell_ids_obj, cell_ids, size, "cell_ids") < 0 ||
        copy_u8_values(errors_obj, errors, size, "errors") < 0) {
        free(cell_ids);
        free(errors);
        return NULL;
    }

    alea_slice_view_t view;
    alea_slice_view_init(&view, ox, oy, oz, nx, ny, nz, ux, uy, uz,
                              u_min, u_max, v_min, v_max);

    int rc = alea_check_grid_overlaps(self->sys, &view, nu, nv,
                                           universe_depth, cell_ids, errors);
    free(cell_ids);

    if (rc < 0) {
        free(errors);
        PyErr_SetString(PyExc_RuntimeError, "Overlap check failed");
        return NULL;
    }

    /* Return updated errors list */
    PyObject* result = PyList_New(size);
    for (Py_ssize_t i = 0; i < size; i++) {
        PyList_SET_ITEM(result, i, PyLong_FromLong(errors[i]));
    }
    free(errors);
    return result;
}

/* ==========================================================================
 * Bounded tile coverage refinement
 * ========================================================================== */

static int tile_stats_dict_set_size_t(PyObject* dict, const char* key,
                                      size_t value) {
    PyObject* item = PyLong_FromSize_t(value);
    if (!item) return -1;
    int rc = PyDict_SetItemString(dict, key, item);
    Py_DECREF(item);
    return rc;
}

static int tile_stats_dict_set_bool(PyObject* dict, const char* key,
                                    int value) {
    PyObject* item = PyBool_FromLong(value);
    if (!item) return -1;
    int rc = PyDict_SetItemString(dict, key, item);
    Py_DECREF(item);
    return rc;
}

static PyObject* build_tile_coverage_stats_dict(
    const alea_tile_coverage_stats_t* stats) {
    PyObject* result = PyDict_New();
    if (!result) return NULL;
    if (tile_stats_dict_set_size_t(result, "tiles", stats->tiles) < 0 ||
        tile_stats_dict_set_size_t(result, "fallback_tiles",
                                   stats->fallback_tiles) < 0 ||
        tile_stats_dict_set_size_t(result, "skipped_tiles",
                                   stats->skipped_tiles) < 0 ||
        tile_stats_dict_set_size_t(result, "exact_fallback_pixels",
                                   stats->exact_fallback_pixels) < 0 ||
        tile_stats_dict_set_size_t(result, "skipped_pixels",
                                   stats->skipped_pixels) < 0 ||
        tile_stats_dict_set_size_t(result, "candidate_total",
                                   stats->candidate_total) < 0 ||
        tile_stats_dict_set_size_t(result, "candidate_max",
                                   stats->candidate_max) < 0 ||
        tile_stats_dict_set_size_t(result, "dedup_candidate_total",
                                   stats->dedup_candidate_total) < 0 ||
        tile_stats_dict_set_size_t(result, "dedup_candidate_max",
                                   stats->dedup_candidate_max) < 0 ||
        tile_stats_dict_set_size_t(result, "candidate_pixel_tests",
                                   stats->candidate_pixel_tests) < 0 ||
        tile_stats_dict_set_size_t(result, "path_groups",
                                   stats->path_groups) < 0 ||
        tile_stats_dict_set_size_t(result, "path_group_pixels_max",
                                   stats->path_group_pixels_max) < 0 ||
        tile_stats_dict_set_size_t(result, "path_group_candidates_max",
                                   stats->path_group_candidates_max) < 0 ||
        tile_stats_dict_set_size_t(result, "bbox_pixel_tests",
                                   stats->bbox_pixel_tests) < 0 ||
        tile_stats_dict_set_size_t(result, "bbox_pixel_rejects",
                                   stats->bbox_pixel_rejects) < 0 ||
        tile_stats_dict_set_size_t(result, "refined_pixels",
                                   stats->refined_pixels) < 0 ||
        tile_stats_dict_set_bool(result, "incomplete",
                                 stats->incomplete) < 0) {
        Py_DECREF(result);
        return NULL;
    }
    return result;
}

static PyObject* PyAleaSystem_refine_grid_coverage_tiles(
        PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    PyObject *origin_obj, *normal_obj, *up_obj;
    PyObject *coverage_obj, *errors_obj, *secondary_obj, *tile_mask_obj = Py_None;
    double u_min, u_max, v_min, v_max;
    int nu, nv, universe_depth = -1, tile_w = 16, tile_h = 16;
    Py_ssize_t max_candidates = 4096, max_evaluated_candidates = 4096,
               max_exact_fallback_pixels = 0;
    int chain_candidates = 0;
    static char* kwlist[] = {
        "origin", "normal", "up", "u_min", "u_max", "v_min", "v_max",
        "nu", "nv", "coverage", "errors", "secondary_cell_ids",
        "universe_depth", "tile_w", "tile_h", "max_candidates",
        "max_evaluated_candidates",
        "max_exact_fallback_pixels", "chain_candidates", "tile_mask", NULL};
    if (!PyArg_ParseTupleAndKeywords(
            args, kwds, "OOOddddiiOOO|iiinnnpO", kwlist,
            &origin_obj, &normal_obj, &up_obj, &u_min, &u_max, &v_min, &v_max,
            &nu, &nv, &coverage_obj, &errors_obj, &secondary_obj,
            &universe_depth, &tile_w, &tile_h, &max_candidates,
            &max_evaluated_candidates,
            &max_exact_fallback_pixels, &chain_candidates, &tile_mask_obj))
        return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    if (nu <= 0 || nv <= 0 || tile_w <= 0 || tile_h <= 0 ||
        max_candidates <= 0 || max_evaluated_candidates < 0 ||
        max_exact_fallback_pixels < 0) {
        PyErr_SetString(PyExc_ValueError, "grid dimensions and limits must be positive");
        return NULL;
    }
    if ((size_t)nu > SIZE_MAX / (size_t)nv ||
        (size_t)nu * (size_t)nv > (size_t)PY_SSIZE_T_MAX) {
        PyErr_SetString(PyExc_OverflowError, "grid is too large");
        return NULL;
    }
    Py_ssize_t size = (Py_ssize_t)((size_t)nu * (size_t)nv);
    double ox, oy, oz, nx, ny, nz, ux, uy, uz;
    if (!PyArg_ParseTuple(origin_obj, "ddd", &ox, &oy, &oz) ||
        !PyArg_ParseTuple(normal_obj, "ddd", &nx, &ny, &nz) ||
        !PyArg_ParseTuple(up_obj, "ddd", &ux, &uy, &uz))
        return NULL;

    uint8_t* coverage = malloc((size_t)size * sizeof(*coverage));
    uint8_t* errors = malloc((size_t)size * sizeof(*errors));
    int* secondary = malloc((size_t)size * sizeof(*secondary));
    uint8_t* provisional = calloc((size_t)size, sizeof(*provisional));
    if (!coverage || !errors || !secondary || !provisional) {
        free(coverage); free(errors); free(secondary); free(provisional);
        return PyErr_NoMemory();
    }
    if (copy_u8_values(coverage_obj, coverage, size, "coverage") < 0 ||
        copy_u8_values(errors_obj, errors, size, "errors") < 0 ||
        copy_int_values(secondary_obj, secondary, size, "secondary_cell_ids") < 0) {
        free(coverage); free(errors); free(secondary); free(provisional);
        return NULL;
    }

    uint8_t* tile_mask = NULL;
    size_t tile_mask_count = 0;
    if (tile_mask_obj != Py_None) {
        size_t tile_cols = ((size_t)nu + (size_t)tile_w - 1) / (size_t)tile_w;
        size_t tile_rows = ((size_t)nv + (size_t)tile_h - 1) / (size_t)tile_h;
        if (tile_cols > SIZE_MAX / tile_rows ||
            tile_cols * tile_rows > (size_t)PY_SSIZE_T_MAX) {
            free(coverage); free(errors); free(secondary); free(provisional);
            PyErr_SetString(PyExc_OverflowError, "tile mask is too large");
            return NULL;
        }
        tile_mask_count = tile_cols * tile_rows;
        tile_mask = malloc(tile_mask_count * sizeof(*tile_mask));
        if (!tile_mask) {
            free(coverage); free(errors); free(secondary); free(provisional);
            return PyErr_NoMemory();
        }
        if (copy_u8_values(tile_mask_obj, tile_mask,
                           (Py_ssize_t)tile_mask_count, "tile_mask") < 0) {
            free(tile_mask); free(coverage); free(errors); free(secondary);
            free(provisional);
            return NULL;
        }
    }

    alea_slice_view_t view;
    alea_slice_view_init(&view, ox, oy, oz, nx, ny, nz, ux, uy, uz,
                         u_min, u_max, v_min, v_max);
    alea_tile_refinement_options_t options;
    alea_tile_refinement_options_init(&options);
    options.max_candidates = (size_t)max_candidates;
    options.max_evaluated_candidates = (size_t)max_evaluated_candidates;
    options.max_exact_fallback_pixels = (size_t)max_exact_fallback_pixels;
    options.use_hier_chain_candidates = chain_candidates != 0;
    options.use_chain_bitset = chain_candidates != 0;
    options.tile_mask = tile_mask;
    options.tile_mask_count = tile_mask_count;
    options.provisional_mask = provisional;
    options.provisional_mask_count = (size_t)size;
    int rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    rc = alea_refine_grid_coverage_tiles_exact_ex(
        self->sys, &view, nu, nv, universe_depth, tile_w, tile_h, &options,
        secondary, coverage, errors);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) {
        free(tile_mask); free(coverage); free(errors); free(secondary); free(provisional);
        return NULL;
    }
    if (rc < 0) {
        free(tile_mask); free(coverage); free(errors); free(secondary); free(provisional);
        PyErr_SetString(PyExc_RuntimeError, "Tile coverage refinement failed");
        return NULL;
    }
    free(tile_mask);

    alea_tile_coverage_stats_t native_stats = alea_tile_coverage_stats_get();
    PyObject* result = PyDict_New();
    PyObject* stats = build_tile_coverage_stats_dict(&native_stats);
    PyObject* coverage_array = grid_array_from_owned_data(coverage, size, NPY_UBYTE);
    PyObject* errors_array = grid_array_from_owned_data(errors, size, NPY_UBYTE);
    PyObject* secondary_array = grid_array_from_owned_data(secondary, size, NPY_INT);
    PyObject* provisional_array = grid_array_from_owned_data(provisional, size, NPY_UBYTE);
    if (!result || !stats || !coverage_array || !errors_array || !secondary_array ||
        !provisional_array) {
        Py_XDECREF(result); Py_XDECREF(stats); Py_XDECREF(coverage_array);
        Py_XDECREF(errors_array); Py_XDECREF(secondary_array); Py_XDECREF(provisional_array);
        return NULL;
    }
    PyDict_SetItemString(result, "coverage", coverage_array);
    PyDict_SetItemString(result, "errors", errors_array);
    PyDict_SetItemString(result, "secondary_cell_ids", secondary_array);
    PyDict_SetItemString(result, "provisional_mask", provisional_array);
    PyDict_SetItemString(result, "stats", stats);
    if (tile_stats_dict_set_bool(result, "incomplete",
                                 native_stats.incomplete) < 0) {
        Py_DECREF(result); Py_DECREF(stats); Py_DECREF(coverage_array);
        Py_DECREF(errors_array); Py_DECREF(secondary_array);
        Py_DECREF(provisional_array);
        return NULL;
    }
    Py_DECREF(stats); Py_DECREF(coverage_array); Py_DECREF(errors_array);
    Py_DECREF(secondary_array); Py_DECREF(provisional_array);
    return result;
}

static PyObject* PyAleaSystem_refine_grid_coverage_paths(
        PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    PyObject *origin_obj, *normal_obj, *up_obj;
    PyObject *primary_obj, *coverage_obj, *errors_obj, *secondary_obj;
    PyObject* tile_mask_obj = Py_None;
    double u_min, u_max, v_min, v_max;
    int nu, nv, universe_depth = -1, tile_w = 16, tile_h = 16;
    Py_ssize_t max_candidates = 4096, max_exact_fallback_pixels = 0;
    int use_path_2d_index = 0;
    static char* kwlist[] = {
        "origin", "normal", "up", "u_min", "u_max", "v_min", "v_max",
        "nu", "nv", "cell_ids", "coverage", "errors", "secondary_cell_ids",
        "universe_depth", "tile_w", "tile_h", "max_candidates",
        "max_exact_fallback_pixels", "use_path_2d_index", "tile_mask", NULL};
    if (!PyArg_ParseTupleAndKeywords(
            args, kwds, "OOOddddiiOOOO|iiinnpO", kwlist,
            &origin_obj, &normal_obj, &up_obj, &u_min, &u_max, &v_min, &v_max,
            &nu, &nv, &primary_obj, &coverage_obj, &errors_obj, &secondary_obj,
            &universe_depth, &tile_w, &tile_h, &max_candidates,
            &max_exact_fallback_pixels, &use_path_2d_index, &tile_mask_obj))
        return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    if (nu <= 0 || nv <= 0 || tile_w <= 0 || tile_h <= 0 ||
        max_candidates <= 0 || max_exact_fallback_pixels < 0) {
        PyErr_SetString(PyExc_ValueError,
                        "grid dimensions and limits must be positive");
        return NULL;
    }
    if ((size_t)nu > SIZE_MAX / (size_t)nv ||
        (size_t)nu * (size_t)nv > (size_t)PY_SSIZE_T_MAX) {
        PyErr_SetString(PyExc_OverflowError, "grid is too large");
        return NULL;
    }
    Py_ssize_t size = (Py_ssize_t)((size_t)nu * (size_t)nv);
    double ox, oy, oz, nx, ny, nz, ux, uy, uz;
    if (!PyArg_ParseTuple(origin_obj, "ddd", &ox, &oy, &oz) ||
        !PyArg_ParseTuple(normal_obj, "ddd", &nx, &ny, &nz) ||
        !PyArg_ParseTuple(up_obj, "ddd", &ux, &uy, &uz))
        return NULL;

    int* primary = malloc((size_t)size * sizeof(*primary));
    uint8_t* coverage = malloc((size_t)size * sizeof(*coverage));
    uint8_t* errors = malloc((size_t)size * sizeof(*errors));
    int* secondary = malloc((size_t)size * sizeof(*secondary));
    uint32_t* path_ids = malloc((size_t)size * sizeof(*path_ids));
    uint8_t* provisional = calloc((size_t)size, sizeof(*provisional));
    if (!primary || !coverage || !errors || !secondary || !path_ids ||
        !provisional) {
        free(primary); free(coverage); free(errors); free(secondary);
        free(path_ids); free(provisional);
        return PyErr_NoMemory();
    }
    if (copy_int_values(primary_obj, primary, size, "cell_ids") < 0 ||
        copy_u8_values(coverage_obj, coverage, size, "coverage") < 0 ||
        copy_u8_values(errors_obj, errors, size, "errors") < 0 ||
        copy_int_values(secondary_obj, secondary, size,
                        "secondary_cell_ids") < 0) {
        free(primary); free(coverage); free(errors); free(secondary);
        free(path_ids); free(provisional);
        return NULL;
    }

    uint8_t* tile_mask = NULL;
    size_t tile_mask_count = 0;
    if (tile_mask_obj != Py_None) {
        size_t tile_cols = ((size_t)nu + (size_t)tile_w - 1) / (size_t)tile_w;
        size_t tile_rows = ((size_t)nv + (size_t)tile_h - 1) / (size_t)tile_h;
        if (tile_cols > SIZE_MAX / tile_rows ||
            tile_cols * tile_rows > (size_t)PY_SSIZE_T_MAX) {
            free(primary); free(coverage); free(errors); free(secondary);
            free(path_ids); free(provisional);
            PyErr_SetString(PyExc_OverflowError, "tile mask is too large");
            return NULL;
        }
        tile_mask_count = tile_cols * tile_rows;
        tile_mask = malloc(tile_mask_count * sizeof(*tile_mask));
        if (!tile_mask) {
            free(primary); free(coverage); free(errors); free(secondary);
            free(path_ids); free(provisional);
            return PyErr_NoMemory();
        }
        if (copy_u8_values(tile_mask_obj, tile_mask,
                           (Py_ssize_t)tile_mask_count, "tile_mask") < 0) {
            free(tile_mask); free(primary); free(coverage); free(errors);
            free(secondary); free(path_ids); free(provisional);
            return NULL;
        }
    }

    alea_slice_view_t view;
    alea_slice_view_init(&view, ox, oy, oz, nx, ny, nz, ux, uy, uz,
                         u_min, u_max, v_min, v_max);
    alea_slice_path_table_t paths = {0};
    alea_tile_refinement_options_t options;
    alea_tile_refinement_options_init(&options);
    options.max_candidates = (size_t)max_candidates;
    options.max_exact_fallback_pixels = (size_t)max_exact_fallback_pixels;
    options.use_path_2d_index = use_path_2d_index != 0;
    options.path_2d_bucket_limit = (size_t)max_candidates;
    options.tile_mask = tile_mask;
    options.tile_mask_count = tile_mask_count;
    options.provisional_mask = provisional;
    options.provisional_mask_count = (size_t)size;
    int grid_rc, refine_rc = -1;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    grid_rc = alea_find_cells_grid_paths_selected(
        self->sys, &view, nu, nv, universe_depth, tile_w, tile_h,
        tile_mask, tile_mask_count, primary, path_ids, &paths);
    if (grid_rc == 0) {
        refine_rc = alea_refine_grid_coverage_paths_exact_ex(
            self->sys, &view, nu, nv, universe_depth, tile_w, tile_h,
            primary, path_ids, &paths, &options,
            secondary, coverage, errors);
    }
    Py_END_ALLOW_THREADS
    int interrupted = restore_sigint(old_sigint);
    size_t path_count = paths.count;
    alea_slice_path_table_free(&paths);
    free(tile_mask);
    free(path_ids);
    free(primary);
    if (interrupted) {
        free(coverage); free(errors); free(secondary); free(provisional);
        return NULL;
    }
    if (grid_rc < 0 || refine_rc < 0) {
        free(coverage); free(errors); free(secondary); free(provisional);
        PyErr_SetString(PyExc_RuntimeError,
                        grid_rc < 0 ? "Path grid query failed" :
                                      "Path coverage refinement failed");
        return NULL;
    }

    alea_tile_coverage_stats_t native_stats = alea_tile_coverage_stats_get();
    PyObject* result = PyDict_New();
    PyObject* stats = build_tile_coverage_stats_dict(&native_stats);
    PyObject* coverage_array =
        grid_array_from_owned_data(coverage, size, NPY_UBYTE);
    PyObject* errors_array =
        grid_array_from_owned_data(errors, size, NPY_UBYTE);
    PyObject* secondary_array =
        grid_array_from_owned_data(secondary, size, NPY_INT);
    PyObject* provisional_array =
        grid_array_from_owned_data(provisional, size, NPY_UBYTE);
    if (!result || !stats || !coverage_array || !errors_array ||
        !secondary_array || !provisional_array) {
        Py_XDECREF(result); Py_XDECREF(stats); Py_XDECREF(coverage_array);
        Py_XDECREF(errors_array); Py_XDECREF(secondary_array);
        Py_XDECREF(provisional_array);
        return NULL;
    }
    PyDict_SetItemString(result, "coverage", coverage_array);
    PyDict_SetItemString(result, "errors", errors_array);
    PyDict_SetItemString(result, "secondary_cell_ids", secondary_array);
    PyDict_SetItemString(result, "provisional_mask", provisional_array);
    PyDict_SetItemString(result, "stats", stats);
    if (tile_stats_dict_set_size_t(result, "path_count", path_count) < 0 ||
        tile_stats_dict_set_bool(result, "incomplete",
                                 native_stats.incomplete) < 0) {
        Py_DECREF(result); Py_DECREF(stats); Py_DECREF(coverage_array);
        Py_DECREF(errors_array); Py_DECREF(secondary_array);
        Py_DECREF(provisional_array);
        return NULL;
    }
    Py_DECREF(stats); Py_DECREF(coverage_array); Py_DECREF(errors_array);
    Py_DECREF(secondary_array); Py_DECREF(provisional_array);
    return result;
}

/* ============================================================================
 * PyAleaSystem Methods - Surface Label Positions
 * ============================================================================ */

static PyObject* PyAleaSystem_find_surface_label_positions(PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    PyObject* origin_obj, *normal_obj, *up_obj;
    PyObject* boundary_ids_obj = Py_None;
    const char* boundary_by = "cell";
    double u_min, u_max, v_min, v_max;
    int width, height, margin = 20;

    static char* kwlist[] = {"origin", "normal", "up", "u_min", "u_max", "v_min", "v_max",
                             "width", "height", "margin", "boundary_ids", "boundary_by", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOddddii|iOs", kwlist,
            &origin_obj, &normal_obj, &up_obj, &u_min, &u_max, &v_min, &v_max,
            &width, &height, &margin, &boundary_ids_obj, &boundary_by)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }
    if (width <= 0 || height <= 0) {
        PyErr_SetString(PyExc_ValueError, "width and height must be positive");
        return NULL;
    }

    /* boundary_by names which rendered grid the caller passed as
     * boundary_ids; the visibility test itself reads that grid directly. */
    if (strcmp(boundary_by, "material") != 0 && strcmp(boundary_by, "cell") != 0) {
        PyErr_SetString(PyExc_ValueError, "boundary_by must be 'cell' or 'material'");
        return NULL;
    }

    /* Curve extraction currently materializes one analytical intersection
     * candidate per system surface before sparse attribution can shortlist
     * it. That is inappropriate for whole-plant models: it can exhaust a
     * notebook kernel before any label is emitted. Keep this public path
     * recoverable until viewport-aware surface culling is available. */
    const size_t max_surface_label_curves = 50000;
    size_t surface_count = alea_surface_count(self->sys);
    if (surface_count > max_surface_label_curves) {
        PyErr_Format(PyExc_MemoryError,
            "surface label extraction refused for %zu surfaces (limit %zu); "
            "use a smaller model/viewport or boundary-map labels",
            surface_count, max_surface_label_curves);
        return NULL;
    }

    double ox, oy, oz, nx, ny, nz, ux, uy, uz;
    if (!PyArg_ParseTuple(origin_obj, "ddd", &ox, &oy, &oz)) return NULL;
    if (!PyArg_ParseTuple(normal_obj, "ddd", &nx, &ny, &nz)) return NULL;
    if (!PyArg_ParseTuple(up_obj, "ddd", &ux, &uy, &uz)) return NULL;

    alea_slice_view_t view;
    alea_slice_view_init(&view, ox, oy, oz, nx, ny, nz, ux, uy, uz,
                              u_min, u_max, v_min, v_max);

    int* boundary_ids = NULL;
    if (boundary_ids_obj != Py_None) {
        Py_ssize_t size = (Py_ssize_t)width * height;
        boundary_ids = malloc((size_t)size * sizeof(int));
        if (!boundary_ids) {
            return PyErr_NoMemory();
        }
        if (copy_int_values(boundary_ids_obj, boundary_ids, size,
                            "boundary_ids") < 0) {
            free(boundary_ids);
            return NULL;
        }
    }

    alea_slice_curves_t* curves = NULL;
    alea_label_position_t* labels = NULL;
    int count = 0;
    int rc = -1;
    /* An explicitly unfiltered analytical request has no rendered boundary
     * grid to verify against.  Do not send it through provenance verification:
     * that path correctly rejects every candidate without an edge map. */
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    curves = alea_get_slice_curves(self->sys, &view);
    if (curves) {
        rc = boundary_ids
            ? alea_find_surface_label_positions_with_provenance(
                  self->sys, &view, curves, boundary_ids,
                  u_min, u_max, v_min, v_max, width, height, margin,
                  strcmp(boundary_by, "material") == 0, &labels, &count)
            : alea_find_surface_label_positions(
                  curves, u_min, u_max, v_min, v_max, width, height, margin,
                  &labels, &count);
    }
    Py_END_ALLOW_THREADS
    int interrupted = restore_sigint(old_sigint);
    int curves_created = curves != NULL;
    free(boundary_ids);
    alea_slice_curves_free(curves);

    if (interrupted) {
        free(labels);
        return NULL;
    }
    if (rc < 0) {
        PyErr_SetString(PyExc_RuntimeError,
                        curves_created ?
                            "Surface label position computation failed" :
                            alea_error());
        return NULL;
    }

    PyObject* result = PyList_New(count);
    for (int i = 0; i < count; i++) {
        PyObject* d = PyDict_New();
        dict_set_new(d, "id", PyLong_FromLong(labels[i].id));
        dict_set_new(d, "px", PyLong_FromLong(labels[i].px));
        dict_set_new(d, "py", PyLong_FromLong(labels[i].py));
        if (labels[i].provenance_group >= 0) {
            PyObject* key = PyUnicode_FromFormat(
                "%d:%d:%d:%d", labels[i].provenance_orientation,
                labels[i].provenance_edge_x, labels[i].provenance_edge_y,
                labels[i].provenance_group);
            if (key) {
                PyDict_SetItemString(d, "coincident_group", key);
                Py_DECREF(key);
            }
        }
        PyList_SET_ITEM(result, i, d);
    }
    free(labels);
    return result;
}

static PyObject* PyAleaSystem_find_surface_labels_on_boundary_map(
        PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    PyObject* origin_obj, *normal_obj, *up_obj, *grid_ids_obj;
    const char* boundary_by = "cell";
    double u_min, u_max, v_min, v_max;
    int width, height, margin = 20, universe_depth = -1;

    static char* kwlist[] = {"origin", "normal", "up", "u_min", "u_max", "v_min", "v_max",
                             "width", "height", "grid_ids", "margin", "boundary_by",
                             "universe_depth", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOddddiiO|isi", kwlist,
            &origin_obj, &normal_obj, &up_obj, &u_min, &u_max, &v_min, &v_max,
            &width, &height, &grid_ids_obj, &margin, &boundary_by,
            &universe_depth)) return NULL;
    if (!self->sys) { PyErr_SetString(PyExc_RuntimeError, "System not initialized"); return NULL; }
    if (width <= 0 || height <= 0) {
        PyErr_SetString(PyExc_ValueError, "width and height must be positive");
        return NULL;
    }
    if (margin < 0) {
        PyErr_SetString(PyExc_ValueError, "margin must not be negative");
        return NULL;
    }

    /* The classifier must reproduce the identity semantics of the grid the
     * caller rendered, or the map would attribute transitions the grid does
     * not show. */
    alea_slice_classify_point_fn classify;
    if (strcmp(boundary_by, "material") == 0) {
        classify = alea_slice_classify_material_at_depth;
    } else if (strcmp(boundary_by, "cell") == 0) {
        classify = alea_slice_classify_cell_at_depth;
    } else {
        PyErr_SetString(PyExc_ValueError, "boundary_by must be 'cell' or 'material'");
        return NULL;
    }

    double ox, oy, oz, nx, ny, nz, ux, uy, uz;
    if (!PyArg_ParseTuple(origin_obj, "ddd", &ox, &oy, &oz)) return NULL;
    if (!PyArg_ParseTuple(normal_obj, "ddd", &nx, &ny, &nz)) return NULL;
    if (!PyArg_ParseTuple(up_obj, "ddd", &ux, &uy, &uz)) return NULL;

    Py_ssize_t size = (Py_ssize_t)width * height;
    int* grid_ids = malloc((size_t)size * sizeof(int));
    if (!grid_ids) return PyErr_NoMemory();
    if (copy_int_values(grid_ids_obj, grid_ids, size, "grid_ids") < 0) {
        free(grid_ids);
        return NULL;
    }

    alea_slice_view_t view;
    alea_slice_view_init(&view, ox, oy, oz, nx, ny, nz, ux, uy, uz,
                              u_min, u_max, v_min, v_max);

    alea_slice_surface_boundary_map_t* map = NULL;
    alea_label_position_t* labels = NULL;
    int count = 0;
    int rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    rc = alea_slice_surface_boundary_map_create(
        self->sys, &view, width, height, grid_ids, classify, &universe_depth, &map);
    if (rc == 0) {
        rc = alea_find_surface_labels_on_boundary_map(map, margin, &labels, &count);
    }
    Py_END_ALLOW_THREADS
    int interrupted = restore_sigint(old_sigint);
    free(grid_ids);

    if (interrupted) {
        free(labels);
        alea_slice_surface_boundary_map_free(map);
        return NULL;
    }
    if (rc != 0) {
        alea_slice_surface_boundary_map_free(map);
        PyErr_SetString(PyExc_RuntimeError, "Surface boundary label computation failed");
        return NULL;
    }

    PyObject* result = PyList_New(count);
    if (!result) {
        free(labels);
        alea_slice_surface_boundary_map_free(map);
        return NULL;
    }
    for (int i = 0; i < count; i++) {
        PyObject* d = PyDict_New();
        dict_set_new(d, "id", PyLong_FromLong(labels[i].id));
        dict_set_new(d, "px", PyLong_FromLong(labels[i].px));
        dict_set_new(d, "py", PyLong_FromLong(labels[i].py));
        dict_set_new(d, "edge_count", PyLong_FromLong(labels[i].pixel_count));
        /* Only an explicit physical crossing group may request presentation
         * stacking.  Sharing the integer label pixel is insufficient: two
         * distinct dense contours can quantize there independently. */
        int found_coincident_group = 0;
        for (int orient = ALEA_SLICE_EDGE_RIGHT;
             orient <= ALEA_SLICE_EDGE_DOWN && !found_coincident_group;
             orient++) {
            size_t groups = alea_slice_surface_boundary_group_count(
                map, labels[i].px, labels[i].py, orient);
            for (size_t group = 0;
                 group < groups && !found_coincident_group; group++) {
                size_t participants = alea_slice_surface_boundary_group_surface_count(
                    map, labels[i].px, labels[i].py, orient, group);
                if (participants < 2) continue;
                for (size_t participant = 0; participant < participants;
                     participant++) {
                    if (alea_slice_surface_boundary_group_surface_id(
                            map, labels[i].px, labels[i].py, orient, group,
                            participant) != labels[i].id)
                        continue;
                    PyObject* key = PyUnicode_FromFormat(
                        "%d:%d:%d:%zu", orient, labels[i].px,
                        labels[i].py, group);
                    if (key) {
                        PyDict_SetItemString(d, "coincident_group", key);
                        Py_DECREF(key);
                    }
                    found_coincident_group = 1;
                    break;
                }
            }
        }
        PyList_SET_ITEM(result, i, d);
    }
    free(labels);
    alea_slice_surface_boundary_map_free(map);
    return result;
}

static PyObject* PyAleaSystem_find_surface_labels_sparse_on_grid(
        PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    PyObject* origin_obj, *normal_obj, *up_obj, *grid_ids_obj;
    const char* boundary_by = "cell";
    double u_min, u_max, v_min, v_max;
    int width, height, margin = 20, universe_depth = -1;
    Py_ssize_t max_queries = 128, max_labels = 256;

    static char* kwlist[] = {
        "origin", "normal", "up", "u_min", "u_max", "v_min", "v_max",
        "width", "height", "grid_ids", "margin", "boundary_by",
        "universe_depth", "max_queries", "max_labels", NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "OOOddddiiO|isinn", kwlist,
            &origin_obj, &normal_obj, &up_obj,
            &u_min, &u_max, &v_min, &v_max,
            &width, &height, &grid_ids_obj, &margin, &boundary_by,
            &universe_depth, &max_queries, &max_labels))
        return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    if (width <= 0 || height <= 0 || margin < 0) {
        PyErr_SetString(PyExc_ValueError,
                        "width and height must be positive and margin non-negative");
        return NULL;
    }
    if (max_queries <= 0 || max_labels <= 0) {
        PyErr_SetString(PyExc_ValueError,
                        "max_queries and max_labels must be positive");
        return NULL;
    }

    alea_slice_classify_point_fn classify;
    if (strcmp(boundary_by, "material") == 0) {
        classify = alea_slice_classify_material_at_depth;
    } else if (strcmp(boundary_by, "cell") == 0) {
        classify = alea_slice_classify_cell_at_depth;
    } else {
        PyErr_SetString(PyExc_ValueError,
                        "boundary_by must be 'cell' or 'material'");
        return NULL;
    }

    double ox, oy, oz, nx, ny, nz, ux, uy, uz;
    if (!PyArg_ParseTuple(origin_obj, "ddd", &ox, &oy, &oz) ||
        !PyArg_ParseTuple(normal_obj, "ddd", &nx, &ny, &nz) ||
        !PyArg_ParseTuple(up_obj, "ddd", &ux, &uy, &uz))
        return NULL;
    if ((size_t)width > SIZE_MAX / (size_t)height ||
        (size_t)width * (size_t)height > (size_t)PY_SSIZE_T_MAX) {
        PyErr_SetString(PyExc_OverflowError, "grid dimensions are too large");
        return NULL;
    }
    Py_ssize_t size = (Py_ssize_t)((size_t)width * (size_t)height);
    int* grid_ids = malloc((size_t)size * sizeof(*grid_ids));
    if (!grid_ids) return PyErr_NoMemory();
    if (copy_int_values(grid_ids_obj, grid_ids, size, "grid_ids") < 0) {
        free(grid_ids);
        return NULL;
    }

    alea_slice_view_t view;
    alea_slice_view_init(&view, ox, oy, oz, nx, ny, nz, ux, uy, uz,
                         u_min, u_max, v_min, v_max);
    alea_label_position_t* labels = NULL;
    int count = 0;
    int rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    rc = alea_find_surface_labels_sparse_on_grid(
        self->sys, &view, width, height, grid_ids, classify, &universe_depth,
        margin, (size_t)max_queries, (size_t)max_labels, &labels, &count);
    Py_END_ALLOW_THREADS
    int interrupted = restore_sigint(old_sigint);
    free(grid_ids);
    if (interrupted) {
        free(labels);
        return NULL;
    }
    if (rc != 0) {
        free(labels);
        PyErr_SetString(PyExc_RuntimeError,
                        "Sparse surface label computation failed");
        return NULL;
    }

    PyObject* result = PyList_New(count);
    if (!result) { free(labels); return NULL; }
    for (int i = 0; i < count; i++) {
        PyObject* d = PyDict_New();
        if (!d) {
            Py_DECREF(result);
            free(labels);
            return NULL;
        }
        dict_set_new(d, "id", PyLong_FromLong(labels[i].id));
        dict_set_new(d, "px", PyLong_FromLong(labels[i].px));
        dict_set_new(d, "py", PyLong_FromLong(labels[i].py));
        dict_set_new(d, "edge_count",
                             PyLong_FromLong(labels[i].pixel_count));
        if (labels[i].provenance_group >= 0) {
            PyObject* key = PyUnicode_FromFormat(
                "%d:%d:%d:%d", labels[i].provenance_orientation,
                labels[i].provenance_edge_x, labels[i].provenance_edge_y,
                labels[i].provenance_group);
            if (key) {
                PyDict_SetItemString(d, "coincident_group", key);
                Py_DECREF(key);
            }
        }
        PyList_SET_ITEM(result, i, d);
    }
    free(labels);
    return result;
}

static PyObject* PyAleaSystem_sparse_surface_label_stats(
        PyAleaSystemObject* self, PyObject* Py_UNUSED(ignored)) {
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    alea_sparse_surface_label_stats_t stats =
        alea_sparse_surface_label_stats_get();
    return Py_BuildValue(
        "{s:n,s:n,s:n,s:n,s:n,s:n,s:n,s:n,s:n,s:n,s:n,s:n,s:n,s:n}",
        "tiles_examined", (Py_ssize_t)stats.tiles_examined,
        "changed_edges", (Py_ssize_t)stats.changed_edges,
        "candidate_edges", (Py_ssize_t)stats.candidate_edges,
        "forward_trace_calls", (Py_ssize_t)stats.forward_trace_calls,
        "reverse_trace_calls", (Py_ssize_t)stats.reverse_trace_calls,
        "accepted_edges", (Py_ssize_t)stats.accepted_edges,
        "observations", (Py_ssize_t)stats.observations,
        "labels", (Py_ssize_t)stats.labels,
        "batch_attempts", (Py_ssize_t)stats.batch_attempts,
        "batch_traces_used", (Py_ssize_t)stats.batch_traces_used,
        "local_provenance_traces_used", (Py_ssize_t)stats.local_provenance_traces_used,
        "local_path_pairs_resolved", (Py_ssize_t)stats.local_path_pairs_resolved,
        "local_path_pairs_same_universe", (Py_ssize_t)stats.local_path_pairs_same_universe,
        "local_path_pairs_same_transform", (Py_ssize_t)stats.local_path_pairs_same_transform,
        "breakpoint_hits", stats.breakpoint_hits,
        "selected_segments", stats.selected_segments);
}
