// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* This file is part of pyalea_binding.c — do NOT compile separately.
 * It is #included from pyalea_binding.c to keep a single translation unit.
 *
 * Contents: Mesh export, mesh sample.
 */

/* ============================================================================
 * Mesh Module
 * ============================================================================ */

static PyObject* PyAleaSystem_mesh_export(PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    const char* filename;
    int nx = 10, ny = 10, nz = 10;
    double x_min = 0, x_max = 0, y_min = 0, y_max = 0, z_min = 0, z_max = 0;
    const char* format_str = "gmsh";
    int void_material_id = 0;
    double auto_pad = 0.01;

    static char* kwlist[] = {"filename", "nx", "ny", "nz",
                             "x_min", "x_max", "y_min", "y_max", "z_min", "z_max",
                             "format", "void_material_id", "auto_pad", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|iiiddddddsid", kwlist,
            &filename, &nx, &ny, &nz,
            &x_min, &x_max, &y_min, &y_max, &z_min, &z_max,
            &format_str, &void_material_id, &auto_pad)) return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);
    cfg.nx = nx; cfg.ny = ny; cfg.nz = nz;
    cfg.x_min = x_min; cfg.x_max = x_max;
    cfg.y_min = y_min; cfg.y_max = y_max;
    cfg.z_min = z_min; cfg.z_max = z_max;
    cfg.void_material_id = void_material_id;
    cfg.auto_pad = auto_pad;

    if (strcmp(format_str, "vtk") == 0) {
        cfg.format = ALEA_MESH_VTK;
    } else {
        cfg.format = ALEA_MESH_GMSH;
    }

    int rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    rc = alea_mesh_export_system(self->sys, &cfg, filename);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) return NULL;

    if (rc < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    Py_RETURN_NONE;
}

static PyObject* PyAleaSystem_mesh_sample(PyAleaSystemObject* self, PyObject* args, PyObject* kwds) {
    int nx = 10, ny = 10, nz = 10;
    double x_min = 0, x_max = 0, y_min = 0, y_max = 0, z_min = 0, z_max = 0;
    int void_material_id = 0;
    double auto_pad = 0.01;
    int as_buffers = 0;

    static char* kwlist[] = {"nx", "ny", "nz",
                             "x_min", "x_max", "y_min", "y_max", "z_min", "z_max",
                             "void_material_id", "auto_pad", "_as_buffers", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|iiiddddddidp", kwlist,
            &nx, &ny, &nz,
            &x_min, &x_max, &y_min, &y_max, &z_min, &z_max,
            &void_material_id, &auto_pad, &as_buffers)) return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }

    alea_mesh_config_t cfg;
    alea_mesh_config_init(&cfg);
    cfg.nx = nx; cfg.ny = ny; cfg.nz = nz;
    cfg.x_min = x_min; cfg.x_max = x_max;
    cfg.y_min = y_min; cfg.y_max = y_max;
    cfg.z_min = z_min; cfg.z_max = z_max;
    cfg.void_material_id = void_material_id;
    cfg.auto_pad = auto_pad;
    /* The legacy entry point returns only these two arrays. */
    cfg.fields = ALEA_MESH_FIELD_MATERIAL_ID | ALEA_MESH_FIELD_CELL_ID;

    alea_mesh_result_t* mesh;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    mesh = alea_mesh_sample(self->sys, &cfg);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) return NULL;

    if (!mesh) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }

    /* Build result dict */
    Py_ssize_t total = (Py_ssize_t)mesh->nx * mesh->ny * mesh->nz;

    PyObject* mat_ids = as_buffers
        ? grid_array_from_owned_data(mesh->material_ids, total, NPY_INT)
        : PyList_New(total);
    PyObject* cell_ids = as_buffers
        ? grid_array_from_owned_data(mesh->cell_ids, total, NPY_INT)
        : PyList_New(total);
    if (!mat_ids || !cell_ids) {
        Py_XDECREF(mat_ids); Py_XDECREF(cell_ids);
        if (as_buffers) {
            mesh->material_ids = NULL;
            mesh->cell_ids = NULL;
        }
        alea_mesh_result_free(mesh);
        return NULL;
    }
    if (!as_buffers) {
        for (Py_ssize_t i = 0; i < total; i++) {
            PyList_SET_ITEM(mat_ids, i, PyLong_FromLong(mesh->material_ids[i]));
            PyList_SET_ITEM(cell_ids, i, PyLong_FromLong(mesh->cell_ids[i]));
        }
    }

    /* Node positions */
    PyObject* x_nodes = PyList_New(mesh->nx + 1);
    PyObject* y_nodes = PyList_New(mesh->ny + 1);
    PyObject* z_nodes = PyList_New(mesh->nz + 1);
    if (!x_nodes || !y_nodes || !z_nodes) {
        Py_XDECREF(mat_ids); Py_XDECREF(cell_ids);
        Py_XDECREF(x_nodes); Py_XDECREF(y_nodes); Py_XDECREF(z_nodes);
        if (as_buffers) {
            mesh->material_ids = NULL;
            mesh->cell_ids = NULL;
        }
        alea_mesh_result_free(mesh);
        return NULL;
    }
    for (int i = 0; i <= mesh->nx; i++)
        PyList_SET_ITEM(x_nodes, i, PyFloat_FromDouble(mesh->x_nodes[i]));
    for (int i = 0; i <= mesh->ny; i++)
        PyList_SET_ITEM(y_nodes, i, PyFloat_FromDouble(mesh->y_nodes[i]));
    for (int i = 0; i <= mesh->nz; i++)
        PyList_SET_ITEM(z_nodes, i, PyFloat_FromDouble(mesh->z_nodes[i]));

    PyObject* result = Py_BuildValue("{s:N, s:N, s:N, s:N, s:N, s:i, s:i, s:i}",
        "material_ids", mat_ids,
        "cell_ids", cell_ids,
        "x_nodes", x_nodes,
        "y_nodes", y_nodes,
        "z_nodes", z_nodes,
        "nx", mesh->nx,
        "ny", mesh->ny,
        "nz", mesh->nz);

    if (as_buffers) {
        /* Ownership moved into the NumPy arrays' base capsules. */
        mesh->material_ids = NULL;
        mesh->cell_ids = NULL;
    }
    alea_mesh_result_free(mesh);
    return result;
}

/* ============================================================================
 * Configurable structured and adaptive mesh APIs
 * ============================================================================ */

typedef struct {
    double *x;
    double *y;
    double *z;
    double *ray_points;
} mesh_node_buffers_t;

static void mesh_node_buffers_free(mesh_node_buffers_t *nodes) {
    if (!nodes) return;
    free(nodes->x);
    free(nodes->y);
    free(nodes->z);
    free(nodes->ray_points);
    memset(nodes, 0, sizeof(*nodes));
}

static int mesh_option_int(PyObject *options, const char *name, int *value) {
    PyObject *item = PyDict_GetItemString(options, name);
    if (!item) return 0;
    long parsed = PyLong_AsLong(item);
    if (parsed == -1 && PyErr_Occurred()) return -1;
    if (parsed < INT_MIN || parsed > INT_MAX) {
        PyErr_Format(PyExc_OverflowError, "%s is outside the C int range", name);
        return -1;
    }
    *value = (int)parsed;
    return 0;
}

static int mesh_option_uint32(PyObject *options, const char *name, uint32_t *value) {
    PyObject *item = PyDict_GetItemString(options, name);
    if (!item) return 0;
    unsigned long parsed = PyLong_AsUnsignedLong(item);
    if (PyErr_Occurred()) return -1;
    if (parsed > UINT32_MAX) {
        PyErr_Format(PyExc_OverflowError, "%s is outside the uint32 range", name);
        return -1;
    }
    *value = (uint32_t)parsed;
    return 0;
}

static int mesh_option_uint64(PyObject *options, const char *name, uint64_t *value) {
    PyObject *item = PyDict_GetItemString(options, name);
    if (!item) return 0;
    unsigned long long parsed = PyLong_AsUnsignedLongLong(item);
    if (PyErr_Occurred()) return -1;
    *value = (uint64_t)parsed;
    return 0;
}

static int mesh_option_size(PyObject *options, const char *name, size_t *value) {
    uint64_t parsed = 0;
    if (!PyDict_GetItemString(options, name)) return 0;
    if (mesh_option_uint64(options, name, &parsed) < 0) return -1;
    if (parsed > SIZE_MAX) {
        PyErr_Format(PyExc_OverflowError, "%s is outside the size_t range", name);
        return -1;
    }
    *value = (size_t)parsed;
    return 0;
}

static int mesh_option_double(PyObject *options, const char *name, double *value) {
    PyObject *item = PyDict_GetItemString(options, name);
    if (!item) return 0;
    double parsed = PyFloat_AsDouble(item);
    if (parsed == -1.0 && PyErr_Occurred()) return -1;
    *value = parsed;
    return 0;
}

static int mesh_option_bool(PyObject *options, const char *name, int *value) {
    PyObject *item = PyDict_GetItemString(options, name);
    if (!item) return 0;
    int parsed = PyObject_IsTrue(item);
    if (parsed < 0) return -1;
    *value = parsed;
    return 0;
}

static int mesh_parse_sampling_mode(PyObject *options,
                                    alea_mesh_sampling_mode_t *mode) {
    PyObject *item = PyDict_GetItemString(options, "sampling_mode");
    if (!item) return 0;
    const char *value = PyUnicode_AsUTF8(item);
    if (!value) return -1;
    if (strcmp(value, "center") == 0) *mode = ALEA_MESH_SAMPLE_CENTER;
    else if (strcmp(value, "corners") == 0) *mode = ALEA_MESH_SAMPLE_CORNERS;
    else if (strcmp(value, "subcell") == 0) *mode = ALEA_MESH_SAMPLE_SUBCELL;
    else if (strcmp(value, "stratified") == 0) *mode = ALEA_MESH_SAMPLE_STRATIFIED;
    else if (strcmp(value, "adaptive") == 0) *mode = ALEA_MESH_SAMPLE_ADAPTIVE;
    else if (strcmp(value, "ray") == 0) *mode = ALEA_MESH_SAMPLE_RAY;
    else {
        PyErr_Format(PyExc_ValueError, "unsupported mesh sampling mode: %s", value);
        return -1;
    }
    return 0;
}

static int mesh_parse_format(PyObject *options, alea_mesh_format_t *format) {
    PyObject *item = PyDict_GetItemString(options, "format");
    if (!item) return 0;
    const char *value = PyUnicode_AsUTF8(item);
    if (!value) return -1;
    if (strcmp(value, "gmsh") == 0) *format = ALEA_MESH_GMSH;
    else if (strcmp(value, "vtk") == 0) *format = ALEA_MESH_VTK;
    else {
        PyErr_Format(PyExc_ValueError, "unsupported mesh format: %s", value);
        return -1;
    }
    return 0;
}

static int mesh_parse_bounds_mode(PyObject *options,
                                  alea_mesh_bounds_mode_t *mode) {
    PyObject *item = PyDict_GetItemString(options, "bounds_mode");
    if (!item) return 0;
    const char *value = PyUnicode_AsUTF8(item);
    if (!value) return -1;
    if (strcmp(value, "legacy") == 0) *mode = ALEA_MESH_BOUNDS_LEGACY;
    else if (strcmp(value, "auto") == 0) *mode = ALEA_MESH_BOUNDS_AUTO;
    else if (strcmp(value, "explicit") == 0) *mode = ALEA_MESH_BOUNDS_EXPLICIT;
    else {
        PyErr_Format(PyExc_ValueError, "unsupported mesh bounds mode: %s", value);
        return -1;
    }
    return 0;
}

static int mesh_parse_ray_origin_mode(PyObject *options,
                                      alea_mesh_ray_origin_mode_t *mode) {
    PyObject *item = PyDict_GetItemString(options, "ray_origin_mode");
    if (!item) return 0;
    const char *value = PyUnicode_AsUTF8(item);
    if (!value) return -1;
    if (strcmp(value, "grid") == 0) *mode = ALEA_MESH_RAY_ORIGINS_GRID;
    else if (strcmp(value, "sobol") == 0) *mode = ALEA_MESH_RAY_ORIGINS_SOBOL;
    else if (strcmp(value, "custom") == 0) *mode = ALEA_MESH_RAY_ORIGINS_CUSTOM;
    else {
        PyErr_Format(PyExc_ValueError,
                     "unsupported ray origin mode: %s", value);
        return -1;
    }
    return 0;
}

static int mesh_parse_nodes(PyObject *options, const char *name, int expected,
                            double **storage, const double **target) {
    PyObject *item = PyDict_GetItemString(options, name);
    if (!item || item == Py_None) return 0;
    PyObject *sequence = PySequence_Fast(item, "mesh nodes must be a sequence");
    if (!sequence) return -1;
    Py_ssize_t count = PySequence_Fast_GET_SIZE(sequence);
    if (count != (Py_ssize_t)expected + 1) {
        PyErr_Format(PyExc_ValueError, "%s must contain %d values", name, expected + 1);
        Py_DECREF(sequence);
        return -1;
    }
    double *values = malloc((size_t)count * sizeof(*values));
    if (!values) {
        Py_DECREF(sequence);
        PyErr_NoMemory();
        return -1;
    }
    for (Py_ssize_t i = 0; i < count; i++) {
        values[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(sequence, i));
        if (values[i] == -1.0 && PyErr_Occurred()) {
            free(values);
            Py_DECREF(sequence);
            return -1;
        }
    }
    Py_DECREF(sequence);
    *storage = values;
    *target = values;
    return 0;
}

static int mesh_parse_ray_points(PyObject *options, mesh_node_buffers_t *buffers,
                                 alea_mesh_config_t *cfg) {
    PyObject *item = PyDict_GetItemString(options, "ray_points");
    if (!item || item == Py_None) return 0;
    PyObject *points = PySequence_Fast(item,
        "ray_points must be a sequence of (u, v) pairs");
    if (!points) return -1;
    Py_ssize_t count = PySequence_Fast_GET_SIZE(points);
    if (count <= 0 || (size_t)count > UINT32_MAX) {
        Py_DECREF(points);
        PyErr_SetString(PyExc_ValueError, "ray_points must not be empty");
        return -1;
    }
    double *values = malloc((size_t)count * 2u * sizeof(*values));
    if (!values) {
        Py_DECREF(points);
        PyErr_NoMemory();
        return -1;
    }
    for (Py_ssize_t i = 0; i < count; i++) {
        PyObject *pair = PySequence_Fast(PySequence_Fast_GET_ITEM(points, i),
                                         "each ray point must be a (u, v) pair");
        if (!pair || PySequence_Fast_GET_SIZE(pair) != 2) {
            Py_XDECREF(pair);
            free(values);
            Py_DECREF(points);
            if (!PyErr_Occurred())
                PyErr_SetString(PyExc_ValueError,
                                "each ray point must contain two values");
            return -1;
        }
        values[2 * i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(pair, 0));
        values[2 * i + 1] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(pair, 1));
        Py_DECREF(pair);
        if (PyErr_Occurred()) {
            free(values);
            Py_DECREF(points);
            return -1;
        }
    }
    Py_DECREF(points);
    buffers->ray_points = values;
    cfg->ray_points = values;
    cfg->ray_point_count = (uint32_t)count;
    return 0;
}

static int mesh_config_from_dict(PyObject *options, alea_mesh_config_t *cfg,
                                 mesh_node_buffers_t *nodes) {
    if (!PyDict_Check(options)) {
        PyErr_SetString(PyExc_TypeError, "mesh options must be a dict");
        return -1;
    }
    if (mesh_option_int(options, "nx", &cfg->nx) < 0 ||
        mesh_option_int(options, "ny", &cfg->ny) < 0 ||
        mesh_option_int(options, "nz", &cfg->nz) < 0 ||
        mesh_option_double(options, "x_min", &cfg->x_min) < 0 ||
        mesh_option_double(options, "x_max", &cfg->x_max) < 0 ||
        mesh_option_double(options, "y_min", &cfg->y_min) < 0 ||
        mesh_option_double(options, "y_max", &cfg->y_max) < 0 ||
        mesh_option_double(options, "z_min", &cfg->z_min) < 0 ||
        mesh_option_double(options, "z_max", &cfg->z_max) < 0 ||
        mesh_option_int(options, "void_material_id", &cfg->void_material_id) < 0 ||
        mesh_option_double(options, "auto_pad", &cfg->auto_pad) < 0 ||
        mesh_parse_sampling_mode(options, &cfg->sampling_mode) < 0 ||
        mesh_option_int(options, "subsamples_per_axis", &cfg->subsamples_per_axis) < 0 ||
        mesh_option_double(options, "mixed_threshold", &cfg->mixed_threshold) < 0 ||
        mesh_option_double(options, "target_error", &cfg->target_error) < 0 ||
        mesh_option_int(options, "max_refine_depth", &cfg->max_refine_depth) < 0 ||
        mesh_option_uint32(options, "max_samples_per_voxel", &cfg->max_samples_per_voxel) < 0 ||
        mesh_option_uint64(options, "max_total_samples", &cfg->max_total_samples) < 0 ||
        mesh_option_uint64(options, "sampling_seed", &cfg->sampling_seed) < 0 ||
        mesh_option_int(options, "workers", &cfg->workers) < 0 ||
        mesh_option_int(options, "ray_grid_u", &cfg->ray_grid_u) < 0 ||
        mesh_option_int(options, "ray_grid_v", &cfg->ray_grid_v) < 0 ||
        mesh_parse_ray_origin_mode(options, &cfg->ray_origin_mode) < 0 ||
        mesh_option_uint32(options, "ray_samples", &cfg->ray_samples) < 0 ||
        mesh_option_uint32(options, "fields", &cfg->fields) < 0 ||
        mesh_parse_bounds_mode(options, &cfg->bounds_mode) < 0 ||
        mesh_parse_format(options, &cfg->format) < 0)
        return -1;

    uint32_t ray_directions = cfg->ray_directions;
    if (mesh_option_uint32(options, "ray_directions", &ray_directions) < 0)
        return -1;
    if (ray_directions > UINT8_MAX) {
        PyErr_SetString(PyExc_ValueError, "ray_directions exceeds uint8 range");
        return -1;
    }
    cfg->ray_directions = (uint8_t)ray_directions;

    if (mesh_parse_nodes(options, "x_nodes", cfg->nx, &nodes->x, &cfg->x_nodes) < 0 ||
        mesh_parse_nodes(options, "y_nodes", cfg->ny, &nodes->y, &cfg->y_nodes) < 0 ||
        mesh_parse_nodes(options, "z_nodes", cfg->nz, &nodes->z, &cfg->z_nodes) < 0 ||
        mesh_parse_ray_points(options, nodes, cfg) < 0)
        return -1;
    return 0;
}

static const char *mesh_sampling_mode_name(alea_mesh_sampling_mode_t mode) {
    switch (mode) {
        case ALEA_MESH_SAMPLE_CENTER: return "center";
        case ALEA_MESH_SAMPLE_CORNERS: return "corners";
        case ALEA_MESH_SAMPLE_SUBCELL: return "subcell";
        case ALEA_MESH_SAMPLE_STRATIFIED: return "stratified";
        case ALEA_MESH_SAMPLE_ADAPTIVE: return "adaptive";
        case ALEA_MESH_SAMPLE_RAY: return "ray";
        default: return "unknown";
    }
}

static const char *mesh_bounds_source_name(alea_mesh_bounds_source_t source) {
    switch (source) {
        case ALEA_MESH_BOUNDS_SOURCE_EXPLICIT: return "explicit";
        case ALEA_MESH_BOUNDS_SOURCE_CUSTOM_NODES: return "custom_nodes";
        case ALEA_MESH_BOUNDS_SOURCE_INFERRED_ROOT_AABB: return "inferred_root_aabb";
        default: return "unknown";
    }
}

static int mesh_dict_take_array(PyObject *dict, const char *name, void **data,
                                Py_ssize_t count, int typenum) {
    if (!*data) return 0;
    PyObject *array = grid_array_from_owned_data(*data, count, typenum);
    if (!array) return -1;
    *data = NULL;
    int rc = PyDict_SetItemString(dict, name, array);
    Py_DECREF(array);
    return rc;
}

static PyObject *mesh_nodes_list(const double *values, int count) {
    PyObject *list = PyList_New(count);
    if (!list) return NULL;
    for (int i = 0; i < count; i++) {
        PyObject *value = PyFloat_FromDouble(values[i]);
        if (!value) {
            Py_DECREF(list);
            return NULL;
        }
        PyList_SET_ITEM(list, i, value);
    }
    return list;
}

static PyObject *mesh_result_to_dict(alea_mesh_result_t *mesh) {
    PyObject *result = PyDict_New();
    if (!result) return NULL;
    Py_ssize_t total = (Py_ssize_t)mesh->nx * mesh->ny * mesh->nz;

#define MESH_SET_VALUE(KEY, EXPR) do { \
    PyObject *_value = (EXPR); \
    if (!_value || PyDict_SetItemString(result, (KEY), _value) < 0) { \
        Py_XDECREF(_value); Py_DECREF(result); return NULL; \
    } \
    Py_DECREF(_value); \
} while (0)

    MESH_SET_VALUE("nx", PyLong_FromLong(mesh->nx));
    MESH_SET_VALUE("ny", PyLong_FromLong(mesh->ny));
    MESH_SET_VALUE("nz", PyLong_FromLong(mesh->nz));
    MESH_SET_VALUE("fields", PyLong_FromUnsignedLong(mesh->fields));
    MESH_SET_VALUE("bounds_source", PyUnicode_FromString(mesh_bounds_source_name(mesh->bounds_source)));
    MESH_SET_VALUE("bounds_padding", PyFloat_FromDouble(mesh->bounds_padding));
    MESH_SET_VALUE("sampling_mode", PyUnicode_FromString(mesh_sampling_mode_name(mesh->sampling_mode)));
    MESH_SET_VALUE("sampling_seed", PyLong_FromUnsignedLongLong(mesh->sampling_seed));
    MESH_SET_VALUE("target_error", PyFloat_FromDouble(mesh->target_error));
    MESH_SET_VALUE("mixed_count", PyLong_FromLong(mesh->mixed_count));
    MESH_SET_VALUE("x_nodes", mesh_nodes_list(mesh->x_nodes, mesh->nx + 1));
    MESH_SET_VALUE("y_nodes", mesh_nodes_list(mesh->y_nodes, mesh->ny + 1));
    MESH_SET_VALUE("z_nodes", mesh_nodes_list(mesh->z_nodes, mesh->nz + 1));

    if (mesh_dict_take_array(result, "material_ids", (void **)&mesh->material_ids,
                             total, NPY_INT) < 0 ||
        mesh_dict_take_array(result, "cell_ids", (void **)&mesh->cell_ids,
                             total, NPY_INT) < 0 ||
        mesh_dict_take_array(result, "mixed_flags", (void **)&mesh->mixed_flags,
                             total, NPY_UINT8) < 0 ||
        mesh_dict_take_array(result, "dominant_fractions", (void **)&mesh->dominant_fractions,
                             total, NPY_DOUBLE) < 0 ||
        mesh_dict_take_array(result, "estimated_errors", (void **)&mesh->estimated_errors,
                             total, NPY_DOUBLE) < 0 ||
        mesh_dict_take_array(result, "sample_counts", (void **)&mesh->sample_counts,
                             total, NPY_UINT32) < 0 ||
        mesh_dict_take_array(result, "tie_flags", (void **)&mesh->tie_flags,
                             total, NPY_UINT8) < 0 ||
        mesh_dict_take_array(result, "refinement_flags", (void **)&mesh->refinement_flags,
                             total, NPY_UINT8) < 0 ||
        mesh_dict_take_array(result, "unique_materials", (void **)&mesh->unique_materials,
                             mesh->num_materials, NPY_INT) < 0) {
        Py_DECREF(result);
        return NULL;
    }

    if (mesh->fraction_spans) {
        uint32_t *offsets = malloc(((size_t)total + 1) * sizeof(*offsets));
        int *materials = malloc(mesh->fraction_count * sizeof(*materials));
        double *fractions = malloc(mesh->fraction_count * sizeof(*fractions));
        if (!offsets || (mesh->fraction_count && (!materials || !fractions))) {
            free(offsets); free(materials); free(fractions);
            Py_DECREF(result);
            PyErr_NoMemory();
            return NULL;
        }
        for (Py_ssize_t i = 0; i < total; i++) offsets[i] = mesh->fraction_spans[i].offset;
        offsets[total] = (uint32_t)mesh->fraction_count;
        for (size_t i = 0; i < mesh->fraction_count; i++) {
            materials[i] = mesh->fractions[i].material_id;
            fractions[i] = mesh->fractions[i].fraction;
        }
        if (mesh_dict_take_array(result, "fraction_offsets", (void **)&offsets,
                                 total + 1, NPY_UINT32) < 0 ||
            mesh_dict_take_array(result, "fraction_material_ids", (void **)&materials,
                                 (Py_ssize_t)mesh->fraction_count, NPY_INT) < 0 ||
            mesh_dict_take_array(result, "fraction_values", (void **)&fractions,
                                 (Py_ssize_t)mesh->fraction_count, NPY_DOUBLE) < 0) {
            free(offsets); free(materials); free(fractions);
            Py_DECREF(result);
            return NULL;
        }
    }
    if (mesh->cell_fraction_spans) {
        uint32_t *offsets = malloc(((size_t)total + 1) * sizeof(*offsets));
        int *cells = malloc(mesh->cell_fraction_count * sizeof(*cells));
        int *materials = malloc(mesh->cell_fraction_count * sizeof(*materials));
        double *fractions = malloc(mesh->cell_fraction_count * sizeof(*fractions));
        if (!offsets || (mesh->cell_fraction_count &&
                         (!cells || !materials || !fractions))) {
            free(offsets); free(cells); free(materials); free(fractions);
            Py_DECREF(result);
            PyErr_NoMemory();
            return NULL;
        }
        for (Py_ssize_t i = 0; i < total; i++)
            offsets[i] = mesh->cell_fraction_spans[i].offset;
        offsets[total] = (uint32_t)mesh->cell_fraction_count;
        for (size_t i = 0; i < mesh->cell_fraction_count; i++) {
            cells[i] = mesh->cell_fractions[i].cell_id;
            materials[i] = mesh->cell_fractions[i].material_id;
            fractions[i] = mesh->cell_fractions[i].fraction;
        }
        if (mesh_dict_take_array(result, "cell_fraction_offsets", (void **)&offsets,
                                 total + 1, NPY_UINT32) < 0 ||
            mesh_dict_take_array(result, "cell_fraction_cell_ids", (void **)&cells,
                                 (Py_ssize_t)mesh->cell_fraction_count, NPY_INT) < 0 ||
            mesh_dict_take_array(result, "cell_fraction_material_ids",
                                 (void **)&materials,
                                 (Py_ssize_t)mesh->cell_fraction_count, NPY_INT) < 0 ||
            mesh_dict_take_array(result, "cell_fraction_values", (void **)&fractions,
                                 (Py_ssize_t)mesh->cell_fraction_count,
                                 NPY_DOUBLE) < 0) {
            free(offsets); free(cells); free(materials); free(fractions);
            Py_DECREF(result);
            return NULL;
        }
    }
#undef MESH_SET_VALUE
    return result;
}

static PyObject *PyAleaSystem_mesh_sample_configured(
        PyAleaSystemObject *self, PyObject *options) {
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    alea_mesh_config_t cfg;
    mesh_node_buffers_t nodes = {0};
    alea_mesh_config_init(&cfg);
    if (mesh_config_from_dict(options, &cfg, &nodes) < 0) {
        mesh_node_buffers_free(&nodes);
        return NULL;
    }
    alea_mesh_result_t *mesh;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    mesh = alea_mesh_sample(self->sys, &cfg);
    Py_END_ALLOW_THREADS
    mesh_node_buffers_free(&nodes);
    if (restore_sigint(old_sigint)) {
        alea_mesh_result_free(mesh);
        return NULL;
    }
    if (!mesh) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }
    PyObject *result = mesh_result_to_dict(mesh);
    alea_mesh_result_free(mesh);
    return result;
}

static uint32_t mesh_export_required_fields(uint32_t export_fields) {
    uint32_t fields = ALEA_MESH_FIELD_MATERIAL_ID | ALEA_MESH_FIELD_CELL_ID;
    if (export_fields & ALEA_MESH_EXPORT_MIXED_FLAG) fields |= ALEA_MESH_FIELD_MIXED_FLAG;
    if (export_fields & ALEA_MESH_EXPORT_DOMINANT_FRACTION) fields |= ALEA_MESH_FIELD_DOMINANT_FRACTION;
    if (export_fields & ALEA_MESH_EXPORT_TIE_FLAG) fields |= ALEA_MESH_FIELD_TIE_FLAG;
    if (export_fields & ALEA_MESH_EXPORT_SAMPLE_COUNT) fields |= ALEA_MESH_FIELD_SAMPLE_COUNT;
    if (export_fields & ALEA_MESH_EXPORT_MATERIAL_FRACTIONS) fields |= ALEA_MESH_FIELD_SAMPLED_FRACTIONS;
    if (export_fields & ALEA_MESH_EXPORT_ESTIMATED_ERROR) fields |= ALEA_MESH_FIELD_ESTIMATED_ERROR;
    if (export_fields & ALEA_MESH_EXPORT_REFINEMENT_FLAG) fields |= ALEA_MESH_FIELD_REFINEMENT_FLAG;
    return fields;
}

static PyObject *PyAleaSystem_mesh_export_configured(
        PyAleaSystemObject *self, PyObject *args) {
    const char *filename;
    PyObject *options;
    if (!PyArg_ParseTuple(args, "sO!:mesh_export_configured", &filename,
                          &PyDict_Type, &options)) return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    alea_mesh_config_t cfg;
    alea_mesh_export_options_t export_options;
    mesh_node_buffers_t nodes = {0};
    alea_mesh_config_init(&cfg);
    alea_mesh_export_options_init(&export_options);
    if (mesh_config_from_dict(options, &cfg, &nodes) < 0 ||
        mesh_option_uint32(options, "export_fields", &export_options.fields) < 0 ||
        mesh_option_int(options, "max_fraction_materials",
                        &export_options.max_fraction_materials) < 0) {
        mesh_node_buffers_free(&nodes);
        return NULL;
    }
    cfg.fields = mesh_export_required_fields(export_options.fields);
    alea_mesh_result_t *mesh;
    int rc = -1;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    mesh = alea_mesh_sample(self->sys, &cfg);
    if (mesh) rc = alea_mesh_export_ex(mesh, cfg.format, filename, &export_options);
    Py_END_ALLOW_THREADS
    mesh_node_buffers_free(&nodes);
    if (restore_sigint(old_sigint)) {
        alea_mesh_result_free(mesh);
        return NULL;
    }
    alea_mesh_result_free(mesh);
    if (rc < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }
    Py_RETURN_NONE;
}

static PyObject *adaptive_result_to_dict(const alea_adaptive_grid_result_t *grid) {
    npy_intp one[1] = {(npy_intp)grid->cell_count};
    npy_intp offsets_dim[1] = {(npy_intp)grid->cell_count + 1};
    npy_intp fraction_dim[1] = {(npy_intp)grid->fraction_count};
    npy_intp cell_fraction_dim[1] = {(npy_intp)grid->cell_fraction_count};
    npy_intp children_dims[2] = {(npy_intp)grid->cell_count, 8};
    npy_intp bounds_dims[2] = {(npy_intp)grid->cell_count, 6};
    PyObject *result = PyDict_New();
    PyObject *ids = NULL, *parents = NULL, *children = NULL, *levels = NULL;
    PyObject *leaves = NULL, *flags = NULL, *bounds = NULL, *materials = NULL;
    PyObject *cells = NULL, *mixed = NULL, *ties = NULL, *refinement = NULL;
    PyObject *dominant = NULL, *errors = NULL, *samples = NULL;
    PyObject *fraction_offsets = NULL, *fraction_materials = NULL;
    PyObject *fraction_values = NULL, *cell_fraction_offsets = NULL;
    PyObject *cell_fraction_cells = NULL, *cell_fraction_materials = NULL;
    PyObject *cell_fraction_values = NULL;
    if (!result ||
        !(ids = PyArray_SimpleNew(1, one, NPY_UINT64)) ||
        !(parents = PyArray_SimpleNew(1, one, NPY_UINT64)) ||
        !(children = PyArray_SimpleNew(2, children_dims, NPY_UINT64)) ||
        !(levels = PyArray_SimpleNew(1, one, NPY_UINT32)) ||
        !(leaves = PyArray_SimpleNew(1, one, NPY_UINT8)) ||
        !(flags = PyArray_SimpleNew(1, one, NPY_UINT8)) ||
        !(bounds = PyArray_SimpleNew(2, bounds_dims, NPY_DOUBLE)) ||
        !(materials = PyArray_SimpleNew(1, one, NPY_INT)) ||
        !(cells = PyArray_SimpleNew(1, one, NPY_INT)) ||
        !(mixed = PyArray_SimpleNew(1, one, NPY_UINT8)) ||
        !(ties = PyArray_SimpleNew(1, one, NPY_UINT8)) ||
        !(refinement = PyArray_SimpleNew(1, one, NPY_UINT8)) ||
        !(dominant = PyArray_SimpleNew(1, one, NPY_DOUBLE)) ||
        !(errors = PyArray_SimpleNew(1, one, NPY_DOUBLE)) ||
        !(samples = PyArray_SimpleNew(1, one, NPY_UINT32)) ||
        !(fraction_offsets = PyArray_SimpleNew(1, offsets_dim, NPY_UINT32)) ||
        !(fraction_materials = PyArray_SimpleNew(1, fraction_dim, NPY_INT)) ||
        !(fraction_values = PyArray_SimpleNew(1, fraction_dim, NPY_DOUBLE)) ||
        !(cell_fraction_offsets = PyArray_SimpleNew(1, offsets_dim, NPY_UINT32)) ||
        !(cell_fraction_cells = PyArray_SimpleNew(1, cell_fraction_dim, NPY_INT)) ||
        !(cell_fraction_materials = PyArray_SimpleNew(1, cell_fraction_dim, NPY_INT)) ||
        !(cell_fraction_values = PyArray_SimpleNew(1, cell_fraction_dim, NPY_DOUBLE)))
        goto error;

    for (size_t i = 0; i < grid->cell_count; i++) {
        const alea_adaptive_grid_cell_t *cell = &grid->cells[i];
        ((uint64_t *)PyArray_DATA((PyArrayObject *)ids))[i] = cell->id;
        ((uint64_t *)PyArray_DATA((PyArrayObject *)parents))[i] = cell->parent_id;
        memcpy(&((uint64_t *)PyArray_DATA((PyArrayObject *)children))[i * 8],
               cell->child_ids, 8 * sizeof(uint64_t));
        ((uint32_t *)PyArray_DATA((PyArrayObject *)levels))[i] = cell->level;
        ((uint8_t *)PyArray_DATA((PyArrayObject *)leaves))[i] = cell->is_leaf;
        ((uint8_t *)PyArray_DATA((PyArrayObject *)flags))[i] = cell->flags;
        double *box = &((double *)PyArray_DATA((PyArrayObject *)bounds))[i * 6];
        box[0] = cell->x_min; box[1] = cell->x_max;
        box[2] = cell->y_min; box[3] = cell->y_max;
        box[4] = cell->z_min; box[5] = cell->z_max;
        ((int *)PyArray_DATA((PyArrayObject *)materials))[i] = cell->material_id;
        ((int *)PyArray_DATA((PyArrayObject *)cells))[i] = cell->cell_id;
        ((uint8_t *)PyArray_DATA((PyArrayObject *)mixed))[i] = cell->mixed;
        ((uint8_t *)PyArray_DATA((PyArrayObject *)ties))[i] = cell->tie_flags;
        ((uint8_t *)PyArray_DATA((PyArrayObject *)refinement))[i] = cell->refinement_flags;
        ((double *)PyArray_DATA((PyArrayObject *)dominant))[i] = cell->dominant_fraction;
        ((double *)PyArray_DATA((PyArrayObject *)errors))[i] = cell->estimated_error;
        ((uint32_t *)PyArray_DATA((PyArrayObject *)samples))[i] = cell->sample_count;
        ((uint32_t *)PyArray_DATA((PyArrayObject *)fraction_offsets))[i] =
            cell->fraction_span.offset;
        ((uint32_t *)PyArray_DATA((PyArrayObject *)cell_fraction_offsets))[i] =
            cell->cell_fraction_span.offset;
    }
    ((uint32_t *)PyArray_DATA((PyArrayObject *)fraction_offsets))[grid->cell_count] =
        (uint32_t)grid->fraction_count;
    ((uint32_t *)PyArray_DATA((PyArrayObject *)cell_fraction_offsets))[grid->cell_count] =
        (uint32_t)grid->cell_fraction_count;
    for (size_t i = 0; i < grid->fraction_count; i++) {
        ((int *)PyArray_DATA((PyArrayObject *)fraction_materials))[i] =
            grid->fractions[i].material_id;
        ((double *)PyArray_DATA((PyArrayObject *)fraction_values))[i] =
            grid->fractions[i].fraction;
    }
    for (size_t i = 0; i < grid->cell_fraction_count; i++) {
        ((int *)PyArray_DATA((PyArrayObject *)cell_fraction_cells))[i] =
            grid->cell_fractions[i].cell_id;
        ((int *)PyArray_DATA((PyArrayObject *)cell_fraction_materials))[i] =
            grid->cell_fractions[i].material_id;
        ((double *)PyArray_DATA((PyArrayObject *)cell_fraction_values))[i] =
            grid->cell_fractions[i].fraction;
    }

#define ADAPTIVE_SET(KEY, VALUE) do { \
    if (PyDict_SetItemString(result, (KEY), (VALUE)) < 0) goto error; \
    Py_CLEAR(VALUE); \
} while (0)
    ADAPTIVE_SET("ids", ids); ADAPTIVE_SET("parent_ids", parents);
    ADAPTIVE_SET("child_ids", children); ADAPTIVE_SET("levels", levels);
    ADAPTIVE_SET("is_leaf", leaves); ADAPTIVE_SET("flags", flags);
    ADAPTIVE_SET("bounds", bounds); ADAPTIVE_SET("material_ids", materials);
    ADAPTIVE_SET("cell_ids", cells); ADAPTIVE_SET("mixed_flags", mixed);
    ADAPTIVE_SET("tie_flags", ties); ADAPTIVE_SET("refinement_flags", refinement);
    ADAPTIVE_SET("dominant_fractions", dominant); ADAPTIVE_SET("estimated_errors", errors);
    ADAPTIVE_SET("sample_counts", samples);
    ADAPTIVE_SET("fraction_offsets", fraction_offsets);
    ADAPTIVE_SET("fraction_material_ids", fraction_materials);
    ADAPTIVE_SET("fraction_values", fraction_values);
    ADAPTIVE_SET("cell_fraction_offsets", cell_fraction_offsets);
    ADAPTIVE_SET("cell_fraction_cell_ids", cell_fraction_cells);
    ADAPTIVE_SET("cell_fraction_material_ids", cell_fraction_materials);
    ADAPTIVE_SET("cell_fraction_values", cell_fraction_values);
#undef ADAPTIVE_SET
    PyObject *value;
#define ADAPTIVE_META(KEY, EXPR) do { \
    value = (EXPR); \
    if (!value || PyDict_SetItemString(result, (KEY), value) < 0) { Py_XDECREF(value); goto error; } \
    Py_DECREF(value); \
} while (0)
    ADAPTIVE_META("cell_count", PyLong_FromSize_t(grid->cell_count));
    ADAPTIVE_META("leaf_count", PyLong_FromSize_t(grid->leaf_count));
    ADAPTIVE_META("root_count", PyLong_FromSize_t(grid->root_count));
    ADAPTIVE_META("max_level", PyLong_FromUnsignedLong(grid->max_level));
    ADAPTIVE_META("balanced", PyBool_FromLong(grid->balanced));
#undef ADAPTIVE_META
    return result;

error:
    Py_XDECREF(ids); Py_XDECREF(parents); Py_XDECREF(children); Py_XDECREF(levels);
    Py_XDECREF(leaves); Py_XDECREF(flags); Py_XDECREF(bounds); Py_XDECREF(materials);
    Py_XDECREF(cells); Py_XDECREF(mixed); Py_XDECREF(ties); Py_XDECREF(refinement);
    Py_XDECREF(dominant); Py_XDECREF(errors); Py_XDECREF(samples); Py_XDECREF(result);
    Py_XDECREF(fraction_offsets); Py_XDECREF(fraction_materials);
    Py_XDECREF(fraction_values); Py_XDECREF(cell_fraction_offsets);
    Py_XDECREF(cell_fraction_cells); Py_XDECREF(cell_fraction_materials);
    Py_XDECREF(cell_fraction_values);
    return NULL;
}

static int adaptive_config_from_dict(PyObject *options,
                                     alea_adaptive_grid_config_t *cfg,
                                     mesh_node_buffers_t *nodes) {
    if (mesh_config_from_dict(options, &cfg->sampling, nodes) < 0 ||
        mesh_option_uint32(options, "max_grid_depth", &cfg->max_grid_depth) < 0 ||
        mesh_option_size(options, "max_cells", &cfg->max_cells) < 0 ||
        mesh_option_bool(options, "refine_mixed", &cfg->refine_mixed) < 0 ||
        mesh_option_bool(options, "refine_high_error", &cfg->refine_high_error) < 0)
        return -1;
    return 0;
}

static PyObject *PyAleaSystem_mesh_sample_adaptive(
        PyAleaSystemObject *self, PyObject *options) {
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    alea_adaptive_grid_config_t cfg;
    mesh_node_buffers_t nodes = {0};
    alea_adaptive_grid_config_init(&cfg);
    if (adaptive_config_from_dict(options, &cfg, &nodes) < 0) {
        mesh_node_buffers_free(&nodes);
        return NULL;
    }
    alea_adaptive_grid_result_t *grid;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    grid = alea_adaptive_grid_sample(self->sys, &cfg);
    Py_END_ALLOW_THREADS
    mesh_node_buffers_free(&nodes);
    if (restore_sigint(old_sigint)) {
        alea_adaptive_grid_result_free(grid);
        return NULL;
    }
    if (!grid) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }
    PyObject *result = adaptive_result_to_dict(grid);
    alea_adaptive_grid_result_free(grid);
    return result;
}

static PyObject *PyAleaSystem_mesh_export_adaptive(
        PyAleaSystemObject *self, PyObject *args) {
    const char *filename;
    PyObject *options;
    if (!PyArg_ParseTuple(args, "sO!:mesh_export_adaptive", &filename,
                          &PyDict_Type, &options)) return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    alea_adaptive_grid_config_t cfg;
    mesh_node_buffers_t nodes = {0};
    alea_adaptive_grid_config_init(&cfg);
    if (adaptive_config_from_dict(options, &cfg, &nodes) < 0) {
        mesh_node_buffers_free(&nodes);
        return NULL;
    }
    alea_adaptive_grid_result_t *grid;
    int rc = -1;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    grid = alea_adaptive_grid_sample(self->sys, &cfg);
    if (grid) rc = alea_adaptive_grid_export(grid, cfg.sampling.format, filename);
    Py_END_ALLOW_THREADS
    mesh_node_buffers_free(&nodes);
    if (restore_sigint(old_sigint)) {
        alea_adaptive_grid_result_free(grid);
        return NULL;
    }
    alea_adaptive_grid_result_free(grid);
    if (rc < 0) {
        PyErr_SetString(PyExc_RuntimeError, alea_error());
        return NULL;
    }
    Py_RETURN_NONE;
}
