// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* This file is included from pyalea_binding.c.  It deliberately wraps the
 * native tile renderer as one C call: invoking scalar Python ray queries for
 * every image pixel would make call-boundary overhead dominate rendering. */

static int render_parse_vec3(PyObject* object, const char* name, double out[3]) {
    if (!object || object == Py_None) return 0;
    PyObject* sequence = PySequence_Fast(object, name);
    if (!sequence) return -1;
    if (PySequence_Fast_GET_SIZE(sequence) != 3) {
        Py_DECREF(sequence);
        PyErr_Format(PyExc_ValueError, "%s must contain exactly three values", name);
        return -1;
    }
    for (Py_ssize_t i = 0; i < 3; ++i) {
        out[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(sequence, i));
        if (PyErr_Occurred()) {
            Py_DECREF(sequence);
            return -1;
        }
    }
    Py_DECREF(sequence);
    return 1;
}

static int render_parse_clips(PyObject* object, render_config_t* config) {
    if (!object || object == Py_None) return 0;
    PyObject* planes = PySequence_Fast(
        object, "clips must be a sequence of (nx, ny, nz, d) planes");
    if (!planes) return -1;
    const Py_ssize_t count = PySequence_Fast_GET_SIZE(planes);
    if (count > RENDER_MAX_CLIPS) {
        Py_DECREF(planes);
        PyErr_Format(PyExc_ValueError, "at most %d clip planes are supported",
                     RENDER_MAX_CLIPS);
        return -1;
    }
    for (Py_ssize_t index = 0; index < count; ++index) {
        PyObject* entry = PySequence_Fast(
            PySequence_Fast_GET_ITEM(planes, index),
            "each clip plane must contain nx, ny, nz, and d");
        if (!entry) {
            Py_DECREF(planes);
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(entry) != 4) {
            Py_DECREF(entry);
            Py_DECREF(planes);
            PyErr_SetString(PyExc_ValueError,
                            "each clip plane must contain four values");
            return -1;
        }
        double values[4];
        for (int component = 0; component < 4; ++component)
            values[component] = PyFloat_AsDouble(
                PySequence_Fast_GET_ITEM(entry, component));
        Py_DECREF(entry);
        if (PyErr_Occurred()) {
            Py_DECREF(planes);
            return -1;
        }
        const double length = sqrt(values[0] * values[0] +
                                   values[1] * values[1] +
                                   values[2] * values[2]);
        if (!isfinite(length) || length <= 1e-15 || !isfinite(values[3])) {
            Py_DECREF(planes);
            PyErr_SetString(PyExc_ValueError,
                            "clip plane values must be finite and the normal non-zero");
            return -1;
        }
        render_clip_plane_t* clip = &config->clips[index];
        clip->normal[0] = values[0] / length;
        clip->normal[1] = values[1] / length;
        clip->normal[2] = values[2] / length;
        clip->d = values[3] / length;
    }
    config->num_clips = (int)count;
    Py_DECREF(planes);
    return 0;
}

static int render_dict_add_array_copy(PyObject* dictionary, const char* key,
                                      const void* values, int ndim,
                                      const npy_intp* dimensions, int typenum) {
    PyArrayObject* array = (PyArrayObject*)PyArray_SimpleNew(ndim, dimensions, typenum);
    if (!array) return -1;
    if (values) memcpy(PyArray_DATA(array), values, PyArray_NBYTES(array));
    int rc = PyDict_SetItemString(dictionary, key, (PyObject*)array);
    Py_DECREF(array);
    return rc;
}

static int render_color_mode(const char* value, render_color_mode_t* out) {
    if (strcmp(value, "material") == 0) *out = RENDER_COLOR_MATERIAL;
    else if (strcmp(value, "cell") == 0) *out = RENDER_COLOR_CELL;
    else if (strcmp(value, "universe") == 0) *out = RENDER_COLOR_UNIVERSE;
    else if (strcmp(value, "density") == 0) *out = RENDER_COLOR_DENSITY;
    else {
        PyErr_SetString(PyExc_ValueError,
                        "color_by must be 'material', 'cell', 'universe', or 'density'");
        return -1;
    }
    return 0;
}

static int render_python_color_id_compare(const void* left, const void* right) {
    const render_color_entry_t* a = left;
    const render_color_entry_t* b = right;
    return (a->id > b->id) - (a->id < b->id);
}

static int render_python_int_compare(const void* left, const void* right) {
    const int a = *(const int*)left;
    const int b = *(const int*)right;
    return (a > b) - (a < b);
}

static int render_parse_filter_mode(
        const char* value, const char* name, render_filter_mode_t* out) {
    if (strcmp(value, "all") == 0) *out = RENDER_FILTER_ALL;
    else if (strcmp(value, "include") == 0) *out = RENDER_FILTER_INCLUDE;
    else if (strcmp(value, "exclude") == 0) *out = RENDER_FILTER_EXCLUDE;
    else {
        PyErr_Format(PyExc_ValueError,
                     "%s must be 'all', 'include', or 'exclude'", name);
        return -1;
    }
    return 0;
}

static int render_parse_id_filter(
        PyObject* object, const char* name, render_id_filter_t* filter) {
    if (!object || object == Py_None) return 0;
    PyObject* values = PySequence_Fast(object, "filter IDs must be a sequence");
    if (!values) return -1;
    const Py_ssize_t count = PySequence_Fast_GET_SIZE(values);
    if (count > 0 && (size_t)count > SIZE_MAX / sizeof(int)) {
        Py_DECREF(values);
        PyErr_NoMemory();
        return -1;
    }
    int* ids = count ? malloc((size_t)count * sizeof(int)) : NULL;
    if (count && !ids) {
        Py_DECREF(values);
        return PyErr_NoMemory(), -1;
    }
    for (Py_ssize_t i = 0; i < count; i++) {
        PyObject* item = PySequence_Fast_GET_ITEM(values, i);
        if (!PyLong_CheckExact(item)) {
            free(ids);
            Py_DECREF(values);
            PyErr_Format(PyExc_TypeError, "%s must contain only integers", name);
            return -1;
        }
        long id = PyLong_AsLong(item);
        if (id < INT_MIN || id > INT_MAX || PyErr_Occurred()) {
            free(ids);
            Py_DECREF(values);
            if (!PyErr_Occurred())
                PyErr_Format(PyExc_OverflowError, "%s ID is outside C int range", name);
            return -1;
        }
        ids[i] = (int)id;
    }
    Py_DECREF(values);
    if (count > 1)
        qsort(ids, (size_t)count, sizeof(int), render_python_int_compare);
    size_t unique = 0;
    for (Py_ssize_t i = 0; i < count; i++)
        if (unique == 0 || ids[i] != ids[unique - 1])
            ids[unique++] = ids[i];
    filter->ids = ids;
    filter->count = unique;
    return 0;
}

static int render_parse_clip_mode(const char* value, render_clip_mode_t* out) {
    if (strcmp(value, "and") == 0) *out = RENDER_CLIP_AND;
    else if (strcmp(value, "or") == 0) *out = RENDER_CLIP_OR;
    else {
        PyErr_SetString(PyExc_ValueError, "clip_mode must be 'and' or 'or'");
        return -1;
    }
    return 0;
}



static int render_parse_custom_colors(PyObject* object, render_config_t* config) {
    if (!object || object == Py_None) return 0;
    PyObject* entries = PySequence_Fast(
        object, "custom_colors must be a sequence of (id, r, g, b) entries");
    if (!entries) return -1;
    Py_ssize_t count = PySequence_Fast_GET_SIZE(entries);
    if (count > INT_MAX) {
        Py_DECREF(entries);
        PyErr_SetString(PyExc_OverflowError, "too many custom color entries");
        return -1;
    }
    if (count > 0) {
        config->custom_colors = calloc((size_t)count, sizeof(render_color_entry_t));
        if (!config->custom_colors) {
            Py_DECREF(entries);
            PyErr_NoMemory();
            return -1;
        }
    }
    config->num_custom_colors = (int)count;
    for (Py_ssize_t index = 0; index < count; ++index) {
        PyObject* entry = PySequence_Fast(
            PySequence_Fast_GET_ITEM(entries, index),
            "each custom color must be an (id, r, g, b) sequence");
        if (!entry) {
            Py_DECREF(entries);
            return -1;
        }
        if (PySequence_Fast_GET_SIZE(entry) != 4) {
            Py_DECREF(entry);
            Py_DECREF(entries);
            PyErr_SetString(PyExc_ValueError,
                            "each custom color must contain id, r, g, and b");
            return -1;
        }
        long identifier = PyLong_AsLong(PySequence_Fast_GET_ITEM(entry, 0));
        double channels[3];
        for (int channel = 0; channel < 3; ++channel)
            channels[channel] = PyFloat_AsDouble(
                PySequence_Fast_GET_ITEM(entry, channel + 1));
        Py_DECREF(entry);
        if (PyErr_Occurred()) {
            Py_DECREF(entries);
            return -1;
        }
        if (identifier < INT_MIN || identifier > INT_MAX ||
            !isfinite(channels[0]) || !isfinite(channels[1]) ||
            !isfinite(channels[2]) || channels[0] < 0.0 || channels[0] > 1.0 ||
            channels[1] < 0.0 || channels[1] > 1.0 ||
            channels[2] < 0.0 || channels[2] > 1.0) {
            Py_DECREF(entries);
            PyErr_SetString(PyExc_ValueError,
                            "custom color ids must fit int and RGB values must be in [0, 1]");
            return -1;
        }
        render_color_entry_t* color = &config->custom_colors[index];
        color->id = (int)identifier;
        color->r = (float)channels[0];
        color->g = (float)channels[1];
        color->b = (float)channels[2];
    }
    Py_DECREF(entries);
    if (config->num_custom_colors > 1)
        qsort(config->custom_colors, (size_t)config->num_custom_colors,
              sizeof(render_color_entry_t), render_python_color_id_compare);
    config->custom_colors_sorted = 1;
    return 0;
}

static int render_mode(const char* value, render_mode_t* out) {
    if (strcmp(value, "solid") == 0) *out = RENDER_MODE_SOLID;
    else if (strcmp(value, "xray") == 0) *out = RENDER_MODE_XRAY;
    else if (strcmp(value, "depth") == 0) *out = RENDER_MODE_DEPTH;
    else if (strcmp(value, "cell_id") == 0) *out = RENDER_MODE_CELLID;
    else if (strcmp(value, "material_id") == 0) *out = RENDER_MODE_MATID;
    else {
        PyErr_SetString(PyExc_ValueError,
                        "mode must be 'solid', 'xray', 'depth', 'cell_id', or 'material_id'");
        return -1;
    }
    return 0;
}

static PyObject* PyAleaSystem_render_3d(PyAleaSystemObject* self,
                                           PyObject* args, PyObject* kwds) {
    int width = 800, height = 600, shadows = 0, edges = 0, aa_samples = 1;
    int auxiliary = 0;
    double fov = RENDER_DEFAULT_FOV, ortho_height = 0.0;
    double xray_density_scale = 0.1;
    const char* color_by = "material";
    const char* mode = "solid";
    const char* clip_mode = "and";
    const char* material_filter_mode = "all";
    const char* cell_filter_mode = "all";
    PyObject *eye = NULL, *target = NULL, *up = NULL, *background = NULL;
    PyObject *custom_colors = NULL, *clips = NULL;
    PyObject *material_ids = NULL, *cell_ids = NULL;
    static char* kwlist[] = {
        "width", "height", "eye", "target", "up", "fov", "ortho_height",
        "color_by", "mode", "background", "shadows", "edges", "aa_samples",
        "auxiliary", "xray_density_scale", "custom_colors",
        "clips", "clip_mode", "material_filter_mode", "material_ids",
        "cell_filter_mode", "cell_ids", NULL
    };
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|iiOOOddssOppipdOOssOsO", kwlist,
                                     &width, &height, &eye, &target, &up, &fov,
                                     &ortho_height, &color_by, &mode, &background,
                                     &shadows, &edges, &aa_samples,
                                     &auxiliary, &xray_density_scale, &custom_colors,
                                     &clips, &clip_mode, &material_filter_mode,
                                     &material_ids, &cell_filter_mode, &cell_ids))
        return NULL;
    if (!self->sys) {
        PyErr_SetString(PyExc_RuntimeError, "System not initialized");
        return NULL;
    }
    if (width < 1 || height < 1 || aa_samples < 1 || aa_samples > RENDER_MAX_AA) {
        PyErr_Format(PyExc_ValueError, "width and height must be positive; aa_samples must be in [1, %d]",
                     RENDER_MAX_AA);
        return NULL;
    }
    if (!isfinite(xray_density_scale) || xray_density_scale < 0.0) {
        PyErr_SetString(PyExc_ValueError,
                        "xray_density_scale must be finite and non-negative");
        return NULL;
    }
    if (ensure_query_acceleration(self) < 0) return NULL;

    render_config_t config;
    render_config_init(&config);
    config.width = width;
    config.height = height;
    config.fov = fov;
    config.ortho_height = ortho_height;
    config.shadows = shadows;
    config.edges = edges;
    config.aa_samples = aa_samples;
    config.aux_output = auxiliary;
    config.xray_density_scale = (float)xray_density_scale;
    if (render_parse_clips(clips, &config) != 0)
        goto failed_config;
    if (render_parse_clip_mode(clip_mode, &config.clip_mode) != 0 ||
        render_parse_filter_mode(material_filter_mode, "material_filter_mode",
                                 &config.material_filter.mode) != 0 ||
        render_parse_filter_mode(cell_filter_mode, "cell_filter_mode",
                                 &config.cell_filter.mode) != 0 ||
        render_parse_id_filter(material_ids, "material_ids",
                               &config.material_filter) != 0 ||
        render_parse_id_filter(cell_ids, "cell_ids",
                               &config.cell_filter) != 0)
        goto failed_config;
    if (render_parse_custom_colors(custom_colors, &config) != 0)
        goto failed_config;
    if (render_color_mode(color_by, &config.color_mode) != 0 ||
        render_mode(mode, &config.render_mode) != 0)
        goto failed_config;

    int present = render_parse_vec3(eye, "eye", config.eye);
    if (present < 0) goto failed_config;
    config.eye_set = present;
    present = render_parse_vec3(target, "target", config.target);
    if (present < 0) goto failed_config;
    config.target_set = present;
    /* Background is stored as float by libalea; parse it through a temporary
     * double triple to avoid writing doubles into the float buffer. */
    if (render_parse_vec3(up, "up", config.up) < 0) goto failed_config;
    if (background && background != Py_None) {
        double rgb[3];
        if (render_parse_vec3(background, "background", rgb) < 0) goto failed_config;
        for (int i = 0; i < 3; ++i) config.background[i] = (float)rgb[i];
    }

    render_camera_t camera;
    render_framebuffer_t* framebuffer = render_framebuffer_create(width, height, auxiliary);
    if (!framebuffer) {
        PyErr_SetString(PyExc_MemoryError, "failed to allocate 3D render framebuffer");
        goto failed_config;
    }
    int rc;
    sighandler_func old_sigint = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    rc = render_camera_setup(&camera, &config, self->sys);
    if (rc == 0) rc = render_scene(self->sys, &config, &camera, framebuffer);
    Py_END_ALLOW_THREADS
    if (restore_sigint(old_sigint)) goto failed;
    if (rc != 0) {
        const char* detail = alea_error();
        PyErr_Format(PyExc_RuntimeError, "3D render failed: %s",
                     detail ? detail : "unknown libalea error");
        goto failed;
    }

    npy_intp image_dims[3] = {height, width, 3};
    npy_intp plane_dims[2] = {height, width};
    size_t pixels = (size_t)width * (size_t)height;
    uint8_t* rgb = malloc(pixels * 3);
    if (!rgb) { PyErr_NoMemory(); goto failed; }
    render_tonemap(framebuffer, rgb);
    PyObject* out = PyDict_New();
    if (!out) { free(rgb); goto failed; }
    if (render_dict_add_array_copy(out, "rgb", rgb, 3, image_dims, NPY_UINT8) != 0) {
        free(rgb); Py_DECREF(out); goto failed;
    }
    free(rgb);
    if (auxiliary &&
        (render_dict_add_array_copy(out, "depth", framebuffer->depth, 2, plane_dims, NPY_FLOAT) != 0 ||
         render_dict_add_array_copy(out, "cell_ids", framebuffer->cell_id, 2, plane_dims, NPY_INT) != 0 ||
         render_dict_add_array_copy(out, "material_ids", framebuffer->material_id, 2, plane_dims, NPY_INT) != 0 ||
         render_dict_add_array_copy(out, "normals", framebuffer->normal, 3, image_dims, NPY_FLOAT) != 0)) {
        Py_DECREF(out);
        goto failed;
    }
    render_framebuffer_free(framebuffer);
    render_config_free(&config);
    return out;

failed:
    render_framebuffer_free(framebuffer);
failed_config:
    render_config_free(&config);
    return NULL;
}
