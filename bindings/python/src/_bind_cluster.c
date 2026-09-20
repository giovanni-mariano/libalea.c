// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

/* Included by pyalea_binding.c; all MPI calls use the initializing thread. */
#ifdef PYALEA_USE_MPI
#include "alea_cluster.h"
#include "alea_cluster_mpi.h"

typedef struct {
    PyObject_HEAD
    alea_cluster_t* cluster;
    int busy;
} PyAleaClusterObject;

static PyTypeObject PyAleaClusterType;
static unsigned long cluster_runtime_thread;
static int cluster_runtime_active;

static int cluster_thread_ok(void) {
    if (!cluster_runtime_active || cluster_runtime_thread != PyThread_get_thread_ident()) {
        PyErr_SetString(PyExc_RuntimeError,
            "cluster calls must run on the thread that initialized the runtime");
        return 0;
    }
    return 1;
}

static int cluster_ready(PyAleaClusterObject* self) {
    if (!cluster_thread_ok()) return 0;
    if (!self->cluster || self->busy) {
        PyErr_SetString(PyExc_RuntimeError, self->busy
            ? "cluster operation already active" : "cluster session is closed");
        return 0;
    }
    return 1;
}

static PyObject* mod_cluster_initialize(PyObject* ignored, PyObject* unused) {
    (void)ignored; (void)unused;
    if (cluster_runtime_active) {
        if (!cluster_thread_ok()) return NULL;
        Py_RETURN_NONE;
    }
    int mpi_initialized = 0, is_main_thread = 1;
    if (MPI_Initialized(&mpi_initialized) != MPI_SUCCESS ||
        (mpi_initialized &&
         (MPI_Is_thread_main(&is_main_thread) != MPI_SUCCESS ||
          !is_main_thread))) {
        PyErr_SetString(PyExc_RuntimeError,
            "MPI was initialized on another thread; cluster requires its MPI main thread");
        return NULL;
    }
    alea_cluster_status_t status = alea_cluster_initialize(NULL, NULL);
    if (status != ALEA_CLUSTER_OK) {
        PyErr_Format(PyExc_RuntimeError, "cluster initialization failed: %s",
                     alea_cluster_status_string(status));
        return NULL;
    }
    cluster_runtime_thread = PyThread_get_thread_ident();
    cluster_runtime_active = 1;
    Py_RETURN_NONE;
}

static PyObject* mod_cluster_finalize(PyObject* ignored, PyObject* unused) {
    (void)ignored; (void)unused;
    if (!cluster_thread_ok()) return NULL;
    alea_cluster_status_t status = alea_cluster_finalize();
    if (status != ALEA_CLUSTER_OK) {
        PyErr_Format(PyExc_RuntimeError, "cluster finalization failed: %s",
                     alea_cluster_status_string(status));
        return NULL;
    }
    cluster_runtime_active = 0;
    Py_RETURN_NONE;
}

static void PyAleaCluster_dealloc(PyAleaClusterObject* self) {
    /* Collective MPI_Comm_free is unsafe from GC. Explicit close is required. */
    if (self->cluster) {
        PyErr_WarnEx(PyExc_ResourceWarning,
            "unclosed cluster session; call close() collectively", 1);
        PyErr_Clear();
    }
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyObject* mod_cluster_create(PyObject* ignored, PyObject* arg) {
    (void)ignored;
    if (!cluster_thread_ok()) return NULL;
    MPI_Comm communicator = MPI_COMM_WORLD;
    if (arg != Py_None) {
        long handle = PyLong_AsLong(arg);
        if (PyErr_Occurred()) return NULL;
        communicator = MPI_Comm_f2c((MPI_Fint)handle);
    }
    alea_cluster_t* cluster = alea_cluster_create_mpi(communicator);
    if (!cluster) {
        PyErr_SetString(PyExc_RuntimeError, "collective cluster creation failed");
        return NULL;
    }
    PyAleaClusterObject* self = (PyAleaClusterObject*)
        PyAleaClusterType.tp_alloc(&PyAleaClusterType, 0);
    alea_cluster_status_t status = alea_cluster_agree(cluster, self
        ? ALEA_CLUSTER_OK : ALEA_CLUSTER_OUT_OF_MEMORY);
    if (status != ALEA_CLUSTER_OK) {
        alea_cluster_destroy(cluster);
        Py_XDECREF(self);
        if (!PyErr_Occurred())
            PyErr_SetString(PyExc_MemoryError,
                "a cluster rank could not allocate its Python context");
        return NULL;
    }
    self->cluster = cluster;
    return (PyObject*)self;
}

static PyObject* PyAleaCluster_close(PyAleaClusterObject* self, PyObject* unused) {
    (void)unused;
    if (!cluster_ready(self)) return NULL;
    alea_cluster_destroy(self->cluster);
    self->cluster = NULL;
    Py_RETURN_NONE;
}

static PyObject* PyAleaCluster_get_rank(PyAleaClusterObject* self, void* unused) {
    (void)unused;
    if (!cluster_ready(self)) return NULL;
    return PyLong_FromLong(alea_cluster_rank(self->cluster));
}
static PyObject* PyAleaCluster_get_size(PyAleaClusterObject* self, void* unused) {
    (void)unused;
    if (!cluster_ready(self)) return NULL;
    return PyLong_FromLong(alea_cluster_size(self->cluster));
}

static PyObject* PyAleaCluster_agree(PyAleaClusterObject* self, PyObject* arg) {
    if (!cluster_ready(self)) return NULL;
    int local = PyObject_IsTrue(arg);
    if (local < 0) return NULL;
    self->busy = 1;
    alea_cluster_status_t status;
    status = alea_cluster_agree(self->cluster, local
        ? ALEA_CLUSTER_INVALID_ARGUMENT : ALEA_CLUSTER_OK);
    self->busy = 0;
    return PyLong_FromLong(status);
}

static PyObject* PyAleaCluster_read(PyAleaClusterObject* self, PyObject* args) {
    int mcnp;
    PyObject* path_obj;
    if (!PyArg_ParseTuple(args, "iO", &mcnp, &path_obj)) return NULL;
    if (!cluster_ready(self)) return NULL;
    const char* path = NULL;
    if (alea_cluster_rank(self->cluster) == 0) {
        if (PyUnicode_Check(path_obj)) path = PyUnicode_AsUTF8(path_obj);
        else if (!PyErr_Occurred())
            PyErr_SetString(PyExc_TypeError, "root path must be a string");
    }
    alea_cluster_status_t prep = PyErr_Occurred()
        ? ALEA_CLUSTER_INVALID_ARGUMENT : ALEA_CLUSTER_OK;
    alea_cluster_status_t status = alea_cluster_agree(self->cluster, prep);
    if (status != ALEA_CLUSTER_OK) {
        if (!PyErr_Occurred())
            PyErr_SetString(PyExc_RuntimeError, "another rank has an invalid input path");
        return NULL;
    }
    char* data = NULL;
    size_t length = 0;
    self->busy = 1;
    status = mcnp ? alea_cluster_read_mcnp_input(self->cluster, path, &data, &length)
                  : alea_cluster_read_file(self->cluster, path, &data, &length);
    self->busy = 0;
    if (status != ALEA_CLUSTER_OK) {
        free(data);
        PyErr_Format(PyExc_RuntimeError, "cluster read failed: %s",
                     alea_cluster_status_string(status));
        return NULL;
    }
    PyObject* out = PyBytes_FromStringAndSize(data, (Py_ssize_t)length);
    free(data);
    status = alea_cluster_agree(self->cluster, out
        ? ALEA_CLUSTER_OK : ALEA_CLUSTER_OUT_OF_MEMORY);
    if (status != ALEA_CLUSTER_OK) {
        Py_XDECREF(out);
        if (!PyErr_Occurred())
            PyErr_SetString(PyExc_MemoryError,
                "a cluster rank could not allocate input bytes");
        return NULL;
    }
    return out;
}

static PyObject* PyAleaCluster_estimate_volumes(PyAleaClusterObject* self,
                                                    PyObject* args) {
    PyObject* obj;
    unsigned long long max_rays, seed, workers, batch_size;
    double target, center_x = 0.0, center_y = 0.0, center_z = 0.0, radius = 0.0;
    if (!PyArg_ParseTuple(args, "OKKKKd|dddd", &obj, &max_rays, &seed,
                          &workers, &batch_size, &target,
                          &center_x, &center_y, &center_z, &radius)) return NULL;
    if (!cluster_ready(self)) return NULL;
    PyAleaSystemObject* sys = PyObject_TypeCheck(obj, &PyAleaSystemType)
        ? (PyAleaSystemObject*)obj : NULL;
    alea_cluster_status_t local = sys && sys->sys && max_rays > 0 &&
        isfinite(center_x) && isfinite(center_y) && isfinite(center_z) &&
        isfinite(radius) && radius >= 0.0
        ? ALEA_CLUSTER_OK : ALEA_CLUSTER_INVALID_ARGUMENT;
    size_t count = local == ALEA_CLUSTER_OK ? alea_volume_path_count(sys->sys) : 0;
    if (count == 0) local = ALEA_CLUSTER_INVALID_ARGUMENT;
    double* volumes = local == ALEA_CLUSTER_OK ? calloc(count, sizeof(double)) : NULL;
    double* errors = local == ALEA_CLUSTER_OK ? calloc(count, sizeof(double)) : NULL;
    alea_volume_path_t* paths = local == ALEA_CLUSTER_OK
        ? calloc(count, sizeof(alea_volume_path_t)) : NULL;
    if (local == ALEA_CLUSTER_OK && (!volumes || !errors || !paths))
        local = ALEA_CLUSTER_OUT_OF_MEMORY;
    alea_cluster_status_t status = alea_cluster_agree(self->cluster, local);
    if (status != ALEA_CLUSTER_OK) {
        free(volumes); free(errors); free(paths);
        PyErr_Format(PyExc_RuntimeError, "cluster volume preparation failed: %s",
                     alea_cluster_status_string(status));
        return NULL;
    }
    alea_volume_estimate_options_t options;
    alea_volume_estimate_options_init(&options);
    options.max_rays = (size_t)max_rays;
    options.seed = (uint64_t)seed;
    options.requested_workers = (size_t)workers;
    options.batch_size = (size_t)batch_size;
    options.target_rel_error = target;
    options.use_sampling_sphere = radius > 0.0;
    options.sampling_center[0] = center_x;
    options.sampling_center[1] = center_y;
    options.sampling_center[2] = center_z;
    options.sampling_radius = radius;
    alea_cluster_volume_stats_t stats;
    self->busy = 1;
    status = alea_cluster_estimate_volumes(self->cluster, sys->sys,
                                            &options, volumes, errors, &stats);
    self->busy = 0;
    if (status != ALEA_CLUSTER_OK) {
        free(volumes); free(errors); free(paths);
        PyErr_Format(PyExc_RuntimeError, "cluster volumes failed: %s",
                     alea_cluster_status_string(status));
        return NULL;
    }
    PyObject* vol_list = PyList_New(count);
    PyObject* err_list = PyList_New(count);
    PyObject* path_list = PyList_New(count);
    size_t got = alea_volume_paths_get(sys->sys, paths, count);
    if (got > count) got = count;
    for (size_t i = 0; vol_list && err_list && path_list && i < count; i++) {
        PyObject* v = PyFloat_FromDouble(volumes[i]);
        PyObject* e = PyFloat_FromDouble(errors[i]);
        PyObject* p = i < got ? volume_path_to_dict(&paths[i]) : Py_NewRef(Py_None);
        if (!v || !e || !p) {
            Py_XDECREF(v); Py_XDECREF(e); Py_XDECREF(p);
            break;
        }
        PyList_SET_ITEM(vol_list, i, v);
        PyList_SET_ITEM(err_list, i, e);
        PyList_SET_ITEM(path_list, i, p);
    }
    free(volumes); free(errors); free(paths);
    PyObject* out = NULL;
    if (!PyErr_Occurred() && vol_list && err_list && path_list)
        out = Py_BuildValue("{s:N,s:N,s:N,s:K,s:K,s:K,s:K,s:d,s:i,s:i}",
        "volumes", vol_list, "rel_errors", err_list, "paths", path_list,
        "rays_completed", (unsigned long long)stats.volume.rays_completed,
        "actual_workers", (unsigned long long)stats.volume.actual_workers,
        "rank_count", (unsigned long long)stats.rank_count,
        "local_rays_completed", (unsigned long long)stats.local_rays_completed,
        "maximum_relative_error", stats.volume.maximum_relative_error,
        "converged", stats.volume.converged, "cancelled", stats.volume.cancelled);
    else { Py_XDECREF(vol_list); Py_XDECREF(err_list); Py_XDECREF(path_list); }
    if (out && (dict_set_new(out, "seed", PyLong_FromUnsignedLongLong(stats.volume.seed)) < 0 ||
                dict_set_new(out, "requested_workers", PyLong_FromSize_t(stats.volume.requested_workers)) < 0 ||
                dict_set_new(out, "local_workers", PyLong_FromSize_t(stats.local_workers)) < 0 ||
                dict_set_new(out, "batch_size", PyLong_FromSize_t(stats.volume.batch_size)) < 0 ||
                dict_set_new(out, "rng_algorithm", PyUnicode_FromString(
                    alea_rng_algorithm_name(stats.volume.rng_algorithm))) < 0 ||
                dict_set_new(out, "rng_address_version", PyLong_FromUnsignedLong(
                    stats.volume.rng_address_version)) < 0))
        Py_CLEAR(out);
    /* Result construction can fail on only one rank. Agree before returning. */
    status = alea_cluster_agree(self->cluster, out
        ? ALEA_CLUSTER_OK : ALEA_CLUSTER_OUT_OF_MEMORY);
    if (status != ALEA_CLUSTER_OK) {
        Py_XDECREF(out);
        if (!PyErr_Occurred())
            PyErr_SetString(PyExc_RuntimeError, "a cluster rank could not build volume results");
        return NULL;
    }
    return out;
}

static PyObject* PyAleaCluster_validate_geometry(PyAleaClusterObject* self,
                                                     PyObject* args) {
    PyObject *obj, *options_obj;
    if (!PyArg_ParseTuple(args, "OO", &obj, &options_obj)) return NULL;
    if (!cluster_ready(self)) return NULL;
    PyAleaSystemObject* sys = PyObject_TypeCheck(obj, &PyAleaSystemType)
        ? (PyAleaSystemObject*)obj : NULL;
    alea_geom_validator_options_t options;
    int parsed = parse_validator_options(options_obj, &options);
    alea_cluster_status_t local = !sys || !sys->sys || parsed < 0
        ? ALEA_CLUSTER_INVALID_ARGUMENT : ALEA_CLUSTER_OK;
    alea_cluster_status_t status = alea_cluster_agree(self->cluster, local);
    if (status != ALEA_CLUSTER_OK) {
        if (!PyErr_Occurred())
            PyErr_SetString(PyExc_RuntimeError, "a cluster rank could not prepare validation");
        return NULL;
    }
    alea_geom_validator_result_t result;
    alea_geom_validator_result_init(&result);
    self->busy = 1;
    status = alea_cluster_validate_geometry(self->cluster, sys->sys, &options,
        alea_cluster_is_root(self->cluster) ? &result : NULL);
    self->busy = 0;
    if (status != ALEA_CLUSTER_OK) {
        alea_geom_validator_result_free(&result);
        PyErr_Format(PyExc_RuntimeError, "cluster validation failed: %s",
                     alea_cluster_status_string(status));
        return NULL;
    }
    PyObject* out;
    if (alea_cluster_is_root(self->cluster)) out = build_result_dict(&result);
    else { Py_INCREF(Py_None); out = Py_None; }
    alea_geom_validator_result_free(&result);
    status = alea_cluster_agree(self->cluster, out
        ? ALEA_CLUSTER_OK : ALEA_CLUSTER_OUT_OF_MEMORY);
    if (status != ALEA_CLUSTER_OK) {
        Py_XDECREF(out);
        if (!PyErr_Occurred())
            PyErr_SetString(PyExc_RuntimeError, "a cluster rank could not build validation results");
        return NULL;
    }
    return out;
}

static PyObject* PyAleaCluster_render_3d(PyAleaClusterObject* self,
                                           PyObject* args) {
    PyObject *obj, *eye_obj, *target_obj, *up_obj, *colors_obj;
    int width, height;
    if (!PyArg_ParseTuple(args, "OiiOOOO", &obj, &width, &height,
                          &eye_obj, &target_obj, &up_obj, &colors_obj)) return NULL;
    if (!cluster_ready(self)) return NULL;
    PyAleaSystemObject* sys = PyObject_TypeCheck(obj, &PyAleaSystemType)
        ? (PyAleaSystemObject*)obj : NULL;
    render_config_t config;
    render_config_init(&config);
    config.width = width;
    config.height = height;
    alea_cluster_status_t local = sys && sys->sys && width > 0 && height > 0 &&
        (size_t)width <= SIZE_MAX / (size_t)height / (3 * sizeof(float))
        ? ALEA_CLUSTER_OK : ALEA_CLUSTER_INVALID_ARGUMENT;
    if (local == ALEA_CLUSTER_OK &&
        (render_parse_vec3(eye_obj, "eye", config.eye) < 0 ||
         render_parse_vec3(target_obj, "target", config.target) < 0 ||
         render_parse_vec3(up_obj, "up", config.up) < 0 ||
         render_parse_custom_colors(colors_obj, &config) < 0))
        local = ALEA_CLUSTER_INVALID_ARGUMENT;
    config.eye_set = eye_obj != Py_None;
    config.target_set = target_obj != Py_None;
    render_framebuffer_t* frame = NULL;
    if (local == ALEA_CLUSTER_OK && alea_cluster_is_root(self->cluster)) {
        frame = render_framebuffer_create(width, height, 0);
        if (!frame) local = ALEA_CLUSTER_OUT_OF_MEMORY;
    }
    alea_cluster_status_t status = alea_cluster_agree(self->cluster, local);
    if (status != ALEA_CLUSTER_OK) {
        render_framebuffer_free(frame);
        render_config_free(&config);
        if (!PyErr_Occurred())
            PyErr_Format(PyExc_RuntimeError, "cluster render preparation failed: %s",
                         alea_cluster_status_string(status));
        return NULL;
    }
    render_camera_t camera;
    int setup;
    self->busy = 1;
    setup = render_camera_setup(&camera, &config, sys->sys);
    status = alea_cluster_agree(self->cluster, setup == 0
        ? ALEA_CLUSTER_OK : ALEA_CLUSTER_COMPUTE_ERROR);
    if (status == ALEA_CLUSTER_OK) {
        status = alea_cluster_render_scene(self->cluster, sys->sys,
                                           &config, &camera, frame);
    }
    self->busy = 0;
    if (status != ALEA_CLUSTER_OK) {
        render_framebuffer_free(frame);
        render_config_free(&config);
        PyErr_Format(PyExc_RuntimeError, "cluster render failed: %s",
                     alea_cluster_status_string(status));
        return NULL;
    }
    if (!frame) {
        render_config_free(&config);
        status = alea_cluster_agree(self->cluster, ALEA_CLUSTER_OK);
        if (status != ALEA_CLUSTER_OK) {
            PyErr_SetString(PyExc_RuntimeError, "root could not build render result");
            return NULL;
        }
        Py_RETURN_NONE;
    }
    npy_intp dimensions[3] = {height, width, 3};
    PyObject* rgb = PyArray_SimpleNew(3, dimensions, NPY_UINT8);
    if (rgb) render_tonemap(frame, PyArray_DATA((PyArrayObject*)rgb));
    render_framebuffer_free(frame);
    render_config_free(&config);
    PyObject* result = rgb ? PyDict_New() : NULL;
    if (result && PyDict_SetItemString(result, "rgb", rgb) < 0) {
        Py_DECREF(result);
        result = NULL;
    }
    Py_XDECREF(rgb);
    status = alea_cluster_agree(self->cluster, result
        ? ALEA_CLUSTER_OK : ALEA_CLUSTER_OUT_OF_MEMORY);
    if (status != ALEA_CLUSTER_OK) {
        Py_XDECREF(result);
        if (!PyErr_Occurred())
            PyErr_SetString(PyExc_RuntimeError, "a cluster rank could not build render result");
        return NULL;
    }
    return result;
}

static PyMethodDef PyAleaCluster_methods[] = {
    {"close", (PyCFunction)PyAleaCluster_close, METH_NOARGS, NULL},
    {"agree", (PyCFunction)PyAleaCluster_agree, METH_O, NULL},
    {"read", (PyCFunction)PyAleaCluster_read, METH_VARARGS, NULL},
    {"estimate_volumes", (PyCFunction)PyAleaCluster_estimate_volumes, METH_VARARGS, NULL},
    {"validate_geometry", (PyCFunction)PyAleaCluster_validate_geometry, METH_VARARGS, NULL},
    {"render_3d", (PyCFunction)PyAleaCluster_render_3d, METH_VARARGS, NULL},
    {NULL}
};
static PyGetSetDef PyAleaCluster_getset[] = {
    {"rank", (getter)PyAleaCluster_get_rank, NULL, NULL, NULL},
    {"size", (getter)PyAleaCluster_get_size, NULL, NULL, NULL},
    {NULL}
};
static PyTypeObject PyAleaClusterType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "pyalea._alea.ClusterContext",
    .tp_basicsize = sizeof(PyAleaClusterObject),
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_dealloc = (destructor)PyAleaCluster_dealloc,
    .tp_methods = PyAleaCluster_methods,
    .tp_getset = PyAleaCluster_getset,
};
#endif
