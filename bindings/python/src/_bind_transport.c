// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/* Included by pyalea_binding.c after _bind_nucdata.c. */

static int tr_number(PyObject* dict, const char* key, double* value) {
    PyObject* item = PyDict_GetItemString(dict, key);
    if (!item) return 0;
    *value = PyFloat_AsDouble(item);
    return PyErr_Occurred() ? -1 : 0;
}

static int tr_uint(PyObject* dict, const char* key, unsigned long long* value) {
    PyObject* item = PyDict_GetItemString(dict, key);
    if (!item) return 0;
    *value = PyLong_AsUnsignedLongLong(item);
    return PyErr_Occurred() ? -1 : 0;
}

static int tr_vec3(PyObject* dict, const char* key, double out[3], int required) {
    PyObject* item = PyDict_GetItemString(dict, key);
    if (!item) {
        if (!required) return 0;
        PyErr_Format(PyExc_ValueError, "missing %s", key);
        return -1;
    }
    PyObject* seq = PySequence_Fast(item, "expected a three-element sequence");
    if (!seq) return -1;
    if (PySequence_Fast_GET_SIZE(seq) != 3) {
        PyErr_Format(PyExc_ValueError, "%s must have three elements", key);
        Py_DECREF(seq);
        return -1;
    }
    for (int i = 0; i < 3; ++i) {
        out[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq, i));
        if (PyErr_Occurred()) { Py_DECREF(seq); return -1; }
    }
    Py_DECREF(seq);
    return 0;
}

static int tr_name(PyObject* dict, const char* key, const char* const* names,
                   int count, int fallback) {
    PyObject* item = PyDict_GetItemString(dict, key);
    if (!item) return fallback;
    const char* name = PyUnicode_AsUTF8(item);
    if (!name) return -1;
    for (int i = 0; i < count; ++i)
        if (strcmp(name, names[i]) == 0) return i;
    PyErr_Format(PyExc_ValueError, "unknown %s: %s", key, name);
    return -1;
}

/* Results are freed after conversion. Copy into NumPy-owned buffers so the
 * arrays remain valid independently of the native transport result. */
static PyObject* tr_array_copy(const void* data, size_t count, int typenum) {
    if (count > (size_t)NPY_MAX_INTP) {
        PyErr_SetString(PyExc_OverflowError, "transport result is too large for NumPy");
        return NULL;
    }
    npy_intp dimensions[1] = {(npy_intp)count};
    PyArrayObject* array = (PyArrayObject*)PyArray_SimpleNew(1, dimensions, typenum);
    if (!array) return NULL;
    if (count > 0) {
        if (!data) {
            Py_DECREF(array);
            PyErr_SetString(PyExc_RuntimeError, "missing transport result data");
            return NULL;
        }
        memcpy(PyArray_DATA(array), data, PyArray_NBYTES(array));
    }
    return (PyObject*)array;
}

static PyObject* tr_floats(const double* data, size_t count) {
    return tr_array_copy(data, count, NPY_FLOAT64);
}

static PyObject* tr_ints(const int* data, size_t count) {
    if (!data) Py_RETURN_NONE;
    return tr_array_copy(data, count, NPY_INT);
}

static const char* const tr_scores[] = {
    "track_length", "collision", "reaction_event", "reaction_rate",
    "heating", "local_deposition"
};
static const char* const tr_domains[] = {"cell", "universe", "mesh"};
static const char* const tr_particles[] = {"neutron", "photon", "all"};

typedef struct {
    PyObject_HEAD
    alea_source_t* source;
} PyAleaSourceObject;

static PyTypeObject PyAleaSourceType;

static int tr_source_spec(PyObject* config, alea_source_spec_t* spec) {
    if (!PyDict_Check(config)) {
        PyErr_SetString(PyExc_TypeError, "source description must be a dict");
        return -1;
    }
    memset(spec, 0, sizeof(*spec));
    spec->weight = 1.0;
    int particle = tr_name(config, "particle", tr_particles, 2, 0);
    if (particle < 0) return -1;
    spec->particle = (alea_nuc_particle_t)particle;
    PyObject* space = PyDict_GetItemString(config, "space");
    PyObject* angle = PyDict_GetItemString(config, "angle");
    if (!space || !angle || !PyDict_Check(space) || !PyDict_Check(angle)) {
        PyErr_SetString(PyExc_ValueError, "space and angle must both be dicts");
        return -1;
    }
    if (PyDict_GetItemString(config, "kind") ||
        PyDict_GetItemString(config, "position") ||
        PyDict_GetItemString(config, "direction") ||
        PyDict_GetItemString(config, "lower") ||
        PyDict_GetItemString(config, "upper")) {
        PyErr_SetString(PyExc_ValueError,
            "put spatial and angular fields inside space and angle");
        return -1;
    }
    static const char* const spaces[] = {"point", "box", "line", "sphere", "cylinder"};
    static const char* const angles[] = {
        "monodirectional", "isotropic", "cone", "cosine", "tabulated_mu", "radial"
    };
    int spatial = tr_name(space, "type", spaces, 5, -1);
    int angular = tr_name(angle, "type", angles, 6, -1);
    if (spatial < 0 || angular < 0) {
        if (!PyErr_Occurred()) PyErr_SetString(PyExc_ValueError, "space and angle require type");
        return -1;
    }
    spec->space = (alea_source_space_t)spatial;
    spec->angle = (alea_source_angle_t)angular;
    if (spatial == 0) {
        if (tr_vec3(space, "position", spec->position, 1) < 0) return -1;
    } else if (spatial == 1) {
        if (tr_vec3(space, "lower", spec->lower, 1) < 0 ||
            tr_vec3(space, "upper", spec->upper, 1) < 0) return -1;
    } else if (spatial == 2) {
        if (tr_vec3(space, "start", spec->start, 1) < 0 ||
            tr_vec3(space, "end", spec->end, 1) < 0) return -1;
    } else {
        if (spatial == 3) {
            if (tr_vec3(space, "center", spec->center, 1) < 0) return -1;
        } else if (tr_vec3(space, "base", spec->base, 1) < 0 ||
                   tr_vec3(space, "axis", spec->axis, 1) < 0) return -1;
        PyObject* radius = PyDict_GetItemString(space, "outer_radius");
        if (!radius) { PyErr_SetString(PyExc_ValueError, "volume source requires outer_radius"); return -1; }
        spec->outer_radius = PyFloat_AsDouble(radius);
        if (PyErr_Occurred() ||
            tr_number(space, "inner_radius", &spec->inner_radius) < 0) return -1;
    }
    if ((angular == 0 || angular == 2 || angular == 3 || angular == 4) &&
        tr_vec3(angle, "direction", spec->direction, 1) < 0)
        return -1;
    if (angular == 2) {
        PyObject* width = PyDict_GetItemString(angle, "half_angle");
        if (!width) { PyErr_SetString(PyExc_ValueError, "cone requires half_angle in radians"); return -1; }
        spec->cone_half_angle = PyFloat_AsDouble(width);
        if (PyErr_Occurred()) return -1;
    } else if (angular == 4) {
        static const char* const interpolation[] = {"histogram", "linear"};
        int mode = tr_name(angle, "interpolation", interpolation, 2, -1);
        if (mode < 0) {
            if (!PyErr_Occurred()) PyErr_SetString(PyExc_ValueError, "tabulated_mu requires interpolation");
            return -1;
        }
        spec->angle_interpolation = (alea_source_pdf_t)mode;
    } else if (angular == 5) {
        if (tr_vec3(angle, "origin", spec->angle_origin, 1) < 0) return -1;
        PyObject* inward = PyDict_GetItemString(angle, "inward");
        if (inward) {
            spec->radial_inward = PyObject_IsTrue(inward);
            if (spec->radial_inward < 0) return -1;
        }
    }
    PyObject* energy = PyDict_GetItemString(config, "energy");
    if (!energy) { PyErr_SetString(PyExc_ValueError, "source requires energy in MeV"); return -1; }
    if (PyDict_Check(energy)) {
        static const char* const types[] = {"mono", "lines"};
        int energy_type = tr_name(energy, "type", types, 2, -1);
        if (energy_type < 0) {
            if (!PyErr_Occurred()) PyErr_SetString(PyExc_ValueError, "energy requires type");
            return -1;
        }
        spec->energy_type = (alea_source_energy_t)energy_type;
        if (energy_type == 0) {
            PyObject* value = PyDict_GetItemString(energy, "value");
            if (!value) { PyErr_SetString(PyExc_ValueError, "mono energy requires value"); return -1; }
            spec->energy = PyFloat_AsDouble(value);
        }
    } else spec->energy = PyFloat_AsDouble(energy);
    if (PyErr_Occurred()) return -1;
    PyObject* time = PyDict_GetItemString(config, "time");
    if (time) {
        if (PyDict_Check(time)) {
            static const char* const types[] = {"constant"};
            if (tr_name(time, "type", types, 1, -1) < 0) {
                if (!PyErr_Occurred()) PyErr_SetString(PyExc_ValueError, "time requires type");
                return -1;
            }
            PyObject* value = PyDict_GetItemString(time, "value");
            if (!value) { PyErr_SetString(PyExc_ValueError, "constant time requires value"); return -1; }
            spec->time = PyFloat_AsDouble(value);
        } else spec->time = PyFloat_AsDouble(time);
        if (PyErr_Occurred()) return -1;
    }
    return tr_number(config, "weight", &spec->weight);
}

static alea_source_t* tr_prepare_source(PyObject* config) {
    alea_source_spec_t spec;
    if (tr_source_spec(config, &spec) < 0) return NULL;
    double *values = NULL, *weights = NULL;
    double *angle_mu = NULL, *angle_pdf = NULL;
    if (spec.angle == ALEA_SOURCE_TABULATED_MU) {
        PyObject* angle = PyDict_GetItemString(config, "angle");
        PyObject* raw_mu = PyDict_GetItemString(angle, "mu");
        PyObject* raw_pdf = PyDict_GetItemString(angle, "pdf");
        if (!raw_mu || !raw_pdf) {
            PyErr_SetString(PyExc_ValueError, "tabulated_mu requires mu and pdf");
            return NULL;
        }
        PyObject* mu_seq = PySequence_Fast(raw_mu, "mu must be a sequence");
        if (!mu_seq) return NULL;
        PyObject* pdf_seq = PySequence_Fast(raw_pdf, "pdf must be a sequence");
        if (!pdf_seq) { Py_DECREF(mu_seq); return NULL; }
        Py_ssize_t count = PySequence_Fast_GET_SIZE(mu_seq);
        if (count < 2 || count != PySequence_Fast_GET_SIZE(pdf_seq) ||
            (size_t)count > UINT32_MAX) {
            PyErr_SetString(PyExc_ValueError, "tabulated_mu needs equal arrays of at least two entries");
            Py_DECREF(mu_seq); Py_DECREF(pdf_seq); return NULL;
        }
        angle_mu = PyMem_New(double, count);
        angle_pdf = PyMem_New(double, count);
        if (!angle_mu || !angle_pdf) {
            PyErr_NoMemory();
            Py_DECREF(mu_seq); Py_DECREF(pdf_seq);
            PyMem_Free(angle_mu); PyMem_Free(angle_pdf); return NULL;
        }
        for (Py_ssize_t i = 0; i < count; ++i) {
            angle_mu[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(mu_seq, i));
            if (!PyErr_Occurred())
                angle_pdf[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(pdf_seq, i));
            if (PyErr_Occurred()) {
                Py_DECREF(mu_seq); Py_DECREF(pdf_seq);
                PyMem_Free(angle_mu); PyMem_Free(angle_pdf); return NULL;
            }
        }
        Py_DECREF(mu_seq); Py_DECREF(pdf_seq);
        spec.angle_mu = angle_mu;
        spec.angle_pdf = angle_pdf;
        spec.angle_count = (size_t)count;
    }
    if (spec.energy_type == ALEA_SOURCE_ENERGY_LINES) {
        PyObject* energy = PyDict_GetItemString(config, "energy");
        PyObject* raw_values = PyDict_GetItemString(energy, "values");
        PyObject* raw_weights = PyDict_GetItemString(energy, "weights");
        if (!raw_values || !raw_weights) {
            PyErr_SetString(PyExc_ValueError, "energy lines require values and weights");
            PyMem_Free(angle_mu); PyMem_Free(angle_pdf);
            return NULL;
        }
        PyObject* value_seq = PySequence_Fast(raw_values, "energy values must be a sequence");
        if (!value_seq) { PyMem_Free(angle_mu); PyMem_Free(angle_pdf); return NULL; }
        PyObject* weight_seq = PySequence_Fast(raw_weights, "energy weights must be a sequence");
        if (!weight_seq) { Py_DECREF(value_seq); PyMem_Free(angle_mu); PyMem_Free(angle_pdf); return NULL; }
        Py_ssize_t count = PySequence_Fast_GET_SIZE(value_seq);
        if (count == 0 || count != PySequence_Fast_GET_SIZE(weight_seq) ||
            (size_t)count > UINT32_MAX) {
            PyErr_SetString(PyExc_ValueError, "energy lines need equal nonempty arrays");
            Py_DECREF(value_seq); Py_DECREF(weight_seq);
            PyMem_Free(angle_mu); PyMem_Free(angle_pdf); return NULL;
        }
        values = PyMem_New(double, count);
        weights = PyMem_New(double, count);
        if (!values || !weights) {
            PyErr_NoMemory();
            Py_DECREF(value_seq); Py_DECREF(weight_seq);
            PyMem_Free(values); PyMem_Free(weights);
            PyMem_Free(angle_mu); PyMem_Free(angle_pdf); return NULL;
        }
        for (Py_ssize_t i = 0; i < count; ++i) {
            values[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(value_seq, i));
            if (!PyErr_Occurred())
                weights[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(weight_seq, i));
            if (PyErr_Occurred()) {
                Py_DECREF(value_seq); Py_DECREF(weight_seq);
                PyMem_Free(values); PyMem_Free(weights);
                PyMem_Free(angle_mu); PyMem_Free(angle_pdf); return NULL;
            }
        }
        Py_DECREF(value_seq); Py_DECREF(weight_seq);
        spec.energy_values = values;
        spec.energy_weights = weights;
        spec.energy_count = (size_t)count;
    }
    alea_source_t* source = NULL;
    alea_error_t err = alea_source_prepare(&spec, &source);
    PyMem_Free(values);
    PyMem_Free(weights);
    PyMem_Free(angle_mu);
    PyMem_Free(angle_pdf);
    if (err != ALEA_OK) PyErr_Format(PyExc_ValueError,
        "invalid source: %s", alea_error_string(err));
    return source;
}

static PyObject* tr_source_new(PyTypeObject* type, PyObject* args, PyObject* kwds) {
    (void)args; (void)kwds;
    PyAleaSourceObject* self = (PyAleaSourceObject*)type->tp_alloc(type, 0);
    if (self) self->source = NULL;
    return (PyObject*)self;
}

static int tr_source_init(PyAleaSourceObject* self, PyObject* args, PyObject* kwds) {
    static char* names[] = {"config", NULL};
    PyObject* config;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!", names, &PyDict_Type, &config)) return -1;
    if (self->source) { PyErr_SetString(PyExc_RuntimeError, "Source is already prepared"); return -1; }
    self->source = tr_prepare_source(config);
    return self->source ? 0 : -1;
}

static void tr_source_dealloc(PyAleaSourceObject* self) {
    alea_source_free(self->source);
    Py_TYPE(self)->tp_free((PyObject*)self);
}

static PyTypeObject PyAleaSourceType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "pyalea._alea.Source",
    .tp_doc = PyDoc_STR("Prepared primary source; does not load nuclear data."),
    .tp_basicsize = sizeof(PyAleaSourceObject),
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_new = tr_source_new,
    .tp_init = (initproc)tr_source_init,
    .tp_dealloc = (destructor)tr_source_dealloc,
};

static PyObject* mod_sample_source(PyObject* module, PyObject* args, PyObject* kwds) {
    (void)module;
    static char* names[] = {"source", "histories", "seed", "history_offset", NULL};
    PyObject* input;
    unsigned int count = 1, offset = 0;
    unsigned long long seed = 1;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|IKI", names,
        &input, &count, &seed, &offset)) return NULL;
    if (count && count - 1 > UINT32_MAX - offset) {
        PyErr_SetString(PyExc_OverflowError, "history range exceeds uint32"); return NULL;
    }
    alea_source_t* owned = NULL;
    const alea_source_t* source;
    if (PyObject_TypeCheck(input, &PyAleaSourceType))
        source = ((PyAleaSourceObject*)input)->source;
    else if (PyDict_Check(input)) source = owned = tr_prepare_source(input);
    else { PyErr_SetString(PyExc_TypeError, "source must be Source or dict"); return NULL; }
    if (!source) return NULL;
    npy_intp vector_dims[2] = {(npy_intp)count, 3};
    npy_intp scalar_dims[1] = {(npy_intp)count};
    PyArrayObject* position = (PyArrayObject*)PyArray_SimpleNew(2, vector_dims, NPY_FLOAT64);
    PyArrayObject* direction = (PyArrayObject*)PyArray_SimpleNew(2, vector_dims, NPY_FLOAT64);
    PyArrayObject* energy = (PyArrayObject*)PyArray_SimpleNew(1, scalar_dims, NPY_FLOAT64);
    PyArrayObject* time = (PyArrayObject*)PyArray_SimpleNew(1, scalar_dims, NPY_FLOAT64);
    PyArrayObject* weight = (PyArrayObject*)PyArray_SimpleNew(1, scalar_dims, NPY_FLOAT64);
    PyArrayObject* particle = (PyArrayObject*)PyArray_SimpleNew(1, scalar_dims, NPY_INT32);
    PyArrayObject* history_id = (PyArrayObject*)PyArray_SimpleNew(1, scalar_dims, NPY_UINT32);
    if (!position || !direction || !energy || !time || !weight ||
        !particle || !history_id) goto preview_fail;
    alea_error_t err = ALEA_OK;
    sighandler_func old = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    for (uint32_t i = 0; i < count; ++i) {
        alea_transport_source_t sample;
        err = alea_source_sample((void*)source, (uint64_t)seed, offset + i, &sample);
        if (err != ALEA_OK) break;
        double* p = (double*)PyArray_DATA(position) + 3*(size_t)i;
        double* d = (double*)PyArray_DATA(direction) + 3*(size_t)i;
        memcpy(p, sample.position, 3*sizeof(double));
        memcpy(d, sample.particle.direction, 3*sizeof(double));
        ((double*)PyArray_DATA(energy))[i] = sample.particle.energy;
        ((double*)PyArray_DATA(time))[i] = sample.particle.time;
        ((double*)PyArray_DATA(weight))[i] = sample.particle.weight;
        ((int32_t*)PyArray_DATA(particle))[i] = (int32_t)sample.particle.type;
        ((uint32_t*)PyArray_DATA(history_id))[i] = offset + i;
    }
    Py_END_ALLOW_THREADS
    int interrupted = restore_sigint(old);
    if (interrupted) goto preview_fail;
    if (err != ALEA_OK) {
        PyErr_Format(PyExc_RuntimeError, "source sampling failed: %s", alea_error_string(err));
        goto preview_fail;
    }
    PyObject* out = PyDict_New();
    if (!out) goto preview_fail;
    if (dict_set_new(out, "position", (PyObject*)position) < 0) { position = NULL; goto preview_dict_fail; }
    position = NULL;
    if (dict_set_new(out, "direction", (PyObject*)direction) < 0) { direction = NULL; goto preview_dict_fail; }
    direction = NULL;
    if (dict_set_new(out, "energy", (PyObject*)energy) < 0) { energy = NULL; goto preview_dict_fail; }
    energy = NULL;
    if (dict_set_new(out, "time", (PyObject*)time) < 0) { time = NULL; goto preview_dict_fail; }
    time = NULL;
    if (dict_set_new(out, "weight", (PyObject*)weight) < 0) { weight = NULL; goto preview_dict_fail; }
    weight = NULL;
    if (dict_set_new(out, "particle", (PyObject*)particle) < 0) { particle = NULL; goto preview_dict_fail; }
    particle = NULL;
    if (dict_set_new(out, "history_id", (PyObject*)history_id) < 0) { history_id = NULL; goto preview_dict_fail; }
    history_id = NULL;
    alea_source_free(owned);
    return out;
preview_dict_fail:
    Py_DECREF(out);
preview_fail:
    Py_XDECREF(position); Py_XDECREF(direction); Py_XDECREF(energy);
    Py_XDECREF(time); Py_XDECREF(weight); Py_XDECREF(particle);
    Py_XDECREF(history_id);
    alea_source_free(owned);
    return NULL;
}

static int tr_tally(alea_tally_plan_t* plan, PyObject* item) {
    if (!PyDict_Check(item)) {
        PyErr_SetString(PyExc_TypeError, "each tally must be a dict");
        return -1;
    }
    alea_tally_spec_t spec = {0};
    int score = tr_name(item, "score", tr_scores, 6, 0);
    int domain = tr_name(item, "domain", tr_domains, 3, 0);
    int particle = tr_name(item, "particle", tr_particles, 3, 2);
    if (score < 0 || domain < 0 || particle < 0) return -1;
    spec.score = (alea_tally_score_t)score;
    spec.domain = (alea_tally_domain_t)domain;
    spec.particle_mask = particle == 2 ? 0 : (1u << particle);
    if (tr_vec3(item, "lower", spec.lower, domain == 2) < 0 ||
        tr_vec3(item, "upper", spec.upper, domain == 2) < 0) return -1;
    if (domain == 2) {
        PyObject* dims = PyDict_GetItemString(item, "dimensions");
        if (!dims) { PyErr_SetString(PyExc_ValueError, "mesh requires dimensions"); return -1; }
        PyObject* seq = PySequence_Fast(dims, "dimensions must be a sequence");
        if (!seq) return -1;
        if (PySequence_Fast_GET_SIZE(seq) != 3) {
            Py_DECREF(seq);
            PyErr_SetString(PyExc_ValueError, "dimensions must have three elements");
            return -1;
        }
        for (int i = 0; i < 3; ++i) {
            unsigned long n = PyLong_AsUnsignedLong(PySequence_Fast_GET_ITEM(seq, i));
            if (PyErr_Occurred()) { Py_DECREF(seq); return -1; }
            if (n > UINT32_MAX) { Py_DECREF(seq); PyErr_SetString(PyExc_OverflowError, "mesh dimension exceeds uint32"); return -1; }
            spec.dimensions[i] = (uint32_t)n;
        }
        Py_DECREF(seq);
    }
    if (tr_number(item, "energy_min", &spec.energy_min) < 0 ||
        tr_number(item, "energy_max", &spec.energy_max) < 0 ||
        tr_number(item, "time_min", &spec.time_min) < 0 ||
        tr_number(item, "time_max", &spec.time_max) < 0) return -1;
    unsigned long long n = 0;
    PyObject* material = PyDict_GetItemString(item, "material_id");
    if (material) {
        long value = PyLong_AsLong(material);
        if (PyErr_Occurred()) return -1;
        if (value < -1 || value > INT_MAX) { PyErr_SetString(PyExc_ValueError, "invalid material_id"); return -1; }
        spec.material_id = (int)value;
    }
    if (tr_uint(item, "reaction_mt", &n) < 0 || n > INT_MAX) goto bad_int;
    spec.reaction_mt = (int)n; n = 0;
    if (tr_uint(item, "nuclide_zaid", &n) < 0 || n > INT_MAX) goto bad_int;
    spec.nuclide_zaid = (int)n;
    PyObject* edges = PyDict_GetItemString(item, "energy_edges");
    PyObject* seq = NULL;
    double* edge_data = NULL;
    if (edges) {
        seq = PySequence_Fast(edges, "energy_edges must be a sequence");
        if (!seq) return -1;
        Py_ssize_t len = PySequence_Fast_GET_SIZE(seq);
        if (len < 2) { PyErr_SetString(PyExc_ValueError, "energy_edges needs at least two entries"); goto done; }
        edge_data = PyMem_New(double, len);
        if (!edge_data) { PyErr_NoMemory(); goto done; }
        for (Py_ssize_t i = 0; i < len; ++i) {
            edge_data[i] = PyFloat_AsDouble(PySequence_Fast_GET_ITEM(seq, i));
            if (PyErr_Occurred()) goto done;
        }
        spec.energy_edges = edge_data;
        spec.energy_group_count = (size_t)len - 1;
    }
    alea_error_t err = alea_tally_plan_add(plan, &spec, NULL);
    if (err != ALEA_OK) PyErr_Format(PyExc_ValueError, "invalid tally: %s", alea_error_string(err));
done:
    PyMem_Free(edge_data);
    Py_XDECREF(seq);
    return PyErr_Occurred() ? -1 : 0;
bad_int:
    if (!PyErr_Occurred()) PyErr_SetString(PyExc_OverflowError, "tally integer exceeds int range");
    return -1;
}

static PyObject* tr_result(const alea_transport_result_t* result) {
    PyObject* out = Py_BuildValue("{s:I,s:K,s:K,s:K,s:K,s:K,s:K,s:K,s:K,s:K,s:K,s:K,s:K}",
        "histories", result->histories,
        "absorbed", result->absorbed, "replaced", result->replaced,
        "leaked", result->leaked, "collisions", result->collisions,
        "boundary_crossings", result->boundary_crossings,
        "reflections", result->reflections, "emitted_neutrons", result->emitted_neutrons,
        "emitted_photons", result->emitted_photons,
        "photon_collisions", result->photon_collisions,
        "photon_absorbed", result->photon_absorbed,
        "photon_replaced", result->photon_replaced,
        "photon_leaked", result->photon_leaked);
    if (!out) return NULL;
    PyObject* tallies = PyList_New((Py_ssize_t)alea_tally_results_count(result->tallies));
    if (!tallies) { Py_DECREF(out); return NULL; }
    for (Py_ssize_t i = 0; i < PyList_GET_SIZE(tallies); ++i) {
        alea_tally_view_t v;
        if (alea_tally_results_view(result->tallies, (size_t)i, &v) != ALEA_OK) goto fail;
        PyObject* t = Py_BuildValue("{s:s,s:s,s:s,s:I,s:n,s:n,s:n,s:i,s:i,s:i,s:d,s:d,s:d,s:d}",
            "score", tr_scores[v.score], "domain", tr_domains[v.domain],
            "particle", v.particle_mask == ALEA_TALLY_NEUTRON ? "neutron" :
                        v.particle_mask == ALEA_TALLY_PHOTON ? "photon" : "all",
            "histories", v.histories, "bin_count", (Py_ssize_t)v.bin_count,
            "spatial_bin_count", (Py_ssize_t)v.spatial_bin_count,
            "energy_group_count", (Py_ssize_t)v.energy_group_count,
            "material_id", v.material_id, "reaction_mt", v.reaction_mt,
            "nuclide_zaid", v.nuclide_zaid,
            "energy_min", v.energy_min, "energy_max", v.energy_max,
            "time_min", v.time_min, "time_max", v.time_max);
        if (!t) goto fail;
        if (dict_set_new(t, "sum", tr_floats(v.sum, v.bin_count)) < 0 ||
            dict_set_new(t, "sum_squared", tr_floats(v.sum_squared, v.bin_count)) < 0 ||
            dict_set_new(t, "bin_ids", tr_ints(v.bin_ids, v.spatial_bin_count)) < 0 ||
            dict_set_new(t, "energy_edges", v.energy_edges ?
                tr_floats(v.energy_edges, v.energy_group_count + 1) : Py_NewRef(Py_None)) < 0) {
            Py_DECREF(t); goto fail;
        }
        if (v.bin_count > (size_t)NPY_MAX_INTP) {
            PyErr_SetString(PyExc_OverflowError, "tally is too large for NumPy");
            Py_DECREF(t); goto fail;
        }
        npy_intp dimensions[1] = {(npy_intp)v.bin_count};
        PyArrayObject* means = (PyArrayObject*)PyArray_SimpleNew(1, dimensions, NPY_FLOAT64);
        if (!means) { Py_DECREF(t); goto fail; }
        PyArrayObject* errors = (PyArrayObject*)PyArray_SimpleNew(1, dimensions, NPY_FLOAT64);
        if (!errors) { Py_DECREF(means); Py_DECREF(t); goto fail; }
        double* mean_data = (double*)PyArray_DATA(means);
        double* error_data = (double*)PyArray_DATA(errors);
        for (size_t k = 0; k < v.bin_count; ++k) {
            double h = (double)v.histories;
            mean_data[k] = v.sum[k] / h;
            double variance = h > 1 ? (v.sum_squared[k] - v.sum[k] * v.sum[k] / h) / (h - 1) : 0;
            error_data[k] = sqrt(fmax(0.0, variance / h));
        }
        if (dict_set_new(t, "mean", (PyObject*)means) < 0) {
            Py_DECREF(errors); Py_DECREF(t); goto fail;
        }
        if (dict_set_new(t, "standard_error", (PyObject*)errors) < 0) {
            Py_DECREF(t); goto fail;
        }
        if (dict_set_new(t, "lower", Py_BuildValue("(ddd)",
                v.lower[0], v.lower[1], v.lower[2])) < 0 ||
            dict_set_new(t, "upper", Py_BuildValue("(ddd)",
                v.upper[0], v.upper[1], v.upper[2])) < 0 ||
            dict_set_new(t, "dimensions", Py_BuildValue("(III)",
                v.dimensions[0], v.dimensions[1], v.dimensions[2])) < 0) {
            Py_DECREF(t); goto fail;
        }
        PyList_SET_ITEM(tallies, i, t);
    }
    if (dict_set_new(out, "tallies", tallies) < 0 ||
        dict_set_new(out, "track_length", tr_floats(result->track_length, result->cell_count)) < 0 ||
        dict_set_new(out, "track_length_squared", tr_floats(result->track_length_squared, result->cell_count)) < 0) {
        Py_DECREF(out); return NULL;
    }
    return out;
fail:
    Py_DECREF(tallies); Py_DECREF(out);
    return NULL;
}

static PyObject* mod_transport_run(PyObject* module, PyObject* args, PyObject* kwds) {
    (void)module;
    static char* names[] = {"system", "xsdir", "config", NULL};
    PyAleaSystemObject* system;
    PyAleaXsDirObject* xsdir;
    PyObject* config;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O!O!O!", names,
        &PyAleaSystemType, &system, &PyAleaXsDirType, &xsdir,
        &PyDict_Type, &config)) return NULL;
    if (!system->sys || !xsdir->xsdir) {
        PyErr_SetString(PyExc_RuntimeError, "system and xsdir must be initialized"); return NULL;
    }
    alea_transport_options_t options = {.histories = 1, .seed = 1,
        .max_events_per_history = 100000, .max_segment_distance = 100.0};
    unsigned long long n = options.histories;
    if (tr_uint(config, "histories", &n) < 0) return NULL;
    if (n == 0 || n > UINT32_MAX) { PyErr_SetString(PyExc_ValueError, "histories must be 1..UINT32_MAX"); return NULL; }
    options.histories = (uint32_t)n;
    n = options.seed;
    if (tr_uint(config, "seed", &n) < 0) return NULL;
    options.seed = (uint64_t)n;
    n = 0; if (tr_uint(config, "history_offset", &n) < 0) return NULL;
    if (n > UINT32_MAX) { PyErr_SetString(PyExc_ValueError, "history_offset exceeds UINT32_MAX"); return NULL; }
    options.history_offset = (uint32_t)n;
    n = options.max_events_per_history;
    if (tr_uint(config, "max_events_per_history", &n) < 0) return NULL;
    if (n > UINT32_MAX) { PyErr_SetString(PyExc_ValueError, "max_events_per_history exceeds UINT32_MAX"); return NULL; }
    options.max_events_per_history = (uint32_t)n;
    n = 0; if (tr_uint(config, "max_pending_particles", &n) < 0) return NULL;
    if (n > SIZE_MAX) { PyErr_SetString(PyExc_ValueError, "max_pending_particles exceeds SIZE_MAX"); return NULL; }
    options.max_pending_particles = (size_t)n;
    if (tr_number(config, "max_segment_distance", &options.max_segment_distance) < 0) return NULL;
    PyObject* coupled_obj = PyDict_GetItemString(config, "coupled");
    int coupled = coupled_obj ? PyObject_IsTrue(coupled_obj) : 0;
    if (coupled < 0) return NULL;
    PyObject* source_config = PyDict_GetItemString(config, "source");
    Py_XINCREF(source_config); /* pin a prepared Source while the GIL is released */
    alea_source_t* owned_source = NULL;
    const alea_source_t* source = NULL;
    if (source_config && PyObject_TypeCheck(source_config, &PyAleaSourceType))
        source = ((PyAleaSourceObject*)source_config)->source;
    else if (source_config && PyDict_Check(source_config))
        source = owned_source = tr_prepare_source(source_config);
    else PyErr_SetString(PyExc_ValueError, "source must be Source or dict");
    if (!source) { Py_XDECREF(source_config); return NULL; }
    alea_tally_plan_t* plan = NULL;
    PyObject* tallies = PyDict_GetItemString(config, "tallies");
    if (tallies) {
        PyObject* seq = PySequence_Fast(tallies, "tallies must be a sequence");
        if (!seq) { alea_source_free(owned_source); Py_DECREF(source_config); return NULL; }
        plan = alea_tally_plan_create(system->sys);
        if (!plan) { Py_DECREF(seq); alea_source_free(owned_source);
            Py_DECREF(source_config); return PyErr_NoMemory(); }
        for (Py_ssize_t i = 0; i < PySequence_Fast_GET_SIZE(seq); ++i) {
            if (tr_tally(plan, PySequence_Fast_GET_ITEM(seq, i)) < 0) {
                Py_DECREF(seq); alea_tally_plan_free(plan);
                alea_source_free(owned_source); Py_DECREF(source_config); return NULL;
            }
        }
        Py_DECREF(seq);
    }
    options.tally_plan = plan;
    alea_nuc_prepare_requirements_t req = {.required_capabilities =
        ALEA_NUC_CAP_CONTINUOUS_NEUTRON |
        (coupled ? ALEA_NUC_CAP_PHOTON_PRODUCTION : 0)};
    uint32_t mask = alea_source_particle_mask(source) |
        (coupled ? ALEA_NUC_BIND_NEUTRON | ALEA_NUC_BIND_PHOTON : 0);
    alea_nuc_cell_bindings_t* bindings = NULL;
    alea_transport_result_t result = {0};
    alea_transport_failure_t failure = {0};
    alea_error_t err;
    int prepared = 0;
    sighandler_func old = install_sigint();
    Py_BEGIN_ALLOW_THREADS
    err = alea_nuc_cell_bindings_prepare(system->sys, xsdir->xsdir, mask,
        mask & ALEA_NUC_BIND_NEUTRON ? &req : NULL, NULL, &bindings);
    if (err == ALEA_OK) {
        prepared = 1;
        err = alea_transport_run_sampled_source(system->sys, bindings,
            alea_source_sample, (void*)source, &options, &result, &failure);
    }
    Py_END_ALLOW_THREADS
    int interrupted = restore_sigint(old);
    PyObject* out = NULL;
    if (!interrupted) {
        if (err == ALEA_OK) out = tr_result(&result);
        else if (!prepared) PyErr_Format(PyExc_RuntimeError,
            "nuclear-data binding failed: %s", alea_error_string(err));
        else PyErr_Format(PyExc_RuntimeError,
            "transport failed: %s (history %u, cell %d, position %.6g %.6g %.6g)",
            alea_error_string(err), failure.history_id, failure.cell_id,
            failure.position[0], failure.position[1], failure.position[2]);
    }
    alea_transport_result_free(&result);
    alea_nuc_cell_bindings_free(bindings);
    alea_tally_plan_free(plan);
    alea_source_free(owned_source);
    Py_DECREF(source_config);
    return out;
}
