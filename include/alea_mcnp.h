// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/**
 * @file alea_mcnp.h
 * @brief Public API for MCNP module
 *
 * Users who want to load/export MCNP files should include this header
 * and link against libalea_mcnp.a.
 *
 * Basic usage:
 *   #include "alea.h"
 *   #include "alea_mcnp.h"
 *
 *   mcnp_model_t* model = mcnp_load("input.inp");
 *   alea_system_t* sys = model->sys;
 *   // ... work with sys (geometry, queries, void generation) ...
 *   mcnp_export(model, "output.inp");
 *   mcnp_model_destroy(model);
 */

#ifndef ALEA_MCNP_H
#define ALEA_MCNP_H

#include "alea.h"
#include <stddef.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * MCNP CELL PARAMETERS
 * ============================================================================ */

/**
 * @brief X-macro table for simple keyword=value cell parameters.
 *
 * Each entry: X(field, type, keyword, export_precision)
 * Adding a new parameter = adding one line here. The struct fields,
 * has_* flags, parser, copy, merge, and export are all auto-generated.
 *
 * Complex params (IMP:*, FILL, TRCL, LAT, U, MAT, RHO) stay manual.
 */
#define MCNP_CELL_SIMPLE_PARAMS(X) \
    X(vol,   double, "VOL",   6)   \
    X(tmp,   double, "TMP",   10)  \
    X(pwt,   double, "PWT",   4)   \
    X(nonu,  int,    "NONU",  0)   \
    X(pd,    double, "PD",    4)   \
    X(elpt,  double, "ELPT",  6)   \
    X(unc,   int,    "UNC",   0)   \
    X(bflcl, int,    "BFLCL", 0)

/**
 * @brief MCNP-specific cell parameters (parallel array to sys->cells)
 */
typedef struct {
    double imp_n, imp_p, imp_e;

    #define X_FIELD(name, type, kw, prec) type name;
    MCNP_CELL_SIMPLE_PARAMS(X_FIELD)
    #undef X_FIELD

    int trcl;
    int trcl_inline;
    int trcl_degrees;
    int trcl_count;
    double trcl_data[13];
    int fill_transform_inline;
    int fill_transform_degrees;
    int fill_transform_count;
    double fill_transform_data[13];
    int like_cell_id;

    /* Flags indicating which parameters were explicitly set */
    unsigned int has_imp_n : 1;
    unsigned int has_imp_p : 1;
    unsigned int has_imp_e : 1;

    #define X_HAS(name, type, kw, prec) unsigned int has_##name : 1;
    MCNP_CELL_SIMPLE_PARAMS(X_HAS)
    #undef X_HAS

    unsigned int has_trcl : 1;
    unsigned int has_mat : 1;
    unsigned int has_rho : 1;
} mcnp_cell_params_t;

/**
 * @brief MCNP export configuration
 */
typedef struct {
    int surface_policy;     /* ALEA_EMIT_MACROBODY (0) */
    int trcl_mode;          /* 0=preserve, 1=bake */
    int transform_mode;     /* 0=original, 1=inline, 2=cards */
    int mcnp_max_col;       /* 80 */
    int mcnp_cont_indent;   /* 5 */
} mcnp_export_config_t;

/* ============================================================================
 * MCNP MODEL
 * ============================================================================ */

/**
 * @brief MCNP model: geometry system + MCNP-specific data
 */
typedef struct {
    alea_system_t* sys;
    int owns_sys;                           /* 1 if destroy should free sys */
    mcnp_cell_params_t* cell_params;        /* parallel array, indexed like sys->cells */
    size_t cell_params_count;
    size_t cell_params_capacity;
    mcnp_export_config_t export_config;
} mcnp_model_t;

/**
 * @brief Default MCNP export configuration
 */
extern const mcnp_export_config_t MCNP_EXPORT_CONFIG_DEFAULT;

/* ============================================================================
 * LOADING
 * ============================================================================ */

/**
 * @brief Load an MCNP file into an mcnp_model_t
 * @param filename Path to MCNP input file
 * @return Model on success, NULL on error
 */
mcnp_model_t* mcnp_load(const char* filename);

/**
 * @brief Load MCNP from a string buffer
 * @param input MCNP input text
 * @param len Length of input (0 = use strlen)
 * @return Model on success, NULL on error
 */
mcnp_model_t* mcnp_load_string(const char* input, size_t len);

/* ============================================================================
 * EXPORT
 * ============================================================================ */

/**
 * @brief Export MCNP model to file
 * @return 0 on success, -1 on error
 */
int mcnp_export(const mcnp_model_t* model, const char* filename);

/**
 * @brief Export MCNP model to stream
 * @return 0 on success, -1 on error
 */
int mcnp_export_stream(const mcnp_model_t* model, FILE* out);

/**
 * @brief Export a bare system to MCNP file (convenience, uses default config)
 * @return 0 on success, -1 on error
 */
int mcnp_export_system(alea_system_t* sys, const char* filename);

/**
 * @brief Export a bare system to MCNP stream (convenience, uses default config)
 * @return 0 on success, -1 on error
 */
int mcnp_export_system_stream(alea_system_t* sys, FILE* out);

/* ============================================================================
 * MODEL MANAGEMENT
 * ============================================================================ */

/**
 * @brief Destroy an MCNP model (frees sys + params)
 */
void mcnp_model_destroy(mcnp_model_t* model);

/**
 * @brief Get cell params by index (bounds-checked)
 * @return Pointer to params, or NULL if out of range
 */
mcnp_cell_params_t* mcnp_cell_params(mcnp_model_t* m, size_t idx);

/**
 * @brief Get cell params by index (const, bounds-checked)
 */
const mcnp_cell_params_t* mcnp_cell_params_const(const mcnp_model_t* m, size_t idx);

/**
 * @brief Ensure cell_params array has room for at least `cap` entries
 * @return 0 on success, -1 on allocation failure
 */
int mcnp_model_reserve_params(mcnp_model_t* model, size_t cap);

/**
 * @brief Append a new cell params entry with defaults (imp=1.0)
 * @return Index of new entry, or -1 on failure
 */
int mcnp_model_add_params(mcnp_model_t* model);

/**
 * @brief Register cell callbacks on the system for this model
 *
 * After calling this, any alea_add_cell() / split / merge / etc.
 * will automatically grow/copy the parallel params array.
 */
void mcnp_model_register_hooks(mcnp_model_t* model);

/**
 * @brief Create a non-owning model wrapper around an existing system
 *
 * The returned model does NOT own the system — mcnp_model_destroy()
 * will free the model and its cell_params but not the system.
 * Cell params are populated with defaults (imp=1.0).
 *
 * @param sys System to wrap (caller retains ownership)
 * @return Model on success, NULL on failure
 */
mcnp_model_t* mcnp_model_wrap(alea_system_t* sys);

#ifdef __cplusplus
}
#endif

#endif /* ALEA_MCNP_H */
