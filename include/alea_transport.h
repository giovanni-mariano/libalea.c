// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

/** @file alea_transport.h
 * Fixed-source continuous-energy neutron and photon transport through geometry.
 * Nuclear data is requested only by this explicit transport operation.
 */
#ifndef ALEA_TRANSPORT_H
#define ALEA_TRANSPORT_H

#include "alea_nucdata.h"
#include "alea_raycast.h"
#include "alea_tally.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Explicit, optional transport binding. Geometry loading, plotting and
 * validation never load nuclear data. Both inputs must outlive the binding.
 * Preparation is all-or-nothing and loads every requested table up front.
 * A geometry mutation invalidates lookups; recreate the binding afterwards.
 */
typedef struct alea_nuc_cell_bindings alea_nuc_cell_bindings_t;
#define ALEA_NUC_BIND_NEUTRON (1u << 0)
#define ALEA_NUC_BIND_PHOTON  (1u << 1)

typedef struct {
    const char* material_law; /* name stored in the geometry material */
    const char* xsdir_table;  /* explicit ACE name, e.g. "lwtr.20t" */
} alea_nuc_thermal_name_map_t;

typedef struct {
    /* Zero disables nearest-table fallback. Applies only outside a bracket;
     * an explicit library suffix can use only that exact table. Kelvin. */
    double nearest_temperature_tolerance;
    /* Zero disables natural-isotope omission. This bound is the sum of
     * omitted natural atom abundances per element; explicit nuclides are
     * never omitted. Remaining isotopes are renormalized within the element. */
    double max_omitted_natural_abundance;
    /* Optional explicit mapping for aliases such as OpenMC S(a,b) names.
     * Borrowed for the duration of preparation only. */
    const alea_nuc_thermal_name_map_t* thermal_name_map;
    size_t thermal_name_map_count;
} alea_nuc_binding_policy_t;

typedef enum {
    ALEA_NUC_BINDING_TEMPERATURE_INTERPOLATED,
    ALEA_NUC_BINDING_TEMPERATURE_NEAREST,
    ALEA_NUC_BINDING_NATURAL_ISOTOPE_OMITTED
} alea_nuc_binding_notice_kind_t;

typedef struct {
    alea_nuc_binding_notice_kind_t kind;
    alea_nuc_table_type_t table_type;
    size_t cell_index;
    int cell_id;
    int material_id;
    int zaid;
    char table[24];            /* lower/nearest table; empty for omission */
    char upper_table[24];      /* upper table for interpolation */
    double requested_kelvin;   /* zero for omission */
    double selected_kelvin;    /* lower/nearest temperature; zero for omission */
    double upper_kelvin;       /* upper temperature for interpolation */
    double upper_fraction;     /* interpolation weight of upper table */
    double omitted_abundance;  /* natural atom abundance; zero otherwise */
} alea_nuc_binding_notice_t;

typedef struct {
    size_t cell_index;
    int cell_id;
    int material_id;
    alea_nuc_particle_t particle;
    alea_nuc_table_type_t table_type;
    int zaid;                 /* input isotope ZAID */
    char table[24];          /* selected table */
    char upper_table[24];    /* second table for interpolation, if any */
    double requested_kelvin; /* zero if the cell did not specify temperature */
    double selected_kelvin;
    double upper_kelvin;
    double upper_fraction;
    double number_density;   /* atoms/barn-cm */
} alea_nuc_binding_selection_t;

/** NULL policy is strict: bracket interpolation is allowed for unsuffixed
 * nuclides; missing isotopes and out-of-range temperatures fail. Explicit
 * suffixes remain fixed. Thermal laws bind only to applicable neutron
 * components at matching table temperatures; aliases require an explicit
 * thermal_name_map entry. A policy enables only its bounded substitutions. */
alea_error_t alea_nuc_cell_bindings_prepare(
    struct alea_system* sys, alea_nuc_xsdir_t* xsdir, uint32_t particles,
    const alea_nuc_prepare_requirements_t* neutron_requirements,
    const alea_nuc_binding_policy_t* policy,
    alea_nuc_cell_bindings_t** output);
void alea_nuc_cell_bindings_free(alea_nuc_cell_bindings_t* bindings);
size_t alea_nuc_cell_bindings_notice_count(const alea_nuc_cell_bindings_t* bindings);
const alea_nuc_binding_notice_t* alea_nuc_cell_bindings_notice(
    const alea_nuc_cell_bindings_t* bindings, size_t index);
size_t alea_nuc_cell_bindings_selection_count(const alea_nuc_cell_bindings_t* bindings);
const alea_nuc_binding_selection_t* alea_nuc_cell_bindings_selection(
    const alea_nuc_cell_bindings_t* bindings, size_t index);

/** O(1) lookup by the terminal geometry cell index. Returns NULL for void
 * cells, unrequested particles, or an invalidated binding. Borrowed result. */
const alea_nuc_prepared_material_t* alea_nuc_cell_bindings_get(
    const alea_nuc_cell_bindings_t* bindings, size_t cell_index,
    alea_nuc_particle_t particle);
bool alea_nuc_cell_bindings_matches_geometry(
    const alea_nuc_cell_bindings_t* bindings, const struct alea_system* sys);

typedef struct {
    double position[3];
    alea_nuc_particle_state_t particle;
} alea_transport_source_t;

/** Sample one independent source particle for a global history ID. The
 * callback must initialize every output field and return an error if it
 * cannot sample. seed and history_id permit reproducible batch splitting.
 * The context is borrowed and the callback runs synchronously. */
typedef alea_error_t (*alea_transport_source_sampler_fn)(
    void* context, uint64_t seed, uint32_t history_id,
    alea_transport_source_t* output);

typedef struct {
    uint32_t histories;
    uint64_t seed;
    uint32_t max_events_per_history; /* shared by all descendants */
    double max_segment_distance; /* finite cm; distance-limit events continue */
    size_t max_pending_particles; /* zero uses a 1024-particle default */
    const alea_tally_plan_t* tally_plan; /* NULL disables configured tallies */
    uint32_t history_offset; /* first global history ID; zero by default */
} alea_transport_options_t;

typedef struct {
    size_t cell_count;
    uint32_t histories;
    uint64_t absorbed;
    uint64_t replaced; /* incident particles consumed by secondary production */
    uint64_t leaked;
    uint64_t collisions;
    uint64_t boundary_crossings;
    uint64_t reflections;
    uint64_t emitted_neutrons;
    uint64_t emitted_photons;
    uint64_t photon_collisions;
    uint64_t photon_absorbed;
    uint64_t photon_replaced;
    uint64_t photon_leaked;
    /* Neutron-only sum and sum of squares of each source history's weighted
     * cell path, retained for compatibility. Use tallies for photon paths.
     * Arrays have cell_count entries and are owned by this result. */
    double* track_length;
    double* track_length_squared;
    alea_tally_results_t* tallies; /* owned result; NULL without tally_plan */
} alea_transport_result_t;

typedef struct {
    uint32_t history_id;
    uint32_t particle_ordinal;
    uint32_t event_index;
    int cell_id;
    double position[3];
    double energy;
    alea_error_t error;
} alea_transport_failure_t;

/** Run one fixed neutron source repeatedly, including emitted neutron
 * descendants from reactions such as (n,2n) and fission when prepared.
 * Initialize output to zero before its first use and free a successful result
 * before reusing it. Bindings must be valid for every history. The bank and
 * event limits apply to each complete source history. Produced photons are
 * tracked when photon bindings are available. White/periodic boundaries
 * and URR sampling return explicit errors.
 * On failure output is zeroed; failure identifies the first incomplete
 * particle when tracking has begun. */
alea_error_t alea_transport_run_fixed_neutron(
    alea_system_t* sys,
    const alea_nuc_cell_bindings_t* bindings,
    const alea_transport_source_t* source,
    const alea_transport_options_t* options,
    alea_transport_result_t* output,
    alea_transport_failure_t* failure);

/** Fixed neutron or photon source with a shared descendant bank and tallies.
 * Material cells visited by each particle type need matching prepared
 * bindings; void cells require none. The neutron-only entry point above is
 * retained for existing callers. */
alea_error_t alea_transport_run_fixed_source(
    alea_system_t* sys,
    const alea_nuc_cell_bindings_t* bindings,
    const alea_transport_source_t* source,
    const alea_transport_options_t* options,
    alea_transport_result_t* output,
    alea_transport_failure_t* failure);

/** Run independently sampled neutron or photon source histories. Each
 * callback result is validated before transport. history_offset permits
 * reproducible nonoverlapping batches, up to UINT32_MAX as the last ID.
 * A callback failure identifies its history and leaves output zeroed. */
alea_error_t alea_transport_run_sampled_source(
    alea_system_t* sys,
    const alea_nuc_cell_bindings_t* bindings,
    alea_transport_source_sampler_fn sampler,
    void* source_context,
    const alea_transport_options_t* options,
    alea_transport_result_t* output,
    alea_transport_failure_t* failure);

void alea_transport_result_free(alea_transport_result_t* result);

#ifdef __cplusplus
}
#endif

#endif
