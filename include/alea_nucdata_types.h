// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_NUCDATA_TYPES_H
#define ALEA_NUCDATA_TYPES_H

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>

#include "alea_types.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================================
 * PARTICLE AND TABLE TYPES
 * ============================================================================ */

typedef enum {
    ALEA_NUC_PARTICLE_NEUTRON = 0,
    ALEA_NUC_PARTICLE_PHOTON,
} alea_nuc_particle_t;

/** ACE table type, determined from ZAID suffix */
typedef enum {
    ALEA_NUC_TABLE_CONTINUOUS_NEUTRON = 0,   /* .XXc */
    ALEA_NUC_TABLE_PHOTOATOMIC,              /* .XXp */
    ALEA_NUC_TABLE_PHOTONUCLEAR,             /* .XXu */
    ALEA_NUC_TABLE_THERMAL_SAB,              /* .XXt */
    ALEA_NUC_TABLE_ELECTRON,                 /* .XXe */
} alea_nuc_table_type_t;

/* Use libalea's alea_error_t (from alea_types.h) for all error codes */

/* ============================================================================
 * XSDIR ENTRY
 * ============================================================================ */

/** One entry from xsdir/xsdata directory file */
typedef struct {
    char zaid[24];          /* e.g. "92235.80c" */
    double awr;             /* atomic weight ratio */
    char filename[512];     /* path to ACE file */
    int access_route;       /* 0 = default */
    int file_type;          /* 1 = ASCII, 2 = binary */
    int address;            /* start line (Type 1) or record number (Type 2) */
    int table_length;       /* number of words in XSS */
    int record_length;      /* 0 for ASCII */
    int num_entries;        /* entries per record, 0 for ASCII */
    double temperature;     /* kT in MeV */
    alea_nuc_table_type_t type;  /* parsed from ZAID suffix */
} alea_nuc_xsdir_entry_t;

/** Cached nuclide entry (nuclide is alea_nuc_nuclide_t*, stored as void*
 *  because the full type is defined later in this header) */
typedef struct {
    char zaid[24];          /* empty string = unused slot */
    void* nuclide;
} alea_nuc_cache_entry_t;

/** xsdir directory */
typedef struct {
    char datapath[512];         /* base path for relative filenames */
    alea_nuc_xsdir_entry_t* entries;
    size_t count;
    size_t capacity;

    /* Nuclide cache — open-addressing hash table, owned by the xsdir.
     * Nuclides loaded via alea_nuc_xsdir_get_nuclide() are cached here
     * and freed when the xsdir is destroyed. */
    alea_nuc_cache_entry_t* cache;
    size_t cache_count;         /* number of occupied slots */
    size_t cache_capacity;      /* total slots (always power of 2) */
} alea_nuc_xsdir_t;

/* ============================================================================
 * ACE TABLE (raw decoded data)
 * ============================================================================ */

#define ALEA_NUC_NXS_SIZE 16
#define ALEA_NUC_JXS_SIZE 32

/** Raw ACE table — header + XSS array */
typedef struct {
    /* Identity */
    char zaid[24];
    double awr;
    double temperature;     /* kT in MeV */
    char date[12];
    char comment[72];

    /* Index arrays (1-based in ACE, stored as-is) */
    int nxs[ALEA_NUC_NXS_SIZE];
    int jxs[ALEA_NUC_JXS_SIZE];

    /* Data array */
    double* xss;
    int xss_length;         /* NXS[1] */

    /* Legacy IZ/AW pairs */
    int iz[16];
    double aw[16];

    /* Table type (parsed from ZAID suffix) */
    alea_nuc_table_type_t type;

    /* Set true if any decode step attempted an out-of-bounds XSS access,
     * signalling corrupt/malformed data. The loader rejects the nuclide
     * rather than fabricating zero-valued physics. */
    bool decode_error;
    bool allocation_error;    /* allocation failed during nested decoding */
} alea_nuc_ace_table_t;

/* ============================================================================
 * NEUTRON REACTION
 * ============================================================================ */

/** Angular distribution type */
typedef enum {
    ALEA_NUC_ANG_ISOTROPIC = 0,
    ALEA_NUC_ANG_EQUIPROBABLE,       /* 32 equiprobable cosine bins */
    ALEA_NUC_ANG_TABULAR,            /* tabulated pdf/cdf */
} alea_nuc_angular_type_t;

/** Angular distribution at one incident energy */
typedef struct {
    alea_nuc_angular_type_t type;
    int interpolation;      /* ACE JJ: 1=histogram, 2=lin-lin */
    int n_cosines;
    double* cosine;         /* NULL for isotropic */
    double* pdf;
    double* cdf;
} alea_nuc_angular_point_t;

/** Angular distribution for a reaction */
typedef struct {
    int n_energies;
    double* energy;             /* incident energies */
    alea_nuc_angular_point_t* data;  /* one per energy point */
} alea_nuc_angular_dist_t;

/** Energy distribution law identifier */
typedef enum {
    ALEA_NUC_ELAW_DISCRETE_PHOTON = 2,  /* discrete secondary photon */
    ALEA_NUC_ELAW_LEVEL = 3,            /* level scattering */
    ALEA_NUC_ELAW_CONT_TABULAR = 4,     /* continuous tabular */
    ALEA_NUC_ELAW_GENERAL_EVAP = 5,     /* general evaporation */
    ALEA_NUC_ELAW_MAXWELL = 7,          /* Maxwell fission spectrum */
    ALEA_NUC_ELAW_EVAPORATION = 9,      /* evaporation spectrum */
    ALEA_NUC_ELAW_WATT = 11,            /* Watt fission spectrum */
    ALEA_NUC_ELAW_KALBACH = 44,         /* Kalbach-Mann */
    ALEA_NUC_ELAW_NBODY = 66,           /* N-body phase space */
    ALEA_NUC_ELAW_CORRELATED = 61,      /* correlated energy-angle */
    ALEA_NUC_ELAW_LAB_ANGLE_ENERGY = 67,/* laboratory angle then energy */
} alea_nuc_energy_law_t;

typedef struct {
    int interpolation;
    int n_points;
    double* energy;
    double* pdf;
    double* cdf;
} alea_nuc_law67_energy_t;

typedef struct {
    int interpolation;
    int n_cosines;
    double* cosine;
    alea_nuc_law67_energy_t* spectrum; /* one conditional spectrum per cosine */
} alea_nuc_law67_incident_t;

/** Energy distribution (can be a chain of laws with probability) */
typedef struct alea_nuc_energy_dist {
    alea_nuc_energy_law_t law;

    /* Applicability (probability of this law vs next) */
    int n_regions;              /* interpolation regions */
    int* nbt;                   /* breakpoints */
    int* interp;                /* interpolation types */
    int n_energies;
    double* energy;             /* incident energies for applicability */
    double* probability;

    /* Law-specific data stored as raw XSS slice for now */
    double* data;
    int data_length;

    /* Level scattering (law 3) */
    double level_A;             /* threshold energy */
    double level_Q;             /* outgoing/incident mass-ratio factor */

    /* Discrete photon (law 2) */
    int discrete_photon_primary;
    double discrete_photon_energy;
    double discrete_photon_awr;

    /* Maxwell/Evaporation/Watt parameters */
    int n_temp_regions;         /* interpolation regions for T(E) / Watt a(E) */
    int* temp_nbt;
    int* temp_interp;
    int n_temp;
    double* temp_energy;        /* incident energy grid */
    double* temp_T;             /* nuclear temperature T(E) */
    double* temp_C;             /* Watt b(E), or restriction energy C(E) */
    int n_watt_b;               /* number of points in Watt b(E) */
    int n_watt_b_regions;       /* interpolation regions for Watt b(E) */
    int* watt_b_nbt;
    int* watt_b_interp;
    double* watt_b_energy;      /* incident energy grid for Watt b(E) */
    double watt_a, watt_b;      /* Watt parameters (if constant) */
    double restriction_energy;  /* signed U; outgoing upper limit is E-U */

    /* General evaporation law 5: equiprobable dimensionless X boundaries */
    int n_general_evap;
    double* general_evap_x;

    /* N-body phase space (law 66) */
    int nbody_particles;
    double nbody_total_mass;
    double nbody_target_awr;
    double nbody_q_value;

    /* Continuous tabular (law 4) and Kalbach-Mann (law 44) */
    struct {
        int n_ein;              /* number of incident energies */
        int n_regions;          /* interpolation regions for incident energy */
        int* nbt;
        int* interp;
        double* ein;            /* incident energy grid */
        int* interpolation;     /* outgoing interpolation code per Ein */
        int* n_discrete;        /* discrete outgoing lines per Ein */
        int* n_eout;            /* number of outgoing energies per Ein */
        double** eout;          /* outgoing energy grids [n_ein][n_eout[i]] */
        double** pdf;           /* probability density [n_ein][n_eout[i]] */
        double** cdf;           /* cumulative distribution [n_ein][n_eout[i]] */
        /* Kalbach-Mann precompound fraction (law 44 only) */
        double** precompound_r; /* R values [n_ein][n_eout[i]], NULL for law 4 */
        double** precompound_a; /* a values [n_ein][n_eout[i]], NULL for law 4 */
        /* Correlated angular (law 61 only): LC locators per (ein, eout) */
        int** ang_lc;           /* angular locators [n_ein][n_eout[i]], NULL for law 4/44 */
        alea_nuc_angular_point_t** correlated_mu; /* decoded law 61 angles */
        int dlw_base;           /* DLW block base for resolving angular locators */
    } tab;

    /* Law 67: incident energy -> equiprobable cosine -> conditional energy. */
    struct {
        int n_ein;
        int n_regions;
        int* nbt;
        int* interp;
        double* ein;
        alea_nuc_law67_incident_t* incident;
    } law67;

    struct alea_nuc_energy_dist* next; /* linked list for multiple laws */
} alea_nuc_energy_dist_t;

/** Single reaction (MT) */
typedef struct {
    int mt;                     /* ENDF MT number */
    double q_value;             /* Q-value (MeV) */
    int ty;                     /* TYR value: yield, sign indicates ang. dist */
    bool center_of_mass;        /* negative TYR: distribution is in CM frame */
    int threshold_index;        /* first energy index (1-based) on main grid */
    int n_energies;             /* number of XS values */
    double* xs;                 /* cross-section array (on sub-grid) */

    alea_nuc_angular_dist_t* angular;
    alea_nuc_energy_dist_t* energy;
} alea_nuc_reaction_t;

/** One neutron-induced photon-production channel from ACE MTRP/SIGP. */
typedef struct {
    int mt;                     /* encoded photon-production MT */
    int parent_mt;              /* incident-neutron reaction MT */
    int mf;                     /* source ENDF file: 12, 13, or 16 */
    bool production_xs;         /* values are production XS rather than yield */
    int threshold_index;        /* 1-based main-grid index for MF=13 */
    int n_regions;
    int* nbt;
    int* interp;
    int n_energies;
    double* energy;             /* incident grid for MF=12/16 */
    double* values;             /* yield or photon-production XS */
    alea_nuc_angular_dist_t* angular;
    alea_nuc_energy_dist_t* spectrum;
} alea_nuc_photon_production_t;

/* ============================================================================
 * FISSION DATA
 * ============================================================================ */

typedef enum {
    ALEA_NUC_NU_POLYNOMIAL = 0,
    ALEA_NUC_NU_TABULAR,
} alea_nuc_nu_type_t;

typedef struct {
    alea_nuc_nu_type_t type;

    /* Polynomial: ν̄ = c0 + c1*E + c2*E² + ... */
    int n_coeffs;
    double* coeffs;

    /* Tabular */
    int n_regions;
    int* nbt;
    int* interp;
    int n_energies;
    double* energy;
    double* nu;
} alea_nuc_nu_bar_t;

typedef struct {
    double decay_rate;          /* s^-1 */
    int n_regions;
    int* nbt;
    int* interp;
    int n_energies;
    double* energy;
    double* probability;
    alea_nuc_energy_dist_t* spectrum;
} alea_nuc_delayed_group_t;

typedef struct {
    alea_nuc_nu_bar_t* total;        /* total ν̄ */
    alea_nuc_nu_bar_t* prompt;       /* prompt ν̄ (NULL if not given) */
    alea_nuc_nu_bar_t* delayed;      /* delayed ν̄ (NULL if not given) */
    int n_delayed_groups;
    alea_nuc_delayed_group_t* delayed_groups;
} alea_nuc_fission_t;

/* ============================================================================
 * URR PROBABILITY TABLES
 * ============================================================================ */

typedef struct {
    int n_energies;
    int n_bands;                /* number of probability bands */
    int interp;                 /* interpolation flag */
    int inelastic_flag;         /* treatment of inelastic */
    int absorption_flag;
    bool multiply_smooth;       /* table values multiply smooth cross sections */
    double* energy;             /* incident energies */
    double* table;              /* n_energies × n_bands × 6 */
} alea_nuc_urr_t;

/* ============================================================================
 * PHOTON DATA (photoatomic .p tables)
 * ============================================================================ */

typedef struct {
    int primary_designator;    /* shell receiving the first vacancy */
    int secondary_designator;  /* 0 for radiative, otherwise second vacancy */
    double energy;             /* tabulated photon/electron energy, MeV */
    double cumulative_probability;
} alea_nuc_atomic_transition_t;

typedef struct {
    int designator;
    double occupancy;
    double binding_energy;              /* MeV */
    double compton_vacancy_probability;
    int n_transitions;
    alea_nuc_atomic_transition_t* transitions;
    double* ln_photoelectric_xs;         /* common photon energy grid */
    size_t max_relaxation_photons;       /* exact cascade capacity bound */
} alea_nuc_atomic_subshell_t;

typedef struct {
    double electron_count;
    double binding_energy;              /* MeV */
    double cumulative_probability;
    int profile_index;
} alea_nuc_compton_shell_t;

typedef struct {
    int interpolation;
    int n_momenta;
    double* momentum;                    /* p_z in atomic units */
    double* pdf;
    double* cdf;
} alea_nuc_compton_profile_t;

typedef struct {
    int n_energies;
    double* energy;
    double* sigma_incoherent;   /* Compton */
    double* sigma_coherent;     /* Rayleigh */
    double* sigma_photoelectric;
    double* sigma_pair;         /* pair + triplet production */
    double* heating;

    /* Pre-stored log-space arrays for fast log-log interpolation */
    double* ln_energy;
    double* ln_sigma_incoherent;
    double* ln_sigma_coherent;
    double* ln_sigma_photoelectric;
    double* ln_sigma_pair;          /* -HUGE_VAL sentinel for below threshold */

    /* Form factors */
    int n_incoherent_ff;
    double* incoherent_momentum;    /* momentum transfer values */
    double* incoherent_ff;          /* S(q,Z) scattering function */

    int n_coherent_ff;
    double* coherent_momentum;
    double* coherent_ff;            /* F(q,Z) form factor */
    double* coherent_ff_cumulative; /* integrated form factor */

    /* Detailed EPR subshell photoelectric and atomic-relaxation data. */
    int epr_format;
    int n_compton_shells;
    alea_nuc_compton_shell_t* compton_shells;
    int n_compton_profiles;
    alea_nuc_compton_profile_t* compton_profiles;
    int n_subshells;
    alea_nuc_atomic_subshell_t* subshells;

    /* Averaged fluorescence representation in older photoatomic tables. */
    int n_fluorescence;
    double* fluorescence_edge;   /* incident-energy edge, MeV */
    double* fluorescence_phi;    /* cumulative edge-jump parameter */
    double* fluorescence_yield;  /* cumulative photon yield */
    double* fluorescence_energy; /* representative photon energy, MeV */
} alea_nuc_photon_data_t;

/* ============================================================================
 * THERMAL SCATTERING DATA (.t tables)
 * ============================================================================ */

typedef enum {
    ALEA_NUC_THERMAL_ELASTIC_NONE = 0,
    ALEA_NUC_THERMAL_ELASTIC_INCOHERENT = 3,
    ALEA_NUC_THERMAL_ELASTIC_COHERENT = 4,
    ALEA_NUC_THERMAL_ELASTIC_MIXED = 5,
} alea_nuc_thermal_elastic_mode_t;

/** Decoded thermal ACE table with discrete or continuous inelastic data. */
typedef struct {
    char zaid[24];
    double awr;
    double temperature;          /* kT in MeV */

    int n_applicable_zaids;
    int applicable_zaids[16];    /* ACE IZ identifiers */

    int n_inelastic_energies;
    double* inelastic_energy;    /* incident energy grid, MeV */
    double* inelastic_xs;        /* barns */
    int n_inelastic_outgoing;
    int n_inelastic_cosines;
    bool inelastic_continuous;   /* ACE IFENG=2 */
    bool inelastic_skewed;
    int n_inelastic_outgoing_total;
    int* inelastic_outgoing_offset; /* [incident + 1], continuous only */
    double* inelastic_energy_out; /* flattened outgoing energies, MeV */
    double* inelastic_pdf;        /* continuous only, MeV^-1 */
    double* inelastic_cdf;        /* continuous only */
    double* inelastic_mu;         /* [flattened outgoing][cosine] */

    alea_nuc_thermal_elastic_mode_t elastic_mode;
    int n_coherent_edges;
    double* coherent_edge;       /* Bragg edges, MeV */
    double* coherent_factor;     /* cumulative structure factor, MeV*b */

    int n_incoherent_energies;
    double* incoherent_energy;   /* MeV */
    double* incoherent_xs;       /* barns */
    int n_incoherent_cosines;
    double* incoherent_mu;       /* [incident][cosine] */
} alea_nuc_thermal_t;

/* ============================================================================
 * NUCLIDE — fully decoded ACE table
 * ============================================================================ */

typedef struct {
    char zaid[24];
    int Z, A, metastable;
    alea_nuc_particle_t particle;
    double awr;
    double temperature;         /* kT in MeV */

    /* Main energy grid */
    int n_energies;
    double* energy;             /* MeV, ascending */
    double* sigma_total;
    double* sigma_abs;
    double* sigma_elastic;
    double* heating;

    /* Elastic angular distribution (decoded from LAND[0]) */
    alea_nuc_angular_dist_t* elastic_angular;

    /* Reactions (non-elastic) */
    int n_reactions;
    alea_nuc_reaction_t* reactions;

    int n_photon_productions;
    alea_nuc_photon_production_t* photon_productions;
    double* total_photon_production_xs; /* GPD, on the main neutron grid */
    int n_photon_yield_multipliers;
    int* photon_yield_multipliers;      /* YP neutron MT identifiers */

    /* Fission (NULL if non-fissile) */
    alea_nuc_fission_t* fission;

    /* URR probability tables (NULL if none) */
    alea_nuc_urr_t* urr;

    /* Photon data (NULL for neutron tables) */
    alea_nuc_photon_data_t* photon;

    /* MT → reaction index lookup table (-1 = absent), size ALEA_NUC_MT_TABLE_SIZE */
#define ALEA_NUC_MT_TABLE_SIZE 1000
    int* mt_to_rxn;

    /* Raw ACE table (kept for lazy decode of angular/energy distributions) */
    alea_nuc_ace_table_t raw;
} alea_nuc_nuclide_t;

/* ============================================================================
 * MATERIAL — composition of nuclides
 * ============================================================================ */

typedef struct {
    alea_nuc_nuclide_t* nuclide;
    double number_density;      /* atoms/barn-cm */
} alea_nuc_mat_component_t;

typedef struct {
    alea_nuc_mat_component_t* components;
    int n_components;
    int capacity;
} alea_nuc_material_t;

/* ============================================================================
 * MULTIGROUP DATA
 * ============================================================================ */

/**
 * Weighting spectrum function for multigroup collapse.
 * @param E   Energy in MeV
 * @param ctx User context pointer (passed through from mg->spectrum_ctx)
 * @return    φ(E), the weighting flux at energy E (arbitrary units)
 */
typedef double (*alea_nuc_spectrum_fn)(double E, void* ctx);

/** Multigroup cross sections and scattering matrix */
typedef struct {
    int n_groups;               /* number of energy groups */
    double* bounds;             /* group boundaries [n_groups+1], descending */
    double* sigma_t;            /* total XS per group */
    double* sigma_a;            /* absorption XS per group */
    double* sigma_s;            /* elastic scattering XS per group */
    double* sigma_f;            /* fission XS per group */
    double* nu_sigma_f;         /* ν̄·σ_f per group */
    double* chi;                /* fission spectrum per group (sums to 1) */
    double* scatter;            /* scattering matrix [g_from * n_groups + g_to] */

    /* User-defined weighting spectrum (NULL = default Maxwellian+1/E+fission) */
    alea_nuc_spectrum_fn spectrum;
    void* spectrum_ctx;
} alea_nuc_multigroup_t;

/* ============================================================================
 * REACTION CLASSIFICATION
 * ============================================================================ */

/** Reaction class for neutron transport */
typedef enum {
    ALEA_NUC_RXN_ABSORPTION = 0,     /* no secondary neutrons (capture, (n,α), etc.) */
    ALEA_NUC_RXN_SCATTER,            /* one neutron out (elastic, inelastic levels) */
    ALEA_NUC_RXN_MULTIPLY,           /* multiple neutrons out ((n,2n), (n,3n), fission) */
} alea_nuc_reaction_class_t;

/* ============================================================================
 * INTERACTION RESULT
 * ============================================================================ */

typedef struct {
    int mt;                     /* reaction MT */
    double energy_out;          /* outgoing energy (MeV) */
    double mu;                  /* scattering cosine */
    int n_secondary;            /* secondary neutrons (from TYR / ν̄) */
    double weight_factor;       /* 1.0 forward, correction for adjoint */
} alea_nuc_interaction_t;

/* ============================================================================
 * PREPARED CONTINUOUS-ENERGY COLLISION PHYSICS
 * ============================================================================ */

typedef struct alea_nuc_prepared_material alea_nuc_prepared_material_t;

typedef enum {
    ALEA_NUC_CAP_STATIONARY_ELASTIC = 1u << 0,
    ALEA_NUC_CAP_ABSORPTION         = 1u << 1,
    ALEA_NUC_CAP_FREE_GAS           = 1u << 2,
    ALEA_NUC_CAP_THERMAL_SAB        = 1u << 3,
    ALEA_NUC_CAP_FISSION            = 1u << 4,
    ALEA_NUC_CAP_PHOTON             = 1u << 5,
    ALEA_NUC_CAP_URR                = 1u << 6,
    ALEA_NUC_CAP_NEUTRON_EMISSION   = 1u << 7,
    ALEA_NUC_CAP_DELAYED_NEUTRON    = 1u << 8,
    ALEA_NUC_CAP_PHOTON_PRODUCTION  = 1u << 9,
} alea_nuc_capability_t;

#define ALEA_NUC_CAP_RESTRICTED_NEUTRON \
    (ALEA_NUC_CAP_STATIONARY_ELASTIC | ALEA_NUC_CAP_ABSORPTION)

#define ALEA_NUC_CAP_CONTINUOUS_NEUTRON \
    (ALEA_NUC_CAP_RESTRICTED_NEUTRON | ALEA_NUC_CAP_NEUTRON_EMISSION)

typedef enum {
    ALEA_NUC_PREP_OK = 0,
    ALEA_NUC_PREP_EMPTY_MATERIAL,
    ALEA_NUC_PREP_UNSUPPORTED_CAPABILITY,
    ALEA_NUC_PREP_UNSUPPORTED_PARTICLE,
    ALEA_NUC_PREP_UNSUPPORTED_URR,
    ALEA_NUC_PREP_UNSUPPORTED_REACTION,
    ALEA_NUC_PREP_INVALID_ENERGY_DISTRIBUTION,
    ALEA_NUC_PREP_INVALID_ANGULAR,
    ALEA_NUC_PREP_INVALID_CROSS_SECTIONS,
    ALEA_NUC_PREP_INVALID_THERMAL_ASSOCIATION,
} alea_nuc_prepare_issue_t;

/** Bind one material component to a thermal scattering table. */
typedef struct {
    int component_index;
    const alea_nuc_thermal_t* thermal;
} alea_nuc_thermal_association_t;

typedef struct {
    uint32_t required_capabilities;
    const alea_nuc_thermal_association_t* thermal_associations;
    size_t n_thermal_associations;
    double thermal_temperature_tolerance; /* absolute kT tolerance, MeV */
} alea_nuc_prepare_requirements_t;

typedef struct {
    uint32_t available_capabilities;
    uint32_t missing_capabilities;
    alea_nuc_prepare_issue_t issue;
    int component_index;    /* -1 if the issue is not component-specific */
    int mt;                 /* 0 if the issue is not reaction-specific */
    int law;                /* 0 if the issue is not distribution-specific */
    char detail[192];
} alea_nuc_capability_report_t;

typedef double (*alea_nuc_random_fn)(void* context);

/** Caller-owned Philox address. Initialize before passing to the RNG callback. */
typedef struct {
    uint64_t seed;
    uint64_t entity_id;
    uint64_t local_draw;
    uint32_t event_index;
    uint32_t domain;
} alea_nuc_rng_t;

typedef enum {
    ALEA_NUC_RNG_FLIGHT = 0,
    ALEA_NUC_RNG_COLLISION,
    ALEA_NUC_RNG_URR
} alea_nuc_rng_domain_t;

typedef struct {
    alea_nuc_particle_t type;
    double energy;          /* MeV */
    double direction[3];    /* unit vector */
    double weight;
    double time;            /* seconds */
} alea_nuc_particle_state_t;

typedef struct {
    alea_nuc_particle_state_t outgoing;
    double mu_cm;
    double mu_lab;
} alea_nuc_free_gas_result_t;

typedef struct {
    bool active;
    double factors[5]; /* total, elastic, fission, capture, heating */
} alea_nuc_urr_sample_t;

typedef struct {
    alea_nuc_urr_sample_t* components;
    size_t capacity;
} alea_nuc_evaluation_workspace_t;

typedef struct {
    const alea_nuc_prepared_material_t* prepared;
    alea_nuc_particle_state_t incident;
    double macro_total;       /* cm^-1 */
    double macro_elastic;     /* cm^-1 */
    double macro_thermal;     /* cm^-1, replaces free-atom elastic in range */
    double macro_absorption;  /* cm^-1 */
    double macro_neutron_emission; /* cm^-1 */
    const alea_nuc_evaluation_workspace_t* workspace;
} alea_nuc_evaluation_t;

typedef enum {
    ALEA_NUC_OUTCOME_ABSORBED = 0,
    ALEA_NUC_OUTCOME_SCATTERED,
    ALEA_NUC_OUTCOME_REPLACED,
} alea_nuc_collision_outcome_t;

typedef struct {
    alea_nuc_particle_state_t* particles;
    size_t capacity;
    size_t count;
} alea_nuc_secondary_buffer_t;

typedef struct {
    alea_nuc_collision_outcome_t outcome;
    int component_index;
    int mt;
    double mu_cm;
    double mu_lab;
    bool deposition_available;
    double local_energy_deposition; /* MeV */
    alea_nuc_particle_state_t outgoing; /* valid when outcome is SCATTERED */
    size_t n_emitted; /* particles appended to the caller's secondary buffer */
} alea_nuc_collision_result_t;

#ifdef __cplusplus
}
#endif

#endif /* ALEA_NUCDATA_TYPES_H */
