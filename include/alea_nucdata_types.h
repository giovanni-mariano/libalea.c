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
    int address;            /* start line (Type 1) or byte offset (Type 2) */
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
    ALEA_NUC_ELAW_LEVEL = 3,            /* level scattering */
    ALEA_NUC_ELAW_CONT_TABULAR = 4,     /* continuous tabular */
    ALEA_NUC_ELAW_GENERAL_EVAP = 5,     /* general evaporation */
    ALEA_NUC_ELAW_MAXWELL = 7,          /* Maxwell fission spectrum */
    ALEA_NUC_ELAW_EVAPORATION = 9,      /* evaporation spectrum */
    ALEA_NUC_ELAW_WATT = 11,            /* Watt fission spectrum */
    ALEA_NUC_ELAW_KALBACH = 44,         /* Kalbach-Mann */
    ALEA_NUC_ELAW_NBODY = 66,           /* N-body phase space */
    ALEA_NUC_ELAW_CORRELATED = 61,      /* correlated energy-angle */
} alea_nuc_energy_law_t;

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
    double level_A;             /* (A+1)/A factor squared */
    double level_Q;

    /* Maxwell/Evaporation/Watt parameters */
    int n_temp;
    double* temp_energy;        /* incident energy grid */
    double* temp_T;             /* nuclear temperature T(E) */
    double* temp_C;             /* restriction energy C(E) */
    double watt_a, watt_b;      /* Watt parameters (if constant) */

    /* Continuous tabular (law 4) and Kalbach-Mann (law 44) */
    struct {
        int n_ein;              /* number of incident energies */
        double* ein;            /* incident energy grid */
        int* n_eout;            /* number of outgoing energies per Ein */
        double** eout;          /* outgoing energy grids [n_ein][n_eout[i]] */
        double** pdf;           /* probability density [n_ein][n_eout[i]] */
        double** cdf;           /* cumulative distribution [n_ein][n_eout[i]] */
        /* Kalbach-Mann precompound fraction (law 44 only) */
        double** precompound_r; /* R values [n_ein][n_eout[i]], NULL for law 4 */
        double** precompound_a; /* a values [n_ein][n_eout[i]], NULL for law 4 */
        /* Correlated angular (law 61 only): LC locators per (ein, eout) */
        int** ang_lc;           /* angular locators [n_ein][n_eout[i]], NULL for law 4/44 */
        int dlw_base;           /* DLW block base for resolving angular locators */
    } tab;

    struct alea_nuc_energy_dist* next; /* linked list for multiple laws */
} alea_nuc_energy_dist_t;

/** Single reaction (MT) */
typedef struct {
    int mt;                     /* ENDF MT number */
    double q_value;             /* Q-value (MeV) */
    int ty;                     /* TYR value: yield, sign indicates ang. dist */
    int threshold_index;        /* first energy index (1-based) on main grid */
    int n_energies;             /* number of XS values */
    double* xs;                 /* cross-section array (on sub-grid) */

    alea_nuc_angular_dist_t* angular;
    alea_nuc_energy_dist_t* energy;
} alea_nuc_reaction_t;

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
    int n_energies;
    double* energy;
    double* nu;
} alea_nuc_nu_bar_t;

typedef struct {
    alea_nuc_nu_bar_t* total;        /* total ν̄ */
    alea_nuc_nu_bar_t* prompt;       /* prompt ν̄ (NULL if not given) */
    alea_nuc_nu_bar_t* delayed;      /* delayed ν̄ (NULL if not given) */
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
    double* energy;             /* incident energies */
    double* table;              /* n_energies × n_bands × 6 */
} alea_nuc_urr_t;

/* ============================================================================
 * PHOTON DATA (photoatomic .p tables)
 * ============================================================================ */

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

    /* Fluorescence */
    int n_fluorescence;
    /* TODO: fluorescence shell data */
} alea_nuc_photon_data_t;

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

#ifdef __cplusplus
}
#endif

#endif /* ALEA_NUCDATA_TYPES_H */
