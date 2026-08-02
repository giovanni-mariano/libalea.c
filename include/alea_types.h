// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_TYPES_H
#define ALEA_TYPES_H

#include <stdint.h>
#include <stdbool.h>
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @file alea_types.h
 * @brief Public type definitions for CSG system
 * 
 * All fundamental types needed by users of the CSG API.
 * This header merges primitive types and CSG node types into one place.
 */


 /* ============================================================================
 * OPAQUE TYPES
 * ============================================================================ */

/** @brief Main CSG system handle (opaque) */
typedef struct alea_system alea_system_t;

/** @brief Export options handle (opaque) */
typedef struct alea_export_options alea_export_options_t;

/** @brief Void generation result handle (opaque) */
typedef struct void_result void_result_t;


/* ============================================================================
 * PUBLIC ID TYPES
 * ============================================================================ */

/**
 * @brief Opaque handle to a CSG tree node
 * 
 * Nodes are stored in an array and accessed by ID.
 * This allows the array to grow without invalidating pointers.
 */
typedef uint32_t alea_node_id_t;

/**
 * @brief Opaque handle to a primitive
 * 
 * Primitives are stored in a deduplicated array.
 */
typedef uint32_t alea_primitive_id_t;

/**
 * @brief Material identifier
 */
typedef uint32_t alea_material_id_t;

/**
 * @brief Surface identifier (for MCNP conversion)
 */
typedef uint32_t alea_surface_id_t;

typedef struct alea_flatten_config alea_flatten_config_t;

// ============================================================================
// PRIMITIVE TYPES
// ============================================================================

/**
 * @brief Types of geometric primitives supported by the CSG system
 */
typedef enum {
    ALEA_PRIMITIVE_PLANE = 1,
    ALEA_PRIMITIVE_SPHERE,
    ALEA_PRIMITIVE_CYLINDER_X,
    ALEA_PRIMITIVE_CYLINDER_Y,
    ALEA_PRIMITIVE_CYLINDER_Z,
    ALEA_PRIMITIVE_CONE_X,
    ALEA_PRIMITIVE_CONE_Y,
    ALEA_PRIMITIVE_CONE_Z,
    ALEA_PRIMITIVE_RPP,          // RPP - axis-aligned box
    ALEA_PRIMITIVE_QUADRIC,
    ALEA_PRIMITIVE_TORUS_X,
    ALEA_PRIMITIVE_TORUS_Y,
    ALEA_PRIMITIVE_TORUS_Z,
    // Macrobodies
    ALEA_PRIMITIVE_RCC,          // Right Circular Cylinder
    ALEA_PRIMITIVE_BOX,          // BOX - general oriented box
    ALEA_PRIMITIVE_SPH,          // SPH - sphere macrobody
    ALEA_PRIMITIVE_TRC,          // Truncated Right Cone
    ALEA_PRIMITIVE_ELL,          // Ellipsoid
    ALEA_PRIMITIVE_REC,          // Right Elliptical Cylinder
    ALEA_PRIMITIVE_WED,          // Wedge
    ALEA_PRIMITIVE_RHP,          // Right Hexagonal Prism
    ALEA_PRIMITIVE_ARB           // Arbitrary Polyhedron
} alea_primitive_type_t;

/**
 * @brief Axis enumeration for axis-aligned primitives
 */
typedef enum {
    ALEA_AXIS_X = 0,
    ALEA_AXIS_Y = 1,
    ALEA_AXIS_Z = 2
} alea_axis_t;

/* ============================================================================
 * EXPORT FORMATS
 * ============================================================================ */

typedef enum {
    ALEA_EXPORT_FORMAT_MCNP = 0,
    ALEA_EXPORT_FORMAT_OPENMC,
    ALEA_EXPORT_FORMAT_SERPENT,
} alea_export_format_t;

typedef enum {
    ALEA_EXPORT_PREFER_MACROBODY,    /* Emit RCC/RPP if target supports */
    ALEA_EXPORT_FORCE_PRIMITIVES,    /* Always decompose to primitive surfaces */
} alea_export_policy_t;

/**
 * @brief Surface emission policy for macrobodies
 * 
 * Controls how macrobodies (RCC, RPP) are emitted during export:
 * - ALEA_EMIT_MACROBODY: Emit as single macrobody surface (MCNP native)
 * - ALEA_EMIT_SURFACES: Decompose to constituent primitive surfaces
 */
typedef enum {
    ALEA_EMIT_MACROBODY,    /* Emit as-is (RCC, RPP) */
    ALEA_EMIT_SURFACES,     /* Decompose to cylinder+planes, 6 planes, etc. */
} alea_surface_emit_policy_t;

/* ============================================================================
 * UNIFIED CONFIGURATION
 * ============================================================================ */

typedef struct {
    /* Tolerance */
    double abs_tol;             /* 1e-14 */
    double rel_tol;             /* 1e-12 */
    double zero_threshold;      /* 1e-14 */

    /* Behavior */
    bool dedup;                 /* true  */
    int log_level;              /* ALEA_LOG_WARN (2) */

    /* Export */
    bool export_materials;      /* true  */
    bool export_transforms;     /* true  */
    int universe_depth;         /* -1 (all) */
    int fill_depth;             /* 0 (no expansion) */

    /* Void generation */
    int void_max_depth;         /* 8    */
    double void_min_size;       /* 0.1  */
    int void_probes_per_axis;   /* 3    */

    /* Void merge */
    double merge_cell_weight;   /* 1.0  */
    double merge_surface_weight;/* 0.1  */
    int merge_max_surfaces;     /* 24   */
    int merge_min_cells;        /* 1    */
    bool merge_use_greedy;      /* false */
    int void_consolidate;       /* 100 (0 = off) */

    /* Flatten */
    int flatten_max_depth;      /* 0 = unlimited */
} alea_config_t;

extern const alea_config_t ALEA_CONFIG_DEFAULT;

// ============================================================================
// CSG OPERATION TYPES

/**
 * @brief CSG operation types (used in packed type_and_flags)
 */
typedef enum {
    ALEA_OP_PRIMITIVE = 0,      // Leaf node (primitive)
    ALEA_OP_UNION,              // Boolean union (OR)
    ALEA_OP_INTERSECTION,       // Boolean intersection (AND)
    ALEA_OP_DIFFERENCE,          // Boolean difference (A - B)
    ALEA_OP_COMPLEMENT          // Boolean complement (NOT A)
} alea_operation_t;


typedef enum {
    ALEA_BOUNDARY_TRANSMISSIVE = 0,
    ALEA_BOUNDARY_REFLECTIVE,    // '*' prefix in MCNP
    ALEA_BOUNDARY_WHITE,         // '+' prefix in MCNP
    ALEA_BOUNDARY_PERIODIC,      // Paired surfaces
    ALEA_BOUNDARY_VACUUM         // Particles escape (killed)
} alea_boundary_type_t;

// ============================================================================
// SPECIAL VALUES
// ============================================================================

#define ALEA_NODE_ID_INVALID      UINT32_MAX
#define ALEA_NODE_ID_LIKE_PENDING (UINT32_MAX - 1)  /* LIKE cell awaiting resolution */
#define ALEA_PRIMITIVE_ID_INVALID UINT32_MAX
#define ALEA_MATERIAL_NONE        0

/** @brief Error codes */
typedef enum {
    /* Success */
    ALEA_OK = 0,

    /* Input validation errors */
    ALEA_ERR_NULL_ARG,           /**< NULL pointer where not allowed */
    ALEA_ERR_INVALID_ID,         /**< Node/primitive/cell ID out of range */
    ALEA_ERR_INVALID_ARG,        /**< Argument value out of valid range */
    ALEA_ERR_INVALID_STATE,      /**< Operation not valid in current state */

    /* Resource errors */
    ALEA_ERR_OUT_OF_MEMORY,      /**< Memory allocation failed */
    ALEA_ERR_FILE_NOT_FOUND,     /**< File does not exist */
    ALEA_ERR_FILE_READ,          /**< Error reading file */
    ALEA_ERR_FILE_WRITE,         /**< Error writing file */

    /* Operation errors */
    ALEA_ERR_PARSE_ERROR,        /**< Syntax/parse error */
    ALEA_ERR_UNSUPPORTED,        /**< Feature not supported */
    ALEA_ERR_UNSUPPORTED_SURFACE,/**< Surface type not supported (legacy alias) */
    ALEA_ERR_EXPORT_FAILED,      /**< Export operation failed */
    ALEA_ERR_NOT_IMPLEMENTED,    /**< Feature not yet implemented */

    /* Interruption */
    ALEA_ERR_INTERRUPTED,        /**< Operation interrupted by user (SIGINT) */

    /* Query results (not strictly errors) */
    ALEA_ERR_NOT_FOUND,          /**< Item not found (valid result, not error) */
    ALEA_ERR_EMPTY,              /**< Collection is empty */
    ALEA_ERR_OVERFLOW            /**< Buffer too small, result truncated */
} alea_error_t;

/* --- Error Functions --- */

/**
 * Get static error description string
 * @param error Error code
 * @return Static string describing the error (do not free)
 */
const char* alea_error_string(alea_error_t error);

/**
 * Set detailed error message (thread-local)
 * @param code Error code
 * @param fmt Printf-style format string
 */
void alea_set_error_detail(alea_error_t code, const char* fmt, ...);

/**
 * Get detailed error message
 * @return Detailed message if set, otherwise static string from error code
 */
const char* alea_get_error_detail(void);

/**
 * Get current error code (thread-local)
 * @return Last error code set via alea_set_error_detail
 */
alea_error_t alea_get_last_error(void);

/**
 * Clear error detail
 */
void alea_clear_error_detail(void);

/* ============================================================================
 * DATA STRUCTURES (Read-only views for queries)
 * ============================================================================ */

/** Resolution flags for alea_cell_hit_t / alea_ray_segment_t.
 *  UNDEFINED_FILL: the answer is a fill or lattice container cell whose
 *  filling universe has no containing cell at the query point (an MCNP
 *  undefined region). The container is the deepest cell that actually
 *  contains the point; the flag marks that the fill chain below it is
 *  unresolved. */
#define ALEA_RESOLVE_UNDEFINED_FILL 0x01

/** @brief Point query result for hierarchical traversal */
typedef struct {
    int cell_id;
    int cell_index;
    int material_id;
    int universe_id;
    int fill_universe;
    int depth;
    double local_x, local_y, local_z;
    uint8_t resolution_flags;  /**< ALEA_RESOLVE_* bits (0 = fully resolved) */
} alea_cell_hit_t;

/** @brief Axis-aligned bounding box */
typedef struct {
    double min_x, max_x;
    double min_y, max_y;
    double min_z, max_z;
} alea_bbox_t;

/** @brief Cell information (read-only view) */
typedef struct {
    int cell_id;            /**< MCNP cell ID */
    int material_id;        /**< Material number (0 = void) */
    double density;         /**< Material density (always positive) */
    bool is_mass_density;   /**< true = g/cm³, false = atoms/b-cm */
    int universe_id;        /**< Universe this cell belongs to */
    int fill_universe;      /**< Universe filling this cell (-1 = none) */
    int fill_transform;     /**< Transform applied to fill (0 = none) */
    alea_node_id_t root;     /**< Root node of CSG tree */
    alea_bbox_t bbox;        /**< Bounding box */
    double temperature;     /**< Temperature in Kelvin (0.0 = not set) */
    bool has_temperature;   /**< true if temperature was explicitly set */
    int lat_type;           /**< Lattice type: 0=none, 1=rect, 2=hex */
    int lat_fill_dims[6];   /**< Lattice dimensions: imin,imax,jmin,jmax,kmin,kmax */
    const int* lat_fill;    /**< Array of universe IDs (NULL if no lattice) */
    size_t lat_fill_count;  /**< Number of elements in lat_fill */
    int lat_fill_repeating; /**< Nonzero when a simple FILL=N repeats indefinitely */
    int lat_fill_zero_element_coords; /**< Fill universe uses the lattice's (0,0,0)-element coordinates */
    double lat_pitch[3];    /**< Element pitch in each dimension */
    double lat_lower_left[3]; /**< Lower-left corner of the lattice */
    const char* comments;     /**< "C" comment lines before cell (NULL if none) */
    const char* inline_comment; /**< Inline "$" comment (NULL if none) */
} alea_cell_info_t;

/** @brief Surface information (read-only view) */
typedef struct {
    int surface_id;         /**< MCNP surface ID */
    alea_primitive_type_t type;
    alea_node_id_t node_id;  /**< Node ID for this surface */
} alea_surface_info_t;

/**
 * @brief Surface entry with pre-allocated sense nodes
 *
 * Each surface gets exactly two nodes:
 * - pos_node: positive sense (outside, +S in MCNP)
 * - neg_node: negative sense (inside, -S in MCNP)
 */
typedef struct {
    int mc_surface_id;            /**< Surface ID (always positive) */
    alea_primitive_id_t primitive_id; /**< Primitive in GLOBAL coordinates */
    alea_node_id_t pos_node;         /**< +S node (sense = +1) */
    alea_node_id_t neg_node;         /**< -S node (sense = -1) */
    int transform_id;               /**< TRn applied (0 = none) */
    bool transform_applied;         /**< true if primitive is transformed */
    alea_boundary_type_t boundary_type;
    int periodic_surface_id;
    alea_node_id_t expanded_pos_node;  /**< Expanded exterior (for macrobodies) */
    alea_node_id_t expanded_neg_node;  /**< Expanded interior (for macrobodies) */
} alea_surface_entry_t;



typedef struct alea_stats alea_stats_t;


// ============================================================================
// PRIMITIVE DATA STRUCTURES
// ============================================================================

/**
 * @brief Plane: ax + by + cz + d = 0
 */
typedef struct {
    double a, b, c, d;
} alea_plane_data_t;

/**
 * @brief Sphere: (x-cx)² + (y-cy)² + (z-cz)² = r²
 */
typedef struct {
    double center_x, center_y, center_z;
    double radius;
} alea_sphere_data_t;

/**
 * @brief Cylinder along X-axis
 */
typedef struct {
    double center_y, center_z;     // Center in YZ plane
    double radius;
} alea_cylinder_x_data_t;

/**
 * @brief Cylinder along Y-axis
 */
typedef struct {
    double center_x, center_z;     // Center in XZ plane
    double radius;
} alea_cylinder_y_data_t;

/**
 * @brief Cylinder along Z-axis
 */
typedef struct {
    double center_x, center_y;     // Center in XY plane
    double radius;
} alea_cylinder_z_data_t;

/**
 * @brief Cone along X-axis
 */
typedef struct {
    double apex_x, apex_y, apex_z;
    double tan_angle_sq;            // tan²(half-angle)
    int sheet_selection;
} alea_cone_x_data_t;

/**
 * @brief Cone along Y-axis
 */
typedef struct {
    double apex_x, apex_y, apex_z;
    double tan_angle_sq;
    int sheet_selection;
} alea_cone_y_data_t;

/**
 * @brief Cone along Z-axis
 */
typedef struct {
    double apex_x, apex_y, apex_z;
    double tan_angle_sq;
    int sheet_selection;
} alea_cone_z_data_t;

/**
 * @brief Axis-aligned box
 */
typedef struct {
    // double x_min, x_max;
    // double y_min, y_max; 
    // double z_min, z_max;

    double min_x, max_x;
    double min_y, max_y;
    double min_z, max_z;
} alea_box_data_t;

/**
 * @brief General quadric surface
 * Ax² + By² + Cz² + Dxy + Eyz + Fxz + Gx + Hy + Iz + J = 0
 */
typedef struct {
    double coeffs[10];  // A, B, C, D, E, F, G, H, I, J
} alea_quadric_data_t;

/**
 * @brief Torus (all orientations use same structure)
 */
typedef struct {
    alea_axis_t axis;                        // Axis of revolution
    double center_x, center_y, center_z;    // Torus center
    double major_radius;                    // Distance from axis to tube center
    double minor_radius;                    // Tube radius
    double axial_semiwidth_B;               // For elliptical cross-section
} alea_torus_data_t;


/**
 * @brief Right Circular Cylinder (arbitrary orientation)
 * 
 * MCNP RCC macrobody: defined by base point, height vector, and radius.
 * The cylinder extends from base to base + height.
 */
typedef struct {
    double base_x, base_y, base_z;          // Base center point
    double height_x, height_y, height_z;    // Height vector (axis direction * length)
    double radius;
} alea_rcc_data_t;

/**
 * @brief General oriented Box (BOX macrobody)
 *
 * MCNP BOX macrobody: defined by corner vertex and 3 edge vectors.
 * The box spans from corner to corner + v1 + v2 + v3.
 */
typedef struct {
    double corner_x, corner_y, corner_z;    // Corner vertex
    double v1_x, v1_y, v1_z;                // First edge vector
    double v2_x, v2_y, v2_z;                // Second edge vector
    double v3_x, v3_y, v3_z;                // Third edge vector
} alea_box_general_data_t;

/**
 * @brief Sphere macrobody (SPH)
 *
 * Same as regular sphere but stored as MCNP macrobody format.
 */
typedef struct {
    double center_x, center_y, center_z;
    double radius;
} alea_sph_data_t;

/**
 * @brief Truncated Right Cone (TRC macrobody)
 *
 * MCNP TRC: defined by base center, height vector, base radius, top radius.
 * A cylinder if both radii are equal, a cone frustum otherwise.
 */
typedef struct {
    double base_x, base_y, base_z;          // Base center point
    double height_x, height_y, height_z;    // Height vector
    double base_radius;                      // Radius at base
    double top_radius;                       // Radius at top
} alea_trc_data_t;

/**
 * @brief Ellipsoid (ELL macrobody)
 *
 * MCNP ELL: defined by two foci (V1, V2) and the length of the major axis.
 * OR by center point and three semi-axes.
 */
typedef struct {
    double v1_x, v1_y, v1_z;    // First focus (or center if semi_a/b/c form)
    double v2_x, v2_y, v2_z;    // Second focus (or semi-axis lengths)
    double major_axis_len;       // Length of major axis (2a)
} alea_ell_data_t;

/**
 * @brief Right Elliptical Cylinder (REC macrobody)
 *
 * MCNP REC: base center, height vector, and two perpendicular semi-axis vectors.
 */
typedef struct {
    double base_x, base_y, base_z;          // Base center point
    double height_x, height_y, height_z;    // Height vector
    double axis1_x, axis1_y, axis1_z;       // First semi-axis vector
    double axis2_x, axis2_y, axis2_z;       // Second semi-axis vector (perpendicular)
} alea_rec_data_t;

/**
 * @brief Wedge (WED macrobody)
 *
 * MCNP WED: vertex and three edge vectors defining triangular prism.
 * The base is a right triangle, extruded along v3.
 */
typedef struct {
    double vertex_x, vertex_y, vertex_z;    // Vertex of right angle
    double v1_x, v1_y, v1_z;                // First leg of right triangle
    double v2_x, v2_y, v2_z;                // Second leg (perpendicular to v1)
    double v3_x, v3_y, v3_z;                // Height/extrusion vector
} alea_wed_data_t;

/**
 * @brief Right Hexagonal Prism (RHP macrobody)
 *
 * MCNP RHP/HEX: base center, height vector, and three vectors to hex vertices.
 * Vertices are at ±v1, ±v2, ±v3 from center, forming regular hexagon cross-section.
 */
typedef struct {
    double base_x, base_y, base_z;          // Base center
    double height_x, height_y, height_z;    // Height vector
    double r1_x, r1_y, r1_z;                // Vector to first vertex pair
    double r2_x, r2_y, r2_z;                // Vector to second vertex pair
    double r3_x, r3_y, r3_z;                // Vector to third vertex pair
} alea_rhp_data_t;

/**
 * @brief Arbitrary Polyhedron (ARB macrobody)
 *
 * MCNP ARB: up to 8 corner vertices and up to 6 four-sided faces.
 * Each face is defined by 4 vertex indices (1-8).
 */
typedef struct {
    double corners[8][3];                    // Up to 8 corner vertices (x,y,z)
    int faces[6][4];                         // Up to 6 faces, each with 4 vertex indices
    int num_corners;                         // Actual number of corners (4-8)
    int num_faces;                           // Actual number of faces (4-6)
} alea_arb_data_t;

/**
 * @brief Union of all primitive data types
 * 
 * Allows storing any primitive in a single memory location.
 * The actual type is determined by the accompanying alea_primitive_type_t.
 */
typedef union {
    alea_plane_data_t plane;
    alea_sphere_data_t sphere;
    alea_cylinder_x_data_t cyl_x;
    alea_cylinder_y_data_t cyl_y;
    alea_cylinder_z_data_t cyl_z;
    alea_cone_x_data_t cone_x;
    alea_cone_y_data_t cone_y;
    alea_cone_z_data_t cone_z;
    alea_box_data_t box;
    alea_quadric_data_t quadric;
    alea_torus_data_t torus;
    alea_rcc_data_t rcc;
    // Macrobodies
    alea_box_general_data_t box_general;
    alea_sph_data_t sph;
    alea_trc_data_t trc;
    alea_ell_data_t ell;
    alea_rec_data_t rec;
    alea_wed_data_t wed;
    alea_rhp_data_t rhp;
    alea_arb_data_t arb;
} alea_primitive_data_t;



// ============================================================================
// CSG NODE TYPES
// ============================================================================





// ============================================================================
// HELPER FUNCTIONS
// ============================================================================

/**
 * @brief Get minimum of two doubles
 */
static inline double alea_min(double a, double b) {
    return a < b ? a : b;
}

/**
 * @brief Get maximum of two doubles
 */
static inline double alea_max(double a, double b) {
    return a > b ? a : b;
}

/**
 * @brief Clamp value to range [min, max]
 */
static inline double alea_clamp(double value, double min, double max) {
    return value < min ? min : (value > max ? max : value);
}

#ifdef __cplusplus
}
#endif

#endif // ALEA_TYPES_H
