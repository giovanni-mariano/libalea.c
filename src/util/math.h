// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_MATH_H
#define ALEA_MATH_H

#include <math.h>

/* ===============================================
 * Common math constants
 * =============================================== */

#ifndef M_PI
#define M_PI 3.141592653589793238462643383279502884L
#endif

#ifndef M_SQRT3
#define M_SQRT3 1.7320508075688772
#endif

/* ===============================================
 * Common macros
 * =============================================== */

#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif

#ifndef MAX
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif

#ifndef CLAMP
#define CLAMP(x, lo, hi) ((x) < (lo) ? (lo) : ((x) > (hi) ? (hi) : (x)))
#endif

#define DEG_TO_RAD(x) ((x) * (M_PI / 180.0))
#define RAD_TO_DEG(x) ((x) * (180.0 / M_PI))

static inline double clamp01(double x) {
    return x < 0.0 ? 0.0 : (x > 1.0 ? 1.0 : x);
}

static inline float clamp01f(float x) {
    return x < 0.0f ? 0.0f : (x > 1.0f ? 1.0f : x);
}

/* ===============================================
 * Types
 * =============================================== */

typedef struct {
    double x, y, z;
} vec3;

/* Column-major 3x3 matrix */
typedef struct {
    double m[3][3]; /* m[col][row] */
} mat3;

/* Column-major 4x4 */
typedef struct {
    double m[4][4];
} mat4;



/*===============================================*/ 
/* Constructors */
static inline vec3 vec3_make(double x, double y, double z) {
    return (vec3){x, y, z};
}

/* Basic ops */
static inline vec3 vec3_add(vec3 a, vec3 b) {
    return (vec3){a.x + b.x, a.y + b.y, a.z + b.z};
}

static inline vec3 vec3_sub(vec3 a, vec3 b) {
    return (vec3){a.x - b.x, a.y - b.y, a.z - b.z};
}

static inline vec3 vec3_scale(vec3 v, double s) {
    return (vec3){v.x * s, v.y * s, v.z * s};
}

static inline double vec3_dot(vec3 a, vec3 b) {
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

static inline vec3 vec3_cross(vec3 a, vec3 b) {
    return (vec3){
        a.y*b.z - a.z*b.y,
        a.z*b.x - a.x*b.z,
        a.x*b.y - a.y*b.x
    };
}

static inline double vec3_length(vec3 v) {
    return sqrt(vec3_dot(v, v));
}

static inline vec3 vec3_normalize(vec3 v) {
    double len = vec3_length(v);
    return (len > 0.0) ? vec3_scale(v, 1.0 / len) : v;
}


static inline mat3 mat3_identity(void) {
    return (mat3){{
        {1,0,0},
        {0,1,0},
        {0,0,1}
    }};
}

static inline vec3 mat3_mul_vec3(mat3 A, vec3 v) {
    return (vec3){
        A.m[0][0]*v.x + A.m[1][0]*v.y + A.m[2][0]*v.z,
        A.m[0][1]*v.x + A.m[1][1]*v.y + A.m[2][1]*v.z,
        A.m[0][2]*v.x + A.m[1][2]*v.y + A.m[2][2]*v.z
    };
}

static inline mat3 mat3_mul(mat3 A, mat3 B) {
    mat3 R;
    for (int c = 0; c < 3; ++c)
        for (int r = 0; r < 3; ++r)
            R.m[c][r] =
                A.m[0][r]*B.m[c][0] +
                A.m[1][r]*B.m[c][1] +
                A.m[2][r]*B.m[c][2];
    return R;
}

static inline mat3 mat3_transpose(mat3 A) {
    mat3 T;
    for (int c = 0; c < 3; ++c)
        for (int r = 0; r < 3; ++r)
            T.m[c][r] = A.m[r][c];
    return T;
}

/* Determinant of a 3x3 column-major matrix.
 * Row i = (A.m[0][i], A.m[1][i], A.m[2][i]) */
static inline double mat3_det(mat3 A) {
    return A.m[0][0] * (A.m[1][1]*A.m[2][2] - A.m[2][1]*A.m[1][2])
         - A.m[1][0] * (A.m[0][1]*A.m[2][2] - A.m[2][1]*A.m[0][2])
         + A.m[2][0] * (A.m[0][1]*A.m[1][2] - A.m[1][1]*A.m[0][2]);
}

/* Solve Ax = b via Cramer's rule.  Returns 0 on success, -1 if singular.
 * A is column-major: row i of A is (A.m[0][i], A.m[1][i], A.m[2][i]). */
static inline int mat3_solve_cramer(mat3 A, vec3 b, vec3 *out) {
    double det = mat3_det(A);
    if (fabs(det) < 1e-14) return -1;
    double inv = 1.0 / det;

    /* Replace column 0 with b */
    out->x = (b.x * (A.m[1][1]*A.m[2][2] - A.m[2][1]*A.m[1][2])
            - A.m[1][0] * (b.y*A.m[2][2] - A.m[2][1]*b.z)
            + A.m[2][0] * (b.y*A.m[1][2] - A.m[1][1]*b.z)) * inv;

    /* Replace column 1 with b */
    out->y = (A.m[0][0] * (b.y*A.m[2][2] - A.m[2][1]*b.z)
            - b.x * (A.m[0][1]*A.m[2][2] - A.m[2][1]*A.m[0][2])
            + A.m[2][0] * (A.m[0][1]*b.z - b.y*A.m[0][2])) * inv;

    /* Replace column 2 with b */
    out->z = (A.m[0][0] * (A.m[1][1]*b.z - b.y*A.m[1][2])
            - A.m[1][0] * (A.m[0][1]*b.z - b.y*A.m[0][2])
            + b.x * (A.m[0][1]*A.m[1][2] - A.m[1][1]*A.m[0][2])) * inv;

    return 0;
}


static inline mat4 mat4_identity(void) {
    mat4 M = {0};
    for (int i = 0; i < 4; ++i) M.m[i][i] = 1.0;
    return M;
}

static inline vec3 mat4_mul_point(mat4 M, vec3 p) {
    double x =
        M.m[0][0]*p.x + M.m[1][0]*p.y + M.m[2][0]*p.z + M.m[3][0];
    double y =
        M.m[0][1]*p.x + M.m[1][1]*p.y + M.m[2][1]*p.z + M.m[3][1];
    double z =
        M.m[0][2]*p.x + M.m[1][2]*p.y + M.m[2][2]*p.z + M.m[3][2];
    return (vec3){x, y, z};
}

static inline vec3 mat4_mul_dir(mat4 M, vec3 v) {
    return (vec3){
        M.m[0][0]*v.x + M.m[1][0]*v.y + M.m[2][0]*v.z,
        M.m[0][1]*v.x + M.m[1][1]*v.y + M.m[2][1]*v.z,
        M.m[0][2]*v.x + M.m[1][2]*v.y + M.m[2][2]*v.z
    };
}


/* ===============================================
 * Array-based vector helpers (for double[3] fields)
 * =============================================== */

static inline double v3arr_dot(const double* a, const double* b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static inline void v3arr_cross(const double* a, const double* b, double* out) {
    out[0] = a[1]*b[2] - a[2]*b[1];
    out[1] = a[2]*b[0] - a[0]*b[2];
    out[2] = a[0]*b[1] - a[1]*b[0];
}

static inline double v3arr_length(const double* v) {
    return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

static inline void v3arr_normalize(double* v) {
    double len = v3arr_length(v);
    if (len > 1e-15) {
        v[0] /= len; v[1] /= len; v[2] /= len;
    }
}

static inline void v3arr_sub(const double* a, const double* b, double* out) {
    out[0] = a[0]-b[0]; out[1] = a[1]-b[1]; out[2] = a[2]-b[2];
}

static inline void v3arr_scale(const double* v, double s, double* out) {
    out[0] = v[0]*s; out[1] = v[1]*s; out[2] = v[2]*s;
}

static inline void v3arr_add(const double* a, const double* b, double* out) {
    out[0] = a[0]+b[0]; out[1] = a[1]+b[1]; out[2] = a[2]+b[2];
}

#endif /* ALEA_MATH_H */
