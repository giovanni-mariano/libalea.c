// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

#ifndef ALEA_COMPUTE_H
#define ALEA_COMPUTE_H

#include "alea_types.h"
#include "alea_render.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Experimental, optional rendering and batch compute module (libalea_compute).
 * Native backends own their device layouts and execution; no GPU portability layer.
 * Rendering and batch point classification support transformed fill hierarchies.
 */
typedef struct alea_compute alea_compute_t;

typedef enum {
    ALEA_COMPUTE_CPU = 0,
    ALEA_COMPUTE_METAL = 1
} alea_compute_backend_t;

typedef enum {
    ALEA_COMPUTE_FP64 = 0, /**< Existing CPU double-precision evaluation. */
    ALEA_COMPUTE_FP32 = 1  /**< Approximate GPU evaluation; opt in explicitly. */
} alea_compute_precision_t;

typedef struct {
    int32_t cell_id;      /**< -1 unless exactly one cell owns the point. */
    int32_t material_id;  /**< -1 unless exactly one cell owns the point. */
    uint32_t owner_count; /**< Terminal occurrences: 0 = gap, 1 = unique, >1 = overlap. */
} alea_compute_result_t;

/** Return whether a backend was built and a suitable device is available.
 * Does not guarantee that a particular geometry or allocation is supported.
 */
int alea_compute_available(alea_compute_backend_t backend);

/** Prepare an owned snapshot; the source may subsequently change or be freed.
 * Recreate the context to pick up geometry changes. No implicit backend or
 * precision fallback. CPU requires FP64; Metal requires FP32.
 * Metal supports planes, spheres, axis-aligned cylinders/cones/tori (including
 * elliptical tubes), general quadrics, RPP boxes and RCC cylinders. CSG is
 * compiled into short-circuit bytecode, with source depth <= 4096. Rendering
 * supports transformed fills (up to 32 universe levels) and rectangular lattices.
 * Hexagonal lattices and lattice-local fill transforms are rejected.
 * FP32 can change ownership near boundaries or erase small features. It is
 * intended for approximate visualization, not transport or verified queries.
 * Returns NULL on failure; details are available through alea_error().
 */
alea_compute_t* alea_compute_create(const alea_system_t* sys,
                                   alea_compute_backend_t backend,
                                   alea_compute_precision_t precision);
void alea_compute_destroy(alea_compute_t* compute);

/** Synchronously classify count interleaved XYZ triples.
 * Counts terminal owning occurrences from universe 0, including terminal void,
 * under all containing fill ancestors. Repeated placements count separately.
 * Containers are not owners; a gap in their contents remains a gap. IDs are -1
 * unless ownership is unique. Supports transformed fills through 32 universe
 * levels and finite/repeating rectangular lattices, including outer fills.
 * Empty/nonpositive element fills are gaps; lattice containers are not owners.
 * Element selection uses floor on tiled axes. Finite singleton Y/Z windows
 * include the upper face, matching core point lookup; untiled Z remains
 * unshifted in center coordinates. Hexagonal lattices and lattice-local fill
 * transforms are rejected. Metal also rejects element indices whose magnitude
 * exceeds 16777215 rather than aliasing adjacent FP32 indices.
 * This samples points, not complete intervals or thin-feature visibility.
 * Metal uses approximate FP32 evaluation, including transformed coordinates.
 * All inputs must be finite; Metal inputs must also fit in float. Input and
 * output storage must not overlap. count == 0 permits NULL arrays.
 * One context must not be used concurrently. Different contexts are independent.
 * Metal buffers are reused between calls; hierarchy buffers upload lazily once.
 * Certified spatial CSG programs reduce expression work inside their domains;
 * original programs remain the fallback. All cells in visited universes remain
 * eligible, including void. ALEA_METAL_POINT_SPATIAL=0 at snapshot creation
 * disables point-program selection for comparison (default: enabled).
 * Hierarchical Metal batches submit up to 4096 points at a time, checking
 * cancellation between submissions; an individual submission has no time bound.
 * Returns 0 on success, -1 on error (output is unspecified on error).
 */
int alea_compute_classify(alea_compute_t* compute, const double* xyz,
                          size_t count, alea_compute_result_t* results);

/** Render the prepared snapshot into an existing framebuffer (linear RGB plus
 * available auxiliary channels). Reuses device allocations across frames.
 * CPU uses render_scene(). Metal is an approximate FP32 solid renderer with
 * shared universe BVHs, ancestor clipping, one user clipping plane,
 * material/cell colors, and up to 4x4 supersampling. Shadows,
 * X-ray rendering, density/universe colors and selection filters are rejected.
 * Edge darkening and depth/ID visualization remain caller-side postprocessing.
 * Config and derived camera must match framebuffer dimensions. Returns 0/-1.
 */
int alea_compute_render(alea_compute_t* compute, const render_config_t* config,
                        const render_camera_t* camera, render_framebuffer_t* framebuffer);

/** Progress callback for tiled Metal ray rendering. Called synchronously on the
 * calling thread before work (completed=0) and after each GPU submission. The
 * completed count may be unchanged while rays are suspended inside CSG. Only
 * completed pixels have been overwritten; seed the framebuffer for previews.
 * Return nonzero to cancel this call, without setting the global interrupt.
 * Do not reenter or destroy the context from the callback. */
typedef int (*alea_compute_render_progress_t)(void* user,
    const render_framebuffer_t* framebuffer, size_t completed, size_t total);
/** Same final pixels as alea_compute_render; visits 64x64 screen tiles, starting
 * at the center. Up to 4096 pixels advance in parallel; unfinished traversal and
 * CSG stacks survive between submissions. Work per ray adapts toward 20 ms
 * per submission (not a latency guarantee). Continuation scratch is at most
 * 64 MiB, independent of image size; image buffers are additional.
 * Returns 0 complete, 1 callback cancellation, -1 error. On cancellation,
 * completed pixels remain available and the context can be reused immediately.
 * Requires a Metal ray context and non-NULL callback; CPU contexts reject
 * this operation. Callbacks can publish copies to a UI and check cancellation. */
int alea_compute_render_progressive(alea_compute_t* compute,
    const render_config_t* config, const render_camera_t* camera,
    render_framebuffer_t* framebuffer, alea_compute_render_progress_t progress,
    void* user);

#ifdef __cplusplus
}
#endif
#endif
