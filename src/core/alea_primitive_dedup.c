// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_system.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

// ============================================================================
// HASH COMPUTATION
// ============================================================================

// FNV-1a hash function
static uint64_t fnv1a_hash(const void* data, size_t len) {
    const uint8_t* bytes = (const uint8_t*)data;
    uint64_t hash = 14695981039346656037ULL;
    
    for (size_t i = 0; i < len; i++) {
        hash ^= bytes[i];
        hash *= 1099511628211ULL;
    }
    
    return hash;
}

// Hash a double using tolerance quantization
static uint64_t hash_double(double value, const alea_config_t* tol) {
    // Quantize to tolerance grid
    if (fabs(value) < tol->zero_threshold) {
        value = 0.0;
    } else {
        // Round to tolerance
        double quant = tol->abs_tol;
        value = round(value / quant) * quant;
    }
    
    uint64_t bits;
    memcpy(&bits, &value, sizeof(double));
    return bits;
}

uint64_t alea_compute_primitive_hash(alea_primitive_type_t type, const alea_primitive_data_t* data,
                                     const alea_config_t* config) {
    const alea_config_t* tol = config ? config : &ALEA_CONFIG_DEFAULT;

    uint64_t hash = fnv1a_hash(&type, sizeof(type));
    
    switch (type) {
        case ALEA_PRIMITIVE_PLANE: {
            const alea_plane_data_t* p = &data->plane;
            hash ^= hash_double(p->a, tol);
            hash ^= hash_double(p->b, tol);
            hash ^= hash_double(p->c, tol);
            hash ^= hash_double(p->d, tol);
            break;
        }
        
        case ALEA_PRIMITIVE_SPHERE: {
            const alea_sphere_data_t* s = &data->sphere;
            hash ^= hash_double(s->center_x, tol);
            hash ^= hash_double(s->center_y, tol);
            hash ^= hash_double(s->center_z, tol);
            hash ^= hash_double(s->radius, tol);
            break;
        }
        
        case ALEA_PRIMITIVE_CYLINDER_Z: {
            const alea_cylinder_z_data_t* c = &data->cyl_z;
            hash ^= hash_double(c->center_x, tol);
            hash ^= hash_double(c->center_y, tol);
            hash ^= hash_double(c->radius, tol);
            break;
        }

        case ALEA_PRIMITIVE_CYLINDER_X: {
            const alea_cylinder_x_data_t* c = &data->cyl_x;
            hash ^= hash_double(c->center_y, tol);
            hash ^= hash_double(c->center_z, tol);
            hash ^= hash_double(c->radius, tol);
            break;
        }

        case ALEA_PRIMITIVE_CYLINDER_Y: {
            const alea_cylinder_y_data_t* c = &data->cyl_y;
            hash ^= hash_double(c->center_x, tol);
            hash ^= hash_double(c->center_z, tol);
            hash ^= hash_double(c->radius, tol);
            break;
        }
        
        case ALEA_PRIMITIVE_CONE_X: {
            const alea_cone_x_data_t* c = &data->cone_x;
            hash ^= hash_double(c->apex_x, tol);
            hash ^= hash_double(c->apex_y, tol);
            hash ^= hash_double(c->apex_z, tol);
            hash ^= hash_double(c->tan_angle_sq, tol);
            hash ^= fnv1a_hash(&c->sheet_selection, sizeof(c->sheet_selection));
            break;
        }

        case ALEA_PRIMITIVE_CONE_Y: {
            const alea_cone_y_data_t* c = &data->cone_y;
            hash ^= hash_double(c->apex_x, tol);
            hash ^= hash_double(c->apex_y, tol);
            hash ^= hash_double(c->apex_z, tol);
            hash ^= hash_double(c->tan_angle_sq, tol);
            hash ^= fnv1a_hash(&c->sheet_selection, sizeof(c->sheet_selection));
            break;
        }

        case ALEA_PRIMITIVE_CONE_Z: {
            const alea_cone_z_data_t* c = &data->cone_z;
            hash ^= hash_double(c->apex_x, tol);
            hash ^= hash_double(c->apex_y, tol);
            hash ^= hash_double(c->apex_z, tol);
            hash ^= hash_double(c->tan_angle_sq, tol);
            hash ^= fnv1a_hash(&c->sheet_selection, sizeof(c->sheet_selection));
            break;
        }
        
        case ALEA_PRIMITIVE_TORUS_X:
        case ALEA_PRIMITIVE_TORUS_Y:
        case ALEA_PRIMITIVE_TORUS_Z: {
            const alea_torus_data_t* t = &data->torus;
            hash ^= fnv1a_hash(&t->axis, sizeof(t->axis));  // Hash the enum
            hash ^= hash_double(t->center_x, tol);
            hash ^= hash_double(t->center_y, tol);
            hash ^= hash_double(t->center_z, tol);
            hash ^= hash_double(t->major_radius, tol);
            hash ^= hash_double(t->minor_radius, tol);
            hash ^= hash_double(t->axial_semiwidth_B, tol);
            break;
        }
        
        case ALEA_PRIMITIVE_RPP: {
            const alea_box_data_t* b = &data->box;
            hash ^= hash_double(b->min_x, tol);
            hash ^= hash_double(b->max_x, tol);
            hash ^= hash_double(b->min_y, tol);
            hash ^= hash_double(b->max_y, tol);
            hash ^= hash_double(b->min_z, tol);
            hash ^= hash_double(b->max_z, tol);
            break;
        }
        
        case ALEA_PRIMITIVE_QUADRIC: {
            const alea_quadric_data_t* q = &data->quadric;
            for (int i = 0; i < 10; i++) {
                hash ^= hash_double(q->coeffs[i], tol);
            }
            break;
        }

        case ALEA_PRIMITIVE_RCC: {
            const alea_rcc_data_t* rcc = &data->rcc;
            hash ^= hash_double(rcc->base_x, tol);
            hash ^= hash_double(rcc->base_y, tol);
            hash ^= hash_double(rcc->base_z, tol);
            hash ^= hash_double(rcc->height_x, tol);
            hash ^= hash_double(rcc->height_y, tol);
            hash ^= hash_double(rcc->height_z, tol);
            hash ^= hash_double(rcc->radius, tol);
            break;
        }

        case ALEA_PRIMITIVE_BOX: {
            const alea_box_general_data_t* bg = &data->box_general;
            hash ^= hash_double(bg->corner_x, tol);
            hash ^= hash_double(bg->corner_y, tol);
            hash ^= hash_double(bg->corner_z, tol);
            hash ^= hash_double(bg->v1_x, tol);
            hash ^= hash_double(bg->v1_y, tol);
            hash ^= hash_double(bg->v1_z, tol);
            hash ^= hash_double(bg->v2_x, tol);
            hash ^= hash_double(bg->v2_y, tol);
            hash ^= hash_double(bg->v2_z, tol);
            hash ^= hash_double(bg->v3_x, tol);
            hash ^= hash_double(bg->v3_y, tol);
            hash ^= hash_double(bg->v3_z, tol);
            break;
        }

        case ALEA_PRIMITIVE_SPH: {
            const alea_sph_data_t* sph = &data->sph;
            hash ^= hash_double(sph->center_x, tol);
            hash ^= hash_double(sph->center_y, tol);
            hash ^= hash_double(sph->center_z, tol);
            hash ^= hash_double(sph->radius, tol);
            break;
        }

        case ALEA_PRIMITIVE_TRC: {
            const alea_trc_data_t* trc = &data->trc;
            hash ^= hash_double(trc->base_x, tol);
            hash ^= hash_double(trc->base_y, tol);
            hash ^= hash_double(trc->base_z, tol);
            hash ^= hash_double(trc->height_x, tol);
            hash ^= hash_double(trc->height_y, tol);
            hash ^= hash_double(trc->height_z, tol);
            hash ^= hash_double(trc->base_radius, tol);
            hash ^= hash_double(trc->top_radius, tol);
            break;
        }

        case ALEA_PRIMITIVE_ELL: {
            const alea_ell_data_t* ell = &data->ell;
            hash ^= hash_double(ell->v1_x, tol);
            hash ^= hash_double(ell->v1_y, tol);
            hash ^= hash_double(ell->v1_z, tol);
            hash ^= hash_double(ell->v2_x, tol);
            hash ^= hash_double(ell->v2_y, tol);
            hash ^= hash_double(ell->v2_z, tol);
            hash ^= hash_double(ell->major_axis_len, tol);
            break;
        }

        case ALEA_PRIMITIVE_REC: {
            const alea_rec_data_t* rec = &data->rec;
            hash ^= hash_double(rec->base_x, tol);
            hash ^= hash_double(rec->base_y, tol);
            hash ^= hash_double(rec->base_z, tol);
            hash ^= hash_double(rec->height_x, tol);
            hash ^= hash_double(rec->height_y, tol);
            hash ^= hash_double(rec->height_z, tol);
            hash ^= hash_double(rec->axis1_x, tol);
            hash ^= hash_double(rec->axis1_y, tol);
            hash ^= hash_double(rec->axis1_z, tol);
            hash ^= hash_double(rec->axis2_x, tol);
            hash ^= hash_double(rec->axis2_y, tol);
            hash ^= hash_double(rec->axis2_z, tol);
            break;
        }

        case ALEA_PRIMITIVE_WED: {
            const alea_wed_data_t* wed = &data->wed;
            hash ^= hash_double(wed->vertex_x, tol);
            hash ^= hash_double(wed->vertex_y, tol);
            hash ^= hash_double(wed->vertex_z, tol);
            hash ^= hash_double(wed->v1_x, tol);
            hash ^= hash_double(wed->v1_y, tol);
            hash ^= hash_double(wed->v1_z, tol);
            hash ^= hash_double(wed->v2_x, tol);
            hash ^= hash_double(wed->v2_y, tol);
            hash ^= hash_double(wed->v2_z, tol);
            hash ^= hash_double(wed->v3_x, tol);
            hash ^= hash_double(wed->v3_y, tol);
            hash ^= hash_double(wed->v3_z, tol);
            break;
        }

        case ALEA_PRIMITIVE_RHP: {
            const alea_rhp_data_t* rhp = &data->rhp;
            hash ^= hash_double(rhp->base_x, tol);
            hash ^= hash_double(rhp->base_y, tol);
            hash ^= hash_double(rhp->base_z, tol);
            hash ^= hash_double(rhp->height_x, tol);
            hash ^= hash_double(rhp->height_y, tol);
            hash ^= hash_double(rhp->height_z, tol);
            hash ^= hash_double(rhp->r1_x, tol);
            hash ^= hash_double(rhp->r1_y, tol);
            hash ^= hash_double(rhp->r1_z, tol);
            hash ^= hash_double(rhp->r2_x, tol);
            hash ^= hash_double(rhp->r2_y, tol);
            hash ^= hash_double(rhp->r2_z, tol);
            hash ^= hash_double(rhp->r3_x, tol);
            hash ^= hash_double(rhp->r3_y, tol);
            hash ^= hash_double(rhp->r3_z, tol);
            break;
        }

        case ALEA_PRIMITIVE_ARB: {
            const alea_arb_data_t* arb = &data->arb;
            // Hash corner count and face count
            hash ^= fnv1a_hash(&arb->num_corners, sizeof(arb->num_corners));
            hash ^= fnv1a_hash(&arb->num_faces, sizeof(arb->num_faces));
            // Hash corners
            for (int i = 0; i < arb->num_corners; i++) {
                hash ^= hash_double(arb->corners[i][0], tol);
                hash ^= hash_double(arb->corners[i][1], tol);
                hash ^= hash_double(arb->corners[i][2], tol);
            }
            // Hash faces (integer indices)
            for (int i = 0; i < arb->num_faces; i++) {
                hash ^= fnv1a_hash(arb->faces[i], sizeof(arb->faces[i]));
            }
            break;
        }
    }

    return hash;
}

// ============================================================================
// HASH TABLE IMPLEMENTATION
// ============================================================================

primitive_hash_table_t* primitive_hash_table_create(void) {
    primitive_hash_table_t* table = malloc(sizeof(primitive_hash_table_t));
    if (!table) return NULL;
    
    table->bucket_count = PRIMITIVE_HASH_INITIAL_SIZE;
    table->entry_count = 0;
    table->buckets = calloc(table->bucket_count, sizeof(primitive_hash_entry_t*));
    
    if (!table->buckets) {
        free(table);
        return NULL;
    }
    
    return table;
}

void primitive_hash_table_destroy(primitive_hash_table_t* table) {
    if (!table) return;
    
    for (size_t i = 0; i < table->bucket_count; i++) {
        primitive_hash_entry_t* entry = table->buckets[i];
        while (entry) {
            primitive_hash_entry_t* next = entry->next;
            free(entry);
            entry = next;
        }
    }
    
    free(table->buckets);
    free(table);
}

static void primitive_hash_table_resize(primitive_hash_table_t* table) {
    size_t new_bucket_count = table->bucket_count * 2;
    primitive_hash_entry_t** new_buckets = calloc(new_bucket_count, sizeof(primitive_hash_entry_t*));
    
    if (!new_buckets) return; // Continue with current size
    
    // Rehash all entries
    for (size_t i = 0; i < table->bucket_count; i++) {
        primitive_hash_entry_t* entry = table->buckets[i];
        while (entry) {
            primitive_hash_entry_t* next = entry->next;
            
            size_t new_bucket = entry->hash % new_bucket_count;
            entry->next = new_buckets[new_bucket];
            new_buckets[new_bucket] = entry;
            
            entry = next;
        }
    }
    
    free(table->buckets);
    table->buckets = new_buckets;
    table->bucket_count = new_bucket_count;
}

uint32_t primitive_hash_table_find(primitive_hash_table_t* table, alea_system_t* sys,
                                   alea_primitive_type_t type, alea_primitive_data_t* data,
                                   uint64_t hash, int8_t* match_inverted) {
    size_t bucket = hash % table->bucket_count;
    primitive_hash_entry_t* entry = table->buckets[bucket];

    while (entry) {
        if (entry->hash == hash) {
            // Hash matches, check full equality
            const alea_primitive_entry_t* prim = &sys->primitives.data[entry->primitive_id];
            if (alea_primitives_equal(type, data, prim->type, &prim->data,
                                     &sys->config, match_inverted)) {
                return entry->primitive_id;
            }
        }
        entry = entry->next;
    }

    return UINT32_MAX; // Not found
}

void primitive_hash_table_insert(primitive_hash_table_t* table, uint32_t primitive_id, uint64_t hash) {
    // Resize if load factor > 0.75
    if (table->entry_count > table->bucket_count * 3 / 4) {
        primitive_hash_table_resize(table);
    }
    
    size_t bucket = hash % table->bucket_count;
    
    primitive_hash_entry_t* entry = malloc(sizeof(primitive_hash_entry_t));
    if (!entry) return;
    
    entry->primitive_id = primitive_id;
    entry->hash = hash;
    entry->next = table->buckets[bucket];
    
    table->buckets[bucket] = entry;
    table->entry_count++;
}