// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
//
// SPDX-License-Identifier: MPL-2.0

#include "alea_xml.h"
#include "alea_xml/region.h"
#include "alea_xml/xml_dom.h"
#include "alea_xml/xml_writer.h"
#include "core/alea_system.h"
#include "util/compat.h"
#include "util/str_builder.h"

#include <errno.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define ALEA_XML_VERSION "1"

static int fail(const char* message) {
    alea_set_error_detail(ALEA_ERR_PARSE_ERROR, "ALEA XML: %s", message);
    return -1;
}

static int fail_element(const alea_xml_dom_element_t* e, const char* message) {
    alea_set_error_detail(ALEA_ERR_PARSE_ERROR, "ALEA XML <%s>: %s",
                          e && e->tag_name ? e->tag_name : "?", message);
    return -1;
}

static int parse_int_strict(const char* text, int* out) {
    char* end = NULL;
    long value;
    if (!text || !*text || !out) return -1;
    errno = 0;
    value = strtol(text, &end, 10);
    while (end && isspace((unsigned char)*end)) end++;
    if (errno || !end || *end || value < INT_MIN || value > INT_MAX) return -1;
    *out = (int)value;
    return 0;
}

static int parse_double_strict(const char* text, double* out) {
    char* end = NULL;
    double value;
    if (!text || !*text || !out) return -1;
    errno = 0;
    value = strtod(text, &end);
    while (end && isspace((unsigned char)*end)) end++;
    if (errno || !end || *end || !isfinite(value)) return -1;
    *out = value;
    return 0;
}

static int parse_bool(const char* text, int default_value, int* out) {
    if (!out) return -1;
    if (!text) { *out = default_value; return 0; }
    if (!strcmp(text, "1") || !strcmp(text, "true")) { *out = 1; return 0; }
    if (!strcmp(text, "0") || !strcmp(text, "false")) { *out = 0; return 0; }
    return -1;
}

static int parse_values(const char* text, double* values, size_t count) {
    const char* p = text;
    if (!p) return -1;
    for (size_t i = 0; i < count; i++) {
        char* end = NULL;
        while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') p++;
        if (!*p) return -1;
        errno = 0;
        values[i] = strtod(p, &end);
        if (errno || end == p || !isfinite(values[i])) return -1;
        p = end;
    }
    while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') p++;
    return *p ? -1 : 0;
}

static int parse_int_values(const char* text, int* values, size_t count) {
    const char* p = text;
    if (!p) return -1;
    for (size_t i = 0; i < count; i++) {
        char* end = NULL;
        long v;
        while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') p++;
        errno = 0;
        v = strtol(p, &end, 10);
        if (errno || end == p || v < INT_MIN || v > INT_MAX) return -1;
        values[i] = (int)v;
        p = end;
    }
    while (*p == ' ' || *p == '\t' || *p == '\r' || *p == '\n') p++;
    return *p ? -1 : 0;
}

static char* copy_string(const char* value) {
    return value && *value ? alea_strdup(value) : NULL;
}

static char* library_internal(const char* value) {
    size_t n;
    char* result;
    if (!value || !*value) return NULL;
    if (value[0] == '.') return alea_strdup(value);
    n = strlen(value);
    result = malloc(n + 2);
    if (!result) return NULL;
    result[0] = '.';
    memcpy(result + 1, value, n + 1);
    return result;
}

static int attr_allowed(const alea_xml_dom_element_t* e,
                        const char* const* names, size_t count) {
    for (size_t i = 0; i < e->attr_count; i++) {
        size_t j;
        for (j = 0; j < count; j++)
            if (!strcmp(e->attrs[i].name, names[j])) break;
        if (j == count) {
            alea_set_error_detail(ALEA_ERR_PARSE_ERROR,
                "ALEA XML <%s>: unknown attribute '%s'",
                e->tag_name, e->attrs[i].name);
            return -1;
        }
    }
    if (e->text_content && *e->text_content)
        return fail_element(e, "text content is not allowed");
    return 0;
}

#define VALIDATE_ATTRS(e, ...) do { \
    const char* const names[] = {__VA_ARGS__}; \
    if (attr_allowed((e), names, sizeof(names) / sizeof(names[0]))) return -1; \
} while (0)

static int validate_structure(const alea_xml_dom_element_t* root) {
    VALIDATE_ATTRS(root, "version", "length_units", "name", "title", "comments");
    int materials = 0, transforms = 0, surfaces = 0, cells = 0;
    for (size_t i = 0; i < root->child_count; i++) {
        alea_xml_dom_element_t* section = root->children[i];
        if (!strcmp(section->tag_name, "materials")) materials++;
        else if (!strcmp(section->tag_name, "transforms")) transforms++;
        else if (!strcmp(section->tag_name, "surfaces")) surfaces++;
        else if (!strcmp(section->tag_name, "cells")) cells++;
        else return fail_element(section, "unknown top-level section");
        if (section->attr_count || (section->text_content && *section->text_content))
            return fail_element(section, "section attributes and text are not allowed");

        for (size_t j = 0; j < section->child_count; j++) {
            alea_xml_dom_element_t* e = section->children[j];
            if (!strcmp(section->tag_name, "materials")) {
                if (!strcmp(e->tag_name, "material")) {
                    VALIDATE_ATTRS(e, "id", "name", "comments", "fraction_basis", "density", "density_units");
                    for (size_t k = 0; k < e->child_count; k++) {
                        alea_xml_dom_element_t* c = e->children[k];
                        if (!strcmp(c->tag_name, "nuclide")) VALIDATE_ATTRS(c, "zaid", "fraction", "library");
                        else if (!strcmp(c->tag_name, "element")) VALIDATE_ATTRS(c, "z", "fraction", "library");
                        else if (!strcmp(c->tag_name, "thermal")) VALIDATE_ATTRS(c, "identifier", "zaid_match");
                        else return fail_element(c, "unknown material child");
                        if (c->child_count) return fail_element(c, "child elements are not allowed");
                    }
                } else if (!strcmp(e->tag_name, "mixture")) {
                    VALIDATE_ATTRS(e, "id", "fraction_basis", "mc_material_id", "name", "comments");
                    for (size_t k = 0; k < e->child_count; k++) {
                        alea_xml_dom_element_t* c = e->children[k];
                        if (strcmp(c->tag_name, "component")) return fail_element(c, "unknown mixture child");
                        VALIDATE_ATTRS(c, "material", "fraction");
                        if (c->child_count) return fail_element(c, "child elements are not allowed");
                    }
                } else return fail_element(e, "unknown materials child");
            } else if (!strcmp(section->tag_name, "transforms")) {
                if (strcmp(e->tag_name, "transform")) return fail_element(e, "unknown transforms child");
                VALIDATE_ATTRS(e, "id", "values", "degrees", "inline");
                if (e->child_count) return fail_element(e, "child elements are not allowed");
            } else if (!strcmp(section->tag_name, "surfaces")) {
                if (strcmp(e->tag_name, "surface")) return fail_element(e, "unknown surfaces child");
                VALIDATE_ATTRS(e, "id", "type", "coeffs", "boundary", "periodic_surface", "transform", "transform_applied");
                if (e->child_count) return fail_element(e, "child elements are not allowed");
            } else {
                if (strcmp(e->tag_name, "cell")) return fail_element(e, "unknown cells child");
                VALIDATE_ATTRS(e, "id", "name", "comments", "inline_comment", "universe", "material", "material_kind", "density", "density_units", "region", "bbox", "temperature", "fill", "fill_transform", "lattice_type", "lattice_dims", "lattice_fill", "lattice_pitch", "lattice_lower_left", "lattice_outer", "lattice_repeating", "lattice_zero_coords");
                for (size_t k = 0; k < e->child_count; k++) {
                    alea_xml_dom_element_t* c = e->children[k];
                    if (!strcmp(c->tag_name, "importance")) VALIDATE_ATTRS(c, "particle", "value");
                    else if (!strcmp(c->tag_name, "parameter")) VALIDATE_ATTRS(c, "name", "value");
                    else return fail_element(c, "unknown cell child");
                    if (c->child_count) return fail_element(c, "child elements are not allowed");
                }
            }
        }
    }
    if (materials > 1 || transforms > 1 || surfaces > 1 || cells != 1)
        return fail("sections must be unique and exactly one <cells> is required");
    return 0;
}

static int surface_type(const char* name, alea_primitive_type_t* type,
                        size_t* coeff_count) {
    static const struct { const char* name; alea_primitive_type_t type; size_t n; } map[] = {
        {"plane", ALEA_PRIMITIVE_PLANE, 4}, {"sphere", ALEA_PRIMITIVE_SPHERE, 4},
        {"cylinder-x", ALEA_PRIMITIVE_CYLINDER_X, 3},
        {"cylinder-y", ALEA_PRIMITIVE_CYLINDER_Y, 3},
        {"cylinder-z", ALEA_PRIMITIVE_CYLINDER_Z, 3},
        {"cone-x", ALEA_PRIMITIVE_CONE_X, 5}, {"cone-y", ALEA_PRIMITIVE_CONE_Y, 5},
        {"cone-z", ALEA_PRIMITIVE_CONE_Z, 5}, {"rpp", ALEA_PRIMITIVE_RPP, 6},
        {"quadric", ALEA_PRIMITIVE_QUADRIC, 10},
        {"torus-x", ALEA_PRIMITIVE_TORUS_X, 6},
        {"torus-y", ALEA_PRIMITIVE_TORUS_Y, 6},
        {"torus-z", ALEA_PRIMITIVE_TORUS_Z, 6}, {"rcc", ALEA_PRIMITIVE_RCC, 7},
        {"box", ALEA_PRIMITIVE_BOX, 12}, {"sph", ALEA_PRIMITIVE_SPH, 4},
        {"trc", ALEA_PRIMITIVE_TRC, 8}, {"ell", ALEA_PRIMITIVE_ELL, 7},
        {"rec", ALEA_PRIMITIVE_REC, 12}, {"wed", ALEA_PRIMITIVE_WED, 12},
        {"rhp", ALEA_PRIMITIVE_RHP, 15}, {"arb", ALEA_PRIMITIVE_ARB, 50}
    };
    for (size_t i = 0; i < sizeof(map) / sizeof(map[0]); i++) {
        if (name && !strcmp(name, map[i].name)) {
            *type = map[i].type; *coeff_count = map[i].n; return 0;
        }
    }
    return -1;
}

static const char* surface_type_name(alea_primitive_type_t type) {
    switch (type) {
        case ALEA_PRIMITIVE_PLANE: return "plane";
        case ALEA_PRIMITIVE_SPHERE: return "sphere";
        case ALEA_PRIMITIVE_CYLINDER_X: return "cylinder-x";
        case ALEA_PRIMITIVE_CYLINDER_Y: return "cylinder-y";
        case ALEA_PRIMITIVE_CYLINDER_Z: return "cylinder-z";
        case ALEA_PRIMITIVE_CONE_X: return "cone-x";
        case ALEA_PRIMITIVE_CONE_Y: return "cone-y";
        case ALEA_PRIMITIVE_CONE_Z: return "cone-z";
        case ALEA_PRIMITIVE_RPP: return "rpp";
        case ALEA_PRIMITIVE_QUADRIC: return "quadric";
        case ALEA_PRIMITIVE_TORUS_X: return "torus-x";
        case ALEA_PRIMITIVE_TORUS_Y: return "torus-y";
        case ALEA_PRIMITIVE_TORUS_Z: return "torus-z";
        case ALEA_PRIMITIVE_RCC: return "rcc";
        case ALEA_PRIMITIVE_BOX: return "box";
        case ALEA_PRIMITIVE_SPH: return "sph";
        case ALEA_PRIMITIVE_TRC: return "trc";
        case ALEA_PRIMITIVE_ELL: return "ell";
        case ALEA_PRIMITIVE_REC: return "rec";
        case ALEA_PRIMITIVE_WED: return "wed";
        case ALEA_PRIMITIVE_RHP: return "rhp";
        case ALEA_PRIMITIVE_ARB: return "arb";
        default: return NULL;
    }
}

static int data_from_values(alea_primitive_type_t type, const double* v,
                            alea_primitive_data_t* d) {
    memset(d, 0, sizeof(*d));
    switch (type) {
        case ALEA_PRIMITIVE_PLANE: memcpy(&d->plane, v, 4 * sizeof(double)); break;
        case ALEA_PRIMITIVE_SPHERE: memcpy(&d->sphere, v, 4 * sizeof(double)); break;
        case ALEA_PRIMITIVE_CYLINDER_X: memcpy(&d->cyl_x, v, 3 * sizeof(double)); break;
        case ALEA_PRIMITIVE_CYLINDER_Y: memcpy(&d->cyl_y, v, 3 * sizeof(double)); break;
        case ALEA_PRIMITIVE_CYLINDER_Z: memcpy(&d->cyl_z, v, 3 * sizeof(double)); break;
        case ALEA_PRIMITIVE_CONE_X:
            if (v[4] != trunc(v[4]) || v[4] < -1 || v[4] > 1) return -1;
            memcpy(&d->cone_x, v, 4 * sizeof(double)); d->cone_x.sheet_selection = (int)v[4]; break;
        case ALEA_PRIMITIVE_CONE_Y:
            if (v[4] != trunc(v[4]) || v[4] < -1 || v[4] > 1) return -1;
            memcpy(&d->cone_y, v, 4 * sizeof(double)); d->cone_y.sheet_selection = (int)v[4]; break;
        case ALEA_PRIMITIVE_CONE_Z:
            if (v[4] != trunc(v[4]) || v[4] < -1 || v[4] > 1) return -1;
            memcpy(&d->cone_z, v, 4 * sizeof(double)); d->cone_z.sheet_selection = (int)v[4]; break;
        case ALEA_PRIMITIVE_RPP: memcpy(&d->box, v, 6 * sizeof(double)); break;
        case ALEA_PRIMITIVE_QUADRIC: memcpy(d->quadric.coeffs, v, 10 * sizeof(double)); break;
        case ALEA_PRIMITIVE_TORUS_X: case ALEA_PRIMITIVE_TORUS_Y: case ALEA_PRIMITIVE_TORUS_Z:
            d->torus.axis = type == ALEA_PRIMITIVE_TORUS_X ? ALEA_AXIS_X :
                            type == ALEA_PRIMITIVE_TORUS_Y ? ALEA_AXIS_Y : ALEA_AXIS_Z;
            d->torus.center_x=v[0]; d->torus.center_y=v[1]; d->torus.center_z=v[2];
            d->torus.major_radius=v[3]; d->torus.minor_radius=v[4]; d->torus.axial_semiwidth_B=v[5]; break;
        case ALEA_PRIMITIVE_RCC: memcpy(&d->rcc, v, 7 * sizeof(double)); break;
        case ALEA_PRIMITIVE_BOX: memcpy(&d->box_general, v, 12 * sizeof(double)); break;
        case ALEA_PRIMITIVE_SPH: memcpy(&d->sph, v, 4 * sizeof(double)); break;
        case ALEA_PRIMITIVE_TRC: memcpy(&d->trc, v, 8 * sizeof(double)); break;
        case ALEA_PRIMITIVE_ELL: memcpy(&d->ell, v, 7 * sizeof(double)); break;
        case ALEA_PRIMITIVE_REC: memcpy(&d->rec, v, 12 * sizeof(double)); break;
        case ALEA_PRIMITIVE_WED: memcpy(&d->wed, v, 12 * sizeof(double)); break;
        case ALEA_PRIMITIVE_RHP: memcpy(&d->rhp, v, 15 * sizeof(double)); break;
        case ALEA_PRIMITIVE_ARB:
            if (v[0] != trunc(v[0]) || v[1] != trunc(v[1])) return -1;
            d->arb.num_corners=(int)v[0]; d->arb.num_faces=(int)v[1];
            if (d->arb.num_corners < 4 || d->arb.num_corners > 8 ||
                d->arb.num_faces < 4 || d->arb.num_faces > 6) return -1;
            memcpy(d->arb.corners, v + 2, 24 * sizeof(double));
            for (int i=0; i<24; i++) {
                if (v[26+i] != trunc(v[26+i]) || v[26+i] < 0 || v[26+i] > 8) return -1;
                ((int*)d->arb.faces)[i]=(int)v[26+i];
            }
            break;
        default: return -1;
    }
    return 0;
}

static size_t values_from_data(alea_primitive_type_t type,
                               const alea_primitive_data_t* d, double* v) {
    size_t n = 0;
    switch (type) {
        case ALEA_PRIMITIVE_PLANE: n=4; memcpy(v,&d->plane,n*sizeof(double)); break;
        case ALEA_PRIMITIVE_SPHERE: n=4; memcpy(v,&d->sphere,n*sizeof(double)); break;
        case ALEA_PRIMITIVE_CYLINDER_X: n=3; memcpy(v,&d->cyl_x,n*sizeof(double)); break;
        case ALEA_PRIMITIVE_CYLINDER_Y: n=3; memcpy(v,&d->cyl_y,n*sizeof(double)); break;
        case ALEA_PRIMITIVE_CYLINDER_Z: n=3; memcpy(v,&d->cyl_z,n*sizeof(double)); break;
        case ALEA_PRIMITIVE_CONE_X: n=5; memcpy(v,&d->cone_x,4*sizeof(double)); v[4]=d->cone_x.sheet_selection; break;
        case ALEA_PRIMITIVE_CONE_Y: n=5; memcpy(v,&d->cone_y,4*sizeof(double)); v[4]=d->cone_y.sheet_selection; break;
        case ALEA_PRIMITIVE_CONE_Z: n=5; memcpy(v,&d->cone_z,4*sizeof(double)); v[4]=d->cone_z.sheet_selection; break;
        case ALEA_PRIMITIVE_RPP: n=6; memcpy(v,&d->box,n*sizeof(double)); break;
        case ALEA_PRIMITIVE_QUADRIC: n=10; memcpy(v,d->quadric.coeffs,n*sizeof(double)); break;
        case ALEA_PRIMITIVE_TORUS_X: case ALEA_PRIMITIVE_TORUS_Y: case ALEA_PRIMITIVE_TORUS_Z:
            n=6; v[0]=d->torus.center_x; v[1]=d->torus.center_y; v[2]=d->torus.center_z;
            v[3]=d->torus.major_radius; v[4]=d->torus.minor_radius; v[5]=d->torus.axial_semiwidth_B; break;
        case ALEA_PRIMITIVE_RCC: n=7; memcpy(v,&d->rcc,n*sizeof(double)); break;
        case ALEA_PRIMITIVE_BOX: n=12; memcpy(v,&d->box_general,n*sizeof(double)); break;
        case ALEA_PRIMITIVE_SPH: n=4; memcpy(v,&d->sph,n*sizeof(double)); break;
        case ALEA_PRIMITIVE_TRC: n=8; memcpy(v,&d->trc,n*sizeof(double)); break;
        case ALEA_PRIMITIVE_ELL: n=7; memcpy(v,&d->ell,n*sizeof(double)); break;
        case ALEA_PRIMITIVE_REC: n=12; memcpy(v,&d->rec,n*sizeof(double)); break;
        case ALEA_PRIMITIVE_WED: n=12; memcpy(v,&d->wed,n*sizeof(double)); break;
        case ALEA_PRIMITIVE_RHP: n=15; memcpy(v,&d->rhp,n*sizeof(double)); break;
        case ALEA_PRIMITIVE_ARB:
            n=50; v[0]=d->arb.num_corners; v[1]=d->arb.num_faces;
            memcpy(v+2,d->arb.corners,24*sizeof(double));
            for (int i=0;i<24;i++) v[26+i]=((const int*)d->arb.faces)[i];
            break;
        default: break;
    }
    return n;
}

static int add_surface(alea_system_t* sys, int id, alea_primitive_type_t type,
                       alea_primitive_data_t* data,
                       alea_xml_region_ctx_t* region) {
    int8_t inverted = 0;
    alea_primitive_id_t prim = alea_get_or_create_primitive(sys, type, data, &inverted);
    alea_surface_entry_t* entry;
    alea_node_id_t pos, neg;
    if (prim == ALEA_PRIMITIVE_ID_INVALID) return -1;
    pos = alea_add_primitive_node(sys, prim, +1, inverted, id);
    neg = alea_add_primitive_node(sys, prim, -1, inverted, id);
    if (pos == ALEA_NODE_ID_INVALID || neg == ALEA_NODE_ID_INVALID) return -1;
    entry = alea_vec_push_uninit(&sys->surfaces, alea_surface_entry_t);
    if (!entry) return -1;
    *entry = (alea_surface_entry_t){id, prim, pos, neg, 0, false,
        ALEA_BOUNDARY_TRANSMISSIVE, 0, ALEA_NODE_ID_INVALID, ALEA_NODE_ID_INVALID};
    if (id >= sys->next_auto_surface_id) sys->next_auto_surface_id = id + 1;
    return alea_xml_region_register_surface(region, id, pos, neg);
}

static alea_boundary_type_t boundary_from_name(const char* name) {
    if (!name || !strcmp(name,"transmissive")) return ALEA_BOUNDARY_TRANSMISSIVE;
    if (!strcmp(name,"reflective")) return ALEA_BOUNDARY_REFLECTIVE;
    if (!strcmp(name,"white")) return ALEA_BOUNDARY_WHITE;
    if (!strcmp(name,"periodic")) return ALEA_BOUNDARY_PERIODIC;
    if (!strcmp(name,"vacuum")) return ALEA_BOUNDARY_VACUUM;
    return (alea_boundary_type_t)-1;
}

static const char* boundary_name(alea_boundary_type_t value) {
    switch (value) {
        case ALEA_BOUNDARY_REFLECTIVE: return "reflective";
        case ALEA_BOUNDARY_WHITE: return "white";
        case ALEA_BOUNDARY_PERIODIC: return "periodic";
        case ALEA_BOUNDARY_VACUUM: return "vacuum";
        default: return "transmissive";
    }
}

static int parse_materials(alea_system_t* sys, const alea_xml_dom_element_t* root) {
    alea_xml_dom_element_t* section = alea_xml_dom_find_child(root, "materials");
    if (!section) return 0;
    for (size_t i=0; i<section->child_count; i++) {
        alea_xml_dom_element_t* e=section->children[i];
        if (strcmp(e->tag_name,"material")) continue;
        int id, index; const char* a=alea_xml_dom_get_attr(e,"id");
        if (parse_int_strict(a,&id) || id <= 0 || alea_find_material_by_id(sys,id)>=0)
            return fail_element(e,"invalid or duplicate id");
        index=alea_add_material(sys,id); if(index<0) return -1;
        alea_material_t* mat=&sys->materials.data[index];
        mat->name=copy_string(alea_xml_dom_get_attr(e,"name"));
        mat->comments=copy_string(alea_xml_dom_get_attr(e,"comments"));
        const char* basis=alea_xml_dom_get_attr(e,"fraction_basis");
        if (basis && strcmp(basis,"atom") && strcmp(basis,"weight")) return fail_element(e,"fraction_basis must be atom or weight");
        mat->is_weight_fraction=basis && !strcmp(basis,"weight");
        const char* density=alea_xml_dom_get_attr(e,"density");
        if (density) {
            double d; const char* units=alea_xml_dom_get_attr(e,"density_units");
            if(parse_double_strict(density,&d)||d<0||!units) return fail_element(e,"invalid density");
            if(!strcmp(units,"g/cm3")) alea_material_set_density(sys,index,d);
            else if(!strcmp(units,"atom/b-cm")) alea_material_set_density(sys,index,-d);
            else return fail_element(e,"density_units must be g/cm3 or atom/b-cm");
        }
        for(size_t j=0;j<e->child_count;j++) {
            alea_xml_dom_element_t* c=e->children[j];
            if(!strcmp(c->tag_name,"nuclide")) {
                int zaid; double fraction; char* library;
                if(parse_int_strict(alea_xml_dom_get_attr(c,"zaid"),&zaid)||zaid<=0||
                   parse_double_strict(alea_xml_dom_get_attr(c,"fraction"),&fraction)||fraction<0)
                    return fail_element(c,"invalid zaid or fraction");
                library=library_internal(alea_xml_dom_get_attr(c,"library"));
                if(!library) return fail_element(c,"each nuclide requires a library extension");
                int rc=alea_material_add_nuclide(sys,index,zaid,library,fraction); free(library);
                if(rc) return -1;
            } else if(!strcmp(c->tag_name,"element")) {
                int z; double fraction; char* library;
                if(parse_int_strict(alea_xml_dom_get_attr(c,"z"),&z)||z<1||z>118||
                   parse_double_strict(alea_xml_dom_get_attr(c,"fraction"),&fraction)||fraction<0)
                    return fail_element(c,"invalid z or fraction");
                library=library_internal(alea_xml_dom_get_attr(c,"library"));
                if(!library) return fail_element(c,"each element requires a library extension");
                int rc=alea_material_add_element(sys,index,z,library,fraction); free(library);
                if(rc) return -1;
            } else if(!strcmp(c->tag_name,"thermal")) {
                const char* identifier=alea_xml_dom_get_attr(c,"identifier"); int zaid=0;
                if(!identifier||!*identifier) return fail_element(c,"identifier is required");
                if(alea_xml_dom_get_attr(c,"zaid_match") && parse_int_strict(alea_xml_dom_get_attr(c,"zaid_match"),&zaid))
                    return fail_element(c,"invalid zaid_match");
                alea_thermal_law_t* law=alea_vec_push_uninit(&mat->thermal_laws,alea_thermal_law_t);
                if(!law) return -1;
                law->identifier=alea_strdup(identifier); law->zaid_match=zaid;
                if(!law->identifier) return -1;
            }
        }
    }
    for (size_t i=0; i<section->child_count; i++) {
        alea_xml_dom_element_t* e=section->children[i];
        if (strcmp(e->tag_name,"mixture")) continue;
        int id, mc_id=0; const char* basis=alea_xml_dom_get_attr(e,"fraction_basis");
        if(parse_int_strict(alea_xml_dom_get_attr(e,"id"),&id)||id<=0||alea_find_mixture_by_id(sys,id)>=0)
            return fail_element(e,"invalid or duplicate id");
        alea_mixture_t* mix=alea_mixture_create(id); if(!mix) return -1;
        if(basis && strcmp(basis,"atom")&&strcmp(basis,"weight")){alea_mixture_destroy(mix);return fail_element(e,"invalid fraction_basis");}
        mix->is_weight_fraction=basis && !strcmp(basis,"weight");
        if(alea_xml_dom_get_attr(e,"mc_material_id")&&parse_int_strict(alea_xml_dom_get_attr(e,"mc_material_id"),&mc_id)){alea_mixture_destroy(mix);return fail_element(e,"invalid mc_material_id");}
        mix->mc_material_id=mc_id; mix->name=copy_string(alea_xml_dom_get_attr(e,"name")); mix->comments=copy_string(alea_xml_dom_get_attr(e,"comments"));
        for(size_t j=0;j<e->child_count;j++) if(!strcmp(e->children[j]->tag_name,"component")) {
            int material; double fraction; alea_xml_dom_element_t* c=e->children[j];
            if(parse_int_strict(alea_xml_dom_get_attr(c,"material"),&material)||alea_find_material_by_id(sys,material)<0||
               parse_double_strict(alea_xml_dom_get_attr(c,"fraction"),&fraction)||fraction<0||
               alea_mixture_add_component(mix,material,fraction)<0){alea_mixture_destroy(mix);return fail_element(c,"invalid component");}
        }
        if(alea_add_mixture(sys,mix)<0){alea_mixture_destroy(mix);return -1;} alea_mixture_destroy(mix);
    }
    return 0;
}

static int parse_transforms(alea_system_t* sys, const alea_xml_dom_element_t* root) {
    alea_xml_dom_element_t* section=alea_xml_dom_find_child(root,"transforms");
    if(!section)return 0;
    for(size_t i=0;i<section->child_count;i++) {
        alea_xml_dom_element_t* e=section->children[i]; if(strcmp(e->tag_name,"transform"))continue;
        int id,degrees=0,inline_flag=0; double v[12]; const char* values=alea_xml_dom_get_attr(e,"values");
        if(parse_int_strict(alea_xml_dom_get_attr(e,"id"),&id)||id<=0||id==INT_MAX||
           parse_bool(alea_xml_dom_get_attr(e,"degrees"),0,&degrees)||parse_bool(alea_xml_dom_get_attr(e,"inline"),0,&inline_flag))
            return fail_element(e,"invalid transform attributes");
        if(alea_get_transform(sys,id))return fail_element(e,"duplicate transform id");
        size_t n=alea_xml_dom_parse_doubles(values,v,12); if((n!=3&&n!=12)||parse_values(values,v,n)) return fail_element(e,"values must contain 3 or 12 numbers");
        if(alea_add_transform(sys,id,v,(int)n,degrees)) return fail_element(e,"could not add transform");
        sys->transforms.data[sys->transforms.count-1].from_inline=inline_flag;
    }
    alea_finalize_transform_ids(sys); return 0;
}

static int parse_surfaces(alea_system_t* sys, const alea_xml_dom_element_t* root,
                          alea_xml_region_ctx_t* region) {
    alea_xml_dom_element_t* section=alea_xml_dom_find_child(root,"surfaces");
    if(!section)return 0;
    for(size_t i=0;i<section->child_count;i++) {
        alea_xml_dom_element_t* e=section->children[i]; if(strcmp(e->tag_name,"surface"))continue;
        int id,transform=0,applied=0,periodic=0; size_t n; double values[50];
        alea_primitive_type_t type; alea_primitive_data_t data;
        if(parse_int_strict(alea_xml_dom_get_attr(e,"id"),&id)||id<=0||id==INT_MAX||
           alea_xml_region_has_surface(region,id))
            return fail_element(e,"invalid or duplicate id");
        if(surface_type(alea_xml_dom_get_attr(e,"type"),&type,&n)||parse_values(alea_xml_dom_get_attr(e,"coeffs"),values,n)||data_from_values(type,values,&data))
            return fail_element(e,"invalid type or coefficients");
        if(add_surface(sys,id,type,&data,region)) return -1;
        alea_surface_entry_t* s=&sys->surfaces.data[sys->surfaces.count-1];
        alea_boundary_type_t boundary=boundary_from_name(alea_xml_dom_get_attr(e,"boundary"));
        if((int)boundary<0)return fail_element(e,"invalid boundary");
        s->boundary_type=boundary;
        if(alea_xml_dom_get_attr(e,"transform")&&parse_int_strict(alea_xml_dom_get_attr(e,"transform"),&transform))return fail_element(e,"invalid transform");
        if(parse_bool(alea_xml_dom_get_attr(e,"transform_applied"),0,&applied))return fail_element(e,"invalid transform_applied");
        if(transform && !alea_get_transform(sys,transform))return fail_element(e,"unknown transform");
        s->transform_id=transform; s->transform_applied=applied;
        if(alea_xml_dom_get_attr(e,"periodic_surface")&&parse_int_strict(alea_xml_dom_get_attr(e,"periodic_surface"),&periodic))return fail_element(e,"invalid periodic_surface");
        s->periodic_surface_id=periodic;
    }
    if(alea_xml_region_finalize_surfaces(region))return fail("duplicate surface id");
    for(size_t i=0;i<sys->surfaces.count;i++) {
        alea_surface_entry_t* surface=&sys->surfaces.data[i];
        if(surface->boundary_type==ALEA_BOUNDARY_PERIODIC && !surface->periodic_surface_id)
            return fail("periodic boundary requires periodic_surface");
        if(surface->periodic_surface_id) {
            int peer_id=surface->periodic_surface_id;
            alea_surface_entry_t* peer=NULL;
            if(peer_id<=0||!alea_xml_region_has_surface(region,peer_id))
                return fail("periodic surface references an unknown surface");
            for(size_t j=0;j<sys->surfaces.count;j++)
                if(sys->surfaces.data[j].mc_surface_id==peer_id){peer=&sys->surfaces.data[j];break;}
            if(!peer||surface->boundary_type!=ALEA_BOUNDARY_PERIODIC||
               peer->boundary_type!=ALEA_BOUNDARY_PERIODIC||
               peer->periodic_surface_id!=surface->mc_surface_id)
                return fail("periodic surfaces must form a reciprocal pair");
        }
    }
    return 0;
}

static int set_parameter(alea_model_cell_metadata_t* m, const char* name, const char* value) {
    double d; int v;
    if(!name||!value)return -1;
    if(!strcmp(name,"volume")){if(parse_double_strict(value,&d))return -1;m->user_volume=d;m->parameter_flags|=ALEA_CELL_PARAM_VOLUME;}
    else if(!strcmp(name,"pwt")){if(parse_double_strict(value,&d))return -1;m->photon_weight=d;m->parameter_flags|=ALEA_CELL_PARAM_PWT;}
    else if(!strcmp(name,"pd")){if(parse_double_strict(value,&d))return -1;m->detector_contribution=d;m->parameter_flags|=ALEA_CELL_PARAM_PD;}
    else if(!strcmp(name,"elpt")){if(parse_double_strict(value,&d))return -1;m->energy_cutoff=d;m->parameter_flags|=ALEA_CELL_PARAM_ELPT;}
    else if(!strcmp(name,"nonu")){if(parse_int_strict(value,&v))return -1;m->fission_turnoff=v;m->parameter_flags|=ALEA_CELL_PARAM_NONU;}
    else if(!strcmp(name,"unc")){if(parse_int_strict(value,&v))return -1;m->uncollided_secondaries=v;m->parameter_flags|=ALEA_CELL_PARAM_UNC;}
    else if(!strcmp(name,"bflcl")){if(parse_int_strict(value,&v))return -1;m->magnetic_field=v;m->parameter_flags|=ALEA_CELL_PARAM_BFLCL;}
    else return -1;
    return 0;
}

static int parse_cells(alea_system_t* sys, const alea_xml_dom_element_t* root,
                       alea_xml_region_ctx_t* region, alea_model_t** model_out) {
    alea_xml_dom_element_t* section=alea_xml_dom_find_child(root,"cells");
    if(!section)return fail("missing <cells>");
    for(size_t i=0;i<section->child_count;i++) {
        alea_xml_dom_element_t* e=section->children[i]; if(strcmp(e->tag_name,"cell"))continue;
        int id,universe=0,material=0,material_index=ALEA_MATERIAL_VOID,mass_density=0; double density=0;
        const char* region_text=alea_xml_dom_get_attr(e,"region"); const char* units=alea_xml_dom_get_attr(e,"density_units");
        if(parse_int_strict(alea_xml_dom_get_attr(e,"id"),&id)||id<=0||!region_text)return fail_element(e,"id and region are required");
        if(alea_xml_dom_get_attr(e,"universe")&&parse_int_strict(alea_xml_dom_get_attr(e,"universe"),&universe))return fail_element(e,"invalid universe");
        if(alea_xml_dom_get_attr(e,"material")) {
            const char* kind=alea_xml_dom_get_attr(e,"material_kind");
            if(parse_int_strict(alea_xml_dom_get_attr(e,"material"),&material)||material<0)return fail_element(e,"invalid material");
            if(material>0) {
                if(kind && !strcmp(kind,"mixture")) {
                    if(alea_find_mixture_by_id(sys,material)<0)return fail_element(e,"unknown mixture");
                    material_index=ALEA_MATERIAL_VOID;
                } else {
                    if(kind && strcmp(kind,"material"))return fail_element(e,"material_kind must be material or mixture");
                    material_index=alea_find_material_by_id(sys,material);
                    if(material_index<0)return fail_element(e,"unknown material");
                }
            } else if(kind) {
                return fail_element(e,"void cells cannot specify material_kind");
            }
        }
        if(material>0) {
            if(parse_double_strict(alea_xml_dom_get_attr(e,"density"),&density)||density<0||!units)return fail_element(e,"material cells require nonnegative density and units");
            if(!strcmp(units,"g/cm3")){density=-density;mass_density=1;} else if(strcmp(units,"atom/b-cm"))return fail_element(e,"invalid density_units");
        } else if(alea_xml_dom_get_attr(e,"density")||units) {
            return fail_element(e,"void cells must omit density and density_units");
        }
        alea_node_id_t node=alea_xml_parse_region(region,region_text);
        if(node==ALEA_NODE_ID_INVALID){alea_set_error_detail(ALEA_ERR_PARSE_ERROR,"ALEA XML cell %d region: %s",id,alea_xml_region_get_error(region));return -1;}
        int index=alea_add_cell_with_id(sys,id,node,material_index,density,universe); if(index<0)return -1;
        alea_cell_entry_t* cell=&sys->cells.data[index];
        if(material>0) cell->is_mass_density=(unsigned int)mass_density;
        if(material_index<0&&material>0){cell->material_id=material;cell->material_index=-1;}
        const char* comments=alea_xml_dom_get_attr(e,"comments"); const char* inline_comment=alea_xml_dom_get_attr(e,"inline_comment");
        if(comments&&alea_cell_set_comment(sys,index,comments))return -1;
        if(inline_comment&&alea_cell_set_inline_comment(sys,index,inline_comment))return -1;
        if(alea_xml_dom_get_attr(e,"temperature")&&alea_cell_set_temperature(sys,index,alea_xml_dom_get_attr_double(e,"temperature",0)))return fail_element(e,"invalid temperature");
        int fill=0,fill_transform=0;
        if(alea_xml_dom_get_attr(e,"fill")&&parse_int_strict(alea_xml_dom_get_attr(e,"fill"),&fill))return fail_element(e,"invalid fill");
        if(alea_xml_dom_get_attr(e,"fill_transform")&&parse_int_strict(alea_xml_dom_get_attr(e,"fill_transform"),&fill_transform))return fail_element(e,"invalid fill_transform");
        if(fill&&alea_set_cell_fill(sys,index,fill,fill_transform))return -1;
        if(alea_xml_dom_get_attr(e,"lattice_type")) {
            if(parse_int_strict(alea_xml_dom_get_attr(e,"lattice_type"),&cell->lat_type)||(cell->lat_type!=1&&cell->lat_type!=2))return fail_element(e,"invalid lattice_type");
            if(parse_int_values(alea_xml_dom_get_attr(e,"lattice_dims"),cell->lat_fill_dims,6))return fail_element(e,"invalid lattice_dims");
            long long nx=(long long)cell->lat_fill_dims[1]-cell->lat_fill_dims[0]+1, ny=(long long)cell->lat_fill_dims[3]-cell->lat_fill_dims[2]+1, nz=(long long)cell->lat_fill_dims[5]-cell->lat_fill_dims[4]+1;
            if(nx<=0||ny<=0||nz<=0||nx*ny*nz>INT_MAX)return fail_element(e,"invalid lattice dimensions");
            cell->lat_fill_count=(size_t)(nx*ny*nz);cell->lat_fill=malloc(cell->lat_fill_count*sizeof(int));if(!cell->lat_fill)return -1;
            if(parse_int_values(alea_xml_dom_get_attr(e,"lattice_fill"),cell->lat_fill,cell->lat_fill_count))return fail_element(e,"invalid lattice_fill");
            if(alea_xml_dom_get_attr(e,"lattice_pitch")&&parse_values(alea_xml_dom_get_attr(e,"lattice_pitch"),cell->lat_pitch,3))return fail_element(e,"invalid lattice_pitch");
            if(alea_xml_dom_get_attr(e,"lattice_lower_left")&&parse_values(alea_xml_dom_get_attr(e,"lattice_lower_left"),cell->lat_lower_left,3))return fail_element(e,"invalid lattice_lower_left");
            cell->lat_outer_universe=-1;if(alea_xml_dom_get_attr(e,"lattice_outer")&&parse_int_strict(alea_xml_dom_get_attr(e,"lattice_outer"),&cell->lat_outer_universe))return fail_element(e,"invalid lattice_outer");
            int b;if(parse_bool(alea_xml_dom_get_attr(e,"lattice_repeating"),0,&b))return fail_element(e,"invalid lattice_repeating");cell->lat_fill_repeating=b;
            if(parse_bool(alea_xml_dom_get_attr(e,"lattice_zero_coords"),0,&b))return fail_element(e,"invalid lattice_zero_coords");
            cell->lat_fill_zero_element_coords=b;sys->has_lattice=true;
        }
        /* bbox is deliberately parsed only for validation; geometry recomputes it. */
        if(alea_xml_dom_get_attr(e,"bbox")){double bbox[6];if(parse_values(alea_xml_dom_get_attr(e,"bbox"),bbox,6)||bbox[0]>bbox[1]||bbox[2]>bbox[3]||bbox[4]>bbox[5])return fail_element(e,"invalid bbox hint");}
    }
    if(alea_validate_cell_ids(sys))return -1;
    alea_model_t* model=alea_model_adopt(sys);if(!model)return -1;
    *model_out=model;
    for(size_t i=0,cell_index=0;i<section->child_count;i++) {
        alea_xml_dom_element_t* e=section->children[i];if(strcmp(e->tag_name,"cell"))continue;
        alea_model_cell_metadata_t* m=alea_model_cell_metadata_mut(model,cell_index++);
        if(alea_model_cell_set_name(model,cell_index-1,alea_xml_dom_get_attr(e,"name"))<0)return -1;
        for(size_t j=0;j<e->child_count;j++) {
            alea_xml_dom_element_t* c=e->children[j];
            if(!strcmp(c->tag_name,"importance")) {
                const char* particle=alea_xml_dom_get_attr(c,"particle");double value;
                if(parse_double_strict(alea_xml_dom_get_attr(c,"value"),&value)||value<0)return fail_element(c,"invalid importance");
                if(particle&&!strcmp(particle,"neutron")){m->importance_neutron=value;m->has_importance_neutron=1;}
                else if(particle&&!strcmp(particle,"photon")){m->importance_photon=value;m->has_importance_photon=1;}
                else if(particle&&!strcmp(particle,"electron")){m->importance_electron=value;m->has_importance_electron=1;}
                else return fail_element(c,"invalid particle");
            } else if(!strcmp(c->tag_name,"parameter")&&set_parameter(m,alea_xml_dom_get_attr(c,"name"),alea_xml_dom_get_attr(c,"value"))) {
                return fail_element(c,"unknown or invalid parameter");
            }
        }
    }
    return 0;
}

static alea_model_t* load_doc(alea_xml_dom_doc_t* doc) {
    alea_system_t* sys=NULL; alea_model_t* model=NULL; arena_t region_arena; int arena_ok=0;
    if(!doc){alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,"ALEA XML: parser allocation failed");return NULL;}
    if(alea_xml_dom_get_error(doc)){alea_set_error_detail(ALEA_ERR_PARSE_ERROR,"ALEA XML: %s",alea_xml_dom_get_error(doc));alea_xml_dom_doc_free(doc);return NULL;}
    alea_xml_dom_element_t* root=doc->root;
    if(!root||strcmp(root->tag_name,"alea")||strcmp(alea_xml_dom_get_attr(root,"version")?alea_xml_dom_get_attr(root,"version"):"",ALEA_XML_VERSION)||
       strcmp(alea_xml_dom_get_attr(root,"length_units")?alea_xml_dom_get_attr(root,"length_units"):"","cm")) {
        fail("root must be <alea version=\"1\" length_units=\"cm\">");alea_xml_dom_doc_free(doc);return NULL;
    }
    if(validate_structure(root)){alea_xml_dom_doc_free(doc);return NULL;}
    sys=alea_create();if(!sys)goto done;
    if(!arena_init(&region_arena))goto done;
    arena_ok=1;
    alea_xml_region_ctx_t region;alea_xml_region_ctx_init(&region,sys,&region_arena);
    if(parse_materials(sys,root)||parse_transforms(sys,root)||parse_surfaces(sys,root,&region))goto done;
    if(parse_cells(sys,root,&region,&model)) {
        if(model) { sys=NULL; alea_model_destroy(model); model=NULL; }
        goto done;
    }
    sys=NULL;
    alea_model_system(model)->source=ALEA_SOURCE_ALEA_XML;
    if(alea_model_set_name(model,alea_xml_dom_get_attr(root,"name"))||alea_model_set_title(model,alea_xml_dom_get_attr(root,"title"))||alea_model_set_comments(model,alea_xml_dom_get_attr(root,"comments"))){alea_model_destroy(model);model=NULL;goto done;}
    if(alea_build_universe_index(alea_model_system(model))){alea_model_destroy(model);model=NULL;}
done:
    if(arena_ok)arena_free(&region_arena);
    if(sys)alea_destroy(sys);
    alea_xml_dom_doc_free(doc);
    return model;
}

alea_model_t* alea_xml_load(const char* filename) {
    if(!filename){alea_set_error_detail(ALEA_ERR_NULL_ARG,"alea_xml_load: filename is NULL");return NULL;}
    return load_doc(alea_xml_dom_parse_file(filename));
}

alea_model_t* alea_xml_load_string(const char* xml,size_t length) {
    if(!xml){alea_set_error_detail(ALEA_ERR_NULL_ARG,"alea_xml_load_string: xml is NULL");return NULL;}
    return load_doc(alea_xml_dom_parse_string(xml,length));
}

static int append_values(char* out,size_t capacity,const double* v,size_t n) {
    size_t used=0;
    for(size_t i=0;i<n;i++){int w=snprintf(out+used,capacity-used,i?" %.17g":"%.17g",v[i]);if(w<0||(size_t)w>=capacity-used)return -1;used+=(size_t)w;}return 0;
}

static int append_int_values(char* out,size_t capacity,const int* v,size_t n) {
    size_t used=0;
    for(size_t i=0;i<n;i++){int w=snprintf(out+used,capacity-used,i?" %d":"%d",v[i]);if(w<0||(size_t)w>=capacity-used)return -1;used+=(size_t)w;}return 0;
}

static int write_attr_optional(alea_xml_writer_t* x,const char* name,const char* value){return !value||!*value||alea_xml_writer_attribute(x,name,value);}

static int node_expr(const alea_system_t* sys,alea_node_id_t id,str_builder_t* sb) {
    if(id>=sys->nodes.count)return -1;
    const alea_node_t* n=&sys->nodes.data[id];
    alea_operation_t op=ALEA_GET_OPERATION(n);
    if(op==ALEA_OP_PRIMITIVE){int sid=n->primitive.mc_surface_id;int flip=(n->primitive.sense>0)!=(n->primitive.inverted!=0);return str_builder_int(sb,flip?sid:-sid)?0:-1;}
    if(op==ALEA_OP_COMPLEMENT){if(!str_builder_puts(sb,"~(" )||node_expr(sys,n->operation.left,sb)||!str_builder_putc(sb,')'))return -1;return 0;}
    if(!str_builder_putc(sb,'(')||node_expr(sys,n->operation.left,sb))return -1;
    if(op==ALEA_OP_UNION){if(!str_builder_puts(sb," | "))return -1;}
    else if(op==ALEA_OP_INTERSECTION){if(!str_builder_putc(sb,' '))return -1;}
    else if(op==ALEA_OP_DIFFERENCE){if(!str_builder_puts(sb," ~(" )||node_expr(sys,n->operation.right,sb)||!str_builder_putc(sb,')')||!str_builder_putc(sb,')'))return -1;return 0;}
    else return -1;
    return node_expr(sys,n->operation.right,sb)||!str_builder_putc(sb,')')?-1:0;
}

#define XOK(expr) do { if (!(expr)) return -1; } while (0)

static int preflight_export(const alea_model_t* model, const alea_system_t* sys) {
    if (model && alea_model_cell_metadata_count(model) != sys->cells.count) {
        alea_set_error_detail(ALEA_ERR_EXPORT_FAILED,
                              "ALEA XML: cell metadata is out of sync");
        return -1;
    }
    for (size_t i = 0; i < sys->materials.count; i++) {
        const alea_material_t* material = &sys->materials.data[i];
        for (size_t j = 0; j < material->nuclides.count; j++) {
            if (!material->nuclides.data[j].library ||
                !material->nuclides.data[j].library[0]) {
                alea_set_error_detail(ALEA_ERR_EXPORT_FAILED,
                    "ALEA XML: nuclide %d in material %d has no library extension",
                    material->nuclides.data[j].zaid, material->material_id);
                return -1;
            }
        }
        for (size_t j = 0; j < material->elements.count; j++) {
            if (!material->elements.data[j].library ||
                !material->elements.data[j].library[0]) {
                alea_set_error_detail(ALEA_ERR_EXPORT_FAILED,
                    "ALEA XML: element %d in material %d has no library extension",
                    material->elements.data[j].atomic_number, material->material_id);
                return -1;
            }
        }
    }
    for (size_t i = 0; i < sys->surfaces.count; i++) {
        const alea_surface_entry_t* surface = &sys->surfaces.data[i];
        if (surface->mc_surface_id <= 0 ||
            !surface_type_name(sys->primitives.data[surface->primitive_id].type)) {
            alea_set_error_detail(ALEA_ERR_EXPORT_FAILED,
                                  "ALEA XML: surface %zu is not serializable", i);
            return -1;
        }
    }
    for (size_t i = 0; i < sys->cells.count; i++) {
        if (sys->cells.data[i].root_node_id >= sys->nodes.count) {
            alea_set_error_detail(ALEA_ERR_EXPORT_FAILED,
                                  "ALEA XML: cell %d has an invalid region",
                                  sys->cells.data[i].mc_cell_id);
            return -1;
        }
    }
    return 0;
}

static int write_model(const alea_model_t* model,const alea_system_t* sys,FILE* stream) {
    if (preflight_export(model, sys)) return -1;
    alea_xml_writer_t x;alea_xml_writer_init_stream_ex(&x,stream,16384,120,2,true);
    XOK(alea_xml_writer_start_document(&x,"1.0","utf-8"));XOK(alea_xml_writer_start_element(&x,"alea"));
    XOK(alea_xml_writer_attribute(&x,"version",ALEA_XML_VERSION));XOK(alea_xml_writer_attribute(&x,"length_units","cm"));
    if(model){XOK(write_attr_optional(&x,"name",alea_model_name(model)));XOK(write_attr_optional(&x,"title",alea_model_title(model)));XOK(write_attr_optional(&x,"comments",alea_model_comments(model)));}
    XOK(alea_xml_writer_end_start_tag(&x,false));
    XOK(alea_xml_writer_start_element(&x,"materials"));XOK(alea_xml_writer_end_start_tag(&x,false));
    for(size_t i=0;i<sys->materials.count;i++) {
        const alea_material_t* m=&sys->materials.data[i];XOK(alea_xml_writer_start_element(&x,"material"));XOK(alea_xml_writer_attribute_i(&x,"id",m->material_id));
        XOK(alea_xml_writer_attribute(&x,"fraction_basis",m->is_weight_fraction?"weight":"atom"));XOK(write_attr_optional(&x,"name",m->name));XOK(write_attr_optional(&x,"comments",m->comments));
        if(m->has_standard_density){XOK(alea_xml_writer_attribute_f(&x,"density",fabs(m->standard_density),17,false));XOK(alea_xml_writer_attribute(&x,"density_units",m->standard_density<0?"atom/b-cm":"g/cm3"));}
        XOK(alea_xml_writer_end_start_tag(&x,false));
        for(size_t j=0;j<m->nuclides.count;j++){const alea_nuclide_t* n=&m->nuclides.data[j];XOK(alea_xml_writer_start_element(&x,"nuclide"));XOK(alea_xml_writer_attribute_i(&x,"zaid",n->zaid));XOK(alea_xml_writer_attribute_f(&x,"fraction",n->fraction,17,false));XOK(alea_xml_writer_attribute(&x,"library",n->library[0]=='.'?n->library+1:n->library));XOK(alea_xml_writer_end_start_tag(&x,true));}
        for(size_t j=0;j<m->elements.count;j++){const alea_element_comp_t* n=&m->elements.data[j];XOK(alea_xml_writer_start_element(&x,"element"));XOK(alea_xml_writer_attribute_i(&x,"z",n->atomic_number));XOK(alea_xml_writer_attribute_f(&x,"fraction",n->fraction,17,false));XOK(alea_xml_writer_attribute(&x,"library",n->library[0]=='.'?n->library+1:n->library));XOK(alea_xml_writer_end_start_tag(&x,true));}
        for(size_t j=0;j<m->thermal_laws.count;j++){XOK(alea_xml_writer_start_element(&x,"thermal"));XOK(alea_xml_writer_attribute(&x,"identifier",m->thermal_laws.data[j].identifier));if(m->thermal_laws.data[j].zaid_match)XOK(alea_xml_writer_attribute_i(&x,"zaid_match",m->thermal_laws.data[j].zaid_match));XOK(alea_xml_writer_end_start_tag(&x,true));}
        XOK(alea_xml_writer_end_element(&x,"material"));
    }
    for(size_t i=0;i<sys->mixtures.count;i++){const alea_mixture_t* m=&sys->mixtures.data[i];XOK(alea_xml_writer_start_element(&x,"mixture"));XOK(alea_xml_writer_attribute_i(&x,"id",m->mixture_id));XOK(alea_xml_writer_attribute(&x,"fraction_basis",m->is_weight_fraction?"weight":"atom"));if(m->mc_material_id)XOK(alea_xml_writer_attribute_i(&x,"mc_material_id",m->mc_material_id));XOK(write_attr_optional(&x,"name",m->name));XOK(write_attr_optional(&x,"comments",m->comments));XOK(alea_xml_writer_end_start_tag(&x,false));for(size_t j=0;j<m->components.count;j++){XOK(alea_xml_writer_start_element(&x,"component"));XOK(alea_xml_writer_attribute_i(&x,"material",m->components.data[j].material_id));XOK(alea_xml_writer_attribute_f(&x,"fraction",m->components.data[j].fraction,17,false));XOK(alea_xml_writer_end_start_tag(&x,true));}XOK(alea_xml_writer_end_element(&x,"mixture"));}
    XOK(alea_xml_writer_end_element(&x,"materials"));
    if(sys->transforms.count){XOK(alea_xml_writer_start_element(&x,"transforms"));XOK(alea_xml_writer_end_start_tag(&x,false));for(size_t i=0;i<sys->transforms.count;i++){char values[512];const alea_transform_t* t=&sys->transforms.data[i];if(append_values(values,sizeof(values),t->data,(size_t)t->value_count))return -1;XOK(alea_xml_writer_start_element(&x,"transform"));XOK(alea_xml_writer_attribute_i(&x,"id",t->transform_id));XOK(alea_xml_writer_attribute(&x,"values",values));if(t->degrees)XOK(alea_xml_writer_attribute(&x,"degrees","true"));if(t->from_inline)XOK(alea_xml_writer_attribute(&x,"inline","true"));XOK(alea_xml_writer_end_start_tag(&x,true));}XOK(alea_xml_writer_end_element(&x,"transforms"));}
    XOK(alea_xml_writer_start_element(&x,"surfaces"));XOK(alea_xml_writer_end_start_tag(&x,false));
    for(size_t i=0;i<sys->surfaces.count;i++){const alea_surface_entry_t* s=&sys->surfaces.data[i];const alea_primitive_entry_t* p=&sys->primitives.data[s->primitive_id];alea_primitive_data_t d;double v[50];char values[1400];const char* type=surface_type_name(p->type);if(!type||!alea_primitive_copy_data(sys,s->primitive_id,&d))return -1;size_t n=values_from_data(p->type,&d,v);if(!n||append_values(values,sizeof(values),v,n))return -1;XOK(alea_xml_writer_start_element(&x,"surface"));XOK(alea_xml_writer_attribute_i(&x,"id",s->mc_surface_id));XOK(alea_xml_writer_attribute(&x,"type",type));XOK(alea_xml_writer_attribute(&x,"coeffs",values));if(s->boundary_type!=ALEA_BOUNDARY_TRANSMISSIVE)XOK(alea_xml_writer_attribute(&x,"boundary",boundary_name(s->boundary_type)));if(s->periodic_surface_id)XOK(alea_xml_writer_attribute_i(&x,"periodic_surface",s->periodic_surface_id));if(s->transform_id)XOK(alea_xml_writer_attribute_i(&x,"transform",s->transform_id));if(s->transform_applied)XOK(alea_xml_writer_attribute(&x,"transform_applied","true"));XOK(alea_xml_writer_end_start_tag(&x,true));}
    XOK(alea_xml_writer_end_element(&x,"surfaces"));
    XOK(alea_xml_writer_start_element(&x,"cells"));XOK(alea_xml_writer_end_start_tag(&x,false));
    for(size_t i=0;i<sys->cells.count;i++){const alea_cell_entry_t* c=&sys->cells.data[i];arena_t a;if(!arena_init(&a))return -1;str_builder_t sb;str_builder_init(&sb,&a,128);if(node_expr(sys,c->root_node_id,&sb)){arena_free(&a);return -1;}const char* region=str_builder_get(&sb);XOK(alea_xml_writer_start_element(&x,"cell"));XOK(alea_xml_writer_attribute_i(&x,"id",c->mc_cell_id));XOK(alea_xml_writer_attribute_i(&x,"universe",c->universe_id));XOK(alea_xml_writer_attribute_i(&x,"material",c->material_id));if(c->material_id>0){if(c->material_index<0&&alea_find_mixture_by_id(sys,c->material_id)>=0)XOK(alea_xml_writer_attribute(&x,"material_kind","mixture"));XOK(alea_xml_writer_attribute_f(&x,"density",c->density,17,false));XOK(alea_xml_writer_attribute(&x,"density_units",c->is_mass_density?"g/cm3":"atom/b-cm"));}XOK(alea_xml_writer_attribute(&x,"region",region));alea_bbox_t bbox=alea_node_bbox_get(&sys->nodes.data[c->root_node_id].bbox);if(isfinite(bbox.min_x)&&isfinite(bbox.max_x)&&isfinite(bbox.min_y)&&isfinite(bbox.max_y)&&isfinite(bbox.min_z)&&isfinite(bbox.max_z)&&bbox.min_x<=bbox.max_x&&bbox.min_y<=bbox.max_y&&bbox.min_z<=bbox.max_z){double bv[6]={bbox.min_x,bbox.max_x,bbox.min_y,bbox.max_y,bbox.min_z,bbox.max_z};char bbox_text[384];if(append_values(bbox_text,sizeof(bbox_text),bv,6)){arena_free(&a);return -1;}XOK(alea_xml_writer_attribute(&x,"bbox",bbox_text));}
        const alea_model_cell_metadata_t* m=model?alea_model_cell_metadata(model,i):NULL;if(m)XOK(write_attr_optional(&x,"name",m->name));XOK(write_attr_optional(&x,"comments",c->comments));XOK(write_attr_optional(&x,"inline_comment",c->inline_comment));if(c->has_temperature)XOK(alea_xml_writer_attribute_f(&x,"temperature",c->temperature,17,false));if(c->fill_universe>0)XOK(alea_xml_writer_attribute_i(&x,"fill",c->fill_universe));if(c->fill_transform)XOK(alea_xml_writer_attribute_i(&x,"fill_transform",c->fill_transform));
        if(c->lat_type){char dims[128],pitch[256],lower[256];str_builder_t fill_sb;str_builder_init(&fill_sb,&a,128);for(size_t j=0;j<c->lat_fill_count;j++){if(j&&!str_builder_putc(&fill_sb,' ')){arena_free(&a);return -1;}if(!str_builder_int(&fill_sb,c->lat_fill[j])){arena_free(&a);return -1;}}const char* fill=str_builder_get(&fill_sb);if(append_int_values(dims,sizeof(dims),c->lat_fill_dims,6)||append_values(pitch,sizeof(pitch),c->lat_pitch,3)||append_values(lower,sizeof(lower),c->lat_lower_left,3)){arena_free(&a);return -1;}XOK(alea_xml_writer_attribute_i(&x,"lattice_type",c->lat_type));XOK(alea_xml_writer_attribute(&x,"lattice_dims",dims));XOK(alea_xml_writer_attribute(&x,"lattice_fill",fill));XOK(alea_xml_writer_attribute(&x,"lattice_pitch",pitch));XOK(alea_xml_writer_attribute(&x,"lattice_lower_left",lower));if(c->lat_outer_universe>=0)XOK(alea_xml_writer_attribute_i(&x,"lattice_outer",c->lat_outer_universe));if(c->lat_fill_repeating)XOK(alea_xml_writer_attribute(&x,"lattice_repeating","true"));if(c->lat_fill_zero_element_coords)XOK(alea_xml_writer_attribute(&x,"lattice_zero_coords","true"));}
        if(m&&(m->has_importance_neutron||m->has_importance_photon||m->has_importance_electron||m->parameter_flags)){XOK(alea_xml_writer_end_start_tag(&x,false));if(m->has_importance_neutron){XOK(alea_xml_writer_start_element(&x,"importance"));XOK(alea_xml_writer_attribute(&x,"particle","neutron"));XOK(alea_xml_writer_attribute_f(&x,"value",m->importance_neutron,17,false));XOK(alea_xml_writer_end_start_tag(&x,true));}if(m->has_importance_photon){XOK(alea_xml_writer_start_element(&x,"importance"));XOK(alea_xml_writer_attribute(&x,"particle","photon"));XOK(alea_xml_writer_attribute_f(&x,"value",m->importance_photon,17,false));XOK(alea_xml_writer_end_start_tag(&x,true));}if(m->has_importance_electron){XOK(alea_xml_writer_start_element(&x,"importance"));XOK(alea_xml_writer_attribute(&x,"particle","electron"));XOK(alea_xml_writer_attribute_f(&x,"value",m->importance_electron,17,false));XOK(alea_xml_writer_end_start_tag(&x,true));}struct {uint32_t flag;const char* name;double value;int integer;} params[]={{ALEA_CELL_PARAM_VOLUME,"volume",m->user_volume,0},{ALEA_CELL_PARAM_PWT,"pwt",m->photon_weight,0},{ALEA_CELL_PARAM_PD,"pd",m->detector_contribution,0},{ALEA_CELL_PARAM_ELPT,"elpt",m->energy_cutoff,0},{ALEA_CELL_PARAM_NONU,"nonu",m->fission_turnoff,1},{ALEA_CELL_PARAM_UNC,"unc",m->uncollided_secondaries,1},{ALEA_CELL_PARAM_BFLCL,"bflcl",m->magnetic_field,1}};for(size_t j=0;j<sizeof(params)/sizeof(params[0]);j++)if(m->parameter_flags&params[j].flag){XOK(alea_xml_writer_start_element(&x,"parameter"));XOK(alea_xml_writer_attribute(&x,"name",params[j].name));if(params[j].integer)XOK(alea_xml_writer_attribute_i(&x,"value",(int)params[j].value));else XOK(alea_xml_writer_attribute_f(&x,"value",params[j].value,17,false));XOK(alea_xml_writer_end_start_tag(&x,true));}XOK(alea_xml_writer_end_element(&x,"cell"));}else XOK(alea_xml_writer_end_start_tag(&x,true));arena_free(&a);}
    XOK(alea_xml_writer_end_element(&x,"cells"));XOK(alea_xml_writer_end_element(&x,"alea"));
    if(!alea_xml_writer_write(&x,stream)){alea_set_error_detail(ALEA_ERR_FILE_WRITE,"ALEA XML: write failed");return -1;}return 0;
}

int alea_xml_export_stream(const alea_model_t* model,FILE* stream){if(!model||!stream){alea_set_error_detail(ALEA_ERR_NULL_ARG,"alea_xml_export_stream: NULL argument");return -1;}return write_model(model,alea_model_system_const(model),stream);}
int alea_xml_export_system_stream(const alea_system_t* sys,FILE* stream){if(!sys||!stream){alea_set_error_detail(ALEA_ERR_NULL_ARG,"alea_xml_export_system_stream: NULL argument");return -1;}return write_model(NULL,sys,stream);}
int alea_xml_export(const alea_model_t* model,const char* filename){if(!model||!filename){alea_set_error_detail(ALEA_ERR_NULL_ARG,"alea_xml_export: NULL argument");return -1;}FILE* f=fopen(filename,"wb");if(!f){alea_set_error_detail(ALEA_ERR_FILE_WRITE,"cannot open %s",filename);return -1;}int rc=alea_xml_export_stream(model,f);if(fclose(f)&&!rc)rc=-1;return rc;}
int alea_xml_export_system(const alea_system_t* sys,const char* filename){if(!sys||!filename){alea_set_error_detail(ALEA_ERR_NULL_ARG,"alea_xml_export_system: NULL argument");return -1;}FILE* f=fopen(filename,"wb");if(!f){alea_set_error_detail(ALEA_ERR_FILE_WRITE,"cannot open %s",filename);return -1;}int rc=alea_xml_export_system_stream(sys,f);if(fclose(f)&&!rc)rc=-1;return rc;}
