// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* Host inventory of the existing Metal preparation format, not GPU execution. */
#include "alea_mcnp.h"
#include "compute/compute_prepare.h"
#include "util/compat.h"
#include <inttypes.h>
#include <stdio.h>
#include <string.h>
#ifndef _WIN32
#include <sys/resource.h>
#endif

static void peak_rss(void) {
#ifndef _WIN32
    struct rusage r;
    if (!getrusage(RUSAGE_SELF,&r)) {
#ifdef __APPLE__
        printf("%" PRIu64,(uint64_t)r.ru_maxrss);
#else
        printf("%" PRIu64,(uint64_t)r.ru_maxrss*1024);
#endif
        return;
    }
#endif
    printf("null");
}

int main(int argc,char **argv) {
    int accelerate=argc==3 && !strcmp(argv[2],"--accelerate");
    if(argc<2 || argc>3 || (argc==3 && !accelerate)) {
        fprintf(stderr,"Usage: %s model.i [--accelerate]\n",argv[0]);return 2;
    }
    double start=alea_monotonic_seconds();
    mcnp_model_t *model=mcnp_load(argv[1]);
    if(!model){fprintf(stderr,"Loading model: %s\n",alea_error());return 1;}
    double loaded=alea_monotonic_seconds();
    compute_scene s={0};
    if(alea_compute_pack_geometry(mcnp_model_system(model),&s,accelerate)) {
        fprintf(stderr,"Packing model: %s\n",alea_error());
        alea_compute_scene_free(&s);mcnp_model_destroy(model);return 1;
    }
    double packed=alea_monotonic_seconds();
    struct {const char *name;size_t count,stride;} rows[]={
        {"instructions",s.node_count,sizeof(*s.nodes)},
        {"primitives",s.primitive_count,sizeof(*s.primitives)},
        {"cells",s.cell_count,sizeof(*s.cells)},
        {"spatial_programs",s.spatial_count,sizeof(*s.spatial)},
        {"point_cells",s.point_cells?s.cell_count:0,sizeof(*s.point_cells)},
        {"point_lattices",s.point_lattice_count,sizeof(*s.point_lattices)},
        {"point_fills",s.point_fill_count,sizeof(*s.point_fills)},
        {"render_bvh",s.render.bvh_count,sizeof(*s.render.bvh)},
        {"transforms",s.render.transform_count,sizeof(*s.render.transforms)},
        {"render_lattices",s.render.lattice_count,sizeof(*s.render.lattices)},
        {"render_fills",s.render.fill_count,sizeof(*s.render.fills)}
    };
    printf("{\"schema\":1,\"execution\":\"host-preparation\","
           "\"layout\":\"existing-metal-wide\",\"accelerated\":%s,"
           "\"load_seconds\":%.9f,\"pack_seconds\":%.9f,\"peak_rss_bytes\":",
           accelerate?"true":"false",loaded-start,packed-loaded);
    peak_rss();printf(",\"buffers\":[");uint64_t total=0;
    for(size_t i=0;i<sizeof(rows)/sizeof(*rows);++i) {
        uint64_t bytes=(uint64_t)rows[i].count*rows[i].stride;total+=bytes;
        printf("%s{\"name\":\"%s\",\"count\":%zu,\"stride\":%zu,\"bytes\":%" PRIu64 "}",
               i?",":"",rows[i].name,rows[i].count,rows[i].stride,bytes);
    }
    printf("],\"buffer_payload_bytes\":%" PRIu64 "}\n",total);
    alea_compute_scene_free(&s);mcnp_model_destroy(model);return 0;
}
