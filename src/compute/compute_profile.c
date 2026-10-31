// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "compute_internal.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

_Static_assert(sizeof(compute_profile)==36,"GPU profile layout");
typedef struct { size_t index; uint32_t instructions; } profile_row;
static int compare_cost(const void* a,const void* b) {
    const profile_row* x=a;const profile_row* y=b;
    if(x->instructions!=y->instructions)return x->instructions<y->instructions?1:-1;
    return (x->index>y->index)-(x->index<y->index);
}
static int finish(FILE* f) { int failed=ferror(f);return fclose(f)!=0 || failed; }
int alea_compute_profile_write(const char* prefix,const compute_profile* rays,
    const compute_profile* costs,const compute_cell* cells,size_t cell_count,
    unsigned width,unsigned height,int grouped) {
    size_t n=(size_t)width*height;
    char* path=malloc(strlen(prefix)+32);
    profile_row* rows=malloc((cell_count?cell_count:1)*sizeof(*rows));
    unsigned char* rgb=malloc(n*3);
    if(!path || !rows || !rgb) {
        free(path);free(rows);free(rgb);
        alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,"Writing Metal profile");return -1;
    }
    sprintf(path,"%s-rays.csv",prefix);FILE* f=fopen(path,"w");if(!f)goto io;
    fprintf(f,"x,y,bvh_visits,cell_calls,instructions,primitive_intersections,interval_steps,hierarchy_steps,parent_checks,parent_cache_hits,rejected_ranges,csg_groups,instruction_unit\n");
    uint32_t max_work=0;
    for(size_t i=0;i<n;++i) {
        const compute_profile* r=&rays[i];
        fprintf(f,"%zu,%zu,%u,%u,%u,%u,%u,%u,%u,%u,%u,%d,%s\n",i%width,i/width,r->bvh,r->cells,r->instructions,r->primitives,r->intervals,r->hierarchy,r->parents,r->cache_hits,r->rejected,grouped,grouped?"bytecode_and_group_ops":"bytecode");
        if(r->instructions>max_work)max_work=r->instructions;
    }
    if(finish(f))goto io;
    for(size_t i=0;i<cell_count;++i)rows[i]=(profile_row){i,costs[i].instructions};
    qsort(rows,cell_count,sizeof(*rows),compare_cost);
    sprintf(path,"%s-cells.csv",prefix);f=fopen(path,"w");if(!f)goto io;
    fprintf(f,"cell_index,cell_id,material_id,program_instructions,cell_calls,instructions,primitive_intersections,interval_steps,is_container,csg_groups,instruction_unit\n");
    for(size_t j=0;j<cell_count;++j) {
        size_t i=rows[j].index;const compute_profile* r=&costs[i];
        if(!r->cells)continue;
        fprintf(f,"%zu,%d,%d,%u,%u,%u,%u,%u,%u,%d,%s\n",i,cells[i].cell_id,cells[i].material_id,cells[i].count,r->cells,r->instructions,r->primitives,r->intervals,cells[i].child!=COMPUTE_NONE || cells[i].lattice!=COMPUTE_NONE,grouped,grouped?"bytecode_and_group_ops":"bytecode");
    }
    if(finish(f))goto io;
    /* Logarithmic black-red-yellow-white heatmap of executed instructions. */
    for(size_t i=0;i<n;++i) {
        double v=max_work?log1p(rays[i].instructions)/log1p(max_work):0;
        for(int a=0;a<3;++a)rgb[3*i+a]=(unsigned char)(255*fmin(1,fmax(0,3*v-a)));
    }
    sprintf(path,"%s-heatmap.png",prefix);
    if(render_write_png(path,rgb,(int)width,(int)height))goto io;
    free(path);free(rows);free(rgb);return 0;
io:
    alea_set_error_detail(ALEA_ERR_INVALID_STATE,"Writing Metal profile '%s' failed",path);
    free(path);free(rows);free(rgb);return -1;
}
