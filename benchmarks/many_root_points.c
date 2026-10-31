// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#define _POSIX_C_SOURCE 200809L
#include "alea.h"
#include "alea_compute.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
static double seconds(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec*1e-9;}
static int order(const void*a,const void*b){double x=*(const double*)a,y=*(const double*)b;return (x>y)-(x<y);}
/* Independently constructed public fixture: a complete regular grid of boxes.
 * Pixel centers avoid faces; floor(x),floor(y) directly determine ownership.
 * Scanner and coarse candidates receive the same 16x16-block input order.
 * The overlap case places every box over the full domain: classification must
 * return 1,024 owners, whereas selection must return the first cell.
 * The hierarchy case translates the same grid into one filled universe, so
 * root filtering cannot remove its child scan. Analytic ownership stays known.
 * Timings include staging and oracle checks; no finer packet filter is retained. */
int main(int argc,char**argv) {
    int broad=argc==2 && !strcmp(argv[1],"overlap");
    int hierarchy=argc==2 && !strcmp(argv[1],"hierarchy");
    if(argc>2 || (argc==2 && !broad && !hierarchy)){fprintf(stderr,"Usage: %s [overlap|hierarchy]\n",argv[0]);return 1;}
    alea_system_t* s=alea_create();if(!s)return 1;
    int material=alea_add_material(s,7);if(material<0)return 1;
    if(hierarchy) {
        int box=alea_box_surface(s,0,1000,1032,0,32,-1,1);
        int parent=alea_add_cell(s,100000,alea_halfspace(s,box,-1),ALEA_MATERIAL_VOID,0,0);
        const double translation[]={1000,0,0};
        if(parent<0 || alea_add_transform(s,1,translation,3,0) || alea_set_fill(s,parent,10,1))return 1;
    }
    for(int y=0;y<32;++y)for(int x=0;x<32;++x) {
        int box=alea_box_surface(s,0,broad?0:x,broad?32:x+1,broad?0:y,broad?32:y+1,-1,1);
        if(box<0 || alea_add_cell(s,y*32+x+1,alea_halfspace(s,box,-1),material,-1,hierarchy?10:0)<0)return 1;
    }
    for(int mode=0;mode<2;++mode) {
        if(setenv("ALEA_METAL_POINT_CANDIDATES",mode?"1":"0",1))return 1;
        alea_compute_t* ctx=alea_compute_create(s,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);
        if(!ctx){fprintf(stderr,"%s\n",alea_error());return 2;}
        for(int selected=0;selected<2;++selected) {
            double times[5],first=0;
            for(int pass=0;pass<6;++pass) {
                double start=seconds();
                for(int tile=0;tile<4;++tile) {
                    double xyz[4096*3];alea_compute_result_t all[4096];alea_compute_selection_t one[4096];
                    for(int i=0;i<4096;++i) {
                        int x=tile%2*64+i/256%4*16+i%16,y=tile/2*64+i/1024*16+i%256/16;
                        xyz[3*i]=(x+.5)/4+(hierarchy?1000:0);xyz[3*i+1]=(y+.5)/4;xyz[3*i+2]=0;
                    }
                    if(selected?alea_compute_select(ctx,xyz,4096,one):alea_compute_classify(ctx,xyz,4096,all)) {
                        fprintf(stderr,"%s\n",alea_error());return 3;
                    }
                    for(int i=0;i<4096;++i) {
                        int expected=broad?1:(int)xyz[3*i+1]*32+(int)(xyz[3*i]-(hierarchy?1000:0))+1;
                        if(selected?(one[i].cell_id!=expected || one[i].material_id!=7 || one[i].depth!=(hierarchy?1:0) || one[i].status!=ALEA_COMPUTE_FOUND):
                            (all[i].cell_id!=(broad?-1:expected) || all[i].material_id!=(broad?-1:7) || all[i].owner_count!=(broad?1024u:1u)))return 4;
                    }
                }
                double elapsed=seconds()-start;
                if(pass)times[pass-1]=elapsed;else first=elapsed;
            }
            qsort(times,5,sizeof(double),order);
            printf("fixture=%s candidates=%d operation=%s first_s=%.6f warm_median_s=%.6f min_s=%.6f max_s=%.6f oracle=PASS\n",broad?"overlap":hierarchy?"hierarchy":"grid",mode,selected?"selected":"all",first,times[2],times[0],times[4]);
        }
        alea_compute_destroy(ctx);
    }
    alea_destroy(s);return 0;
}
