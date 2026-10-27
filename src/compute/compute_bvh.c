// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "compute_internal.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Rebuild after conservative content refitting. Each universe retains its
 * original root offset: a binary tree with N leaves always has 2*N-1 nodes. */
typedef struct {compute_render_bvh leaf;float center[3],cost;} item;
typedef struct {float box[6];double cost;size_t count;} bin;
static void clear_box(float b[6]) {for(int a=0;a<3;++a){b[2*a]=INFINITY;b[2*a+1]=-INFINITY;}}
static void join(float b[6],const float x[6]) {
    for(int a=0;a<3;++a){b[2*a]=fminf(b[2*a],x[2*a]);b[2*a+1]=fmaxf(b[2*a+1],x[2*a+1]);}
}
static double area(const float b[6]) {
    double x=(double)b[1]-b[0],y=(double)b[3]-b[2],z=(double)b[5]-b[4];
    if(x<0 || y<0 || z<0)return 0;
    double value=2*(x*y+x*z+y*z);return isfinite(value)?value:INFINITY;
}
static int cmp_x(const void* a,const void* b) {float x=((const item*)a)->center[0],y=((const item*)b)->center[0];return (x>y)-(x<y);}
static int cmp_y(const void* a,const void* b) {float x=((const item*)a)->center[1],y=((const item*)b)->center[1];return (x>y)-(x<y);}
static int cmp_z(const void* a,const void* b) {float x=((const item*)a)->center[2],y=((const item*)b)->center[2];return (x>y)-(x<y);}
static unsigned bucket(float center,double low,double extent) {
    double k=16*((double)center-low)/extent;
    return k<=0?0:k>=15?15:(unsigned)k;
}
static size_t partition(item* v,size_t count,unsigned depth) {
    float low[3],high[3];for(int a=0;a<3;++a)low[a]=high[a]=v[0].center[a];
    for(size_t i=1;i<count;++i) {
        if((i&1023)==0 && alea_interrupted())return 0;
        for(int a=0;a<3;++a){low[a]=fminf(low[a],v[i].center[a]);high[a]=fmaxf(high[a],v[i].center[a]);}
    }
    double best=INFINITY;unsigned axis=0,split=0;int found=0;
    for(unsigned a=1;a<3;++a)if((double)high[a]-low[a]>(double)high[axis]-low[axis])axis=a;
    for(unsigned a=0;a<3;++a) {
        double extent=(double)high[a]-low[a];if(!(extent>0) || depth>=32)continue;
        bin bins[16]={0};for(int j=0;j<16;++j)clear_box(bins[j].box);
        for(size_t i=0;i<count;++i) {
            if((i&1023)==1023 && alea_interrupted())return 0;
            bin* b=&bins[bucket(v[i].center[a],low[a],extent)];join(b->box,v[i].leaf.bounds);b->cost+=v[i].cost;++b->count;
        }
        float left[6],right[16][6];clear_box(left);float box[6];clear_box(box);
        double costs[16],cost=0;size_t sizes[16],size=0;
        for(int j=15;j>=0;--j){join(box,bins[j].box);memcpy(right[j],box,sizeof(box));cost+=bins[j].cost;size+=bins[j].count;costs[j]=cost;sizes[j]=size;}
        cost=0;size=0;
        for(unsigned j=0;j<15;++j) {
            join(left,bins[j].box);cost+=bins[j].cost;size+=bins[j].count;
            if(!size || !sizes[j+1])continue;
            double score=area(left)*cost+area(right[j+1])*costs[j+1];
            if(score<best){best=score;axis=a;split=j;found=1;}
        }
    }
    if(found) {
        size_t lo=0,hi=count,visited=0;double extent=(double)high[axis]-low[axis];
        while(lo<hi) {
            if((visited++&1023)==1023 && alea_interrupted())return 0;
            if(bucket(v[lo].center[axis],low[axis],extent)<=split)++lo;
            else {item swap=v[lo];v[lo]=v[--hi];v[hi]=swap;}
        }
        if(lo && lo<count)return lo;
    }
    qsort(v,count,sizeof(*v),axis==0?cmp_x:axis==1?cmp_y:cmp_z);
    return count/2;
}
static int build(item* v,size_t count,unsigned depth,compute_render_bvh* out,size_t* cursor) {
    if((*cursor&1023)==0 && alea_interrupted())return -1;
    size_t id=(*cursor)++;
    if(count==1){out[id]=v[0].leaf;out[id].escape=(uint32_t)*cursor;return 0;}
    size_t middle=partition(v,count,depth);
    if(!middle || build(v,middle,depth+1,out,cursor))return -1;
    size_t right=*cursor;
    if(build(v+middle,count-middle,depth+1,out,cursor))return -1;
    compute_render_bvh* n=&out[id];*n=out[id+1];n->cell=COMPUTE_NONE;n->escape=(uint32_t)*cursor;
    join(n->bounds,out[right].bounds);join(n->fallback,out[right].fallback);n->guard=fminf(n->guard,out[right].guard);
    return 0;
}
int alea_compute_bvh_rebuild(compute_scene* s) {
    const char* mode=getenv("ALEA_METAL_BVH");
    if(!mode || !strcmp(mode,"0"))return 0;
    int weighted=!strcmp(mode,"weighted");
    if(!weighted && strcmp(mode,"sah")){alea_set_error_detail(ALEA_ERR_INVALID_ARG,"ALEA_METAL_BVH must be 0, sah, or weighted");return -1;}
    compute_render_scene* r=&s->render;if(!r->bvh_count)return 0;
    item* items=malloc(((r->bvh_count+1)/2)*sizeof(*items));
    compute_render_bvh* out=malloc(r->bvh_count*sizeof(*out));
    if(!items || !out){free(items);free(out);alea_set_error_detail(ALEA_ERR_OUT_OF_MEMORY,"Rebuilding render BVHs");return -1;}
    for(size_t root=0;root<r->bvh_count;) {
        if(alea_interrupted())goto interrupted;
        size_t end=r->bvh[root].escape,count=0;
        for(size_t i=root;i<end;++i) {
            if(((i-root)&1023)==0 && alea_interrupted())goto interrupted;
            if(r->bvh[i].cell==COMPUTE_NONE)continue;
            item* v=&items[count++];v->leaf=r->bvh[i];compute_cell c=s->cells[v->leaf.cell];
            v->cost=weighted?fmaxf(1,fminf((float)c.count,262144)):1;
            for(int a=0;a<3;++a){double center=.5*((double)v->leaf.bounds[2*a]+v->leaf.bounds[2*a+1]);v->center[a]=isfinite(center)?(float)center:0;}
        }
        if(!count){free(items);free(out);alea_set_error_detail(ALEA_ERR_INVALID_STATE,"Empty render BVH during rebuild");return -1;}
        size_t cursor=root;if(build(items,count,0,out,&cursor))goto interrupted;
        root=end;
    }
    free(items);free(r->bvh);r->bvh=out;return 0;
interrupted:
    free(items);free(out);alea_set_error_detail(ALEA_ERR_INTERRUPTED,"Rebuilding render BVHs interrupted");return -1;
}
