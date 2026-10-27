// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#define _POSIX_C_SOURCE 200809L
#include "alea.h"
#include "alea_mcnp.h"
#include "alea_compute.h"
#include "core/alea_cell.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include <signal.h>
#include <string.h>
#include <sys/stat.h>
#include <inttypes.h>

static double seconds(void) {struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+t.tv_nsec*1e-9;}
static void interrupted(int sig) {(void)sig;alea_interrupt();}
static int integer(const char* s,int max) {
    char* end;errno=0;long value=strtol(s,&end,10);
    return errno || end==s || *end || value<1 || value>max?-1:(int)value;
}
typedef struct {double start,first,last,gap;size_t done,calls;} bench_progress;
static int double_order(const void* a,const void* b) {double x=*(const double*)a,y=*(const double*)b;return (x>y)-(x<y);}
static int report_progress(void* user,const render_framebuffer_t* fb,size_t done,size_t total) {
    (void)fb;(void)total;bench_progress* p=user;double now=seconds();
    if(done && !p->done)p->first=now-p->start;
    if(p->calls)p->gap=fmax(p->gap,now-p->last);
    p->last=now;p->done=done;++p->calls;return 0;
}
/* Small fixed workloads also exercise more pixels than the resident pool. */
static alea_system_t* synthetic_scene(const char* name) {
    alea_system_t* sys=alea_create();if(!sys)return NULL;int mat=alea_add_material(sys,7);
    int mixed=!strcmp(name,"synthetic:mixed"),reject=!strcmp(name,"synthetic:reject");
    alea_node_id_t parts=ALEA_NODE_ID_INVALID;
    if(mixed) {
        for(int i=0;i<256;++i) {
            int sphere=alea_sphere_surface(sys,0,-8+(i%4)*.9,(i/4%8-3.5)*.9,(i/32-3.5)*.9,.32);
            alea_node_id_t part=alea_halfspace(sys,sphere,-1);parts=parts==ALEA_NODE_ID_INVALID?part:alea_union(sys,parts,part);
        }
        alea_add_cell(sys,42,parts,mat,-1,0);
        for(int i=0;i<64;++i) {
            int sphere=alea_sphere_surface(sys,0,3+(i%4)*1.5,(i/4%4-1.5)*1.5,(i/16-1.5)*1.5,.55);
            alea_add_cell(sys,100+i,alea_halfspace(sys,sphere,-1),mat,-1,0);
        }
    } else {
        int box=alea_box_surface(sys,0,-12,12,reject?0:-12,reject?.1:12,-12,12);
        alea_node_id_t parent=alea_halfspace(sys,box,-1);
        if(!reject)for(int i=0;i<128;++i) {
            int sphere=alea_sphere_surface(sys,0,(i%8-3.5)*2.5,(i/8%4-1.5)*3.5,(i/32-1.5)*3.5,.6);
            parent=alea_difference(sys,parent,alea_halfspace(sys,sphere,-1));
        }
        int cell=alea_add_cell(sys,10,parent,0,0,0);alea_set_cell_fill(sys,cell,1,0);
        for(int i=0;i<(reject?128:8);++i) {
            int sphere=alea_sphere_surface(sys,0,reject?(i%8-3.5)*2.7:(i&1?7:-7),reject?(i/8%4-1)*4:(i&2?7:-7),reject?(i/32-1.5)*3:(i&4?7:-7),.8);
            alea_add_cell(sys,100+i,alea_halfspace(sys,sphere,-1),mat,-1,1);
        }
    }
    if(alea_prepare_query_acceleration(sys)){alea_destroy(sys);return NULL;}return sys;
}
static uint64_t executable_identity(const char* path) {
    FILE* f=fopen(path,"rb");if(!f)return 0;
    uint64_t hash=UINT64_C(14695981039346656037);unsigned char bytes[4096];size_t n;
    while((n=fread(bytes,1,sizeof(bytes),f)))for(size_t i=0;i<n;++i){hash^=bytes[i];hash*=UINT64_C(1099511628211);}
    int failed=ferror(f);fclose(f);return failed?0:hash;
}
int main(int argc,char**argv) {
    if(argc!=3 && argc!=5 && argc!=6 && argc!=7) {
        fprintf(stderr,"Usage: %s MODEL auto|bounds|current|lazy|sah|weighted|refill|balanced|all [WIDTH HEIGHT [FRAMES [VIEW]]]\n",argv[0]);
        fprintf(stderr,"ITER cameras: VIEW 1=cutaway (default), 2=opposite cutaway, 3=exterior, 4=component close-up. Default 16x12, two frames. Prefix a mode with p- for progressive comparisons. ALEA_BENCH_PROFILE enables instrumented CSV output.\n");return 1;
    }
    int width=argc>=5?integer(argv[3],4096):16,height=argc>=5?integer(argv[4],4096):12;
    int frames=argc>=6?integer(argv[5],100):2,view=argc==7?integer(argv[6],4):1;
    if(width<0 || height<0 || frames<0 || view<0)return 1;
    int progressive=!strncmp(argv[2],"p-",2);const char* mode_name=argv[2]+(progressive?2:0);
    int experimental=!strcmp(mode_name,"current") || !strcmp(mode_name,"lazy") || !strcmp(mode_name,"sah") ||
        !strcmp(mode_name,"weighted") || !strcmp(mode_name,"refill") || !strcmp(mode_name,"balanced") || !strcmp(mode_name,"all");
    int automatic=!strcmp(mode_name,"auto"),bounds_only=!strcmp(mode_name,"bounds");
    if(!experimental && !automatic && !bounds_only){fprintf(stderr,"Unknown benchmark mode: %s\n",mode_name);return 1;}
    const char* label=experimental?mode_name:automatic?"Automatic":"Bounds";
    setbuf(stdout,NULL);signal(SIGINT,interrupted);signal(SIGTERM,interrupted);
    unsetenv("ALEA_METAL_PROFILE");setenv("ALEA_METAL_SPATIAL",experimental?"1":"0",1);
    const char* flags[]={"ALEA_METAL_LAZY","ALEA_METAL_BVH","ALEA_METAL_REFILL"};
    for(int i=0;i<3;++i)setenv(flags[i],"0",1);
    setenv("ALEA_BENCH_TRACE","1",1);
    printf("Executable fingerprint fnv1a64=%016" PRIx64 "\n",executable_identity(argv[0]));
    struct stat info;if(!stat(argv[1],&info))printf("Model %s (%lld bytes, mtime %lld)\n",argv[1],(long long)info.st_size,(long long)info.st_mtime);
    printf("Build %s %s, run identity %s, mode %s, image %dx%d, view %d, AA1, progressive %d\n",__DATE__,__TIME__,getenv("ALEA_BENCH_ID")?getenv("ALEA_BENCH_ID"):"unspecified",mode_name,width,height,view,progressive);
    int synthetic=!strcmp(argv[1],"synthetic:sparse") || !strcmp(argv[1],"synthetic:reject") || !strcmp(argv[1],"synthetic:mixed");
    double start=seconds();mcnp_model_t* model=synthetic?NULL:mcnp_load(argv[1]);
    alea_system_t* sys=synthetic?synthetic_scene(argv[1]):model?model->sys:NULL;
    if(!sys){fprintf(stderr,"%s\n",alea_error());return 1;}
    printf("Load %.3f s\n",seconds()-start);
    render_config_t cfg;render_config_init(&cfg);cfg.width=width;cfg.height=height;cfg.log_level=0;
    cfg.eye_set=cfg.target_set=1;cfg.eye[0]=1800;cfg.eye[1]=-2200;cfg.eye[2]=1400;
    cfg.target[0]=cfg.target[1]=cfg.target[2]=0;
    cfg.num_clips=1;cfg.clips[0]=(render_clip_plane_t){{0,1,0},0};
    if(view==2){cfg.eye[0]=-1800;cfg.eye[1]=2200;cfg.eye[2]=1000;cfg.clips[0].normal[1]=-1;}
    if(view==3)cfg.num_clips=0;
    if(view==4){cfg.eye[0]=300;cfg.eye[1]=80;cfg.eye[2]=-120;cfg.target[0]=300;cfg.target[1]=300;cfg.target[2]=-200;cfg.num_clips=0;}
    if(synthetic){cfg.eye[0]=0;cfg.eye[1]=-35;cfg.eye[2]=12;cfg.target[0]=cfg.target[1]=cfg.target[2]=0;cfg.num_clips=0;}
    printf("Camera eye %.9g,%.9g,%.9g target %.9g,%.9g,%.9g, clips %d\n",cfg.eye[0],cfg.eye[1],cfg.eye[2],cfg.target[0],cfg.target[1],cfg.target[2],cfg.num_clips);
    render_camera_t camera;
    render_framebuffer_t* image[2]={NULL,NULL};alea_compute_t* gpu[2]={NULL,NULL};double times[2][100];int status=1;
    if(render_camera_setup(&camera,&cfg,sys))goto done;
    for(int mode=0;mode<2;++mode) {
        if(mode) {
            if(experimental) {
                int all=!strcmp(mode_name,"all");
                if(all || !strcmp(mode_name,"balanced") || !strcmp(mode_name,"lazy"))setenv(flags[0],"1",1);
                if(all || !strcmp(mode_name,"weighted"))setenv(flags[1],"weighted",1);
                else if(!strcmp(mode_name,"sah"))setenv(flags[1],"sah",1);
                if(all || !strcmp(mode_name,"balanced") || !strcmp(mode_name,"refill"))setenv(flags[2],"1",1);
            } else if(automatic || bounds_only){if(setenv("ALEA_METAL_SPATIAL",automatic?"1":"bounds",1))goto done;}
        }
        const char* prefix=getenv("ALEA_BENCH_PROFILE");
        if(prefix) {
            char path[4096];int length=snprintf(path,sizeof(path),"%s-%s",prefix,mode?"accelerated":"baseline");
            if(length<0 || (size_t)length>=sizeof(path) || setenv("ALEA_METAL_PROFILE",path,1))goto done;
        }
        image[mode]=render_framebuffer_create(width,height,1);if(!image[mode])goto done;
        start=seconds();gpu[mode]=alea_compute_create(sys,ALEA_COMPUTE_METAL,ALEA_COMPUTE_FP32);
        if(!gpu[mode])goto done;
        printf("%s prepare %.3f s\n",mode?label:"Baseline",seconds()-start);
        printf("%s flags: spatial=%s lazy=%s bvh=%s refill=%s\n",mode?label:"Baseline",getenv("ALEA_METAL_SPATIAL"),getenv(flags[0]),getenv(flags[1]),getenv(flags[2]));
    }
    for(int frame=-1;frame<frames;++frame)for(int order=0;order<2;++order) {
        int mode=frame>=0?(order+(frame&1))%2:order;
        bench_progress progress={.start=seconds()};start=progress.start;
        if(progressive?alea_compute_render_progressive(gpu[mode],&cfg,&camera,image[mode],report_progress,&progress):alea_compute_render(gpu[mode],&cfg,&camera,image[mode]))goto done;
        double elapsed=seconds()-start;
        printf("%s %s %d %.3f s",mode?label:"Baseline",frame<0?"warmup":"frame",frame,elapsed);
        if(progressive)printf(", first %.3f s, callback gap %.3f s, %zu callbacks",progress.first,progress.gap,progress.calls);
        putchar('\n');if(frame>=0)times[mode][frame]=elapsed;
    }
    for(int mode=0;mode<2;++mode) {
        qsort(times[mode],frames,sizeof(double),double_order);
        double median=(times[mode][frames/2]+times[mode][(frames-1)/2])*.5;
        printf("%s median %.3f s, range %.3f..%.3f s (%d warm frames)\n",mode?label:"Baseline",median,times[mode][0],times[mode][frames-1],frames);
    }
    size_t differences=0,depth_differences=0,hits[2]={0,0};double depth=0,color=0,normal=0;
    for(size_t i=0;i<(size_t)width*height;++i) {
        for(int mode=0;mode<2;++mode)hits[mode]+=image[mode]->cell_id[i]>=0;
        if(image[0]->cell_id[i]!=image[1]->cell_id[i] || image[0]->material_id[i]!=image[1]->material_id[i]) {
            if(differences<12)printf("Pixel %zu baseline cell %d depth %.9g, experiment cell %d depth %.9g\n",i,
                image[0]->cell_id[i],image[0]->depth[i],image[1]->cell_id[i],image[1]->depth[i]);
            ++differences;
        } else if(image[0]->cell_id[i]>=0) {
            double delta=fabs(image[0]->depth[i]-image[1]->depth[i]);depth=fmax(depth,delta);depth_differences+=delta!=0;
        }
        for(int k=0;k<3;++k) {
            color=fmax(color,fabs(image[0]->color[3*i+k]-image[1]->color[3*i+k]));
            normal=fmax(normal,fabs(image[0]->normal[3*i+k]-image[1]->normal[3*i+k]));
        }
    }
    printf("Compared %d pixels: ID differences %zu, changed matching depths %zu, max depth %.9g, RGB %.9g, normal %.9g\n",
        width*height,differences,depth_differences,depth,color,normal);
    printf("Hit pixels: reference %zu, experiment %zu\n",hits[0],hits[1]);
    puts("Agreement with the baseline is a regression check, not a CPU/analytic accuracy proof.");
    status=differences || depth || color || normal || (synthetic && !hits[0])?2:0;
done:
    if(status==1)fprintf(stderr,"%s\n",alea_error());
    alea_compute_destroy(gpu[0]);alea_compute_destroy(gpu[1]);render_framebuffer_free(image[0]);render_framebuffer_free(image[1]);
    render_config_free(&cfg);if(model)mcnp_model_destroy(model);else alea_destroy(sys);return status;
}
