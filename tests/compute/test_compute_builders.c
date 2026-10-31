// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "alea_test.h"
#include "alea.h"
#include "compute/compute_internal.h"
/* Deterministic interruption at internal checkpoints, without sleeps or
 * production failure hooks. These implementations replace archive symbols. */
static size_t checkpoints,interrupt_at;
static int builder_interrupted(void) {++checkpoints;return interrupt_at && checkpoints>=interrupt_at;}
#define alea_interrupted builder_interrupted
#include "../../src/compute/compute_bvh.c"
#include "../../src/compute/compute_point.c"
#include "../../src/compute/compute_groups.c"
#undef alea_interrupted

static alea_system_t* grid_scene(void) {
    alea_system_t* sys=alea_create();int mat=alea_add_material(sys,1);
    for(int i=0;i<2048;++i) {
        int sphere=alea_sphere_surface(sys,0,3*(i%16),3*(i/16%16),3*(i/256),.5);
        alea_add_cell(sys,i+1,alea_halfspace(sys,sphere,-1),mat,-1,0);
    }
    return sys;
}
TEST(bvh_interruption_preserves_scene_and_allows_retry) {
    alea_system_t* sys=grid_scene();
    ASSERT_EQ(alea_setenv("ALEA_METAL_BVH","weighted",1),0);
    for(unsigned stage=0;stage<3;++stage) {
        compute_scene s;ASSERT_EQ(alea_compute_pack_geometry(sys,&s,0),0);
        size_t bytes=s.render.bvh_count*sizeof(*s.render.bvh);
        compute_render_bvh* original=s.render.bvh,*copy=malloc(bytes);ASSERT_NOT_NULL(copy);memcpy(copy,original,bytes);
        checkpoints=0;interrupt_at=stage==0?2:stage==1?8:20;
        int result=alea_compute_bvh_rebuild(&s);interrupt_at=0;
        ASSERT_EQ(result,-1);ASSERT_EQ(alea_error_code(),ALEA_ERR_INTERRUPTED);
        ASSERT_TRUE(s.render.bvh==original);ASSERT_EQ(memcmp(copy,original,bytes),0);
        ASSERT_EQ(alea_compute_bvh_rebuild(&s),0);
        ASSERT_EQ(s.render.bvh[s.render.root].escape,s.render.bvh_count);
        free(copy);alea_compute_scene_free(&s);
    }
    ASSERT_EQ(alea_unsetenv("ALEA_METAL_BVH"),0);alea_destroy(sys);
}

TEST(point_program_certification_preserves_errors_and_retries) {
    compute_primitive primitive={.type=COMPUTE_SPHERE,.p={0,0,0,1e20}};
    compute_node node={COMPUTE_LEAF,0};
    compute_cell cell={.root=0,.count=1,.spatial=0,.lattice=COMPUTE_NONE};
    compute_spatial_program region={.bounds={-1,1,-1,1,-1,1},.child=COMPUTE_NONE};
    metal_point_cell link={.spatial=COMPUTE_NONE};
    compute_scene s={.cells=&cell,.cell_count=1,.primitives=&primitive,.nodes=&node,
        .point_cells=&link,.spatial=&region,.spatial_count=1};
    /* The original sphere overflows even in this small region. */
    ASSERT_EQ(alea_compute_point_programs_pack(&s),0);ASSERT_EQ(link.spatial,COMPUTE_NONE);
    primitive.p[3]=1;
    checkpoints=0;interrupt_at=1;
    ASSERT_EQ(alea_compute_point_programs_pack(&s),-1);ASSERT_EQ(link.spatial,COMPUTE_NONE);
    interrupt_at=0;
    ASSERT_EQ(alea_compute_point_programs_pack(&s),0);ASSERT_EQ(link.spatial,0);
    ASSERT_EQ(alea_setenv("ALEA_METAL_POINT_SPATIAL","invalid",1),0);
    ASSERT_EQ(alea_compute_point_programs_pack(&s),-1);ASSERT_EQ(alea_error_code(),ALEA_ERR_INVALID_ARG);
    ASSERT_EQ(alea_unsetenv("ALEA_METAL_POINT_SPATIAL"),0);
}
TEST(point_candidates_keep_unknown_boundaries_and_retry) {
    alea_system_t* sys=grid_scene();compute_scene s;
    ASSERT_EQ(alea_compute_pack_geometry(sys,&s,0),0);alea_destroy(sys);
    int8_t* signs=malloc(s.primitive_count);uint32_t* next=malloc(s.cell_count*sizeof(*next));
    ASSERT_NOT_NULL(signs);ASSERT_NOT_NULL(next);
    float bounds[]={-.1,.1,-.1,.1,-.1,.1};uint32_t root;compute_candidate_stats stats={0};
    ASSERT_EQ(alea_compute_point_candidates(&s,bounds,signs,next,&root,&stats),0);
    ASSERT_EQ(stats.cells,2048);ASSERT_EQ(stats.rejected,2047);ASSERT_EQ(root,0);ASSERT_EQ(next[0],COMPUTE_NONE);
    bounds[0]=bounds[1]=.5f;stats=(compute_candidate_stats){0};
    ASSERT_EQ(alea_compute_point_candidates(&s,bounds,signs,next,&root,&stats),0);ASSERT_EQ(root,0);
    bounds[0]=bounds[1]=1e30f;stats=(compute_candidate_stats){0};
    ASSERT_EQ(alea_compute_point_candidates(&s,bounds,signs,next,&root,&stats),0);ASSERT_EQ(stats.rejected,0);
    checkpoints=0;interrupt_at=2;
    ASSERT_EQ(alea_compute_point_candidates(&s,bounds,signs,next,&root,&stats),-1);
    ASSERT_EQ(alea_error_code(),ALEA_ERR_INTERRUPTED);interrupt_at=0;
    bounds[0]=-.1f;bounds[1]=.1f;stats=(compute_candidate_stats){0};
    ASSERT_EQ(alea_compute_point_candidates(&s,bounds,signs,next,&root,&stats),0);ASSERT_EQ(stats.rejected,2047);
    free(signs);free(next);alea_compute_scene_free(&s);
}
TEST(point_candidates_preserve_nonfinite_operands_and_budget_fallback) {
    compute_primitive primitives[2]={
        {.type=COMPUTE_SPHERE,.p={0,0,0,1e20}},
        {.type=COMPUTE_BOX,.p={10,11,10,11,10,11}}};
    compute_node nodes[]={{COMPUTE_LEAF,0},{COMPUTE_SKIP_INTERSECTION,2},
                          {COMPUTE_LEAF,2},{COMPUTE_INTERSECTION,0}};
    compute_cell cell={.root=0,.count=4,.lattice=COMPUTE_NONE};
    compute_scene scene={.nodes=nodes,.primitives=primitives,.primitive_count=2,.cells=&cell,.cell_count=1};
    const float bounds[]={-1,1,-1,1,-1,1};int8_t signs[2];uint32_t next,root;
    compute_candidate_stats stats={0};
    ASSERT_EQ(alea_compute_point_candidates(&scene,bounds,signs,&next,&root,&stats),0);
    ASSERT_EQ(root,0);ASSERT_EQ(stats.rejected,0);
    primitives[0].p[3]=2;stats=(compute_candidate_stats){0};
    ASSERT_EQ(alea_compute_point_candidates(&scene,bounds,signs,&next,&root,&stats),0);
    ASSERT_EQ(root,COMPUTE_NONE);ASSERT_EQ(stats.rejected,1);
    cell.count=POINT_PROOF_BUDGET+1;stats=(compute_candidate_stats){0};
    ASSERT_EQ(alea_compute_point_candidates(&scene,bounds,signs,&next,&root,&stats),0);
    ASSERT_EQ(root,0);ASSERT_EQ(stats.instructions,0);ASSERT_EQ(stats.budget_fallbacks,1);
}
TEST(compact_code_exact_in_place_and_wide_fallback) {
    compute_node original[]={{COMPUTE_LEAF,0x0fffffffu},{COMPUTE_SKIP_UNION,12345},
        {COMPUTE_CONSTANT,1},{COMPUTE_COMPLEMENT,0},{15,0}};
    enum {N=sizeof(original)/sizeof(*original)};
    compute_node nodes[N];memcpy(nodes,original,sizeof(nodes));
    uint32_t words[N];
    ASSERT_EQ(alea_compute_code_compact(nodes,N,words),1);
    compute_scene scene={.compact_nodes=words};
    for(size_t i=0;i<N;++i) {
        compute_node decoded=alea_compute_node_at(&scene,i);
        ASSERT_EQ(decoded.op,original[i].op);ASSERT_EQ(decoded.arg,original[i].arg);
    }
    ASSERT_EQ(memcmp(nodes,original,sizeof(nodes)),0);
    ASSERT_EQ(alea_compute_code_compact(nodes,N,(uint32_t*)nodes),1);
    ASSERT_EQ(memcmp(nodes,words,sizeof(words)),0);
    for(unsigned invalid=0;invalid<2;++invalid) {
        memcpy(nodes,original,sizeof(nodes));
        if(invalid)nodes[N-1].op=16;else nodes[N-1].arg=0x10000000u;
        compute_node copy[N];memcpy(copy,nodes,sizeof(copy));
        ASSERT_EQ(alea_compute_code_compact(nodes,N,(uint32_t*)nodes),0);
        ASSERT_EQ(memcmp(nodes,copy,sizeof(nodes)),0);
    }
    memcpy(nodes,original,sizeof(nodes));alea_interrupt();
    int result=alea_compute_code_compact(nodes,N,(uint32_t*)nodes);
    alea_clear_interrupt();
    ASSERT_EQ(result,-1);ASSERT_EQ(alea_error_code(),ALEA_ERR_INTERRUPTED);
    ASSERT_EQ(memcmp(nodes,original,sizeof(nodes)),0);
    ASSERT_EQ(alea_compute_code_compact(nodes,N,(uint32_t*)nodes),1);
    ASSERT_EQ(memcmp(nodes,words,sizeof(words)),0);
    ASSERT_EQ(alea_compute_code_compact(NULL,0,NULL),1);
}

typedef struct { float value;unsigned trace,error; } jump_result;
static jump_result run_jump_test(const compute_node* code,size_t count,const float* input) {
    float stack[COMPUTE_STACK];unsigned top=0; jump_result r={0};
    for(size_t pc=0;pc<count;++pc) {
        compute_node n=code[pc];
        if(n.op==COMPUTE_LEAF) {
            r.trace=r.trace*7+n.arg+1;float v=input[n.arg];
            if(!isfinite(v)){r.error=n.arg+1;return r;}stack[top++]=v;
        } else if(n.op==COMPUTE_SKIP_UNION || n.op==COMPUTE_SKIP_INTERSECTION) {
            if(n.op==COMPUTE_SKIP_UNION?stack[top-1]<=0:stack[top-1]>0)pc+=n.arg;
        } else {
            float right=stack[--top],left=stack[top-1];
            stack[top-1]=n.op==COMPUTE_UNION?fminf(left,right):fmaxf(left,right);
        }
    }
    r.value=stack[0];return r;
}
TEST(jump_threading_preserves_trace_errors_boundaries_and_encoding) {
    float values[]={-1,0,-0.0f,1,NAN,INFINITY};
    for(int op=0;op<2;++op) {
        compute_node original[]={{COMPUTE_LEAF,0},{op?COMPUTE_SKIP_UNION:COMPUTE_SKIP_INTERSECTION,2},
            {COMPUTE_LEAF,1},{op?COMPUTE_UNION:COMPUTE_INTERSECTION,0},
            {op?COMPUTE_SKIP_UNION:COMPUTE_SKIP_INTERSECTION,2},{COMPUTE_LEAF,2},
            {op?COMPUTE_UNION:COMPUTE_INTERSECTION,0}};
        compute_node optimized[7];memcpy(optimized,original,sizeof(original));
        compute_cell cells[2]={{.root=0,.count=7},{.root=0,.count=7}};
        compute_scene s={.nodes=optimized,.node_count=7,.cells=cells,.cell_count=2};size_t rewritten;
        ASSERT_EQ(alea_compute_thread_jumps(&s,&rewritten),0);ASSERT_EQ(rewritten,1);ASSERT_EQ(optimized[1].arg,5);
        ASSERT_EQ(alea_compute_thread_jumps(&s,&rewritten),0);ASSERT_EQ(rewritten,0);
        uint32_t compact[7];ASSERT_EQ(alea_compute_code_compact(optimized,7,compact),1);
        compute_scene packed={.compact_nodes=compact};
        for(int i=0;i<7;++i){compute_node n=alea_compute_node_at(&packed,i);ASSERT_EQ(n.op,optimized[i].op);ASSERT_EQ(n.arg,optimized[i].arg);}
        for(unsigned a=0;a<6;++a)for(unsigned b=0;b<6;++b)for(unsigned c=0;c<6;++c) {
            float in[]={values[a],values[b],values[c]};jump_result x=run_jump_test(original,7,in),y=run_jump_test(optimized,7,in);
            ASSERT_EQ(x.trace,y.trace);ASSERT_EQ(x.error,y.error);ASSERT_EQ(memcmp(&x.value,&y.value,sizeof(float)),0);
        }
    }
}
TEST(jump_threading_rejects_invalid_ranges_and_skips_aliases) {
    compute_node code[]={{COMPUTE_LEAF,0},{COMPUTE_SKIP_INTERSECTION,2},{COMPUTE_LEAF,1},
        {COMPUTE_INTERSECTION,0},{COMPUTE_SKIP_INTERSECTION,2},{COMPUTE_LEAF,2},{COMPUTE_INTERSECTION,0}};
    compute_node copy[7];memcpy(copy,code,sizeof(code));
    compute_cell cells[2]={{.root=0,.count=7},{.root=0,.count=4}};
    compute_scene s={.nodes=code,.node_count=7,.cells=cells,.cell_count=2};size_t changed;
    ASSERT_EQ(alea_compute_thread_jumps(&s,&changed),0);ASSERT_EQ(changed,0);ASSERT_EQ(memcmp(copy,code,sizeof(code)),0);
    cells[1].count=8;ASSERT_EQ(alea_compute_thread_jumps(&s,&changed),-1);ASSERT_EQ(memcmp(copy,code,sizeof(code)),0);
    s.cell_count=1;code[1].arg=100;ASSERT_EQ(alea_compute_thread_jumps(&s,&changed),-1);code[1].arg=2;
    alea_interrupt();int rc=alea_compute_thread_jumps(&s,&changed);alea_clear_interrupt();ASSERT_EQ(rc,-1);
    ASSERT_EQ(memcmp(copy,code,sizeof(code)),0);
    code[4].op=COMPUTE_SKIP_UNION;ASSERT_EQ(alea_compute_thread_jumps(&s,&changed),0);ASSERT_EQ(changed,0);
    code[1].op=code[4].op=COMPUTE_SKIP_REVERSE_DIFFERENCE;
    ASSERT_EQ(alea_compute_thread_jumps(&s,&changed),0);ASSERT_EQ(changed,0);
}
#include "../../benchmarks/region_program_audit.h"
TEST(region_and_recognizer_checks_structure_and_compact_encoding) {
    compute_node nodes[]={{COMPUTE_LEAF,1},{COMPUTE_SKIP_INTERSECTION,5},
        {COMPUTE_LEAF,2},{COMPUTE_SKIP_INTERSECTION,2},{COMPUTE_LEAF,0},
        {COMPUTE_INTERSECTION,0},{COMPUTE_INTERSECTION,0}};
    compute_scene s={.nodes=nodes,.node_count=7,.primitive_count=2};
    compute_cell c={.root=0,.count=7};uint32_t leaves;
    ASSERT_EQ(region_and_program(&s,c,&leaves),1);ASSERT_EQ(leaves,3);
    uint32_t compact[7];ASSERT_EQ(alea_compute_code_compact(nodes,7,compact),1);
    s.compact_nodes=compact;ASSERT_EQ(region_and_program(&s,c,&leaves),1);ASSERT_EQ(leaves,3);s.compact_nodes=NULL;
    nodes[1].arg=4;ASSERT_EQ(region_and_program(&s,c,&leaves),0);ASSERT_EQ(leaves,0);
    nodes[1].arg=100;ASSERT_EQ(region_and_program(&s,c,&leaves),-1);nodes[1].arg=5;
    nodes[2].arg=4;ASSERT_EQ(region_and_program(&s,c,&leaves),-1);nodes[2].arg=2;
    nodes[5].op=COMPUTE_UNION;ASSERT_EQ(region_and_program(&s,c,&leaves),0);nodes[5].op=COMPUTE_INTERSECTION;
    c.count=8;ASSERT_EQ(region_and_program(&s,c,&leaves),-1);
    c.count=1;ASSERT_EQ(region_and_program(&s,c,&leaves),1);ASSERT_EQ(leaves,1);
    c.count=2;ASSERT_EQ(region_and_program(&s,c,&leaves),-1);
    c.count=0;ASSERT_EQ(region_and_program(&s,c,&leaves),-1);
}
#include "../../benchmarks/region_group_audit.h"
typedef struct {region_group groups[8];size_t count;int fail;} group_test_output;
static int test_group_emit(void* user,region_group g) {
    group_test_output* out=user;
    if(out->fail || out->count==8)return -1;
    out->groups[out->count++]=g;return 0;
}
TEST(region_groups_preserve_mixed_structure_and_ranges) {
    compute_node nodes[]={
        {COMPUTE_LEAF,0},{COMPUTE_SKIP_INTERSECTION,2},{COMPUTE_LEAF,1},{COMPUTE_INTERSECTION,0},
        {COMPUTE_SKIP_UNION,5},
        {COMPUTE_LEAF,2},{COMPUTE_SKIP_INTERSECTION,2},{COMPUTE_LEAF,3},{COMPUTE_INTERSECTION,0},
        {COMPUTE_UNION,0},{COMPUTE_COMPLEMENT,0}};
    compute_scene s={.nodes=nodes,.node_count=11,.primitive_count=2};
    compute_cell c={.root=0,.count=11};group_test_output out={0};
    ASSERT_EQ(region_group_program(&s,c,test_group_emit,&out),0);ASSERT_EQ(out.count,2);
    ASSERT_EQ(out.groups[0].begin,0);ASSERT_EQ(out.groups[0].end,4);
    ASSERT_EQ(out.groups[1].begin,5);ASSERT_EQ(out.groups[1].end,9);
    uint32_t compact[11];ASSERT_EQ(alea_compute_code_compact(nodes,11,compact),1);
    s.compact_nodes=compact;out.count=0;ASSERT_EQ(region_group_program(&s,c,test_group_emit,&out),0);ASSERT_EQ(out.count,2);s.compact_nodes=NULL;
    out.fail=1;ASSERT_EQ(region_group_program(&s,c,test_group_emit,&out),-1);out.fail=0;
    nodes[4].arg=4;ASSERT_EQ(region_group_program(&s,c,test_group_emit,&out),-1);nodes[4].arg=5;
    c.root=5;c.count=4;out.count=0;
    ASSERT_EQ(region_group_program(&s,c,test_group_emit,&out),0);ASSERT_EQ(out.count,1);
    ASSERT_EQ(out.groups[0].begin,5);ASSERT_EQ(out.groups[0].end,9);
    c.count=1;out.count=0;ASSERT_EQ(region_group_program(&s,c,test_group_emit,&out),0);ASSERT_EQ(out.count,0);
    nodes[5].op=COMPUTE_CONSTANT;nodes[5].arg=0;c.count=4;out.count=0;
    ASSERT_EQ(region_group_program(&s,c,test_group_emit,&out),0);ASSERT_EQ(out.count,0);
    nodes[5].op=COMPUTE_LEAF;nodes[5].arg=4;ASSERT_EQ(region_group_program(&s,c,test_group_emit,&out),-1);
    alea_interrupt();ASSERT_EQ(region_group_program(&s,c,test_group_emit,&out),-1);alea_clear_interrupt();
}
TEST(group_plan_aliases_budget_interruption_and_cleanup) {
    compute_node code[]={{COMPUTE_LEAF,0},{COMPUTE_SKIP_INTERSECTION,5},
        {COMPUTE_LEAF,2},{COMPUTE_SKIP_INTERSECTION,2},{COMPUTE_LEAF,4},
        {COMPUTE_INTERSECTION,0},{COMPUTE_INTERSECTION,0}};
    compute_node original[7];memcpy(original,code,sizeof(code));
    compute_cell cells[]={{.root=0,.count=7},{.root=2,.count=4},{.root=0,.count=1},{.root=0,.count=7}};
    compute_scene scene={.nodes=code,.node_count=7,.cells=cells,.cell_count=4,.primitive_count=3};
    compute_group_plan plan;
    ASSERT_EQ(alea_compute_groups_plan(&scene,1,&plan),0);ASSERT_NULL(plan.programs);ASSERT_NULL(plan.spans);
    ASSERT_EQ(alea_compute_groups_plan(&scene,78,&plan),0);ASSERT_NULL(plan.programs);ASSERT_NULL(plan.spans);
    ASSERT_EQ(alea_compute_groups_plan(&scene,65536,&plan),1);
    ASSERT_EQ(plan.program_count,3);ASSERT_EQ(plan.span_count,2);
    ASSERT_EQ(plan.programs[0].count,1);ASSERT_EQ(plan.programs[0].groups,0);
    ASSERT_EQ(plan.programs[1].count,7);ASSERT_EQ(plan.programs[1].groups,1);
    ASSERT_EQ(plan.programs[2].root,2);ASSERT_EQ(plan.programs[2].groups,1);
    ASSERT_EQ(plan.spans[0].begin,0);ASSERT_EQ(plan.spans[0].end,7);
    ASSERT_EQ(plan.spans[1].begin,2);ASSERT_EQ(plan.spans[1].end,6);
    alea_compute_groups_free(&plan);alea_compute_groups_free(&plan);
    const unsigned stops[]={1,4,6};
    for(unsigned i=0;i<3;++i) {
        checkpoints=0;interrupt_at=stops[i];int rc=alea_compute_groups_plan(&scene,65536,&plan);interrupt_at=0;
        ASSERT_EQ(rc,-1);ASSERT_EQ(alea_error_code(),ALEA_ERR_INTERRUPTED);
        ASSERT_NULL(plan.programs);ASSERT_NULL(plan.spans);
    }
    ASSERT_EQ(memcmp(original,code,sizeof(code)),0);
    code[1].arg=4;ASSERT_EQ(alea_compute_groups_plan(&scene,65536,&plan),-1);
    ASSERT_EQ(alea_error_code(),ALEA_ERR_INVALID_STATE);ASSERT_NULL(plan.programs);ASSERT_NULL(plan.spans);code[1].arg=5;
    scene.cell_count=1;cells[0].count=1;
    ASSERT_EQ(alea_compute_groups_plan(&scene,65536,&plan),0);ASSERT_NULL(plan.programs);
}
TEST_MAIN()
