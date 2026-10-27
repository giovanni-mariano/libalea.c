// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* Persistent per-pixel continuations. The host queries this Metal-only layout
 * rather than duplicating SIMD alignment rules in C. A bounded pool is resident.
 * Each work unit is one CSG opcode or one bounded traversal/state transition. */
enum { RESUME_BEGIN, RESUME_WALK, RESUME_INTERVAL, RESUME_VALIDATE, RESUME_SHADE, RESUME_DONE };
enum { RESUME_PARENT=1, RESUME_SOLID=2, RESUME_GUARD=3 };
struct resume_expression {
    int values[COMPUTE_STACK];float until[COMPUTE_STACK];uint boundaries[COMPUTE_STACK];
    uint pc,sp,at_cursor;compute_cell program;float end;
};
struct resume_interval {
    uint id,mode,phase,steps,entering,surface,boundary;bool found,hit;
    float cursor,end,enter,exit;
};
struct resume_ray {
    render_frame stack[COMPUTE_HIER_STACK];uint deferred[COMPUTE_BVH_STACK];
    int top;uint deferred_count,steps,chosen,stage,sample;
    float begin,best;float3 world_d,clip,chosen_normal,sum;
    compute_render_pixel out,center;int error;
    compute_profile cost,interval_cost;
    resume_expression expression;resume_interval interval;
    uint guard_frame,candidate_id,candidate_surface;
    float candidate_t,candidate_end;float3 candidate_normal;
};
kernel void alea_resume_layout(device uint* size [[buffer(0)]],uint i [[thread_position_in_grid]]) {
    if(i==0)*size=sizeof(resume_ray);
}
static int resume_expression_step(device const compute_node* code,device const compute_primitive* primitives,
        thread resume_expression& e,float3 o,float3 d,float cursor,
        thread bool& error,thread compute_profile& cost) {
    float end=e.end;
    compute_node n=code[e.program.root+e.pc];
    if(profile_enabled)cost.instructions=profile_sum(cost.instructions,1);
    if(n.op==COMPUTE_CONSTANT) {
        e.values[e.sp]=n.arg?1:-1;e.until[e.sp]=end;e.boundaries[e.sp]=COMPUTE_NONE;++e.sp;
    } else if(n.op==COMPUTE_LEAF) {
        uint id=n.arg>>1;float4 roots;float crossing=end;uint boundary=COMPUTE_NONE;
        if(profile_enabled)cost.primitives=profile_sum(cost.primitives,1);
        uint count=primitive_roots(primitives[id],o,d,roots,error);if(error)return -1;
        for(uint j=0;j<count;++j) {
            if(roots[j]==cursor)e.at_cursor=id;
            if(roots[j]>cursor && roots[j]<crossing){crossing=roots[j];boundary=id;}
        }
        float probe=cursor+(crossing-cursor)*.5f;
        /* Do not evaluate astronomical points just because this primitive
         * has no further crossing. A finite interior probe is sufficient. */
        probe=min(probe,cursor+max(1.0f,abs(cursor)*.001f));
        if(probe==cursor)probe=crossing;
        float v=primitive_value(primitives[id],o+probe*d);
        if(!isfinite(v)){error=true;return -1;}
        int sign=(v>0)-(v<0);
        e.values[e.sp]=(n.arg&1)?-sign:sign;e.until[e.sp]=crossing;e.boundaries[e.sp]=boundary;++e.sp;
    } else if(n.op==COMPUTE_SKIP_UNION || n.op==COMPUTE_SKIP_INTERSECTION || n.op==COMPUTE_SKIP_REVERSE_DIFFERENCE) {
        int v=e.values[e.sp-1];
        if((n.op==COMPUTE_SKIP_UNION && v<=0) || (n.op==COMPUTE_SKIP_INTERSECTION && v>0) ||
           (n.op==COMPUTE_SKIP_REVERSE_DIFFERENCE && v<0)) {
            if(n.op==COMPUTE_SKIP_REVERSE_DIFFERENCE)e.values[e.sp-1]=-v;
            e.pc+=n.arg;
        }
    } else if(n.op==COMPUTE_COMPLEMENT)e.values[e.sp-1]=-e.values[e.sp-1];
    else {
        uint b=--e.sp,a=e.sp-1;
        e.values[a]=n.op==COMPUTE_UNION?min(e.values[a],e.values[b]):
            n.op==COMPUTE_INTERSECTION?max(e.values[a],e.values[b]):
            n.op==COMPUTE_DIFFERENCE?max(e.values[a],-e.values[b]):max(e.values[b],-e.values[a]);
        if(e.until[b]<e.until[a]){e.until[a]=e.until[b];e.boundaries[a]=e.boundaries[b];}
    }
    ++e.pc;return 0;
}
static void resume_start_interval(thread resume_ray& s,uint id,uint mode,float begin,float end) {
    s.interval={};s.interval.id=id;s.interval.mode=mode;s.interval.cursor=begin;s.interval.end=end;
    s.interval.entering=COMPUTE_NONE;s.interval.surface=COMPUTE_NONE;s.interval_cost={};s.stage=RESUME_INTERVAL;
    s.expression.end=begin; /* Invalidate the previous cell's spatial selection. */
}
/* Return 0 suspended, 1 finished, -1 failed. mode 1 asks for the whole parent
 * interval; mode 2 stops at the first occupied point of a material cell. */
static int resume_interval_step(thread resume_ray& s,device const compute_node* code,
        device const compute_primitive* primitives,device const compute_cell* cells,
        device const compute_spatial_program* spatial) {
    thread resume_interval& v=s.interval;thread resume_expression& e=s.expression;
    thread render_frame& f=s.stack[v.mode==RESUME_GUARD?s.guard_frame-1:s.top];
    if(v.phase==0) {
        if(v.cursor>=v.end){v.hit=v.found;if(v.found)v.exit=v.end;return 1;}
        if(++v.steps>100000)return -1;
        if(profile_enabled)s.interval_cost.intervals=profile_sum(s.interval_cost.intervals,1);
        if(v.cursor>=e.end) {
            e.end=v.end;e.program=spatial_segment(cells[v.id],f.o,f.d,v.cursor,e.end,spatial);
        }
        e.pc=e.sp=0;e.at_cursor=COMPUTE_NONE;v.phase=1;return 0;
    }
    if(e.pc<e.program.count) {
        bool error=false;
        return resume_expression_step(code,primitives,e,f.o,f.d,v.cursor,error,s.interval_cost);
    }
    bool inside=e.values[0]<=0;float next=e.until[0];uint boundary=e.boundaries[0];
    if(v.mode==RESUME_GUARD) {
        if(next<=v.cursor)return -1;
        v.hit=inside;v.enter=v.cursor;v.exit=next;v.surface=e.at_cursor;v.boundary=boundary;return 1;
    }
    if(v.entering==COMPUTE_NONE)v.entering=e.at_cursor;
    if(inside && !v.found) {
        v.found=true;v.enter=v.cursor;v.surface=v.entering;
        if(v.mode==RESUME_SOLID){v.exit=next;v.hit=true;return 1;}
    } else if(!inside && v.found){v.exit=v.cursor;v.hit=true;return 1;}
    if(next<=v.cursor)return -1;
    v.cursor=next;v.entering=boundary;v.phase=0;return 0;
}
static void resume_profile_interval(thread resume_ray& s,device compute_profile* cell_costs) {
    if(!profile_enabled)return;
    s.cost.cells=profile_sum(s.cost.cells,1);
    s.cost.instructions=profile_sum(s.cost.instructions,s.interval_cost.instructions);
    s.cost.primitives=profile_sum(s.cost.primitives,s.interval_cost.primitives);
    s.cost.intervals=profile_sum(s.cost.intervals,s.interval_cost.intervals);
    device atomic_uint* dst=(device atomic_uint*)(cell_costs+s.interval.id);
    profile_atomic_add(dst+1,1);profile_atomic_add(dst+2,s.interval_cost.instructions);
    profile_atomic_add(dst+3,s.interval_cost.primitives);profile_atomic_add(dst+4,s.interval_cost.intervals);
}
static void resume_descend(thread resume_ray& s,compute_cell c,uint child,float lo,float hi,
        float3 normal,float3 translation,device const compute_transform* transforms,uint guard,bool boundary) {
    thread render_frame& f=s.stack[s.top];f.cursor=hi;
    if(child==COMPUTE_NONE)return;
    if(s.top+1>=COMPUTE_HIER_STACK){s.out.error=4;s.stage=RESUME_SHADE;return;}
    render_frame next={f.o,f.d,normal,f.world_to_local,lo,hi,lo,hi,child,s.deferred_count,COMPUTE_NONE};
    next.guard=guard;
    next.entry_owner=boundary?uint(s.top+1):f.entry_owner;next.guard_surface=next.guard_boundary=COMPUTE_NONE;
    if(c.lattice!=COMPUTE_NONE)next.o-=translation;
    else {
        device const float* m=transforms[c.transform].inverse;
        float3x3 linear(float3(m[0],m[4],m[8]),float3(m[1],m[5],m[9]),float3(m[2],m[6],m[10]));
        next.o=linear*f.o+float3(m[3],m[7],m[11]);next.d=linear*f.d;
        next.world_to_local=linear*f.world_to_local;
    }
    if(!all(isfinite(next.o)) || !all(isfinite(next.d)) ||
       !all(isfinite(next.world_to_local[0])) || !all(isfinite(next.world_to_local[1])) ||
       !all(isfinite(next.world_to_local[2]))) {s.out.error=1;s.stage=RESUME_SHADE;return;}
    s.stack[++s.top]=next;
}
static void resume_begin(thread resume_ray& s,constant compute_render_params& p,float2 pixel) {
    float3 world_o(p.eye[0],p.eye[1],p.eye[2]),world_d(p.forward[0],p.forward[1],p.forward[2]);
    float3 right(p.right[0],p.right[1],p.right[2]),up(p.up[0],p.up[1],p.up[2]);
    float u=(2*pixel.x/p.width-1)*p.half_width,v=(1-2*pixel.y/p.height)*p.half_height;
    if(p.perspective)world_d=normalize(world_d+u*right+v*up);else world_o+=u*right+v*up;
    s.out={};thread compute_render_pixel& out=s.out;out.cell_id=-1;out.depth=1e30f;
    for(uint a=0;a<3;++a)out.color[a]=p.background[a];
    if(!all(isfinite(world_o)) || !all(isfinite(world_d))){out.error=1;s.stage=RESUME_SHADE;return;}
    if(p.bvh_root==COMPUTE_NONE){s.stage=RESUME_SHADE;return;}
    float begin=0,end=1e15f;float3 clip(p.clip[0],p.clip[1],p.clip[2]);
    if(p.clip_enabled) {
        float den=dot(clip,world_d),num=dot(clip,world_o)+p.clip[3];
        if(!isfinite(den) || !isfinite(num)){out.error=1;s.stage=RESUME_SHADE;return;}
        if(den==0){if(num<=0){s.stage=RESUME_SHADE;return;}}
        else if(den>0)begin=max(begin,-num/den);else end=min(end,-num/den);
        if(end<=begin){s.stage=RESUME_SHADE;return;}
    }
    s.world_d=world_d;s.clip=clip;s.begin=begin;s.best=end;s.chosen=COMPUTE_NONE;s.chosen_normal=float3(0);
    s.top=0;s.deferred_count=0;s.steps=0;
    s.stack[0]={world_o,world_d,float3(0),float3x3(1),begin,end,begin,end,p.bvh_root,0,COMPUTE_NONE};
    s.stack[0].guard=COMPUTE_NONE;
    s.stage=RESUME_WALK;
}
static void resume_walk(thread resume_ray& s,device const compute_primitive* primitives,
        device const compute_cell* cells,device const compute_render_bvh* bvh,
        device const compute_transform* transforms,device const compute_lattice* lattices,device const uint* fills) {
    if(s.top<0){s.stage=RESUME_SHADE;return;}
    if(profile_enabled)s.cost.hierarchy=profile_sum(s.cost.hierarchy,1);
    if(++s.steps>1000000){s.out.error=2;s.stage=RESUME_SHADE;return;}
    thread render_frame& f=s.stack[s.top];float limit=min(s.best,f.end);
    if(f.begin>=limit){s.deferred_count=f.deferred_start;--s.top;return;}
    if(f.pending!=COMPUTE_NONE) {
        limit=min(limit,f.pending_end);uint id=f.pending;compute_cell c=cells[id];
        if(f.cursor>=limit){f.pending=COMPUTE_NONE;return;}
        if(c.lattice==COMPUTE_NONE) {
            if(defer_parent(c))resume_descend(s,c,c.child,f.cursor,limit,f.entry_normal,float3(0),transforms,id,false);
            else {
                if(profile_enabled)s.cost.parents=profile_sum(s.cost.parents,1);
                resume_start_interval(s,id,RESUME_PARENT,f.cursor,limit);
            }
            return;
        }
        float lo=f.cursor,hi=limit;uint child=c.child;float3 translation,local_normal;bool error=false;
        if(!lattice_segment(lattices[c.lattice],fills,f.o,f.d,lo,limit,hi,child,translation,local_normal,error)) {
            if(error){s.out.error=3;s.stage=RESUME_SHADE;}else f.pending=COMPUTE_NONE;return;
        }
        float3 normal=lo>f.begin?to_world_normal(f,local_normal):f.entry_normal;
        resume_descend(s,c,child,lo,hi,normal,translation,transforms,COMPUTE_NONE,lo>f.begin);return;
    }
    if(f.node==COMPUTE_NONE) {
        if(s.deferred_count==f.deferred_start){--s.top;return;}
        f.node=s.deferred[--s.deferred_count];
    }
    if(profile_enabled)s.cost.bvh=profile_sum(s.cost.bvh,1);
    uint node=f.node;device const compute_render_bvh& b=bvh[node];
    f.node=COMPUTE_NONE;
    float boxlo=f.begin,boxhi=limit;
    if(!render_box(b,f.o,f.d,boxlo,boxhi))return;
    if(b.cell==COMPUTE_NONE) {
        uint left=node+1,right=bvh[left].escape;
        float ll=f.begin,lh=limit,rl=f.begin,rh=limit;
        bool hit_left=render_box(bvh[left],f.o,f.d,ll,lh),hit_right=render_box(bvh[right],f.o,f.d,rl,rh);
        if(hit_left && hit_right) {
            if(s.deferred_count==COMPUTE_BVH_STACK){s.out.error=5;s.stage=RESUME_SHADE;return;}
            f.node=ll<=rl?left:right;s.deferred[s.deferred_count++]=ll<=rl?right:left;
        } else if(hit_left)f.node=left;else if(hit_right)f.node=right;
        return;
    }
    uint id=b.cell;compute_cell c=cells[id];
    if(c.lattice==COMPUTE_NONE && c.transform!=COMPUTE_NONE && c.child==COMPUTE_NONE)return;
    if(c.lattice!=COMPUTE_NONE && lattices[c.lattice].repeating && fills[lattices[c.lattice].offset]==COMPUTE_NONE)return;
    if(c.lattice!=COMPUTE_NONE) {
        device const compute_lattice& l=lattices[c.lattice];
        if(!l.repeating && l.outer==COMPUTE_NONE) {
            /* BVH boxes are outward rounded; clip against the actual
             * finite grid before stepping its last element. */
            bool miss=false;
            for(uint a=0;a<3;++a)if(l.pitch[a]>0) {
                float low=l.lower[a],high=low+float(l.dims[2*a+1]-l.dims[2*a]+1)*l.pitch[a];
                if(f.d[a]==0){if(f.o[a]<low || f.o[a]>high)miss=true;}
                else {
                    float x=(low-f.o[a])/f.d[a],y=(high-f.o[a])/f.d[a];
                    boxlo=max(boxlo,min(x,y));boxhi=min(boxhi,max(x,y));
                }
            }
            if(miss || boxlo>=boxhi)return;
        }
    }
    if(c.transform!=COMPUTE_NONE || c.lattice!=COMPUTE_NONE) {
        f.pending=id;
        float padding=max(.001f,max(abs(boxlo),abs(boxhi))*1e-5f);
        f.pending_end=c.lattice!=COMPUTE_NONE?boxhi:min(limit,boxhi+padding);
        /* Only a lattice grid box is itself a real boundary. Ordinary
         * cell boxes are conservative and must not replace CSG crossings. */
        f.cursor=c.lattice!=COMPUTE_NONE?boxlo:max(f.begin,boxlo-padding);
        return;
    }
    float padding=max(.001f,max(abs(boxlo),abs(boxhi))*1e-5f);
    s.candidate_end=min(limit,boxhi+padding);
    resume_start_interval(s,id,RESUME_SOLID,max(f.begin,boxlo-padding),s.candidate_end);
}
/* Continue after one cached or freshly evaluated ancestor truth segment. */
static void resume_guard_result(thread resume_ray& s,device const compute_primitive* primitives) {
    thread render_frame& guard=s.stack[s.guard_frame];
    if(!guard.guard_inside) {
        if(profile_enabled)s.cost.rejected=profile_sum(s.cost.rejected,1);
        if(guard.guard_hi<=s.candidate_t){s.out.error=1;s.stage=RESUME_SHADE;return;}
        resume_start_interval(s,s.candidate_id,RESUME_SOLID,guard.guard_hi,s.candidate_end);return;
    }
    if(s.candidate_surface==COMPUTE_NONE && s.candidate_t==guard.guard_lo && guard.guard_surface!=COMPUTE_NONE &&
       (s.candidate_t!=s.stack[s.top].begin || s.guard_frame>=s.stack[s.top].entry_owner)) {
        thread render_frame& owner=s.stack[s.guard_frame-1];
        s.candidate_normal=to_world_normal(owner,surface_normal(primitives[guard.guard_surface],owner.o+s.candidate_t*owner.d));
    }
    ++s.guard_frame;s.stage=RESUME_VALIDATE;
}
static void resume_validate(thread resume_ray& s,device const compute_primitive* primitives) {
    if(s.guard_frame>uint(s.top)) {
        if(s.candidate_t<s.best){s.best=s.candidate_t;s.chosen=s.candidate_id;s.chosen_normal=s.candidate_normal;}
        s.stage=RESUME_WALK;return;
    }
    thread render_frame& guard=s.stack[s.guard_frame];
    if(guard.guard==COMPUTE_NONE){++s.guard_frame;return;}
    if(profile_enabled)s.cost.parents=profile_sum(s.cost.parents,1);
    if(s.candidate_t>=guard.guard_lo && s.candidate_t<guard.guard_hi) {
        if(profile_enabled)s.cost.cache_hits=profile_sum(s.cost.cache_hits,1);
        resume_guard_result(s,primitives);return;
    }
    resume_start_interval(s,guard.guard,RESUME_GUARD,s.candidate_t,s.candidate_end);
}
static void resume_finish_interval(thread resume_ray& s,device const compute_primitive* primitives,
        device const compute_cell* cells,device const compute_transform* transforms) {
    thread resume_interval& v=s.interval;thread render_frame& f=s.stack[s.top];s.stage=RESUME_WALK;
    if(v.mode==RESUME_GUARD) {
        thread render_frame& guard=s.stack[s.guard_frame];
        uint entering=guard.guard_hi==v.enter && guard.guard_boundary!=COMPUTE_NONE?guard.guard_boundary:v.surface;
        guard.guard_lo=v.enter;guard.guard_hi=v.exit;guard.guard_inside=v.hit;
        guard.guard_surface=entering;guard.guard_boundary=v.boundary;
        resume_guard_result(s,primitives);
    } else if(v.mode==RESUME_PARENT) {
        if(!v.hit){f.pending=COMPUTE_NONE;return;}
        float3 normal=v.surface==COMPUTE_NONE?f.entry_normal:
            to_world_normal(f,surface_normal(primitives[v.surface],f.o+v.enter*f.d));
        compute_cell c=cells[v.id];resume_descend(s,c,c.child,v.enter,v.exit,normal,float3(0),transforms,COMPUTE_NONE,v.surface!=COMPUTE_NONE);
    } else if(v.hit && v.enter<s.best) {
        float3 normal=v.surface==COMPUTE_NONE?f.entry_normal:
            to_world_normal(f,surface_normal(primitives[v.surface],f.o+v.enter*f.d));
        if(lazy_threshold) {
            s.candidate_t=v.enter;s.candidate_id=v.id;s.candidate_surface=v.surface;
            s.candidate_normal=normal;s.guard_frame=1;s.stage=RESUME_VALIDATE;
        } else {s.best=v.enter;s.chosen=v.id;s.chosen_normal=normal;}
    }
}
static void resume_shade(thread resume_ray& s,device const compute_cell* cells,
        device const float* colors,constant compute_render_params& p) {
    if(s.out.error || s.chosen==COMPUTE_NONE)return;
    thread compute_render_pixel& out=s.out;
    float3 normal=s.chosen_normal,base(colors[3*s.chosen],colors[3*s.chosen+1],colors[3*s.chosen+2]);
    bool cut=p.clip_enabled && s.begin>0 && s.best==s.begin;
    if(cut)normal=-normalize(s.clip);if(dot(normal,s.world_d)>0)normal=-normal;
    float3 color;
    if(cut)color=base*p.tint;
    else {
        float3 light=p.light_set?float3(p.light[0],p.light[1],p.light[2]):s.world_d;
        float ndl=max(-dot(normal,light),0.0f);
        float shine=pow(max(-dot(2*ndl*normal+light,s.world_d),0.0f),p.shininess);
        color=base*(p.ambient+p.diffuse*ndl)+p.specular*shine;
    }
    if(!all(isfinite(color)) || !all(isfinite(normal))){out.error=1;return;}
    color=clamp(color,0.0f,1.0f);
    for(uint a=0;a<3;++a){out.color[a]=color[a];out.normal[a]=normal[a];}
    out.depth=s.best;out.cell_id=cells[s.chosen].cell_id;out.material_id=cells[s.chosen].material_id;return;
}
kernel void alea_render_resume(device const compute_node* code [[buffer(0)]],
        device const compute_cell* cells [[buffer(1)]],device const compute_render_bvh* bvh [[buffer(2)]],
        device const float* colors [[buffer(3)]],device compute_render_pixel* pixels [[buffer(4)]],
        constant compute_render_params& p [[buffer(5)]],device const compute_primitive* primitives [[buffer(6)]],
        device const compute_transform* transforms [[buffer(7)]],device const compute_lattice* lattices [[buffer(8)]],
        device const uint* fills [[buffer(9)]],constant compute_render_tile& tile [[buffer(10)]],
        device compute_profile* ray_costs [[buffer(11)]],device compute_profile* cell_costs [[buffer(12)]],
        device const compute_spatial_program* spatial [[buffer(15)]],
        device resume_ray* states [[buffer(16)]],device uint* done [[buffer(17)]],
        constant compute_resume_params& dispatch [[buffer(18)]],
        device const uint* slot_pixels [[buffer(20)]],
        uint index [[thread_position_in_grid]]) {
    bool reset=dispatch.reset;
    if(refill_enabled) {
        if(done[index]==2)return;
        reset=done[index]==3;
    } else if(!reset && done[index])return;
    resume_ray s;
    if(reset){s={};done[index]=0;}else s=states[index];
    uint pixel=refill_enabled?slot_pixels[index]:(tile.y+index/tile.width)*p.width+tile.x+index%tile.width;
    float2 origin(pixel%p.width,pixel/p.width);
    for(uint work=0;work<dispatch.work && s.stage!=RESUME_DONE;++work) {
        if(s.stage==RESUME_BEGIN) {
            // Clear chosen before early camera/clip rejection of an AA sample.
            s.chosen=COMPUTE_NONE;
            float2 sample(float(s.sample%p.aa),float(s.sample/p.aa));
            resume_begin(s,p,origin+(sample+.5f)/float(p.aa));
        } else if(s.stage==RESUME_WALK)resume_walk(s,primitives,cells,bvh,transforms,lattices,fills);
        else if(s.stage==RESUME_VALIDATE)resume_validate(s,primitives);
        else if(s.stage==RESUME_INTERVAL) {
            int result=resume_interval_step(s,code,primitives,cells,spatial);
            if(result) {
                resume_profile_interval(s,cell_costs);
                if(result<0){s.out.error=1;s.stage=RESUME_SHADE;}
                else resume_finish_interval(s,primitives,cells,transforms);
            }
        } else {
            resume_shade(s,cells,colors,p);
            s.sum+=float3(s.out.color[0],s.out.color[1],s.out.color[2]);s.error=max(s.error,s.out.error);
            if(s.sample%p.aa==p.aa/2 && s.sample/p.aa==p.aa/2)s.center=s.out;
            ++s.sample;s.stage=s.sample==p.aa*p.aa?RESUME_DONE:RESUME_BEGIN;
        }
    }
    if(s.stage==RESUME_DONE) {
        for(uint a=0;a<3;++a)s.center.color[a]=s.sum[a]/float(p.aa*p.aa);
        s.center.error=s.error;pixels[pixel]=s.center;
        if(profile_enabled)ray_costs[pixel]=s.cost;
        done[index]=1;
    } else states[index]=s;
}
