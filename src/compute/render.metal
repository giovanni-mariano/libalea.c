// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "intersect.metal"

constant bool profile_enabled [[function_constant(0)]];
constant uint lazy_threshold [[function_constant(2)]];
constant bool refill_enabled [[function_constant(4)]];
static uint profile_sum(uint a,uint b) { return a+min(b,0xffffffffu-a); }
static void profile_atomic_add(device atomic_uint* p,uint value) {
    uint old=atomic_load_explicit(p,memory_order_relaxed);
    while(!atomic_compare_exchange_weak_explicit(p,&old,profile_sum(old,value),
          memory_order_relaxed,memory_order_relaxed)) {}
}


/* Evaluate a CSG program over the open ray interval beginning at cursor.
 * Each value carries the next crossing up to which its truth is constant.
 * Short-circuiting also skips surface intersections: a failed intersection
 * remains false until that operand changes; a satisfied union remains true.
 * This avoids sorting every surface in large cell expressions. */
static bool ray_expression(device const compute_node* code,device const compute_primitive* primitives,
        compute_cell cell,float3 o,float3 d,float cursor,float end,
        thread float& next,thread uint& surface,thread uint& at_cursor,thread bool& error,thread compute_profile& cost) {
    int values[COMPUTE_STACK];float until[COMPUTE_STACK];uint boundaries[COMPUTE_STACK];uint sp=0;
    for(uint i=0;i<cell.count;++i) {
        compute_node n=code_at(code,cell.root+i);
        if(profile_enabled)cost.instructions=profile_sum(cost.instructions,1);
        if(n.op==COMPUTE_CONSTANT) {
            values[sp]=n.arg?1:-1;until[sp]=end;boundaries[sp]=COMPUTE_NONE;++sp;
        } else if(n.op==COMPUTE_LEAF) {
            uint id=n.arg>>1;float4 roots;float crossing=end;uint boundary=COMPUTE_NONE;
            if(profile_enabled)cost.primitives=profile_sum(cost.primitives,1);
            uint count=primitive_roots(primitives[id],o,d,roots,error);if(error)return false;
            for(uint j=0;j<count;++j) {
                if(roots[j]==cursor)at_cursor=id;
                if(roots[j]>cursor && roots[j]<crossing){crossing=roots[j];boundary=id;}
            }
            float probe=cursor+(crossing-cursor)*.5f;
            /* Do not evaluate astronomical points just because this primitive
             * has no further crossing. A finite interior probe is sufficient. */
            probe=min(probe,cursor+max(1.0f,abs(cursor)*.001f));
            if(probe==cursor)probe=crossing;
            float v=primitive_value(primitives[id],o+probe*d);
            if(!isfinite(v)){error=true;return false;}
            int sign=(v>0)-(v<0);values[sp]=(n.arg&1)?-sign:sign;until[sp]=crossing;boundaries[sp]=boundary;++sp;
        } else if(n.op==COMPUTE_SKIP_UNION || n.op==COMPUTE_SKIP_INTERSECTION || n.op==COMPUTE_SKIP_REVERSE_DIFFERENCE) {
            int v=values[sp-1];
            if((n.op==COMPUTE_SKIP_UNION && v<=0) || (n.op==COMPUTE_SKIP_INTERSECTION && v>0) ||
               (n.op==COMPUTE_SKIP_REVERSE_DIFFERENCE && v<0)) {
                if(n.op==COMPUTE_SKIP_REVERSE_DIFFERENCE)values[sp-1]=-v;
                i+=n.arg;
            }
        } else if(n.op==COMPUTE_COMPLEMENT)values[sp-1]=-values[sp-1];
        else {
            uint b=--sp,a=sp-1;
            values[a]=n.op==COMPUTE_UNION?min(values[a],values[b]):
                n.op==COMPUTE_INTERSECTION?max(values[a],values[b]):
                n.op==COMPUTE_DIFFERENCE?max(values[a],-values[b]):max(values[b],-values[a]);
            if(until[b]<until[a]){until[a]=until[b];boundaries[a]=boundaries[b];}
        }
    }
    next=until[0];surface=boundaries[0];return values[0]<=0;
}

/* Spatial proofs overlap their logical boxes. Evaluate only ray segments
 * covered by the proof, and retain the original program for poorly conditioned
 * coordinates. Spatial split planes never provide a hit or shading normal. */
static compute_cell spatial_segment(compute_cell cell,float3 o,float3 d,float cursor,
        thread float& end,device const compute_spatial_program* spatial) {
    if(cell.spatial==COMPUTE_NONE)return cell;
    uint id=cell.spatial;float lo=cursor,hi=end;
    if(!ray_box(spatial[id].bounds,o,d,lo,hi) || hi<=cursor)return cell;
    if(lo>cursor){end=min(end,lo);return cell;}
    float next=hi;
    while(spatial[id].child!=COMPUTE_NONE) {
        uint child=spatial[id].child,a=spatial[id].axis;
        float split=spatial[child+1].bounds[2*a];
        bool right=o[a]>=split;
        if(d[a]!=0) {
            float crossing=(split-o[a])/d[a];
            right=d[a]>0?crossing<=cursor:crossing>cursor;
            if(crossing>cursor)next=min(next,crossing);
        }
        id=child+uint(right);
    }
    if(next<=cursor || !isfinite(next))return cell;
    float3 lower,upper,margin;
    for(uint a=0;a<3;++a) {
        lower[a]=spatial[id].bounds[2*a];upper[a]=spatial[id].bounds[2*a+1];
        margin[a]=max((upper[a]-lower[a])*.001f,max(abs(lower[a]),abs(upper[a]))*1e-5f);
    }
    float3 start=o+cursor*d,finish=o+next*d;
    float3 error=(abs(o)+max(abs(cursor),abs(next))*abs(d))*1e-6f;
    if(any(error>margin*.5f) || any(min(start,finish)-error<lower-margin*.5f) ||
        any(max(start,finish)+error>upper+margin*.5f))return cell;
    cell.root=spatial[id].root;cell.count=spatial[id].count;end=next;return cell;
}

static bool cell_interval_impl(device const compute_node* code,device const compute_primitive* primitives,
        compute_cell cell,
        float3 o,float3 d,float begin,float end,bool want_exit,
        thread float& enter,thread float& exit,thread uint& surface,thread bool& error,
        thread compute_profile& cost,device const compute_spatial_program* spatial) {
    float cursor=begin,segment_end=begin;uint entering=COMPUTE_NONE;bool found=false;uint steps=0;
    compute_cell program=cell;
    while(cursor<end) {
        if(++steps>100000){error=true;return false;}
        if(profile_enabled)cost.intervals=profile_sum(cost.intervals,1);
        float next;uint boundary,at_cursor=COMPUTE_NONE;
        /* The proof covers the whole segment. Reuse it across primitive
         * crossings, and choose again only after leaving the spatial region. */
        if(cursor>=segment_end) {
            segment_end=end;program=spatial_segment(cell,o,d,cursor,segment_end,spatial);
        }
        bool inside=ray_expression(code,primitives,program,o,d,cursor,segment_end,next,boundary,at_cursor,error,cost);
        if(error)return false;
        if(entering==COMPUTE_NONE)entering=at_cursor;
        if(inside && !found) {
            found=true;enter=cursor;surface=entering;
            if(!want_exit){exit=next;return true;}
        } else if(!inside && found){exit=cursor;return true;}
        if(next<=cursor){error=true;return false;}
        cursor=next;entering=boundary;
    }
    if(found)exit=end;
    return found;
}

static bool cell_interval(device const compute_node* code,device const compute_primitive* primitives,
        compute_cell cell,float3 o,float3 d,float begin,float end,bool want_exit,
        thread float& enter,thread float& exit,thread uint& surface,thread bool& error,
        thread compute_profile& ray,device compute_profile* cells,uint id,
        device const compute_spatial_program* spatial) {
    compute_profile cost={};
    bool hit=cell_interval_impl(code,primitives,cell,o,d,begin,end,want_exit,enter,exit,surface,error,cost,spatial);
    if(profile_enabled) {
        ray.cells=profile_sum(ray.cells,1);
        ray.instructions=profile_sum(ray.instructions,cost.instructions);
        ray.primitives=profile_sum(ray.primitives,cost.primitives);
        ray.intervals=profile_sum(ray.intervals,cost.intervals);
        device atomic_uint* dst=(device atomic_uint*)(cells+id);
        profile_atomic_add(dst+1,1);
        profile_atomic_add(dst+2,cost.instructions);
        profile_atomic_add(dst+3,cost.primitives);
        profile_atomic_add(dst+4,cost.intervals);
    }
    return hit;
}

/* A distant ray can lose small coordinates when transformed or reconstructed.
 * Keep the original hierarchy boxes in that regime. Child coordinate budgets
 * are propagated upward by the host compiler. */
static bool render_box(device const compute_render_bvh& b,float3 o,float3 d,thread float& lo,thread float& hi) {
    if(isinf(b.guard))return ray_box(b.bounds,o,d,lo,hi);
    float scale=max(abs(o.x),max(abs(o.y),abs(o.z))),extent=0;
    for(uint a=0;a<3;++a)extent=max(extent,max(abs(b.bounds[2*a]),abs(b.bounds[2*a+1])));
    bool safe=4*(scale+extent+1)<b.guard;
    return ray_box(safe?b.bounds:b.fallback,o,d,lo,hi);
}

struct render_frame {
    float3 o,d,entry_normal;
    float3x3 world_to_local;
    float begin,end,cursor,pending_end;
    uint node,deferred_start,pending;
    uint guard,guard_inside,guard_surface,guard_boundary,entry_owner;
    float guard_lo,guard_hi;
};
static bool defer_parent(compute_cell c) {
    return lazy_threshold && c.count>=lazy_threshold && !(c.flags&COMPUTE_CELL_LATTICE_SUBTREE) &&
        (lazy_threshold==1 || c.count/4>=(c.flags>>1));
}
/* A cache covers one open truth segment of this placed parent's local ray.
 * Check both bounds: near-first BVH traversal may revisit earlier distances. */
static bool parent_truth(device const compute_node* code,device const compute_primitive* primitives,
        device const compute_cell* cells,thread render_frame& cache,thread const render_frame& owner,
        float cursor,float limit,thread bool& error,thread compute_profile& cost,
        device compute_profile* cell_costs,device const compute_spatial_program* spatial) {
    if(profile_enabled)cost.parents=profile_sum(cost.parents,1);
    if(cursor>=cache.guard_lo && cursor<cache.guard_hi) {
        if(profile_enabled)cost.cache_hits=profile_sum(cost.cache_hits,1);
        return cache.guard_inside;
    }
    float end=limit;compute_cell program=spatial_segment(cells[cache.guard],owner.o,owner.d,cursor,end,spatial);
    uint surface,at_cursor=COMPUTE_NONE;float next;compute_profile work={};
    bool inside=ray_expression(code,primitives,program,owner.o,owner.d,cursor,end,next,surface,at_cursor,error,work);
    if(error || next<=cursor){error=true;return false;}
    uint entering=cache.guard_hi==cursor && cache.guard_boundary!=COMPUTE_NONE?cache.guard_boundary:at_cursor;
    cache.guard_lo=cursor;cache.guard_hi=next;cache.guard_inside=inside;
    cache.guard_surface=entering;cache.guard_boundary=surface;
    if(profile_enabled) {
        cost.cells=profile_sum(cost.cells,1);cost.intervals=profile_sum(cost.intervals,1);
        cost.instructions=profile_sum(cost.instructions,work.instructions);
        cost.primitives=profile_sum(cost.primitives,work.primitives);
        device atomic_uint* dst=(device atomic_uint*)(cell_costs+cache.guard);
        profile_atomic_add(dst+1,1);profile_atomic_add(dst+2,work.instructions);
        profile_atomic_add(dst+3,work.primitives);profile_atomic_add(dst+4,1);
    }
    return inside;
}
static float3 to_world_normal(thread const render_frame& f,float3 normal) {
    float3 v=transpose(f.world_to_local)*normal;
    return length_squared(v)>0?normalize(v):float3(0);
}

/* One lattice segment at a time. Repeating grids stay procedural: a ray only
 * visits the elements within its enclosing CSG interval, never an expanded
 * array of millions of placements. */
static float adjacent_float(float x,float direction) {
    if(direction==0)return x;
    if(x==0)return copysign(as_type<float>(1u),direction);
    uint bits=as_type<uint>(x);
    return as_type<float>((x>0)==(direction>0)?bits+1:bits-1);
}
static bool lattice_segment(device const compute_lattice& l,device const uint* fills,
        float3 o,float3 d,float begin,float end,thread float& next,thread uint& child,
        thread float3& translation,thread float3& normal,thread bool& error) {
    int3 logical(0),relative(0);float3 position=o+begin*d;
    float3 low,high;
    for(uint a=0;a<3;++a) {
        int length=l.dims[2*a+1]-l.dims[2*a]+1;
        bool tiled=a==0 || length>1 || (l.repeating && (a!=2 || l.pitch[a]>0));
        float shifted=adjacent_float(position[a],d[a]);
        float index=tiled?floor((shifted-l.lower[a])/l.pitch[a]):0;
        if(!isfinite(index) || abs(index)>2147483000.0f){error=true;return false;}
        relative[a]=int(index);
        if(l.pitch[a]>0) {
            low[a]=l.lower[a]+float(relative[a])*l.pitch[a];high[a]=low[a]+l.pitch[a];
            if(tiled && d[a]!=0) {
                /* Select the element after an exact ray/grid crossing using
                 * the same t calculation as the DDA. Reconstructing o+t*d
                 * alone may round back into the previous element. */
                float leave=((d[a]>0?high[a]:low[a])-o[a])/d[a];
                if(leave<=begin) {
                    relative[a]+=d[a]>0?1:-1;
                    low[a]=l.lower[a]+float(relative[a])*l.pitch[a];high[a]=low[a]+l.pitch[a];
                }
            }
        } else {low[a]=-INFINITY;high[a]=INFINITY;}
        if((l.dims[2*a]>0 && relative[a]>2147483647-l.dims[2*a]) ||
           (l.dims[2*a]<0 && relative[a]<(-2147483647-1)-l.dims[2*a])){error=true;return false;}
        logical[a]=relative[a]+l.dims[2*a];
        translation[a]=l.zero_coords?float(logical[a])*l.pitch[a]:
            a==2 && !tiled?0:l.lower[a]+(float(relative[a])+.5f)*l.pitch[a];
    }
    bool outside=false;
    for(uint a=0;a<3;++a)outside|=logical[a]<l.dims[2*a] || logical[a]>l.dims[2*a+1];
    if(l.repeating)child=fills[l.offset];
    else if(outside)child=l.outer;
    else {
        uint nj=uint(l.dims[3]-l.dims[2]+1),nk=uint(l.dims[5]-l.dims[4]+1);
        child=fills[l.offset+(uint(relative.x)*nj+uint(relative.y))*nk+uint(relative.z)];
    }
    next=end;normal=float3(0);
    float latest=-INFINITY;
    for(uint a=0;a<3;++a)if(d[a]!=0) {
        float entry=((d[a]>0?low[a]:high[a])-o[a])/d[a];
        float leave=((d[a]>0?high[a]:low[a])-o[a])/d[a];
        if(entry>latest){latest=entry;normal=float3(0);normal[a]=d[a]>0?-1:1;}
        if(leave>begin)next=min(next,leave);
        else if(isfinite(leave)) {
            /* Rounding at large coordinates can put the reconstructed point
             * one ULP behind a grid face; never silently skip an element. */
            error=true;return false;
        }
    }
    return next>begin;
}

static compute_render_pixel trace_pixel(device const compute_node* code,
        device const compute_primitive* primitives,device const compute_cell* cells,
        device const compute_render_bvh* bvh,device const compute_transform* transforms,
        device const compute_lattice* lattices,device const uint* fills,
        device const float* colors,constant compute_render_params& p,float2 pixel,
        thread compute_profile& cost,device compute_profile* cell_costs,
        device const compute_spatial_program* spatial) {
    float3 world_o(p.eye[0],p.eye[1],p.eye[2]),world_d(p.forward[0],p.forward[1],p.forward[2]);
    float3 right(p.right[0],p.right[1],p.right[2]),up(p.up[0],p.up[1],p.up[2]);
    float u=(2*pixel.x/p.width-1)*p.half_width,v=(1-2*pixel.y/p.height)*p.half_height;
    if(p.perspective)world_d=normalize(world_d+u*right+v*up);else world_o+=u*right+v*up;
    compute_render_pixel out={};out.cell_id=-1;out.depth=1e30f;
    for(uint a=0;a<3;++a)out.color[a]=p.background[a];
    if(!all(isfinite(world_o)) || !all(isfinite(world_d))){out.error=1;return out;}
    if(p.bvh_root==COMPUTE_NONE)return out;
    float begin=0,end=1e15f;float3 clip(p.clip[0],p.clip[1],p.clip[2]);
    if(p.clip_enabled) {
        float den=dot(clip,world_d),num=dot(clip,world_o)+p.clip[3];
        if(!isfinite(den) || !isfinite(num)){out.error=1;return out;}
        if(den==0){if(num<=0)return out;}
        else if(den>0)begin=max(begin,-num/den);else end=min(end,-num/den);
        if(end<=begin)return out;
    }
    render_frame stack[COMPUTE_HIER_STACK];int top=0;
    uint deferred[COMPUTE_BVH_STACK],deferred_count=0;
    stack[0]={world_o,world_d,float3(0),float3x3(1),begin,end,begin,end,p.bvh_root,0,COMPUTE_NONE};
    stack[0].guard=COMPUTE_NONE;
    float best=end;uint chosen=COMPUTE_NONE;float3 chosen_normal(0);bool error=false;
    uint steps=0;
    while(top>=0) {
        if(profile_enabled)cost.hierarchy=profile_sum(cost.hierarchy,1);
        if(++steps>1000000){out.error=2;return out;}
        thread render_frame& f=stack[top];
        float limit=min(best,f.end);
        if(f.begin>=limit){deferred_count=f.deferred_start;--top;continue;}
        if(f.pending!=COMPUTE_NONE) {
            limit=min(limit,f.pending_end);
            uint id=f.pending;compute_cell c=cells[id];
            if(f.cursor>=limit){f.pending=COMPUTE_NONE;continue;}
            float lo=f.cursor,hi=limit;uint surface=COMPUTE_NONE,child=c.child;
            float3 normal=f.entry_normal,translation(0);bool lattice=c.lattice!=COMPUTE_NONE;
            bool lazy=!lattice && defer_parent(c);
            if(lattice) {
                float3 local_normal;
                if(!lattice_segment(lattices[c.lattice],fills,f.o,f.d,lo,limit,hi,child,translation,local_normal,error)) {
                    if(error){out.error=3;return out;}f.pending=COMPUTE_NONE;continue;
                }
                if(lo>f.begin)normal=to_world_normal(f,local_normal);
            } else if(!lazy) {
                if(profile_enabled)cost.parents=profile_sum(cost.parents,1);
                if(!cell_interval(code,primitives,c,f.o,f.d,lo,limit,true,lo,hi,surface,error,cost,cell_costs,id,spatial)) {
                    if(error){out.error=1;return out;}f.pending=COMPUTE_NONE;continue;
                }
                if(surface!=COMPUTE_NONE)normal=to_world_normal(f,surface_normal(primitives[surface],f.o+lo*f.d));
            }
            f.cursor=hi;
            if(child==COMPUTE_NONE)continue;
            if(top+1>=COMPUTE_HIER_STACK){out.error=4;return out;}
            render_frame next={f.o,f.d,normal,f.world_to_local,lo,hi,lo,hi,child,deferred_count,COMPUTE_NONE};
            next.guard=lazy?id:COMPUTE_NONE;
            next.entry_owner=(lattice?lo>f.begin:surface!=COMPUTE_NONE)?uint(top+1):f.entry_owner;
            next.guard_surface=next.guard_boundary=COMPUTE_NONE;
            if(lattice)next.o-=translation;
            else {
                device const float* m=transforms[c.transform].inverse;
                float3x3 linear(float3(m[0],m[4],m[8]),float3(m[1],m[5],m[9]),float3(m[2],m[6],m[10]));
                next.o=linear*f.o+float3(m[3],m[7],m[11]);next.d=linear*f.d;
                next.world_to_local=linear*f.world_to_local;
            }
            if(!all(isfinite(next.o)) || !all(isfinite(next.d)) ||
               !all(isfinite(next.world_to_local[0])) || !all(isfinite(next.world_to_local[1])) ||
               !all(isfinite(next.world_to_local[2]))) {out.error=1;return out;}
            stack[++top]=next;continue;
        }
        if(f.node==COMPUTE_NONE) {
            if(deferred_count==f.deferred_start){--top;continue;}
            f.node=deferred[--deferred_count];
        }
        if(profile_enabled)cost.bvh=profile_sum(cost.bvh,1);
        uint node=f.node;device const compute_render_bvh& b=bvh[node];
        f.node=COMPUTE_NONE;
        float boxlo=f.begin,boxhi=limit;
        if(!render_box(b,f.o,f.d,boxlo,boxhi))continue;
        if(b.cell==COMPUTE_NONE) {
            uint left=node+1,right=bvh[left].escape;
            float ll=f.begin,lh=limit,rl=f.begin,rh=limit;
            bool hit_left=render_box(bvh[left],f.o,f.d,ll,lh),hit_right=render_box(bvh[right],f.o,f.d,rl,rh);
            if(hit_left && hit_right) {
                if(deferred_count==COMPUTE_BVH_STACK){out.error=5;return out;}
                f.node=ll<=rl?left:right;deferred[deferred_count++]=ll<=rl?right:left;
            } else if(hit_left)f.node=left;else if(hit_right)f.node=right;
            continue;
        }
        uint id=b.cell;compute_cell c=cells[id];
        if(c.lattice==COMPUTE_NONE && c.transform!=COMPUTE_NONE && c.child==COMPUTE_NONE)continue;
        if(c.lattice!=COMPUTE_NONE && lattices[c.lattice].repeating && fills[lattices[c.lattice].offset]==COMPUTE_NONE)continue;
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
                if(miss || boxlo>=boxhi)continue;
            }
        }
        if(c.transform!=COMPUTE_NONE || c.lattice!=COMPUTE_NONE) {
            f.pending=id;
            float padding=max(.001f,max(abs(boxlo),abs(boxhi))*1e-5f);
            f.pending_end=c.lattice!=COMPUTE_NONE?boxhi:min(limit,boxhi+padding);
            /* Only a lattice grid box is itself a real boundary. Ordinary
             * cell boxes are conservative and must not replace CSG crossings. */
            f.cursor=c.lattice!=COMPUTE_NONE?boxlo:max(f.begin,boxlo-padding);
            continue;
        }
        float lo,hi;uint surface;
        /* Search only the conservative box interval, with padding so the box
         * never substitutes for a physical surface or its shading normal. */
        float padding=max(.001f,max(abs(boxlo),abs(boxhi))*1e-5f);
        float cell_begin=max(f.begin,boxlo-padding),cell_end=min(limit,boxhi+padding);
        while(cell_begin<cell_end && cell_interval(code,primitives,c,f.o,f.d,cell_begin,cell_end,false,lo,hi,surface,error,cost,cell_costs,id,spatial)) {
            float3 normal=surface==COMPUTE_NONE?f.entry_normal:
                to_world_normal(f,surface_normal(primitives[surface],f.o+lo*f.d));
            bool accepted=true;
            if(lazy_threshold)for(int a=1;a<=top;++a) {
                thread render_frame& guard=stack[a];if(guard.guard==COMPUTE_NONE)continue;
                if(!parent_truth(code,primitives,cells,guard,stack[a-1],lo,cell_end,error,cost,cell_costs,spatial)) {
                    if(profile_enabled)cost.rejected=profile_sum(cost.rejected,1);
                    cell_begin=guard.guard_hi;accepted=false;break;
                }
                if(surface==COMPUTE_NONE && lo==guard.guard_lo && guard.guard_surface!=COMPUTE_NONE &&
                   (lo!=f.begin || uint(a)>=f.entry_owner))
                    normal=to_world_normal(stack[a-1],surface_normal(primitives[guard.guard_surface],stack[a-1].o+lo*stack[a-1].d));
            }
            if(error)break;
            if(accepted){if(lo<best){best=lo;chosen=id;chosen_normal=normal;}break;}
        }
        if(error){out.error=1;return out;}
    }
    if(chosen==COMPUTE_NONE)return out;
    float3 normal=chosen_normal,base(colors[3*chosen],colors[3*chosen+1],colors[3*chosen+2]);
    bool cut=p.clip_enabled && begin>0 && best==begin;
    if(cut)normal=-normalize(clip);if(dot(normal,world_d)>0)normal=-normal;
    float3 color;
    if(cut)color=base*p.tint;
    else {
        float3 light=p.light_set?float3(p.light[0],p.light[1],p.light[2]):world_d;
        float ndl=max(-dot(normal,light),0.0f);
        float shine=pow(max(-dot(2*ndl*normal+light,world_d),0.0f),p.shininess);
        color=base*(p.ambient+p.diffuse*ndl)+p.specular*shine;
    }
    if(!all(isfinite(color)) || !all(isfinite(normal))){out.error=1;return out;}
    color=clamp(color,0.0f,1.0f);
    for(uint a=0;a<3;++a){out.color[a]=color[a];out.normal[a]=normal[a];}
    out.depth=best;out.cell_id=cells[chosen].cell_id;out.material_id=cells[chosen].material_id;return out;
}

kernel void alea_render_image(device const compute_node* code [[buffer(0)]],
        device const compute_cell* cells [[buffer(1)]],device const compute_render_bvh* bvh [[buffer(2)]],
        device const float* colors [[buffer(3)]],device compute_render_pixel* pixels [[buffer(4)]],
        constant compute_render_params& p [[buffer(5)]],device const compute_primitive* primitives [[buffer(6)]],
        device const compute_transform* transforms [[buffer(7)]],device const compute_lattice* lattices [[buffer(8)]],
        device const uint* fills [[buffer(9)]],constant compute_render_tile& tile [[buffer(10)]],
        device compute_profile* ray_costs [[buffer(11)]],
        device compute_profile* cell_costs [[buffer(12)]],
        device const compute_spatial_program* spatial [[buffer(15)]],
        uint index [[thread_position_in_grid]]) {
    index+=tile.offset;
    index=(tile.y+index/tile.width)*p.width+tile.x+index%tile.width;
    if(index>=p.width*p.height)return;
    compute_profile cost={};
    float2 origin(index%p.width,index/p.width);compute_render_pixel output={};float3 sum(0);int error=0;
    for(uint y=0;y<p.aa;++y)for(uint x=0;x<p.aa;++x) {
        compute_render_pixel sample=trace_pixel(code,primitives,cells,bvh,transforms,lattices,fills,colors,p,
            origin+(float2(x,y)+.5f)/float(p.aa),cost,cell_costs,spatial);
        sum+=float3(sample.color[0],sample.color[1],sample.color[2]);error=max(error,sample.error);
        if(x==p.aa/2 && y==p.aa/2)output=sample;
    }
    for(uint a=0;a<3;++a)output.color[a]=sum[a]/float(p.aa*p.aa);
    output.error=error;pixels[index]=output;
    if(profile_enabled)ray_costs[index]=cost;
}

#include "resume.metal"
