// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* Ordinary fills propagate a finite working region into their child's local
 * coordinates. A query beyond a working region can always use the BLAS. */
static alea_bbox_t hierarchy_box(interval_box_t b) {
    return (alea_bbox_t){.min_x=b.lo[0],.max_x=b.hi[0],.min_y=b.lo[1],
                         .max_y=b.hi[1],.min_z=b.lo[2],.max_z=b.hi[2]};
}

static bool merge_fill_box(interval_box_t* region, interval_box_t parent,
                           const alea_matrix_t* transform) {
    interval_box_t child={{INFINITY,INFINITY,INFINITY},{-INFINITY,-INFINITY,-INFINITY}};
    for (unsigned i=0;i<8;i++) {
        double p[]={i&1?parent.hi[0]:parent.lo[0],i&2?parent.hi[1]:parent.lo[1],
                    i&4?parent.hi[2]:parent.lo[2]};
        alea_matrix_transform_point_inverse(transform,p,p+1,p+2);
        for(unsigned k=0;k<3;k++) {
            if(!isfinite(p[k]))return false;
            if(p[k]<child.lo[k])child.lo[k]=p[k];
            if(p[k]>child.hi[k])child.hi[k]=p[k];
        }
    }
    bool changed=false;
    for(unsigned k=0;k<3;k++) {
        if(child.lo[k]<region->lo[k]){region->lo[k]=child.lo[k];changed=true;}
        if(child.hi[k]>region->hi[k]){region->hi[k]=child.hi[k];changed=true;}
    }
    return changed;
}

static int choose_domains(interval_index_t* idx) {
    const alea_system_t* sys=idx->sys;
    const alea_universe_t* world=alea_get_universe(sys,0);
    size_t n=sys->universes.count;
    interval_box_t* reach=calloc(n+1,sizeof(*reach));
    bool* changed=calloc(n+1,sizeof(*changed));
    if(!reach || !changed){free(reach);free(changed);return -1;}
    for(size_t i=0;i<n;i++) {
        idx->domains[i]=idx->domain;
        reach[i]=(interval_box_t){{INFINITY,INFINITY,INFINITY},{-INFINITY,-INFINITY,-INFINITY}};
    }
    if(world) {
        size_t root=world-sys->universes.data;
        reach[root]=interval_box_from(idx->domain);
        changed[root]=true;
        for(unsigned pass=0;pass<MAX_FILL_RAYCAST_DEPTH;pass++) {
            bool work=false;
            for(size_t i=0;i<n;i++) {
                if(!changed[i])continue;
                changed[i]=false;
                const alea_universe_t* u=sys->universes.data+i;
                for(size_t j=0;j<u->cell_indices.count;j++) {
                    uint32_t id=u->cell_indices.data[j];
                    const alea_cell_entry_t* c=sys->cells.data+id;
                    if(c->fill_universe<=0 || c->lat_type)continue;
                    const alea_universe_t* target=alea_get_universe(sys,c->fill_universe);
                    const alea_matrix_t* matrix=alea_hier_spatial_get_cell_fill_matrix(sys,id);
                    if(!target || !matrix || !matrix->has_inverse)continue;
                    size_t child=target-sys->universes.data;
                    if(child==root || child==i)continue;
                    interval_box_t clip=bounded_cell_domain(idx,id,hierarchy_box(reach[i]));
                    if(!box_empty(clip) && merge_fill_box(reach+child,clip,matrix)) {
                        changed[child]=true;
                        work=true;
                    }
                }
            }
            if(!work)break;
        }
        for(size_t i=0;i<n;i++) {
            if(i==root || box_empty(reach[i]))continue;
            bool finite=true;
            for(unsigned k=0;k<3;k++) {
                double pad=1e-6*(1+fabs(reach[i].lo[k])+fabs(reach[i].hi[k]))+
                           0.005*(reach[i].hi[k]-reach[i].lo[k]);
                reach[i].lo[k]-=pad;reach[i].hi[k]+=pad;
                if(!isfinite(reach[i].hi[k]-reach[i].lo[k]))finite=false;
            }
            if(finite)idx->domains[i]=hierarchy_box(reach[i]);
        }
    }
    free(reach);free(changed);return 0;
}
