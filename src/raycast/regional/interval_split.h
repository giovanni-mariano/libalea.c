// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* Score interior cuts using both the volume crossed by points and the
 * surface area seen by rays. Sampling affects tree shape only. Each child
 * still receives every cell that passes the conservative overlap test. */
static double select_split(const interval_index_t* idx,alea_bbox_t box,const uint32_t* ids,
                           size_t count,unsigned* axis) {
    interval_box_t b=interval_box_from(box);
    double width[3],longest=0;
    *axis=0;
    for(unsigned k=0;k<3;k++) {
        width[k]=b.hi[k]-b.lo[k];
        if(width[k]>longest){longest=width[k];*axis=k;}
    }
    double split=.5*b.lo[*axis]+.5*b.hi[*axis];
    if(count<24 || !isfinite(longest) || longest<=0)return split;
    double x=width[0]/longest,y=width[1]/longest,z=width[2]/longest;
    double area=x*y+y*z+z*x;
    if(!(area>0))return split;
    double face[]={y*z,z*x,x*y};
    size_t samples=count<64?count:64;
    double best=INFINITY;
    for(unsigned k=0;k<3;k++) {
        if(width[k]<=0)continue;
        unsigned starts[9]={0},ends[9]={0};
        double scale=8/width[k];
        if(!isfinite(scale))continue;
        for(size_t j=0;j<samples;j++) {
            const interval_box_t* bound=idx->bounds+ids[j*(count-1)/(samples-1)];
            double lo=(bound->lo[k]-b.lo[k])*scale,hi=(bound->hi[k]-b.lo[k])*scale;
            unsigned first=lo<=0?0:lo>=8?8:(unsigned)ceil(lo);
            unsigned last=hi<=0?0:hi>=8?8:(unsigned)hi;
            starts[first]++;ends[last]++;
        }
        unsigned left=starts[0],right=(unsigned)samples;
        double cross=face[k]/area;
        for(unsigned cut=1;cut<8;cut++) {
            left+=starts[cut];right-=ends[cut-1];
            double f=cut/8.0;
            double weight=.5*cross;
            double cost=left*(f+weight*(1-f))+right*((1-f)+weight*f);
            /* Resolve near-ties towards a balanced cut on a longer axis. */
            cost+=.02*samples*(fabs(f-.5)+.25*(1-width[k]/longest));
            if(cost<best){best=cost;*axis=k;split=b.lo[k]+f*width[k];}
        }
    }
    return split;
}
