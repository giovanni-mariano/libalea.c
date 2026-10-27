// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* Analytic surface crossings. Ray directions need not be normalized: affine
 * fill transforms preserve the world ray parameter throughout the hierarchy. */
static bool ray_box(device const float* b,float3 o,float3 d,thread float& lo,thread float& hi) {
    for(uint a=0;a<3;++a) {
        if(b[2*a]>b[2*a+1])return false;
        if(d[a]==0){if(o[a]<b[2*a] || o[a]>b[2*a+1])return false;}
        else {float x=(b[2*a]-o[a])/d[a],y=(b[2*a+1]-o[a])/d[a];
            lo=max(lo,min(x,y));hi=min(hi,max(x,y));if(lo>hi)return false;}
    }
    return true;
}
static uint quadratic(float a,float b,float c,thread float4& roots,thread bool& error) {
    float scale=max(abs(a),max(abs(b),abs(c)));
    if(!isfinite(scale)){error=true;return 0;}
    if(scale==0)return 0;
    a/=scale;b/=scale;c/=scale;
    if(a==0){if(b==0)return 0;roots[0]=-c/b;return isfinite(roots[0])?1:0;}
    float disc=fma(b,b,-4*a*c);
    if(disc<0)return 0;
    float q=-.5f*(b+copysign(sqrt(disc),b));
    roots[0]=q==0?-b/(2*a):q/a;roots[1]=q==0?roots[0]:c/q;
    if(roots[0]>roots[1]){float t=roots[0];roots[0]=roots[1];roots[1]=t;}
    return 2;
}
static uint cubic_stationary(float a,float b,float c,thread float3& r) {
    /* Monic cubic x^3+a*x^2+b*x+c. Used only on a centered, scaled torus. */
    float p=b-a*a/3, q=(2*a*a*a)/27-a*b/3+c;
    float disc=q*q/4+p*p*p/27;
    if(disc>0) {
        float d=sqrt(disc),x=-q/2+d,y=-q/2-d;
        r[0]=copysign(pow(abs(x),1.0f/3),x)+copysign(pow(abs(y),1.0f/3),y)-a/3;return 1;
    }
    if(p>=0){r[0]=-a/3;return 1;}
    float radius=2*sqrt(-p/3);
    float angle=acos(clamp((3*q/(2*p))*sqrt(-3/p),-1.0f,1.0f))/3;
    r[0]=radius*cos(angle)-a/3;
    r[1]=radius*cos(angle-2.09439510239f)-a/3;
    r[2]=radius*cos(angle-4.18879020479f)-a/3;
    if(r[0]>r[1]){float t=r[0];r[0]=r[1];r[1]=t;}
    if(r[1]>r[2]){float t=r[1];r[1]=r[2];r[2]=t;}
    if(r[0]>r[1]){float t=r[0];r[0]=r[1];r[1]=t;}
    return 3;
}
static float torus_distance(float3 p,float R,float r) {
    float radial=length(p.xy)-R;
    return radial*radial+p.z*p.z-r*r;
}
static uint torus_roots(device const compute_primitive& n,float3 o,float3 d,
                         thread float4& roots,thread bool& error) {
    device const float* p=n.p;
    o-=float3(p[0],p[1],p[2]);
    uint axis=n.type-COMPUTE_TORUS_X;
    if(axis==0){o=o.zyx;d=d.zyx;}else if(axis==1){o=o.xzy;d=d.xzy;}
    o.z*=p[4]/p[5];d.z*=p[4]/p[5];
    float speed=length(d),scale=p[3]+p[4];
    if(!isfinite(speed) || !isfinite(scale) || speed==0 || scale==0){error=true;return 0;}
    d/=speed;
    float shift=-dot(o,d);
    o=(o+shift*d)/scale;
    float R=p[3]/scale,r=p[4]/scale;
    float oo=dot(o,o),od=dot(o,d),dd=dot(d,d);
    float sphere_disc=od*od+dd*(1-oo);
    if(!isfinite(sphere_disc)){error=true;return 0;}
    if(sphere_disc<0)return 0;
    float extent=sqrt(sphere_disc);
    float lower=(-od-extent)/dd,upper=(-od+extent)/dd;
    float k=oo+R*R-r*r;
    float c4=dd*dd,c3=4*dd*od;
    float c2=2*dd*k+4*od*od-4*R*R*dot(d.xy,d.xy);
    float c1=4*od*k-8*R*R*dot(o.xy,d.xy);
    float3 critical;
    uint nc=cubic_stationary(3*c3/(4*c4),c2/(2*c4),c1/(4*c4),critical);
    float ends[5];uint count=0;ends[count++]=lower;
    for(uint j=0;j<nc;++j)if(critical[j]>lower && critical[j]<upper)ends[count++]=critical[j];
    ends[count++]=upper;
    uint nr=0;
    for(uint j=0;j+1<count && nr<4;++j) {
        float a=ends[j],b=ends[j+1];
        float fa=torus_distance(o+a*d,R,r),fb=torus_distance(o+b*d,R,r);
        if(!isfinite(fa) || !isfinite(fb)){error=true;return 0;}
        if(fa==0){float t=(shift+a*scale)/speed;if(!nr || roots[nr-1]!=t)roots[nr++]=t;}
        if(nr==4)break;
        if((fa<0 && fb>0) || (fa>0 && fb<0)) {
            /* Quartic critical points isolate the crossings. Evaluate the
             * unsquared torus equation during bisection to retain thin tubes
             * and reject the extraneous branch of spindle tori. */
            for(uint it=0;it<30;++it) {
                float mid=a+(b-a)*.5f;if(mid==a || mid==b)break;
                float fm=torus_distance(o+mid*d,R,r);
                if((fm<0)==(fa<0)){a=mid;fa=fm;}else b=mid;
            }
            roots[nr++]=(shift+(a+(b-a)*.5f)*scale)/speed;
        }
        if(j+2==count && fb==0 && nr<4)roots[nr++]=(shift+ends[j+1]*scale)/speed;
    }
    return nr;
}

static uint primitive_roots(device const compute_primitive& n,float3 o,float3 d,
                             thread float4& t,thread bool& error) {
    device const float* p=n.p;
    if(n.type>=COMPUTE_TORUS_X && n.type<=COMPUTE_TORUS_Z)return torus_roots(n,o,d,t,error);
    if(n.type==COMPUTE_PLANE) {
        float3 normal(p[0],p[1],p[2]);float den=dot(normal,d),num=dot(normal,o)+p[3];
        if(!isfinite(den) || !isfinite(num)){error=true;return 0;}
        if(den==0)return 0;t[0]=-num/den;return isfinite(t[0])?1:0;
    }
    if(n.type==COMPUTE_BOX) {
        float lo=-INFINITY,hi=INFINITY;bool hit=ray_box(p,o,d,lo,hi);t[0]=lo;t[1]=hi;return hit?2:0;
    }
    if(n.type==COMPUTE_QUADRIC) {
        float3 q=o-float3(p[10],p[11],p[12]);
        /* Center the ray at the quadric's local origin before substitution. */
        float shift=-dot(q,d)/dot(d,d);q+=shift*d;
        float a=d.x*(p[0]*d.x+p[3]*d.y+p[5]*d.z)+d.y*(p[1]*d.y+p[4]*d.z)+p[2]*d.z*d.z;
        float b=2*(p[0]*q.x*d.x+p[1]*q.y*d.y+p[2]*q.z*d.z)+p[3]*(q.x*d.y+q.y*d.x)+
            p[4]*(q.y*d.z+q.z*d.y)+p[5]*(q.x*d.z+q.z*d.x)+p[6]*d.x+p[7]*d.y+p[8]*d.z;
        float c=q.x*(p[0]*q.x+p[3]*q.y+p[5]*q.z+p[6])+q.y*(p[1]*q.y+p[4]*q.z+p[7])+q.z*(p[2]*q.z+p[8])+p[9];
        uint count=quadratic(a,b,c,t,error);for(uint j=0;j<count;++j)t[j]+=shift;return count;
    }
    if(n.type>=COMPUTE_CONE_X && n.type<=COMPUTE_CONE_Z) {
        float3 q=o-float3(p[0],p[1],p[2]);uint axis=n.type-COMPUTE_CONE_X;
        float shift=-dot(q,d)/dot(d,d);q+=shift*d;
        float qa=q[axis],da=d[axis];q[axis]=0;float3 rd=d;rd[axis]=0;
        uint count=quadratic(dot(rd,rd)-p[3]*da*da,2*(dot(q,rd)-p[3]*qa*da),dot(q,q)-p[3]*qa*qa,t,error);
        uint nr=0;for(uint j=0;j<count;++j)if(p[4]==0 || p[4]*(qa+t[j]*da)>=0)t[nr++]=t[j]+shift;
        return nr;
    }
    if(n.type==COMPUTE_RCC) {
        float3 q=o-float3(p[0],p[1],p[2]),axis(p[3],p[4],p[5]);
        float h=dot(q,axis),v=dot(d,axis);float3 radial=q-h*axis,rd=d-v*axis;
        float4 side;uint count=0,nr=0;float a=dot(rd,rd);
        if(a>0) {
            float center=-dot(radial,rd)/a;float3 closest=radial+center*rd;
            float disc=(p[7]*p[7]-dot(closest,closest))/a;
            if(!isfinite(disc)){error=true;return 0;}
            if(disc>=0){float span=sqrt(disc);side[0]=center-span;side[1]=center+span;count=2;}
        }
        for(uint j=0;j<count;++j)if(h+side[j]*v>=0 && h+side[j]*v<=p[6])t[nr++]=side[j];
        if(v!=0)for(uint j=0;j<2;++j){float at=(j*p[6]-h)/v;float3 r=radial+at*rd;if(dot(r,r)<=p[7]*p[7] && nr<4)t[nr++]=at;}
        return nr;
    }
    float3 q;float radius;
    if(n.type==COMPUTE_SPHERE){q=o-float3(p[0],p[1],p[2]);radius=p[3];}
    else {
        uint axis=n.type-COMPUTE_CYL_X;
        q=axis==0?float3(0,o.y-p[0],o.z-p[1]):axis==1?float3(o.x-p[0],0,o.z-p[1]):float3(o.x-p[0],o.y-p[1],0);
        d[axis]=0;radius=p[2];
    }
    float a=dot(d,d),b=dot(q,d);if(!isfinite(a) || !isfinite(b)){error=true;return 0;}if(a==0)return 0;
    float center=-b/a;float3 closest=q+center*d;float disc=(radius*radius-dot(closest,closest))/a;
    if(!isfinite(disc)){error=true;return 0;}if(disc<0)return 0;
    float span=sqrt(disc);t[0]=center-span;t[1]=center+span;return 2;
}

static float3 surface_normal(device const compute_primitive& n,float3 q) {
    device const float* p=n.p;float3 v(0);
    switch(n.type) {
    case COMPUTE_PLANE:v=float3(p[0],p[1],p[2]);break;
    case COMPUTE_SPHERE:v=q-float3(p[0],p[1],p[2]);break;
    case COMPUTE_CYL_X:v=float3(0,q.y-p[0],q.z-p[1]);break;
    case COMPUTE_CYL_Y:v=float3(q.x-p[0],0,q.z-p[1]);break;
    case COMPUTE_CYL_Z:v=float3(q.x-p[0],q.y-p[1],0);break;
    case COMPUTE_BOX: {
        float closest=INFINITY;
        for(uint i=0;i<6;++i)if(abs(q[i/2]-p[i])<closest){closest=abs(q[i/2]-p[i]);v=float3(0);v[i/2]=(i&1)?1:-1;}
        break;
    }
    case COMPUTE_QUADRIC:
        q-=float3(p[10],p[11],p[12]);
        v=float3(2*p[0]*q.x+p[3]*q.y+p[5]*q.z+p[6],2*p[1]*q.y+p[3]*q.x+p[4]*q.z+p[7],2*p[2]*q.z+p[4]*q.y+p[5]*q.x+p[8]);break;
    case COMPUTE_CONE_X:case COMPUTE_CONE_Y:case COMPUTE_CONE_Z:
        v=q-float3(p[0],p[1],p[2]);v[n.type-COMPUTE_CONE_X]*=-p[3];break;
    case COMPUTE_TORUS_X:case COMPUTE_TORUS_Y:case COMPUTE_TORUS_Z: {
        q-=float3(p[0],p[1],p[2]);uint axis=n.type-COMPUTE_TORUS_X;float axial=q[axis];q[axis]=0;
        float rho=length(q);v=rho>0?q*((rho-p[3])/rho):float3(0);
        v[axis]=axial*(p[4]/p[5])*(p[4]/p[5]);break;
    }
    case COMPUTE_RCC: {
        q-=float3(p[0],p[1],p[2]);float3 axis(p[3],p[4],p[5]);float h=dot(q,axis);v=q-h*axis;
        float side=abs(length(v)-p[7]);if(abs(h)<side)v=-axis;else if(abs(h-p[6])<side)v=axis;break;
    }
    }
    return length_squared(v)>0?normalize(v):float3(0);
}
