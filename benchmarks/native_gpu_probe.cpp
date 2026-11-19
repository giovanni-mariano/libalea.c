// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
// Standalone feasibility probe. It does not register a libalea compute backend.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <cmath>
#include <vector>
#include <stdexcept>
#if defined(PROBE_CUDA)
#include <cuda_runtime.h>
#define GPU(name) cuda##name
#define BACKEND "cuda"
#elif defined(PROBE_HIP)
#include <hip/hip_runtime.h>
#define GPU(name) hip##name
#define BACKEND "hip"
#else
#define BACKEND "host-validation"
#endif
#include "native_probe_kernel.h"

static double now() {
    return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
#if defined(PROBE_CUDA) || defined(PROBE_HIP)
static void check(GPU(Error_t) result) {
    if (result != GPU(Success)) throw std::runtime_error(GPU(GetErrorString)(result));
}
template<class T> struct Buffer {
    T *ptr = nullptr;
    explicit Buffer(size_t n) { check(GPU(Malloc)((void**)&ptr,n*sizeof(T))); }
    ~Buffer() { if (ptr) GPU(Free)(ptr); }
    Buffer(const Buffer&) = delete;
    Buffer& operator=(const Buffer&) = delete;
};
__global__ static void points_kernel(const ProbePrimitive *p, uint32_t n,
                                    const ProbePoint *q, int *out, uint32_t count) {
    uint32_t i=blockIdx.x*blockDim.x+threadIdx.x;
    if (i<count) out[i]=probe_contains(p,n,q[i]);
}
__global__ static void rays_kernel(const ProbePrimitive *p, uint32_t n,
                                  const ProbeRay *r, ProbeState *s,
                                  uint32_t count, uint32_t budget) {
    uint32_t i=blockIdx.x*blockDim.x+threadIdx.x;
    if (i<count) probe_slice(p,n,r[i],s+i,budget);
}
#endif

static void json_string(const char *s) {
    putchar('"');
    for (; *s; ++s) {
        unsigned char c=(unsigned char)*s;
        if (c=='"' || c=='\\') printf("\\%c",c);
        else if (c<32) printf("\\u%04x",c);
        else putchar(c);
    }
    putchar('"');
}

/* Independent analytic expectations for sphere radius 2 clipped to -1<=x<=1.
 * Rays use t>=0; tangencies count as singleton intervals in THIS probe. */
static const ProbeRay rays[] = {
    {{-3,0,0},{1,0,0}}, {{0,0,0},{1,0,0}}, {{0,3,0},{0,-1,0}},
    {{2,3,0},{0,-1,0}}, {{-3,0,0},{-1,0,0}}, {{0,2,0},{1,0,0}},
    {{0,0,0},{0,0,0}}, {{0,0,-3},{0,0,2}}
};
static const int expected_status[] = {1,1,1,2,2,1,3,1};
static const float expected_enter[] = {2,0,1,0,0,0,0,.5f};
static const float expected_exit[] = {4,1,5,0,0,0,0,2.5f};
static const ProbePoint points[] = {{0,0,0},{1,0,0},{1.25f,0,0},{0,2,0},
                                   {0,2.25f,0},{-1,0,0},{-1.25f,0,0},{0,0,-1}};
static const int expected_point[] = {1,1,0,1,0,1,0,1};

static void continuation_checks() {
    const ProbePrimitive p[]={{0,0,0,0,2},{1,1,0,0,-1},{1,-1,0,0,-1}};
    ProbeState s={0,INFINITY,0,PROBE_ACTIVE};
    probe_slice(p,3,rays[0],&s,0);
    if (s.cursor!=0 || s.status!=PROBE_ACTIVE || s.enter!=0 || !std::isinf(s.exit))
        throw std::runtime_error("zero budget changed continuation");
    probe_slice(p,3,rays[0],&s,1);
    if (s.cursor!=1 || s.status!=PROBE_ACTIVE || s.enter!=1 || s.exit!=5)
        throw std::runtime_error("slice did not stop at work budget");
    // Resume after an arbitrary pause, then verify completed work stays stable.
    probe_slice(p,3,rays[0],&s,2);
    probe_slice(p,3,rays[0],&s,2);
    if (s.cursor!=3 || s.status!=PROBE_HIT || s.enter!=2 || s.exit!=4)
        throw std::runtime_error("resume/completed continuation mismatch");
}

static void run(uint32_t count, uint32_t budget, uint32_t repeats) {
    double begin=now();
    std::vector<ProbePrimitive> p;
    p.push_back({0,0,0,0,2});
    for (uint32_t i=0;i<repeats;++i) {
        p.push_back({1,1,0,0,-1}); p.push_back({1,-1,0,0,-1});
    }
    uint32_t n=(uint32_t)p.size(), slices=(n+budget-1)/budget;
    std::vector<ProbePoint> q(count);
    std::vector<ProbeRay> r(count);
    std::vector<ProbeState> s(count, {0,INFINITY,0,PROBE_ACTIVE});
    std::vector<int> out(count);
    for (uint32_t i=0;i<count;++i) { q[i]=points[i%8]; r[i]=rays[i%8]; }
    double prepared=now(), uploaded=prepared, evaluated;
    size_t bytes=p.size()*sizeof(ProbePrimitive)+(size_t)count*
        (sizeof(ProbePoint)+sizeof(ProbeRay)+sizeof(ProbeState)+sizeof(int));
#if defined(PROBE_CUDA) || defined(PROBE_HIP)
    Buffer<ProbePrimitive> dp(n); Buffer<ProbePoint> dq(count);
    Buffer<ProbeRay> dr(count); Buffer<ProbeState> ds(count); Buffer<int> dout(count);
    check(GPU(Memcpy)(dp.ptr,p.data(),p.size()*sizeof(p[0]),GPU(MemcpyHostToDevice)));
    check(GPU(Memcpy)(dq.ptr,q.data(),q.size()*sizeof(q[0]),GPU(MemcpyHostToDevice)));
    check(GPU(Memcpy)(dr.ptr,r.data(),r.size()*sizeof(r[0]),GPU(MemcpyHostToDevice)));
    check(GPU(Memcpy)(ds.ptr,s.data(),s.size()*sizeof(s[0]),GPU(MemcpyHostToDevice)));
    uploaded=now();
    points_kernel<<<(count+127)/128,128>>>(dp.ptr,n,dq.ptr,dout.ptr,count);
    check(GPU(GetLastError)());
    for (uint32_t i=0;i<slices;++i) {
        rays_kernel<<<(count+127)/128,128>>>(dp.ptr,n,dr.ptr,ds.ptr,count,budget);
        check(GPU(GetLastError)());
    }
    check(GPU(DeviceSynchronize)()); evaluated=now();
    check(GPU(Memcpy)(out.data(),dout.ptr,out.size()*sizeof(out[0]),GPU(MemcpyDeviceToHost)));
    check(GPU(Memcpy)(s.data(),ds.ptr,s.size()*sizeof(s[0]),GPU(MemcpyDeviceToHost)));
#else
    for (uint32_t i=0;i<count;++i) out[i]=probe_contains(p.data(),n,q[i]);
    for (uint32_t j=0;j<slices;++j)
        for (uint32_t i=0;i<count;++i) probe_slice(p.data(),n,r[i],&s[i],budget);
    evaluated=now();
#endif
    double read=now();
    for (uint32_t i=0;i<count;++i) {
        uint32_t k=i%8;
        if (out[i]!=expected_point[k] || s[i].status!=expected_status[k] ||
            (s[i].status==PROBE_HIT &&
             (s[i].cursor!=n || !std::isfinite(s[i].enter) || !std::isfinite(s[i].exit) ||
              fabsf(s[i].enter-expected_enter[k])>1e-5f ||
              fabsf(s[i].exit-expected_exit[k])>1e-5f)))
            throw std::runtime_error("analytic point/ray fixture mismatch");
    }
    printf("{\"queries\":%u,\"primitives\":%u,\"slice_budget\":%u,\"slices\":%u,"
           "\"buffer_payload_bytes\":%zu,\"prepare_seconds\":%.9f,"
           "\"allocate_upload_seconds\":%.9f,\"submit_and_wait_seconds\":%.9f,"
           "\"readback_seconds\":%.9f,\"verified\":true}",
           count,n,budget,slices,bytes,prepared-begin,uploaded-prepared,evaluated-uploaded,read-evaluated);
}

int main(int argc, char **argv) {
    bool test=argc==2 && !strcmp(argv[1],"--self-test");
    if (argc>1 && !test) { fprintf(stderr,"Usage: %s [--self-test]\n",argv[0]); return 2; }
    try {
        continuation_checks();
        const char *device="CPU (no GPU execution)";
        double begin=now();
#if defined(PROBE_CUDA) || defined(PROBE_HIP)
        int devices=0; check(GPU(GetDeviceCount)(&devices));
        if (!devices) throw std::runtime_error("no compatible GPU available");
        check(GPU(SetDevice)(0));
#if defined(PROBE_CUDA)
        cudaDeviceProp props;
#else
        hipDeviceProp_t props;
#endif
        check(GPU(GetDeviceProperties)(&props,0)); device=props.name;
        check(GPU(Free)(nullptr)); // include lazy context initialization
#endif
        printf("{\"schema\":1,\"backend\":\"%s\",\"device\":",BACKEND);
        json_string(device);
        printf(",\"device_init_seconds\":%.9f,\"fixture\":\"clipped-sphere-intersection\","
               "\"timing\":\"host-wall-clock\",\"runs\":[",now()-begin);
        if (test) {
            const uint32_t budgets[]={1,7,128,8192};
            run(16,1,1);
            for (unsigned i=0;i<4;++i) { putchar(','); run(16,budgets[i],2049); }
        } else {
            for (unsigned i=0;i<4;++i) { if (i) putchar(','); run(8192,128,2049); }
        }
        puts("]}"); return 0;
    } catch (const std::exception &e) {
        fprintf(stderr,"%s probe failed: %s\n",BACKEND,e.what()); return 1;
    }
}
