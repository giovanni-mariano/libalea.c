// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
/* Exercise the packaged evaluator with independent diagnostic leaf traces. */
#define COMPUTE_RAY_GROUP_TRACE_LEAF(arg) if(ray_audit_active){ray_add(audit_cost.v[2]);audit_cost.trace=audit_cost.trace*31+arg+1;}
#define RAY_PROBE_BUFFERS device const uint* ray_maps [[buffer(6)]],device const uint* ray_refs [[buffer(7)]],
#define RAY_PROBE_ARGUMENTS ,ray_maps,ray_refs
#include "ray_group_work.metal"
