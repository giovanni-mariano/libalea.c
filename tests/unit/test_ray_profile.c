// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0
#include "alea_test.h"
#ifndef ALEA_RAY_PROFILE
#define ALEA_RAY_PROFILE
#endif
#include "util/ray_profile.h"
static double fake_time;
static double profile_clock(void) {return fake_time;}
#define alea_monotonic_seconds profile_clock
#include "../../src/util/ray_profile.c"
#undef alea_monotonic_seconds

RAY_PROFILE_FUNCTION(int, profile_early_return, (int fail), (fail), RAY_RECOVERY) {
    fake_time += 2;
    if (fail) return -1;
    return 0;
}

TEST(ray_profile_wrapper_balances_early_return) {
    fake_time = 0;
    ASSERT_EQ(alea_ray_profile_start(0), 0);
    ASSERT_EQ(profile_early_return(1), -1);
    ASSERT_EQ(alea_ray_profile.stage, RAY_OTHER);
    ASSERT_EQ(alea_ray_profile.seconds[RAY_RECOVERY], 2);
    ASSERT_EQ(profile_early_return(0), 0);
    alea_ray_profile_stop();
    ASSERT_EQ(alea_ray_profile.seconds[RAY_RECOVERY], 4);
    ASSERT_EQ(alea_ray_profile.calls[RAY_RECOVERY], 2);
    ASSERT_EQ(profile_early_return(1), -1);
    ASSERT_EQ(alea_ray_profile.seconds[RAY_RECOVERY], 4);
    alea_ray_profile_free();
}

TEST(ray_profile_nested_exclusive_accounting) {
    fake_time=0;ASSERT_EQ(alea_ray_profile_start(0),0);
    fake_time=1;int root=alea_ray_profile_enter(RAY_HIERARCHY);
    fake_time=3;int csg=alea_ray_profile_enter(RAY_MEMBERSHIP);
    fake_time=4;int nested=alea_ray_profile_enter(RAY_MEMBERSHIP);
    fake_time=5;alea_ray_profile_leave(nested);
    fake_time=6;alea_ray_profile_leave(csg);
    fake_time=8;alea_ray_profile_leave(root);
    fake_time=9;alea_ray_profile_stop();
    ASSERT_EQ(alea_ray_profile.seconds[RAY_OTHER],2);
    ASSERT_EQ(alea_ray_profile.seconds[RAY_HIERARCHY],4);
    ASSERT_EQ(alea_ray_profile.seconds[RAY_MEMBERSHIP],3);
    ASSERT_EQ(alea_ray_profile.calls[RAY_MEMBERSHIP],2);
    ASSERT_EQ(alea_ray_profile.stage,RAY_OTHER);
    fake_time=20;ASSERT_EQ(alea_ray_profile_enter(RAY_RECOVERY),-1);
    alea_ray_profile_leave(-1);alea_ray_profile_stop();
    ASSERT_EQ(alea_ray_profile.seconds[RAY_OTHER],2);
    alea_ray_profile_free();alea_ray_profile_free();
}
TEST(ray_profile_attributes_recursive_descent_to_each_universe_once) {
    fake_time = 0;
    ASSERT_EQ(alea_ray_profile_start(2), 0);
    ray_profile_universe* parent = alea_ray_profile.universes;
    ray_profile_universe* child = parent + 1;
    alea_ray_profile_switch_universe(parent);
    int previous = alea_ray_profile_enter(RAY_MEMBERSHIP);
    fake_time = 1;
    alea_ray_profile_leave(previous);
    fake_time = 2;
    alea_ray_profile_switch_universe(child);
    previous = alea_ray_profile_enter(RAY_MEMBERSHIP);
    fake_time = 5;
    alea_ray_profile_leave(previous);
    fake_time = 6;
    alea_ray_profile_switch_universe(parent);
    fake_time = 8;
    alea_ray_profile_switch_universe(NULL);
    ASSERT_EQ(parent->lookup_seconds, 4);
    ASSERT_EQ(child->lookup_seconds, 4);
    ASSERT_EQ(parent->membership_seconds, 1);
    ASSERT_EQ(child->membership_seconds, 3);
    ASSERT_EQ(parent->membership_calls, 1);
    ASSERT_EQ(child->membership_calls, 1);
    alea_ray_profile_free();
}

TEST(ray_profile_root_aggregation_and_bounded_collisions) {
    ASSERT_EQ(alea_ray_profile_start(0),0);
    alea_ray_profile_root_record(0,.25);alea_ray_profile_root_record(0,.75);
    ASSERT_EQ(alea_ray_profile.roots[0].calls,2);
    ASSERT_EQ(alea_ray_profile.roots[0].seconds,1);
    /* Multiples of the table size deliberately share the initial bucket. */
    for(unsigned i=1;i<=32;++i)alea_ray_profile_root_record(i*RAY_PROFILE_ROOTS,1);
    ASSERT_EQ(alea_ray_profile.dropped_roots,1);
    ASSERT_EQ(alea_ray_profile_start(0),0);
    ASSERT_EQ(alea_ray_profile.dropped_roots,0);
    ASSERT_EQ(alea_ray_profile.roots[0].calls,0);
    alea_ray_profile_free();
}

TEST(ray_profile_list_identity_counts_reset_and_saturation) {
    ASSERT_EQ(alea_ray_profile_start(0), 0);
    /* The address hash repeats every 64 KiB. Use real pointers into one
     * allocation so this also exercises collisions under the sanitizers. */
    const size_t stride = RAY_PROFILE_LISTS * 16 / sizeof(uint32_t);
    uint32_t *cards = calloc(33 * stride, sizeof(*cards));
    ASSERT_NOT_NULL(cards);
    alea_ray_profile_list_record(cards, 10, 8);
    alea_ray_profile_list_record(cards, 10, 7);
    unsigned bucket = ((uintptr_t)cards >> 4) & (RAY_PROFILE_LISTS - 1);
    ASSERT_EQ(alea_ray_profile.lists[bucket].calls, 2);
    ASSERT_EQ(alea_ray_profile.lists[bucket].tests, 15);
    for (unsigned i = 1; i <= 32; ++i)
        alea_ray_profile_list_record(cards + i * stride, 1, 1);
    ASSERT_EQ(alea_ray_profile.dropped_lists, 1);
    alea_ray_profile_stop();
    alea_ray_profile_list_record(cards, 10, 100);
    ASSERT_EQ(alea_ray_profile.lists[bucket].tests, 15);
    ASSERT_EQ(alea_ray_profile_start(0), 0);
    alea_ray_profile_list_record(cards, 10, 5);
    alea_ray_profile_list_record(cards, 11, 6);
    ASSERT_EQ(alea_ray_profile.lists[bucket].tests, 5);
    ASSERT_EQ(alea_ray_profile.lists[(bucket + 1) & (RAY_PROFILE_LISTS - 1)].tests, 6);
    ASSERT_EQ(alea_ray_profile.dropped_lists, 0);
    free(cards);
    alea_ray_profile_free();
}
TEST_MAIN()
