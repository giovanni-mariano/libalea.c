// SPDX-FileCopyrightText: 2026 Giovanni MARIANO
// SPDX-License-Identifier: MPL-2.0

#include "alea_test.h"
#include "alea.h"
#include "alea_mcnp.h"
#include "alea_raycast.h"
#include "raycast/ray_intersect.h"
#include "primitives/primitive_eval.h"
#include "core/alea_system.h"
#include "slice/curve_intersect.h"

static alea_primitive_type_t torus_type(int axis) {
    return (alea_primitive_type_t)(ALEA_PRIMITIVE_TORUS_X + axis);
}

TEST(elliptical_torus_axial_crossings_match_both_semiwidths) {
    const double widths[] = {.25, 2, 1e-6, 1e9};
    for (int axis = 0; axis < 3; axis++) {
        int radial = (axis + 1) % 3;
        for (size_t w = 0; w < 4; w++) {
            double b = widths[w];
            alea_primitive_data_t data = {0};
            data.torus = (alea_torus_data_t){.axis=axis, .major_radius=5,
                .minor_radius=1, .axial_semiwidth_B=b};
            double origin[3] = {0}, direction[3] = {0};
            origin[radial] = 5; origin[axis] = 3*b; direction[axis] = -1;
            alea_ray_t ray;
            ASSERT_EQ(alea_ray_init(&ray, origin[0], origin[1], origin[2],
                                   direction[0], direction[1], direction[2]), 0);
            double roots[4];
            int count = ray_intersect_torus(&ray, &data.torus, roots);
            ASSERT_EQ(count, 2);
            ASSERT_NEAR(roots[0], 2*b, 1e-8*(1+b));
            ASSERT_NEAR(roots[1], 4*b, 1e-8*(1+b));
            for (int j = 0; j < count; j++) {
                double x, y, z;
                alea_ray_point_at(&ray, roots[j], &x, &y, &z);
                ASSERT_NEAR(alea_primitive_eval(torus_type(axis), &data, x,y,z), 0, 1e-8);
                ASSERT(ray_torus_window_may_hit(&ray, &data.torus,
                        nextafter(roots[j], 0), nextafter(roots[j], INFINITY)));
            }
        }
    }
}

TEST(elliptical_torus_radial_roots_tangency_and_miss) {
    alea_torus_data_t torus = {.axis=ALEA_AXIS_Z, .major_radius=5,
                              .minor_radius=1, .axial_semiwidth_B=2};
    alea_ray_t ray;
    double roots[4];
    ASSERT_EQ(alea_ray_init(&ray,-10,0,1,1,0,0), 0);
    ASSERT_EQ(ray_intersect_torus(&ray,&torus,roots), 4);
    double half = sqrt(.75);
    const double expected[] = {5-half,5+half,15-half,15+half};
    for (int j=0;j<4;j++) ASSERT_NEAR(roots[j],expected[j],1e-8);
    ASSERT_EQ(alea_ray_init(&ray,-10,0,2,1,0,0), 0);
    ASSERT_EQ(ray_intersect_torus(&ray,&torus,roots), 4);
    ASSERT_NEAR(roots[0],5,1e-7); ASSERT_NEAR(roots[1],5,1e-7);
    ASSERT_NEAR(roots[2],15,1e-7); ASSERT_NEAR(roots[3],15,1e-7);
    ASSERT_EQ(alea_ray_init(&ray,-10,0,2.01,1,0,0), 0);
    ASSERT_EQ(ray_intersect_torus(&ray,&torus,roots), 0);
}

TEST(elliptical_torus_normals_follow_the_evaluator_gradient) {
    for (int axis=0;axis<3;axis++) for (int inward=0;inward<2;inward++) {
        int radial=(axis+1)%3;
        alea_primitive_data_t data={0};
        data.torus=(alea_torus_data_t){.axis=axis,.major_radius=5,
            .minor_radius=1,.axial_semiwidth_B=2};
        double p[3]={0}, expected[3]={0};
        double sign=inward?-1:1;
        p[radial]=5+sign*.6; p[axis]=1.6;
        expected[radial]=sign*.6; expected[axis]=.4;
        double length=hypot(.6,.4), normal[3];
        primitive_normal_at(torus_type(axis),&data,p[0],p[1],p[2],
                             &normal[0],&normal[1],&normal[2]);
        for (int k=0;k<3;k++) {
            ASSERT_NEAR(normal[k],expected[k]/length,1e-10);
            double lo[3]={p[0],p[1],p[2]}, hi[3]={p[0],p[1],p[2]};
            lo[k]-=1e-5; hi[k]+=1e-5;
            double derivative=(alea_primitive_eval(torus_type(axis),&data,hi[0],hi[1],hi[2])-
                alea_primitive_eval(torus_type(axis),&data,lo[0],lo[1],lo[2]))/2e-5;
            ASSERT_NEAR(normal[k],derivative/(2*length),1e-8);
        }
    }
}

TEST(mcnp_elliptical_torus_tracks_have_analytic_material_lengths) {
    mcnp_model_t* model=mcnp_load("tests/data/mcnp_elliptical_torus.mcnp");
    ASSERT_NOT_NULL(model);
    const double origins[]={5,0,3, 5,0,-3};
    const double directions[]={0,0,-1, 0,0,1};
    alea_raycast_batch_result_t* tracks=alea_raycast_batch_result_create();
    ASSERT_NOT_NULL(tracks);
    for (int rich=0;rich<2;rich++) {
        alea_raycast_batch_options_t options={.struct_size=sizeof(options),
            .fields=ALEA_RAY_BATCH_MATERIAL | (rich?ALEA_RAY_BATCH_SURFACES:0),
            .projected_depth=-1};
        ASSERT_EQ(alea_raycast_hier_batch(model->sys,origins,directions,2,6,&options,tracks),0);
        ASSERT_EQ(alea_raycast_batch_segment_count(tracks),6);
        const double* enter=alea_raycast_batch_t_enter(tracks);
        const double* exit=alea_raycast_batch_t_exit(tracks);
        const int32_t* cells=alea_raycast_batch_cell_ids(tracks);
        const int32_t* materials=alea_raycast_batch_material_ids(tracks);
        const double expected_enter[]={0,1,5}, expected_exit[]={1,5,6};
        const int expected_cells[]={2,1,2};
        for (int i=0;i<6;i++) {
            ASSERT_NEAR(enter[i],expected_enter[i%3],1e-8);
            ASSERT_NEAR(exit[i],expected_exit[i%3],1e-8);
            ASSERT_EQ(cells[i],expected_cells[i%3]);
            ASSERT_EQ(materials[i],expected_cells[i%3]);
        }
    }
    alea_raycast_batch_result_destroy(tracks);
    mcnp_model_destroy(model);
}

TEST(elliptical_torus_oblique_roots_and_normals_match_point_classification) {
    uint64_t seed=793;
    size_t checked=0;
    for (int axis=0;axis<3;axis++) for (int i=0;i<400;i++) {
        alea_primitive_data_t data={0};
        data.torus=(alea_torus_data_t){.axis=axis,.center_x=1e6,
            .center_y=-2e6,.center_z=3e6,.major_radius=5,
            .minor_radius=1,.axial_semiwidth_B=(i&1)?.5:2};
        double values[6];
        for (int k=0;k<6;k++) {
            seed=seed*UINT64_C(6364136223846793005)+1;
            values[k]=2.0*(double)(seed>>11)/9007199254740992.0-1;
        }
        alea_ray_t ray;
        ASSERT_EQ(alea_ray_init(&ray,1e6+8*values[0],-2e6+8*values[1],
            3e6+8*values[2],values[3],values[4],values[5]),0);
        double roots[4];
        int count=ray_intersect_torus(&ray,&data.torus,roots);
        for (int j=0;j<count;j++) {
            double x,y,z;
            alea_ray_point_at(&ray,roots[j],&x,&y,&z);
            ASSERT_NEAR(alea_primitive_eval(torus_type(axis),&data,x,y,z),0,1e-7);
            ASSERT(ray_torus_window_may_hit(&ray,&data.torus,
                nextafter(roots[j],0),nextafter(roots[j],INFINITY)));
            checked++;
        }
    }
    ASSERT(checked>100);
}

TEST(elliptical_torus_zero_major_radius_spindle_and_unbounded_axial_cases) {
    alea_torus_data_t torus={.axis=ALEA_AXIS_Z,.major_radius=0,
        .minor_radius=1,.axial_semiwidth_B=2};
    alea_ray_t ray;
    double roots[4];
    ASSERT_EQ(alea_ray_init(&ray,0,0,3,0,0,-1),0);
    ASSERT_EQ(ray_intersect_torus(&ray,&torus,roots),2);
    ASSERT_NEAR(roots[0],1,1e-10); ASSERT_NEAR(roots[1],5,1e-10);
    torus.major_radius=.5;
    ASSERT_EQ(alea_ray_init(&ray,-3,0,0,1,0,0),0);
    ASSERT_EQ(ray_intersect_torus(&ray,&torus,roots),2);
    ASSERT_NEAR(roots[0],1.5,1e-10); ASSERT_NEAR(roots[1],4.5,1e-10);
    torus.major_radius=5; torus.axial_semiwidth_B=0;
    ASSERT_EQ(alea_ray_init(&ray,-10,0,100,1,0,0),0);
    ASSERT_EQ(ray_intersect_torus(&ray,&torus,roots),4);
    const double expected[]={4,6,14,16};
    for (int j=0;j<4;j++) ASSERT_NEAR(roots[j],expected[j],1e-8);
    ASSERT_EQ(alea_ray_init(&ray,5,0,100,0,0,1),0);
    ASSERT_EQ(ray_intersect_torus(&ray,&torus,roots),0);
}

TEST(transformed_elliptical_torus_navigation_preserves_physical_lengths) {
    mcnp_model_t* model=mcnp_load("tests/data/mcnp_elliptical_torus.mcnp");
    ASSERT_NOT_NULL(model);
    alea_system_t* sys=model->sys;
    ASSERT_EQ(alea_cell_set_universe(sys,0,1),0);
    ASSERT_EQ(alea_cell_set_universe(sys,1,1),0);
    int box=alea_box_surface(sys,2,-100,100,-100,100,-100,100);
    ASSERT(box>=0);
    int parent=alea_add_cell(sys,10,alea_halfspace(sys,box,-1),ALEA_MATERIAL_VOID,0,0);
    /* Rotate local Z onto world Y and translate the occurrence. */
    const double transform[]={10,20,30, 1,0,0, 0,0,-1, 0,1,0};
    ASSERT_EQ(alea_add_transform(sys,1,transform,12,0),0);
    ASSERT_EQ(alea_set_fill(sys,parent,1,1),0);
    alea_ray_navigator_t* nav=alea_ray_navigator_create(sys);
    ASSERT_NOT_NULL(nav);
    ASSERT_EQ(alea_ray_navigator_set_validation_mode(nav,ALEA_NAV_VALIDATE_STRICT),0);
    const double origin[]={15,23,30}, direction[]={0,-1,0};
    alea_nav_location_t location;
    ASSERT_EQ(alea_ray_navigator_restart(nav,origin,direction,&location),0);
    ASSERT_EQ(location.cell_id,2);
    alea_nav_event_t event;
    ASSERT_EQ(alea_ray_navigator_advance(nav,INFINITY,6,&event),0);
    ASSERT_EQ(event.kind,ALEA_NAV_BOUNDARY);
    ASSERT_NEAR(event.distance,1,1e-8);
    ASSERT_EQ(event.after.cell_id,1);
    ASSERT_EQ(alea_ray_navigator_advance(nav,INFINITY,5,&event),0);
    ASSERT_EQ(event.kind,ALEA_NAV_BOUNDARY);
    ASSERT_NEAR(event.distance,4,1e-8);
    ASSERT_EQ(event.after.cell_id,2);
    alea_ray_navigator_destroy(nav);
    mcnp_model_destroy(model);
}

TEST(elliptical_torus_slice_curves_match_analytic_sections) {
    for (int axis=0;axis<3;axis++) {
        alea_primitive_data_t data={0};
        data.torus=(alea_torus_data_t){.axis=axis,.major_radius=5,
            .minor_radius=1,.axial_semiwidth_B=2};
        int other=(axis+1)%3;
        double normal[3]={0}, up[3]={0}; normal[other]=1; up[axis]=1;
        alea_slice_plane_t plane;
        alea_slice_plane_init(&plane,0,0,0,normal[0],normal[1],normal[2],up[0],up[1],up[2]);
        alea_curve_2d_t curve;
        ASSERT(alea_intersect_primitive_plane(torus_type(axis),&data,&plane,&curve));
        double roots[4];
        ASSERT_EQ(alea_curve_scanline_intersect(&curve,1,roots,4),4);
        double half=sqrt(.75), expected[]={-5-half,-5+half,5-half,5+half};
        double u_min,u_max,v_min,v_max;
        alea_curve_bbox(&curve,&u_min,&u_max,&v_min,&v_max);
        ASSERT(u_min<=-6 && u_max>=6 && v_min<=-2 && v_max>=2);
        for (int j=0;j<4;j++) {
            ASSERT_NEAR(roots[j],expected[j],1e-8);
            double x,y,z;
            alea_plane_to_3d(&plane,roots[j],1,&x,&y,&z);
            ASSERT_NEAR(alea_primitive_eval(torus_type(axis),&data,x,y,z),0,1e-8);
        }
        /* A perpendicular plane at axial offset 1.5 still cuts the ellipse. */
        alea_slice_plane_init_axis(&plane,axis,1.5);
        ASSERT(alea_intersect_primitive_plane(torus_type(axis),&data,&plane,&curve));
        ASSERT_EQ(alea_curve_scanline_intersect(&curve,0,roots,4),4);
        for (int j=0;j<4;j++) for (int k=j+1;k<4;k++)
            if (roots[k]<roots[j]) { double tmp=roots[j]; roots[j]=roots[k]; roots[k]=tmp; }
        half=sqrt(1-1.5*1.5/4);
        const double perpendicular[]={-5-half,-5+half,5-half,5+half};
        for (int j=0;j<4;j++) ASSERT_NEAR(roots[j],perpendicular[j],1e-8);
        alea_slice_plane_init_axis(&plane,axis,2.1);
        ASSERT(!alea_intersect_primitive_plane(torus_type(axis),&data,&plane,&curve));
    }
}

TEST(elliptical_torus_normals_remain_unit_for_very_wide_tubes) {
    alea_primitive_data_t data={0};
    data.torus=(alea_torus_data_t){.axis=ALEA_AXIS_Z,.major_radius=5,
        .minor_radius=1,.axial_semiwidth_B=1e15};
    for (int sign=-1;sign<=1;sign+=2) {
        double nx,ny,nz;
        primitive_normal_at(ALEA_PRIMITIVE_TORUS_Z,&data,5,0,sign*1e15,&nx,&ny,&nz);
        ASSERT_NEAR(nx,0,1e-12); ASSERT_NEAR(ny,0,1e-12); ASSERT_NEAR(nz,sign,1e-12);
    }
}

TEST(elliptical_torus_thin_ring_keeps_axial_crossings) {
    alea_torus_data_t torus={.axis=ALEA_AXIS_Z,.major_radius=1e6,
        .minor_radius=1,.axial_semiwidth_B=2};
    alea_ray_t ray;
    ASSERT_EQ(alea_ray_init(&ray,1e6,0,6,0,0,-1),0);
    double roots[4];
    ASSERT_EQ(ray_intersect_torus(&ray,&torus,roots),2);
    ASSERT_NEAR(roots[0],4,1e-8); ASSERT_NEAR(roots[1],8,1e-8);
}

TEST(elliptical_slice_scanlines_keep_original_coordinates_after_scaling) {
    alea_primitive_data_t data={0};
    data.torus=(alea_torus_data_t){.axis=ALEA_AXIS_Z,.major_radius=5,
        .minor_radius=1,.axial_semiwidth_B=1e6};
    alea_slice_plane_t plane;
    /* U is axial, V radial, forcing transformed direction normalization. */
    alea_slice_plane_init(&plane,0,0,0,0,1,0,1,0,0);
    alea_curve_2d_t curve;
    ASSERT(alea_intersect_primitive_plane(ALEA_PRIMITIVE_TORUS_Z,&data,&plane,&curve));
    double roots[4];
    ASSERT_EQ(alea_curve_scanline_intersect(&curve,5,roots,4),2);
    ASSERT_NEAR(roots[0],-1e6,1e-3); ASSERT_NEAR(roots[1],1e6,1e-3);
}

TEST(elliptical_torus_raster_matches_independent_pixel_formula) {
    mcnp_model_t* model=mcnp_load("tests/data/mcnp_elliptical_torus.mcnp");
    ASSERT_NOT_NULL(model);
    enum { SIZE=32, COUNT=SIZE*SIZE };
    int32_t cells[COUNT], materials[COUNT];
    double densities[COUNT]; uint8_t flags[COUNT];
    alea_slice_view_t view;
    alea_slice_view_axis(&view,2,1.5,-8,8,-8,8);
    alea_slice_raster_t raster;
    alea_slice_raster_init(&raster);
    raster.nu=raster.nv=SIZE;
    raster.fields=ALEA_SLICE_RASTER_CELL_ID|ALEA_SLICE_RASTER_MATERIAL_ID|
                  ALEA_SLICE_RASTER_DENSITY|ALEA_SLICE_RASTER_RESOLUTION_FLAGS;
    raster.cell_ids=cells; raster.material_ids=materials;
    raster.densities=densities; raster.resolution_flags=flags;
    alea_raycast_batch_result_t* tracks=alea_raycast_batch_result_create();
    ASSERT_NOT_NULL(tracks);
    for (int fused=0;fused<2;fused++) {
        if (fused) {
            alea_slice_raster_options_t options;
            alea_slice_raster_options_init(&options);
            ASSERT_EQ(alea_trace_ray_slice_raster(model->sys,&view,&options,&raster),0);
        } else {
            alea_raycast_batch_options_t options={.struct_size=sizeof(options),
                .fields=ALEA_RAY_BATCH_MATERIAL|ALEA_RAY_BATCH_DENSITY|
                        ALEA_RAY_BATCH_RESOLUTION_FLAGS,.projected_depth=-1};
            ASSERT_EQ(alea_trace_ray_slice_compact(model->sys,&view,SIZE,&options,tracks),0);
            ASSERT_EQ(alea_rasterize_ray_slice_compact(&view,tracks,&raster),0);
        }
        for (int row=0;row<SIZE;row++) for (int col=0;col<SIZE;col++) {
            double x=-8+(col+.5)*16/SIZE, y=-8+(row+.5)*16/SIZE;
            double radial=hypot(x,y)-5;
            int expected=radial*radial+1.5*1.5/4<1 ? 1 : 2;
            int i=row*SIZE+col;
            ASSERT_EQ(cells[i],expected); ASSERT_EQ(materials[i],expected);
            ASSERT_EQ(densities[i],expected==1?1:2); ASSERT_EQ(flags[i],0);
        }
    }
    alea_raycast_batch_result_destroy(tracks);
    mcnp_model_destroy(model);
}

TEST_MAIN()
