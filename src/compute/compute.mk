# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0

USE_METAL ?= 0
ifeq ($(USE_METAL),1)
  ifneq ($(UNAME_S),Darwin)
    $(error USE_METAL=1 requires macOS)
  endif
  COMPUTE_VARIANT = metal
  COMPUTE_CPPFLAGS = -DALEA_USE_METAL=1
  COMPUTE_LDLIBS = -framework Foundation -framework Metal
else ifeq ($(USE_METAL),0)
  COMPUTE_VARIANT = cpu
else
  $(error USE_METAL must be 0 or 1)
endif

COMPUTE_BUILD_DIR = $(BUILD_DIR)/compute/$(COMPUTE_VARIANT)
COMPUTE_BIN_DIR = $(BIN_DIR)/compute/$(COMPUTE_VARIANT)
LIB_COMPUTE = $(COMPUTE_BIN_DIR)/libalea_compute.a
COMPUTE_OBJS = $(COMPUTE_BUILD_DIR)/compute.o
ifeq ($(USE_METAL),1)
  COMPUTE_OBJS += $(COMPUTE_BUILD_DIR)/compute_scene.o $(COMPUTE_BUILD_DIR)/compute_groups.o $(COMPUTE_BUILD_DIR)/compute_render_scene.o $(COMPUTE_BUILD_DIR)/compute_bvh.o $(COMPUTE_BUILD_DIR)/compute_spatial.o $(COMPUTE_BUILD_DIR)/compute_local.o $(COMPUTE_BUILD_DIR)/compute_point.o $(COMPUTE_BUILD_DIR)/compute_profile.o $(COMPUTE_BUILD_DIR)/compute_metal.o
endif

.PHONY: compute test-compute benchmark-compute install-compute
compute: lib-core $(LIB_COMPUTE)

$(COMPUTE_BUILD_DIR) $(COMPUTE_BIN_DIR):
	@$(MKDIR_P) $@

$(COMPUTE_BUILD_DIR)/%.o: $(SRC_DIR)/compute/%.c | $(COMPUTE_BUILD_DIR)
	@echo "CC  $< ($(COMPUTE_VARIANT))"
	@$(CC) $(CFLAGS) $(DEPFLAGS) $(COMPUTE_CPPFLAGS) $(INCLUDES) -c $< -o $@

$(COMPUTE_BUILD_DIR)/compute.air: $(SRC_DIR)/compute/compute.metal $(SRC_DIR)/compute/classify.metal $(SRC_DIR)/compute/select.metal $(SRC_DIR)/compute/render.metal $(SRC_DIR)/compute/resume.metal $(SRC_DIR)/compute/intersect.metal $(SRC_DIR)/compute/compute_layout.h $(SRC_DIR)/compute/compute_root_reuse.h $(SRC_DIR)/compute/groups.metal $(SRC_DIR)/compute/ray_groups.metal | $(COMPUTE_BUILD_DIR)
	@echo "METAL $<"
	@xcrun --sdk macosx metal -std=macos-metal2.3 -fno-fast-math \
		-fmodules-cache-path=$(abspath $(COMPUTE_BUILD_DIR))/module-cache -c $< -o $@

$(COMPUTE_BUILD_DIR)/compute.metallib: $(COMPUTE_BUILD_DIR)/compute.air
	@xcrun --sdk macosx metallib $< -o $@

$(COMPUTE_BUILD_DIR)/compute_metallib.h: $(COMPUTE_BUILD_DIR)/compute.metallib
	@xxd -i -n alea_compute_metallib $< > $@

$(COMPUTE_BUILD_DIR)/compute_metal.o: $(SRC_DIR)/compute/compute_metal.m $(SRC_DIR)/compute/compute_groups_metal.h $(COMPUTE_BUILD_DIR)/compute_metallib.h
	@echo "OBJC $<"
	@xcrun --sdk macosx clang $(CFLAGS) $(DEPFLAGS) $(COMPUTE_CPPFLAGS) \
		$(INCLUDES) -I$(COMPUTE_BUILD_DIR) -fno-objc-arc -c $< -o $@

$(LIB_COMPUTE): $(COMPUTE_OBJS) $(SRC_DIR)/compute/compute.mk | $(COMPUTE_BIN_DIR)
	@$(RM) $@.tmp
	@$(AR) rcs $@.tmp $(COMPUTE_OBJS)
	@mv $@.tmp $@

$(COMPUTE_BIN_DIR)/test_compute$(EXEEXT): $(TEST_DIR)/compute/test_compute.c $(LIB_COMPUTE) $(LIB_CORE)
	@$(CC) $(CFLAGS) $(COMPUTE_CPPFLAGS) $(INCLUDES) -I$(TEST_DIR) $< \
		$(LIB_COMPUTE) $(LIB_CORE) $(LDFLAGS) $(COMPUTE_LDLIBS) -o $@

test-compute: check-public-headers compute $(COMPUTE_BIN_DIR)/test_compute$(EXEEXT)
	@$(COMPUTE_BIN_DIR)/test_compute$(EXEEXT)
ifeq ($(USE_METAL),1)
	@$(COMPUTE_BIN_DIR)/test_compute_builders
	@$(COMPUTE_BIN_DIR)/test_metal_resources

test-compute: $(COMPUTE_BIN_DIR)/test_compute_builders $(COMPUTE_BIN_DIR)/test_metal_resources

$(COMPUTE_BIN_DIR)/test_compute_builders: $(TEST_DIR)/compute/test_compute_builders.c benchmarks/region_program_audit.h benchmarks/region_group_audit.h $(SRC_DIR)/compute/compute_bvh.c $(SRC_DIR)/compute/compute_spatial.c $(SRC_DIR)/compute/compute_point.c $(LIB_COMPUTE) $(LIB_CORE)
	@$(CC) $(CFLAGS) $(COMPUTE_CPPFLAGS) $(INCLUDES) -I$(TEST_DIR) $< \
		$(LIB_COMPUTE) $(LIB_CORE) $(LDFLAGS) $(COMPUTE_LDLIBS) -o $@

$(COMPUTE_BUILD_DIR)/group_eval_probe.air: $(TEST_DIR)/compute/group_eval_probe.metal $(COMPUTE_BUILD_DIR)/compute.air
	@xcrun --sdk macosx metal -std=macos-metal2.3 -fno-fast-math -fmodules-cache-path=$(abspath $(COMPUTE_BUILD_DIR))/module-cache -c $< -o $@

$(COMPUTE_BUILD_DIR)/group_eval_probe.metallib: $(COMPUTE_BUILD_DIR)/group_eval_probe.air
	@xcrun --sdk macosx metallib $< -o $@

$(COMPUTE_BUILD_DIR)/group_eval_probe_metallib.h: $(COMPUTE_BUILD_DIR)/group_eval_probe.metallib
	@xxd -i -n alea_group_eval_probe_metallib $< > $@

$(COMPUTE_BIN_DIR)/test_metal_resources: $(TEST_DIR)/compute/test_metal_resources.m $(TEST_DIR)/compute/test_group_trace.h $(COMPUTE_BUILD_DIR)/group_eval_probe_metallib.h $(COMPUTE_BUILD_DIR)/ray_eval_probe_metallib.h $(TEST_DIR)/compute/test_ray_contract.h $(SRC_DIR)/compute/compute_metal.m $(SRC_DIR)/compute/compute_groups_metal.h $(COMPUTE_BUILD_DIR)/compute_metallib.h $(LIB_COMPUTE) $(LIB_CORE)
	@xcrun --sdk macosx clang $(CFLAGS) $(COMPUTE_CPPFLAGS) $(INCLUDES) -I$(TEST_DIR) -I$(COMPUTE_BUILD_DIR) -fno-objc-arc $< \
		$(LIB_COMPUTE) $(LIB_CORE) $(LDFLAGS) $(COMPUTE_LDLIBS) -o $@
endif

$(COMPUTE_BIN_DIR)/benchmark_compute$(EXEEXT): benchmarks/compute.c $(LIB_COMPUTE) $(LIB_CORE)
	@$(CC) $(CFLAGS) $(COMPUTE_CPPFLAGS) $(INCLUDES) $< \
		$(LIB_COMPUTE) $(LIB_CORE) $(LDFLAGS) $(COMPUTE_LDLIBS) -o $@

benchmark-compute: compute $(COMPUTE_BIN_DIR)/benchmark_compute$(EXEEXT)
	@$(COMPUTE_BIN_DIR)/benchmark_compute$(EXEEXT)

install-compute: compute
	@$(MKDIR_P) "$(DESTDIR)$(LIBDIR)" "$(DESTDIR)$(INCLUDEDIR)"
	@$(INSTALL_DATA) $(LIB_COMPUTE) "$(DESTDIR)$(LIBDIR)/"
	@$(INSTALL_DATA) $(INCLUDE_DIR)/alea_compute.h "$(DESTDIR)$(INCLUDEDIR)/"

-include $(COMPUTE_OBJS:.o=.d)

ifeq ($(USE_METAL),1)
$(COMPUTE_BIN_DIR)/benchmark_gpu_geometry$(EXEEXT): benchmarks/gpu_geometry.c benchmarks/gpu_geometry_fixture.h $(LIB_COMPUTE) $(LIB_CORE)
	@$(CC) $(CFLAGS) $(COMPUTE_CPPFLAGS) $(INCLUDES) $< \
		$(LIB_COMPUTE) $(LIB_CORE) $(LDFLAGS) $(COMPUTE_LDLIBS) -o $@

$(COMPUTE_BUILD_DIR)/point_work_audit.air: benchmarks/point_work_audit.metal benchmarks/point_work_audit_layout.h $(COMPUTE_BUILD_DIR)/compute.air
	@xcrun --sdk macosx metal -std=macos-metal2.3 -fno-fast-math \
		-fmodules-cache-path=$(abspath $(COMPUTE_BUILD_DIR))/module-cache -c $< -o $@

$(COMPUTE_BUILD_DIR)/point_work_audit.metallib: $(COMPUTE_BUILD_DIR)/point_work_audit.air
	@xcrun --sdk macosx metallib $< -o $@

$(COMPUTE_BUILD_DIR)/point_work_audit_metallib.h: $(COMPUTE_BUILD_DIR)/point_work_audit.metallib
	@xxd -i -n alea_point_work_audit_metallib $< > $@

$(COMPUTE_BIN_DIR)/benchmark_metal_point_work$(EXEEXT): benchmarks/metal_point_work.m benchmarks/point_work_audit_layout.h benchmarks/gpu_geometry_fixture.h $(SRC_DIR)/compute/compute_metal.m $(SRC_DIR)/compute/compute_groups_metal.h $(COMPUTE_BUILD_DIR)/compute_metallib.h $(COMPUTE_BUILD_DIR)/point_work_audit_metallib.h $(LIB_COMPUTE) $(LIB_CORE)
	@xcrun --sdk macosx clang $(CFLAGS) $(COMPUTE_CPPFLAGS) $(INCLUDES) -I$(COMPUTE_BUILD_DIR) -fno-objc-arc $< \
		$(LIB_COMPUTE) $(LIB_CORE) $(LDFLAGS) $(COMPUTE_LDLIBS) -o $@

$(COMPUTE_BUILD_DIR)/region_work_audit.air: benchmarks/region_work_audit.metal benchmarks/point_work_audit_layout.h $(COMPUTE_BUILD_DIR)/compute.air
	@xcrun --sdk macosx metal -std=macos-metal2.3 -fno-fast-math -fmodules-cache-path=$(abspath $(COMPUTE_BUILD_DIR))/module-cache -c $< -o $@

$(COMPUTE_BUILD_DIR)/region_work_audit.metallib: $(COMPUTE_BUILD_DIR)/region_work_audit.air
	@xcrun --sdk macosx metallib $< -o $@

$(COMPUTE_BUILD_DIR)/region_work_audit_metallib.h: $(COMPUTE_BUILD_DIR)/region_work_audit.metallib
	@xxd -i -n alea_region_work_audit_metallib $< > $@

$(COMPUTE_BIN_DIR)/benchmark_metal_region_work$(EXEEXT): benchmarks/metal_region_work.m benchmarks/region_program_audit.h benchmarks/point_work_audit_layout.h benchmarks/gpu_geometry_fixture.h $(SRC_DIR)/compute/compute_metal.m $(SRC_DIR)/compute/compute_groups_metal.h $(COMPUTE_BUILD_DIR)/compute_metallib.h $(COMPUTE_BUILD_DIR)/region_work_audit_metallib.h $(LIB_COMPUTE) $(LIB_CORE)
	@xcrun --sdk macosx clang $(CFLAGS) $(COMPUTE_CPPFLAGS) $(INCLUDES) -I$(COMPUTE_BUILD_DIR) -fno-objc-arc $< $(LIB_COMPUTE) $(LIB_CORE) $(LDFLAGS) $(COMPUTE_LDLIBS) -o $@

$(COMPUTE_BUILD_DIR)/group_work_audit.air: benchmarks/group_work_audit.metal benchmarks/region_work_audit.metal benchmarks/point_work_audit_layout.h $(COMPUTE_BUILD_DIR)/compute.air
	@xcrun --sdk macosx metal -std=macos-metal2.3 -fno-fast-math -fmodules-cache-path=$(abspath $(COMPUTE_BUILD_DIR))/module-cache -c $< -o $@

$(COMPUTE_BUILD_DIR)/group_work_audit.metallib: $(COMPUTE_BUILD_DIR)/group_work_audit.air
	@xcrun --sdk macosx metallib $< -o $@

$(COMPUTE_BUILD_DIR)/group_work_audit_metallib.h: $(COMPUTE_BUILD_DIR)/group_work_audit.metallib
	@xxd -i -n alea_group_work_audit_metallib $< > $@

$(COMPUTE_BIN_DIR)/benchmark_metal_group_work$(EXEEXT): benchmarks/metal_group_work.m benchmarks/metal_region_work.m benchmarks/region_group_audit.h benchmarks/region_group_storage.h benchmarks/region_program_audit.h benchmarks/point_work_audit_layout.h benchmarks/gpu_geometry_fixture.h $(SRC_DIR)/compute/compute_metal.m $(SRC_DIR)/compute/compute_groups_metal.h $(COMPUTE_BUILD_DIR)/compute_metallib.h $(COMPUTE_BUILD_DIR)/group_work_audit_metallib.h $(LIB_COMPUTE) $(LIB_CORE)
	@xcrun --sdk macosx clang $(CFLAGS) $(COMPUTE_CPPFLAGS) $(INCLUDES) -I$(COMPUTE_BUILD_DIR) -fno-objc-arc $< $(LIB_COMPUTE) $(LIB_CORE) $(LDFLAGS) $(COMPUTE_LDLIBS) -o $@

$(COMPUTE_BUILD_DIR)/group_execute.air: benchmarks/group_execute.metal $(COMPUTE_BUILD_DIR)/compute.air
	@xcrun --sdk macosx metal -std=macos-metal2.3 -fno-fast-math -fmodules-cache-path=$(abspath $(COMPUTE_BUILD_DIR))/module-cache -c $< -o $@

$(COMPUTE_BUILD_DIR)/group_execute.metallib: $(COMPUTE_BUILD_DIR)/group_execute.air
	@xcrun --sdk macosx metallib $< -o $@

$(COMPUTE_BUILD_DIR)/group_execute_metallib.h: $(COMPUTE_BUILD_DIR)/group_execute.metallib
	@xxd -i -n alea_group_execute_metallib $< > $@

$(COMPUTE_BIN_DIR)/benchmark_metal_group_execute$(EXEEXT): benchmarks/metal_group_execute.m benchmarks/point_work_audit_layout.h benchmarks/metal_region_work.m benchmarks/region_group_storage.h benchmarks/region_group_audit.h benchmarks/region_program_audit.h benchmarks/gpu_geometry_fixture.h $(SRC_DIR)/compute/compute_metal.m $(SRC_DIR)/compute/compute_groups_metal.h $(COMPUTE_BUILD_DIR)/compute_metallib.h $(COMPUTE_BUILD_DIR)/group_execute_metallib.h $(LIB_COMPUTE) $(LIB_CORE)
	@xcrun --sdk macosx clang $(CFLAGS) $(COMPUTE_CPPFLAGS) $(INCLUDES) -I$(COMPUTE_BUILD_DIR) -fno-objc-arc $< $(LIB_COMPUTE) $(LIB_CORE) $(LDFLAGS) $(COMPUTE_LDLIBS) -o $@

$(COMPUTE_BUILD_DIR)/bounds_preview.air: benchmarks/bounds_preview.metal benchmarks/bounds_preview_layout.h | $(COMPUTE_BUILD_DIR)
	@xcrun --sdk macosx metal -std=macos-metal2.3 -fno-fast-math \
		-fmodules-cache-path=$(abspath $(COMPUTE_BUILD_DIR))/module-cache -c $< -o $@

$(COMPUTE_BUILD_DIR)/bounds_preview.metallib: $(COMPUTE_BUILD_DIR)/bounds_preview.air
	@xcrun --sdk macosx metallib $< -o $@

$(COMPUTE_BUILD_DIR)/bounds_preview_metallib.h: $(COMPUTE_BUILD_DIR)/bounds_preview.metallib
	@xxd -i -n alea_bounds_preview_metallib $< > $@

$(COMPUTE_BIN_DIR)/benchmark_metal_bounds_preview$(EXEEXT): benchmarks/metal_bounds_preview.m benchmarks/bounds_preview_layout.h benchmarks/gpu_geometry_fixture.h $(SRC_DIR)/compute/compute_metal.m $(SRC_DIR)/compute/compute_groups_metal.h $(COMPUTE_BUILD_DIR)/compute_metallib.h $(COMPUTE_BUILD_DIR)/bounds_preview_metallib.h $(LIB_COMPUTE) $(LIB_CORE)
	@xcrun --sdk macosx clang $(CFLAGS) $(COMPUTE_CPPFLAGS) $(INCLUDES) -I$(COMPUTE_BUILD_DIR) -fno-objc-arc $< \
		$(LIB_COMPUTE) $(LIB_CORE) $(LDFLAGS) $(COMPUTE_LDLIBS) -o $@

$(COMPUTE_BUILD_DIR)/reprojection.air: benchmarks/reprojection.metal benchmarks/reprojection_layout.h | $(COMPUTE_BUILD_DIR)
	@xcrun --sdk macosx metal -std=macos-metal2.3 -fno-fast-math \
		-fmodules-cache-path=$(abspath $(COMPUTE_BUILD_DIR))/module-cache -c $< -o $@

$(COMPUTE_BUILD_DIR)/reprojection.metallib: $(COMPUTE_BUILD_DIR)/reprojection.air
	@xcrun --sdk macosx metallib $< -o $@

$(COMPUTE_BUILD_DIR)/reprojection_metallib.h: $(COMPUTE_BUILD_DIR)/reprojection.metallib
	@xxd -i -n alea_reprojection_metallib $< > $@

$(COMPUTE_BIN_DIR)/benchmark_metal_reprojection$(EXEEXT): benchmarks/metal_reprojection.m benchmarks/reprojection_layout.h benchmarks/gpu_geometry_fixture.h $(COMPUTE_BUILD_DIR)/reprojection_metallib.h $(LIB_COMPUTE) $(LIB_CORE)
	@xcrun --sdk macosx clang $(CFLAGS) $(COMPUTE_CPPFLAGS) $(INCLUDES) -I$(COMPUTE_BUILD_DIR) -fno-objc-arc $< \
		$(LIB_COMPUTE) $(LIB_CORE) $(LDFLAGS) $(COMPUTE_LDLIBS) -o $@

$(COMPUTE_BIN_DIR)/benchmark_metal_query_scaling$(EXEEXT): benchmarks/metal_query_scaling.m benchmarks/gpu_geometry_fixture.h $(SRC_DIR)/compute/compute_metal.m $(SRC_DIR)/compute/compute_groups_metal.h $(COMPUTE_BUILD_DIR)/compute_metallib.h $(LIB_COMPUTE) $(LIB_CORE)
	@xcrun --sdk macosx clang $(CFLAGS) $(COMPUTE_CPPFLAGS) $(INCLUDES) -I$(COMPUTE_BUILD_DIR) -fno-objc-arc $< \
		$(LIB_COMPUTE) $(LIB_CORE) $(LDFLAGS) $(COMPUTE_LDLIBS) -o $@

$(COMPUTE_BIN_DIR)/benchmark_many_root_points$(EXEEXT): benchmarks/many_root_points.c $(LIB_COMPUTE) $(LIB_CORE)
	@$(CC) $(CFLAGS) $(COMPUTE_CPPFLAGS) $(INCLUDES) $< \
		$(LIB_COMPUTE) $(LIB_CORE) $(LDFLAGS) $(COMPUTE_LDLIBS) -o $@

.PHONY: benchmark-many-root-points
benchmark-many-root-points: compute $(COMPUTE_BIN_DIR)/benchmark_many_root_points$(EXEEXT)
	@$(COMPUTE_BIN_DIR)/benchmark_many_root_points$(EXEEXT)
	@$(COMPUTE_BIN_DIR)/benchmark_many_root_points$(EXEEXT) overlap

.PHONY: benchmark-gpu-geometry
benchmark-gpu-geometry: compute $(COMPUTE_BIN_DIR)/benchmark_gpu_geometry$(EXEEXT)
	@$(COMPUTE_BIN_DIR)/benchmark_gpu_geometry$(EXEEXT)

$(COMPUTE_BIN_DIR)/benchmark_render$(EXEEXT): benchmarks/render_compute.c $(LIB_COMPUTE) $(LIB_CORE)
	@$(CC) $(CFLAGS) $(COMPUTE_CPPFLAGS) $(INCLUDES) $< \
		$(LIB_COMPUTE) $(LIB_CORE) $(LDFLAGS) $(COMPUTE_LDLIBS) -o $@

.PHONY: benchmark-render-compute
benchmark-render-compute: compute $(COMPUTE_BIN_DIR)/benchmark_render$(EXEEXT)
	@$(COMPUTE_BIN_DIR)/benchmark_render$(EXEEXT)

$(COMPUTE_BIN_DIR)/benchmark_render_model$(EXEEXT): benchmarks/render_model.c $(LIB_COMPUTE) $(LIB_MCNP) $(LIB_CORE)
	@$(CC) $(CFLAGS) $(COMPUTE_CPPFLAGS) $(INCLUDES) $< \
		$(LIB_COMPUTE) $(LIB_MCNP) $(LIB_CORE) $(LDFLAGS) $(COMPUTE_LDLIBS) -o $@

$(COMPUTE_BIN_DIR)/benchmark_render_progressive$(EXEEXT): benchmarks/render_progressive.c $(LIB_COMPUTE) $(LIB_MCNP) $(LIB_CORE)
	@$(CC) $(CFLAGS) $(COMPUTE_CPPFLAGS) $(INCLUDES) $< \
		$(LIB_COMPUTE) $(LIB_MCNP) $(LIB_CORE) $(LDFLAGS) $(COMPUTE_LDLIBS) -o $@

$(COMPUTE_BIN_DIR)/benchmark_metal_counters$(EXEEXT): benchmarks/metal_counters.m | $(COMPUTE_BIN_DIR)
	@xcrun --sdk macosx clang $(CFLAGS) -fno-objc-arc $< $(COMPUTE_LDLIBS) -o $@

.PHONY: benchmark-metal-counters
benchmark-metal-counters: $(COMPUTE_BIN_DIR)/benchmark_metal_counters$(EXEEXT)
	@$(COMPUTE_BIN_DIR)/benchmark_metal_counters$(EXEEXT)

endif

ifeq ($(USE_METAL),1)
$(COMPUTE_BUILD_DIR)/ray_group_work.air: benchmarks/ray_group_work.metal $(TEST_DIR)/compute/ray_contract_kernels.metal $(COMPUTE_BUILD_DIR)/compute.air
	@xcrun --sdk macosx metal -std=macos-metal2.3 -fno-fast-math -fmodules-cache-path=$(abspath $(COMPUTE_BUILD_DIR))/module-cache -c $< -o $@

$(COMPUTE_BUILD_DIR)/ray_group_work.metallib: $(COMPUTE_BUILD_DIR)/ray_group_work.air
	@xcrun --sdk macosx metallib $< -o $@

$(COMPUTE_BUILD_DIR)/ray_group_work_metallib.h: $(COMPUTE_BUILD_DIR)/ray_group_work.metallib
	@xxd -i -n alea_ray_group_work_metallib $< > $@

$(COMPUTE_BIN_DIR)/benchmark_metal_ray_groups$(EXEEXT): benchmarks/metal_ray_groups.m benchmarks/ray_contract_oracle.h $(TEST_DIR)/compute/test_ray_contract.h benchmarks/gpu_geometry_fixture.h $(SRC_DIR)/compute/compute_metal.m $(SRC_DIR)/compute/compute_groups_metal.h $(COMPUTE_BUILD_DIR)/compute_metallib.h $(COMPUTE_BUILD_DIR)/ray_group_work_metallib.h $(LIB_COMPUTE) $(LIB_MCNP) $(LIB_CORE)
	@xcrun --sdk macosx clang $(CFLAGS) $(COMPUTE_CPPFLAGS) $(INCLUDES) -I$(COMPUTE_BUILD_DIR) -fno-objc-arc $< $(LIB_COMPUTE) $(LIB_MCNP) $(LIB_CORE) $(LDFLAGS) $(COMPUTE_LDLIBS) -o $@
endif

ifeq ($(USE_METAL),1)
$(COMPUTE_BUILD_DIR)/ray_group_execute.air: benchmarks/ray_group_execute.metal benchmarks/ray_group_work.metal $(TEST_DIR)/compute/ray_contract_kernels.metal $(COMPUTE_BUILD_DIR)/compute.air
	@xcrun --sdk macosx metal -std=macos-metal2.3 -fno-fast-math -fmodules-cache-path=$(abspath $(COMPUTE_BUILD_DIR))/module-cache -c $< -o $@

$(COMPUTE_BUILD_DIR)/ray_group_execute.metallib: $(COMPUTE_BUILD_DIR)/ray_group_execute.air
	@xcrun --sdk macosx metallib $< -o $@

$(COMPUTE_BUILD_DIR)/ray_group_execute_metallib.h: $(COMPUTE_BUILD_DIR)/ray_group_execute.metallib
	@xxd -i -n alea_ray_group_execute_metallib $< > $@

$(COMPUTE_BIN_DIR)/benchmark_metal_ray_execute$(EXEEXT): benchmarks/metal_ray_execute.m benchmarks/metal_ray_groups.m benchmarks/ray_execute_driver.h benchmarks/ray_lattice_fixture.h benchmarks/ray_contract_oracle.h $(TEST_DIR)/compute/test_ray_contract.h benchmarks/gpu_geometry_fixture.h $(SRC_DIR)/compute/compute_metal.m $(SRC_DIR)/compute/compute_groups_metal.h $(COMPUTE_BUILD_DIR)/compute_metallib.h $(COMPUTE_BUILD_DIR)/ray_group_execute_metallib.h $(LIB_COMPUTE) $(LIB_MCNP) $(LIB_CORE)
	@xcrun --sdk macosx clang $(CFLAGS) $(COMPUTE_CPPFLAGS) $(INCLUDES) -I$(COMPUTE_BUILD_DIR) -fno-objc-arc $< $(LIB_COMPUTE) $(LIB_MCNP) $(LIB_CORE) $(LDFLAGS) $(COMPUTE_LDLIBS) -o $@
endif

ifeq ($(USE_METAL),1)
$(COMPUTE_BUILD_DIR)/ray_eval_probe.air: $(TEST_DIR)/compute/ray_eval_probe.metal $(TEST_DIR)/compute/ray_contract_kernels.metal $(COMPUTE_BUILD_DIR)/compute.air
	@xcrun --sdk macosx metal -std=macos-metal2.3 -fno-fast-math -fmodules-cache-path=$(abspath $(COMPUTE_BUILD_DIR))/module-cache -c $< -o $@

$(COMPUTE_BUILD_DIR)/ray_eval_probe.metallib: $(COMPUTE_BUILD_DIR)/ray_eval_probe.air
	@xcrun --sdk macosx metallib $< -o $@

$(COMPUTE_BUILD_DIR)/ray_eval_probe_metallib.h: $(COMPUTE_BUILD_DIR)/ray_eval_probe.metallib
	@xxd -i -n alea_ray_eval_probe_metallib $< > $@
endif
