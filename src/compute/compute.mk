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
  COMPUTE_OBJS += $(COMPUTE_BUILD_DIR)/compute_scene.o $(COMPUTE_BUILD_DIR)/compute_render_scene.o $(COMPUTE_BUILD_DIR)/compute_bvh.o $(COMPUTE_BUILD_DIR)/compute_spatial.o $(COMPUTE_BUILD_DIR)/compute_local.o $(COMPUTE_BUILD_DIR)/compute_point.o $(COMPUTE_BUILD_DIR)/compute_profile.o $(COMPUTE_BUILD_DIR)/compute_metal.o
endif

.PHONY: compute test-compute benchmark-compute install-compute
compute: lib-core $(LIB_COMPUTE)

$(COMPUTE_BUILD_DIR) $(COMPUTE_BIN_DIR):
	@$(MKDIR_P) $@

$(COMPUTE_BUILD_DIR)/%.o: $(SRC_DIR)/compute/%.c | $(COMPUTE_BUILD_DIR)
	@echo "CC  $< ($(COMPUTE_VARIANT))"
	@$(CC) $(CFLAGS) $(DEPFLAGS) $(COMPUTE_CPPFLAGS) $(INCLUDES) -c $< -o $@

$(COMPUTE_BUILD_DIR)/compute.air: $(SRC_DIR)/compute/compute.metal $(SRC_DIR)/compute/classify.metal $(SRC_DIR)/compute/render.metal $(SRC_DIR)/compute/resume.metal $(SRC_DIR)/compute/intersect.metal $(SRC_DIR)/compute/compute_layout.h | $(COMPUTE_BUILD_DIR)
	@echo "METAL $<"
	@xcrun --sdk macosx metal -std=macos-metal2.3 -fno-fast-math \
		-fmodules-cache-path=$(abspath $(COMPUTE_BUILD_DIR))/module-cache -c $< -o $@

$(COMPUTE_BUILD_DIR)/compute.metallib: $(COMPUTE_BUILD_DIR)/compute.air
	@xcrun --sdk macosx metallib $< -o $@

$(COMPUTE_BUILD_DIR)/compute_metallib.h: $(COMPUTE_BUILD_DIR)/compute.metallib
	@xxd -i -n alea_compute_metallib $< > $@

$(COMPUTE_BUILD_DIR)/compute_metal.o: $(SRC_DIR)/compute/compute_metal.m $(COMPUTE_BUILD_DIR)/compute_metallib.h
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

$(COMPUTE_BIN_DIR)/test_compute_builders: $(TEST_DIR)/compute/test_compute_builders.c $(SRC_DIR)/compute/compute_bvh.c $(SRC_DIR)/compute/compute_spatial.c $(SRC_DIR)/compute/compute_point.c $(LIB_COMPUTE) $(LIB_CORE)
	@$(CC) $(CFLAGS) $(COMPUTE_CPPFLAGS) $(INCLUDES) -I$(TEST_DIR) $< \
		$(LIB_COMPUTE) $(LIB_CORE) $(LDFLAGS) $(COMPUTE_LDLIBS) -o $@

$(COMPUTE_BIN_DIR)/test_metal_resources: $(TEST_DIR)/compute/test_metal_resources.m $(SRC_DIR)/compute/compute_metal.m $(COMPUTE_BUILD_DIR)/compute_metallib.h $(LIB_COMPUTE) $(LIB_CORE)
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
$(COMPUTE_BIN_DIR)/benchmark_gpu_geometry$(EXEEXT): benchmarks/gpu_geometry.c $(LIB_COMPUTE) $(LIB_CORE)
	@$(CC) $(CFLAGS) $(COMPUTE_CPPFLAGS) $(INCLUDES) $< \
		$(LIB_COMPUTE) $(LIB_CORE) $(LDFLAGS) $(COMPUTE_LDLIBS) -o $@

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

endif
