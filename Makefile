# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

# Alea Project Makefile

# ============================================================================
# Compiler Configuration
# ============================================================================

CC = gcc
CFLAGS = -Wall -Wextra -g -std=c11 -fPIC
LDFLAGS = -lm

# Automatic dependency generation
DEPFLAGS = -MMD -MP

# Optional OpenMP support (set USE_OPENMP=1 to enable)
ifdef USE_OPENMP
  # Detect macOS (Apple clang needs different flags)
  UNAME_S := $(shell uname -s)
  ifeq ($(UNAME_S),Darwin)
    # macOS: use Homebrew libomp with -Xpreprocessor
    LIBOMP_PREFIX := $(shell brew --prefix libomp 2>/dev/null || echo "/usr/local/opt/libomp")
    CFLAGS += -Xpreprocessor -fopenmp -I$(LIBOMP_PREFIX)/include
    LDFLAGS += -L$(LIBOMP_PREFIX)/lib -lomp
  else
    # Linux/Windows: check if using clang (needs -lomp) or gcc (uses -fopenmp)
    IS_CLANG := $(shell $(CC) --version 2>/dev/null | grep -qi clang && echo 1)
    ifeq ($(IS_CLANG),1)
      CFLAGS += -fopenmp
      LDFLAGS += -fopenmp -lomp
    else
      CFLAGS += -fopenmp
      LDFLAGS += -fopenmp
    endif
  endif
endif

# Release build (set RELEASE=1)
ifdef RELEASE
    ifdef PORTABLE
        CFLAGS += -O3 -DNDEBUG
    else
        CFLAGS += -O3 -march=native -DNDEBUG
    endif
else
    CFLAGS += -O0
endif

# ============================================================================
# Directory Structure
# ============================================================================

# Source directories
SRC_DIR = src
CORE_DIR = $(SRC_DIR)/core
UTIL_DIR = $(SRC_DIR)/util
PRIMITIVES_DIR = $(SRC_DIR)/primitives
MCNP_PARSER_DIR = $(SRC_DIR)/mcnp/parser
MCNP_CONV_DIR = $(SRC_DIR)/mcnp/conversion
MCNP_EXPO_DIR = $(SRC_DIR)/mcnp/exporter
MCNP_MODEL_DIR = $(SRC_DIR)/mcnp

OPENMC_DIR = $(SRC_DIR)/openmc
SERPENT_DIR = $(SRC_DIR)/serpent
NUCDATA_DIR = $(SRC_DIR)/nucdata
RAYCAST_DIR = $(SRC_DIR)/raycast
SLICE_DIR = $(SRC_DIR)/slice
RENDER_DIR = $(SRC_DIR)/render
MESH_DIR = $(SRC_DIR)/mesh

# Lua (vendored)
LUA_DIR = vendor/lua
LUA_BIND_DIR = $(SRC_DIR)/lua_bind
LINENOISE_DIR = vendor/linenoise

# Build directories
BUILD_DIR = build
BIN_DIR = bin
INCLUDE_DIR = include

# Include paths
INCLUDES = -I$(INCLUDE_DIR) -I$(SRC_DIR)

# Test directory
TEST_DIR = tests
UNIT_TEST_DIR = $(TEST_DIR)/unit
INTEGRATION_TEST_DIR = $(TEST_DIR)/integration

# ============================================================================
# Source Files by Layer (NEW Flattened Architecture)
# ============================================================================


# Core system 
CORE_SRCS = \
	$(CORE_DIR)/alea_system.c \
	$(CORE_DIR)/alea_export.c \
	$(CORE_DIR)/alea_tolerance.c \
	$(CORE_DIR)/alea_primitive_dedup.c \
	$(CORE_DIR)/alea_materials.c \
	$(CORE_DIR)/alea_elements_data.c \
	$(CORE_DIR)/alea_error.c \
	$(CORE_DIR)/alea_public_api.c \
	$(CORE_DIR)/alea_eval.c \
	$(CORE_DIR)/alea_ops.c \
	$(CORE_DIR)/alea_simplify.c \
	$(CORE_DIR)/alea_universe.c \
	$(CORE_DIR)/alea_cell_complement.c \
	$(CORE_DIR)/alea_void.c \
	$(CORE_DIR)/alea_macrobody.c \
	$(CORE_DIR)/alea_spatial.c

# Memory management and utilities
UTIL_SRCS = \
	$(UTIL_DIR)/arena.c \
	$(UTIL_DIR)/str_builder.c \
	$(UTIL_DIR)/alea_log.c \
	$(UTIL_DIR)/poly_solve.c \
	$(UTIL_DIR)/compat.c \
	$(UTIL_DIR)/alea_svg.c

# Primitives
PRIMITIVES_SRCS = \
	$(PRIMITIVES_DIR)/primitive_eval.c \
	$(PRIMITIVES_DIR)/primitive_create.c \
	$(PRIMITIVES_DIR)/primitive_desc.c \
	$(PRIMITIVES_DIR)/bbox.c


	

# MCNP parser
MCNP_PARSER_SRCS = \
	$(MCNP_PARSER_DIR)/mcnp_parser.c \
	$(MCNP_PARSER_DIR)/mcnp_lexer.c \

# MCNP geometry parser
MCNP_GEOM_SRCS = \
	$(MCNP_PARSER_DIR)/geom_parser.c \
	$(MCNP_PARSER_DIR)/geom_lexer.c

# MCNP conversion (NEW simplified)
MCNP_CONV_SRCS = \
	$(MCNP_CONV_DIR)/mcnp_conversion.c \
	$(MCNP_CONV_DIR)/surface_conv.c \
	$(MCNP_CONV_DIR)/cell_conv.c

MCNP_EXPO_SRCS = \
	$(MCNP_EXPO_DIR)/mcnp_str.c \
	$(MCNP_EXPO_DIR)/mcnp_export.c

MCNP_MODEL_SRCS = \
	$(MCNP_MODEL_DIR)/mcnp_model.c

# Nuclear data module
NUCDATA_SRCS = \
	$(NUCDATA_DIR)/context.c \
	$(NUCDATA_DIR)/ace_reader.c \
	$(NUCDATA_DIR)/xsdir.c \
	$(NUCDATA_DIR)/xs_decode.c \
	$(NUCDATA_DIR)/lookup.c \
	$(NUCDATA_DIR)/material.c \
	$(NUCDATA_DIR)/reaction.c \
	$(NUCDATA_DIR)/angular.c \
	$(NUCDATA_DIR)/energy_dist.c \
	$(NUCDATA_DIR)/doppler.c \
	$(NUCDATA_DIR)/multigroup.c \
	$(NUCDATA_DIR)/material_bridge.c

OPENMC_EXPO_SRCS = \
	$(OPENMC_DIR)/openmc_xml.c \
	$(OPENMC_DIR)/openmc_export.c

SERPENT_EXPO_SRCS = \
	$(SERPENT_DIR)/serpent_export.c

# OpenMC parser and conversion
OPENMC_PARSE_SRCS = \
	$(OPENMC_DIR)/openmc_parse.c \
	$(OPENMC_DIR)/openmc_region.c \
	$(OPENMC_DIR)/openmc_conversion.c

# Raycast module
RAYCAST_SRCS = \
	$(RAYCAST_DIR)/raycast.c \
	$(RAYCAST_DIR)/ray_intersect.c \
	$(RAYCAST_DIR)/bvh.c \
	$(RAYCAST_DIR)/raycast_api.c

# Slice module (analytical intersection + vector export)
SLICE_SRCS = \
	$(SLICE_DIR)/curve_intersect.c \
	$(SLICE_DIR)/slice_api.c

# Render module (3D batch renderer)
RENDER_SRCS = \
	$(RENDER_DIR)/render3d.c

# Mesh export module
MESH_SRCS = \
	$(MESH_DIR)/mesh_export.c

# Core library sources (geometry engine + raycast/slice/render/mesh)
CORE_LIB_SRCS = \
	$(CORE_SRCS) \
	$(UTIL_SRCS) \
	$(PRIMITIVES_SRCS) \
	$(RAYCAST_SRCS) \
	$(SLICE_SRCS) \
	$(RENDER_SRCS) \
	$(MESH_SRCS)

# MCNP module sources (parser + conversion + exporter + model)
MCNP_MODULE_SRCS = $(MCNP_PARSER_SRCS) $(MCNP_GEOM_SRCS) $(MCNP_CONV_SRCS) $(MCNP_EXPO_SRCS) $(MCNP_MODEL_SRCS)

# OpenMC model
OPENMC_MODEL_SRCS = $(OPENMC_DIR)/openmc_model.c

# OpenMC module sources (parser + conversion + exporter + model)
OPENMC_MODULE_SRCS = $(OPENMC_EXPO_SRCS) $(OPENMC_PARSE_SRCS) $(OPENMC_MODEL_SRCS)

# Serpent module sources (exporter only)
SERPENT_MODULE_SRCS = $(SERPENT_EXPO_SRCS)

# Lua 5.4 (vendored) - exclude standalone binaries
LUA_SRCS = $(filter-out $(LUA_DIR)/lua.c $(LUA_DIR)/luac.c $(LUA_DIR)/onelua.c, $(wildcard $(LUA_DIR)/*.c))
LUA_OBJS = $(patsubst $(LUA_DIR)/%.c,$(BUILD_DIR)/lua/%.o,$(LUA_SRCS))

# Lua bindings
LUA_BIND_SRCS = $(filter-out $(LUA_BIND_DIR)/lua_main.c, $(wildcard $(LUA_BIND_DIR)/*.c))
LUA_BIND_OBJS = $(LUA_BIND_SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)

# ============================================================================
# Object Files
# ============================================================================

CORE_OBJS = $(CORE_SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)
UTIL_OBJS = $(UTIL_SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)
PRIMITIVES_OBJS = $(PRIMITIVES_SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)
MCNP_PARSER_OBJS = $(MCNP_PARSER_SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)
MCNP_GEOM_OBJS = $(MCNP_GEOM_SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)
MCNP_CONV_OBJS = $(MCNP_CONV_SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)
MCNP_EXPO_OBJS = $(MCNP_EXPO_SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)
MCNP_MODEL_OBJS = $(MCNP_MODEL_SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)
NUCDATA_OBJS = $(NUCDATA_SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)
OPENMC_EXPO_OBJS = $(OPENMC_EXPO_SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)
SERPENT_EXPO_OBJS = $(SERPENT_EXPO_SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)
OPENMC_PARSE_OBJS = $(OPENMC_PARSE_SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)
OPENMC_MODEL_OBJS = $(OPENMC_MODEL_SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)
RAYCAST_OBJS = $(RAYCAST_SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)
SLICE_OBJS = $(SLICE_SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)
RENDER_OBJS = $(RENDER_SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)
MESH_OBJS = $(MESH_SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)

# Core library objects (geometry engine + raycast/slice/render/mesh)
CORE_LIB_OBJS = $(CORE_OBJS) $(UTIL_OBJS) $(PRIMITIVES_OBJS) $(RAYCAST_OBJS) $(SLICE_OBJS) $(RENDER_OBJS) $(MESH_OBJS)

# MCNP module objects
MCNP_MODULE_OBJS = $(MCNP_PARSER_OBJS) $(MCNP_GEOM_OBJS) $(MCNP_CONV_OBJS) $(MCNP_EXPO_OBJS) $(MCNP_MODEL_OBJS)

# OpenMC module objects
OPENMC_MODULE_OBJS = $(OPENMC_EXPO_OBJS) $(OPENMC_PARSE_OBJS) $(OPENMC_MODEL_OBJS)

# Serpent module objects
SERPENT_MODULE_OBJS = $(SERPENT_EXPO_OBJS)

# All objects (for full library)
ALL_OBJS = $(CORE_LIB_OBJS) $(NUCDATA_OBJS) $(MCNP_MODULE_OBJS) $(OPENMC_MODULE_OBJS) $(SERPENT_MODULE_OBJS)

# Dependency files (generated by -MMD)
ALL_DEPS = $(ALL_OBJS:.o=.d)

# ============================================================================
# Libraries
# ============================================================================

# Core library (geometry engine + raycast/slice/render/mesh)
LIB_CORE = $(BIN_DIR)/libalea.a

# Format modules (optional)
LIB_MCNP = $(BIN_DIR)/libalea_mcnp.a
LIB_OPENMC = $(BIN_DIR)/libalea_openmc.a
LIB_SERPENT = $(BIN_DIR)/libalea_serpent.a
LIB_NUCDATA = $(BIN_DIR)/libalea_nucdata.a

# Full library (everything)
LIB = $(BIN_DIR)/libalea_full.a

# ============================================================================
# Test Sources
# ============================================================================

# Unit tests
UNIT_TEST_SRCS = $(wildcard $(UNIT_TEST_DIR)/*.c)
UNIT_TEST_BINS = $(UNIT_TEST_SRCS:$(UNIT_TEST_DIR)/%.c=$(BIN_DIR)/tests/unit/%)

# Integration tests
INTEGRATION_TEST_SRCS = $(wildcard $(INTEGRATION_TEST_DIR)/*.c)
INTEGRATION_TEST_BINS = $(INTEGRATION_TEST_SRCS:$(INTEGRATION_TEST_DIR)/%.c=$(BIN_DIR)/tests/integration/%)

# All tests
ALL_TEST_BINS = $(UNIT_TEST_BINS) $(INTEGRATION_TEST_BINS)

# ============================================================================
# Main Targets
# ============================================================================

.PHONY: all clean full lib-core modules tests structure help test cli test-lua tools

# Default target: core library only
all: lib-core

# Build core library (recommended)
lib-core: structure $(LIB_CORE)

# Build optional format modules
modules: structure $(LIB_MCNP) $(LIB_OPENMC) $(LIB_SERPENT) $(LIB_NUCDATA)

# Build everything (core + modules + full archive)
full: lib-core modules $(LIB)

# CLI binary
ALEA_CLI = $(BIN_DIR)/alea
cli: lib-core modules $(ALEA_CLI)

# Build tools (mc_convert, mc_plotter)
tools: lib-core modules
	$(MAKE) -C tools

# Build all tests (requires core + modules)
tests: lib-core modules $(ALL_TEST_BINS)

# ============================================================================
# Directory Structure Creation
# ============================================================================

# Directory creation rules (order-only prerequisites for parallel safety)
BUILD_DIRS = $(BUILD_DIR)/core $(BUILD_DIR)/util $(BUILD_DIR)/primitives \
	$(BUILD_DIR)/nucdata \
	$(BUILD_DIR)/mcnp/parser $(BUILD_DIR)/mcnp/geometry $(BUILD_DIR)/mcnp/conversion \
	$(BUILD_DIR)/mcnp/exporter $(BUILD_DIR)/mcnp $(BUILD_DIR)/openmc $(BUILD_DIR)/serpent \
	$(BUILD_DIR)/raycast $(BUILD_DIR)/slice $(BUILD_DIR)/render $(BUILD_DIR)/mesh \
	$(BUILD_DIR)/lua $(BUILD_DIR)/lua_bind $(BUILD_DIR)/linenoise \
	$(BIN_DIR) $(BIN_DIR)/tests/unit $(BIN_DIR)/tests/integration

$(BUILD_DIRS):
	@mkdir -p $@

structure: $(BUILD_DIRS)

# ============================================================================
# Library Build Rules
# ============================================================================

# Core library (default)
$(LIB_CORE): $(CORE_LIB_OBJS) | $(BIN_DIR)
	@echo "AR  $@"
	@ar rcs $@ $^
	@echo "✓ Built library: $@"

# MCNP module (optional)
$(LIB_MCNP): $(MCNP_MODULE_OBJS) | $(BIN_DIR)
	@echo "AR  $@"
	@ar rcs $@ $^
	@echo "✓ Built MCNP module: $@"

# Nuclear data module (optional)
$(LIB_NUCDATA): $(NUCDATA_OBJS) | $(BIN_DIR)
	@echo "AR  $@"
	@ar rcs $@ $^
	@echo "✓ Built nucdata module: $@"

# OpenMC module (optional)
$(LIB_OPENMC): $(OPENMC_MODULE_OBJS) | $(BIN_DIR)
	@echo "AR  $@"
	@ar rcs $@ $^
	@echo "✓ Built OpenMC module: $@"

# Serpent module (optional)
$(LIB_SERPENT): $(SERPENT_MODULE_OBJS) | $(BIN_DIR)
	@echo "AR  $@"
	@ar rcs $@ $^
	@echo "✓ Built Serpent module: $@"

# Full library (core + all format modules)
$(LIB): $(ALL_OBJS) | $(BIN_DIR)
	@echo "AR  $@"
	@ar rcs $@ $^
	@echo "✓ Built full library: $@"

# ============================================================================
# Object File Rules by Layer
# ============================================================================



# Core system
$(BUILD_DIR)/core/%.o: $(CORE_DIR)/%.c | $(BUILD_DIR)/core
	@echo "CC  $<"
	@$(CC) $(CFLAGS) $(DEPFLAGS) $(INCLUDES) -c $< -o $@

# Utilities
$(BUILD_DIR)/util/%.o: $(UTIL_DIR)/%.c | $(BUILD_DIR)/util
	@echo "CC  $<"
	@$(CC) $(CFLAGS) $(DEPFLAGS) $(INCLUDES) -c $< -o $@

# Primitives
$(BUILD_DIR)/primitives/%.o: $(PRIMITIVES_DIR)/%.c | $(BUILD_DIR)/primitives
	@echo "CC  $<"
	@$(CC) $(CFLAGS) $(DEPFLAGS) $(INCLUDES) -c $< -o $@

# MCNP parser
$(BUILD_DIR)/mcnp/parser/%.o: $(MCNP_PARSER_DIR)/%.c | $(BUILD_DIR)/mcnp/parser
	@echo "CC  $<"
	@$(CC) $(CFLAGS) $(DEPFLAGS) $(INCLUDES) -c $< -o $@

# MCNP conversion
$(BUILD_DIR)/mcnp/conversion/%.o: $(MCNP_CONV_DIR)/%.c | $(BUILD_DIR)/mcnp/conversion
	@echo "CC  $<"
	@$(CC) $(CFLAGS) $(DEPFLAGS) $(INCLUDES) -c $< -o $@

# MCNP exporter
$(BUILD_DIR)/mcnp/exporter/%.o: $(MCNP_EXPO_DIR)/%.c | $(BUILD_DIR)/mcnp/exporter
	@echo "CC  $<"
	@$(CC) $(CFLAGS) $(DEPFLAGS) $(INCLUDES) -c $< -o $@

# MCNP model
$(BUILD_DIR)/mcnp/%.o: $(MCNP_MODEL_DIR)/%.c | $(BUILD_DIR)/mcnp
	@echo "CC  $<"
	@$(CC) $(CFLAGS) $(DEPFLAGS) $(INCLUDES) -c $< -o $@

# OpenMC
$(BUILD_DIR)/openmc/%.o: $(OPENMC_DIR)/%.c | $(BUILD_DIR)/openmc
	@echo "CC  $<"
	@$(CC) $(CFLAGS) $(DEPFLAGS) $(INCLUDES) -c $< -o $@

# Serpent
$(BUILD_DIR)/serpent/%.o: $(SERPENT_DIR)/%.c | $(BUILD_DIR)/serpent
	@echo "CC  $<"
	@$(CC) $(CFLAGS) $(DEPFLAGS) $(INCLUDES) -c $< -o $@

# Nuclear data
$(BUILD_DIR)/nucdata/%.o: $(NUCDATA_DIR)/%.c | $(BUILD_DIR)/nucdata
	@echo "CC  $<"
	@$(CC) $(CFLAGS) $(DEPFLAGS) $(INCLUDES) -c $< -o $@

# Raycast
$(BUILD_DIR)/raycast/%.o: $(RAYCAST_DIR)/%.c | $(BUILD_DIR)/raycast
	@echo "CC  $<"
	@$(CC) $(CFLAGS) $(DEPFLAGS) $(INCLUDES) -c $< -o $@

# Slice (analytical intersection + vector export)
$(BUILD_DIR)/slice/%.o: $(SLICE_DIR)/%.c | $(BUILD_DIR)/slice
	@echo "CC  $<"
	@$(CC) $(CFLAGS) $(DEPFLAGS) $(INCLUDES) -c $< -o $@

# Render (3D batch renderer)
$(BUILD_DIR)/render/%.o: $(RENDER_DIR)/%.c | $(BUILD_DIR)/render
	@echo "CC  $<"
	@$(CC) $(CFLAGS) $(DEPFLAGS) $(INCLUDES) -c $< -o $@

# Mesh export
$(BUILD_DIR)/mesh/%.o: $(MESH_DIR)/%.c | $(BUILD_DIR)/mesh
	@echo "CC  $<"
	@$(CC) $(CFLAGS) $(DEPFLAGS) $(INCLUDES) -c $< -o $@

# Lua 5.4 (vendored) - suppress warnings with -w
# Use LUA_USE_POSIX on Unix, LUA_USE_WINDOWS on Windows
LUA_UNAME_S := $(shell uname -s)
ifeq ($(OS),Windows_NT)
  LUA_PLAT_FLAGS = -DLUA_USE_WINDOWS
else ifneq ($(findstring MINGW,$(LUA_UNAME_S)),)
  LUA_PLAT_FLAGS = -DLUA_USE_WINDOWS
else ifneq ($(findstring MSYS,$(LUA_UNAME_S)),)
  LUA_PLAT_FLAGS = -DLUA_USE_WINDOWS
else
  LUA_PLAT_FLAGS = -DLUA_USE_POSIX
endif

$(BUILD_DIR)/lua/%.o: $(LUA_DIR)/%.c | $(BUILD_DIR)/lua
	@echo "CC  $< (lua)"
	@$(CC) -std=gnu11 -O2 -fPIC $(LUA_PLAT_FLAGS) -w -MMD -MP -c $< -o $@

# Lua bindings
$(BUILD_DIR)/lua_bind/%.o: $(LUA_BIND_DIR)/%.c | $(BUILD_DIR)/lua_bind
	@echo "CC  $<"
	@$(CC) $(CFLAGS) $(DEPFLAGS) $(INCLUDES) -I$(LUA_DIR) -I$(LINENOISE_DIR) -c $< -o $@

# Linenoise (POSIX only - provides line editing in REPL)
# Windows uses fgets fallback in lua_main.c
ifeq ($(OS),Windows_NT)
  LINENOISE_OBJ =
  LINENOISE_INC =
  CLI_LDFLAGS = $(LDFLAGS)
else
  UNAME_S := $(shell uname -s)
  ifneq ($(findstring MINGW,$(UNAME_S)),)
    LINENOISE_OBJ =
    LINENOISE_INC =
    CLI_LDFLAGS = $(LDFLAGS)
  else ifneq ($(findstring MSYS,$(UNAME_S)),)
    LINENOISE_OBJ =
    LINENOISE_INC =
    CLI_LDFLAGS = $(LDFLAGS)
  else
    LINENOISE_OBJ = $(BUILD_DIR)/linenoise/linenoise.o
    LINENOISE_INC = -I$(LINENOISE_DIR)
    CLI_LDFLAGS = $(LDFLAGS) -ldl
  endif
endif

$(BUILD_DIR)/linenoise/linenoise.o: $(LINENOISE_DIR)/linenoise.c | $(BUILD_DIR)/linenoise
	@echo "CC  $< (linenoise)"
	@$(CC) -std=gnu11 -O2 -fPIC -w -MMD -MP -c $< -o $@

# CLI binary (lua_main.c + bindings + linenoise + all alea libs + lua)
# Platform-specific whole-archive flags for format modules
CLI_UNAME_S := $(shell uname -s)
ifeq ($(CLI_UNAME_S),Darwin)
  CLI_LIBS = -Wl,-force_load,$(LIB_MCNP) -Wl,-force_load,$(LIB_OPENMC) -Wl,-force_load,$(LIB_SERPENT) -Wl,-force_load,$(LIB_NUCDATA) $(LIB_CORE)
else
  CLI_LIBS = -Wl,--whole-archive $(LIB_MCNP) $(LIB_OPENMC) $(LIB_SERPENT) $(LIB_NUCDATA) -Wl,--no-whole-archive $(LIB_CORE)
endif

$(ALEA_CLI): $(LUA_BIND_DIR)/lua_main.c $(LUA_BIND_OBJS) $(LINENOISE_OBJ) $(LUA_OBJS) $(LIB_CORE) $(LIB_MCNP) $(LIB_OPENMC) $(LIB_SERPENT) $(LIB_NUCDATA) | $(BIN_DIR)
	@echo "LD  $@"
	@$(CC) $(CFLAGS) $(INCLUDES) -I$(LUA_DIR) $(LINENOISE_INC) $< $(LUA_BIND_OBJS) $(LINENOISE_OBJ) $(LUA_OBJS) \
		$(CLI_LIBS) $(CLI_LDFLAGS) -o $@

# ============================================================================
# Test Build Rules
# ============================================================================

# Test include path (adds tests/ for alea_test.h)
# Strip -DNDEBUG from test builds so assert() remains active
TEST_INCLUDES = $(INCLUDES) -I$(TEST_DIR)
TEST_CFLAGS = $(filter-out -DNDEBUG,$(CFLAGS))

# Platform-specific whole-archive flags for format modules
TEST_UNAME_S := $(shell uname -s)
ifeq ($(TEST_UNAME_S),Darwin)
  TEST_LIBS = -Wl,-force_load,$(LIB_MCNP) -Wl,-force_load,$(LIB_OPENMC) -Wl,-force_load,$(LIB_SERPENT) -Wl,-force_load,$(LIB_NUCDATA) $(LIB_CORE)
else
  TEST_LIBS = -Wl,--whole-archive $(LIB_MCNP) $(LIB_OPENMC) $(LIB_SERPENT) $(LIB_NUCDATA) -Wl,--no-whole-archive $(LIB_CORE)
endif

# Unit tests
$(BIN_DIR)/tests/unit/%: $(UNIT_TEST_DIR)/%.c $(LIB_CORE) $(LIB_MCNP) $(LIB_OPENMC) $(LIB_SERPENT) $(LIB_NUCDATA) | $(BIN_DIR)/tests/unit
	@echo "LD  $@"
	@$(CC) $(TEST_CFLAGS) $(TEST_INCLUDES) $< $(TEST_LIBS) $(LDFLAGS) -o $@

# Integration tests
$(BIN_DIR)/tests/integration/%: $(INTEGRATION_TEST_DIR)/%.c $(LIB_CORE) $(LIB_MCNP) $(LIB_OPENMC) $(LIB_SERPENT) $(LIB_NUCDATA) | $(BIN_DIR)/tests/integration
	@echo "LD  $@"
	@$(CC) $(TEST_CFLAGS) $(TEST_INCLUDES) $< $(TEST_LIBS) $(LDFLAGS) -o $@

# ============================================================================
# Run Tests
# ============================================================================

test: tests
	@echo ""
	@echo "=== Running Unit Tests ==="
	@for test in $(UNIT_TEST_BINS); do \
		if [ -f $$test ]; then \
			echo ""; \
			echo "Running $$test..."; \
			LD_LIBRARY_PATH=$(BIN_DIR) $$test || exit 1; \
		fi \
	done
	@echo ""
	@echo "=== Running Integration Tests ==="
	@for test in $(INTEGRATION_TEST_BINS); do \
		if [ -f $$test ]; then \
			echo ""; \
			echo "Running $$test..."; \
			LD_LIBRARY_PATH=$(BIN_DIR) $$test || exit 1; \
		fi \
	done
	@echo ""
	@echo "✓ All tests passed!"

# Run only unit tests
.PHONY: test-unit
test-unit: $(UNIT_TEST_BINS)
	@echo ""
	@echo "=== Running Unit Tests ==="
	@for test in $(UNIT_TEST_BINS); do \
		if [ -f $$test ]; then \
			echo ""; \
			echo "Running $$test..."; \
			LD_LIBRARY_PATH=$(BIN_DIR) $$test || exit 1; \
		fi \
	done
	@echo ""
	@echo "✓ All unit tests passed!"

# Run only integration tests
.PHONY: test-integration
test-integration: $(INTEGRATION_TEST_BINS)
	@echo ""
	@echo "=== Running Integration Tests ==="
	@for test in $(INTEGRATION_TEST_BINS); do \
		if [ -f $$test ]; then \
			echo ""; \
			echo "Running $$test..."; \
			LD_LIBRARY_PATH=$(BIN_DIR) $$test || exit 1; \
		fi \
	done
	@echo ""
	@echo "✓ All integration tests passed!"


# Run Lua tests
test-lua: cli
	@echo ""
	@echo "=== Running Lua Tests ==="
	@for test in tests/lua/test_*.lua; do \
		if [ -f $$test ]; then \
			echo ""; \
			echo "Running $$test..."; \
			$(ALEA_CLI) $$test || exit 1; \
		fi \
	done
	@echo ""
	@echo "✓ All Lua tests passed!"

# ============================================================================
# Valgrind (memory leak and error detection)
# ============================================================================

VALGRIND = valgrind --leak-check=full --errors-for-leak-kinds=definite,indirect,possible \
	--track-origins=yes --error-exitcode=1 \
	--suppressions=valgrind.supp

.PHONY: test-valgrind
test-valgrind: tests cli
	@echo ""
	@echo "=== Running Unit Tests under Valgrind ==="
	@for test in $(UNIT_TEST_BINS); do \
		if [ -f $$test ]; then \
			echo ""; \
			echo "Running $$test..."; \
			LD_LIBRARY_PATH=$(BIN_DIR) $(VALGRIND) $$test || exit 1; \
		fi \
	done
	@echo ""
	@echo "=== Running Integration Tests under Valgrind ==="
	@for test in $(INTEGRATION_TEST_BINS); do \
		if [ -f $$test ]; then \
			echo ""; \
			echo "Running $$test..."; \
			LD_LIBRARY_PATH=$(BIN_DIR) $(VALGRIND) $$test || exit 1; \
		fi \
	done
	@echo ""
	@echo "=== Running Lua Tests under Valgrind ==="
	@for test in tests/lua/test_*.lua; do \
		if [ -f $$test ]; then \
			echo ""; \
			echo "Running $$test..."; \
			$(VALGRIND) $(ALEA_CLI) $$test || exit 1; \
		fi \
	done
	@echo ""
	@echo "All tests passed under Valgrind (no leaks, no errors)!"

# ============================================================================
# Fuzzing (requires clang with libFuzzer)
# ============================================================================

FUZZ_DIR = tests/fuzz
FUZZ_CC = clang
FUZZ_CFLAGS = -g -O1 -fsanitize=fuzzer,address,undefined -fno-omit-frame-pointer

# MCNP parser fuzzer
$(BIN_DIR)/fuzz_mcnp: $(FUZZ_DIR)/fuzz_mcnp_parser.c $(LIB) | $(BIN_DIR)
	@echo "FUZZ  $@"
	@$(FUZZ_CC) $(FUZZ_CFLAGS) $(INCLUDES) $< $(LIB) $(LDFLAGS) -o $@

# OpenMC parser fuzzer
$(BIN_DIR)/fuzz_openmc: $(FUZZ_DIR)/fuzz_openmc_parser.c $(LIB) | $(BIN_DIR)
	@echo "FUZZ  $@"
	@$(FUZZ_CC) $(FUZZ_CFLAGS) $(INCLUDES) $< $(LIB) $(LDFLAGS) -o $@

# Legacy alias
$(BIN_DIR)/fuzz_parser: $(BIN_DIR)/fuzz_mcnp
	@ln -sf fuzz_mcnp $@

# Build all fuzzers
fuzz-build: $(BIN_DIR)/fuzz_mcnp $(BIN_DIR)/fuzz_openmc
	@echo "✓ Fuzzers built:"
	@echo "    $(BIN_DIR)/fuzz_mcnp    - MCNP parser fuzzer"
	@echo "    $(BIN_DIR)/fuzz_openmc  - OpenMC parser fuzzer"
	@echo ""
	@echo "Run with:"
	@echo "    make fuzz-mcnp    - Fuzz MCNP parser"
	@echo "    make fuzz-openmc  - Fuzz OpenMC parser"

# Fuzz MCNP parser
fuzz-mcnp: $(BIN_DIR)/fuzz_mcnp
	@echo "Starting MCNP fuzzer (Ctrl+C to stop)..."
	@echo "Crashes saved to fuzz/crashes/"
	@mkdir -p $(FUZZ_DIR)/crashes
	@$(BIN_DIR)/fuzz_mcnp $(FUZZ_DIR)/corpus/ -artifact_prefix=$(FUZZ_DIR)/crashes/ -max_len=65536

# Fuzz OpenMC parser
fuzz-openmc: $(BIN_DIR)/fuzz_openmc
	@echo "Starting OpenMC fuzzer (Ctrl+C to stop)..."
	@echo "Crashes saved to fuzz/crashes/"
	@mkdir -p $(FUZZ_DIR)/crashes
	@$(BIN_DIR)/fuzz_openmc $(FUZZ_DIR)/corpus_openmc/ -artifact_prefix=$(FUZZ_DIR)/crashes/ -max_len=65536

# Default fuzz target (MCNP)
fuzz: fuzz-mcnp

# ============================================================================
# Clean Targets
# ============================================================================

clean:
	@echo "Cleaning build artifacts..."
	@rm -rf $(BUILD_DIR) $(BIN_DIR)
	@$(MAKE) -C tools clean 2>/dev/null || true
	@echo "✓ Clean complete"

distclean: clean
	@echo "Deep cleaning..."
	@find . -name "*.o" -delete
	@find . -name "*.a" -delete
	@find . -name "*~" -delete
	@echo "✓ Distclean complete"

# ============================================================================
# Show Project Structure
# ============================================================================

.PHONY: tree
tree:
	@echo "Project structure:"
	@tree -I 'obj|bin|__pycache__|*.pyc' || ls -R

# ============================================================================
# Help
# ============================================================================

help:
	@echo "Alea Project Makefile"
	@echo ""
	@echo "Main Targets:"
	@echo "  all              - Build library (default)"
	@echo "  full             - Build full library (core + MCNP + OpenMC)"
	@echo "  tests            - Build all tests"
	@echo "  test             - Build and run all tests"
	@echo "  test-unit        - Run only unit tests"
	@echo "  test-integration - Run only integration tests"
	@echo "  clean            - Remove build artifacts"
	@echo "  distclean        - Deep clean (all generated files)"
	@echo "  tree             - Show project structure"
	@echo "  help             - Show this help message"
	@echo ""
	@echo "Directory Layout:"
	@echo "  src/             - Source code"
	@echo "    core/          - Core CSG system"
	@echo "    util/          - Memory management"
	@echo "    primitives/    - Geometric primitives"
	@echo "    csg/           - CSG operations"
	@echo "    mcnp/          - MCNP parsing and conversion"
	@echo "  include/         - Public headers"
	@echo "  obj/             - Compiled object files"
	@echo "  bin/             - Executables and libraries"
	@echo "  tests/           - Test code"
	@echo "    unit_/         - Unit tests"
	@echo "    integration/   - Integration tests"
	@echo ""
	@echo "Architecture:"
	@echo "  NEW: Flattened single alea_system_t structure"
	@echo "  - Replaces old node_pool, primitive_pool, alea_world"
	@echo "  - Direct array access for performance"
	@echo "  - Built-in deduplication with hash table"
	@echo "  - Automatic growth (no manual sizing)"

# ============================================================================
# Automatic Dependencies
# ============================================================================

# Include generated dependency files (ignore if they don't exist yet)
-include $(ALL_DEPS)

.PHONY: all full tests test test-unit test-integration test-valgrind clean distclean tree help structure
