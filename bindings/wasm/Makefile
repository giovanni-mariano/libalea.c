# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

LIBALEA_ROOT := ../..
OUT_DIR := dist
SINGLE_BUILD_DIR := build/wasm-single-simd
SINGLE_BIN_DIR := bin/wasm-single-simd
OPENMP_BUILD_DIR := build/wasm-openmp-simd
OPENMP_BIN_DIR := bin/wasm-openmp-simd

EMCC ?= emcc
EMAR ?= emar
NATIVE_CC ?= cc
PTHREAD_POOL_SIZE ?= 8
WASM_INITIAL_MEMORY ?= 268435456
WASM_MAXIMUM_MEMORY ?= 2147483648
LARGE_INPUT_MB ?= 64
WASM_SIMD_FLAGS ?= -msimd128

EXPORTED_FUNCTIONS := [\
"_alea_wasm_init",\
"_alea_wasm_resize",\
"_alea_wasm_load_mcnp",\
"_alea_wasm_render",\
"_alea_wasm_pixels",\
"_alea_wasm_pixel_bytes",\
"_alea_wasm_width",\
"_alea_wasm_height",\
"_alea_wasm_cell_count",\
"_alea_wasm_surface_count",\
"_alea_wasm_center_x",\
"_alea_wasm_center_y",\
"_alea_wasm_center_z",\
"_alea_wasm_radius",\
"_alea_wasm_set_target",\
"_alea_wasm_openmp_enabled",\
"_alea_wasm_openmp_max_threads",\
"_alea_wasm_last_error",\
"_alea_wasm_destroy",\
"_malloc",\
"_free"\
]

COMMON_LINK_FLAGS := -O3 -std=c11 --no-entry \
	-sWASM=1 -sMODULARIZE=1 -sEXPORT_ES6=1 -sENVIRONMENT=web,worker,node \
	-sSTACK_SIZE=2097152 \
	-sALLOW_MEMORY_GROWTH=1 -sABORTING_MALLOC=0 \
	-sMAXIMUM_MEMORY=$(WASM_MAXIMUM_MEMORY) \
	-sEXPORTED_RUNTIME_METHODS=HEAPU8 \
	-sEXPORTED_FUNCTIONS='$(EXPORTED_FUNCTIONS)'

.PHONY: all single openmp check-openmp-toolchain demo test test-openmp test-large-openmp test-node test-openmp-node test-native single-lib openmp-lib clean help

all: single

single: $(OUT_DIR)/alea.js

openmp: check-openmp-toolchain $(OUT_DIR)/alea-openmp.js

check-openmp-toolchain:
	@version=`$(EMCC) --version | sed -n '1s/.* \([0-9][0-9.]*\) (.*/\1/p'`; \
	if [ -z "$$version" ] || ! printf '%s\n%s\n' 6.0.3 "$$version" | sort -V -C; then \
		echo "OpenMP WASM requires Emscripten 6.0.3 or newer (found $$version)" >&2; \
		exit 1; \
	fi

demo: single
	@echo "Serve this directory with: python3 serve.py"
	@echo "Then open: http://localhost:8000/web/"

test: test-native test-node

test-openmp: test-native test-openmp-node

test-large-openmp: openmp
	node tests/test_wasm.mjs alea-openmp.js openmp $(LARGE_INPUT_MB)

test-node: single
	node tests/test_wasm.mjs

test-openmp-node: openmp
	node tests/test_wasm.mjs alea-openmp.js openmp

single-lib:
	$(MAKE) -C $(LIBALEA_ROOT) AR=$(EMAR) \
		CC="$(EMCC) $(WASM_SIMD_FLAGS)" \
		BUILD_DIR=$(SINGLE_BUILD_DIR) BIN_DIR=$(SINGLE_BIN_DIR) \
		RELEASE=1 PORTABLE=1 lib-core $(SINGLE_BIN_DIR)/libalea_mcnp.a

openmp-lib:
	$(MAKE) -C $(LIBALEA_ROOT) AR=$(EMAR) \
		CC="$(EMCC) $(WASM_SIMD_FLAGS)" \
		BUILD_DIR=$(OPENMP_BUILD_DIR) BIN_DIR=$(OPENMP_BIN_DIR) \
		RELEASE=1 PORTABLE=1 USE_OPENMP=1 lib-core $(OPENMP_BIN_DIR)/libalea_mcnp.a

$(OUT_DIR):
	mkdir -p $@

$(OUT_DIR)/alea.js: Makefile alea_wasm.c alea_wasm.h single-lib | $(OUT_DIR)
	$(EMCC) $(WASM_SIMD_FLAGS) $(COMMON_LINK_FLAGS) \
		-I$(LIBALEA_ROOT)/include -I$(LIBALEA_ROOT)/src \
		alea_wasm.c $(LIBALEA_ROOT)/$(SINGLE_BIN_DIR)/libalea_mcnp.a \
		$(LIBALEA_ROOT)/$(SINGLE_BIN_DIR)/libalea.a -lm -o $@

$(OUT_DIR)/alea-openmp.js: Makefile alea_wasm.c alea_wasm.h openmp-lib | $(OUT_DIR)
	$(EMCC) $(WASM_SIMD_FLAGS) $(COMMON_LINK_FLAGS) -fopenmp -pthread \
		-sPTHREAD_POOL_SIZE=$(PTHREAD_POOL_SIZE) \
		-sINITIAL_MEMORY=$(WASM_INITIAL_MEMORY) \
		-I$(LIBALEA_ROOT)/include -I$(LIBALEA_ROOT)/src \
		alea_wasm.c $(LIBALEA_ROOT)/$(OPENMP_BIN_DIR)/libalea_mcnp.a \
		$(LIBALEA_ROOT)/$(OPENMP_BIN_DIR)/libalea.a -lm -o $@

test-native:
	$(MAKE) -C $(LIBALEA_ROOT) lib-core bin/libalea_mcnp.a
	$(NATIVE_CC) -O2 -std=c11 -Wall -Wextra -Werror -fopenmp \
		-I. -I$(LIBALEA_ROOT)/include -I$(LIBALEA_ROOT)/src \
		alea_wasm.c tests/test_wasm_binding.c \
		$(LIBALEA_ROOT)/bin/libalea_mcnp.a $(LIBALEA_ROOT)/bin/libalea.a \
		-lm -o tests/test_wasm_binding
	./tests/test_wasm_binding

clean:
	rm -rf $(OUT_DIR) tests/test_wasm_binding \
		$(LIBALEA_ROOT)/$(SINGLE_BUILD_DIR) $(LIBALEA_ROOT)/$(SINGLE_BIN_DIR) \
		$(LIBALEA_ROOT)/$(OPENMP_BUILD_DIR) $(LIBALEA_ROOT)/$(OPENMP_BIN_DIR)

help:
	@echo "Targets:"
	@echo "  single       Build the portable single-threaded WASM module"
	@echo "  openmp       Build the OpenMP/Web Worker WASM module"
	@echo "  demo         Build single-threaded WASM and print serving instructions"
	@echo "  test         Run native-facade and real WASM smoke tests"
	@echo "  test-openmp  Run native-facade and OpenMP WASM smoke tests"
	@echo "  test-large-openmp  Parse a large synthetic MCNP input (default: 64 MiB)"
	@echo "  test-node    Exercise the generated WASM module under Node.js"
	@echo "  test-openmp-node  Exercise the generated OpenMP module under Node.js"
	@echo "  test-native  Exercise the binding facade with the native library"
	@echo "  clean        Remove binding-specific generated files"
