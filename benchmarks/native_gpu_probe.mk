# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0
# Standalone: make -f benchmarks/native_gpu_probe.mk [BACKEND=host|cuda|hip] test
# No GPU dependency discovery occurs in ordinary libalea builds.
BACKEND ?= host
OUT ?= build/native-gpu-probe/$(BACKEND)
CXX ?= c++
NVCC ?= nvcc
HIPCC ?= hipcc
PROBE_FLAGS ?= -O2 -std=c++11
ifeq ($(BACKEND),host)
PROBE_COMPILER = $(CXX)
PROBE_BACKEND_FLAGS = -Wall -Wextra -Werror
else ifeq ($(BACKEND),cuda)
PROBE_COMPILER = $(NVCC)
PROBE_BACKEND_FLAGS = -x cu -DPROBE_CUDA --fmad=false
else ifeq ($(BACKEND),hip)
PROBE_COMPILER = $(HIPCC)
PROBE_BACKEND_FLAGS = -x hip -DPROBE_HIP -ffp-contract=off
else
$(error BACKEND must be host, cuda or hip)
endif

.PHONY: all test run
all: $(OUT)/probe
$(OUT)/probe: benchmarks/native_gpu_probe.cpp benchmarks/native_probe_kernel.h benchmarks/native_gpu_probe.mk
	@mkdir -p $(OUT)
	$(PROBE_COMPILER) $(PROBE_FLAGS) $(PROBE_BACKEND_FLAGS) $< -o $@
test: $(OUT)/probe
	$(OUT)/probe --self-test
run: $(OUT)/probe
	$(OUT)/probe
