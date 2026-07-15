# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

# Tools Makefile

CC = gcc
CFLAGS ?= -Wall -Wextra -g -std=c11
# Match the library build's include roots: public headers (e.g. alea_mcnp.h)
# pull in src-side utility headers such as util/alea_vec.h.
INCLUDES = -I../include -I../src
BIN_DIR = ../bin
LDFLAGS ?= -lm

ifdef WINDOWS_GNU
  EXEEXT ?= .exe
else
  EXEEXT ?=
endif

# Detect platform
UNAME_S := $(shell uname -s)

# OpenMP support (pass USE_OPENMP=1 to enable)
ifdef USE_OPENMP
  ifeq ($(UNAME_S),Darwin)
    LIBOMP_PREFIX ?= $(shell brew --prefix libomp 2>/dev/null || echo "/usr/local/opt/libomp")
    OPENMP_FLAGS = -Xpreprocessor -fopenmp -I$(LIBOMP_PREFIX)/include
    OPENMP_LDFLAGS = -L$(LIBOMP_PREFIX)/lib -lomp
  else
    IS_CLANG := $(shell $(CC) --version 2>/dev/null | grep -qi clang && echo 1)
    OPENMP_FLAGS = -fopenmp
    OPENMP_LDFLAGS = -fopenmp
    ifdef LIBOMP_PREFIX
      OPENMP_FLAGS += -I$(LIBOMP_PREFIX)/include
      OPENMP_LDFLAGS += -L$(LIBOMP_PREFIX)/lib
    endif
    ifeq ($(IS_CLANG),1)
      OPENMP_LDFLAGS += -lomp
    endif
  endif
else
  OPENMP_FLAGS =
  OPENMP_LDFLAGS =
endif

# Release build (set RELEASE=1)
ifdef RELEASE
  ifdef PORTABLE
    CFLAGS += -O3 -DNDEBUG
  else
    CFLAGS += -O3 -march=native -DNDEBUG
  endif
endif

# Platform-specific whole-archive flags for format modules
ifeq ($(UNAME_S),Darwin)
  LIBS_MCNP = -Wl,-force_load,$(BIN_DIR)/libalea_mcnp.a
  LIBS_OPENMC = -Wl,-force_load,$(BIN_DIR)/libalea_openmc.a
  LIBS_SERPENT = -Wl,-force_load,$(BIN_DIR)/libalea_serpent.a
else
  LIBS_MCNP = -Wl,--whole-archive $(BIN_DIR)/libalea_mcnp.a -Wl,--no-whole-archive
  LIBS_OPENMC = -Wl,--whole-archive $(BIN_DIR)/libalea_openmc.a -Wl,--no-whole-archive
  LIBS_SERPENT = -Wl,--whole-archive $(BIN_DIR)/libalea_serpent.a -Wl,--no-whole-archive
endif

# All libraries
LIBS_ALL = $(LIBS_MCNP) $(LIBS_OPENMC) $(LIBS_SERPENT) $(BIN_DIR)/libalea.a
LDFLAGS += $(OPENMP_LDFLAGS)

# Windows: link psapi for GetProcessMemoryInfo in large_model_probe
ifdef WINDOWS_GNU
  PROBE_LIBS = -lpsapi
else ifeq ($(OS),Windows_NT)
  PROBE_LIBS = -lpsapi
else
  PROBE_LIBS =
endif

TOOLS = $(BIN_DIR)/mc_convert$(EXEEXT) $(BIN_DIR)/mc_plotter$(EXEEXT) $(BIN_DIR)/mc_plotter_experimental$(EXEEXT) $(BIN_DIR)/nuc_plot$(EXEEXT) $(BIN_DIR)/large_model_probe$(EXEEXT)

all: $(TOOLS)

$(BIN_DIR)/mc_convert$(EXEEXT): mc_convert.c
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $< $(LIBS_ALL) $(LDFLAGS)

$(BIN_DIR)/mc_plotter$(EXEEXT): mc_plotter.c
	$(CC) $(CFLAGS) $(OPENMP_FLAGS) -O2 $(INCLUDES) -o $@ $< $(LIBS_ALL) $(LDFLAGS)

$(BIN_DIR)/mc_plotter_experimental$(EXEEXT): mc_plotter_experimental.c
	$(CC) $(CFLAGS) $(OPENMP_FLAGS) -O2 $(INCLUDES) -o $@ $< $(LIBS_ALL) $(LDFLAGS)

$(BIN_DIR)/nuc_plot$(EXEEXT): nuc_plot.c
	$(CC) $(CFLAGS) $(INCLUDES) -I../src -o $@ $< $(BIN_DIR)/libalea_nucdata.a $(BIN_DIR)/libalea.a $(LDFLAGS)

$(BIN_DIR)/large_model_probe$(EXEEXT): large_model_probe.c
	$(CC) $(CFLAGS) $(INCLUDES) -I../src -o $@ $< $(LIBS_ALL) $(LDFLAGS) $(PROBE_LIBS)

clean:
	rm -f $(TOOLS)

.PHONY: all clean
