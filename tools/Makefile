# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

# Tools Makefile

CC = gcc
CFLAGS = -Wall -Wextra -g -std=c11
INCLUDES = -I../include
BIN_DIR = ../bin

# Detect platform
UNAME_S := $(shell uname -s)

# OpenMP support (pass USE_OPENMP=1 to enable)
ifdef USE_OPENMP
  ifeq ($(UNAME_S),Darwin)
    LIBOMP_PREFIX := $(shell brew --prefix libomp 2>/dev/null || echo "/usr/local/opt/libomp")
    OPENMP_FLAGS = -Xpreprocessor -fopenmp -I$(LIBOMP_PREFIX)/include
    OPENMP_LDFLAGS = -L$(LIBOMP_PREFIX)/lib -lomp
  else
    OPENMP_FLAGS = -fopenmp
    OPENMP_LDFLAGS = -fopenmp
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

# Platform-specific whole-archive flags
ifeq ($(UNAME_S),Darwin)
  LIBS_MCNP = -Wl,-force_load,$(BIN_DIR)/libalea_mcnp.a
  LIBS_OPENMC = -Wl,-force_load,$(BIN_DIR)/libalea_openmc.a
  LIBS_RAYCAST = -Wl,-force_load,$(BIN_DIR)/libalea_raycast.a
else
  LIBS_MCNP = -Wl,--whole-archive $(BIN_DIR)/libalea_mcnp.a -Wl,--no-whole-archive
  LIBS_OPENMC = -Wl,--whole-archive $(BIN_DIR)/libalea_openmc.a -Wl,--no-whole-archive
  LIBS_RAYCAST = -Wl,--whole-archive $(BIN_DIR)/libalea_raycast.a -Wl,--no-whole-archive
endif

LIBS_CORE = $(BIN_DIR)/libalea.a $(OPENMP_LDFLAGS) -lm

TOOLS = $(BIN_DIR)/mc_convert $(BIN_DIR)/mc_plotter

all: $(TOOLS)

$(BIN_DIR)/mc_convert: mc_convert.c
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $< $(LIBS_MCNP) $(LIBS_OPENMC) $(LIBS_CORE)

$(BIN_DIR)/mc_plotter: mc_plotter.c
	$(CC) $(CFLAGS) $(OPENMP_FLAGS) -O2 $(INCLUDES) -o $@ $< \
		$(BIN_DIR)/libalea_slice.a $(LIBS_RAYCAST) \
		$(LIBS_MCNP) $(LIBS_OPENMC) $(LIBS_CORE)

clean:
	rm -f $(TOOLS)

.PHONY: all clean
