# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

# Tools Makefile

CC = gcc
CFLAGS ?= -Wall -Wextra -g -std=c11
# Public tools compile against the installed header surface. Tools that inspect
# internals add -I../src on their own rule below.
INCLUDES = -I../include
BIN_DIR = ../bin
LDFLAGS ?= -lm

ifdef WINDOWS_GNU
  EXEEXT ?= .exe
else
  EXEEXT ?=
endif

# Platform-specific linker spelling below.
UNAME_S := $(shell uname -s)

# Threaded library support (pass USE_TINYPAR=1 when linking a threaded build)
ifdef USE_TINYPAR
  THREAD_FLAGS = -pthread
else
  THREAD_FLAGS =
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
LDFLAGS += $(THREAD_FLAGS)

# Windows: link psapi for GetProcessMemoryInfo in large_model_probe
ifdef WINDOWS_GNU
  PROBE_LIBS = -lpsapi
else ifeq ($(OS),Windows_NT)
  PROBE_LIBS = -lpsapi
else
  PROBE_LIBS =
endif

TOOLS = $(BIN_DIR)/mc_convert$(EXEEXT) $(BIN_DIR)/mc_plotter$(EXEEXT) $(BIN_DIR)/nuc_plot$(EXEEXT) $(BIN_DIR)/large_model_probe$(EXEEXT)

all: $(TOOLS)

$(BIN_DIR)/mc_convert$(EXEEXT): mc_convert.c
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $< $(LIBS_ALL) $(LDFLAGS)

$(BIN_DIR)/mc_plotter$(EXEEXT): mc_plotter.c
	$(CC) $(CFLAGS) $(THREAD_FLAGS) -O2 $(INCLUDES) -o $@ $< $(LIBS_ALL) $(LDFLAGS)

$(BIN_DIR)/nuc_plot$(EXEEXT): nuc_plot.c
	$(CC) $(CFLAGS) $(INCLUDES) -I../src -o $@ $< $(BIN_DIR)/libalea_nucdata.a $(BIN_DIR)/libalea.a $(LDFLAGS)

$(BIN_DIR)/large_model_probe$(EXEEXT): large_model_probe.c
	$(CC) $(CFLAGS) $(INCLUDES) -I../src -o $@ $< $(LIBS_ALL) $(LDFLAGS) $(PROBE_LIBS)

clean:
	rm -f $(TOOLS)

.PHONY: all clean
