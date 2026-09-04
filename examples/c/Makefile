# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

# Examples Makefile

CC = gcc
CFLAGS = -Wall -Wextra -g -std=c11
INCLUDES = -I../../include

# Threaded library support (pass USE_TINYPAR=1 when linking a threaded build)
ifdef USE_TINYPAR
  THREAD_FLAGS = -pthread
else
  THREAD_FLAGS =
endif

# Libraries
LIBS_CORE = ../../bin/libalea.a $(THREAD_FLAGS) -lm
LIBS_MCNP = -Wl,--whole-archive ../../bin/libalea_mcnp.a -Wl,--no-whole-archive
LIBS_OPENMC = -Wl,--whole-archive ../../bin/libalea_openmc.a -Wl,--no-whole-archive
LIBS_ALL = $(LIBS_MCNP) $(LIBS_OPENMC) $(LIBS_CORE)
LIBS_NUCDATA = ../../bin/libalea_nucdata.a $(LIBS_CORE)

EXAMPLES = basic mcnp_roundtrip mcnp_volume render3d void_demo

all: $(EXAMPLES)

basic: basic.c
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $< $(LIBS_ALL)

mcnp_roundtrip: mcnp_roundtrip.c
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $< $(LIBS_ALL)

mcnp_volume: mcnp_volume.c
	$(CC) $(CFLAGS) $(THREAD_FLAGS) $(INCLUDES) -o $@ $< $(LIBS_ALL)

render3d: render3d.c
	$(CC) $(CFLAGS) $(THREAD_FLAGS) -O2 $(INCLUDES) -o $@ $< $(LIBS_ALL)

void_demo: void_demo.c
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $< $(LIBS_ALL)

clean:
	rm -f $(EXAMPLES) *.inp *.ppm *.bmp *.png

.PHONY: all clean
