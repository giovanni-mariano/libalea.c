# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0

# SDL stays outside the library and scripting bindings. Opt in for viewers.
USE_SDL ?= 0
FETCH_SDL ?= 1
PKG_CONFIG ?= pkg-config
SDL2_CONFIG ?= sdl2-config
PYTHON ?= python3
SDL_JOBS ?= 4
SDL_TOOLS_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
SDL_LOCAL_DIR := $(abspath $(SDL_TOOLS_DIR)/../build/deps/SDL2-2.32.10)
ifeq ($(USE_SDL),1)
  ifeq ($(shell $(PKG_CONFIG) --exists sdl2 2>/dev/null && echo yes),yes)
    SDL_CFLAGS = -DALEA_USE_SDL $(shell $(PKG_CONFIG) --cflags sdl2)
    SDL_LIBS = $(shell $(PKG_CONFIG) --libs sdl2)
  else ifeq ($(shell $(SDL2_CONFIG) --version >/dev/null 2>&1 && echo yes),yes)
    SDL_CFLAGS = -DALEA_USE_SDL $(shell $(SDL2_CONFIG) --cflags)
    SDL_LIBS = $(shell $(SDL2_CONFIG) --libs)
  else ifeq ($(FETCH_SDL),1)
    SDL_LOCAL_CONFIG := $(SDL_LOCAL_DIR)/install/bin/sdl2-config
    SDL_CFLAGS = -DALEA_USE_SDL $(shell "$(SDL_LOCAL_CONFIG)" --cflags)
    SDL_LIBS = $(shell "$(SDL_LOCAL_CONFIG)" --static-libs)

    .PHONY: sdl-local
    sdl-local:
	$(PYTHON) "$(SDL_TOOLS_DIR)fetch_sdl.py" --directory "$(SDL_LOCAL_DIR)" --compiler "$(CC)" --jobs "$(SDL_JOBS)"

    viewer-config: sdl-local
  else
    $(error USE_SDL=1 with FETCH_SDL=0 requires SDL2 development files discoverable through pkg-config or sdl2-config)
  endif
endif

# Relink viewers when switching USE_SDL without requiring a clean build.
.PHONY: viewer-config
viewer-config:
