# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
# SPDX-License-Identifier: MPL-2.0

VIEW_PROJECT_ROOT := $(abspath $(dir $(lastword $(MAKEFILE_LIST)))/..)
USE_METAL ?= 0
ifeq ($(USE_METAL),1)
  ifneq ($(UNAME_S),Darwin)
    $(error USE_METAL=1 requires macOS)
  endif
  VIEW_COMPUTE_CFLAGS = -DALEA_VIEW_METAL
  VIEW_COMPUTE_LIBS = $(BIN_DIR)/compute/metal/libalea_compute.a -framework Foundation -framework Metal
  VIEW_COMPUTE_DEP = viewer-compute
.PHONY: viewer-compute
viewer-compute:
	$(MAKE) -C "$(VIEW_PROJECT_ROOT)" compute USE_METAL=1
endif
