# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

SHELL := /bin/sh

ROOT := $(abspath ../..)
PYTHON ?= python3
CONFIG := $(PYTHON) scripts/config.py value

PYTHON_TAG := $(shell $(CONFIG) python-tag)
PLATFORM_TAG ?= $(shell $(CONFIG) platform)
EXT_SUFFIX := $(shell $(CONFIG) ext-suffix)
PYTHON_INCLUDE := $(shell $(CONFIG) python-include)
NUMPY_INCLUDE := $(shell $(CONFIG) numpy-include)
PYTHON_CC := $(shell $(CONFIG) cc)

ifeq ($(origin CC),default)
CC := $(PYTHON_CC)
endif

BUILD_ROOT ?= build/$(PYTHON_TAG)-$(PLATFORM_TAG)
LIB_BUILD_DIR := $(abspath $(BUILD_ROOT)/libalea-obj)
LIB_BIN_DIR := $(abspath $(BUILD_ROOT)/libalea-bin)
LIBALEA := $(LIB_BIN_DIR)/libalea_full.a
PACKAGE_ROOT := $(abspath $(BUILD_ROOT)/package)
PACKAGE_DIR := $(PACKAGE_ROOT)/pyalea
PACKAGE_SOURCES := package/pyalea/__init__.py package/pyalea/__init__.pyi package/pyalea/py.typed
OBJECT := $(abspath $(BUILD_ROOT)/pyalea_binding.o)
MODULE := $(PACKAGE_DIR)/_alea$(EXT_SUFFIX)
DIST_DIR ?= $(abspath dist)

BINDING_SOURCES := $(wildcard src/*.c)
PUBLIC_HEADERS := $(wildcard $(ROOT)/include/*.h)
CPPFLAGS += -I$(ROOT)/include -I$(PYTHON_INCLUDE) -I$(NUMPY_INCLUDE)
BINDING_CFLAGS ?= -O2 -std=c11 -fPIC -pthread

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
MODULE_LDFLAGS := -bundle -undefined dynamic_lookup -Wl,-force_load,$(LIBALEA) -lm -pthread
else
MODULE_LDFLAGS := -shared -Wl,--whole-archive $(LIBALEA) -Wl,--no-whole-archive -lm -pthread
endif

.PHONY: all clean help package test test-runtime test-valgrind FORCE

all: $(MODULE)

$(LIBALEA): FORCE
	$(MAKE) -C $(ROOT) RELEASE=1 PORTABLE=1 TINYPAR_BACKEND=native \
		BUILD_DIR=$(LIB_BUILD_DIR) BIN_DIR=$(LIB_BIN_DIR) full

$(OBJECT): $(BINDING_SOURCES) $(PUBLIC_HEADERS) | $(BUILD_ROOT)
	$(CC) $(CPPFLAGS) $(BINDING_CFLAGS) -c src/pyalea_binding.c -o $@

$(MODULE): $(OBJECT) $(LIBALEA) $(PACKAGE_SOURCES) | $(PACKAGE_DIR)
	cp package/pyalea/__init__.py $(PACKAGE_DIR)/__init__.py
	cp package/pyalea/__init__.pyi $(PACKAGE_DIR)/__init__.pyi
	cp package/pyalea/py.typed $(PACKAGE_DIR)/py.typed
	$(CC) -o $@ $(OBJECT) $(MODULE_LDFLAGS)

$(BUILD_ROOT) $(PACKAGE_DIR):
	mkdir -p $@

test: all
	PYTHONPATH=$(PACKAGE_ROOT) $(PYTHON) -m pytest -p no:cacheprovider tests -v

# Test the already-built extension against the active NumPy runtime. This
# target deliberately has no build prerequisites: compatibility jobs must fail
# rather than silently recompiling against the runtime under test.
test-runtime:
	@test -f "$(MODULE)" || { echo "Missing prebuilt extension: $(MODULE)"; exit 1; }
	PYTHONPATH=$(PACKAGE_ROOT) $(PYTHON) -m pytest -p no:cacheprovider tests -v

test-valgrind: all
	PYTHONPATH=$(PACKAGE_ROOT) PYTHONMALLOC=malloc valgrind \
		--leak-check=full --show-leak-kinds=definite --num-callers=30 \
		--log-file=$(BUILD_ROOT)/valgrind.log \
		$(PYTHON) tests/leak_probe.py
	$(PYTHON) scripts/check_valgrind.py $(BUILD_ROOT)/valgrind.log

package: test
	$(PYTHON) scripts/package.py \
		--package-dir $(PACKAGE_DIR) --output-dir $(DIST_DIR) --root $(ROOT) \
		--platform-tag $(PLATFORM_TAG)
	PYTHONPATH=$(PACKAGE_ROOT) $(PYTHON) -m pytest -p no:cacheprovider tests -v

clean:
	rm -rf build dist .pytest_cache tests/__pycache__ scripts/__pycache__ package/pyalea/__pycache__

help:
	@echo "pyAlea native binding"
	@echo "  make                 Build pyalea._alea for PYTHON"
	@echo "  make test            Run standalone binding tests"
	@echo "  make test-runtime    Test the existing binary without rebuilding"
	@echo "  make test-valgrind   Check repeated binding calls for native leaks"
	@echo "  make package         Create a release archive and checksum"
	@echo "  make clean           Remove bindings/python build outputs"
	@echo "Variables: PYTHON, CC, BUILD_ROOT, DIST_DIR, BINDING_CFLAGS"

FORCE:
