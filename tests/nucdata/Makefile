# SPDX-FileCopyrightText: 2026 Giovanni MARIANO
#
# SPDX-License-Identifier: MPL-2.0

# Nuclear data download & test Makefile

CC = gcc
CFLAGS = -Wall -Wextra -g -std=c11
INCLUDES = -I../../include -I../../tests
LIBS = ../../bin/libalea_nucdata.a ../../bin/libalea.a -lm

# Download URLs
FENDL32c_URL = https://www-nds.iaea.org/fendl/data/neutron/fendl-FENDL-3.2c-neutron-ace.zip
ENDFB80_URL = https://nucleardata.lanl.gov/lib/Lib80x.tgz
EPRDATA14_URL = https://nucleardata.lanl.gov/lib/eprdata14.tgz
MCPLIB84_URL = https://nucleardata.lanl.gov/lib/mcplib84.tgz

# Data directories
FENDL32c_DIR = fendl-FENDL-3.2c-neutron-ace
ENDFB80_DIR = Lib80x
EPRDATA14_DIR = eprdata14
MCPLIB84_DIR = mcplib84

# Test binaries
TESTS = test_endfb80 test_fendl32c test_eprdata14 test_photon_xs
BENCH = bench_load

# ============================================================================
# Targets
# ============================================================================

.PHONY: all test clean data data-fendl32c data-endfb80 data-eprdata14 data-mcplib84 help

all: $(TESTS)

help:
	@echo "Nuclear data test targets:"
	@echo ""
	@echo "  make data           - Download all nuclear data libraries"
	@echo "  make data-fendl32c  - Download FENDL-3.2c ACE"
	@echo "  make data-endfb80   - Download ENDF/B-VIII.0 Lib80x"
	@echo "  make data-eprdata14 - Download eprdata14"
	@echo "  make data-mcplib84  - Download mcplib84"
	@echo ""
	@echo "  make all            - Build all test binaries"
	@echo "  make test           - Build and run all tests"
	@echo "  make clean          - Remove test binaries"

# ============================================================================
# Data download
# ============================================================================

data: data-fendl32c data-eprdata14 data-mcplib84

data-fendl32c: $(FENDL32c_DIR)
data-endfb80: $(ENDFB80_DIR)
data-eprdata14: $(EPRDATA14_DIR)
data-mcplib84: $(MCPLIB84_DIR)

$(FENDL32c_DIR):
	@echo "Downloading FENDL-3.2c neutron ACE..."
	curl -fSL -o fendl-ace.zip $(FENDL32c_URL)
	unzip -q fendl-ace.zip
	rm -f fendl-ace.zip
	@echo "Done: $(FENDL32c_DIR)/"

$(ENDFB80_DIR):
	@echo "Downloading ENDF/B-VIII.0 Lib80x ACE (~7 GB)..."
	curl -fSL -o Lib80x.tgz $(ENDFB80_URL)
	tar xzf Lib80x.tgz
	rm -f Lib80x.tgz
	@echo "Done: $(ENDFB80_DIR)/"

$(EPRDATA14_DIR):
	@echo "Downloading eprdata14..."
	curl -fSL -o eprdata14.tgz $(EPRDATA14_URL)
	tar xzf eprdata14.tgz
	rm -f eprdata14.tgz
	@echo "Done: $(EPRDATA14_DIR)/"

$(MCPLIB84_DIR):
	@echo "Downloading mcplib84..."
	curl -fSL -o mcplib84.tgz $(MCPLIB84_URL)
	tar xzf mcplib84.tgz
	rm -f mcplib84.tgz
	@echo "Done: $(MCPLIB84_DIR)/"

# ============================================================================
# Test binaries
# ============================================================================

bench_load: bench_load.c $(LIBS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $< $(LIBS)

test_fendl32c: test_fendl32c.c $(LIBS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $< $(LIBS)

test_endfb80: test_endfb80.c $(LIBS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $< $(LIBS)

test_eprdata14: test_eprdata14.c $(LIBS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $< $(LIBS)

test_photon_xs: test_photon_xs.c $(LIBS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $< $(LIBS)

# ============================================================================
# Run tests
# ============================================================================

test: $(TESTS)
	@echo ""
	@echo "=== Nuclear Data Tests ==="
	@failed=0; \
	for t in $(TESTS); do \
		echo ""; \
		echo "Running $$t..."; \
		./$$t || failed=$$((failed + 1)); \
	done; \
	echo ""; \
	if [ $$failed -eq 0 ]; then \
		echo "All nuclear data tests passed!"; \
	else \
		echo "$$failed test(s) failed"; \
		exit 1; \
	fi

clean:
	rm -f $(TESTS) $(BENCH)
