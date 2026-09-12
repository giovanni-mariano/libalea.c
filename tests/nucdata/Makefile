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
EPRDATA12_URL = https://nucleardata.lanl.gov/lib/eprdata12.zip
EPRDATA12_SHA512 = e39855622b08f3d6bb930c664b1020cef30add6dbfde654ea5c556f0ee8cedf81a459c3beb1805d6fba670c21d87acdee612f64caf294f4b2e014b965bce8a7c
EPRDATA14_URL = https://nucleardata.lanl.gov/lib/eprdata14.zip
EPRDATA14_SHA512 = 701445be628757d1bbd5de5773c882792c44b30a476dd2cd47615d77e43c39bbcd0c0cca675dedc850660a5343c57266af2927042f7dbd8fd005ee9607829a66
EPRDATA25_URL = https://nucleardata.lanl.gov/lib/eprdata25.zip
EPRDATA25_SHA512 = 25f0b63f27c534d8376e285634e3ac0dd6a3e9a670b2888eb19d9ec8b1e08988da17e6b5b4f3acbde21e48e137a39e3fea57a1a2fcabd912304a82e385a4d561
MCPLIB84_URL = https://nucleardata.lanl.gov/lib/mcplib84.tgz
ENDFB71_TSL_URL = https://www.nndc.bnl.gov/endf-b7.1/aceFiles/ENDF-B-VII.1-tsl.tar.gz
ENDFB71_TSL_SHA256 = 601db93498b51aa5a44800ee7ae7d439ac2bfb3ad36d7d1eb7e4a7be28968bc7
NJOY_TEST74_URL = https://raw.githubusercontent.com/njoy/NJOY2016/2016.79/tests/74/referenceTape71
NJOY_TEST74_SHA256 = 624147b06fc517cc7a5227ceac834a3febf107d00c1c5e86ecdddcfea51cd2d2
NJOY_TEST07_URL = https://raw.githubusercontent.com/njoy/NJOY2016/2016.79/tests/07/referenceTape26
NJOY_TEST07_SHA256 = 4459a2c249f064888aa7d38ecffbcae39ebea4170678ae412d0b71a4592d84c8
NJOY_TEST08_URL = https://raw.githubusercontent.com/njoy/NJOY2016/2016.79/tests/08/referenceTape25
NJOY_TEST08_SHA256 = b969e801a59924f526d4f032540dcb3a9ad912214bf4bf4a7c47ec00bdf4c56b
NJOY_TEST59_URL = https://raw.githubusercontent.com/njoy/NJOY2016/2016.79/tests/59/referenceTape27
NJOY_TEST59_SHA256 = c3a6e20ad375227a65b861c9bdee6e70fbde34f23b484721d7054e92cd7b4127
NJOY_TEST25_URL = https://raw.githubusercontent.com/njoy/NJOY2016/2016.79/tests/25/referenceTape70
NJOY_TEST25_SHA256 = 0e2bcca8e078e55d135391904a3cf851e3170f6b6d8907fb7e38964f75996e98

# Data directories
FENDL32c_DIR = fendl-FENDL-3.2c-neutron-ace
ENDFB80_DIR = Lib80x
EPRDATA12_DIR = eprdata12
EPRDATA14_DIR = eprdata14
EPRDATA25_DIR = eprdata25
MCPLIB84_DIR = mcplib84
ENDFB71_TSL_DIR = endfb71-tsl
NJOY_TEST74_DIR = njoy-test74
NJOY_TEST07_DIR = njoy-test07
NJOY_TEST08_DIR = njoy-test08
NJOY_TEST59_DIR = njoy-test59
NJOY_TEST25_DIR = njoy-test25

# Test binaries
TESTS = test_endfb80 test_fendl32c test_eprdata12 test_eprdata14 test_eprdata25 test_photon_xs test_nucdata_integration test_thermal test_photon_production test_production_photon_angles test_production_urr test_urr_matrix test_temperature_mix test_photoatomic_relaxation
BENCH = bench_load

# ============================================================================
# Targets
# ============================================================================

.PHONY: all test test-inventory-required test-photon-production-required test-production-photon-angles-required test-production-urr-required test-urr-matrix-required test-photoatomic-relaxation-required test-temperature-mix-required test-thermal-required clean data data-fendl32c data-endfb80 data-eprdata12 data-eprdata14 data-eprdata25 data-mcplib84 data-endfb71-tsl data-njoy-test74 data-njoy-test07 data-njoy-test08 data-njoy-test25 data-njoy-test59 help

all: $(TESTS)

help:
	@echo "Nuclear data test targets:"
	@echo ""
	@echo "  make data           - Download all nuclear data libraries"
	@echo "  make data-fendl32c  - Download FENDL-3.2c ACE"
	@echo "  make data-endfb80   - Download ENDF/B-VIII.0 Lib80x"
	@echo "  make data-eprdata12 - Download eprdata12"
	@echo "  make data-eprdata14 - Download eprdata14"
	@echo "  make data-eprdata25 - Download eprdata25"
	@echo "  make data-mcplib84  - Download mcplib84"
	@echo "  make data-endfb71-tsl - Download ENDF/B-VII.1 thermal ACE"
	@echo "  make data-njoy-test74 - Download NJOY IFENG=2 thermal ACE"
	@echo "  make data-njoy-test07 - Download NJOY photon-production ACE"
	@echo "  make data-njoy-test08 - Download NJOY MF=16 photon-production ACE"
	@echo "  make data-njoy-test59 - Download NJOY photoatomic relaxation ACE"
	@echo "  make data-njoy-test25 - Download NJOY light-water thermal ACE"
	@echo "  make test-photon-production-required - Download and run pinned photon tests"
	@echo "  make test-production-photon-angles-required - Run production anisotropy benchmark"
	@echo "  make test-production-urr-required - Run production U-238 URR benchmark"
	@echo "  make test-urr-matrix-required - Scan production URR representations"
	@echo "  make test-photoatomic-relaxation-required - Validate fluorescence metadata"
	@echo "  make test-temperature-mix-required - Run real resonance-table mixing tests"
	@echo "  make test-thermal-required - Download and run pinned continuous thermal tests"
	@echo "  make test-inventory-required - Inventory representative pinned tables"
	@echo ""
	@echo "  make all            - Build all test binaries"
	@echo "  make test           - Build and run all tests"
	@echo "  make clean          - Remove test binaries"

# ============================================================================
# Data download
# ============================================================================

data: data-fendl32c data-eprdata12 data-eprdata14 data-eprdata25 data-mcplib84 data-endfb71-tsl data-njoy-test74 data-njoy-test07 data-njoy-test08 data-njoy-test25 data-njoy-test59

data-fendl32c: $(FENDL32c_DIR)
data-endfb80: $(ENDFB80_DIR)
data-eprdata12: $(EPRDATA12_DIR)
data-eprdata14: $(EPRDATA14_DIR)
data-eprdata25: $(EPRDATA25_DIR)
data-mcplib84: $(MCPLIB84_DIR)
data-endfb71-tsl: $(ENDFB71_TSL_DIR)
data-njoy-test74: $(NJOY_TEST74_DIR)/xsdir
data-njoy-test07: $(NJOY_TEST07_DIR)/xsdir
data-njoy-test08: $(NJOY_TEST08_DIR)/xsdir
data-njoy-test59: $(NJOY_TEST59_DIR)/xsdir
data-njoy-test25: $(NJOY_TEST25_DIR)/xsdir

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

$(EPRDATA12_DIR):
	@echo "Downloading eprdata12..."
	curl -fSL -o eprdata12.zip $(EPRDATA12_URL)
	echo "$(EPRDATA12_SHA512)  eprdata12.zip" | sha512sum -c -
	unzip -q eprdata12.zip
	rm -f eprdata12.zip
	@echo "Done: $(EPRDATA12_DIR)/"

$(EPRDATA14_DIR):
	@echo "Downloading eprdata14..."
	curl -fSL -o eprdata14.zip $(EPRDATA14_URL)
	echo "$(EPRDATA14_SHA512)  eprdata14.zip" | sha512sum -c -
	unzip -q eprdata14.zip
	rm -f eprdata14.zip
	@echo "Done: $(EPRDATA14_DIR)/"

$(EPRDATA25_DIR):
	@echo "Downloading eprdata25..."
	curl -fSL -o eprdata25.zip $(EPRDATA25_URL)
	echo "$(EPRDATA25_SHA512)  eprdata25.zip" | sha512sum -c -
	unzip -q eprdata25.zip
	rm -f eprdata25.zip
	@echo "Done: $(EPRDATA25_DIR)/"

$(MCPLIB84_DIR):
	@echo "Downloading mcplib84..."
	curl -fSL -o mcplib84.tgz $(MCPLIB84_URL)
	tar xzf mcplib84.tgz
	rm -f mcplib84.tgz
	@echo "Done: $(MCPLIB84_DIR)/"

$(ENDFB71_TSL_DIR):
	@echo "Downloading ENDF/B-VII.1 thermal ACE..."
	curl -fSL -o endfb71-tsl.tar.gz $(ENDFB71_TSL_URL)
	echo "$(ENDFB71_TSL_SHA256)  endfb71-tsl.tar.gz" | sha256sum -c -
	mkdir -p $(ENDFB71_TSL_DIR)
	tar xzf endfb71-tsl.tar.gz -C $(ENDFB71_TSL_DIR)
	rm -f endfb71-tsl.tar.gz
	@echo "Done: $(ENDFB71_TSL_DIR)/"

$(NJOY_TEST74_DIR)/xsdir:
	@echo "Downloading NJOY2016 test 74 continuous thermal ACE..."
	mkdir -p $(NJOY_TEST74_DIR)
	curl -fSL -o $(NJOY_TEST74_DIR)/hzrh.10t $(NJOY_TEST74_URL)
	echo "$(NJOY_TEST74_SHA256)  $(NJOY_TEST74_DIR)/hzrh.10t" | sha256sum -c -
	printf '%s\n' 'directory' 'hzrh.10t 0.999167 hzrh.10t 0 1 1 2295813 0 0 2.551E-08' > $(NJOY_TEST74_DIR)/xsdir
	@echo "Done: $(NJOY_TEST74_DIR)/"

$(NJOY_TEST07_DIR)/xsdir:
	@echo "Downloading NJOY2016 test 07 U-235 photon-production ACE..."
	mkdir -p $(NJOY_TEST07_DIR)
	curl -fSL -o $(NJOY_TEST07_DIR)/u235.00c $(NJOY_TEST07_URL)
	echo "$(NJOY_TEST07_SHA256)  $(NJOY_TEST07_DIR)/u235.00c" | sha256sum -c -
	printf '%s\n' 'directory' '92235.00c 233.025000 u235.00c 0 1 1 56381 0 0 2.5852E-08' > $(NJOY_TEST07_DIR)/xsdir
	@echo "Done: $(NJOY_TEST07_DIR)/"

$(NJOY_TEST08_DIR)/xsdir:
	@echo "Downloading NJOY2016 test 08 Ni-61 MF=16 photon-production ACE..."
	mkdir -p $(NJOY_TEST08_DIR)
	curl -fSL -o $(NJOY_TEST08_DIR)/ni61.00c $(NJOY_TEST08_URL)
	echo "$(NJOY_TEST08_SHA256)  $(NJOY_TEST08_DIR)/ni61.00c" | sha256sum -c -
	printf '%s\n' 'directory' '28061.00c 60.408000 ni61.00c 0 1 1 89227 0 0 2.5852E-08' > $(NJOY_TEST08_DIR)/xsdir
	@echo "Done: $(NJOY_TEST08_DIR)/"

$(NJOY_TEST59_DIR)/xsdir:
	@echo "Downloading NJOY2016 test 59 photoatomic ACE..."
	mkdir -p $(NJOY_TEST59_DIR)
	curl -fSL -o $(NJOY_TEST59_DIR)/u.31p $(NJOY_TEST59_URL)
	echo "$(NJOY_TEST59_SHA256)  $(NJOY_TEST59_DIR)/u.31p" | sha256sum -c -
	printf '%s\n' 'directory' '92000.31p 235.984000 u.31p 0 1 1 71807 0 0 0.0' > $(NJOY_TEST59_DIR)/xsdir
	@echo "Done: $(NJOY_TEST59_DIR)/"

$(NJOY_TEST25_DIR)/xsdir:
	@echo "Downloading NJOY2016 test 25 light-water thermal ACE..."
	mkdir -p $(NJOY_TEST25_DIR)
	curl -fSL -o $(NJOY_TEST25_DIR)/lwtr.10t $(NJOY_TEST25_URL)
	echo "$(NJOY_TEST25_SHA256)  $(NJOY_TEST25_DIR)/lwtr.10t" | sha256sum -c -
	printf '%s\n' 'directory' 'lwtr.10t 0.999167 lwtr.10t 0 1 1 49375 0 0 2.5852E-08' > $(NJOY_TEST25_DIR)/xsdir
	@echo "Done: $(NJOY_TEST25_DIR)/"

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

test_eprdata12: test_eprdata14.c $(LIBS)
	$(CC) $(CFLAGS) $(INCLUDES) \
		-DEPRDATA_XSDIR='"eprdata12/xsdir"' \
		-DEPRDATA_SUFFIX='".12p"' -DEPRDATA_LABEL='"eprdata12"' \
		-DEPRDATA_FORMAT=1 -o $@ $< $(LIBS)

test_eprdata25: test_eprdata14.c $(LIBS)
	$(CC) $(CFLAGS) $(INCLUDES) \
		-DEPRDATA_XSDIR='"eprdata25/xsdir"' \
		-DEPRDATA_SUFFIX='".25p"' -DEPRDATA_LABEL='"eprdata25"' \
		-o $@ $< $(LIBS)

test_photon_xs: test_photon_xs.c $(LIBS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $< $(LIBS)

test_nucdata_integration: test_nucdata_integration.c $(LIBS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $< $(LIBS)

test_thermal: test_thermal.c $(LIBS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $< $(LIBS)

test_photon_production: test_photon_production.c $(LIBS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $< $(LIBS)

test_production_photon_angles: test_production_photon_angles.c $(LIBS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $< $(LIBS)

test_production_urr: test_production_urr.c $(LIBS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $< $(LIBS)

test_urr_matrix: test_urr_matrix.c $(LIBS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $< $(LIBS)

test_temperature_mix: test_temperature_mix.c $(LIBS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $< $(LIBS)

test_photoatomic_relaxation: test_photoatomic_relaxation.c $(LIBS)
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
		echo "No nuclear data test failures."; \
	else \
		echo "$$failed test(s) failed"; \
		exit 1; \
	fi

test-photon-production-required: data-njoy-test07 data-njoy-test08 test_photon_production ../../bin/alea
	./test_photon_production
	../../bin/alea test_photon_lua.lua

test-production-photon-angles-required: test_production_photon_angles
	@test -n "$(ALEA_PRODUCTION_NEUTRON_XSDIR)" || \
		(echo "ALEA_PRODUCTION_NEUTRON_XSDIR is required"; exit 1)
	ALEA_PRODUCTION_NEUTRON_XSDIR="$(ALEA_PRODUCTION_NEUTRON_XSDIR)" \
		./test_production_photon_angles

test-production-urr-required: test_production_urr
	@test -n "$(ALEA_ENDFB80_XSDIR)" || \
		(echo "ALEA_ENDFB80_XSDIR is required"; exit 1)
	ALEA_ENDFB80_XSDIR="$(ALEA_ENDFB80_XSDIR)" ./test_production_urr

test-urr-matrix-required: test_urr_matrix
	@test -n "$(ALEA_ENDFB80_XSDIR)" || \
		(echo "ALEA_ENDFB80_XSDIR is required"; exit 1)
	@test -n "$(ALEA_ENDFB71_XSDIR)" || \
		(echo "ALEA_ENDFB71_XSDIR is required"; exit 1)
	ALEA_ENDFB80_XSDIR="$(ALEA_ENDFB80_XSDIR)" \
	ALEA_ENDFB71_XSDIR="$(ALEA_ENDFB71_XSDIR)" ./test_urr_matrix

test-temperature-mix-required: data-njoy-test07 test_temperature_mix
	./test_temperature_mix

test-photoatomic-relaxation-required: data-njoy-test59 test_photoatomic_relaxation
	./test_photoatomic_relaxation

../../bin/nuc_inventory: ../../tools/nuc_inventory.c $(LIBS)
	$(MAKE) -C ../../tools ../bin/nuc_inventory CC=$(CC)

test-inventory-required: data-njoy-test07 data-njoy-test08 data-njoy-test25 data-njoy-test59 ../../bin/nuc_inventory
	./test_inventory.sh

../../bin/alea:
	$(MAKE) -C ../.. cli

test-thermal-required: data-endfb71-tsl data-njoy-test74 data-njoy-test25 test_thermal
	ALEA_NUCDATA_REQUIRED=1 ./test_thermal endfb71_graphite
	ALEA_NUCDATA_REQUIRED=1 ./test_thermal njoy74
	ALEA_NUCDATA_REQUIRED=1 ./test_thermal njoy25

clean:
	rm -f $(TESTS) $(BENCH)
