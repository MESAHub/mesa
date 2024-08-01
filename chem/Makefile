include ../make/defaults-module.mk

# Build

MODULE_NAME := chem
SRCS := public/chem_def.f90 \
        public/chem_lib.f90 \
        private/chem_isos_io.f90 \
        private/lodders_mod.f90 \
        private/nuclide_set_mod.f90
SRCS_CHECK := test/src/test_chem.f
INTERNAL_DEPENDS_ON := const utils math
EXTERNAL_DEPENDS_ON :=
BINTYPE := static-lib

# Testing

CHECK_RESULTS_GOLDEN := test/test_output

# Install

MODULES := chem_def.mod chem_lib.mod
INSTALL_INCLUDES :=
INSTALL_COMMANDS += install-chem-data

include $(MAKE_DIR)/Makefile

install-chem-data:
	mkdir -p ../data/chem_data && cp data/* ../data/chem_data
