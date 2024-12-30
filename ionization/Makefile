include ../make/defaults-module.mk

# Build

MODULE_NAME := ionization
SRCS := public/ionization_def.f90 \
        public/ionization_lib.f90 \
        private/ion_tables_load.f90 \
        private/ion_tables_eval.f90 \
        private/ion_table_plot.f90 \
        private/mod_ionization.f90
SRCS_CHECK := test/src/test_ionization.f
INTERNAL_DEPENDS_ON := const utils math interp_2d chem
EXTERNAL_DEPENDS_ON :=
BINTYPE := static-lib

# Testing

CHECK_RESULTS_GOLDEN := test/test_output

# Install

MODULES := ionization_def.mod ionization_lib.mod
INSTALL_INCLUDES :=
INSTALL_COMMANDS += install-ionization-data

include $(MAKE_DIR)/Makefile

install-ionization-data:
	if [ ! -r ../data/ionozation_data ]; then \
		tar -xJf test/ionization_data.tar.xz -C ../data; \
	fi
