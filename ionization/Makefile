include ../make/defaults-module.mk

# Build

MODULE_NAME := ionization
SRCS := public/ionization_def.f90 \
        public/ionization_lib.f90 \
        private/ion_tables_load.f90 \
        private/ion_tables_eval.f90 \
        private/ion_table_plot.f90 \
        private/mod_ionization.f90
SRCS_CHECK := test/src/test_ionization.f90 test/src/test_ionization_support.f90
INTERNAL_DEPENDS_ON := const utils math interp_2d chem
EXTERNAL_DEPENDS_ON :=
BINTYPE := lib

# Testing

CHECK_RESULTS_GOLDEN := test/test_output

# Install

MODULES := ionization_def.mod ionization_lib.mod
INSTALL_INCLUDES :=
EXTRACT_DATA_COMMANDS += ../data/ionization_data/.extract-completed

include $(MAKE_DIR)/Makefile

../data/ionization_data/.extract-completed: test/ionization_data.tar.xz
	rm -rf ../data/ionization_data
	tar -xJf test/ionization_data.tar.xz -C ../data
	touch $@
