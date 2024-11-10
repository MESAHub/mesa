include ../make/defaults-module.mk

# Build

MODULE_NAME := colors
SRCS := public/colors_def.f90 \
        public/colors_lib.f90 \
        private/mod_colors.f90
SRCS_CHECK := test/src/test_colors.f90
INTERNAL_DEPENDS_ON := const utils math
EXTERNAL_DEPENDS_ON :=
BINTYPE := static-lib

# Testing

CHECK_RESULTS_GOLDEN := test/test_output

# Install

MODULES := colors_def.mod colors_lib.mod
INSTALL_INCLUDES :=
INSTALL_COMMANDS += install-colors-data

include $(MAKE_DIR)/Makefile

install-colors-data:
	mkdir -p ../data/colors_data && cp data/lcb98cor.dat data/blackbody_johnson.dat ../data/colors_data
