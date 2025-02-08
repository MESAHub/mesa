include ../make/defaults-module.mk

# Build

MODULE_NAME := turb
SRCS := public/turb.f90 \
        private/mlt.f90 \
        private/tdc_support.f90 \
        private/tdc.f90 \
        private/thermohaline.f90 \
        private/semiconvection.f90
SRCS_CHECK := test/src/test_turb.f90
INTERNAL_DEPENDS_ON := const utils auto_diff num chem
EXTERNAL_DEPENDS_ON :=
BINTYPE := static-lib

# Testing

CHECK_RESULTS_GOLDEN := test/test_output

# Install

MODULES := turb.mod
INSTALL_INCLUDES :=

include $(MAKE_DIR)/Makefile
