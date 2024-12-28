include ../make/defaults-module.mk

# Build

MODULE_NAME := neu
SRCS := public/neu_def.f90 \
        public/neu_lib.f90 \
        private/mod_neu.f90
SRCS_CHECK := test/src/test_neu.f90
INTERNAL_DEPENDS_ON := const utils math
EXTERNAL_DEPENDS_ON :=
BINTYPE := static-lib

# Testing

CHECK_RESULTS_GOLDEN := test/test_output

# Install

MODULES := neu_def.mod neu_lib.mod
INSTALL_INCLUDES :=

include $(MAKE_DIR)/Makefile
