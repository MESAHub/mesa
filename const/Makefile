include ../make/defaults-module.mk

# Build

MODULE_NAME := const
SRCS := public/const_def.f90 public/const_lib.f90
SRCS_CHECK := test/src/test_const.f90
INTERNAL_DEPENDS_ON :=
EXTERNAL_DEPENDS_ON :=
BINTYPE := static-lib

# Testing

CHECK_RESULTS_GOLDEN := test/test_output

# Install

MODULES := const_def.mod const_lib.mod
INSTALL_INCLUDES :=

include $(MAKE_DIR)/Makefile
