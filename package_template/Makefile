include ../make/defaults-module.mk

# Build

MODULE_NAME := pkg
SRCS := public/pkg_def.f90 public/pkg_lib.f90 private/pkg_mod.f90
SRCS_CHECK := test/src/test_pkg.f90
INTERNAL_DEPENDS_ON :=
EXTERNAL_DEPENDS_ON :=
BINTYPE := static-lib

# Testing

CHECK_RESULTS_GOLDEN := test/test_output

# Install

MODULES := pkg_def.mod pkg_lib.mod
INSTALL_INCLUDES :=

include $(MAKE_DIR)/Makefile
