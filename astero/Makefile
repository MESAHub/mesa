include ../make/defaults-module.mk

# Build

MODULE_NAME := astero
SRCS := \
    public/astero_def.f90 \
    public/astero_lib.f90 \
    private/astero_support.f90 \
    private/extras_support.f90 \
    private/astero_run_support.f90

INTERNAL_DEPENDS_ON := utils star

ifeq ($(WITH_PGSTAR),yes)
SRCS += private/pgstar_astero_plots.f90
EXTERNAL_DEPENDS_ON += pgplot
else
SRCS += private/pgstar_astero_plots_stub.f90
endif

ifeq ($(WITH_GYRE),yes)
SRCS += private/gyre_support.f90
INTERNAL_DEPENDS_ON += gyre
else
SRCS += private/gyre_support_stub.f90
endif

ifeq ($(WITH_ADIPLS),yes)
SRCS += private/adipls_support.f90 private/adipls_support_procs.f90
INTERNAL_DEPENDS_ON += adipls
else
SRCS += private/adipls_support_stub.f90 private/adipls_support_procs_stub.f90
endif

SRCS_CHECK = test/src/test_astero.f
BINTYPE := static-lib
INCLUDE_DIRS := -Iprivate -Ipublic -Idefaults

# Testing

CHECK_RESULTS_GOLDEN := test/test_output

# Install

MODULES := astero_def.mod astero_lib.mod
INSTALL_INCLUDES := public/set_star_astero_procs.inc

include $(MAKE_DIR)/Makefile
